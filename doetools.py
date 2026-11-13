"""Small toolbox for designed experiments with several sizes of experimental unit.

Everything is done with projections (numpy only), so the classical multi-stratum
ANOVA and the residual checks need no R and no specialised package.

strata_anova     ANOVA table with one error line per stratum (balanced designs)
stratum_tests    simulation-calibrated residual checks for one error stratum
tukey_nonadditivity, boxcox_interval, neighbour_pairs, make_groupings
"""
import numpy as np
import pandas as pd
import patsy
from scipy import stats

NSIM = 2000


# ----------------------------------------------------------------- linear algebra
def colspace(X, tol=1e-9):
    """Orthonormal basis of the column space of X (rank-revealing)."""
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    return U[:, : int((s > tol * s.max()).sum())]


def residual_maker(X):
    """M = I - H and the residual degrees of freedom."""
    Q = colspace(X)
    return np.eye(X.shape[0]) - Q @ Q.T, X.shape[0] - Q.shape[1]


def unit_ids(df, cols):
    return df[cols].astype(str).agg("|".join, axis=1).values


def dummies(ids):
    return pd.get_dummies(ids).values.astype(float)


# ----------------------------------------------------------------- multi-stratum ANOVA
def strata_anova(df, response, fixed, chain, labels=None):
    """ANOVA for a balanced design with nested error strata.

    chain  list of column lists, coarse -> fine, e.g. [["block"], ["block","gen"], ["block","gen","trt"]].
           The observations themselves form the last (residual) stratum.
    Each treatment term is tested in the coarsest stratum whose units it is constant on.
    Returns a DataFrame: stratum, term, df, SS, MS, F, p.
    """
    y = df[response].astype(float).values
    n = len(y)
    dm = patsy.dmatrix(fixed, df)
    Xf = np.asarray(dm)
    slices = {k: v for k, v in dm.design_info.term_name_slices.items() if k != "Intercept"}
    ids = [unit_ids(df, s) for s in chain]
    D = [dummies(i) for i in ids]
    one = np.ones((n, 1))
    labels = labels or [":".join(s) for s in chain] + ["residual"]

    def ss_rank(X):
        Q = colspace(X)
        return float(((Q.T @ y) ** 2).sum()), Q.shape[1]

    def constant_within(cols, unit):
        g = pd.DataFrame(cols).groupby(unit)
        return bool(((g.max() - g.min()).abs().values < 1e-12).all())

    assign = {}
    for name, sl in slices.items():
        assign[name] = len(chain)
        for k in range(1, len(chain)):
            if constant_within(Xf[:, sl], ids[k]):
                assign[name] = k
                break
    rows = []
    for k in range(len(chain) + 1):
        X = D[k - 1] if k >= 1 else one
        ss0, r0 = ss_rank(X)
        block_rows = []
        for name, sl in slices.items():
            if assign[name] != k:
                continue
            X = np.hstack([X, Xf[:, sl]])
            ss1, r1 = ss_rank(X)
            block_rows.append([labels[k], name, r1 - r0, ss1 - ss0])
            ss0, r0 = ss1, r1
        if k < len(chain):
            ss1, r1 = ss_rank(np.hstack([X, D[k]]))
        else:
            ss1, r1 = float(y @ y), n
        err_df, err_ss = r1 - r0, ss1 - ss0
        for row in block_rows:
            ms = row[3] / row[2]
            F = ms / (err_ss / err_df)
            rows.append(row + [ms, F, float(stats.f.sf(F, row[2], err_df))])
        rows.append([labels[k], "Blocks" if k == 0 else f"Error ({k})", err_df, err_ss, err_ss / err_df, np.nan, np.nan])
    return pd.DataFrame(rows, columns=["stratum", "term", "df", "SS", "MS", "F", "p"])


# ----------------------------------------------------------------- residual statistics (vectorised over simulations)
def bf_stat(R, codes, k):
    """Brown-Forsythe (median Levene) F for each row of R."""
    Z = np.empty_like(R)
    for g in range(k):
        idx = codes == g
        Z[:, idx] = np.abs(R[:, idx] - np.median(R[:, idx], axis=1, keepdims=True))
    grand = Z.mean(1)
    ssb = np.zeros(R.shape[0])
    ssw = np.zeros(R.shape[0])
    for g in range(k):
        idx = codes == g
        zg = Z[:, idx]
        m = zg.mean(1)
        ssb += idx.sum() * (m - grand) ** 2
        ssw += ((zg - m[:, None]) ** 2).sum(1)
    # rounding makes structurally-zero statistics (mirror-image residual pairs) tie at 0 -> p = 1
    return np.round((ssb / (k - 1)) / (ssw / (R.shape[1] - k)), 9)


def var_ratio(r, codes, k):
    v = np.array([np.mean(r[codes == g] ** 2) for g in range(k)])
    return float(v.max() / max(v.min(), 1e-300))


def spearman_rows(A, B):
    ra = stats.rankdata(A, axis=1)
    rb = stats.rankdata(B, axis=1)
    ra -= ra.mean(1, keepdims=True)
    rb -= rb.mean(1, keepdims=True)
    return (ra * rb).sum(1) / np.sqrt((ra**2).sum(1) * (rb**2).sum(1))


def p_upper(obs, sims):
    return float((1 + np.sum(sims >= obs)) / (len(sims) + 1))


def p_lower(obs, sims):
    return float((1 + np.sum(sims <= obs)) / (len(sims) + 1))


def stratum_tests(r, M, df_res, fitted, groupings, pairs, rng, nsim=NSIM):
    """Simulation-calibrated checks for one error stratum.

    r = M y are the stratum residuals.  Under the model r = M e with e iid normal, and every
    statistic below is scale free, so pushing nsim standard-normal vectors through M gives the
    exact null distribution of each statistic for this particular design.
    """
    n = len(r)
    out = {"n_units": int(n), "df": int(df_res)}
    E = rng.standard_normal((nsim, n))
    R = E @ M
    sigma = np.sqrt(r @ r / df_res)
    out["sigma"] = float(sigma)
    # normality
    w_obs = stats.shapiro(r)[0]
    w_sim = np.array([stats.shapiro(x)[0] for x in R])
    out["shapiro_W"] = float(w_obs)
    out["p_normal"] = p_lower(w_obs, w_sim)
    # outliers
    lev = np.sqrt(np.clip(np.diag(M), 1e-12, None))
    t_obs = np.abs(r / lev) / sigma
    t_sim = np.abs(R / lev) / np.sqrt((R**2).sum(1, keepdims=True) / df_res)
    out["max_abs_student"] = float(t_obs.max())
    out["p_outlier"] = p_upper(t_obs.max(), t_sim.max(1))
    # equal variance, family-wise over all groupings (min-p)
    if groupings:
        obs = np.empty(len(groupings))
        sim = np.empty((nsim, len(groupings)))
        ratios = {}
        for j, (name, codes, k) in enumerate(groupings):
            obs[j] = bf_stat(r[None, :], codes, k)[0]
            sim[:, j] = bf_stat(R, codes, k)
            ratios[name] = var_ratio(r, codes, k)
        p_each = np.array([p_upper(obs[j], sim[:, j]) for j in range(len(groupings))])
        rank = (nsim - stats.rankdata(sim, axis=0, method="min") + 2) / (nsim + 1)
        out["p_var_family"] = p_lower(p_each.min(), rank.min(1))
        worst = int(np.argmin(p_each))
        out["var_worst_grouping"] = groupings[worst][0]
        out["var_worst_p"] = float(p_each[worst])
        out["var_worst_ratio"] = ratios[groupings[worst][0]]
        out["var_each"] = {g[0]: {"p": float(p), "max_min_var_ratio": ratios[g[0]]} for g, p in zip(groupings, p_each)}
    # mean-variance trend
    if fitted is not None:
        rho_obs = spearman_rows(np.abs(r)[None, :], fitted[None, :])[0]
        rho_sim = spearman_rows(np.abs(R), fitted[None, :] + sigma * (E - R))
        out["rho_abs_resid_vs_fitted"] = float(rho_obs)
        out["p_trend"] = p_upper(abs(rho_obs), np.abs(rho_sim))
    # spatial independence (neighbouring plots)
    if pairs is not None and len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        s_obs = (r[i] * r[j]).sum() / (r @ r)
        s_sim = (R[:, i] * R[:, j]).sum(1) / (R**2).sum(1)
        out["neighbour_stat"] = float(s_obs)
        out["neighbour_z"] = float((s_obs - s_sim.mean()) / s_sim.std())
        out["p_spatial"] = float(min(1.0, 2 * min(p_upper(s_obs, s_sim), p_lower(s_obs, s_sim))))
    return out


def tukey_nonadditivity(y, X):
    """Tukey's one-degree-of-freedom test: does fitted^2 explain the residuals?"""
    M, df = residual_maker(X)
    r = M @ y
    z = M @ ((y - r) ** 2)
    if z @ z < 1e-10 * max(1.0, (y @ y) ** 2):
        return None
    ss1 = (r @ z) ** 2 / (z @ z)
    F = ss1 / ((r @ r - ss1) / (df - 1))
    return float(stats.f.sf(F, 1, df - 1))


def boxcox_interval(y, X):
    """Box-Cox lambda (profile likelihood on a grid) and its 95% interval."""
    if np.any(y <= 0):
        return None
    M, _ = residual_maker(X)
    y = y / np.exp(np.mean(np.log(y)))  # geometric mean 1: the Jacobian term vanishes
    lam = np.round(np.arange(-2, 2.0001, 0.05), 2)
    ll = []
    for l in lam:
        z = M @ (np.log(y) if abs(l) < 1e-9 else (y**l - 1) / l)
        ll.append(-0.5 * len(y) * np.log(z @ z))
    ll = np.array(ll)
    keep = lam[ll >= ll.max() - stats.chi2.ppf(0.95, 1) / 2]
    return {"lambda_hat": float(lam[ll.argmax()]), "ci": [float(keep.min()), float(keep.max())]}


def neighbour_pairs(df):
    """Index pairs of plots that touch (same row & adjacent column, or same column & adjacent row)."""
    pos = {(r, c): i for i, (r, c) in enumerate(zip(df["row"], df["col"]))}
    out = []
    for (r, c), i in pos.items():
        for key in ((r, c + 1), (r + 1, c)):
            if key in pos:
                out.append((i, pos[key]))
    return np.array(out)


def make_groupings(df, factors, min_cell=4):
    g = []
    for f in factors:
        codes, lv = pd.factorize(df[f])
        if len(lv) > 1 and np.bincount(codes).min() >= min_cell:
            g.append((f, codes, len(lv)))
    return g
