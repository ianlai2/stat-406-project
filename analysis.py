"""STAT 406 Project Part 2 - analysis.

Data: Goulden (1939) wheat root-rot seed-treatment experiment (agridat::goulden.splitsplit).
Design: split-split-plot in 4 blocks, 160 plots.
    whole plot     (half a block, 8 in total)      variety: Marquis, Mindum
    subplot        (1/10 of a whole plot, 80)      seed dust (5) x application method (dry, wet)
    sub-subplot    (half a subplot, 160)           soil: inoculated, uninoculated
Response: plot yield.

Run:  python analysis.py
Writes every table and figure to outputs/ and the full console log to outputs/results.txt.

Sections of the log: 1 data, 2 ANOVA, 3 which levels differ, 4 assumption checks, 5 validation against an
equal-variance REML mixed model, 6 comparison method (REML mixed models with unequal variances, fitted by
unequal_variance.R), 7 how robust each significant effect is.  Section 6 runs the R script when Rscript is
installed and otherwise reads the results saved in outputs/.
"""
import itertools
import pathlib
import re
import shutil
import subprocess
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy import stats

from doetools import (boxcox_interval, dummies, make_groupings, neighbour_pairs, residual_maker,
                      strata_anova, stratum_tests, tukey_nonadditivity, unit_ids)

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)
SEED, NSIM = 406, 5000

# chart style: two validated categorical hues, neutral ink, recessive grid
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "lines.linewidth": 2, "lines.markersize": 7, "legend.frameon": False, "figure.dpi": 110,
})


class Tee:
    def __init__(self, path):
        self.f, self.out = open(path, "w", encoding="utf-8"), sys.stdout

    def write(self, s):
        self.f.write(s)
        self.out.write(s)

    def flush(self):
        self.f.flush()
        self.out.flush()


def header(title):
    print("\n" + "=" * 96 + f"\n{title}\n" + "=" * 96)


def fmt_p(p):
    return "<0.0001" if p < 1e-4 else f"{p:.4f}"


# ======================================================================================= 1. data
def load():
    d = pd.read_csv(HERE / "data" / "goulden_splitsplit.csv")
    d["trt"] = d["trt"].astype(str)
    d["gen"] = pd.Categorical(d["gen"], ["Marquis", "Mindum"])
    d["dust"] = pd.Categorical(d["dust"], ["Check", "CaCO3", "Ceresan", "DuBay", "Semesan"])
    d["dry"] = pd.Categorical(d["dry"], ["Dry", "Wet"])
    d["inoc"] = pd.Categorical(d["inoc"], ["Uninoc", "Inoc"])
    d["wp"] = unit_ids(d, ["block", "gen"])          # whole plot id
    d["sp"] = unit_ids(d, ["block", "gen", "trt"])   # subplot id
    return d.sort_values(["row", "col"]).reset_index(drop=True)


def describe(d):
    header("1. DATA AND DESIGN")
    print(f"observations: {len(d)}   missing values: {int(d.isna().sum().sum())}   yield range: {d['yield'].min()}-{d['yield'].max()}")
    print(f"blocks: {d.block.nunique()}   whole plots: {d.wp.nunique()}   subplots: {d.sp.nunique()}   sub-subplots: {len(d)}")
    cell = d.groupby(["gen", "dust", "dry", "inoc"], observed=True).size()
    print(f"treatment combinations: {len(cell)} (2 x 5 x 2 x 2), replicates per combination: min {cell.min()}, max {cell.max()}  -> balanced")
    ok = (d.groupby("sp")["col"].agg(lambda c: c.max() - c.min()) == 1).all() and (d.groupby("wp")["row"].nunique() == 1).all()
    print(f"layout check - each whole plot is one field row, each subplot is two adjacent columns: {bool(ok)}")
    print("\nmean yield by factor level")
    for f in ["gen", "dust", "dry", "inoc"]:
        print(f"  {f:5s} " + "   ".join(f"{k}: {v:.2f}" for k, v in d.groupby(f, observed=True)["yield"].mean().items()))


# ======================================================================================= 2. ANOVA
FIXED = "gen*dust*dry*inoc"
CHAIN = [["block"], ["block", "gen"], ["block", "gen", "trt"]]
STRATA = ["Blocks", "Whole plots", "Subplots", "Sub-subplots"]
# sums of squares printed in Goulden (1939), p. 157.  He pooled the 4-factor interaction (7.7, 4 df)
# into Error (3), so his Error (3) is 931.1 on 64 df; with that term fitted it is 923.4 on 60 df.
GOULDEN_SS = {"Blocks": 989.6, "gen": 3638.6, "Error (1)": 647.6, "dust": 987.6, "dry": 117.3, "dust:dry": 46.2,
              "gen:dust": 146.7, "gen:dry": 91.5, "gen:dust:dry": 148.1, "Error (2)": 1059.1, "inoc": 965.3,
              "gen:inoc": 0.3, "dust:inoc": 379.8, "dry:inoc": 68.9, "dust:dry:inoc": 25.8, "gen:dust:inoc": 119.4,
              "gen:dry:inoc": 3.9}


def anova(d):
    header("2. SPLIT-SPLIT-PLOT ANOVA (one error term per size of experimental unit)")
    a = strata_anova(d, "yield", FIXED, CHAIN, STRATA)
    a["Goulden_SS"] = a["term"].map(GOULDEN_SS)
    # Holm adjustment over the 15 treatment F tests in the table
    tested = a["p"].notna()
    a["p_holm"] = np.nan
    a.loc[tested, "p_holm"] = multipletests(a.loc[tested, "p"], method="holm")[1]
    show = a.copy()
    show["F"] = show["F"].map(lambda v: "" if np.isnan(v) else f"{v:.2f}")
    show["p"] = a["p"].map(lambda v: "" if np.isnan(v) else fmt_p(v))
    show["p_holm"] = a["p_holm"].map(lambda v: "" if np.isnan(v) else fmt_p(v))
    show["sig"] = a["p"].map(lambda v: "" if np.isnan(v) else ("***" if v < .001 else "**" if v < .01 else "*" if v < .05 else ""))
    show["Goulden_SS"] = show["Goulden_SS"].map(lambda v: "" if np.isnan(v) else f"{v:.1f}")
    print(show.to_string(index=False, formatters={"SS": "{:.1f}".format, "MS": "{:.2f}".format}))
    chk = a.dropna(subset=["Goulden_SS"])
    worst = (chk["SS"] - chk["Goulden_SS"]).abs().max()
    e3 = a.loc[a.term.isin(["Error (3)", "gen:dust:dry:inoc"]), "SS"].sum()
    print(f"\ncheck against the published table: largest |SS - Goulden| over {len(chk)} lines = {worst:.2f}; "
          f"Error (3) + 4-factor interaction = {e3:.1f} (Goulden: 931.1)")
    assert worst < 0.15 and abs(e3 - 931.1) < 0.15, "ANOVA does not reproduce the published table"
    held = a.loc[a["p_holm"] < 0.05, "term"].tolist()
    print(f"terms with p < 0.05: {int((a['p'] < 0.05).sum())} of {int(tested.sum())}; "
          f"still below 0.05 after the Holm adjustment: {', '.join(held)}")
    a.to_csv(OUT / "anova_table.csv", index=False)
    ms = a.set_index("term")["MS"]
    return a, {"blk": ms["Blocks"], 1: ms["Error (1)"], 2: ms["Error (2)"], 3: ms["Error (3)"]}, {1: 3, 2: 54, 3: 60}


# ======================================================================================= 5. validation
def mixed_model(d, a, ms):
    header("5. VALIDATION: LINEAR MIXED MODEL FITTED BY REML WITH EQUAL VARIANCES (statsmodels MixedLM)\n"
           "Balanced data, so agreement with the ANOVA is expected. Section 6 relaxes the equal-variance assumption.")
    f = "Q('yield') ~ C(gen, Sum)*C(dust, Sum)*C(dry, Sum)*C(inoc, Sum)"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        md = smf.mixedlm(f, d, groups="block", re_formula="1",
                         vc_formula={"wholeplot": "0 + C(gen)", "subplot": "0 + C(sp)"})
        # BFGS reaches the REML optimum; the default L-BFGS stops early on this model (lower likelihood)
        fit = md.fit(reml=True, method="bfgs", maxiter=5000)
    vc = dict(zip(fit.model.exog_vc.names, fit.vcomp))
    reml = {"block": float(fit.cov_re.iloc[0, 0]), "whole plot": vc["wholeplot"], "subplot": vc["subplot"], "residual": fit.scale}
    mom = {"block": (ms["blk"] - ms[1]) / 40, "whole plot": (ms[1] - ms[2]) / 20, "subplot": (ms[2] - ms[3]) / 2, "residual": ms[3]}
    print(f"converged: {fit.converged}")
    print("\nvariance components        REML     ANOVA (method of moments)")
    for k in reml:
        print(f"  {k:12s}        {reml[k]:8.2f}     {mom[k]:8.2f}")
    # Wald F tests for each term (sum-to-zero coding, so these are the usual factorial effects)
    b, V = fit.fe_params.values, fit.cov_params().iloc[: len(fit.fe_params), : len(fit.fe_params)].values
    term_of = [re.sub(r"C\((\w+), Sum\)\[S\.[^\]]+\]", r"\1", nm) for nm in fit.model.exog_names]
    rows = []
    aa = a.set_index("term")
    for term in dict.fromkeys(term_of):
        if term == "Intercept":
            continue
        idx = [i for i, t in enumerate(term_of) if t == term]
        W = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx])) / len(idx)
        rows.append([term, len(idx), W, aa.loc[term, "F"]])
    cmp_ = pd.DataFrame(rows, columns=["term", "df", "F_mixed_model", "F_anova"])
    cmp_["ratio"] = cmp_["F_mixed_model"] / cmp_["F_anova"]
    print("\nfixed-effect tests: mixed-model Wald F next to the ANOVA F")
    print(cmp_.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))
    print(f"\nlargest relative difference between the two sets of F statistics: {(cmp_['ratio'] - 1).abs().max():.2%}")
    cmp_.to_csv(OUT / "mixed_model_vs_anova.csv", index=False)
    pd.DataFrame({"REML": reml, "ANOVA": mom}).to_csv(OUT / "variance_components.csv")
    return reml


# ======================================================================================= 4. assumptions
def stratum_residuals(d):
    """Residuals for each error stratum, with the matrix that forms them."""
    y = d["yield"].astype(float).values
    Xf = np.asarray(patsy.dmatrix(FIXED, d))
    res = {}
    # sub-subplots: remove treatment effects and subplot means
    X3 = np.hstack([Xf, dummies(d["sp"].values)])
    M3, df3 = residual_maker(X3)
    res[3] = dict(r=M3 @ y, M=M3, df=df3, fitted=y - M3 @ y, frame=d, X=X3, y=y)
    # subplots: subplot means, remove treatment effects and whole-plot means
    for level, unit, parent in ((2, "sp", "wp"), (1, "wp", "block")):
        g = d.groupby(unit, sort=True)
        ybar = g["yield"].mean()
        Xbar = pd.DataFrame(Xf).groupby(d[unit].values).mean().loc[ybar.index].values
        first = g.first().loc[ybar.index].reset_index()
        Xk = np.hstack([Xbar, dummies(first[parent].values)])
        Mk, dfk = residual_maker(Xk)
        rk = Mk @ ybar.values
        res[level] = dict(r=rk, M=Mk, df=dfk, fitted=ybar.values - rk, frame=first, X=Xk, y=ybar.values)
    return res


def assumptions(d, aov):
    header("4. ASSUMPTION CHECKS (each error stratum separately; p-values calibrated with "
           f"{NSIM} simulations)")
    rng = np.random.default_rng(SEED)
    res = stratum_residuals(d)
    for a, b in itertools.combinations(["gen", "dust", "dry", "inoc"], 2):
        for k in (2, 3):
            fr = res[k]["frame"]
            fr[f"{a}:{b}"] = fr[a].astype(str) + "|" + fr[b].astype(str)
    pairs = [f"{a}:{b}" for a, b in itertools.combinations(["gen", "dust", "dry", "inoc"], 2)]
    g3 = make_groupings(res[3]["frame"], ["gen", "dust", "dry", "inoc", "block"] + pairs)
    g2 = make_groupings(res[2]["frame"], ["gen", "dust", "dry", "block", "gen:dust", "gen:dry", "dust:dry"])
    t3 = stratum_tests(res[3]["r"], res[3]["M"], res[3]["df"], res[3]["fitted"], g3, neighbour_pairs(d), rng, NSIM)
    t2 = stratum_tests(res[2]["r"], res[2]["M"], res[2]["df"], res[2]["fitted"], g2, None, rng, NSIM)
    t3["p_nonadditive"] = tukey_nonadditivity(res[3]["y"], res[3]["X"])
    t2["p_nonadditive"] = tukey_nonadditivity(res[2]["y"], res[2]["X"])
    bc = boxcox_interval(res[3]["y"], res[3]["X"])
    rows = []
    for name, t in (("Subplot error (54 df)", t2), ("Sub-subplot error (60 df)", t3)):
        rows += [
            [name, "Normality", f"Shapiro-Wilk W = {t['shapiro_W']:.3f}", t["p_normal"]],
            [name, "Equal variance", f"Brown-Forsythe over {len(t['var_each'])} groupings; worst: {t['var_worst_grouping']} "
                                     f"(max/min variance {t['var_worst_ratio']:.1f})", t["p_var_family"]],
            [name, "Variance vs mean", f"Spearman(|resid|, fitted) = {t['rho_abs_resid_vs_fitted']:+.3f}", t["p_trend"]],
            [name, "Outliers", f"largest |studentized residual| = {t['max_abs_student']:.2f}", t["p_outlier"]],
            [name, "Additivity", "Tukey one-df test", t["p_nonadditive"]],
        ]
    rows.append(["Sub-subplot error (60 df)", "Independence (field layout)",
                 f"neighbouring-plot residual correlation statistic z = {t3['neighbour_z']:+.2f}", t3["p_spatial"]])
    tab = pd.DataFrame(rows, columns=["stratum", "assumption", "statistic", "p"])
    tab["verdict"] = tab["p"].map(lambda p: "ok" if p >= 0.05 else ("borderline" if p >= 0.01 else "VIOLATED"))
    print(tab.to_string(index=False, formatters={"p": "{:.3f}".format}))
    tab.to_csv(OUT / "assumption_checks.csv", index=False)
    print("\nequal-variance detail, sub-subplot stratum (p for each grouping, max/min residual variance):")
    for k, v in t3["var_each"].items():
        print(f"  {k:12s} p = {v['p']:.3f}   ratio = {v['max_min_var_ratio']:.2f}")
    print("\nequal-variance detail, subplot stratum:")
    for k, v in t2["var_each"].items():
        print(f"  {k:12s} p = {v['p']:.3f}   ratio = {v['max_min_var_ratio']:.2f}")
    print(f"\nBox-Cox (sub-subplot stratum): lambda-hat = {bc['lambda_hat']}, 95% interval {bc['ci']}"
          "  (1 = no transformation)")
    # sensitivity: the interval sits just above 1, so repeat the ANOVA on yield^2 and compare decisions
    a2 = strata_anova(d.assign(y2=d["yield"].astype(float) ** 2), "y2", FIXED, CHAIN, STRATA)
    both = aov.merge(a2, on="term", suffixes=("", "_sq")).dropna(subset=["p"])
    flips = both[(both["p"] < 0.05) != (both["p_sq"] < 0.05)]
    kept = int(((both["p"] < 0.05) & (both["p_sq"] < 0.05)).sum())
    print(f"Sensitivity to scale: ANOVA repeated on yield^2 -> {kept} of the {int((both['p'] < 0.05).sum())} terms "
          "significant at 5% on the raw scale stay significant; "
          + ("no other term changes." if flips.empty else "decision changes: " + ", ".join(
              f"{t} (p = {p:.3f} raw, {q:.3f} squared)" for t, p, q in zip(flips.term, flips.p, flips.p_sq))))
    both[["stratum", "term", "df", "F", "p", "F_sq", "p_sq"]].to_csv(OUT / "sensitivity_yield_squared.csv", index=False)
    r1 = res[1]
    print(f"\nWhole-plot error has only {r1['df']} df (8 whole plots): too few for a formal test. Whole-plot residuals:")
    print("  " + "  ".join(f"{u}: {v:+.2f}" for u, v in zip(r1["frame"]["wp"], r1["r"])))
    print("Independence otherwise rests on the randomisation at each level, which the mixed model mirrors with\n"
          "random effects for blocks, whole plots and subplots.")
    return res, {2: t2, 3: t3}, rng


# ======================================================================================= 3. inference
def ci_row(label, est, se, df, mult=1):
    t = est / se
    p = min(1.0, 2 * stats.t.sf(abs(t), df) * mult)
    h = stats.t.ppf(1 - 0.025 / mult, df) * se
    return [label, est, se, df, est - h, est + h, p]


def letters(means, sig):
    """Compact letter display for means sorted high -> low; sig[i][j] True when i and j differ."""
    k = len(means)
    runs = []
    for i in range(k):
        j = i
        while j + 1 < k and not any(sig[a][j + 1] for a in range(i, j + 1)):
            j += 1
        if not any(a <= i and j <= b for a, b in runs):
            runs.append((i, j))
    out = [""] * k
    for n, (a, b) in enumerate(runs):
        for i in range(a, b + 1):
            out[i] += "abcdefgh"[n]
    return out


def inference(d, ms, dfe):
    header("3. WHICH LEVELS DIFFER (each comparison uses the error term of its own stratum)")
    y = "yield"
    cols = ["comparison", "estimate", "SE", "df", "lower95", "upper95", "p"]
    m = lambda f: d.groupby(f, observed=True)[y].mean()  # noqa: E731
    out = []
    # --- main effects with one df
    g = m("gen")
    out.append(ci_row("Variety: Mindum - Marquis", g["Mindum"] - g["Marquis"], np.sqrt(2 * ms[1] / 80), dfe[1]))
    w = m("dry")
    out.append(ci_row("Method: Dry - Wet", w["Dry"] - w["Wet"], np.sqrt(2 * ms[2] / 80), dfe[2]))
    i = m("inoc")
    out.append(ci_row("Soil: Inoculated - Uninoculated", i["Inoc"] - i["Uninoc"], np.sqrt(2 * ms[3] / 80), dfe[3]))
    main = pd.DataFrame(out, columns=cols)
    print("\n3a. One-df main effects")
    print(main.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))
    # the variety F test has 3 error df, so show what the four blocks support without assuming normality
    vb = d.groupby(["block", "gen"], observed=True)[y].mean().unstack()
    vdiff = (vb["Mindum"] - vb["Marquis"]).values
    flips = np.array(list(itertools.product([-1, 1], repeat=len(vdiff))))
    p_flip = float(np.mean(np.abs((flips * vdiff).mean(1)) >= abs(vdiff.mean()) - 1e-9))
    print("  Variety by block (Mindum - Marquis): " + "   ".join(f"{b}: {v:+.2f}" for b, v in zip(vb.index, vdiff)))
    print(f"  Exact sign-flip test over the 4 blocks: p = {p_flip:.3f} (smallest possible value with 4 blocks: {2 / len(flips):.3f})")

    # --- dust main effect: Tukey HSD on the subplot error
    dm = m("dust").sort_values(ascending=False)
    se_mean = np.sqrt(ms[2] / 32)
    qcrit = stats.studentized_range.ppf(0.95, 5, dfe[2])
    rows, sig = [], [[False] * 5 for _ in range(5)]
    for (a, ma), (b, mb) in itertools.combinations(dm.items(), 2):
        q = abs(ma - mb) / se_mean
        p = float(stats.studentized_range.sf(q, 5, dfe[2]))
        h = qcrit * se_mean
        rows.append([f"{a} - {b}", ma - mb, se_mean * np.sqrt(2), dfe[2], ma - mb - h, ma - mb + h, p])
        ia, ib = list(dm.index).index(a), list(dm.index).index(b)
        sig[ia][ib] = sig[ib][ia] = p < 0.05
    tuk = pd.DataFrame(rows, columns=cols)
    lets = letters(dm.values, sig)
    print("\n3b. Seed dust main effect: Tukey HSD (subplot error, 54 df); means sharing a letter do not differ")
    print("  " + "   ".join(f"{k}: {v:.2f} {l}" for (k, v), l in zip(dm.items(), lets)))
    print(tuk.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))
    chk = d[d.dust == "Check"][y].mean() - d[d.dust != "Check"][y].mean()
    c = ci_row("Check - mean of the four dusts", chk, np.sqrt(ms[2] * (1 / 32 + 1 / 128)), dfe[2])
    print(f"  contrast  {c[0]}: {c[1]:.2f} (95% CI {c[4]:.2f} to {c[5]:.2f}, p = {fmt_p(c[6])})")

    # --- dust x inoculation: effect of inoculation within each dust (sub-subplot error), Bonferroni x5
    cell = d.groupby(["dust", "inoc"], observed=True)[y].mean().unstack()
    rows = [ci_row(f"{k}: Inoc - Uninoc", cell.loc[k, "Inoc"] - cell.loc[k, "Uninoc"], np.sqrt(2 * ms[3] / 16), dfe[3], 5)
            for k in cell.index]
    di = pd.DataFrame(rows, columns=cols)
    print("\n3c. Dust x inoculation: inoculation effect within each dust (sub-subplot error, Bonferroni for 5)")
    print(di.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))
    # dusts compared within one soil condition need both errors: Var = (MS2 + MS3)/16, Satterthwaite df
    v = (ms[2] + ms[3]) / 16
    dfs = (ms[2] + ms[3]) ** 2 / (ms[2] ** 2 / dfe[2] + ms[3] ** 2 / dfe[3])
    rows = [ci_row(f"{s}: Check - {k}", cell.loc["Check", s] - cell.loc[k, s], np.sqrt(v), dfs, 8)
            for s in ["Inoc", "Uninoc"] for k in cell.index if k != "Check"]
    dc = pd.DataFrame(rows, columns=cols)
    print(f"\n3d. Check vs each dust within a soil condition (SE from MS2 and MS3, Satterthwaite df = {dfs:.1f}, Bonferroni for 8)")
    print(dc.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))

    # --- method x inoculation, variety x method
    cell2 = d.groupby(["dry", "inoc"], observed=True)[y].mean().unstack()
    mi = pd.DataFrame([ci_row(f"{k}: Inoc - Uninoc", cell2.loc[k, "Inoc"] - cell2.loc[k, "Uninoc"], np.sqrt(2 * ms[3] / 40), dfe[3], 2)
                       for k in cell2.index], columns=cols)
    print("\n3e. Method x inoculation: inoculation effect within each method (sub-subplot error, Bonferroni for 2)")
    print(mi.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))
    cell3 = d.groupby(["gen", "dry"], observed=True)[y].mean().unstack()
    gm = pd.DataFrame([ci_row(f"{k}: Dry - Wet", cell3.loc[k, "Dry"] - cell3.loc[k, "Wet"], np.sqrt(2 * ms[2] / 40), dfe[2], 2)
                       for k in cell3.index], columns=cols)
    print("\n3f. Variety x method: dry-vs-wet effect within each variety (subplot error, Bonferroni for 2)")
    print(gm.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))

    # --- built-in negative control: "Dry" and "Wet" check plots received no dust at all
    cc = d[d.dust == "Check"].groupby("dry", observed=True)[y].mean()
    nc = ci_row("Check plots: Dry - Wet", cc["Dry"] - cc["Wet"], np.sqrt(2 * ms[2] / 16), dfe[2])
    print("\n3g. Method effect on untreated (Check) seed, where no dust was applied")
    print(pd.DataFrame([nc], columns=cols).to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))

    # --- method effects among the plots that received a dust (error terms pooled from the full model)
    ds = d[d.dust != "Check"]
    n_cell = len(ds) / 4
    w2 = ds.groupby("dry", observed=True)[y].mean()
    c1 = ds.groupby(["gen", "dry"], observed=True)[y].mean().unstack()
    c2 = ds.groupby(["dry", "inoc"], observed=True)[y].mean().unstack()
    gm_int = (c1.loc["Mindum", "Dry"] - c1.loc["Mindum", "Wet"]) - (c1.loc["Marquis", "Dry"] - c1.loc["Marquis", "Wet"])
    mi_int = (c2.loc["Wet", "Inoc"] - c2.loc["Wet", "Uninoc"]) - (c2.loc["Dry", "Inoc"] - c2.loc["Dry", "Uninoc"])
    dusted = pd.DataFrame([
        ci_row("Method: Dry - Wet", w2["Dry"] - w2["Wet"], np.sqrt(ms[2] / n_cell), dfe[2]),
        ci_row("Variety x method: Dry-Wet gap, Mindum - Marquis", gm_int, np.sqrt(4 * ms[2] / n_cell), dfe[2]),
        ci_row("Method x inoculation: Inoc-Uninoc gap, Wet - Dry", mi_int, np.sqrt(4 * ms[3] / n_cell), dfe[3]),
    ], columns=cols)
    print("\n3h. Method effects among the four dusts only (Check plots excluded, because they received no dust to apply)")
    print(dusted.to_string(index=False, float_format="{:.3f}".format, formatters={"p": fmt_p}))

    allc = pd.concat([main.assign(family="main effects"), tuk.assign(family="dust Tukey HSD"),
                      di.assign(family="inoculation within dust"), dc.assign(family="check vs dust within soil"),
                      mi.assign(family="inoculation within method"), gm.assign(family="method within variety"),
                      pd.DataFrame([nc], columns=cols).assign(family="negative control"),
                      dusted.assign(family="method effects, dusted plots only")])
    allc.to_csv(OUT / "comparisons.csv", index=False)
    return dict(p_flip=p_flip, dust_means=dm, dust_letters=lets, se_dust=se_mean, dust_inoc=cell, method_inoc=cell2, gen_method=cell3)


# ======================================================================================= 6. figures
def fig_layout(d):
    fig, ax = plt.subplots(figsize=(11, 4.6))
    tint = {"Inoc": "#b7d3f6", "Uninoc": "#f0efec"}
    for _, r in d.iterrows():
        ax.add_patch(Rectangle((r.col - 1, r.row - 1), 1, 1, facecolor=tint[r.inoc], edgecolor=SURFACE, linewidth=0.8))
        ax.text(r.col - 0.5, r.row - 0.2, "I" if r.inoc == "Inoc" else "U", ha="center", va="center", fontsize=7, color=INK2)
    for (_, g), sub in d.groupby(["sp", "wp"]):
        c0, r0 = sub.col.min() - 1, sub.row.iloc[0] - 1
        ax.add_patch(Rectangle((c0, r0), 2, 1, fill=False, edgecolor=INK2, linewidth=1.1))
        s = sub.iloc[0]
        ax.text(c0 + 1, r0 + 0.4, f"{s.dust}\n{s.dry}", ha="center", va="center", fontsize=6.5, color=INK, linespacing=0.95)
    for blk, sub in d.groupby("block"):
        r0, r1 = sub.row.min() - 1, sub.row.max()
        ax.add_patch(Rectangle((0, r0), 20, r1 - r0, fill=False, edgecolor=INK, linewidth=2.2))
        ax.text(20.25, (r0 + r1) / 2, f"Block {blk[1:]}", va="center", fontsize=9, color=INK)
    for wp, sub in d.groupby("wp"):
        ax.text(-0.25, sub.row.iloc[0] - 0.5, sub.gen.iloc[0], ha="right", va="center", fontsize=9, color=INK)
    ax.set_xlim(-2.6, 22.4)
    ax.set_ylim(8.05, -0.05)
    ax.axis("off")
    ax.set_title("Field layout: three sizes of experimental unit", loc="left")
    ax.text(0, 8.45, "Row = whole plot (variety)   |   outlined pair = subplot (dust x method)   |   "
                     "single cell = sub-subplot (I = inoculated soil, U = uninoculated)", fontsize=8.5, color=INK2)
    fig.savefig(OUT / "fig1_field_layout.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _two_series(ax, x, ya, yb, la, lb, half, xlabels, ylabel):
    xs = np.arange(len(x))
    for yv, col, lab, off, mk in ((ya, BLUE, la, -0.04, "o"), (yb, ORANGE, lb, 0.04, "s")):
        ax.errorbar(xs + off, yv, yerr=half, color=col, marker=mk, capsize=3, elinewidth=1.2, label=lab,
                    markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax.annotate(lab, (xs[-1] + off, yv[-1]), xytext=(9, 0), textcoords="offset points", va="center", fontsize=9, color=INK2)
    ax.set_xticks(xs, xlabels)
    ax.set_xlim(-0.35, len(x) - 1 + 0.95)
    ax.set_ylabel(ylabel)
    ax.grid(axis="x", visible=False)


def fig_effects(d, inf, ms, dfe):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), gridspec_kw={"width_ratios": [2.3, 1, 1]})
    t3, t2 = stats.t.ppf(0.975, dfe[3]), stats.t.ppf(0.975, dfe[2])
    c = inf["dust_inoc"].loc[inf["dust_means"].index]
    _two_series(axes[0], c.index, c["Inoc"].values, c["Uninoc"].values, "Inoculated", "Uninoculated",
                t3 * np.sqrt(ms[3] / 16), list(c.index), "Mean yield")
    axes[0].set_title("Seed dust x soil inoculation (p = 0.0003)")
    c = inf["method_inoc"]
    _two_series(axes[1], c.index, c["Inoc"].values, c["Uninoc"].values, "Inoculated", "Uninoculated",
                t3 * np.sqrt(ms[3] / 40), list(c.index), "")
    axes[1].set_title("Method x inoculation (p = 0.039)")
    c = inf["gen_method"].T
    _two_series(axes[2], c.index, c["Mindum"].values, c["Marquis"].values, "Mindum", "Marquis",
                t2 * np.sqrt(ms[2] / 40), list(c.index), "")
    axes[2].set_title("Variety x method (p = 0.035)")
    # a variety difference within one method uses whole-plot and subplot error: Var = (MS1 + MS2)/40, Satterthwaite df
    df_var = (ms[1] + ms[2]) ** 2 / (ms[1] ** 2 / dfe[1] + ms[2] ** 2 / dfe[2])
    margin_var = stats.t.ppf(0.975, df_var) * np.sqrt((ms[1] + ms[2]) / 40)
    axes[2].text(0.5, 0.42, f"Variety gap at one method:\n95% margin ±{margin_var:.1f}", transform=axes[2].transAxes,
                 ha="center", va="center", fontsize=8.5, color=INK2)
    for ax in axes:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)
    fig.text(0.01, -0.09, "Bars: 95% margins from the error term for one comparison only. Panels 1-2: inoculated vs. uninoculated "
             "(sub-subplot error).\nPanel 3: dry vs. wet within a variety (subplot error). "
             "Panel 3 bars do not apply to comparing the two varieties.", fontsize=8.5, color=INK2)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_interactions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    dm, lets = inf["dust_means"], inf["dust_letters"]
    half = stats.studentized_range.ppf(0.95, 5, dfe[2]) * inf["se_dust"] / 2  # non-overlap <=> Tukey significant
    ys = np.arange(len(dm))[::-1]
    ax.errorbar(dm.values, ys, xerr=half, fmt="o", color=BLUE, capsize=3, elinewidth=1.4, markeredgecolor=SURFACE, markeredgewidth=1.2)
    for yv, mv, l in zip(ys, dm.values, lets):
        ax.text(mv + half + 0.25, yv, f"{mv:.1f}  {l}", va="center", fontsize=9, color=INK2)
    ax.set_yticks(ys, dm.index)
    ax.set_xlabel("Mean yield (averaged over variety, method and soil)")
    ax.set_title("Seed dust main effect: Tukey HSD")
    ax.grid(axis="y", visible=False)
    ax.set_xlim(dm.min() - half - 1, dm.max() + half + 3)
    fig.text(0.01, -0.04, "Bars are half the Tukey HSD width: two dusts differ at the 5% level when their bars do not overlap.\n"
             "Dusts sharing a letter are not significantly different.", fontsize=8.5, color=INK2)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_dust_means.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_diagnostics(d, res, rng):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6))
    names = {2: "Subplot error (54 df)", 3: "Sub-subplot error (60 df)"}
    for row, k in enumerate((2, 3)):
        r, M, df = res[k]["r"], res[k]["M"], res[k]["df"]
        s = np.sqrt(r @ r / df)
        lev = np.sqrt(np.clip(np.diag(M), 1e-12, None))
        t = r / lev / s
        n = len(t)
        q = stats.norm.ppf((np.arange(1, n + 1) - 0.375) / (n + 0.25))
        R = rng.standard_normal((1000, n)) @ M
        T = np.sort(R / lev / np.sqrt((R**2).sum(1, keepdims=True) / df), axis=1)
        lo, hi = np.percentile(T, [2.5, 97.5], axis=0)
        ax = axes[row, 0]
        ax.fill_between(q, lo, hi, color="#cde2fb", linewidth=0, label="95% simulation envelope")
        ax.plot(q, q, color=AXIS, linewidth=1)
        ax.scatter(q, np.sort(t), s=16, color=BLUE, edgecolor=SURFACE, linewidth=0.5, zorder=3)
        ax.set_title(f"{names[k]}: normal Q-Q")
        ax.set_xlabel("Normal quantile")
        ax.set_ylabel("Studentized residual")
        ax.legend(loc="upper left", fontsize=8)
        ax = axes[row, 1]
        ax.axhline(0, color=AXIS, linewidth=1)
        ax.scatter(res[k]["fitted"], t, s=16, color=BLUE, edgecolor=SURFACE, linewidth=0.5)
        ax.set_title("Residual vs fitted")
        ax.set_xlabel("Fitted value")
        ax = axes[row, 2]
        fr = res[k]["frame"]
        grp = fr["dust"].astype(str).values
        order = ["Check", "CaCO3", "Ceresan", "DuBay", "Semesan"]
        jit = rng.uniform(-0.14, 0.14, n)
        ax.axhline(0, color=AXIS, linewidth=1)
        ax.scatter([order.index(g) for g in grp] + jit, t, s=16, color=BLUE, edgecolor=SURFACE, linewidth=0.5)
        ax.set_xticks(range(5), order)
        ax.set_title("Residual spread by seed dust")
        ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_residual_diagnostics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # field map of sub-subplot residuals (diverging: two hues around a neutral midpoint)
    cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f3a583", "#b8451a"])
    grid = np.full((8, 20), np.nan)
    grid[d["row"].values - 1, d["col"].values - 1] = res[3]["r"]
    fig, ax = plt.subplots(figsize=(10, 3.6))
    v = np.nanmax(np.abs(grid))
    im = ax.imshow(grid, cmap=cmap, vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(np.arange(-0.5, 20, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.grid(which="major", visible=False)
    ax.tick_params(which="both", length=0)
    ax.set_xticks(range(0, 20, 2), [str(i) for i in range(1, 21, 2)])
    ax.set_yticks(range(8), [str(i) for i in range(1, 9)])
    ax.set_xlabel("Field column")
    ax.set_ylabel("Field row (whole plot)")
    ax.set_title("Sub-subplot residuals in field position (no spatial pattern expected)")
    fig.colorbar(im, ax=ax, label="Residual", fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_residual_field_map.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ======================================================================================= 6. comparison method
def unequal_variance():
    header("6. COMPARISON METHOD: REML MIXED MODELS WITH UNEQUAL RESIDUAL VARIANCES (nlme, unequal_variance.R)\n"
           "Follows up the closest assumption check in section 4: residual variance by block, and by seed dust.")
    rscript = shutil.which("Rscript")
    ran = False
    if rscript:
        run = subprocess.run([rscript, "unequal_variance.R"], cwd=HERE, capture_output=True, text=True)
        ran = run.returncode == 0
        if not ran:
            print("Rscript failed, so the saved results are shown instead:\n" + run.stderr.strip()[-400:])
    txt, csv = OUT / "unequal_variance.txt", OUT / "unequal_variance.csv"
    if not (txt.exists() and csv.exists()):
        print("Rscript is not installed and outputs/unequal_variance.* are missing: sections 6 and 7 skipped.")
        return None
    print("source: " + ("unequal_variance.R, run just now" if ran else
                        "saved outputs/unequal_variance.txt (Rscript not available; run Rscript unequal_variance.R to refresh)"))
    print()
    print(txt.read_text(encoding="utf-8").rstrip())
    return pd.read_csv(csv)


# ======================================================================================= 7. robustness
def robustness(a, uv, p_flip):
    header("7. HOW ROBUST IS EACH SIGNIFICANT EFFECT\n"
           "Every term with p < 0.05 in section 2, re-tested with the Holm adjustment and under both unequal-variance models.")
    t = a.dropna(subset=["p"]).merge(uv, on="term")
    sig = t[t["p"] < 0.05].copy()

    def reading(r):
        if r["p_holm"] < 0.05 and r["p_largest"] < 0.05:
            return "solid"
        if r["p_largest"] < 0.05:
            return "holds under every variance model, not after Holm"
        return "suggestive only"

    sig["reading"] = sig.apply(reading, axis=1)
    rob = sig[["term", "p", "p_holm", "p_variance_by_block", "p_variance_by_dust", "reading"]].rename(
        columns={"p": "p_anova"})
    order = {"solid": 0, "holds under every variance model, not after Holm": 1, "suggestive only": 2}
    rob = rob.sort_values("reading", key=lambda c: c.map(order), kind="stable")
    pc = ["p_anova", "p_holm", "p_variance_by_block", "p_variance_by_dust"]
    print(rob.to_string(index=False, formatters={c: fmt_p for c in pc}))
    rob.to_csv(OUT / "robustness.csv", index=False)
    counts = rob["reading"].value_counts()
    print("\nsummary: " + "; ".join(f"{int(counts.get(k, 0))} {k}" for k in order))
    print(f"variety: only 3 error df; the exact sign-flip test over the 4 blocks gives p = {p_flip:.3f}, "
          "so its F test leans on normality")
    other = t[(t["p"] >= 0.05) & ((t["p_variance_by_block"] < 0.05) | (t["p_variance_by_dust"] < 0.05))]
    if len(other):
        print("\nterms not significant in the ANOVA that fall below 0.05 under one unequal-variance model:")
        print(other[["term", "p", "p_variance_by_block", "p_variance_by_dust"]].rename(columns={"p": "p_anova"})
              .to_string(index=False, formatters={c: fmt_p for c in pc}))
        print("These are not treated as findings: none holds under both variance models, the by-block variance "
              "ratios are\nimprecise (section 6), and these F tests treat estimated variances as known.")
    return rob


# ======================================================================================= main
def main():
    sys.stdout = Tee(OUT / "results.txt")
    d = load()
    describe(d)
    a, ms, dfe = anova(d)
    inf = inference(d, ms, dfe)
    res, tests, rng = assumptions(d, a)
    mixed_model(d, a, ms)
    uv = unequal_variance()
    if uv is not None:
        robustness(a, uv, inf["p_flip"])
    fig_layout(d)
    fig_effects(d, inf, ms, dfe)
    fig_diagnostics(d, res, rng)
    header("FILES WRITTEN")
    for p in sorted(OUT.iterdir()):
        print("  outputs/" + p.name)


if __name__ == "__main__":
    main()
