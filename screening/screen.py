"""Uniform assumption screening for candidate STAT 406 Part 2 datasets.

For every candidate the same checks are run on the residuals of the model the
design calls for (treatment factorial + blocks + one error term per size of
experimental unit):

  normality       Shapiro-Wilk W on the residuals of each error stratum
  equal variance  Brown-Forsythe (median Levene) F over every treatment factor,
                  block, and every pair of treatment factors; one family-wise
                  p-value (min-p) per stratum
  mean-variance   Spearman correlation of |residual| with fitted value
  outliers        largest absolute studentized residual
  additivity      Tukey one-degree-of-freedom test (fitted^2 term)
  independence    Moran-type lag-1 statistic over neighbouring plots when the
                  field layout (row, col) is known; otherwise by design only
  scale           Box-Cox profile likelihood interval for lambda

Residuals from a fitted model are neither independent nor identically
distributed, so textbook p-values are off.  Every p-value except Tukey's is
therefore calibrated by simulation: normal errors are pushed through the same
residual-forming matrix (2,000 draws) and the observed statistic is ranked
against them.  All statistics are scale free, so this is exact under the model.

Usage:  python screening/screen.py            (all candidates)
        python screening/screen.py goulden.splitsplit
"""
import itertools
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import patsy

HERE = pathlib.Path(__file__).parent
SEED = 406


sys.path.insert(0, str(HERE.parent))
from fetch import fetch  # noqa: E402
from doetools import (boxcox_interval, dummies, make_groupings, neighbour_pairs, residual_maker,  # noqa: E402
                      stratum_tests, tukey_nonadditivity, unit_ids)


def screen(spec, transform="identity", seed=SEED):
    rng = np.random.default_rng(seed)
    if not (HERE / "raw" / f"{spec['item']}.csv").exists():
        fetch(spec["item"])  # downloads the CSV and its documentation page from Rdatasets
    df = pd.read_csv(HERE / "raw" / f"{spec['item']}.csv").drop(columns="rownames", errors="ignore")
    if "prepare" in spec:
        df = spec["prepare"](df)
    df = df.dropna(subset=[spec["response"]]).reset_index(drop=True)
    for f in spec["factors"] + spec.get("block_factors", []):
        df[f] = df[f].astype(str)
    y = df[spec["response"]].astype(float).values
    y = {"identity": y, "sqrt": np.sqrt(y), "log": np.log(y) if y.min() > 0 else np.log1p(y)}[transform]
    Xf = np.asarray(patsy.dmatrix(spec["fixed"], df))
    chain = spec.get("strata", [])  # nested, coarse -> fine
    extra = spec.get("crossed_strata", [])
    ids = [unit_ids(df, s) for s in chain]
    D = [dummies(i) for i in ids] + [dummies(unit_ids(df, s)) for s in extra]
    X = np.hstack([Xf] + D)
    M, df_res = residual_maker(X)
    r = M @ y
    res = {
        "item": spec["item"], "transform": transform, "n": int(len(y)),
        "n_predictors": len(spec["factors"]) + len(spec.get("covariates", [])),
        "strata": {},
    }
    # lowest stratum
    pair_factors = [a + ":" + b for a, b in itertools.combinations(spec["factors"], 2)]
    for a, b in itertools.combinations(spec["factors"], 2):
        df[a + ":" + b] = df[a] + "|" + df[b]
    groupings = make_groupings(df, spec["factors"] + spec.get("block_factors", []) + pair_factors)
    pairs = neighbour_pairs(df) if spec.get("spatial") and {"row", "col"} <= set(df.columns) else None
    low = stratum_tests(r, M, df_res, y - r, groupings, pairs, rng)
    low["p_nonadditive"] = tukey_nonadditivity(y, X)
    res["strata"]["residual"] = low
    res["boxcox"] = boxcox_interval(df[spec["response"]].astype(float).values, X) if transform == "identity" else None
    res["zeros"] = int((df[spec["response"]] == 0).sum())
    # upper strata of the nested chain: regress unit means on fixed effects + parent units
    for k in range(len(chain) - 1, 0, -1):
        u = pd.Series(ids[k])
        ybar = pd.Series(y).groupby(u.values).mean()
        Xbar = pd.DataFrame(Xf).groupby(u.values).mean().loc[ybar.index].values
        parent = pd.Series(ids[k - 1]).groupby(u.values).first().loc[ybar.index].values
        Xk = np.hstack([Xbar, dummies(parent)])
        Mk, dfk = residual_maker(Xk)
        if dfk < 6:
            res["strata"][":".join(chain[k])] = {"n_units": int(len(ybar)), "df": int(dfk), "note": "too few df to test"}
            continue
        first = df.groupby(u.values).first().loc[ybar.index]
        const = [f for f in spec["factors"] + spec.get("block_factors", []) + pair_factors
                 if (df.groupby(u.values)[f].nunique() == 1).all()]
        gk = make_groupings(first, const)
        rk = Mk @ ybar.values
        res["strata"][":".join(chain[k])] = stratum_tests(rk, Mk, dfk, ybar.values - rk, gk, None, rng)
    return res


def verdict(res):
    ps = []
    for name, s in res["strata"].items():
        for key in ("p_normal", "p_var_family", "p_trend", "p_outlier", "p_spatial", "p_nonadditive"):
            if s.get(key) is not None:
                ps.append((s[key], f"{name}:{key[2:]}"))
    fails = [n for p, n in ps if p < 0.01]
    warns = [n for p, n in ps if 0.01 <= p < 0.05]
    grade = "FAIL" if fails else ("PASS" if len(warns) <= 1 else "WARN")
    return grade, fails, warns


if __name__ == "__main__":
    from candidates import CANDIDATES

    wanted = sys.argv[1:] or list(CANDIDATES)
    results = {}
    for item in wanted:
        spec = CANDIDATES[item]
        for tr in ("identity", "sqrt", "log"):
            res = screen(spec, tr)
            grade, fails, warns = verdict(res)
            res.update(grade=grade, fails=fails, warns=warns)
            results[f"{item}|{tr}"] = res
            low = res["strata"]["residual"]
            print(f"{item:28s} {tr:8s} {grade:4s} n={res['n']:3d} df={low['df']:3d} "
                  f"norm={low['p_normal']:.3f} var={low.get('p_var_family', float('nan')):.3f} "
                  f"(worst {low.get('var_worst_grouping', '-')}, ratio {low.get('var_worst_ratio', float('nan')):.1f}) "
                  f"trend={low.get('p_trend', float('nan')):.3f} outl={low['p_outlier']:.3f} "
                  f"nonadd={low['p_nonadditive'] if low['p_nonadditive'] is None else round(low['p_nonadditive'], 3)} "
                  f"spatial={low.get('p_spatial', 'n/a')}")
            for name, s in res["strata"].items():
                if name != "residual":
                    print(f"{'':28s} stratum {name}: df={s['df']} " + (s.get("note") or
                          f"norm={s['p_normal']:.3f} var={s.get('p_var_family', float('nan')):.3f} outl={s['p_outlier']:.3f} trend={s.get('p_trend', float('nan')):.3f}"))
            if fails or warns:
                print(f"{'':28s} fails={fails} warns={warns}")
            if tr == "identity" and res["boxcox"]:
                print(f"{'':28s} Box-Cox lambda={res['boxcox']['lambda_hat']} CI={res['boxcox']['ci']}")
    out = HERE / "results.json"
    old = json.loads(out.read_text()) if out.exists() and sys.argv[1:] else {}
    old.update(results)
    out.write_text(json.dumps(old, indent=1))
