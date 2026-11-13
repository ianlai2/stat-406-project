"""Turn screening/results.json into the table in SCREENING.md.  Run after screen.py."""
import json
import pathlib

import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from candidates import CANDIDATES  # noqa: E402

HERE = pathlib.Path(__file__).parent
res = json.loads((HERE / "results.json").read_text())

NOTES = {
    "chakravertti.factorial": "Baseline (the dataset rejected earlier): 19 of 27 yields are zero in one cell; no transformation rescues it.",
    "goulden.splitsplit": "Selected. Three sizes of experimental unit, 4 treatment factors, balanced, published reference ANOVA.",
    "gregory.cotton": "Passes only after a log transform; values are treatment means and there is a single error term.",
    "gomez.fractionalfactorial": "Runner-up. Cleanest residuals of all, 6 factors, but one size of experimental unit and 2-level factors only.",
    "chinloy.fractionalfactorial": "Unreplicated; residuals of neighbouring plots are correlated.",
    "cochran.factorial": "Passes, but only 32 plots and 14 error df.",
    "besag.triticale": "Neighbouring-plot residuals strongly correlated (the reason the dataset is famous); nitrogen has an unbalanced 18/36 split in the file.",
    "harvey.lsmeans": "Observational and unbalanced (65 steers); mean-variance trend flagged.",
    "gomez.splitsplit": "Passes on the raw scale, but only 3 treatment factors.",
    "gomez.stripsplitplot": "Passes, but only 3 treatment factors.",
    "cox.stripsplit": "Passes, but only 3 treatment factors.",
    "mead.cowpea.maize": "Passes, but only 3 treatment factors and one error term.",
    "woodman.pig": "Passes, but 30 pigs and only 2 treatment factors plus a covariate.",
    "cochran.eelworms": "Flags on variance, additivity and neighbouring plots; 2 treatment factors plus a covariate.",
}


def cell(r):
    low = r["strata"]["residual"]
    f = lambda k: "n/a" if low.get(k) is None else f"{low[k]:.3f}"  # noqa: E731
    return f"{f('p_normal')} | {f('p_var_family')} | {f('p_trend')} | {f('p_outlier')} | {f('p_nonadditive')} | {f('p_spatial')}"


lines = [
    "| Dataset | Design | n | Treatment predictors | Raw-scale result | Normality | Equal variance | Variance vs mean | Outliers | Additivity | Neighbouring plots | Other scales | Note |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
]
for item, spec in CANDIDATES.items():
    r = res[f"{item}|identity"]
    upper = [f"{k} ({s['df']} df): " + (s.get("note") or f"normality {s['p_normal']:.2f}, variance {s.get('p_var_family', float('nan')):.2f}")
             for k, s in r["strata"].items() if k != "residual"]
    other = ", ".join(f"{t} {res[f'{item}|{t}']['grade']}" for t in ("sqrt", "log"))
    npred = f"{len(spec['factors'])}" + (f" + {len(spec['covariates'])} covariate(s)" if spec.get("covariates") else "")
    lines.append(f"| `{item}` | {spec['design']} | {r['n']} | {npred} | **{r['grade']}** | {cell(r)} | {other} | {NOTES[item]}"
                 + (" Upper strata: " + "; ".join(upper) + "." if upper else "") + " |")
table = "\n".join(lines)
(HERE / "screening_table.md").write_text(table + "\n", encoding="utf-8")
print(table)
