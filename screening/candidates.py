"""Model specifications for every screened dataset (all from the R package agridat, via Rdatasets).

fixed          treatment model (patsy formula)
strata         nested error strata, coarse -> fine; whatever is left is the residual stratum
crossed_strata extra unit sets for strip-plot designs (added as dummies, not tested separately)
factors        treatment factors that count toward the "four predictors" requirement
covariates     numerical predictors
block_factors  blocking variables (do not count as predictors here)
spatial        True when row/col field positions are available for an independence check
"""


def _letters(df, letters):
    """Rebuild 0/1 factor columns from treatment labels such as 'abd' or '(1)'."""
    for L in letters:
        df[L.upper()] = df["trt"].astype(str).apply(lambda t: "1" if L in t else "0")
    return df


def _eelworm(df):
    df["trt"] = df["fumigant"] + df["dose"].astype(str)
    return df


def _pig(df):
    df["gain"] = df["weight2"] - df["weight1"]
    return df


def _triticale(df):
    df["rate"] = df["rate"].astype(str)
    df["nitro"] = df["nitro"].astype(str)
    return df


CANDIDATES = {
    # ---- baseline: the dataset ChatGPT rejected (pipeline should reject it too)
    "chakravertti.factorial": dict(
        item="chakravertti.factorial", response="yield",
        fixed="C(date)*C(gen)*C(spacing)*C(seeds)",
        strata=[["block"], ["block", "date"], ["block", "date", "gen"]],
        factors=["date", "gen", "spacing", "seeds"], block_factors=["block"],
        design="split-split-plot, 3 blocks"),
    # ---- tier 1: four or more treatment predictors
    "goulden.splitsplit": dict(
        item="goulden.splitsplit", response="yield",
        fixed="C(gen)*C(dust)*C(dry)*C(inoc)",
        strata=[["block"], ["block", "gen"], ["block", "gen", "trt"]],
        factors=["gen", "dust", "dry", "inoc"], block_factors=["block"], spatial=True,
        design="split-split-plot, 4 blocks"),
    "gregory.cotton": dict(
        item="gregory.cotton", response="yield",
        fixed="C(nitrogen)*C(date)*C(water)*C(spacing)",
        strata=[["year"]],
        factors=["nitrogen", "date", "water", "spacing"], block_factors=["year"],
        design="4-factor factorial, years as blocks (treatment means)"),
    "gomez.fractionalfactorial": dict(
        item="gomez.fractionalfactorial", response="yield", prepare=lambda d: _letters(d, "abcdef"),
        fixed="(A+B+C+D+E+F)**2",
        strata=[["rep"], ["rep", "block"]],
        factors=list("ABCDEF"), block_factors=["rep"], spatial=True,
        design="1/2 of 2^6 in incomplete blocks, 2 reps"),
    "chinloy.fractionalfactorial": dict(
        item="chinloy.fractionalfactorial", response="yield",
        fixed="(C(n)+C(p)+C(k)+C(b)+C(m))**2",
        strata=[["block"]],
        factors=["n", "p", "k", "b", "m"], block_factors=["block"], spatial=True,
        design="1/3 of 3^5 in 9 blocks, unreplicated"),
    "cochran.factorial": dict(
        item="cochran.factorial", response="yield", prepare=lambda d: _letters(d, "dnpk"),
        fixed="(D+N+P+K)**3",
        strata=[["rep"], ["rep", "block"]],
        factors=list("DNPK"), block_factors=["rep"],
        design="2^4 with DNPK confounded with blocks, 2 reps"),
    "besag.triticale": dict(
        item="besag.triticale", response="yield", prepare=_triticale,
        fixed="(C(gen)+C(rate)+C(nitro)+C(regulator))**2",
        strata=[],
        factors=["gen", "rate", "nitro", "regulator"], spatial=True,
        design="4-factor factorial, completely randomized, unreplicated"),
    "harvey.lsmeans": dict(
        item="harvey.lsmeans", response="adg",
        fixed="C(line)+C(damage)+weanage+weight",
        strata=[["sire"]],
        factors=["line", "damage"], covariates=["weanage", "weight"], block_factors=["sire"],
        design="unbalanced; sires nested in lines (random), 2 covariates"),
    # ---- tier 2: three treatment predictors (a fourth only if blocks or a covariate count)
    "gomez.splitsplit": dict(
        item="gomez.splitsplit", response="yield",
        fixed="C(nitro)*C(management)*C(gen)",
        strata=[["rep"], ["rep", "nitro"], ["rep", "nitro", "management"]],
        factors=["nitro", "management", "gen"], block_factors=["rep"], spatial=True,
        design="split-split-plot, 3 reps"),
    "gomez.stripsplitplot": dict(
        item="gomez.stripsplitplot", response="yield",
        fixed="C(nitro)*C(gen)*C(planting)",
        strata=[["rep"]],
        crossed_strata=[["rep", "nitro"], ["rep", "gen"], ["rep", "nitro", "gen"]],
        factors=["nitro", "gen", "planting"], block_factors=["rep"], spatial=True,
        design="strip-split-plot, 3 reps"),
    "cox.stripsplit": dict(
        item="cox.stripsplit", response="yield",
        fixed="C(soil)*C(fert)*C(calcium)",
        strata=[["rep"]],
        crossed_strata=[["rep", "fert"], ["rep", "soil"], ["rep", "fert", "calcium"], ["rep", "fert", "soil"]],
        factors=["soil", "fert", "calcium"], block_factors=["rep"],
        design="strip-split-plot, 4 reps"),
    "mead.cowpea.maize": dict(
        item="mead.cowpea.maize", response="myield",
        fixed="C(nitro)*C(cowpea)*C(maize)",
        strata=[["block"]],
        factors=["nitro", "cowpea", "maize"], block_factors=["block"],
        design="3-factor factorial in 3 blocks"),
    "woodman.pig": dict(
        item="woodman.pig", response="gain", prepare=_pig,
        fixed="C(treatment)*C(sex)+w0",
        strata=[["pen"]],
        factors=["treatment", "sex"], covariates=["w0"], block_factors=["pen"],
        design="RCB with covariate, n=30"),
    "cochran.eelworms": dict(
        item="cochran.eelworms", response="grain", prepare=_eelworm,
        fixed="C(trt)+np.log(initial)",
        strata=[["block"]],
        factors=["fumigant", "dose"], covariates=["initial"], block_factors=["block"], spatial=True,
        design="RCB with covariate, n=48"),
}
