# Dataset screening

Goal: replace the rice experiment (`chakravertti.factorial`), which fails equal variance and normality on every
scale, with a dataset that meets the Part 2 requirements and whose model assumptions actually hold.

## Method

Candidates are designed experiments from the R package agridat (CSV files from Rdatasets) with a numerical
response, complete data, and at least three treatment predictors. Each gets the model its design calls for:
the full treatment factorial (main effects and two-factor interactions where the design is unreplicated or
fractional), blocks, and one error term per size of experimental unit. The same checks are then run on the
residuals, on the raw scale and after square-root and log transforms:

| Check | Statistic |
|---|---|
| Normality | Shapiro-Wilk W |
| Equal variance | Brown-Forsythe F over every treatment factor, block and pair of treatment factors, combined into one family-wise p-value |
| Variance vs mean | Spearman correlation of absolute residual with fitted value |
| Outliers | Largest absolute studentized residual |
| Additivity | Tukey one-degree-of-freedom test |
| Neighbouring plots | Lag-1 residual correlation over adjacent field plots, when row and column positions are recorded |

All p-values except Tukey's are calibrated by simulation (2,000 normal error vectors pushed through the same
residual-forming matrix), because residuals from a fitted model are not independent. Grades: **FAIL** if any
p < 0.01, **WARN** if two or more are between 0.01 and 0.05, otherwise **PASS**. Upper error strata are tested
the same way when they have at least 6 df.

The pipeline reproduces the earlier finding on the rice data (largest-to-smallest residual variance 8,020 to 1
across planting date x variety), which is the sanity check that it detects a real violation.

## Results (raw scale; p-values are for the lowest error stratum)

| Dataset | Design | n | Treatment predictors | Raw-scale result | Normality | Equal variance | Variance vs mean | Outliers | Additivity | Neighbouring plots | Other scales | Note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `chakravertti.factorial` | split-split-plot, 3 blocks | 405 | 4 | **FAIL** | 0.000 | 0.000 | 0.000 | 0.000 | 0.132 | n/a | sqrt FAIL, log FAIL | Baseline (the dataset rejected earlier): 19 of 27 yields are zero in one cell; no transformation rescues it. Upper strata: block:date:gen (20 df): normality 0.71, variance 0.68; block:date (8 df): normality 0.18, variance 0.46. |
| `goulden.splitsplit` | split-split-plot, 4 blocks | 160 | 4 | **PASS** | 0.485 | 0.068 | 0.955 | 0.165 | 0.298 | 0.792 | sqrt PASS, log WARN | Selected. Three sizes of experimental unit, 4 treatment factors, balanced, published reference ANOVA. Upper strata: block:gen:trt (54 df): normality 0.60, variance 0.24; block:gen (3 df): too few df to test. |
| `gregory.cotton` | 4-factor factorial, years as blocks (treatment means) | 144 | 4 | **FAIL** | 0.628 | 0.009 | 0.441 | 0.725 | 0.000 | n/a | sqrt FAIL, log PASS | Passes only after a log transform; values are treatment means and there is a single error term. |
| `gomez.fractionalfactorial` | 1/2 of 2^6 in incomplete blocks, 2 reps | 64 | 6 | **PASS** | 0.590 | 0.989 | 0.907 | 0.635 | 0.975 | 0.929 | sqrt PASS, log PASS | Runner-up. Cleanest residuals of all, 6 factors, but one size of experimental unit and 2-level factors only. Upper strata: rep:block (2 df): too few df to test. |
| `chinloy.fractionalfactorial` | 1/3 of 3^5 in 9 blocks, unreplicated | 81 | 5 | **WARN** | 0.659 | 0.614 | 0.045 | 0.072 | 0.090 | 0.015 | sqrt WARN, log FAIL | Unreplicated; residuals of neighbouring plots are correlated. |
| `cochran.factorial` | 2^4 with DNPK confounded with blocks, 2 reps | 32 | 4 | **PASS** | 0.329 | 0.349 | 0.345 | 0.631 | 0.452 | n/a | sqrt PASS, log PASS | Passes, but only 32 plots and 14 error df. Upper strata: rep:block (2 df): too few df to test. |
| `besag.triticale` | 4-factor factorial, completely randomized, unreplicated | 54 | 4 | **FAIL** | 0.606 | 0.494 | 0.224 | 0.326 | 0.565 | 0.001 | sqrt FAIL, log FAIL | Neighbouring-plot residuals strongly correlated (the reason the dataset is famous); nitrogen has an unbalanced 18/36 split in the file. |
| `harvey.lsmeans` | unbalanced; sires nested in lines (random), 2 covariates | 65 | 2 + 2 covariate(s) | **PASS** | 0.615 | 0.363 | 0.020 | 0.797 | 0.074 | n/a | sqrt WARN, log WARN | Observational and unbalanced (65 steers); mean-variance trend flagged. |
| `gomez.splitsplit` | split-split-plot, 3 reps | 135 | 3 | **PASS** | 0.922 | 0.085 | 0.063 | 0.477 | 0.693 | 0.863 | sqrt FAIL, log FAIL | Passes on the raw scale, but only 3 treatment factors. Upper strata: rep:nitro:management (20 df): normality 0.21, variance 0.85; rep:nitro (8 df): normality 0.50, variance 0.53. |
| `gomez.stripsplitplot` | strip-split-plot, 3 reps | 108 | 3 | **PASS** | 0.750 | 0.708 | 0.462 | 0.207 | 0.838 | 0.701 | sqrt PASS, log FAIL | Passes, but only 3 treatment factors. |
| `cox.stripsplit` | strip-split-plot, 4 reps | 96 | 3 | **PASS** | 0.198 | 0.133 | 0.477 | 0.211 | 0.812 | n/a | sqrt PASS, log PASS | Passes, but only 3 treatment factors. |
| `mead.cowpea.maize` | 3-factor factorial in 3 blocks | 72 | 3 | **PASS** | 0.335 | 0.076 | 0.906 | 0.256 | 0.462 | n/a | sqrt FAIL, log FAIL | Passes, but only 3 treatment factors and one error term. |
| `woodman.pig` | RCB with covariate, n=30 | 30 | 2 + 1 covariate(s) | **PASS** | 0.139 | 0.310 | 0.359 | 0.072 | 0.709 | n/a | sqrt PASS, log PASS | Passes, but 30 pigs and only 2 treatment factors plus a covariate. |
| `cochran.eelworms` | RCB with covariate, n=48 | 48 | 2 + 1 covariate(s) | **WARN** | 0.342 | 0.048 | 0.066 | 0.459 | 0.018 | 0.054 | sqrt PASS, log PASS | Flags on variance, additivity and neighbouring plots; 2 treatment factors plus a covariate. |

## Decision

`goulden.splitsplit` is the choice.

1. It is the only candidate with four treatment factors **and** more than one size of experimental unit that passes
   on the raw scale. The other multi-stratum designs that pass (one split-split-plot, two strip-split-plots) have three treatment factors,
   so they would meet the four-predictor rule only by counting blocks.
2. Both of its testable error strata pass every check (subplot stratum: normality 0.60, variance 0.24).
   Its whole-plot stratum has 3 df and cannot be tested, which is stated in the README.
3. It is balanced with 160 plots, has seven significant effects including three interactions, and comes with a
   published ANOVA table to validate against.

Runner-up: `gomez.fractionalfactorial` has the cleanest residuals of all and six factors, but a single size of
experimental unit and only two-level factors. `gregory.cotton` is usable after a log transform, but its values
are treatment means with one error term.

Not screened: datasets with a binary, count or ordinal response, datasets with missing plots
(`archbold.apple`, `holshouser.splitstrip`), and variety trials whose only real factor is genotype.

## Reproduce

```bash
python screening/screen.py
```

```bash
python screening/summarize.py
```

`screen.py` downloads any missing CSV into `screening/raw/` first. Full numbers are in `screening/results.json`
and `screening/screen_log.txt`.
