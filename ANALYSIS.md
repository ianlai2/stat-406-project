# Full analysis of the model

Every number here comes from `python analysis.py`; the section numbers in brackets point to
[outputs/results.txt](outputs/results.txt). The [README](README.md) has the dataset, the design and the figures.

**Summary.** The model is a four-factor split-split-plot with three error terms. It supports three firm
conclusions: seed dust matters, soil inoculation matters, and the two interact. Variety and application method
hold with caveats, and two interactions involving method are suggestive only.

## The model

$$
\begin{aligned}
y_{ijklm} = \mu &+ \rho_i + \alpha_j + \eta_{ij} \\
&+ \beta_k + \gamma_l + (\alpha\beta)_{jk} + (\alpha\gamma)_{jl} + (\beta\gamma)_{kl} + (\alpha\beta\gamma)_{jkl} + \varepsilon_{ijkl} \\
&+ \delta_m + (\text{all interactions with } \delta) + e_{ijklm}
\end{aligned}
$$

- **Fixed effects:** variety $\alpha$ (2 levels), dust $\beta$ (5), method $\gamma$ (2), inoculation $\delta$ (2), and all their interactions.
- **Random effects:** block $\rho_i$, whole-plot error $\eta_{ij} \sim N(0,\sigma_1^2)$, subplot error
  $\varepsilon_{ijkl} \sim N(0,\sigma_2^2)$, sub-subplot error $e_{ijklm} \sim N(0,\sigma_3^2)$, all independent.

There are three errors because randomization happened at three sizes of unit. Variety was assigned to 8 half-blocks,
the dust and method combination to 80 subplots, and inoculation to 160 sub-subplots. Each effect is tested against
the variation among the units it was randomized to.

## Expected mean squares and variance components

$$
\begin{aligned}
E[MS_{E3}] &= \sigma_3^2 \\
E[MS_{E2}] &= \sigma_3^2 + 2\sigma_2^2 \\
E[MS_{E1}] &= \sigma_3^2 + 2\sigma_2^2 + 20\sigma_1^2 \\
E[MS_{\text{Block}}] &= \sigma_3^2 + 2\sigma_2^2 + 20\sigma_1^2 + 40\sigma_b^2
\end{aligned}
$$

Solving from the bottom up gives the variance components (section 5):

$$
\begin{aligned}
\hat\sigma_3^2 &= 15.39 \\
\hat\sigma_2^2 &= \frac{19.62 - 15.39}{2} = 2.11 \\
\hat\sigma_1^2 &= \frac{215.87 - 19.62}{20} = 9.81 \\
\hat\sigma_b^2 &= \frac{329.84 - 215.87}{40} = 2.85
\end{aligned}
$$

Most of the noise is at the smallest unit (15.39) and at the whole-plot level (9.81). Subplots add little beyond that.

## The tests (section 2)

| Stratum | Error MS (df) | Significant terms | F | p |
|---|---|---|---|---|
| Whole plot | 215.87 (3) | Variety | 16.86 | 0.026 |
| Subplot | 19.62 (54) | Dust | 12.59 | <0.0001 |
| | | Method | 5.98 | 0.018 |
| | | Variety x method | 4.66 | 0.035 |
| Sub-subplot | 15.39 (60) | Inoculation | 62.72 | <0.0001 |
| | | Dust x inoculation | 6.17 | 0.0003 |
| | | Method x inoculation | 4.48 | 0.039 |

The other eight interactions are not significant (p from 0.12 to 0.97).

**Worked example, variety.** The test and its confidence interval both use the whole-plot error:

$$
F = \frac{MS_{\text{Variety}}}{MS_{E1}} = \frac{3638.56}{215.87} = 16.86 \quad \text{on } (1, 3) \text{ df}
$$

$$
SE = \sqrt{\frac{2 \cdot MS_{E1}}{80}} = \sqrt{\frac{2(215.87)}{80}} = 2.32
$$

$$
9.54 \pm t_{0.975,\,3}\,(2.32) = 9.54 \pm 3.182(2.32) = (2.1,\ 16.9)
$$

The same 9.5-unit difference tested against the sub-subplot error would look overwhelmingly significant, and that
would be wrong. Variety has the largest effect in the experiment and the weakest test, because there are only 8 whole plots.

## What the effects are (section 3)

- **Dust:** untreated Check seed yields most (71.2), then CaCO3 (68.8), Semesan (65.8), Ceresan (65.1), DuBay (64.8).
  By Tukey HSD, Check beats the last three and CaCO3 beats the last two.
- **Inoculation:** inoculated soil yields 4.9 more on average (95% CI 3.7 to 6.2).
- **Dust x inoculation:** the main story. The inoculation gain is 10.3 for Check and 6.3 for CaCO3, but only 2.3 to 3.4
  for the three mercury dusts. So the dust main effect is mostly a statement about inoculated soil: Check beats every
  dust there and none of them on uninoculated soil.
- **Variety:** Mindum beats Marquis by 9.5, and leads in all four blocks.
- **Method:** dry application beats wet by 1.7. Among plots that actually received a dust, it is 1.5 (p = 0.061).
- **Method interactions:** dry beats wet by 3.2 for Mindum and 0.2 for Marquis; the inoculation gain is 6.2 with wet
  application and 3.6 with dry.

## Do the assumptions hold (section 4)

- **Normality, outliers, additivity, spatial independence:** all pass in both testable strata (smallest p = 0.18).
- **Equal variance:** the weak point. Across blocks in the sub-subplot stratum the max/min variance ratio is 3.9
  (p = 0.011 alone, 0.066 after allowing for 11 groupings).
- **Whole-plot stratum:** with 3 df it cannot be checked at all. An exact sign-flip test on the four block differences
  gives p = 0.125, the smallest value possible, so variety's 0.026 depends on normality.
- **Box-Cox:** $\hat\lambda = 1.85$, interval 1.05 to 2.65, just excluding 1. Squaring the yield keeps all seven effects significant.

## The comparison method (sections 5 and 6)

The equal-variance REML mixed model reproduces the ANOVA F statistics to within 0.01%, as it must for balanced data,
so it validates the code but is not a real comparison. The real comparison lets the residual variance differ by block
or by dust. Both fit better than one variance (likelihood-ratio p = 0.035 and 0.028).

The by-dust model is well estimated. The by-block model is not: every variance ratio's interval contains 1, so it is
used only as a stress test.

## How far to trust each effect (section 7)

| Effect | ANOVA p | Holm p | Variance by block | Variance by dust | Verdict |
|---|---|---|---|---|---|
| Dust | <0.0001 | <0.0001 | <0.0001 | <0.0001 | Solid |
| Inoculation | <0.0001 | <0.0001 | <0.0001 | <0.0001 | Solid |
| Dust x inoculation | 0.0003 | 0.004 | <0.0001 | 0.0002 | Solid |
| Variety | 0.026 | 0.29 | 0.028 | 0.030 | Holds with caveats |
| Method | 0.018 | 0.21 | 0.034 | 0.026 | Holds with caveats |
| Variety x method | 0.035 | 0.35 | 0.12 | 0.064 | Suggestive |
| Method x inoculation | 0.039 | 0.35 | 0.22 | 0.013 | Suggestive |

## Limits

- **The direction is odd.** In a root-rot study one would expect inoculation to hurt and seed treatment to help; both
  go the other way. The numbers match Goulden's printed totals, so it is not a coding error, but the model cannot explain it.
- **Method means nothing for Check.** Check plots received no dust, so "dry" and "wet" are only labels there. The
  dusted-only contrasts are the cleaner read.
- **Blocks are treated as additive.** The Tukey one-df test supports this (p = 0.30 and 0.31), but with 4 blocks it has little power.
- **The yield units are unknown,** so effect sizes can only be read relative to the mean of 67.1.
