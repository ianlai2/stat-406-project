# Sensitivity analysis: does relaxing the equal-variance assumption change the conclusions?
# Fits the split-split-plot mixed model by REML with (a) one residual variance, (b) a residual variance
# for each block, and (c) a residual variance for each seed dust.
# Run from the repository root:  Rscript unequal_variance.R
# Writes outputs/unequal_variance.txt and outputs/unequal_variance.csv.  Needs only nlme, which ships with R.
library(nlme)

dir.create("outputs", showWarnings = FALSE)
d <- read.csv("data/goulden_splitsplit.csv", stringsAsFactors = TRUE)
d$trt <- factor(d$trt)

ctl <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 100)
equal <- lme(yield ~ gen * dust * dry * inoc, random = ~ 1 | block / gen / trt,
             data = d, method = "REML", control = ctl)
by_block <- update(equal, weights = varIdent(form = ~ 1 | block))
by_dust <- update(equal, weights = varIdent(form = ~ 1 | dust))

a0 <- anova(equal)
a1 <- anova(by_block)
a2 <- anova(by_dust)
tab <- data.frame(term = rownames(a0), numDF = a0$numDF, denDF = a0$denDF,
                  p_equal_variance = a0$`p-value`,
                  p_variance_by_block = a1$`p-value`,
                  p_variance_by_dust = a2$`p-value`)[-1, ]
tab$p_largest <- pmax(tab$p_equal_variance, tab$p_variance_by_block, tab$p_variance_by_dust)
tab$below_05_in_all_three <- tab$p_largest < 0.05
write.csv(tab, "outputs/unequal_variance.csv", row.names = FALSE)

sink("outputs/unequal_variance.txt")
cat("==== Likelihood-ratio tests against one residual variance (REML, same fixed effects) ====\n")
print(anova(equal, by_block))
print(anova(equal, by_dust))
cat("\n==== Residual SD relative to the reference level ====\n")
sd_block <- coef(by_block$modelStruct$varStruct, unconstrained = FALSE)
sd_dust <- coef(by_dust$modelStruct$varStruct, unconstrained = FALSE)
cat("by block (reference ", setdiff(levels(d$block), names(sd_block)), "):\n", sep = "")
print(round(sd_block, 2))
cat("by dust (reference ", setdiff(levels(d$dust), names(sd_dust)), "):\n", sep = "")
print(round(sd_dust, 2))
cat("\n==== Approximate 95% intervals for those SD ratios (1 = same variance as the reference) ====\n")
cat("by block:\n")
print(round(intervals(by_block, which = "var-cov")$varStruct, 2))
cat("by dust:\n")
print(round(intervals(by_dust, which = "var-cov")$varStruct, 2))
cat("\n==== Variance components under the three variance models (Residual = reference level) ====\n")
vc <- sapply(list(equal = equal, by_block = by_block, by_dust = by_dust), function(m) {
  v <- VarCorr(m)
  setNames(as.numeric(v[c(2, 4, 6, 7), "Variance"]), c("block", "whole plot", "subplot", "residual"))
})
print(round(vc, 2))
cat("\nNote: the by-block SD ratios are imprecise (every interval contains 1), and the F tests treat the\n",
    "estimated variances as known. Use the last table to see which conclusions do not depend on the\n",
    "variance model, not to claim new effects.\n", sep = "")
cat("\n==== Fixed-effect p-values under the three variance models ====\n")
show <- tab
show[, 4:7] <- round(show[, 4:7], 4)
print(show, row.names = FALSE)
sink()
