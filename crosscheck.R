# Optional cross-check of analysis.py in base R (no extra packages: stats::aov and nlme::lme).
# Run from the repository root:  Rscript crosscheck.R
d <- read.csv("data/goulden_splitsplit.csv", stringsAsFactors = TRUE)
d$trt <- factor(d$trt)

cat("==== aov with Error(block/gen/trt) ====\n")
print(summary(aov(yield ~ gen * dust * dry * inoc + Error(block / gen / trt), data = d)))

cat("\n==== nlme::lme, REML, random intercepts for block, whole plot, subplot ====\n")
library(nlme)
m <- lme(yield ~ gen * dust * dry * inoc, random = ~ 1 | block / gen / trt, data = d, method = "REML")
print(VarCorr(m))
print(anova(m))
