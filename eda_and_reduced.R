# Run from the submission folder with Rscript eda_and_reduced.R.
library(nlme)
dir.create("outputs", showWarnings = FALSE)
d <- read.csv("data/goulden_splitsplit.csv", stringsAsFactors = TRUE)
d$trt <- factor(d$trt)
print(summary(d$yield))
print(table(table(d$gen, d$dust, d$dry, d$inoc)))
wide <- reshape(d[c("block","gen","trt","dust","inoc","yield")],
                idvar = c("block","gen","trt","dust"), timevar = "inoc", direction = "wide")
wide$gap <- wide$yield.Inoc - wide$yield.Uninoc
stopifnot(nrow(wide) == 80, !anyNA(wide$gap))
print(aggregate(gap ~ dust, wide, mean))
png("outputs/paired_eda.png", 1400, 575, res = 135)
par(mar = c(3.1, 4.1, .6, .6), las = 1)
dusts <- c("Check","CaCO3","Ceresan","DuBay","Semesan")
plot(0, type = "n", xlim = c(.55,5.45), ylim = range(c(0,wide$gap)) + c(-2,2),
     xaxt = "n", xlab = "", ylab = "Paired yield difference", bty = "l")
abline(h = 0, col = "gray70", lty = 2)
axis(1, 1:5, c("Untreated", dusts[-1]))
for (i in seq_along(dusts)) {
  v <- sort(wide$gap[wide$dust == dusts[i]])
  points(i + seq(-.17,.17,length.out=length(v)), v, pch=16, col="#808A93", cex=.8)
  segments(i-.22, mean(v), i+.22, mean(v), col="#155894", lwd=4)
}
dev.off()
# Retain all three random intercepts in both fits. Only fixed interactions change.
full <- lme(yield ~ gen*dust*dry*inoc, random = ~1 | block/gen/trt, data=d, method="REML")
reduced <- update(full, fixed = yield ~ gen*dry + dry*inoc + dust*inoc)
for (fit in list(full, reduced)) {
  print(VarCorr(fit))
  print(anova(fit))
}
