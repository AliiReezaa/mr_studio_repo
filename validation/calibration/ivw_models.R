# IVW effects models with few variants (Round 30). Compares the false positive rate at theta = 0 of
#   mr_ivw_mre : multiplicative random effects with no lower limit on the residual SE (TwoSampleMR::mr_ivw_mre)
#   mr_ivw     : the same, with the residual SE not allowed below 1 (TwoSampleMR::mr_ivw, the TwoSampleMR default; MR Studio's default from Round 30)
#   mr_burgess : MendelianRandomization::mr_ivw defaults (fixed effect below 4 variants, otherwise as mr_ivw)
#   mr_ivw_mre_t : mr_ivw_mre with a t distribution on k - 1 degrees of freedom
#   fixed      : fixed effect (TwoSampleMR::mr_ivw_fe)
# Same data generation as run_calibration.R (R21 cells). The closed forms are checked against MendelianRandomization::mr_ivw first.
# Run from the repository root (a few minutes):  Rscript --vanilla validation/calibration/ivw_models.R
suppressMessages(library(MendelianRandomization))
N <- 1e5
sim_pair <- function(mu, theta, alpha_sd_in_se = 0, n_invalid = 0) {
  k <- length(mu); f <- runif(k, 0.05, 0.5); s <- 1 / sqrt(2 * f * (1 - f) * N)
  g <- mu * s * sample(c(-1, 1), k, TRUE)
  a <- numeric(k); if (n_invalid > 0) a[sample.int(k, n_invalid)] <- rnorm(n_invalid, 0, alpha_sd_in_se) * s
  list(bx = rnorm(k, g, s), by = rnorm(k, theta * g + a, s), sx = s, sy = s)
}
ivw <- function(d) { w <- 1/d$sy^2; b <- sum(w*d$bx*d$by)/sum(w*d$bx^2); sefe <- 1/sqrt(sum(w*d$bx^2)); k <- length(d$bx)
  rse <- sqrt(sum(w*(d$by-b*d$bx)^2)/(k-1))
  c(mr_ivw_mre = 2*pnorm(-abs(b/(sefe*rse))), mr_ivw = 2*pnorm(-abs(b/(sefe*max(rse,1)))), mr_burgess = 2*pnorm(-abs(b/(sefe*(if (k<4) 1 else max(rse,1))))), mr_ivw_mre_t = 2*pt(-abs(b/(sefe*rse)), k-1), fixed = 2*pnorm(-abs(b/sefe))) }
# agreement with MendelianRandomization::mr_ivw defaults
set.seed(1); dif <- 0
for (i in 1:200) { k <- sample(3:12,1); d <- sim_pair(runif(k,6,12), 0, 3, max(1,round(.2*k))); m <- mr_ivw(mr_input(d$bx,d$sx,d$by,d$sy)); dif <- max(dif, abs(m@Pvalue - ivw(d)[["mr_burgess"]])) }
cat(sprintf("closed form vs MendelianRandomization::mr_ivw default, largest P difference over 200 data sets: %.1e\n", dif))
set.seed(29); R <- 20000; out <- NULL
for (sc in c("valid","balanced")) for (k in c(3,4,5,6,8,10,15,20,30)) {
  p <- t(replicate(R, { d <- if (sc=="valid") sim_pair(runif(k,6,12),0) else sim_pair(runif(k,6,12),0,3,max(1,round(.2*k))); ivw(d) }))
  out <- rbind(out, data.frame(scenario=sc, k=k, t(colMeans(p < 0.05)))) }
print(out, digits=3, row.names=FALSE); utils::write.csv(out, file.path("validation", "calibration", "results_ivw_models.csv"), row.names = FALSE); cat("reps per cell", R, "; MC SE near 0.05 =", round(sqrt(.05*.95/R),4), "\n")
