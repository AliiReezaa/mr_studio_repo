# Round 35: the parts of R21 (instrument count) and R23 (instrument strength) that run_calibration.R did not cover.
#   A  R21 under directional pleiotropy (InSIDE holds): 20% of the instruments (at least 1) have a positive pleiotropic effect.
#   B  R21 for MR-PRESSO (left out of Round 29 as too slow): global test and reported estimate at 4 to 20 instruments.
#   C  R23 with sample overlap between the exposure and outcome GWAS, where weak instruments bias toward the observational
#      association instead of the null (Burgess 2016, PMID 27625185).
# Same data generation, methods and settings as run_calibration.R (TwoSampleMR 0.7.13, RadialMR 1.2.5, MRPRESSO 1.0; nboot 1000,
# phi 1, NbDistribution 1000, SignifThreshold 0.05). IVW is TwoSampleMR::mr_ivw, the default since Round 30.
# Run from the repository root:  Rscript --vanilla validation/calibration/counts_strength_r35.R        (about 2 hours on 11 cores)
# Quick check of the script:     Rscript --vanilla validation/calibration/counts_strength_r35.R quick
suppressMessages({ library(TwoSampleMR); library(RadialMR); library(MRPRESSO); library(parallel) })
quick <- identical(commandArgs(TRUE)[1], "quick")
out_dir <- if (quick) tempdir() else file.path("validation", "calibration")
cores <- max(1L, detectCores() - 1L)
REPS <- if (quick) 6L else 500L; REPS_P <- if (quick) 4L else 200L
N <- 100000; PA <- default_parameters(); z <- qnorm(0.975)
RNGkind("L'Ecuyer-CMRG")

# Exposure and outcome errors correlate by rho_e = overlap fraction x observational correlation (equal N, standardised traits).
sim_pair <- function(mu, theta, alpha = NULL, rho_e = 0) {
  k <- length(mu); f <- runif(k, 0.05, 0.5); s <- 1 / sqrt(2 * f * (1 - f) * N)
  g <- mu * s * sample(c(-1, 1), k, TRUE); a <- if (is.null(alpha)) numeric(k) else alpha * s
  e1 <- rnorm(k); e2 <- rho_e * e1 + sqrt(1 - rho_e^2) * rnorm(k)
  list(bx = g + e1 * s, by = theta * g + a + e2 * s, sx = s, sy = s, g = g)
}
fit <- function(d, methods) {
  k <- length(d$bx); res <- list(); add <- function(m, b, se, p, crit) res[[m]] <<- c(b = b, se = se, p = p, crit = crit)
  if ("ivw" %in% methods) { r <- mr_ivw(d$bx, d$by, d$sx, d$sy, PA); add("ivw", r$b, r$se, r$pval, z) }
  if ("radial" %in% methods) { r <- summary(lm(I(d$by / d$sy) ~ -1 + I(d$bx / d$sy)))$coef; add("radial", r[1, 1], r[1, 2], r[1, 4], qt(0.975, k - 1)) }
  if ("egger" %in% methods) { r <- mr_egger_regression(d$bx, d$by, d$sx, d$sy, PA); add("egger", r$b, r$se, r$pval, qt(0.975, r$Q_df)) }
  if ("radegger" %in% methods && k >= 4) { r <- suppressWarnings(egger_radial(format_radial(d$bx, d$by, d$sx, d$sy, paste0("v", seq_len(k))), 0.05, 1L, FALSE))
    add("radegger", r$coef[2, 1], r$coef[2, 2], r$coef[2, 4], qt(0.975, r$df)) }
  if ("wm" %in% methods) { r <- mr_weighted_median(d$bx, d$by, d$sx, d$sy, PA); add("wm", r$b, r$se, r$pval, z) }
  if ("smed" %in% methods) { r <- mr_simple_median(d$bx, d$by, d$sx, d$sy, PA); add("smed", r$b, r$se, r$pval, z) }
  if ("mode" %in% methods) { r <- suppressWarnings(mr_mode(data.frame(beta.exposure = d$bx, beta.outcome = d$by, se.exposure = d$sx, se.outcome = d$sy), PA, mode_method = "all"))
    w <- r[r$method == "Weighted mode", ]; s <- r[r$method == "Simple mode", ]; add("mode", w$b, w$se, w$pval, z); add("smode", s$b, s$se, s$pval, z) }
  res
}
summ <- function(fits, theta, extra) { ms <- unique(unlist(lapply(fits, names)))
  do.call(rbind, lapply(ms, function(m) { x <- do.call(rbind, lapply(fits, `[[`, m)); x <- x[is.finite(x[, "b"]) & is.finite(x[, "se"]), , drop = FALSE]; rj <- mean(x[, "p"] < 0.05)
    data.frame(extra, method = m, reps = nrow(x), mean_est = mean(x[, "b"]), emp_sd = sd(x[, "b"]), mean_se = mean(x[, "se"]), reject_p05 = rj,
      coverage95 = mean(abs(x[, "b"] - theta) < x[, "crit"] * x[, "se"]), mc_se_reject = sqrt(rj * (1 - rj) / nrow(x))) })) }
cell <- function(seed, reps, gen, methods) { set.seed(seed); sd <- sample.int(.Machine$integer.max, reps)
  mclapply(sd, function(s) { set.seed(s); fit(gen(), methods) }, mc.cores = cores, mc.set.seed = FALSE) }
n_inv <- function(k) max(1L, round(0.2 * k))

# ---- A. R21 under directional pleiotropy -----------------------------------------------------------------------------
# theta = 0, mu ~ U(6, 12) as in Round 29; invalid instruments get alpha = |N(0, 3)| outcome SEs (all positive, mean 2.4 SE).
t0 <- Sys.time(); A <- NULL
for (k in c(3, 4, 5, 6, 8, 10, 15, 20, 30)) {
  gen <- function() { a <- numeric(k); j <- sample.int(k, n_inv(k)); a[j] <- abs(rnorm(length(j), 0, 3)); sim_pair(runif(k, 6, 12), 0, a) }
  A <- rbind(A, summ(cell(3500 + k, REPS, gen, c("ivw", "radial", "egger", "radegger", "wm", "smed", "mode")), 0, data.frame(k = k, scenario = "directional")))
  cat("A k =", k, format(Sys.time() - t0), "\n") }
write.csv(A, file.path(out_dir, "results_R21_directional.csv"), row.names = FALSE)

# ---- B. MR-PRESSO under R21 ------------------------------------------------------------------------------------------
# Scenarios as Round 29 (valid; balanced: 20% invalid with alpha ~ N(0, 3) SE). The reported estimate follows the adapter:
# the outlier-corrected row when outliers were removed, else the raw row; interval on n - outliers - 1 (or n - 1) df.
pn <- function(x) suppressWarnings(as.numeric(sub("^<", "", as.character(x))))
presso_one <- function(d) { k <- length(d$bx)
  f <- tryCatch(suppressWarnings(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", data = data.frame(bx = d$bx, by = d$by, sx = d$sx, sy = d$sy),
         OUTLIERtest = TRUE, DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 1000, seed = sample.int(1e6, 1))), error = function(e) NULL)
  if (is.null(f)) return(c(ok = 0, glob = NA, nout = NA, b = NA, se = NA, p = NA, crit = NA))
  main <- f[["Main MR results"]]; dt <- f[["MR-PRESSO results"]][["Distortion Test"]]
  no <- if (!is.null(dt) && is.numeric(dt[["Outliers Indices"]])) length(dt[["Outliers Indices"]]) else 0L
  row <- if (no > 0 && is.finite(main[2, "Causal Estimate"])) 2 else 1; df <- if (row == 2) k - no - 1 else k - 1
  c(ok = 1, glob = pn(f[["MR-PRESSO results"]][["Global Test"]]$Pvalue) < 0.05, nout = no, b = main[row, "Causal Estimate"], se = main[row, "Sd"],
    p = pn(main[row, "P-value"]), crit = qt(0.975, max(df, 1))) }
B <- NULL
for (scen in c("valid", "balanced")) for (k in c(4, 6, 8, 10, 15, 20)) {
  set.seed(3510 + k + 100 * (scen == "balanced"))
  r <- do.call(rbind, mclapply(seq_len(REPS_P), function(i) { a <- numeric(k); if (scen == "balanced") { j <- sample.int(k, n_inv(k)); a[j] <- rnorm(length(j), 0, 3) }
    presso_one(sim_pair(runif(k, 6, 12), 0, a)) }, mc.cores = cores, mc.set.seed = TRUE))
  r <- r[r[, "ok"] == 1 & is.finite(r[, "b"]) & is.finite(r[, "p"]), , drop = FALSE]; rj <- mean(r[, "p"] < 0.05)
  B <- rbind(B, data.frame(k = k, scenario = scen, reps = nrow(r), global_reject = mean(r[, "glob"]), mean_outliers = mean(r[, "nout"]), mean_est = mean(r[, "b"]),
    emp_sd = sd(r[, "b"]), reject_p05 = rj, coverage95 = mean(abs(r[, "b"]) < r[, "crit"] * r[, "se"]), mc_se_reject = sqrt(rj * (1 - rj) / nrow(r))))
  cat("B", scen, "k =", k, format(Sys.time() - t0), "\n") }
write.csv(B, file.path(out_dir, "results_R21_presso.csv"), row.names = FALSE)

# ---- C. R23 with sample overlap ---------------------------------------------------------------------------------------
# 30 instruments of fixed target F, theta = 0, observational correlation 0.3 from confounding, overlap 0, 50% or 100%.
# Expected IVW bias about overlap x 0.3 / mean F (Burgess 2016); bias_vs_obs below is the mean estimate divided by 0.3.
C <- NULL; r_obs <- 0.3
for (ov in c(0, 0.5, 1)) for (Ft in c(3, 5, 10, 20, 50)) {
  gen <- function() sim_pair(rep(sqrt(Ft - 1), 30), 0, NULL, ov * r_obs)
  x <- summ(cell(3600 + Ft + 1000 * ov, REPS, gen, c("ivw", "egger", "wm")), 0, data.frame(overlap = ov, target_f = Ft))
  x$bias_vs_obs <- x$mean_est / r_obs; x$expected_ivw <- ov * r_obs / Ft / r_obs; C <- rbind(C, x)
  cat("C overlap", ov, "F", Ft, format(Sys.time() - t0), "\n") }
write.csv(C, file.path(out_dir, "results_R23_overlap.csv"), row.names = FALSE)

pr <- function(x) print(format(x, digits = 3), row.names = FALSE)
pr(reshape(A[, c("k", "method", "reject_p05")], idvar = "k", timevar = "method", direction = "wide"))
pr(B); pr(reshape(C[, c("overlap", "target_f", "method", "reject_p05", "bias_vs_obs")], idvar = c("overlap", "target_f"), timevar = "method", direction = "wide"))
cat("Total time:", format(Sys.time() - t0), "\n")
