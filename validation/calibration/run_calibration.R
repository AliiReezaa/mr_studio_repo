# Threshold calibration by simulation (Round 29, open item 4 / FR-104).
# Calibrates the three heuristic thresholds that are active in the v1 feasibility engine (prototype/engine/feasibility_rules.csv):
#   R21  at least 10 instruments for MR-Egger, the median, mode and radial methods          (rule of thumb)
#   R23  mean F statistic of at least 10                                                    (Burgess 2011, PMID 21414999)
#   R24  I2GX of at least 0.9 for MR-Egger and radial MR-Egger                              (Bowden 2016, PMID 27616674)
# Every estimate is computed with the functions and default settings the MR Studio adapters call (TwoSampleMR 0.7.13, RadialMR 1.2.5:
# nboot = 1000, phi = 1, first-order weights), and every P value and interval is formed as the adapters form it.
# The data are simulated summary statistics for two non-overlapping samples (exposure and outcome GWAS of 100,000 people each).
# MR-PRESSO is not included (about a minute per data set with its default 1000 simulated data sets); its R21 threshold stays uncalibrated.
# Run from the repository root (about two hours on 11 cores):  Rscript --vanilla validation/calibration/run_calibration.R
# Quick check of the script (a few repetitions):              Rscript --vanilla validation/calibration/run_calibration.R quick
suppressMessages({ library(TwoSampleMR); library(RadialMR); library(parallel) })
quick <- identical(commandArgs(TRUE)[1], "quick")
out_dir <- if (quick) tempdir() else file.path("validation", "calibration")
cores <- max(1L, detectCores() - 1L)
REPS <- if (quick) 8L else 500L        # R21 and R23: Monte Carlo SE of a 5% rejection rate with 500 repetitions is 1.0 percentage point
REPS_I2 <- if (quick) 20L else 2000L   # R24: fast methods only
N <- 100000                            # people in each GWAS
PA <- default_parameters()             # nboot = 1000, the adapters' default; phi = 1 is mr_mode's default
RNGkind("L'Ecuyer-CMRG")

# ---- Data generation ----------------------------------------------------------------------------------------------
# Variant j has allele frequency f_j ~ U(0.05, 0.5); with standardised traits its standard error is 1 / sqrt(2 f (1 - f) N).
# mu_j is its true strength (true effect / SE), so its expected F statistic is mu_j^2 + 1.
# Outcome effect: theta * gamma_j + alpha_j (alpha_j = 0 for a valid instrument) plus sampling noise.
sim_pair <- function(mu, theta, alpha_sd_in_se = 0, n_invalid = 0) {
  k <- length(mu); f <- runif(k, 0.05, 0.5); s <- 1 / sqrt(2 * f * (1 - f) * N)
  g <- mu * s * sample(c(-1, 1), k, TRUE)
  a <- numeric(k); if (n_invalid > 0) a[sample.int(k, n_invalid)] <- rnorm(n_invalid, 0, alpha_sd_in_se) * s
  list(bx = rnorm(k, g, s), by = rnorm(k, theta * g + a, s), sx = s, sy = s, g = g)
}

# ---- The methods, as the adapters run them --------------------------------------------------------------------------
# Radial IVW (first order, weights 1) has the IVW point estimate and the residual-scaled IVW SE with a t interval on k - 1 df;
# it is computed directly here because RadialMR::ivw_radial also iterates exact estimates (3 s per call). check_radial() below
# confirms on 50 data sets that the numbers equal RadialMR's first row, which is the row the adapter reports.
fit_all <- function(d, methods) {
  k <- length(d$bx); z <- qnorm(0.975); res <- list()
  add <- function(m, b, se, p, crit) res[[m]] <<- c(b = b, se = se, p = p, crit = crit)
  if ("ivw_mre" %in% methods) { r <- mr_ivw_mre(d$bx, d$by, d$sx, d$sy, PA); add("ivw_mre", r$b, r$se, r$pval, z) }
  if ("radial" %in% methods) { r <- summary(lm(I(d$by / d$sy) ~ -1 + I(d$bx / d$sy)))$coef; add("radial", r[1, 1], r[1, 2], r[1, 4], qt(0.975, k - 1)) }
  if ("egger" %in% methods && k >= 3) { r <- mr_egger_regression(d$bx, d$by, d$sx, d$sy, PA); add("egger", r$b, r$se, r$pval, qt(0.975, r$Q_df)) }
  if ("radegger" %in% methods && k >= 4) { r <- suppressWarnings(egger_radial(format_radial(d$bx, d$by, d$sx, d$sy, paste0("v", seq_len(k))), 0.05, 1L, FALSE))
    add("radegger", r$coef[2, 1], r$coef[2, 2], r$coef[2, 4], qt(0.975, r$df)) }
  if ("wm" %in% methods) { r <- mr_weighted_median(d$bx, d$by, d$sx, d$sy, PA); add("wm", r$b, r$se, r$pval, z) }
  if ("smed" %in% methods) { r <- mr_simple_median(d$bx, d$by, d$sx, d$sy, PA); add("smed", r$b, r$se, r$pval, z) }
  if (any(c("mode", "smode") %in% methods)) {   # mr_mode computes the simple and the weighted mode in one bootstrap
    r <- suppressWarnings(mr_mode(data.frame(beta.exposure = d$bx, beta.outcome = d$by, se.exposure = d$sx, se.outcome = d$sy), PA, mode_method = "all"))
    w <- r[r$method == "Weighted mode", ]; s <- r[r$method == "Simple mode", ]
    add("mode", w$b, w$se, w$pval, z); add("smode", s$b, s$se, s$pval, z) }
  res
}
check_radial <- function() { set.seed(29); dif <- 0
  for (i in 1:50) { d <- sim_pair(runif(20, 4, 10), 0.2); a <- fit_all(d, "radial")$radial
    r <- suppressWarnings(ivw_radial(format_radial(d$bx, d$by, d$sx, d$sy, paste0("v", 1:20)), 0.05, 1L, 1e-4, FALSE))$coef[1, ]
    dif <- max(dif, abs(a[["b"]] - r[[1]]), abs(a[["se"]] - r[[2]]), abs(a[["p"]] - r[[4]])) }
  dif }

# Summary of one cell: each method's mean estimate, relative bias, empirical SD, mean reported SE, rejection rate of the null
# (its own P < 0.05) and coverage of its 95% interval around the true theta.
summarise <- function(fits, theta, extra = list()) {
  ms <- unique(unlist(lapply(fits, names)))
  do.call(rbind, lapply(ms, function(m) { x <- do.call(rbind, lapply(fits, `[[`, m)); x <- x[is.finite(x[, "b"]) & is.finite(x[, "se"]), , drop = FALSE]
    data.frame(extra, method = m, reps = nrow(x), mean_est = mean(x[, "b"]), rel_bias = if (theta != 0) mean(x[, "b"]) / theta - 1 else NA_real_,
      emp_sd = sd(x[, "b"]), mean_se = mean(x[, "se"]), reject_p05 = mean(x[, "p"] < 0.05), coverage95 = mean(abs(x[, "b"] - theta) < x[, "crit"] * x[, "se"]),
      mc_se_reject = sqrt(mean(x[, "p"] < 0.05) * (1 - mean(x[, "p"] < 0.05)) / nrow(x)), stringsAsFactors = FALSE) }))
}
run_cell <- function(seed, reps, gen, methods) { set.seed(seed); seeds <- sample.int(.Machine$integer.max, reps)
  mclapply(seeds, function(s) { set.seed(s); d <- gen(); c(fit_all(d, methods), list(.meta = c(mean_f = mean((d$bx / d$sx)^2), true_f = mean((d$g / d$sx)^2) + 1,
    i2gx = suppressWarnings(Isq(abs(d$bx), d$sx)), k = length(d$bx)))) }, mc.cores = cores, mc.set.seed = FALSE) }
meta <- function(fits) colMeans(do.call(rbind, lapply(fits, `[[`, ".meta")))
strip <- function(fits) lapply(fits, function(f) f[names(f) != ".meta"])

t0 <- Sys.time()
cat(sprintf("Radial IVW computed directly vs RadialMR::ivw_radial, largest difference over 50 data sets: %.2e\n", check_radial()))
ALL <- c("ivw_mre", "egger", "wm", "smed", "mode", "smode", "radial", "radegger")

# ---- R21: number of instruments -------------------------------------------------------------------------------------
# Strong, valid instruments (mu ~ U(6, 12): expected F 37 to 145), theta = 0, so every rejection is a false positive.
# Scenario "valid": all instruments valid. Scenario "balanced": 20% of instruments (rounded, at least 1) have a direct effect on the
# outcome with mean 0 and SD 3 outcome standard errors (balanced pleiotropy, InSIDE holds).
K <- c(3, 4, 5, 6, 8, 10, 15, 20, 30)
r21 <- list(); cell <- 0L
for (sc in c("valid", "balanced")) for (k in K) { cell <- cell + 1L
  ni <- if (sc == "valid") 0 else max(1, round(0.2 * k))
  f <- run_cell(2100 + cell, REPS, function() sim_pair(runif(k, 6, 12), 0, 3, ni), ALL)
  r21[[cell]] <- summarise(strip(f), 0, list(scenario = sc, k = k, n_invalid = ni))
  cat(sprintf("R21 %-8s k = %2d done (%s)\n", sc, k, format(round(Sys.time() - t0, 1)))) }
r21 <- do.call(rbind, r21)
utils::write.csv(r21, file.path(out_dir, "results_R21_instruments.csv"), row.names = FALSE)

# ---- R23: instrument strength (mean F) ------------------------------------------------------------------------------
# Part A: 30 valid instruments of set strength, without selection: mu_j ~ U(0.5 m, 1.5 m), with m chosen so that the expected F is
# 2, 3, 5, 10, 20 or 50. theta = 0.2 (bias and coverage) and theta = 0 (false positives; the mode methods, 7 s per data set, are left
# out at theta = 0).
# Part B: instruments chosen as MR Studio chooses them, by P < 5e-8 in the same exposure sample: 2,000 candidate variants with
# true strength mu ~ exponential (mean 1, 1.5, 2 or 3), of which those with |z| > 5.45 are kept (data sets with fewer than 3 are redrawn).
# The observed F of a kept variant is then at least 29.7 whatever its true strength (winner's curse). Fast methods only.
r23 <- list(); cell <- 0L
for (Fe in c(2, 3, 5, 10, 20, 50)) for (th in c(0.2, 0)) { cell <- cell + 1L; m <- sqrt((Fe - 1) / (13 / 12))
  f <- run_cell(2300 + cell, REPS, function() sim_pair(runif(30, 0.5 * m, 1.5 * m), th), if (th == 0) setdiff(ALL, c("mode", "smode")) else ALL); mt <- meta(f)
  r23[[cell]] <- summarise(strip(f), th, list(part = "A_fixed_strength", target_f = Fe, theta = th, mean_obs_f = unname(mt["mean_f"]), mean_true_f = unname(mt["true_f"]), mean_k = 30))
  cat(sprintf("R23 A F = %2d theta = %.1f done (%s)\n", Fe, th, format(round(Sys.time() - t0, 1)))) }
zsel <- qnorm(1 - 5e-8 / 2)
gen_sel <- function(mm, th) function() { repeat { mu <- rexp(2000, 1 / mm); d <- sim_pair(mu, th); keep <- abs(d$bx / d$sx) > zsel
  if (sum(keep) >= 3) return(lapply(d, `[`, keep)) } }
for (mm in c(1, 1.5, 2, 3)) for (th in c(0.2, 0)) { cell <- cell + 1L
  f <- run_cell(2300 + cell, REPS, gen_sel(mm, th), c("ivw_mre", "egger", "wm", "smed", "radial")); mt <- meta(f)
  r23[[cell]] <- summarise(strip(f), th, list(part = "B_selected_p5e-8", target_f = NA_real_, theta = th, mean_obs_f = unname(mt["mean_f"]), mean_true_f = unname(mt["true_f"]), mean_k = unname(mt["k"])))
  cat(sprintf("R23 B mean mu = %.1f theta = %.1f done (%s)\n", mm, th, format(round(Sys.time() - t0, 1)))) }
r23 <- do.call(rbind, r23)
utils::write.csv(r23, file.path(out_dir, "results_R23_strength.csv"), row.names = FALSE)

# ---- R24: I2GX and the dilution of MR-Egger ---------------------------------------------------------------------------
# 30 valid instruments, theta = 0.2, mu_j = |10 + s * N(0, 1)|: the spread s sets I2GX (about s^2 / (s^2 + 1)).
# Dilution = 1 - mean Egger slope / theta. Each data set is also kept on its own, so the dilution can be read by observed I2GX.
r24 <- list(); r24d <- list(); cell <- 0L
for (s in c(0.5, 1, 1.5, 2, 3, 4, 6)) { cell <- cell + 1L
  f <- run_cell(2400 + cell, REPS_I2, function() sim_pair(abs(10 + s * rnorm(30)), 0.2), c("ivw_mre", "egger", "radegger")); mt <- meta(f)
  r24[[cell]] <- summarise(strip(f), 0.2, list(spread = s, mean_i2gx = unname(mt["i2gx"])))
  r24d[[cell]] <- data.frame(spread = s, i2gx = vapply(f, function(x) x$.meta[["i2gx"]], 0), egger = vapply(f, function(x) x$egger[["b"]], 0), radegger = vapply(f, function(x) x$radegger[["b"]], 0))
  cat(sprintf("R24 spread = %.1f done (%s)\n", s, format(round(Sys.time() - t0, 1)))) }
r24 <- do.call(rbind, r24); r24d <- do.call(rbind, r24d)
br <- c(-Inf, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, Inf); r24d$bin <- cut(r24d$i2gx, br, right = FALSE)
r24b <- do.call(rbind, lapply(split(r24d, r24d$bin), function(x) if (nrow(x)) data.frame(i2gx_bin = as.character(x$bin[1]), n = nrow(x), mean_i2gx = mean(x$i2gx),
  egger_dilution = 1 - mean(x$egger) / 0.2, radegger_dilution = 1 - mean(x$radegger) / 0.2)))
utils::write.csv(r24, file.path(out_dir, "results_R24_i2gx.csv"), row.names = FALSE)
utils::write.csv(r24b, file.path(out_dir, "results_R24_i2gx_bins.csv"), row.names = FALSE)

cat(sprintf("\nFinished in %s. %s, TwoSampleMR %s, RadialMR %s, %d cores.\n", format(round(Sys.time() - t0, 1)), R.version.string, packageVersion("TwoSampleMR"), packageVersion("RadialMR"), cores))
pr <- function(x) print(format(x, digits = 3), row.names = FALSE)
cat("\nR21, false positive rate (P < 0.05) at theta = 0:\n"); pr(reshape(r21[, c("scenario", "k", "method", "reject_p05")], idvar = c("scenario", "k"), timevar = "method", direction = "wide"))
cat("\nR23, relative bias at theta = 0.2:\n"); pr(reshape(r23[r23$theta == 0.2, c("part", "target_f", "mean_obs_f", "mean_true_f", "method", "rel_bias")], idvar = c("part", "target_f", "mean_obs_f", "mean_true_f"), timevar = "method", direction = "wide"))
cat("\nR24, dilution of MR-Egger by observed I2GX:\n"); pr(r24b)
