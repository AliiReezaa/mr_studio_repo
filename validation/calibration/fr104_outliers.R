# Calibration of the outlier-module rules of FR-90 to FR-99 (requirements_spec.md section H, validation item FR-104), Round 32.
# The module itself is NOT BUILT. This simulation measures, before it is built, what its [design] thresholds would do:
#   detectors   radial IVW and radial MR-Egger outliers (RadialMR 1.2.5, first-order weights) at the package level alpha = 0.05
#               and at the Bonferroni level 0.05 / k; MR-PRESSO (MRPRESSO 1.0, NbDistribution 1000, SignifThreshold 0.05) in a smaller run
#   trigger     FR-90: Cochran's Q P < 0.05 (IVW) or MR-Egger intercept P < 0.05
#   selection   FR-94 to FR-99: candidate sets per detector and their consensus, cap on the share removed (10%, 20% or 30%),
#               floor of 5 instruments, steps 1 to 3 of the selection rule, IVW (TwoSampleMR::mr_ivw, the app's default) on the chosen set
# Measures: false removal rate (valid variants flagged), power (invalid variants flagged), trigger rate, how often the cap or the floor
# decides, and bias, coverage and false positive rate of IVW on all instruments and on the selected set.
#
# Run from the repository root:  Rscript --vanilla validation/calibration/fr104_outliers.R          (about 30 minutes on 11 cores)
# Quick check of the script:      Rscript --vanilla validation/calibration/fr104_outliers.R quick
suppressMessages({ library(TwoSampleMR); library(RadialMR); library(MRPRESSO); library(parallel) })
quick <- identical(commandArgs(TRUE)[1], "quick")
out_dir <- if (quick) tempdir() else file.path("validation", "calibration")
cores <- max(1L, detectCores() - 1L)
REPS <- if (quick) 10L else 1000L          # fast detectors; Monte Carlo SE of a 5% rate is 0.7 percentage points
REPS_P <- if (quick) 2L else c(`10` = 100L, `30` = 60L)   # MR-PRESSO (18 s per data set with 10 variants, 60 s with 30)
N <- 100000
PA <- default_parameters()
RNGkind("L'Ecuyer-CMRG")

# ---- Data: two non-overlapping GWAS of N people, standardised traits, true strength mu_j ~ U(4, 10) (F about 17 to 101) ----
# Pleiotropic effect alpha_j is in units of the variant's outcome SE.
#   valid        all alpha_j = 0
#   hetero       valid but heterogeneous: each variant has its own causal effect theta_j ~ N(theta, 0.05^2) (biology, not pleiotropy)
#   balanced     30% invalid, alpha_j ~ N(0, 3^2)
#   directional  30% invalid, alpha_j ~ N(3, 1)
#   few_strong   5% invalid (at least 1), alpha_j = +-8
#   many_strong  30% invalid, alpha_j = +-8
sim <- function(k, theta, scen) {
  f <- runif(k, 0.05, 0.5); s <- 1 / sqrt(2 * f * (1 - f) * N); g <- runif(k, 4, 10) * s * sample(c(-1, 1), k, TRUE)
  th <- if (scen == "hetero") rnorm(k, theta, 0.05) else rep(theta, k)
  ni <- switch(scen, valid = 0L, hetero = 0L, balanced = round(0.3 * k), directional = round(0.3 * k), few_strong = max(1L, round(0.05 * k)), many_strong = round(0.3 * k))
  inv <- if (ni) sample.int(k, ni) else integer(0); a <- numeric(k)
  if (ni) a[inv] <- switch(scen, balanced = rnorm(ni, 0, 3), directional = rnorm(ni, 3, 1), few_strong = , many_strong = 8 * sample(c(-1, 1), ni, TRUE)) * s[inv]
  list(bx = rnorm(k, g, s), by = rnorm(k, th * g + a, s), sx = s, sy = s, invalid = seq_len(k) %in% inv)
}

# ---- Detectors, first-order weights, as RadialMR defines them (checked against RadialMR in check_radial()) ----
radial_ivw_q <- function(d) { w <- d$bx^2 / d$sy^2; b <- d$by / d$bx; bi <- sum(w * b) / sum(w); w * (b - bi)^2 }
radial_egger_q <- function(d) { w <- d$bx^2 / d$sy^2; b <- d$by / d$bx; x <- sqrt(w); y <- b * sqrt(w)   # the ratio b is the same for either allele
  fit <- lm(y ~ x); stats::residuals(fit)^2 }
flag <- function(q, alpha) stats::pchisq(q, 1, lower.tail = FALSE) < alpha
check_radial <- function() { set.seed(32); dif <- 0L
  for (i in 1:30) { d <- sim(20, 0.2, "few_strong"); rf <- format_radial(d$bx, d$by, d$sx, d$sy, paste0("v", 1:20))
    ri <- suppressWarnings(ivw_radial(rf, 0.05, 1L, 1e-4, FALSE)); re <- suppressWarnings(egger_radial(rf, 0.05, 1L, FALSE))
    oi <- if (is.data.frame(ri$outliers)) as.character(ri$outliers$SNP) else character(0); oe <- if (is.data.frame(re$outliers)) as.character(re$outliers$SNP) else character(0)
    dif <- dif + length(setdiff(union(oi, paste0("v", which(flag(radial_ivw_q(d), 0.05)))), intersect(oi, paste0("v", which(flag(radial_ivw_q(d), 0.05)))))) +
      length(setdiff(union(oe, paste0("v", which(flag(radial_egger_q(d), 0.05)))), intersect(oe, paste0("v", which(flag(radial_egger_q(d), 0.05)))))) }
  dif }

ivw <- function(d, keep) { r <- mr_ivw(d$bx[keep], d$by[keep], d$sx[keep], d$sy[keep], PA); c(b = r$b, se = r$se, p = r$pval) }
signals <- function(d, keep) { k <- sum(keep); if (k < 3) return(c(q = 1, e = 1, Q = NA))
  qi <- mr_ivw_mre(d$bx[keep], d$by[keep], d$sx[keep], d$sy[keep], PA); eg <- mr_egger_regression(d$bx[keep], d$by[keep], d$sx[keep], d$sy[keep], PA)
  c(q = stats::pchisq(qi$Q, qi$Q_df, lower.tail = FALSE), e = eg$pval_i, Q = qi$Q) }

# FR-94 to FR-99 on one data set. sets: named list of logical "flagged" vectors (one per detector).
select_set <- function(d, sets, cap, floor_n = 5L) {
  k <- length(d$bx); s0 <- signals(d, rep(TRUE, k)); trig <- s0[["q"]] < 0.05 || s0[["e"]] < 0.05
  none <- list(keep = rep(TRUE, k), trig = trig, removed = 0L, why = if (trig) "no eligible improvement" else "no trigger", capped = FALSE)
  if (!trig) return(none)
  cand <- sets; if (length(sets) >= 2) { agree <- Reduce(`+`, lapply(sets, as.integer)); cand$consensus <- agree >= 2 }
  cand <- Filter(function(f) any(f), cand); if (!length(cand)) return(none)
  agree <- Reduce(`+`, lapply(sets, as.integer))
  ev <- lapply(names(cand), function(nm) { f <- cand[[nm]]; nr <- sum(f); keep <- !f
    eligible <- nr <= floor(cap * k) && sum(keep) >= floor_n
    s <- if (eligible) signals(d, keep) else c(q = NA, e = NA, Q = NA)
    list(nm = nm, keep = keep, nr = nr, eligible = eligible, over_cap = nr > floor(cap * k), clears = eligible && s[["q"]] >= 0.05 && s[["e"]] >= 0.05,
         agreement = mean(agree[f]), excess = if (eligible) (s0[["Q"]] - s[["Q"]]) - nr else NA_real_) })
  capped <- any(vapply(ev, `[[`, TRUE, "over_cap")); el <- Filter(function(e) e$eligible, ev)
  if (!length(el)) { none$capped <- capped; return(none) }
  cl <- Filter(function(e) e$clears, el)
  pick <- if (length(cl)) { o <- order(vapply(cl, `[[`, 0, "nr"), -vapply(cl, `[[`, 0, "agreement"), -vapply(cl, `[[`, 0, "excess")); cl[[o[1]]] }
          else { ex <- vapply(el, `[[`, 0, "excess"); if (any(ex > 0)) el[[which.max(ex)]] else NULL }
  if (is.null(pick)) { none$capped <- capped; return(none) }
  list(keep = pick$keep, trig = TRUE, removed = pick$nr, why = if (pick$clears) "step 2" else "step 3", capped = capped)
}

one_rep <- function(k, theta, scen) {
  d <- sim(k, theta, scen); qi <- radial_ivw_q(d); qe <- radial_egger_q(d)
  det <- list(rivw_05 = flag(qi, 0.05), rivw_bonf = flag(qi, 0.05 / k), regg_05 = flag(qe, 0.05), regg_bonf = flag(qe, 0.05 / k))
  rates <- unlist(lapply(names(det), function(n) c(setNames(if (any(!d$invalid)) mean(det[[n]][!d$invalid]) else NA_real_, paste0(n, "_false")),
                                                   setNames(if (any(d$invalid)) mean(det[[n]][d$invalid]) else NA_real_, paste0(n, "_power")))))
  full <- ivw(d, rep(TRUE, k))
  rules <- list(r05_cap20 = list(c("rivw_05", "regg_05"), 0.2), bonf_cap10 = list(c("rivw_bonf", "regg_bonf"), 0.1),
                bonf_cap20 = list(c("rivw_bonf", "regg_bonf"), 0.2), bonf_cap30 = list(c("rivw_bonf", "regg_bonf"), 0.3))
  sel <- unlist(lapply(names(rules), function(rn) { ru <- rules[[rn]]; s <- select_set(d, det[ru[[1]]], ru[[2]])
    e <- if (sum(s$keep) >= 3) ivw(d, s$keep) else c(b = NA, se = NA, p = NA)
    setNames(c(e[["b"]], e[["se"]], e[["p"]], s$trig, s$removed, s$removed > 0, s$capped, if (any(d$invalid)) mean(!s$keep[d$invalid]) else NA, if (any(!d$invalid)) mean(!s$keep[!d$invalid]) else NA),
             paste0(rn, c("_b", "_se", "_p", "_trig", "_removed", "_any", "_capped", "_inv_removed", "_valid_removed"))) }))
  c(full_b = full[["b"]], full_se = full[["se"]], full_p = full[["p"]], rates, sel)
}

cells <- expand.grid(scen = c("valid", "hetero", "balanced", "directional", "few_strong", "many_strong"), k = c(10L, 30L, 100L, 200L), theta = c(0, 0.2), stringsAsFactors = FALSE)
cat("radial detector check, differing flags in 30 data sets:", check_radial(), "\n")
set.seed(3201)
res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) { ce <- cells[i, ]
  r <- do.call(rbind, mclapply(seq_len(REPS), function(j) one_rep(ce$k, ce$theta, ce$scen), mc.cores = cores, mc.set.seed = TRUE))
  m <- colMeans(r, na.rm = TRUE); z <- qnorm(0.975)
  cov <- function(b, se) mean(abs(r[, b] - ce$theta) <= z * r[, se], na.rm = TRUE)
  rej <- function(p) mean(r[, p] < 0.05, na.rm = TRUE)
  bias <- function(b) if (ce$theta == 0) mean(r[, b], na.rm = TRUE) else mean(r[, b], na.rm = TRUE) / ce$theta - 1
  out <- data.frame(scenario = ce$scen, k = ce$k, theta = ce$theta, reps = REPS,
    full_bias = bias("full_b"), full_coverage = cov("full_b", "full_se"), full_reject = rej("full_p"))
  for (dn in c("rivw_05", "rivw_bonf", "regg_05", "regg_bonf")) { out[[paste0(dn, "_false")]] <- m[[paste0(dn, "_false")]]; out[[paste0(dn, "_power")]] <- m[[paste0(dn, "_power")]] }
  for (rn in c("r05_cap20", "bonf_cap10", "bonf_cap20", "bonf_cap30")) {
    out[[paste0(rn, "_trigger")]] <- m[[paste0(rn, "_trig")]]; out[[paste0(rn, "_removal_rate")]] <- m[[paste0(rn, "_any")]]
    out[[paste0(rn, "_mean_removed")]] <- m[[paste0(rn, "_removed")]]; out[[paste0(rn, "_cap_decided")]] <- m[[paste0(rn, "_capped")]]
    out[[paste0(rn, "_invalid_removed")]] <- m[[paste0(rn, "_inv_removed")]]; out[[paste0(rn, "_valid_removed")]] <- m[[paste0(rn, "_valid_removed")]]
    out[[paste0(rn, "_bias")]] <- bias(paste0(rn, "_b")); out[[paste0(rn, "_coverage")]] <- cov(paste0(rn, "_b"), paste0(rn, "_se")); out[[paste0(rn, "_reject")]] <- rej(paste0(rn, "_p")) }
  cat(sprintf("%-11s k=%3d theta=%.1f done\n", ce$scen, ce$k, ce$theta)); out }))
utils::write.csv(res, file.path(out_dir, "results_fr104_outliers.csv"), row.names = FALSE)

# ---- MR-PRESSO, smaller run: false removal and power of its outlier test at 10 and 30 variants ----
pcells <- expand.grid(scen = c("valid", "balanced", "few_strong"), k = c(10L, 30L), stringsAsFactors = FALSE)
set.seed(3202)
pres <- do.call(rbind, lapply(seq_len(nrow(pcells)), function(i) { ce <- pcells[i, ]; reps <- if (quick) REPS_P else REPS_P[[as.character(ce$k)]]
  r <- do.call(rbind, mclapply(seq_len(reps), function(j) { d <- sim(ce$k, 0.2, ce$scen)
    pd <- data.frame(bx = d$bx, by = d$by, sx = d$sx, sy = d$sy)
    fit <- tryCatch(suppressWarnings(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", data = pd, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    SignifThreshold = 0.05, NbDistribution = 1000, seed = sample.int(1e6, 1))), error = function(e) NULL)
    if (is.null(fit)) return(c(ok = 0, glob = NA, fl = NA, pw = NA))
    gp <- suppressWarnings(as.numeric(sub("^<", "", as.character(fit$`MR-PRESSO results`$`Global Test`$Pvalue))))
    idx <- fit$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`; fl <- rep(FALSE, ce$k); if (is.numeric(idx)) fl[idx] <- TRUE
    c(ok = 1, glob = gp < 0.05, fl = if (any(!d$invalid)) mean(fl[!d$invalid]) else NA, pw = if (any(d$invalid)) mean(fl[d$invalid]) else NA) }, mc.cores = cores, mc.set.seed = TRUE))
  cat(sprintf("presso %-11s k=%3d done\n", ce$scen, ce$k))
  data.frame(scenario = ce$scen, k = ce$k, theta = 0.2, reps = nrow(r), ran = mean(r[, "ok"]), global_p_below_0.05 = mean(r[, "glob"], na.rm = TRUE),
             false_removal = mean(r[, "fl"], na.rm = TRUE), power = mean(r[, "pw"], na.rm = TRUE)) }))
utils::write.csv(pres, file.path(out_dir, "results_fr104_presso.csv"), row.names = FALSE)
print(res[, 1:8]); print(pres)
