# Validation study 12 (Round 62, redone in Round 70e): BWMR (Zhao et al. 2020, Bioinformatics 36:1501-1508, doi:10.1093/bioinformatics/btz749,
# PMID 31593215). MR Studio runs its own implementation of the published algorithm (BWMR_MRS_SRC in contracts/r/adapters/reference_adapters_v11.R,
# written from the paper and its supplement, arXiv 1811.10223v2), because the authors' package carries no licence and cannot be distributed.
# The package (GitHub jiazhao97/BWMR, commit c00e818, with the authors' R >= 4.2 fix from updates/BWMR_updated.R applied at run time) is
# used here only as the reference to compare against. Data: the package's ExampleData (LDL cholesterol on dyslipidemia, 76 variants), the
# Waterworth et al. lipid example of MendelianRandomization (28 variants) and 20 simulated data sets of 60 variants with about 10% outliers.
#   A. The vignette's printed result (effect 0.6177372, SE 0.0256347) through the app path (study pair, script, replay).
#   B. The variance: the linear-response formula of the supplement built in full ((3N + 4)-square matrices) equals MR Studio's exact
#      arrow-shaped solve.
#   C. Agreement with the package. The package stops when the evidence lower bound changes by about 1e-6 (relative); MR Studio iterates
#      to convergence (1e-12). With the package's stopping rule, MR Studio agrees with the package; at full convergence the remaining
#      difference is reported. Where the two reach different optima, the evidence lower bound of each is reported.
# Criteria (set after a first comparison, see README): A and C within 0.01 SE for the effect and 1% for the SE; B relative 1e-8.
#   Rscript --vanilla validation/zhao2020_bwmr/run_validation.R          (part C needs the package BWMR, installed for validation only)
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "zhao2020_bwmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("bwmrval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
eval(parse(text = BWMR_MRS_SRC))   # the function every study script contains
if (!requireNamespace("BWMR", quietly = TRUE)) stop("Part C compares with the package BWMR (GitHub jiazhao97/BWMR, commit c00e818); install it to run this validation.")
e0 <- new.env(); utils::data("ExampleData", package = "BWMR", envir = e0); x <- e0$ExampleData
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)

# A. vignette, through the app
d <- data.frame(variant_id = if ("SNP" %in% names(x)) x$SNP else sprintf("v%02d", seq_len(nrow(x))), ea = "A", oa = "G", beta_x = x$beta.exposure, se_x = x$se.exposure, eaf_x = NA_real_,
                beta_y = x$beta.outcome, se_y = x$se.outcome, eaf_y = NA_real_, p_x = NA_real_, p_y = NA_real_, n_x = NA_real_, n_y = NA_real_, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "pair_iv", data = d, exposure_id = "LDL-c", outcome_id = "dyslipidemia", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("BWMR"); r <- studio_run(h, "LDL-c", "dyslipidemia", inp, "bwmr", idx = 1L, params = list(ci_level = 0.95)); pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl)
add("A", "effect (vignette)", "0.6177372", sprintf("%.7f", e$estimate), "within 0.01 SE", abs(e$estimate - 0.6177372) / 0.0256347 < 0.01)
add("A", "standard error (vignette)", "0.0256347", sprintf("%.7f", e$se), "within 1%", abs(e$se / 0.0256347 - 1) < 0.01)
add("A (reported)", "P value (vignette)", "2.642837e-128", format(signif(e$p, 7)), "reported", TRUE)
add("A", "replay of the study script", "identical", pl$replay$status, "identical", identical(pl$replay$status, "identical"))

# data sets for B and C
sets <- list(package_example = data.frame(bx = x$beta.exposure, sx = x$se.exposure, by = x$beta.outcome, sy = x$se.outcome),
             waterworth_lipids = data.frame(bx = MendelianRandomization::ldlc, sx = MendelianRandomization::ldlcse, by = MendelianRandomization::chdlodds, sy = MendelianRandomization::chdloddsse))
for (s in 1:20) { set.seed(s); n <- 60; sx <- runif(n, 0.01, 0.03); bx <- rnorm(n, 0.08, 0.04) + rnorm(n, 0, sx); sy <- runif(n, 0.01, 0.04)
  o <- rbinom(n, 1, 0.1); by <- 0.3 * bx + rnorm(n, 0, sy) + rnorm(n, 0, 0.005) + o * rnorm(n, 0.08, 0.03); sets[[sprintf("simulated_%02d", s)]] <- data.frame(bx, sx, by, sy) }

# B. dense linear-response formula of the supplement, -[(V H1 - I)^-1 V]_11, at MR Studio's converged solution
dense_var <- function(g, G, sx, sy, f) {
  N <- length(g); K <- 3 * N + 4; v <- sy^2 + f$tau2; mb <- f$beta; sb <- f$se_mfvb^2; pw <- f$weights; a <- f$a; b <- f$b
  s2g <- f$s2_gamma; mg <- f$mu_gamma   # the same variational state as the arrow solve
  V <- H1 <- matrix(0, K, K); V[1:2, 1:2] <- c(sb, 2 * mb * sb, 2 * mb * sb, 2 * sb^2 + 4 * mb^2 * sb)
  for (j in seq_len(N)) { k <- 3 + 3 * (j - 1)
    V[k:(k + 1), k:(k + 1)] <- c(s2g[j], 2 * mg[j] * s2g[j], 2 * mg[j] * s2g[j], 2 * s2g[j]^2 + 4 * mg[j]^2 * s2g[j]); V[k + 2, k + 2] <- pw[j] * (1 - pw[j])
    H1[1, k] <- H1[k, 1] <- pw[j] * G[j] / v[j]; H1[1, k + 2] <- H1[k + 2, 1] <- mg[j] * G[j] / v[j]
    H1[2, k + 1] <- H1[k + 1, 2] <- -pw[j] / (2 * v[j]); H1[2, k + 2] <- H1[k + 2, 2] <- -(mg[j]^2 + s2g[j]) / (2 * v[j])
    H1[k, k + 2] <- H1[k + 2, k] <- mb * G[j] / v[j]; H1[k + 1, k + 2] <- H1[k + 2, k + 1] <- -(mb^2 + sb) / (2 * v[j])
    H1[K - 1, k + 2] <- H1[k + 2, K - 1] <- 1; H1[K, k + 2] <- H1[k + 2, K] <- -1 }
  V[(K - 1):K, (K - 1):K] <- c(trigamma(a) - trigamma(a + b), -trigamma(a + b), -trigamma(a + b), trigamma(b) - trigamma(a + b))
  -(solve(V %*% H1 - diag(K)) %*% V)[1, 1] }
dB <- max(vapply(sets, function(s) { f <- bwmr_mrs(s$bx, s$by, s$sx, s$sy); abs(dense_var(s$bx, s$by, s$sx, s$sy, f) / f$se_beta^2 - 1) }, 0))
add("B", "variance: full linear-response formula vs the arrow solve, largest relative difference (22 data sets)", "0", sprintf("%.1e", dB), "< 1e-8", dB < 1e-8)

# C. against the package
suppressPackageStartupMessages(library(ggplot2))
src <- deparse(BWMR::BWMR); src <- sub('class(Sigma_hat) == "try-error"', 'inherits(Sigma_hat, "try-error")', src, fixed = TRUE)
bw <- eval(parse(text = src)); environment(bw) <- asNamespace("BWMR")
cmp <- do.call(rbind, lapply(names(sets), function(k) { s <- sets[[k]]
  invisible(utils::capture.output(p <- suppressWarnings(bw(s$bx, s$by, s$sx, s$sy))))
  f <- bwmr_mrs(s$bx, s$by, s$sx, s$sy); f6 <- bwmr_mrs(s$bx, s$by, s$sx, s$sy, tol = 1e-6)
  fp <- bwmr_mrs(s$bx, s$by, s$sx, s$sy, init = list(beta = p$beta, weights = p$weights, tau2 = p$tau^2))   # from the package's answer
  same <- abs(fp$beta - f$beta) / f$se_beta < 0.01
  data.frame(data = k, variants = nrow(s), package_beta = p$beta, package_se = p$se_beta, mrstudio_beta = f$beta, mrstudio_se = f$se_beta, iterations = f$iterations,
             same_optimum = same, elbo_mrstudio = f$elbo, elbo_from_package_answer = fp$elbo,
             stop1e6_dbeta_in_se = abs(f6$beta - p$beta) / p$se_beta, stop1e6_rel_se = abs(f6$se_beta / p$se_beta - 1),
             full_dbeta_in_se = abs(f$beta - p$beta) / p$se_beta, full_rel_se = abs(f$se_beta / p$se_beta - 1)) }))
utils::write.csv(cmp, file.path(vdir, "package_comparison.csv"), row.names = FALSE)
so <- cmp[cmp$same_optimum, ]; dif <- cmp[!cmp$same_optimum, ]
add("C", sprintf("package's stopping rule: largest effect difference, %d of %d data sets with the same optimum", nrow(so), nrow(cmp)), "0", sprintf("%.4f SE", max(so$stop1e6_dbeta_in_se)), "< 0.01 SE", max(so$stop1e6_dbeta_in_se) < 0.01)
add("C", "package's stopping rule: largest SE difference, same data sets", "0", sprintf("%.2f%%", 100 * max(so$stop1e6_rel_se)), "< 1%", max(so$stop1e6_rel_se) < 0.01)
add("C (reported)", "full convergence (MR Studio's default): largest effect difference, same data sets", "0", sprintf("%.4f SE (%s)", max(so$full_dbeta_in_se), so$data[which.max(so$full_dbeta_in_se)]), "reported", TRUE)
add("C (reported)", "full convergence: largest SE difference, same data sets", "0", sprintf("%.2f%% (%s)", 100 * max(so$full_rel_se), so$data[which.max(so$full_rel_se)]), "reported", TRUE)
if (nrow(dif)) add("C (reported)", sprintf("different optimum: %s", paste(dif$data, collapse = ", ")), "package's answer",
                   paste(sprintf("%s: MR Studio %.4f (SE %.4f), evidence lower bound %.3f; package %.4f (SE %.4f), lower bound at its optimum %.3f", dif$data, dif$mrstudio_beta, dif$mrstudio_se, dif$elbo_mrstudio, dif$package_beta, dif$package_se, dif$elbo_from_package_answer), collapse = "; "),
                   "MR Studio's optimum has the higher evidence lower bound", all(dif$elbo_mrstudio >= dif$elbo_from_package_answer))
res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
print(res[, c("part", "quantity", "reference", "mr_studio", "pass")], row.names = FALSE); ok <- all(res$pass); cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
