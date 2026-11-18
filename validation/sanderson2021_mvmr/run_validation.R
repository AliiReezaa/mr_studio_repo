# Validation study 6: the MVMR package tutorial of Sanderson, Spiller and Bowden 2021 (Stat Med 40:5434-5452, doi:10.1002/sim.9133,
# PMID 34338327). Data: rawdat_mvmr bundled with the MVMR package (145 variants; LDL-c, HDL-c and triglycerides from GLGC on systolic
# blood pressure from UK Biobank). Reference values are computed with the MVMR package itself (strength_mvmr, pleiotropy_mvmr, ivw_mvmr),
# the authors' own implementation. Run from the repository root (needs the MVMR package, not part of the app runtime):
#   Rscript --vanilla validation/sanderson2021_mvmr/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "sanderson2021_mvmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mvval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
if (!requireNamespace("MVMR", quietly = TRUE)) stop("install the MVMR package (remotes::install_github(\"WSpiller/MVMR\")) to run this validation")
suppressMessages(library(MVMR)); raw <- MVMR::rawdat_mvmr; ex <- c("LDL-c", "HDL-c", "Triglycerides")

# The authors' implementation (tutorial steps 2 to 5, zero covariance as in the tutorial's first call).
F <- MVMR::format_mvmr(BXGs = raw[, 1:3], BYG = raw[, 7], seBXGs = raw[, 4:6], seBYG = raw[, 8], RSID = raw[, 9])
invisible(utils::capture.output(suppressWarnings({ ref_f <- unlist(MVMR::strength_mvmr(F, gencov = 0)); ref_q <- MVMR::pleiotropy_mvmr(F, gencov = 0) }), ref_iv <- MVMR::ivw_mvmr(F)))
ref_iv <- ref_iv[, 1:2]

# MR Studio, run as a study pair on the same data (exposure-outcome effect alleles are already aligned in rawdat_mvmr).
d <- data.frame(variant_id = raw$SNP, ea = "A", oa = "G", beta_x_1 = raw$LDL_beta, se_x_1 = raw$LDL_se, beta_x_2 = raw$HDL_beta, se_x_2 = raw$HDL_se,
                beta_x_3 = raw$Trg_beta, se_x_3 = raw$Trg_se, beta_y = raw$SBP_beta, se_y = raw$SBP_se, mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "multi_iv", data = d, exposure_id = paste(ex, collapse = "|"), outcome_id = "SBP", outcome_scale = "beta", exposure_scale = "beta")
h <- studio_new("Sanderson 2021 MVMR tutorial")
r <- studio_run(h, "Lipids", "SBP", inp, c("mvivw", "mvegger", "mvmedian", "mvlasso"), idx = 1L, params = list(ci_level = 0.95, mv_model = "default", mv_iter = 1000L))
pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl); iv <- e[e$method == "mvivw" & e$estimand == "direct_effect", , drop = FALSE]
dg <- Filter(function(g) g$method == "mvivw", pl$diagnostics); val <- function(n) vapply(Filter(function(g) g$name == n, dg), function(g) as.numeric(g$value), 0)
cf <- val("mv_cond_f"); qa <- Filter(function(g) g$name == "mv_q_a", dg)[[1]]

res <- data.frame(quantity = c(sprintf("conditional F, %s", ex), "Q_A", "Q_A df", sprintf("MVMR-IVW estimate, %s", ex), sprintf("MVMR-IVW SE, %s", ex)),
                  mvmr_package = c(as.numeric(ref_f), ref_q$Qstat, nrow(raw) - 3, ref_iv[, 1], ref_iv[, 2]),
                  mr_studio = c(cf, qa$value, qa$df, iv$estimate, iv$se), stringsAsFactors = FALSE)
res$covariance <- "zero"
# Part 2 (Round 53): the tutorial's phenotypic correlation matrix (LDL-HDL -0.1, LDL-TG -0.05, HDL-TG 0.2), as phenocov_mvmr builds the
# covariances from it, against strength_mvmr and pleiotropy_mvmr with gencov = Xcovmat. MR Studio takes the same three correlations as settings.
pcor <- matrix(c(1, -0.1, -0.05, -0.1, 1, 0.2, -0.05, 0.2, 1), 3, 3); Xcov <- MVMR::phenocov_mvmr(pcor, F[, 7:9])
invisible(utils::capture.output({ ref_f2 <- unlist(MVMR::strength_mvmr(F, gencov = Xcov)); ref_q2 <- MVMR::pleiotropy_mvmr(F, gencov = Xcov) }))
r2 <- studio_run(h, "Lipids, correlated", "SBP", inp, "mvivw", idx = 2L, params = list(ci_level = 0.95, mv_model = "default", mv_r_1_2 = -0.1, mv_r_1_3 = -0.05, mv_r_2_3 = 0.2))
pl2 <- pair_payload(h$dir, r2$pair_folder); dg <- Filter(function(g) g$method == "mvivw", pl2$diagnostics)
cf2 <- val("mv_cond_f"); qa2 <- Filter(function(g) g$name == "mv_q_a", dg)[[1]]
res2 <- data.frame(quantity = c(sprintf("conditional F, %s", ex), "Q_A"), mvmr_package = c(as.numeric(ref_f2), ref_q2$Qstat), mr_studio = c(cf2, qa2$value), covariance = "tutorial correlation matrix", stringsAsFactors = FALSE)
blocked <- ADAPTER_MVIVW_V11$check(inp, list(mv_r_1_2 = 0.95, mv_r_1_3 = -0.95, mv_r_2_3 = 0.95))$status == "block"   # not a valid correlation matrix
res <- rbind(res, res2)
res$abs_diff <- abs(res$mr_studio - res$mvmr_package); res$rel_diff <- res$abs_diff / abs(res$mvmr_package)
# Criterion: the same formulas on the same data, so agreement to numerical precision (relative difference below 1e-8). The MVMR-IVW
# standard errors are compared with the package's ivw_mvmr, which uses the multiplicative random-effects (residual-scaled) SEs, the
# same as mr_mvivw's default model for this heterogeneous fit.
ok <- all(res$rel_diff < 1e-8) && blocked && identical(pl$replay$status, "identical") && identical(pl2$replay$status, "identical") && all(vapply(r$results, function(x) identical(x$status, "ok"), TRUE))
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
print(format(res, digits = 10), row.names = FALSE); cat(sprintf("\ninvalid correlation matrix blocked: %s", blocked)); cat(sprintf("\nreplay: %s\nVALIDATION %s\n", pl$replay$status, if (ok) "PASSED" else "FAILED"))
if (!ok) quit(status = 1)
