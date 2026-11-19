# Validation study 11 (Round 61): MR-GRIP, Dudbridge et al. 2025, PLoS Genet 21:e1011967 (doi:10.1371/journal.pgen.1011967, PMID 41468517),
# Table 6: BMI on coronary artery disease, 97 GIANT BMI variants (96 with CARDIoGRAMplusC4D effects). Published: IVW OR 1.44 (1.26 to 1.63),
# MR-GRIP 1.54 (1.25 to 1.89), MR-GRIP adjusted for weak instruments 1.56 (1.07 to 2.26). Data from OpenGWAS: BMI, Locke et al. 2015 GIANT
# (ieu-a-2); CAD, CARDIoGRAMplusC4D (ieu-a-7). Opt-in, with the optional cache (no repeated queries):
#   MRSTUDIO_VALIDATE_OPENGWAS=1 MRSTUDIO_VALIDATION_CACHE=<folder> OPENGWAS_JWT=<token> Rscript --vanilla validation/dudbridge2025_grip/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "dudbridge2025_grip"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("gripval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
# Part A (offline): MR Studio against TwoSampleMR::mr_grip on the Waterworth lipid example (LDL-c on CHD, 28 variants), both versions.
suppressMessages(library(MendelianRandomization))
dw <- data.frame(variant_id = sprintf("w%02d", seq_along(ldlc)), ea = "A", oa = "G", beta_x = ldlc, se_x = ldlcse, eaf_x = NA_real_, beta_y = chdlodds, se_y = chdloddsse, eaf_y = NA_real_,
                 p_x = NA_real_, p_y = NA_real_, n_x = 17723, n_y = 17723, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
iw <- list(kind = "pair_iv", data = dw, exposure_id = "LDL-c", outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta"); ha <- studio_new("GRIP A")
ga <- TwoSampleMR::mr_grip(ldlc, chdlodds, ldlcse, chdloddsse, TwoSampleMR::default_parameters())
A <- t(vapply(c("standard", "weak"), function(v) { r <- studio_run(ha, "LDL-c", "CHD", iw, "grip", idx = if (v == "standard") 1L else 2L, params = list(ci_level = 0.95, grip_variant = v)); pl <- pair_payload(ha$dir, r$pair_folder); e <- .est_df(pl)
  b <- e$estimate[e$estimand == "causal_effect"]; i <- e$estimate[e$estimand == "pleiotropy_intercept"]; c(studio = b, package = if (v == "standard") ga$b else ga$b.adj, diff = abs(b - (if (v == "standard") ga$b else ga$b.adj)) + abs(i - ga$b_i), replay = identical(pl$replay$status, "identical")) }, numeric(4)))
print(A); okA <- all(A[, "diff"] < 1e-10) && all(A[, "replay"] == 1); utils::write.csv(data.frame(version = rownames(A), A), file.path(vdir, "results_partA.csv"), row.names = FALSE)
cat(sprintf("Part A (exact against TwoSampleMR): %s\n", if (okA) "PASSED" else "FAILED")); if (!okA) quit(status = 1)
# Part C (offline, exact data): the paper's Tables 5 and 6 from its own input data, S1 Data Tables A (urate on CHD, 31 variants) and B
# (BMI on CAD, 96 variants), saved here as s1_table_a.csv and s1_table_b.csv (PLoS Genet, CC BY 4.0, Dudbridge et al. 2025).
# Criterion: odds ratios and limits equal to the 2 decimals printed, standard errors to the 3 significant digits printed.
paperC <- list(table_a = data.frame(method = c("IVW", "MR-GRIP", "MR-GRIP weak"), or = c(1.11, 1.04, 1.03), lo = c(1.03, 0.96, 0.95), hi = c(1.20, 1.13, 1.12), se = c(0.0401, 0.0417, 0.0421)),
               table_b = data.frame(method = c("IVW", "MR-GRIP", "MR-GRIP weak"), or = c(1.44, 1.54, 1.56), lo = c(1.26, 1.25, 1.07), hi = c(1.63, 1.89, 2.26), se = c(0.0648, 0.105, 0.190)))
hc <- studio_new("GRIP C"); kk <- 0L; C <- list()
for (tb in names(paperC)) { x <- utils::read.csv(file.path(vdir, sprintf("s1_%s.csv", tb)), stringsAsFactors = FALSE)
  dd <- data.frame(variant_id = x$SNP, ea = "A", oa = "G", beta_x = x$beta.exposure, se_x = x$se.exposure, eaf_x = NA_real_, beta_y = x$beta.outcome, se_y = x$se.outcome, eaf_y = NA_real_,
                   p_x = NA_real_, p_y = NA_real_, n_x = NA_real_, n_y = NA_real_, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  ii <- list(kind = "pair_iv", data = dd, exposure_id = tb, outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta")
  kk <- kk + 1L; ra <- studio_run(hc, tb, "CHD", ii, c("ivw_mre", "grip"), idx = kk, params = list(ci_level = 0.95)); kk <- kk + 1L; rb <- studio_run(hc, tb, "CHD", ii, "grip", idx = kk, params = list(ci_level = 0.95, grip_variant = "weak"))
  ea <- .est_df(pair_payload(hc$dir, ra$pair_folder)); eb <- .est_df(pair_payload(hc$dir, rb$pair_folder)); pk <- function(e, m) e[e$method == m & e$estimand == "causal_effect", ]
  got <- rbind(pk(ea, "ivw_mre"), pk(ea, "grip"), pk(eb, "grip")); pp <- paperC[[tb]]
  C[[tb]] <- data.frame(table = if (tb == "table_a") "5 (urate on CHD)" else "6 (BMI on CAD)", method = pp$method, paper = sprintf("%.2f (%.2f to %.2f), %s", pp$or, pp$lo, pp$hi, format(pp$se)),
    mr_studio = sprintf("%.2f (%.2f to %.2f), %s", exp(got$estimate), exp(got$ci_lower), exp(got$ci_upper), signif(got$se, 3)),
    match = abs(round(exp(got$estimate), 2) - pp$or) < 1e-9 & abs(round(exp(got$ci_lower), 2) - pp$lo) < 1e-9 & abs(round(exp(got$ci_upper), 2) - pp$hi) < 1e-9 & abs(signif(got$se, 3) - pp$se) < 1e-9,
    replay = c(pair_payload(hc$dir, ra$pair_folder)$replay$status, pair_payload(hc$dir, ra$pair_folder)$replay$status, pair_payload(hc$dir, rb$pair_folder)$replay$status), stringsAsFactors = FALSE) }
C <- do.call(rbind, C); print(C, row.names = FALSE); utils::write.csv(C, file.path(vdir, "results_partC.csv"), row.names = FALSE)
okC <- all(C$match) && all(C$replay == "identical"); cat(sprintf("Part C (paper Tables 5 and 6 from S1 Data): %s\n", if (okC) "PASSED" else "FAILED"))
# Part B: the paper's Table 6 from OpenGWAS (a diagnostic of the public instrument set; Part C is the reproduction).
jwt <- Sys.getenv("OPENGWAS_JWT"); cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE")
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_OPENGWAS"), "1") || !nzchar(jwt)) { cat("Part B skipped: set MRSTUDIO_VALIDATE_OPENGWAS=1 and OPENGWAS_JWT.\n"); cat(sprintf("VALIDATION %s\n", if (okC) "PASSED" else "FAILED")); quit(status = if (okC) 0 else 1) }
options(ieugwasr_api = "https://api.opengwas.io/api/")
cached <- function(key, expr) { if (!nzchar(cache)) return(expr); f <- file.path(cache, paste0("grip_", key, ".rds")); if (file.exists(f)) return(readRDS(f)); dir.create(cache, showWarnings = FALSE, recursive = TRUE); v <- expr; saveRDS(v, f); v }
ins <- cached("ins_ieu-a-2", suppressMessages(TwoSampleMR::extract_instruments("ieu-a-2"))); o <- cached("out_ieu-a-7", suppressMessages(TwoSampleMR::extract_outcome_data(ins$SNP, "ieu-a-7", proxies = TRUE)))
H <- suppressMessages(TwoSampleMR::harmonise_data(ins, o, action = 2)); H <- H[H$mr_keep, ]
d <- data.frame(variant_id = H$SNP, ea = H$effect_allele.exposure, oa = H$other_allele.exposure, beta_x = H$beta.exposure, se_x = H$se.exposure, eaf_x = H$eaf.exposure, beta_y = H$beta.outcome, se_y = H$se.outcome,
                eaf_y = H$eaf.outcome, p_x = H$pval.exposure, p_y = H$pval.outcome, n_x = H$samplesize.exposure, n_y = H$samplesize.outcome, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "pair_iv", data = d, exposure_id = "BMI", outcome_id = "CAD", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("GRIP"); r1 <- studio_run(h, "BMI", "CAD", inp, c("ivw_mre", "grip"), idx = 1L, params = list(ci_level = 0.95, ivw_model = "mre"))
r2 <- studio_run(h, "BMI", "CAD", inp, "grip", idx = 2L, params = list(ci_level = 0.95, grip_variant = "weak"))
e1 <- .est_df(pair_payload(h$dir, r1$pair_folder)); e2 <- .est_df(pair_payload(h$dir, r2$pair_folder))
g <- TwoSampleMR::mr_grip(d$beta_x, d$beta_y, d$se_x, d$se_y, TwoSampleMR::default_parameters()); iv <- TwoSampleMR::mr_ivw_mre(d$beta_x, d$beta_y, d$se_x, d$se_y)
pick <- function(e, m) e[e$method == m & e$estimand == "causal_effect", ]
res <- data.frame(method = c("IVW", "MR-GRIP", "MR-GRIP weak"), paper_or = c(1.44, 1.54, 1.56), paper_lo = c(1.26, 1.25, 1.07), paper_hi = c(1.63, 1.89, 2.26), paper_se = c(0.0648, 0.105, 0.190),
                  or = exp(c(pick(e1, "ivw_mre")$estimate, pick(e1, "grip")$estimate, pick(e2, "grip")$estimate)), se = c(pick(e1, "ivw_mre")$se, pick(e1, "grip")$se, pick(e2, "grip")$se),
                  package_b = c(iv$b, g$b, g$b.adj), studio_b = c(pick(e1, "ivw_mre")$estimate, pick(e1, "grip")$estimate, pick(e2, "grip")$estimate))
res$within_paper_ci <- res$or >= res$paper_lo & res$or <= res$paper_hi; res$diff_package <- abs(res$studio_b - res$package_b)
reps <- c(pair_payload(h$dir, r1$pair_folder)$replay$status, pair_payload(h$dir, r2$pair_folder)$replay$status)
ok <- okC && all(res$diff_package < 1e-10) && all(res$within_paper_ci) && all(reps == "identical")
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); print(res[, c("method", "paper_or", "paper_lo", "paper_hi", "or", "se", "paper_se", "within_paper_ci", "diff_package")], digits = 3, row.names = FALSE)
cat(sprintf("\nvariants: %d (paper 96); replay: %s\nVALIDATION %s\n", nrow(d), paste(reps, collapse = ", "), if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
