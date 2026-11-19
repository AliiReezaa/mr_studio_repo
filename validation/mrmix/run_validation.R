# Validation study 10 (Round 60): MRMix (Qi and Chatterjee 2019, Nat Commun 10:1941, doi:10.1038/s41467-019-09432-2), on real GWAS data:
# the education instruments of validation study 7 (Lee et al. 2018, ieu-a-1239) with BMI (Yengo et al. 2018, ieu-b-40, continuous) and with
# CHD (CARDIoGRAMplusC4D, ieu-a-7, binary), which carry allele frequencies and sample sizes. It reads the OpenGWAS data saved by study 7
# (MRSTUDIO_VALIDATION_CACHE) and sends no query. Criteria: (1) MR Studio's standardised estimate equals TwoSampleMR::run_mrmix, the
# package's reference wrapper, on the same harmonised data; (2) the reported estimate is that value divided by the median implied SD of
# the exposure (times that of the outcome for a continuous outcome); (3) replay identical.
#   MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/mrmix/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "mrmix"))
cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE"); if (!nzchar(cache) || !length(Sys.glob(file.path(cache, "carter2019_*.rds")))) { cat("Skipped: needs the OpenGWAS data saved by validation/carter2019_mediation (MRSTUDIO_VALIDATION_CACHE).\n"); quit(status = 0) }
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mixval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
ie <- readRDS(Sys.glob(file.path(cache, "carter2019_ins_ieu-a-1239_r2_0.001.rds")))
h <- studio_new("MRMix"); out <- list(); k <- 0L
for (oc in c("ieu-b-40", "ieu-a-7")) { k <- k + 1L
  o <- readRDS(Sys.glob(file.path(cache, sprintf("carter2019_out_%s_*.rds", oc)))[1]); H <- suppressMessages(TwoSampleMR::harmonise_data(ie, o, action = 2)); H <- H[H$mr_keep, ]
  bin <- oc == "ieu-a-7"; H$units.exposure <- "SD"; H$units.outcome <- if (bin) "log odds" else "SD"
  d <- data.frame(variant_id = H$SNP, ea = H$effect_allele.exposure, oa = H$other_allele.exposure, beta_x = H$beta.exposure, se_x = H$se.exposure, eaf_x = H$eaf.exposure, beta_y = H$beta.outcome, se_y = H$se.outcome,
                  eaf_y = H$eaf.outcome, p_x = H$pval.exposure, p_y = H$pval.outcome, n_x = H$samplesize.exposure, n_y = H$samplesize.outcome, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  inp <- list(kind = "pair_iv", data = d, exposure_id = "education", outcome_id = oc, outcome_scale = if (bin) "log_or" else "beta", exposure_scale = "beta")
  r <- studio_run(h, "Education", oc, inp, "mrmix", idx = k, params = list(ci_level = 0.95)); pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl); dg <- function(n) Filter(function(g) g$name == n, pl$diagnostics)[[1]]$value
  # TwoSampleMR::run_mrmix treats "SD" units as already standardised; with units unset it standardises like MR Studio. Use the generic unit.
  H2 <- H; H2$units.exposure <- "units"; H2$units.outcome <- if (bin) "log odds" else "units"; ref <- suppressWarnings(suppressMessages(TwoSampleMR::run_mrmix(H2)))[[1]]
  maf <- ifelse(is.finite(d$eaf_x), d$eaf_x, d$eaf_y); keep <- is.finite(maf) & maf > 0 & maf < 1; g <- sqrt(2 * maf[keep] * (1 - maf[keep]))
  k_exp <- (if (bin) 1 else stats::median(d$se_y[keep] * sqrt(d$n_y[keep]) * g)) / stats::median(d$se_x[keep] * sqrt(d$n_x[keep]) * g)
  out[[k]] <- data.frame(outcome = if (bin) "CHD (binary)" else "BMI (continuous)", n_variants = sum(keep), theta_std_studio = dg("mrmix_theta_std"), theta_std_twosamplemr = ref$theta,
                         estimate = e$estimate, expected = ref$theta * k_exp, se = e$se, replay = pl$replay$status, stringsAsFactors = FALSE) }
res <- do.call(rbind, out); res$diff_std <- abs(res$theta_std_studio - res$theta_std_twosamplemr); res$diff_scale <- abs(res$estimate - res$expected)
ok <- all(res$diff_std < 1e-10) && all(res$diff_scale < 1e-10) && all(res$replay == "identical")
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); print(res, digits = 5, row.names = FALSE); cat(sprintf("\nVALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
