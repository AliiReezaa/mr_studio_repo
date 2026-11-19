# Validation study 12 (Round 62): BWMR (Zhao et al. 2020, Bioinformatics 36:1501-1508, doi:10.1093/bioinformatics/btz749, PMID 31593215).
# The package vignette (BWMR_package.pdf, GitHub jiazhao97/BWMR at the pinned commit) fits its ExampleData (LDL cholesterol on dyslipidemia,
# from the paper's real-data analysis, bundled with the package) and prints: effect 0.6177372, standard error 0.0256347, p-value 2.642837e-128.
# Criterion: MR Studio, run as a study pair with a replayed script, reproduces the printed values to the printed digits.
#   Rscript --vanilla validation/zhao2020_bwmr/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "zhao2020_bwmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("bwmrval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
e0 <- new.env(); utils::data("ExampleData", package = "BWMR", envir = e0); x <- e0$ExampleData
d <- data.frame(variant_id = if ("SNP" %in% names(x)) x$SNP else sprintf("v%02d", seq_len(nrow(x))), ea = "A", oa = "G", beta_x = x$beta.exposure, se_x = x$se.exposure, eaf_x = NA_real_,
                beta_y = x$beta.outcome, se_y = x$se.outcome, eaf_y = NA_real_, p_x = NA_real_, p_y = NA_real_, n_x = NA_real_, n_y = NA_real_, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "pair_iv", data = d, exposure_id = "LDL-c", outcome_id = "dyslipidemia", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("BWMR"); r <- studio_run(h, "LDL-c", "dyslipidemia", inp, "bwmr", idx = 1L, params = list(ci_level = 0.95)); pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl)
res <- data.frame(quantity = c("effect", "standard error", "p-value"), vignette = c(0.6177372, 0.0256347, 2.642837e-128), mr_studio = c(e$estimate, e$se, e$p))
res$match <- c(abs(signif(res$mr_studio[1:2], 7) - res$vignette[1:2]) < 1e-12, abs(signif(res$mr_studio[3], 7) / res$vignette[3] - 1) < 1e-9)
ok <- all(res$match) && identical(pl$replay$status, "identical"); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
print(res, digits = 8, row.names = FALSE); cat(sprintf("variants: %d; replay: %s\nVALIDATION %s\n", nrow(d), pl$replay$status, if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
