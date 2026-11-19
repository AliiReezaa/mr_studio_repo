# Validation study 9 (Round 59): TwoSampleMR estimators and diagnostics added in MR Studio, on the published Waterworth et al. lipid
# example (28 variants, LDL-c on CHD, as bundled with MendelianRandomization and used by Burgess and Thompson 2015 and Patel et al. 2023).
# Methods: MR-Egger with bootstrap SE, penalised weighted median, simple and weighted mode (NOME), Rucker model selection (point,
# bootstrap, jackknife) and the sign concordance test. Criterion: MR Studio, run as a study pair with replayed scripts, equals TwoSampleMR
# called directly with the same data, settings and seed.
#   Rscript --vanilla validation/twosamplemr_methods/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "twosamplemr_methods"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("tsval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f))
  library(MendelianRandomization); library(TwoSampleMR) })
d <- data.frame(variant_id = sprintf("waterworth_%02d", seq_along(ldlc)), ea = "A", oa = "G", beta_x = ldlc, se_x = ldlcse, eaf_x = NA_real_, beta_y = chdlodds, se_y = chdloddsse, eaf_y = NA_real_,
                p_x = 2 * pnorm(-abs(ldlc / ldlcse)), p_y = 2 * pnorm(-abs(chdlodds / chdloddsse)), n_x = 17723, n_y = 17723, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "pair_iv", data = d, exposure_id = "LDL-c", outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("TwoSampleMR methods"); res <- list(); k <- 0L
seed_of <- function(r, m) as.integer(sub("^seed <- ([0-9]+).*", "\\1", grep("^seed <- ", readLines(Sys.glob(file.path(h$dir, "03_pairs", r$pair_folder, "02_scripts", sprintf("*_%s.R", m)))), value = TRUE)[1]))
tsd <- data.frame(id.exposure = "LDL-c", id.outcome = "CHD", exposure = "LDL-c", outcome = "CHD", SNP = d$variant_id, beta.exposure = ldlc, beta.outcome = chdlodds, se.exposure = ldlcse, se.outcome = chdloddsse, mr_keep = TRUE)
P <- function(...) { p <- default_parameters(); a <- list(...); p[names(a)] <- a; p }
ref <- function(m, seed, extra = list()) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  switch(m, egger_boot = mr_egger_regression_bootstrap(ldlc, chdlodds, ldlcse, chdloddsse, P(nboot = 1000))$b, pwm = mr_penalised_weighted_median(ldlc, chdlodds, ldlcse, chdloddsse, P(nboot = 1000))$b,
         smode_nome = mr_simple_mode_nome(ldlc, chdlodds, ldlcse, chdloddsse, P(nboot = 1000))$b, wmode_nome = mr_weighted_mode_nome(ldlc, chdlodds, ldlcse, chdloddsse, P(nboot = 1000))$b,
         rucker = { r1 <- suppressMessages(mr_rucker(tsd, P(nboot = 1000)))[[1]]; if (extra$v == "point") r1$selected$Estimate[1] else {
           rr <- if (extra$v == "bootstrap") suppressMessages(mr_rucker_bootstrap(tsd, P(nboot = 1000)))$res else suppressMessages(mr_rucker_jackknife(tsd, P(nboot = 1000)))[[1]]$res
           rr$Estimate[grepl("Rucker median", rr$Method)][1] } },
         sign = mr_sign(ldlc, chdlodds, ldlcse, chdloddsse)$b) }
runs <- list(list(m = c("egger_boot", "pwm", "smode_nome", "wmode_nome", "sign"), p = list()), list(m = "rucker", p = list(rk_variant = "point")),
             list(m = "rucker", p = list(rk_variant = "bootstrap")), list(m = "rucker", p = list(rk_variant = "jackknife")))
for (ru in runs) { k <- k + 1L; r <- studio_run(h, "LDL-c", "CHD", inp, ru$m, idx = k, params = c(list(ci_level = 0.95), ru$p)); pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl)
  for (m in ru$m) { s <- seed_of(r, m); got <- if (m == "sign") Filter(function(g) g$name == "sign_concordance", pl$diagnostics)[[1]]$value else e$estimate[e$method == m & e$estimand == "causal_effect"]
    want <- ref(m, s, list(v = ru$p$rk_variant %||% "point")); res[[length(res) + 1L]] <- data.frame(method = if (m == "rucker") paste0("rucker (", ru$p$rk_variant, ")") else m, mr_studio = got, twosamplemr = want,
      diff = abs(got - want), replay = pl$replay$status, stringsAsFactors = FALSE) } }
res <- do.call(rbind, res); ok <- all(res$diff < 1e-10) && all(res$replay == "identical")
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); print(res, digits = 6, row.names = FALSE)
cat(sprintf("\nVALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
