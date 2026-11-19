# Validation study 8 (Round 58): the MendelianRandomization estimators added in MR Studio, on the worked example of the package paper
# (Patel et al. 2023, Wellcome Open Res 8:449, doi:10.12688/wellcomeopenres.19995.2, PMID 37915953): 28 Waterworth et al. variants,
# LDL-c (and HDL-c, triglycerides) on CHD, n = 17,723, as bundled with MendelianRandomization 0.10.0. The paper prints the calls
# (for example mr_cML(mr_input(ldlc, ldlcse, chdlodds, chdloddsse), n = 17723)); its printed MVMR-IVW table is reproduced too.
# Criterion: MR Studio, run as a study pair with a replayed script, equals the package called as the paper shows (same seed).
#   Rscript --vanilla validation/patel2023_mr_methods/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "patel2023_mr_methods"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mrval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir)); for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)); library(MendelianRandomization) })
N <- 17723
d <- data.frame(variant_id = sprintf("waterworth_%02d", seq_along(ldlc)), ea = "A", oa = "G", beta_x = ldlc, se_x = ldlcse, eaf_x = NA_real_, beta_y = chdlodds, se_y = chdloddsse, eaf_y = NA_real_,
                p_x = 2 * pnorm(-abs(ldlc / ldlcse)), p_y = 2 * pnorm(-abs(chdlodds / chdloddsse)), n_x = N, n_y = N, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "pair_iv", data = d, exposure_id = "LDL-c", outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("Patel 2023 methods"); ms <- c("maxlik", "mbe", "cml", "robust_ivw", "pivw")
r <- studio_run(h, "LDL-c", "CHD", inp, ms, idx = 1L, params = list(ci_level = 0.95))
pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl); seed <- as.integer(sub("^seed <- ([0-9]+).*", "\\1", grep("^seed <- ", readLines(Sys.glob(file.path(h$dir, "03_pairs", r$pair_folder, "02_scripts", "*_cml.R"))), value = TRUE)[1]))
mi <- mr_input(ldlc, ldlcse, chdlodds, chdloddsse); sd <- function(f) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); f }
ref <- list(maxlik = mr_maxlik(mi), mbe = sd(mr_mbe(mi, seed = seed)), cml = sd(mr_cML(mi, n = N, random_seed = seed)), robust_ivw = mr_ivw(mi, robust = TRUE, penalized = TRUE), pivw = sd(mr_pivw(mi)))
res <- do.call(rbind, lapply(ms, function(m) { x <- e[e$method == m, ]; f <- ref[[m]]
  data.frame(method = m, mr_studio = x$estimate, package = f@Estimate, se_studio = x$se, se_package = f@StdError, diff = max(abs(x$estimate - f@Estimate), abs(x$se - f@StdError)), stringsAsFactors = FALSE) }))
# MVMR: MVMR-cML as the paper calls it, and the paper's printed MVMR-IVW table (estimates 1.925, -0.590, 0.723; SE 0.439, 0.555, 0.230).
dm <- data.frame(variant_id = d$variant_id, ea = "A", oa = "G", beta_x_1 = ldlc, se_x_1 = ldlcse, n_x_1 = N, beta_x_2 = hdlc, se_x_2 = hdlcse, n_x_2 = N, beta_x_3 = trig, se_x_3 = trigse, n_x_3 = N,
                 beta_y = chdlodds, se_y = chdloddsse, n_y = N, mr_keep = TRUE, stringsAsFactors = FALSE)
im <- list(kind = "multi_iv", data = dm, exposure_id = "LDL-c|HDL-c|TG", outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta")
r2 <- studio_run(h, "Lipids", "CHD", im, c("mvivw", "mvcml"), idx = 2L, params = list(ci_level = 0.95, mv_model = "default")); pl2 <- pair_payload(h$dir, r2$pair_folder); e2 <- .est_df(pl2)
seed2 <- as.integer(sub("^seed <- ([0-9]+).*", "\\1", grep("^seed <- ", readLines(Sys.glob(file.path(h$dir, "03_pairs", r2$pair_folder, "02_scripts", "*_mvcml.R"))), value = TRUE)[1]))
mvi <- mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse), by = chdlodds, byse = chdloddsse)
fc <- { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed2); mr_mvcML(mvi, n = N, seed = seed2) }
c2 <- e2[e2$method == "mvcml" & e2$estimand == "direct_effect", ]; iv <- e2[e2$method == "mvivw" & e2$estimand == "direct_effect", ]
res <- rbind(res, data.frame(method = sprintf("mvcml[%d]", 1:3), mr_studio = c2$estimate, package = fc@Estimate, se_studio = c2$se, se_package = fc@StdError,
                             diff = pmax(abs(c2$estimate - fc@Estimate), abs(c2$se - fc@StdError)), stringsAsFactors = FALSE),
             data.frame(method = sprintf("mvivw[%d] vs paper", 1:3), mr_studio = round(iv$estimate, 3), package = c(1.925, -0.590, 0.723), se_studio = round(iv$se, 3), se_package = c(0.439, 0.555, 0.230),
                        diff = pmax(abs(round(iv$estimate, 3) - c(1.925, -0.590, 0.723)), abs(round(iv$se, 3) - c(0.439, 0.555, 0.230))), stringsAsFactors = FALSE))
ok <- all(res$diff < 1e-8) && identical(pl$replay$status, "identical") && identical(pl2$replay$status, "identical")
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); print(res, digits = 6, row.names = FALSE)
cat(sprintf("\nreplay: %s, %s\nVALIDATION %s\n", pl$replay$status, pl2$replay$status, if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
