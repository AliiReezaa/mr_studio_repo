# Validation study 13 (Round 63): BiDirMRcML and BiDirCDcML (Xue and Pan 2022, PLoS Genet 18:e1010205, doi:10.1371/journal.pgen.1010205,
# PMID 35576237), GitHub xue-hr/BiDirectCausal (MIT) with MRcML and MRCD, pinned. The paper's applications use genome-wide data for 16 traits,
# which are not bundled; this study checks (1) that MR Studio's across-pairs script equals the package called directly on the same union of
# instruments, for all four versions and both directions, and (2) that it recovers a known two-way effect (simulated traits B: 0.25 forward,
# 0.10 back), with replay. Rscript --vanilla validation/xue2022_bidir_cml/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "xue2022_bidir_cml"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("bdval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R", "mediation_core.R", "bidir_core.R")) source(file.path(args_dir, f)) })
h <- studio_new("BiDir"); a <- pair_input("expB", "outB"); b <- pair_input("outB", "expB")
f1 <- studio_run(h, a$exposure_label, a$outcome_label, a$input, "ivw_mre", idx = 1L)$pair_folder; f2 <- studio_run(h, b$exposure_label, b$outcome_label, b$input, "ivw_mre", idx = 2L)$pair_folder
pr <- list(list(ex = "expB", out = "outB", folder = f1, exName = "B", outName = "outcome B", rev = FALSE), list(ex = "outB", out = "expB", folder = f2, exName = "outcome B", outName = "B", rev = TRUE))
res <- bidir_cml_run(h, pr)[[1]]
# direct call, as the package documents it, on the same union
x <- a$input$data; y <- b$input$data
u <- rbind(data.frame(id = x$variant_id, bX = x$beta_x, sX = x$se_x, bY = x$beta_y, sY = x$se_y, nX = x$n_x, nY = x$n_y), data.frame(id = y$variant_id, bX = y$beta_y, sX = y$se_y, bY = y$beta_x, sY = y$se_x, nX = y$n_y, nY = y$n_x))
u <- u[!duplicated(u$id), ]; seed <- .seed_of(paste("bidircml", f1)); suppressPackageStartupMessages({ library(MRcML); library(MRCD) })
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); d1 <- BiDirectCausal::BiDirMRcML(u$bX, u$bY, u$sX, u$sY, median(u$nX), median(u$nY), random.seed = seed)
set.seed(seed); d2 <- BiDirectCausal::BiDirCDcML(u$bX, u$bY, u$sX, u$sY, median(u$nX), median(u$nY), random.seed = seed)
res$direct <- mapply(function(m, v, dd) (if (m == "BiDirMRcML") d1 else d2)[[sprintf("%s.est.%s", dd, v)]], res$method, res$version, res$direction)
res$diff <- abs(res$estimate - res$direct)
g <- function(dd) res[res$method == "BiDirMRcML" & res$version == "S.DP" & res$direction == dd, ]
ok <- all(res$diff < 1e-12) && g("XtoY")$ci_lower <= 0.25 && 0.25 <= g("XtoY")$ci_upper && g("YtoX")$ci_lower <= 0.10 && 0.10 <= g("YtoX")$ci_upper
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); print(res[, c("method", "version", "direction", "estimate", "se", "direct", "diff")], digits = 4, row.names = FALSE)
cat(sprintf("\nvariants: %d\nVALIDATION %s\n", res$variants[1], if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
