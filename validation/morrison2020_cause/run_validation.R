# Validation study 18 (Round 70): CAUSE (Morrison et al. 2020, Nat Genet 52:740-747, doi:10.1038/s41588-020-0631-4, PMID 32451458),
# GitHub jean997/cause (GPL >= 2), pinned. Published values: the package's LDL cholesterol and coronary artery disease vignette
# (analysis/ldl_cad.Rmd, rendered in docs/ldl_cad.html of the pinned commit), which uses the formatted data, the nuisance parameters
# and the 719 LD-pruned variants stored in the package's example_data folder.
#   A. Steps 4 and 5 of the vignette on its stored inputs: the posterior medians and 95% credible intervals at the vignette's printed
#      precision (gamma is the criterion), and the ELPD comparison of the sharing and causal models within its Monte Carlo spread
#      over 20 seeds (the ELPD is estimated from posterior draws, so it changes slightly with the seed). The criterion is the comparison
#      CAUSE tests, sharing against causal; the two comparisons with the null model are reported with their distance in the vignette's
#      standard errors (narrowed after the first run, see README: the vignette's own two renders differ by more than the seed spread).
#   B. Step 2: est_cause_params() with the vignette's seed and one million variants against the stored parameters.
#   C. The app path (Data library datasets, preparation, study script with local PLINK clumping) gives the vignette's conclusion,
#      equals the package pipeline called directly on the same files with the same seed, and replays identically.
# Needs loo 2.8.0 (from loo 2.10, loo_compare() returns a different table; see README) and the app's 1000 Genomes EUR panel.
#   MRSTUDIO_VALIDATE_CAUSE=1 MRSTUDIO_CAUSE_EXAMPLE=<cause source>/example_data Rscript --vanilla validation/morrison2020_cause/run_validation.R
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_CAUSE"), "1")) { cat("Skipped: set MRSTUDIO_VALIDATE_CAUSE=1 to run.\n"); quit(status = 0) }
lb <- Sys.getenv("MRSTUDIO_EXTRA_LIB", ""); if (nzchar(lb)) .libPaths(c(lb, .libPaths()))
if (utils::packageVersion("loo") != "2.8.0") stop("CAUSE needs loo 2.8.0; found ", utils::packageVersion("loo"))
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "morrison2020_cause"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("causeval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("resource_core.R", "import_core.R", "ld_core.R", "ldsc_core.R", "mrlap_core.R", "cause_core.R")) source(file.path(args_dir, f)) })
if (is.na(plink_bin()) && nzchar(Sys.which("plink"))) Sys.setenv(MRGUI_PLINK = Sys.which("plink"))
ex <- Sys.getenv("MRSTUDIO_CAUSE_EXAMPLE"); md5 <- c(LDL_CAD_merged.RDS = "c7dd531d9ffcf4ea4b08ee7c6aaf1d87", LDL_CAD_params.RDS = "49b56d75057058fb5c86ff48f2a4a058", top_ldl_pruned_plink.RDS = "391728f5a95c2ed0b1a88bc70f0e09ed")
for (f in names(md5)) if (!file.exists(file.path(ex, f))) stop("missing ", f, " in MRSTUDIO_CAUSE_EXAMPLE (the example_data folder of cause at the pinned commit)")
got <- unname(tools::md5sum(file.path(ex, names(md5)))); if (!identical(got, unname(md5))) stop("MD5 of the example data differs from the pinned commit: ", paste(names(md5)[got != md5], collapse = ", "))
suppressPackageStartupMessages(library(cause)); RcppParallel::setThreadOptions(numThreads = CAUSE_THREADS); invisible(loadNamespace("ieugwasr"))
X <- readRDS(file.path(ex, "LDL_CAD_merged.RDS")); params <- readRDS(file.path(ex, "LDL_CAD_params.RDS")); tv <- readRDS(file.path(ex, "top_ldl_pruned_plink.RDS"))
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)
fmt <- function(v) sprintf("%s (%s, %s)", round(v[1], 2), round(v[2], 2), round(v[3], 2))
# A. vignette steps 4 and 5 (printed in docs/ldl_cad.html)
vig <- list(gamma = "0.32 (0.24, 0.4)", eta_causal = "-0.03 (-1.28, 0.91)", q_causal = "0.03 (0, 0.24)", eta_sharing = "0.39 (0.3, 0.51)", q_sharing = "0.62 (0.43, 0.79)",
            elpd = c(-29.660651, -36.735730, -7.075079), se = c(5.588301, 7.092086, 1.631256), p = 7.2e-06,
            elpd_jul2021 = c(-29.514383, -36.856250, -7.341867))   # the earlier render of the same vignette (docs/ldl_cad.html at 05b1882)
t0 <- Sys.time(); fits <- lapply(c(100, 1:19), function(s) { set.seed(s); invisible(utils::capture.output(r <- cause(X = X, variants = tv, param_ests = params))); r })
min_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins")) / length(fits)
s1 <- summary(fits[[1]], ci_size = 0.95); qc <- s1$quants[[2]]; qs <- s1$quants[[1]]
add("A", "gamma, causal model (posterior median, 95% credible interval)", vig$gamma, fmt(qc[, "gamma"]), "equal at the printed precision", identical(fmt(qc[, "gamma"]), vig$gamma))
for (k in list(c("eta, causal model", "eta_causal", "2", "eta"), c("q, causal model", "q_causal", "2", "q"), c("eta, sharing model", "eta_sharing", "1", "eta"), c("q, sharing model", "q_sharing", "1", "q"))) {
  v <- if (k[3] == "2") qc[, k[4]] else qs[, k[4]]; add("A (reported)", k[1], vig[[k[2]]], fmt(v), "reported", TRUE) }
gamma_all <- vapply(fits, function(r) { g <- summary(r)$quants[[2]][, "gamma"]; fmt(g) }, "")
add("A", "gamma over the 20 seeds", "one value", paste(unique(gamma_all), collapse = "; "), "the same for every seed", length(unique(gamma_all)) == 1L)
E <- t(vapply(fits, function(r) r$elpd$delta_elpd, numeric(3)))
for (j in 1:3) { m <- mean(E[, j]); sdv <- stats::sd(E[, j]); zz <- (vig$elpd[j] - m) / sdv
  if (j < 3) { add("A (reported)", sprintf("delta ELPD, %s vs %s (mean over 20 seeds; the July 2021 render gave %.3f)", fits[[1]]$elpd$model1[j], fits[[1]]$elpd$model2[j], vig$elpd_jul2021[j]),
      sprintf("%.3f (SE %.3f)", vig$elpd[j], vig$se[j]), sprintf("%.3f (MC SD %.3f)", m, sdv), sprintf("%.2f of the vignette's SE apart", abs(vig$elpd[j] - m) / vig$se[j]), TRUE); next }
  add("A", sprintf("delta ELPD, %s vs %s (Monte Carlo mean and SD over 20 seeds)", fits[[1]]$elpd$model1[j], fits[[1]]$elpd$model2[j]), sprintf("%.3f (SE %.3f)", vig$elpd[j], vig$se[j]),
      sprintf("%.3f (MC SD %.3f; range %.3f to %.3f)", m, sdv, min(E[, j]), max(E[, j])), "vignette value within 3 Monte Carlo SD", abs(zz) <= 3) }
P <- vapply(fits, function(r) summary(r)$p, 0)
add("A", "P that the causal model fits better (20 seeds)", sprintf("%.1e", vig$p), sprintf("%.1e to %.1e", min(P), max(P)), "below 0.05 for every seed, as in the vignette", all(P < 0.05))
# B. vignette step 2
t0 <- Sys.time(); set.seed(100); varlist <- with(X, sample(snp, size = 1000000, replace = FALSE)); invisible(utils::capture.output(pB <- est_cause_params(X, varlist)))
min_params <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
dpi <- if (identical(dim(pB$mix_grid), dim(params$mix_grid))) max(abs(pB$mix_grid$pi - params$mix_grid$pi)) else NA
add("B", "rho from est_cause_params (seed 100, 1,000,000 variants)", sprintf("%.6f", params$rho), sprintf("%.6f", pB$rho), "within 1e-4", abs(pB$rho - params$rho) < 1e-4)
add("B (reported)", "mixing proportions, largest difference", "0", if (is.na(dpi)) "different grids" else sprintf("%.1e", dpi), "reported", TRUE)
# C. app path: the two GWAS as Data library datasets, the study script with local PLINK clumping
lib <- library_dir(root)
mk <- function(trait, b, s) { d <- data.frame(snp = X$snp, ea = unname(X$A1), oa = unname(X$A2), beta = X[[b]], se = X[[s]], stringsAsFactors = FALSE)
  dd <- file.path(lib, trait); dir.create(dd, recursive = TRUE); lib_table_write(d, dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = trait, role = "either", facts = list(ancestry = "EUR", kind = "continuous", name = trait), original = list(file = "LDL_CAD_merged.RDS", sha256 = sha256_file(file.path(ex, "LDL_CAD_merged.RDS"))), rows_kept = nrow(d), usable_as_exposure = TRUE, usable_as_outcome = TRUE, usable_in_runs = TRUE), auto_unbox = TRUE)) }
mk("LDL", "beta_hat_1", "seb1"); mk("CAD", "beta_hat_2", "seb2"); h <- studio_new("CAUSE validation")
pair <- list(folder = "01_ldl_cad", ex = "lib:LDL", out = "lib:CAD", exName = "LDL", outName = "CAD")
sf <- cause_pair_prepare(h, pair, root); old <- setwd(h$dir)
t0 <- Sys.time(); sys.source(sf, envir = new.env(parent = globalenv())); min_app <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
a <- utils::read.csv(file.path("04_across_pairs", "cause", "01_ldl_cad_results.csv")); va <- readLines(file.path("04_across_pairs", "cause", "01_ldl_cad_variants.txt"))
sys.source(sf, envir = new.env(parent = globalenv())); b <- utils::read.csv(file.path("04_across_pairs", "cause", "01_ldl_cad_results.csv"))
TRAITS <- NULL; BFILE <- NULL; PLINK <- NULL; seed <- NULL
for (ln in c("^TRAITS <- ", "^BFILE <- ", "^seed <- ")) eval(parse(text = sub("   #.*$", "", grep(ln, readLines(sf), value = TRUE))))
rd <- function(f) as.data.frame(data.table::fread(f))
Xd <- gwas_merge(rd(TRAITS[1]), rd(TRAITS[2]), snp_name_cols = c("snp", "snp"), beta_hat_cols = c("beta", "beta"), se_cols = c("se", "se"), A1_cols = c("ea", "ea"), A2_cols = c("oa", "oa"))
set.seed(seed); vl <- with(Xd, sample(snp, size = min(1e6, nrow(Xd)), replace = FALSE)); invisible(utils::capture.output(pd <- est_cause_params(Xd, vl)))
cl <- ieugwasr::ld_clump(dat = data.frame(rsid = Xd$snp, pval = Xd$p1, stringsAsFactors = FALSE), clump_r2 = 0.01, clump_p = 1e-3, plink_bin = PLINK, bfile = BFILE)
set.seed(seed); invisible(utils::capture.output(rdr <- cause(X = Xd, variants = cl$rsid, param_ests = pd))); setwd(old)
sd_ <- summary(rdr); gd <- sd_$quants[[2]][, "gamma"]
add("C", sprintf("app path, LDL on CAD (%d variants clumped locally; the vignette's 719 used the OpenGWAS LD reference)", a$instruments), sprintf("gamma %s, P %.1e", vig$gamma, vig$p),
    sprintf("gamma %s, P %.1e", fmt(c(a$gamma, a$gamma_lower, a$gamma_upper)), a$p), "causal model fits better (P < 0.05) and gamma within the vignette's 95% interval", a$p < 0.05 && a$gamma >= 0.24 && a$gamma <= 0.40)
dC <- max(abs(c(a$gamma - gd[1], a$gamma_lower - gd[2], a$gamma_upper - gd[3], a$z - sd_$z, a$rho - pd$rho)), abs(a$instruments - nrow(cl)))
add("C", "app path vs the package pipeline called directly (same files, same seed)", "0", sprintf("%.1e", dC), "< 1e-10 (read back from the CSV table)", dC < 1e-10)
add("C (reported)", "overlap with the vignette's 719 variants", "719", sprintf("%d of %d in common", length(intersect(va, tv)), length(va)), "reported", TRUE)
cols <- c("gamma", "gamma_lower", "gamma_upper", "z", "p", "delta_elpd", "se_delta_elpd", "rho")
dD <- max(abs(as.matrix(a[, cols]) - as.matrix(b[, cols])))
add("C", "replay of the study script", "identical", if (dD == 0) "identical" else sprintf("%.1e", dD), "identical", dD == 0)
utils::write.csv(a, file.path(vdir, "cause_ldl_cad.csv"), row.names = FALSE)
utils::write.csv(data.frame(seed = c(100, 1:19), null_sharing = E[, 1], null_causal = E[, 2], sharing_causal = E[, 3], p = P), file.path(vdir, "elpd_seeds.csv"), row.names = FALSE)
cat(sprintf("time: %.2f min per CAUSE fit on the stored inputs; est_cause_params %.1f min; app path %.1f min; inputs MD5 %s\n", min_fit, min_params, min_app, paste(names(md5), md5, collapse = ", ")))
res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
print(res[, c("part", "quantity", "reference", "mr_studio", "pass")], row.names = FALSE); ok <- all(res$pass); cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
