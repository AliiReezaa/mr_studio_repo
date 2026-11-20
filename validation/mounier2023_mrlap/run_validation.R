# Validation study 16 (Round 68): MRlap (Mounier and Kutalik 2023, Genet Epidemiol 47:314-331, doi:10.1002/gepi.22522, PMID 37036286),
# GitHub n-mounier/MRlap (GPL-2), pinned. Data: the package's own Example A (simulated BMI and SBP, about 100,000 people each,
# 1.15 million variants, no sample overlap), bundled with the package together with its stored result (data/A.RDS).
#   A. MRlap called directly equals the stored result: every deterministic quantity exactly; the corrected effect's standard error,
#      which MRlap simulates (random draws, no seed in the example), within three Monte Carlo standard deviations (five seeds).
#   B. The app path (Data library datasets, preparation, study script) gives the stored result.
#   C. The app path equals MRlap called directly on the same prepared files with the same seed (to the precision of its CSV table);
#   D. the script replays identically.
# Needs the LDSC authors' LD scores and HapMap3 list (validation/ldsc_neale_ukb fetches them):
#   MRSTUDIO_VALIDATE_MRLAP=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/mounier2023_mrlap/run_validation.R
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_MRLAP"), "1")) { cat("Skipped: set MRSTUDIO_VALIDATE_MRLAP=1 to run.\n"); quit(status = 0) }
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "mounier2023_mrlap"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mrlapval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE", file.path(tmp, "cache")); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"), MRGUI_LDSC_DIR = file.path(cache, "ldsc"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("resource_core.R", "import_core.R", "ld_core.R", "ldsc_core.R", "mrlap_core.R")) source(file.path(args_dir, f)) })
fetch <- function(url, f, md5) { p <- file.path(cache, f); if (!(file.exists(p) && identical(unname(tools::md5sum(p)), md5))) utils::download.file(url, p, mode = "wb", quiet = TRUE)
  if (!identical(unname(tools::md5sum(p)), md5)) stop("MD5 differs for ", f); p }
z <- lapply(c("hm3", "w_ld"), function(w) fetch(LDSC_REF[[w]]$url, basename(dirname(LDSC_REF[[w]]$url)), LDSC_REF[[w]]$md5)); invisible(ldsc_ref_fetch(files = list(hm3 = z[[1]], w_ld = z[[2]])))
LD <- ldsc_ld_folder("official"); HM3 <- file.path(ldsc_ref_dir(), "w_hm3.snplist"); invisible(loadNamespace("ieugwasr"))
ex <- system.file("data", c("BMI_Data.tsv.gz", "SBP_Data.tsv.gz"), package = "MRlap"); ref <- readRDS(system.file("data", "A.RDS", package = "MRlap"))
run_direct <- function(f1, f2, seed) { set.seed(seed); invisible(utils::capture.output(r <- MRlap::MRlap(exposure = f1, exposure_name = "exposure", outcome = f2, outcome_name = "outcome", ld = LD, hm3 = HM3, verbose = FALSE))); r }
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)
det <- function(r) c(observed_effect = r$MRcorrection$observed_effect, observed_effect_se = r$MRcorrection$observed_effect_se, m_IVs = r$MRcorrection$m_IVs, corrected_effect = r$MRcorrection$corrected_effect,
                     unlist(r$LDSC), unlist(r$GeneticArchitecture))
# A. direct, five seeds
runs <- lapply(1:5, function(s) run_direct(ex[1], ex[2], 1000L + s))
dA <- max(abs(det(runs[[1]]) - det(ref))); same_iv <- identical(runs[[1]]$MRcorrection$IVs, ref$MRcorrection$IVs)
add("A", "deterministic quantities vs the stored Example A (effects, instruments, all LDSC and architecture values)", "0", sprintf("%.1e", dA), "< 1e-8, same instruments", dA < 1e-8 && same_iv)
se <- vapply(runs, function(r) r$MRcorrection$corrected_effect_se, 0); zz <- abs(ref$MRcorrection$corrected_effect_se - mean(se)) / stats::sd(se)
add("A", "corrected effect SE (simulated)", sprintf("%.5f", ref$MRcorrection$corrected_effect_se), sprintf("%.5f (sd %.5f over 5 seeds)", mean(se), stats::sd(se)), "within 3 Monte Carlo sd", zz < 3)
# B, C, D. the app path
lib <- library_dir(root); mk <- function(f, trait) { x <- as.data.frame(data.table::fread(f))
  d <- data.frame(snp = x$rsid, chr = x$chr, pos = x$pos, ea = x$alt, oa = x$ref, eaf = x$af, beta = x$beta, se = x$se, n = x$N, stringsAsFactors = FALSE)
  dd <- file.path(lib, trait); dir.create(dd, recursive = TRUE); lib_table_write(d, dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = trait, role = "either", facts = list(ancestry = "EUR", kind = "continuous", name = trait), original = list(file = basename(f), sha256 = sha256_file(f)), rows_kept = nrow(d), usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE), auto_unbox = TRUE)) }
mk(ex[1], "BMI"); mk(ex[2], "SBP")
h <- studio_new("MRlap validation"); pair <- list(folder = "01_bmi_sbp", ex = "lib:BMI", out = "lib:SBP", exName = "BMI", outName = "SBP")
sf <- mrlap_pair_prepare(h, pair, root); old <- setwd(h$dir); sys.source(sf, envir = new.env(parent = globalenv())); app <- utils::read.csv(file.path("04_across_pairs", "mrlap", "01_bmi_sbp_results.csv"))
sys.source(sf, envir = new.env(parent = globalenv())); rep <- utils::read.csv(file.path("04_across_pairs", "mrlap", "01_bmi_sbp_results.csv"))
seed <- as.integer(sub(";.*", "", sub("^seed <- ([0-9]+)L.*", "\\1", grep("^seed <- ", readLines(sf), value = TRUE))))
TRAITS <- NULL; eval(parse(text = sub("   #.*$", "", grep("^TRAITS <- ", readLines(sf), value = TRUE)))); dr <- run_direct(TRAITS[1], TRAITS[2], seed); setwd(old)
g <- function(q) app[app$quantity == q, ]
bB <- c(abs(g("observed_effect")$estimate - ref$MRcorrection$observed_effect), abs(g("corrected_effect")$estimate - ref$MRcorrection$corrected_effect), abs(g("observed_effect")$instruments - ref$MRcorrection$m_IVs), abs(g("observed_effect")$rg - ref$LDSC$rg))
add("B", "app path vs the stored Example A (observed and corrected effect, instruments, rg)", "0", sprintf("%.1e", max(bB)), "< 1e-6", max(bB) < 1e-6)
cC <- max(abs(c(g("observed_effect")$estimate - dr$MRcorrection$observed_effect, g("corrected_effect")$estimate - dr$MRcorrection$corrected_effect, g("corrected_effect")$se - dr$MRcorrection$corrected_effect_se, g("difference")$p - dr$MRcorrection$p_difference)))
add("C", "app path vs MRlap called directly on the same files, same seed (effects, simulated SE, difference test)", "0", sprintf("%.1e", cC), "< 1e-10 (the app's values are read back from its CSV table)", cC < 1e-10)
dD <- max(abs(as.matrix(app[, c("estimate", "se", "p")]) - as.matrix(rep[, c("estimate", "se", "p")])), na.rm = TRUE)
add("D", "replay of the study script", "identical", if (dD == 0) "identical" else sprintf("%.1e", dD), "identical", dD == 0)
res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); utils::write.csv(app, file.path(vdir, "mrlap_example_a.csv"), row.names = FALSE)
print(res[, c("part", "reference", "mr_studio", "pass")], row.names = FALSE); ok <- all(res$pass); cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
