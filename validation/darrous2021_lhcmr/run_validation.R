# Validation study 15 (Round 65): LHC-MR (Darrous et al. 2021, Nat Commun 12:7274, doi:10.1038/s41467-021-26970-w, PMID 34907193),
# GitHub LizaDarrous/lhcMR (MIT), pinned. Pair BMI (X) and standing height (Y), Neale lab UK Biobank round 2, both sexes, as in the
# paper's Supplementary Data 2 (CC BY 4.0). The paper used the full Neale files (about 4.6 million variants of the authors' LD file);
# this study uses the Neale lab's LDSC-format files (HapMap3 variants, about 1 million), so its standard errors are larger.
#   A. The app path (Data library datasets, preparation, study script calling lhcMR) replays identically with the same seed.
#   B. The estimates agree with the paper's Supplementary Data 2 within two combined standard errors.
# Runs only on request (about 600 MB of downloads; the fits take one to two hours on 8 cores):
#   MRSTUDIO_VALIDATE_LHCMR=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/darrous2021_lhcmr/run_validation.R
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_LHCMR"), "1")) { cat("Skipped: set MRSTUDIO_VALIDATE_LHCMR=1 to run.\n"); quit(status = 0) }
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "darrous2021_lhcmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("lhcval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE", file.path(tmp, "cache")); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"), MRGUI_LDSC_DIR = file.path(cache, "ldsc"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir)); for (f in c("import_core.R", "ld_core.R", "ldsc_core.R", "lhcmr_core.R")) source(file.path(args_dir, f)) })
get <- function(url, f, md5) { p <- file.path(cache, f); if (!(file.exists(p) && identical(unname(tools::md5sum(p)), md5))) utils::download.file(url, p, mode = "wb", quiet = TRUE)
  if (!identical(unname(tools::md5sum(p)), md5)) stop("MD5 differs for ", f); p }
nf <- list(bmi = get("https://www.dropbox.com/s/iosq5liwda3gvca/21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "0b5b1157ed364fc6bac0cde6bd81076d"),
           height = get("https://www.dropbox.com/s/clygq1yvw3v6c96/50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "1e1fcd5828070f1385d13f6a515bb624"))
z <- lapply(c("hm3", "w_ld"), function(w) get(LDSC_REF[[w]]$url, basename(dirname(LDSC_REF[[w]]$url)), LDSC_REF[[w]]$md5)); invisible(ldsc_ref_fetch(files = list(hm3 = z[[1]], w_ld = z[[2]])))
lf <- lapply(c("ld", "rho"), function(w) get(LHCMR_REF[[w]]$url, LHCMR_REF[[w]]$file, LHCMR_REF[[w]]$md5)); invisible(lhcmr_ref_fetch(files = list(ld = lf[[1]], rho = lf[[2]])))
sd2 <- get("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-26970-w/MediaObjects/41467_2021_26970_MOESM5_ESM.xlsx", "41467_2021_26970_MOESM5_ESM.xlsx", "a9122a7299c670d3c62ebba1aead0c2b")

# Standardised effects of inverse-normal transformed traits: beta = Z / sqrt(N), SE = 1 / sqrt(N)
lib <- library_dir(root); mk <- function(f, trait) { x <- utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
  d <- data.frame(snp = x$SNP, ea = x$A1, oa = x$A2, beta = x$Z / sqrt(x$N), se = 1 / sqrt(x$N), n = x$N, stringsAsFactors = FALSE)
  dd <- file.path(lib, trait); dir.create(dd, recursive = TRUE); lib_table_write(d, dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = trait, role = "either", facts = list(ancestry = "EUR", kind = "continuous", name = trait), original = list(file = basename(f), sha256 = sha256_file(f)), usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE), auto_unbox = TRUE)) }
mk(nf$bmi, "BMI"); mk(nf$height, "Height")
h <- studio_new("LHC-MR validation"); pair <- list(folder = "01_bmi_height", ex = "lib:BMI", out = "lib:Height", exName = "BMI", outName = "Height")
t0 <- Sys.time(); app <- lhcmr_pair_run(h, pair); cat("app path:", format(Sys.time() - t0), "\n")
# A. replay (the script again from the study folder) and lhcMR directly on the prepared inputs, same seed
old <- setwd(h$dir); file.rename(file.path("04_across_pairs", "lhcmr", "01_bmi_height_results.csv"), file.path(tmp, "first.csv"))
sys.source(file.path("04_across_pairs", "lhcmr", "01_bmi_height.R"), envir = new.env(parent = globalenv())); rep <- utils::read.csv(file.path("04_across_pairs", "lhcmr", "01_bmi_height_results.csv")); setwd(old)
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)
dR <- max(abs(app$estimate - rep$estimate), abs(app$se - rep$se), na.rm = TRUE)
add("A", "replay of the study script (same seed)", "identical", if (dR == 0) "identical" else sprintf("%.1e", dR), "identical", dR == 0)
# B. paper
P <- as.data.frame(readxl::read_excel(sd2)); P <- P[P$Pair == "BMI-SHeight", ]
for (p in c("axy", "ayx", "tX", "tY", "h2X", "h2Y", "iXY")) { a <- app[app$parameter == p, ]; ps <- as.numeric(P[[paste0(p, "_JKse")]]); pe <- as.numeric(P[[p]])
  zz <- (a$estimate - pe) / sqrt(a$se^2 + ps^2); add("B", c(axy = "BMI on height", ayx = "height on BMI", tX = "confounder on BMI", tY = "confounder on height", h2X = "direct h2 BMI", h2Y = "direct h2 height", iXY = "cross-trait intercept")[[p]],
    sprintf("%.3f (%.3f)", pe, ps), sprintf("%.3f (%.3f)", a$estimate, a$se), "within 2 combined SE", is.finite(zz) && abs(zz) < 2) }
res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE); utils::write.csv(app, file.path(vdir, "lhcmr_bmi_height.csv"), row.names = FALSE)
print(res, row.names = FALSE); ok <- all(res$pass); cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
