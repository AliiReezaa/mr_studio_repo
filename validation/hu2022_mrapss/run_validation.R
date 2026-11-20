# Validation study 17 (Round 69): MR-APSS (Hu et al. 2022, PNAS 119:e2106858119, doi:10.1073/pnas.2106858119, PMID 35787050),
# GitHub YangLabHKUST/MR-APSS (MIT), pinned. Published values: the authors' result file Traits_MRAPSS.MRres (repository
# YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce), rows Height_UKB -> BMI and BMI -> Height_UKB (UK Biobank, Watanabe et al. 2019).
# Data here: Neale lab UK Biobank BMI (21001_irnt) and standing height (50_irnt), LDSC format (HapMap3, as MR-APSS uses), effects
# Z / sqrt(N) with SE 1 / sqrt(N); the LDSC authors' European LD scores; the app's 1000 Genomes EUR panel and PLINK for clumping.
#   A. The package's bundled example (BMI -> T2D, precomputed background) with the pinned version, recorded with V1.0.0 for comparison.
#   B. The app path (Data library datasets, preparation, study script) in both directions, against the published values (2 combined SE).
#   C. The app path equals the package pipeline called directly on the same prepared files with the same seed; D. the script replays.
#   MRSTUDIO_VALIDATE_MRAPSS=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/hu2022_mrapss/run_validation.R
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_MRAPSS"), "1")) { cat("Skipped: set MRSTUDIO_VALIDATE_MRAPSS=1 to run.\n"); quit(status = 0) }
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "hu2022_mrapss"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("apssval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE", file.path(tmp, "cache")); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"), MRGUI_LDSC_DIR = file.path(cache, "ldsc"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("resource_core.R", "import_core.R", "ld_core.R", "ldsc_core.R", "mrlap_core.R", "mrapss_core.R")) source(file.path(args_dir, f)) })
if (is.na(plink_bin()) && nzchar(Sys.which("plink"))) Sys.setenv(MRGUI_PLINK = Sys.which("plink"))
fetch <- function(url, f, md5) { p <- file.path(cache, f); if (!(file.exists(p) && identical(unname(tools::md5sum(p)), md5))) utils::download.file(url, p, mode = "wb", quiet = TRUE)
  if (!identical(unname(tools::md5sum(p)), md5)) stop("MD5 differs for ", f); p }
z <- lapply(c("hm3", "w_ld"), function(w) fetch(LDSC_REF[[w]]$url, basename(dirname(LDSC_REF[[w]]$url)), LDSC_REF[[w]]$md5)); invisible(ldsc_ref_fetch(files = list(hm3 = z[[1]], w_ld = z[[2]])))
nf <- list(bmi = fetch("https://www.dropbox.com/s/iosq5liwda3gvca/21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "0b5b1157ed364fc6bac0cde6bd81076d"),
           height = fetch("https://www.dropbox.com/s/clygq1yvw3v6c96/50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "1e1fcd5828070f1385d13f6a515bb624"))
pub <- utils::read.table(fetch("https://raw.githubusercontent.com/YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce/master/Traits_MRAPSS.MRres", "Traits_MRAPSS.MRres", "d616f1c9335fa0f0abf5858bb56de188"), header = TRUE, stringsAsFactors = FALSE)
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)
# A. bundled example with the pinned version
suppressPackageStartupMessages(library(MRAPSS)); e <- new.env(); utils::data("C", "Omega", "MRdat", package = "MRAPSS", envir = e)
invisible(utils::capture.output(ra <- MRAPSS(e$MRdat, exposure = "BMI", outcome = "T2D", C = e$C, Omega = e$Omega, Cor.SelectionBias = TRUE)))
add("A (recorded)", "bundled example BMI -> T2D (1,227 instruments), pinned version", "V1.0.0 on the same example: 0.209 (0.033)", sprintf("%.3f (%.3f)", ra$beta, ra$beta.se), "within 2 combined SE of V1.0.0", abs(ra$beta - 0.2090028) / sqrt(ra$beta.se^2 + 0.03299989^2) < 2)
# B, C, D. app path, both directions
lib <- library_dir(root); mk <- function(f, trait) { x <- utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
  d <- data.frame(snp = x$SNP, ea = x$A1, oa = x$A2, beta = x$Z / sqrt(x$N), se = 1 / sqrt(x$N), n = x$N, stringsAsFactors = FALSE)
  dd <- file.path(lib, trait); dir.create(dd, recursive = TRUE); lib_table_write(d, dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = trait, role = "either", facts = list(ancestry = "EUR", kind = "continuous", name = trait), original = list(file = basename(f), sha256 = sha256_file(f)), rows_kept = nrow(d), usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE), auto_unbox = TRUE)) }
mk(nf$bmi, "BMI"); mk(nf$height, "Height"); h <- studio_new("MR-APSS validation")
run_pair <- function(ex, out, folder) { pair <- list(folder = folder, ex = paste0("lib:", ex), out = paste0("lib:", out), exName = ex, outName = out)
  sf <- mrapss_pair_prepare(h, pair, root); old <- setwd(h$dir); on.exit(setwd(old))
  sys.source(sf, envir = new.env(parent = globalenv())); a <- utils::read.csv(file.path("04_across_pairs", "mrapss", paste0(folder, "_results.csv")))
  sys.source(sf, envir = new.env(parent = globalenv())); b <- utils::read.csv(file.path("04_across_pairs", "mrapss", paste0(folder, "_results.csv")))
  TRAITS <- NULL; BFILE <- NULL; PLINK <- NULL; LD_SCORES <- NULL; seed <- NULL; IV_THRESHOLD <- NULL
  for (ln in c("^TRAITS <- ", "^LD_SCORES <- ", "^BFILE <- ", "^seed <- ")) eval(parse(text = sub("   #.*$", "", grep(ln, readLines(sf), value = TRUE))))
  invisible(loadNamespace("ieugwasr")); set.seed(seed)
  invisible(utils::capture.output({ x <- format_data(as.data.frame(data.table::fread(TRAITS[1]))); y <- format_data(as.data.frame(data.table::fread(TRAITS[2])))
    pa <- est_paras(dat1 = x, dat2 = y, trait1.name = "exposure", trait2.name = "outcome", ldscore.dir = LD_SCORES)
    iv <- clump(pa$dat, IV.Threshold = IV_THRESHOLD, SNP_col = "SNP", pval_col = "pval.exp", clump_kb = 1000, clump_r2 = 0.001, bfile = BFILE, plink_bin = PLINK)
    dr <- MRAPSS(iv, exposure = "exposure", outcome = "outcome", C = pa$C, Omega = pa$Omega, Cor.SelectionBias = TRUE) }))
  list(app = a, rep = b, direct = dr, n_direct = nrow(iv)) }
for (pp in list(c("Height", "BMI", "Height_UKB", "BMI"), c("BMI", "Height", "BMI", "Height_UKB"))) {
  r <- run_pair(pp[1], pp[2], paste0("0", match(pp[1], c("Height", "BMI")), "_", tolower(pp[1]), "_", tolower(pp[2]))); P <- pub[pub$exposure == pp[3] & pub$outcome == pp[4], ][1, ]
  zz <- (r$app$estimate - P$beta.hat) / sqrt(r$app$se^2 + P$se^2)
  add("B", sprintf("%s -> %s vs the paper (%s instruments there)", pp[1], pp[2], P$nsnp), sprintf("%.4f (%.4f)", P$beta.hat, P$se), sprintf("%.4f (%.4f), %d instruments", r$app$estimate, r$app$se, r$app$instruments), "within 2 combined SE", abs(zz) < 2)
  dC <- max(abs(c(r$app$estimate - r$direct$beta, r$app$se - r$direct$beta.se)), abs(r$app$instruments - r$n_direct))
  add("C", sprintf("%s -> %s: app path vs the package pipeline called directly (same files, same seed)", pp[1], pp[2]), "0", sprintf("%.1e", dC), "< 1e-10 (read back from the CSV table)", dC < 1e-10)
  dD <- max(abs(as.matrix(r$app[, c("estimate", "se", "p")]) - as.matrix(r$rep[, c("estimate", "se", "p")])))
  add("D", sprintf("%s -> %s: replay of the study script", pp[1], pp[2]), "identical", if (dD == 0) "identical" else sprintf("%.1e", dD), "identical", dD == 0)
  utils::write.csv(r$app, file.path(vdir, sprintf("mrapss_%s_%s.csv", tolower(pp[1]), tolower(pp[2]))), row.names = FALSE) }
res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
print(res[, c("part", "reference", "mr_studio", "pass")], row.names = FALSE); ok <- all(res$pass); cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
