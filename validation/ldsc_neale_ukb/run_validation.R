# Validation study 14 (Round 64): LD scores and cross-trait LD score regression (Bulik-Sullivan et al. 2015, Nat Genet 47:291-295 and
# 47:1236-1241; PMIDs 25642630, 26414676) with GenomicSEM::ldsc (Grotzinger et al. 2019, PMID 30962613), pinned.
# Data: the Neale lab UK Biobank round 2 summary statistics in LDSC format for standing height (50_irnt) and BMI (21001_irnt), both
# sexes, and the lab's published LDSC table (ukb31063_h2_all.02Oct2019); the LDSC authors' European LD scores and HapMap3 list from
# their CC BY copies on Zenodo; the 1000 Genomes Phase 3 EUR panel of the app. Every download is checked against its MD5 below.
#   A. LD scores computed by MR Studio equal the corrected r2 sums from PLINK --r2 (chromosome 22, 200 variants).
#   B. LD scores computed from the 1000 Genomes Phase 3 panel against the official scores (1000 Genomes Phase 1).
#   C. The app's path (Data library datasets with alleles swapped for half the variants, preparation, study script, replay) equals
#      GenomicSEM::ldsc run directly on the Neale files (restricted
#      to the HapMap3 list with matching alleles, as munge_sumstats --merge-alleles does).
#   D. Against the Neale lab's published values (lambda GC, intercept, SNP heritability; their h2 is from partitioned LDSC).
#   E. Results with the computed LD scores against those with the official ones.
# Downloads about 130 MB and computes LD scores for 22 chromosomes (10 to 30 minutes), so it runs only on request:
#   MRSTUDIO_VALIDATE_LDSC=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/ldsc_neale_ukb/run_validation.R
# A folder of LD scores computed before (MRSTUDIO_LDSC_COMPUTED=<folder>) skips the computation.
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_LDSC"), "1")) { cat("Skipped: set MRSTUDIO_VALIDATE_LDSC=1 to run (downloads about 130 MB).\n"); quit(status = 0) }
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "ldsc_neale_ukb"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("ldscval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
cache <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE", file.path(tmp, "cache")); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"), MRGUI_LDSC_DIR = file.path(cache, "ldsc"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "ld_core.R", "ldsc_core.R")) source(file.path(args_dir, f)) })
if (is.na(plink_bin()) && nzchar(Sys.which("plink"))) Sys.setenv(MRGUI_PLINK = Sys.which("plink"))
get <- function(url, f, md5) { p <- file.path(cache, f); if (!(file.exists(p) && identical(unname(tools::md5sum(p)), md5))) utils::download.file(url, p, mode = "wb", quiet = TRUE)
  if (!identical(unname(tools::md5sum(p)), md5)) stop("MD5 differs for ", f); p }
NEALE <- list(height = c("https://www.dropbox.com/s/clygq1yvw3v6c96/50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "50_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "1e1fcd5828070f1385d13f6a515bb624"),
              bmi = c("https://www.dropbox.com/s/iosq5liwda3gvca/21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1", "21001_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz", "0b5b1157ed364fc6bac0cde6bd81076d"),
              h2 = c("https://www.dropbox.com/s/ipeqyhrpdqav5uh/ukb31063_h2_all.02Oct2019.tsv.gz?dl=1", "ukb31063_h2_all.02Oct2019.tsv.gz", "dc75cdbd8343c934feb4631050257f7f"))
nf <- lapply(NEALE, function(x) get(x[1], x[2], x[3]))
z <- lapply(c("hm3", "w_ld"), function(w) get(LDSC_REF[[w]]$url, basename(dirname(LDSC_REF[[w]]$url)), LDSC_REF[[w]]$md5))
invisible(ldsc_ref_fetch(files = list(hm3 = z[[1]], w_ld = z[[2]]))); off <- ldsc_ld_folder("official")
rd <- function(f) utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
out <- list(); add <- function(part, quantity, reference, mr_studio, criterion, pass) out[[length(out) + 1L]] <<- data.frame(part, quantity, reference, mr_studio, criterion, pass, stringsAsFactors = FALSE)

# A. PLINK
pd <- file.path(tmp, "plinkA"); dir.create(pd); writeLines(utils::read.table(file.path(ldsc_ref_dir(), "w_hm3.snplist"), header = TRUE)$SNP, file.path(pd, "h.txt"))
system2(plink_bin(), c("--bfile", shQuote(file.path(ld_panel_dir(), "EUR")), "--chr 22 --extract", file.path(pd, "h.txt"), "--make-bed --out", file.path(pd, "c"), "--silent"))
dA <- ldsc_l2_bfile(file.path(pd, "c")); set.seed(64); pick <- sort(sample(nrow(dA), 200)); writeLines(dA$SNP[pick], file.path(pd, "p.txt"))
system2(plink_bin(), c("--bfile", file.path(pd, "c"), "--r2 --ld-snp-list", file.path(pd, "p.txt"), "--ld-window 999999 --ld-window-kb 1000 --ld-window-r2 0 --out", file.path(pd, "r"), "--silent"))
r <- utils::read.table(file.path(pd, "r.ld"), header = TRUE); nA <- length(readLines(file.path(pd, "c.fam"))); r <- r[abs(r$BP_B - r$BP_A) <= 1e6 & r$SNP_B %in% dA$SNP, ]
s <- tapply(r$R2 - (1 - r$R2) / (nA - 2), r$SNP_A, sum)   # PLINK lists each variant with itself (r2 = 1), as the LD score includes it
dif <- max(abs(dA$L2[pick] - as.numeric(s[dA$SNP[pick]])))
add("A", "largest difference from PLINK r2 sums (200 variants, chr22, 1 Mb)", "0", sprintf("%.1e", dif), "< 1e-4 (PLINK prints 6 digits)", dif < 1e-4)

# B. computed against official
cmp <- Sys.getenv("MRSTUDIO_LDSC_COMPUTED", "")
if (!nzchar(cmp)) { ldsc_ldscores_compute("EUR", cm = ldsc_official_cm()); cmp <- ldsc_ld_folder("EUR") } else { dir.create(dirname(ldsc_ld_folder("EUR")), recursive = TRUE, showWarnings = FALSE); file.symlink(normalizePath(cmp), ldsc_ld_folder("EUR")) }
B <- do.call(rbind, lapply(1:22, function(c) merge(rd(file.path(off, paste0(c, ".l2.ldscore.gz")))[, c("SNP", "L2")], rd(file.path(cmp, paste0(c, ".l2.ldscore.gz")))[, c("SNP", "L2")], by = "SNP", suffixes = c("_off", "_new"))))
nOff <- sum(vapply(1:22, function(c) nrow(rd(file.path(off, paste0(c, ".l2.ldscore.gz")))), 0)); cB <- stats::cor(B$L2_off, B$L2_new); rt <- mean(B$L2_new) / mean(B$L2_off); cov <- nrow(B) / nOff
add("B", "correlation of computed and official LD scores (genome)", "1", sprintf("%.4f", cB), "> 0.98", cB > 0.98)
add("B", "mean computed / mean official LD score", sprintf("%.3f (variants in the panel / official)", cov), sprintf("%.3f", rt), "within 0.03 of the panel's share of HapMap3 variants", abs(rt - cov) < 0.03)

# C. app path against GenomicSEM directly
lib <- library_dir(root); set.seed(1); mk <- function(f, trait) { x <- utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE); sw <- stats::runif(nrow(x)) < 0.5
  d <- data.frame(snp = x$SNP, ea = ifelse(sw, x$A2, x$A1), oa = ifelse(sw, x$A1, x$A2), beta = ifelse(sw, -x$Z, x$Z), se = 1, n = x$N, stringsAsFactors = FALSE)
  dd <- file.path(lib, trait); dir.create(dd, recursive = TRUE); lib_table_write(d, dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = trait, role = "either", facts = list(ancestry = "EUR", kind = "continuous", name = trait), original = list(file = basename(f), sha256 = sha256_file(f)), usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE), auto_unbox = TRUE)) }
mk(nf$height, "Height"); mk(nf$bmi, "BMI")
h <- studio_new("LDSC validation"); pair <- list(folder = "01_height_bmi", ex = "lib:Height", out = "lib:BMI", exName = "Height", outName = "BMI")
app <- ldsc_pair_run(h, pair, "official")
old <- setwd(h$dir); rep <- { sys.source(file.path("04_across_pairs", "ldsc", "01_height_bmi.R"), envir = new.env(parent = globalenv())); utils::read.csv(file.path("04_across_pairs", "ldsc", "01_height_bmi_results.csv")) }; setwd(old)
# Reference input, filtered independently of the app as munge_sumstats --merge-alleles does: variants not in w_hm3.snplist or with other alleles
# removed (the official LD score files hold more variants than the list, so GenomicSEM on the raw files would also use those)
hl <- utils::read.table(file.path(ldsc_ref_dir(), "w_hm3.snplist"), header = TRUE, stringsAsFactors = FALSE); dl <- file.path(tmp, "direct"); dir.create(dl); nbad <- c()
ref <- vapply(c(nf$height, nf$bmi), function(f) { x <- utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE); i <- match(x$SNP, hl$SNP)
  bad <- is.na(i) | !((x$A1 == hl$A1[i] & x$A2 == hl$A2[i]) | (x$A1 == hl$A2[i] & x$A2 == hl$A1[i])); nbad <<- c(nbad, sum(bad)); o <- file.path(dl, basename(f)); utils::write.table(x[!bad, ], gzfile(o), sep = "\t", quote = FALSE, row.names = FALSE); o }, "")
old <- setwd(dl); invisible(utils::capture.output(g <- GenomicSEM::ldsc(traits = ref, sample.prev = c(NA, NA), population.prev = c(NA, NA), ld = off, wld = off, trait.names = c("X", "Y"), stand = TRUE, ldsc.log = "direct"))); setwd(old)
dir_v <- c(g$S[1, 1], g$S[2, 2], g$S_Stand[2, 1], g$S[2, 1], g$I[1, 1], g$I[2, 2], g$I[2, 1]); dC <- max(abs(app$estimate - dir_v)); dR <- max(abs(app$estimate - rep$estimate), abs(app$se - rep$se))
add("C", "largest difference, app path vs GenomicSEM::ldsc on the Neale files (h2, rg, covariance, intercepts)", "0", sprintf("%.1e", dC), "< 1e-10", dC < 1e-10)
add("C", "variants not in w_hm3.snplist or with other alleles, removed (height, BMI)", "munge_sumstats rule", paste(nbad, collapse = ", "), "reported", TRUE)
add("C", "h2 and rg estimates, app path vs direct", "0", sprintf("%.1e", max(abs(app$estimate[1:4] - dir_v[1:4]))), "< 1e-10", max(abs(app$estimate[1:4] - dir_v[1:4])) < 1e-10)
add("C", "replay of the study script", "identical", if (dR == 0) "identical" else sprintf("%.1e", dR), "identical", dR == 0)

# D. Neale lab published values
H <- utils::read.delim(gzfile(nf$h2), stringsAsFactors = FALSE); H <- H[H$sex == "both_sexes" & H$phenotype %in% c("50_irnt", "21001_irnt"), ]
lgA <- readLines(file.path(h$dir, "04_across_pairs", "ldsc", "01_height_bmi_ldsc.log")); lam <- as.numeric(sub(".*: ", "", grep("^Lambda GC", lgA, value = TRUE)))
for (k in 1:2) { ph <- c("50_irnt", "21001_irnt")[k]; nm <- c("Height", "BMI")[k]; P <- H[H$phenotype == ph, ]; a <- app[app$trait == nm, ]
  add("D", paste("lambda GC,", nm), sprintf("%.4f", P$lambdaGC), sprintf("%.4f", lam[c(1, length(lam))[k]]), "within 0.005", abs(lam[c(1, length(lam))[k]] - P$lambdaGC) < 0.005)
  ai <- a[a$quantity == "intercept", ]; ah <- a[a$quantity == "h2", ]
  zi <- (ai$estimate - P$intercept) / sqrt(ai$se^2 + P$intercept_se^2); zh <- (ah$estimate - P$h2_observed) / sqrt(ah$se^2 + P$h2_observed_se^2)
  add("D", paste("intercept,", nm), sprintf("%.3f (%.3f)", P$intercept, P$intercept_se), sprintf("%.3f (%.3f)", ai$estimate, ai$se), "within 2 combined SE", abs(zi) < 2)
  add("D", paste("SNP heritability,", nm), sprintf("%.3f (%.3f), partitioned", P$h2_observed, P$h2_observed_se), sprintf("%.3f (%.3f)", ah$estimate, ah$se), "within 2 combined SE", abs(zh) < 2) }

# E. computed LD scores
appE <- ldsc_pair_run(h, modifyList(pair, list(folder = "02_height_bmi_eur")), "EUR")
for (q in c("h2", "rg")) for (t in unique(app$trait[app$quantity == q])) { a <- app[app$quantity == q & app$trait == t, ]; e <- appE[appE$quantity == q & appE$trait == t, ]; zz <- abs(a$estimate - e$estimate) / a$se
  add("E", paste(q, if (t == "both") "" else t, "with computed vs official LD scores"), sprintf("%.3f (%.3f)", a$estimate, a$se), sprintf("%.3f (%.3f)", e$estimate, e$se), "within 1 SE of the official-score estimate", zz < 1) }

res <- do.call(rbind, out); utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
utils::write.csv(cbind(official = app[, c("quantity", "trait", "estimate", "se", "p")], computed_eur = appE[, c("estimate", "se")]), file.path(vdir, "ldsc_height_bmi.csv"), row.names = FALSE)
print(res[, c("part", "quantity", "reference", "mr_studio", "pass")], row.names = FALSE); ok <- all(res$pass)
cat(sprintf("VALIDATION %s\n", if (ok) "PASSED" else "FAILED")); if (!ok) quit(status = 1)
