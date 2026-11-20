# ---- Round 64: LD scores and cross-trait LD score regression (Bulik-Sullivan et al. 2015, Nat Genet 47:291 and 47:1236) ----------
# LD score regression needs LD scores of the HapMap3 variants. Two sources:
#   1. the official European LD scores of the LDSC authors (eur_w_ld_chr: 1000 Genomes Europeans, 1 cM windows, sum of r2 over
#      HapMap3 variants), from the copy on Zenodo (CC BY 4.0) since the Broad Institute links went offline; pinned by MD5;
#   2. LD scores computed here from the 1000 Genomes Phase 3 panel the app already uses for clumping (any of the five
#      super-populations), with the estimator of the LDSC software (bulik/ldsc 2fdeeb3, ldscore/ldscore.py): genotypes centred
#      with missing values set to the mean, scaled by their standard deviation, r2 corrected for sample size as r2 - (1 - r2)/(n - 2),
#      summed over the HapMap3 variants within the window (both sides, the variant itself included).
# The regression itself is GenomicSEM::ldsc (Grotzinger et al. 2019, Nat Hum Behav 3:513), which reimplements LDSC in R.
LDSC_REF <- list(
  w_ld = list(url = "https://zenodo.org/api/records/8182036/files/eur_w_ld_chr.tar.gz/content", md5 = "e2f16343c4cfaa76caa7d0c03d26b489",
              bytes = 33357890, doi = "10.5281/zenodo.8182036", what = "LD scores of the LDSC authors (eur_w_ld_chr), 1000 Genomes Europeans"),
  hm3 = list(url = "https://zenodo.org/api/records/7773502/files/w_hm3.snplist.gz/content", md5 = "153ecc2bcfa740afafe656e6a384d769",
             bytes = 5221397, doi = "10.5281/zenodo.7773502", what = "HapMap3 variant list of the LDSC authors (w_hm3.snplist)"))
LDSC_PIN <- "2fdeeb3b44379408794154993dbd6101b8946b7e"   # bulik/ldsc commit whose estimator is reproduced

ldsc_ref_dir <- function() { d <- Sys.getenv("MRGUI_LDSC_DIR", ""); if (nzchar(d)) d else file.path(tools::R_user_dir("MRStudio", "data"), "ldsc") }
.ldsc_md5 <- function(f) unname(tools::md5sum(f))

#' Download and check the official files (each is checked against its pinned MD5 before use). Returns the status.
ldsc_ref_fetch <- function(what = c("hm3", "w_ld"), dir = ldsc_ref_dir(), files = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (w in what) {
    r <- LDSC_REF[[w]]; f <- files[[w]] %||% file.path(dir, basename(dirname(r$url)))
    if (is.null(files[[w]]) && !(file.exists(f) && identical(.ldsc_md5(f), r$md5))) utils::download.file(r$url, f, mode = "wb", quiet = TRUE)
    if (!identical(.ldsc_md5(f), r$md5)) stop(sprintf("The downloaded %s is not the pinned file (MD5 differs).", r$what))
    if (w == "hm3") { x <- utils::read.table(gzfile(f), header = TRUE, stringsAsFactors = FALSE); utils::write.table(x, file.path(dir, "w_hm3.snplist"), sep = "\t", quote = FALSE, row.names = FALSE) }
    else utils::untar(f, exdir = dir)
    writeLines(c(paste("source:", r$url), paste("DOI:", r$doi), paste("MD5:", r$md5), "licence: CC BY 4.0", paste("installed:", format(Sys.Date()))), file.path(dir, paste0("SOURCE_", w, ".txt")))
  }
  ldsc_ref_status(dir)
}

.ldsc_complete <- function(d) dir.exists(d) && all(file.exists(file.path(d, paste0(1:22, ".l2.ldscore.gz")))) && all(file.exists(file.path(d, paste0(1:22, ".l2.M_5_50"))))
ldsc_ref_status <- function(dir = ldsc_ref_dir()) {
  pops <- if (dir.exists(file.path(dir, "computed"))) Filter(function(p) .ldsc_complete(file.path(dir, "computed", p)), list.files(file.path(dir, "computed"))) else character()
  list(dir = dir, hm3 = file.exists(file.path(dir, "w_hm3.snplist")), official = .ldsc_complete(file.path(dir, "eur_w_ld_chr")), computed = as.list(pops))
}
#' Folder of the LD scores for a source: "official" (European) or a super-population code of the computed scores.
ldsc_ld_folder <- function(source = "official", dir = ldsc_ref_dir()) if (identical(source, "official")) file.path(dir, "eur_w_ld_chr") else file.path(dir, "computed", source)

# PLINK .bed decoding: 2 bits per person, first person in the lowest bits; codes 00, 01, 10, 11 = 2, NA, 1, 0 copies of the A1 allele.
.BED_LUT <- local({ v <- c(2, NA, 1, 0); t(vapply(0:255, function(b) v[bitwAnd(bitwShiftR(b, c(0L, 2L, 4L, 6L)), 3L) + 1L], numeric(4))) })
.bed_cols <- function(raw, j, n) { g <- .BED_LUT[as.integer(raw[, j, drop = FALSE]) + 1L, , drop = FALSE]
  matrix(t(g), ncol = length(j))[seq_len(n), , drop = FALSE] }   # people x variants (byte rows interleave 4 people each)
.ldsc_std <- function(g) { m <- colMeans(g, na.rm = TRUE); i <- which(is.na(g), arr.ind = TRUE); if (nrow(i)) g[i] <- m[i[, 2]]
  g <- sweep(g, 2, m); s <- sqrt(colMeans(g^2)); s[s == 0] <- 1; sweep(g, 2, s, "/") }

#' LD scores of one PLINK fileset (variants already restricted to HapMap3), as in ldsc --l2: returns CHR SNP BP CM MAF L2.
#' `coord` is the position used for the window (cM with `wind` = 1, or base pairs with `wind` = 1e6).
ldsc_l2_bfile <- function(prefix, coord = NULL, wind = 1e6, block = 1000L) {
  bim <- utils::read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  n <- length(readLines(paste0(prefix, ".fam"))); nb <- ceiling(n / 4); m <- nrow(bim)
  con <- file(paste0(prefix, ".bed"), "rb"); sig <- readBin(con, "raw", 3); raw <- matrix(readBin(con, "raw", nb * m), nb, m); close(con)
  if (!identical(sig, as.raw(c(0x6c, 0x1b, 0x01)))) stop("Not a variant-major PLINK .bed file: ", prefix)
  # SNP filter of ldsc (__filter_snps_maf__): keep polymorphic variants with fewer than n heterozygous-or-missing people
  keep <- logical(m); frq <- numeric(m)
  for (s in split(seq_len(m), ceiling(seq_len(m) / 5000))) { g <- .bed_cols(raw, s, n)
    f <- colMeans(g, na.rm = TRUE) / 2; hm <- colSums(is.na(g) | g == 1); frq[s] <- f; keep[s] <- pmin(f, 1 - f) > 0 & hm < n }
  bim <- bim[keep, ]; raw <- raw[, keep, drop = FALSE]; frq <- frq[keep]; m <- nrow(bim)
  x <- if (is.null(coord)) bim$BP else coord[keep]
  lo <- findInterval(x - wind, x, left.open = TRUE) + 1L; hi <- findInterval(x + wind, x)   # window: |x_j - x_i| <= wind
  l2 <- numeric(m)
  for (a in split(seq_len(m), ceiling(seq_len(m) / block))) {
    b <- seq.int(lo[a[1]], hi[a[length(a)]]); za <- .ldsc_std(.bed_cols(raw, a, n)); zb <- .ldsc_std(.bed_cols(raw, b, n))
    r2 <- (crossprod(za, zb) / n)^2; u <- r2 - (1 - r2) / (n - 2)
    u[abs(outer(x[a], x[b], "-")) > wind] <- 0; l2[a] <- rowSums(u)
  }
  data.frame(CHR = bim$CHR, SNP = bim$SNP, BP = bim$BP, CM = if (is.null(coord)) bim$CM else x, MAF = pmin(frq, 1 - frq), L2 = l2, stringsAsFactors = FALSE)
}

#' Compute LD scores in LDSC format (<chr>.l2.ldscore.gz, <chr>.l2.M, <chr>.l2.M_5_50) from the 1000 Genomes panel of `pop`.
#' Window 1 cM when a genetic map is given (`cm`: named vector rsID -> cM, e.g. the map positions of the official files), else 1 Mb.
ldsc_ldscores_compute <- function(pop = "EUR", out = ldsc_ld_folder(pop), hm3 = file.path(ldsc_ref_dir(), "w_hm3.snplist"), panel = ld_panel_dir(),
                                  bin = plink_bin(), cm = NULL, chroms = 1:22, progress = function(chr, done) NULL) {
  if (!file.exists(hm3)) stop("The HapMap3 variant list is missing; download it first.")
  if (is.na(bin)) stop("PLINK is not available.")
  bf <- file.path(panel, pop); if (!file.exists(paste0(bf, ".bed"))) stop("The 1000 Genomes panel for ", pop, " is not installed.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE); tmp <- tempfile("ldsc"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  snps <- utils::read.table(hm3, header = TRUE, stringsAsFactors = FALSE)$SNP; ex <- file.path(tmp, "hm3.txt"); writeLines(snps, ex)
  for (k in seq_along(chroms)) { ch <- chroms[k]; pre <- file.path(tmp, paste0("c", ch))
    st <- system2(bin, c("--bfile", shQuote(bf), "--chr", ch, "--extract", shQuote(ex), "--make-bed", "--out", shQuote(pre), "--silent"), stdout = FALSE, stderr = FALSE)
    if (st != 0 || !file.exists(paste0(pre, ".bed"))) stop("PLINK could not extract chromosome ", ch, ".")
    bim <- utils::read.table(paste0(pre, ".bim"), stringsAsFactors = FALSE)
    if (anyDuplicated(bim$V2)) { writeLines(bim$V2[!duplicated(bim$V2) & !duplicated(bim$V2, fromLast = TRUE)], file.path(tmp, "u.txt"))
      system2(bin, c("--bfile", shQuote(pre), "--extract", shQuote(file.path(tmp, "u.txt")), "--make-bed", "--out", shQuote(paste0(pre, "u")), "--silent"), stdout = FALSE, stderr = FALSE); pre <- paste0(pre, "u"); bim <- utils::read.table(paste0(pre, ".bim"), stringsAsFactors = FALSE) }
    coord <- if (!is.null(cm)) unname(cm[bim$V2]) else NULL
    if (!is.null(coord)) { if (anyNA(coord)) { ok <- !is.na(coord); writeLines(bim$V2[ok], file.path(tmp, "m.txt"))
        system2(bin, c("--bfile", shQuote(pre), "--extract", shQuote(file.path(tmp, "m.txt")), "--make-bed", "--out", shQuote(paste0(pre, "m")), "--silent"), stdout = FALSE, stderr = FALSE)
        pre <- paste0(pre, "m"); bim <- utils::read.table(paste0(pre, ".bim"), stringsAsFactors = FALSE); coord <- unname(cm[bim$V2]) }
      if (is.unsorted(coord)) stop("The genetic map is not in the order of the positions on chromosome ", ch, ".") }
    d <- ldsc_l2_bfile(pre, coord, wind = if (is.null(coord)) 1e6 else 1)
    if (ch == 6) d <- d[!(d$BP >= 25e6 & d$BP <= 34e6), ]   # HLA rows left out after computing, as in the official 6.l2 files (other scores unchanged)
    con <- gzfile(file.path(out, paste0(ch, ".l2.ldscore.gz")), "w"); utils::write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
    writeLines(as.character(nrow(d)), file.path(out, paste0(ch, ".l2.M"))); writeLines(as.character(sum(d$MAF > 0.05)), file.path(out, paste0(ch, ".l2.M_5_50")))
    progress(ch, k / length(chroms))
  }
  writeLines(c(sprintf("LD scores computed by MR Studio from the 1000 Genomes Phase 3 panel (%s) with the estimator of ldsc %s.", pop, LDSC_PIN),
               "Sum of the sample-size corrected r2 over HapMap3 variants (w_hm3.snplist), the variant itself included. HLA (chr6 25-34 Mb) left out.",
               if (is.null(cm)) "Window: 1 Mb on each side." else "Window: 1 cM on each side (genetic map positions of the official LD score files).",
               paste("computed:", format(Sys.Date()))), file.path(out, "README"))
  invisible(out)
}

# ---- Summary statistics for LDSC (as munge_sumstats.py with --merge-alleles w_hm3.snplist) -----------------------------------------
# Kept: HapMap3 variants whose two alleles are those of the list (strand-ambiguous A/T and C/G never are), MAF >= 0.01 when the
# allele frequency is known, N >= 90th percentile / 1.5 when N varies by variant. Z = beta / SE, signed for the A1 allele of the list.
PIN_GENOMICSEM <- "da95d431a4709693d07eae55e172ade55526948b"   # GitHub GenomicSEM/GenomicSEM, version 0.0.5
LDSC_MIN_SNPS <- 200000L   # fewer HapMap3 variants than this means the file is not genome-wide (LDSC needs the whole genome)

ldsc_munge <- function(d, hm3 = file.path(ldsc_ref_dir(), "w_hm3.snplist")) {
  h <- if (is.data.frame(hm3)) hm3 else utils::read.table(hm3, header = TRUE, stringsAsFactors = FALSE)
  led <- list(rows = nrow(d))
  if (!("snp" %in% names(d))) stop("LD score regression needs rsIDs; this dataset identifies variants by position only.")
  if (!("n" %in% names(d)) || all(is.na(d$n))) stop("LD score regression needs the sample size; give it when importing the dataset.")
  i <- match(d$snp, h$SNP); d <- d[!is.na(i), , drop = FALSE]; h <- h[i[!is.na(i)], , drop = FALSE]; led$hapmap3 <- nrow(d)
  same <- d$ea == h$A1 & d$oa == h$A2; swap <- d$ea == h$A2 & d$oa == h$A1
  ok <- same | swap; d <- d[ok, , drop = FALSE]; h <- h[ok, , drop = FALSE]; same <- same[ok]; led$alleles_match <- nrow(d)
  if ("eaf" %in% names(d) && any(is.finite(d$eaf))) { k <- !is.finite(d$eaf) | pmin(d$eaf, 1 - d$eaf) >= 0.01; d <- d[k, , drop = FALSE]; h <- h[k, , drop = FALSE]; same <- same[k] }
  led$maf <- nrow(d)
  if (length(unique(stats::na.omit(d$n))) > 1L) { k <- is.finite(d$n) & d$n >= stats::quantile(d$n, 0.9, na.rm = TRUE) / 1.5; d <- d[k, , drop = FALSE]; h <- h[k, , drop = FALSE]; same <- same[k] }
  z <- d$beta / d$se; k <- is.finite(z) & is.finite(d$n)
  out <- data.frame(SNP = h$SNP, A1 = h$A1, A2 = h$A2, N = d$n, Z = ifelse(same, z, -z), stringsAsFactors = FALSE)[k, , drop = FALSE]
  led$kept <- nrow(out); attr(out, "ledger") <- led; out
}

.ldsc_slug <- function(trait) slugify(trait, 60L)

#' Write the script for one pair: LDSC of the exposure and the outcome (heritability of each, genetic correlation, intercepts).
ldsc_script <- function(folder, files, labels, ld, ld_what, sprev, pprev) {
  q <- function(s) encodeString(s, quote = "\""); v <- function(x) paste0("c(", paste(ifelse(is.na(x), "NA", format(x, digits = 15)), collapse = ", "), ")")
  c("# MR Studio: LD score regression (Bulik-Sullivan et al. 2015, Nat Genet 47:291 and 47:1236) with GenomicSEM::ldsc.",
    sprintf("# Package: GenomicSEM commit %s (GitHub GenomicSEM/GenomicSEM).", PIN_GENOMICSEM),
    "# Run from the study folder. The summary statistics in 02_datasets/ldsc were prepared by MR Studio from the Data library files",
    "# (HapMap3 variants, alleles matched to w_hm3.snplist, Z signed for its A1 allele); each has a record of the variants removed.",
    sprintf("PAIR <- %s", q(folder)), sprintf("TRAITS <- c(%s)", paste(vapply(files, q, ""), collapse = ", ")), sprintf("NAMES <- c(%s)", paste(vapply(labels, q, ""), collapse = ", ")),
    sprintf("LD_SCORES <- %s   # %s", q(ld), ld_what),
    sprintf("SAMPLE_PREV <- %s; POPULATION_PREV <- %s   # case proportion and population prevalence of a binary trait (liability scale), NA otherwise", v(sprev), v(pprev)),
    "if (!requireNamespace(\"GenomicSEM\", quietly = TRUE)) stop(\"Package GenomicSEM is not installed (GitHub GenomicSEM/GenomicSEM).\")",
    "od <- file.path(\"04_across_pairs\", \"ldsc\"); dir.create(od, recursive = TRUE, showWarnings = FALSE)",
    "r <- GenomicSEM::ldsc(traits = TRAITS, sample.prev = SAMPLE_PREV, population.prev = POPULATION_PREV, ld = LD_SCORES, wld = LD_SCORES,",
    "                      trait.names = c(\"X\", \"Y\"), stand = TRUE, ldsc.log = file.path(od, PAIR))",
    "lg <- readLines(file.path(od, paste0(PAIR, \"_ldsc.log\")))   # intercept standard errors are printed in the log only (4 digits)",
    "num <- function(pat) { x <- regmatches(lg, regexpr(pat, lg)); m <- regmatches(x, gregexpr(\"-?[0-9.]+(e-?[0-9]+)?\", x)); t(vapply(m, function(y) as.numeric(utils::tail(y, 2)), numeric(2))) }",
    "ic <- num(\"^Intercept: .*\"); cx <- num(\"^Cross trait Intercept: .*\"); mc <- as.numeric(sub(\".*: \", \"\", grep(\"^Mean Chi\\\\^2\", lg, value = TRUE)))",
    "se <- sqrt(diag(r$V)); ses <- sqrt(diag(r$V_Stand)); p <- function(b, s) 2 * stats::pnorm(-abs(b / s))",
    "out <- data.frame(quantity = c(\"h2\", \"h2\", \"rg\", \"genetic_covariance\", \"intercept\", \"intercept\", \"cross_trait_intercept\"),",
    "  trait = c(NAMES, \"both\", \"both\", NAMES, \"both\"),",
    "  estimate = c(r$S[1, 1], r$S[2, 2], r$S_Stand[2, 1], r$S[2, 1], r$I[1, 1], r$I[2, 2], r$I[2, 1]),",
    "  se = c(se[1], se[3], ses[2], se[2], ic[1, 2], ic[2, 2], cx[1, 2]), stringsAsFactors = FALSE)",
    "out$p <- c(p(out$estimate[1:4], out$se[1:4]), p(out$estimate[5:6] - 1, out$se[5:6]), p(out$estimate[7], out$se[7]))   # intercepts tested against 1 (own) and 0 (cross-trait)",
    "out$ci_lower <- out$estimate - 1.96 * out$se; out$ci_upper <- out$estimate + 1.96 * out$se",
    "out$scale <- c(ifelse(is.na(POPULATION_PREV), \"observed\", \"liability\"), \"\", \"\", \"\", \"\", \"\")",
    "out$mean_chi2 <- c(mc[1], mc[length(mc)], NA, NA, mc[1], mc[length(mc)], NA); out$snps <- r$m",
    "utils::write.csv(out, file.path(od, paste0(PAIR, \"_results.csv\")), row.names = FALSE)",
    "out")
}

#' Prepare both datasets of a pair, write and run its LDSC script. `source` = "official" or a super-population of computed scores.
ldsc_pair_run <- function(h, pair, source = NULL, pprev = c(NA_real_, NA_real_), root = studies_root(), dir = ldsc_ref_dir()) {
  if (isTRUE(pair$mv) || is.null(pair$folder)) stop("LD score regression is for a pair with one exposure.")
  ids <- full_gwas_ids(c(pair$ex, pair$out), root); if (!all(startsWith(ids, "lib:"))) stop("LD score regression needs the full summary statistics of both traits, imported into the Data library.")
  rows <- vapply(ids, function(i) tryCatch(as.numeric(library_dataset(sub("^lib:", "", i), root)$recipe$rows_kept %||% NA), error = function(e) NA_real_), 0)
  if (exists("mem_check")) { mc <- mem_check(ldsc_memory_gb(rows), "LD score regression"); if (!mc$ok) stop(mc$msg, call. = FALSE) }   # Round 66
  ds <- lapply(ids, function(i) library_dataset(sub("^lib:", "", i), root))
  anc <- vapply(ds, function(x) toupper(x$recipe$facts$ancestry %||% "EUR"), "")
  if (!identical(anc[1], anc[2])) stop(sprintf("The two GWAS are from different ancestries (%s and %s); LD score regression needs one LD reference for both.", anc[1], anc[2]))
  source <- source %||% if (identical(anc[1], "EUR") && ldsc_ref_status(dir)$official) "official" else anc[1]
  ld <- ldsc_ld_folder(source, dir); if (!.ldsc_complete(ld)) stop(if (identical(source, "official")) "The European LD scores are not installed yet." else sprintf("LD scores for %s have not been computed yet.", source))
  hm3 <- utils::read.table(file.path(dir, "w_hm3.snplist"), header = TRUE, stringsAsFactors = FALSE)
  rel <- file.path("02_datasets", "ldsc"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE); files <- character(2); led <- list()
  for (k in 1:2) { tr <- ds[[k]]$recipe$trait; f <- file.path(rel, paste0(.ldsc_slug(tr), ".sumstats.gz")); files[k] <- f
    if (!file.exists(file.path(h$dir, f))) { m <- ldsc_munge(ds[[k]]$data, hm3); led[[tr]] <- attr(m, "ledger")
      if (nrow(m) < LDSC_MIN_SNPS) stop(sprintf("\"%s\" has %s HapMap3 variants after preparation; LD score regression needs genome-wide summary statistics (at least %s).", tr, format(nrow(m), big.mark = ","), format(LDSC_MIN_SNPS, big.mark = ",")))
      con <- gzfile(file.path(h$dir, f), "w"); utils::write.table(m, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
      atomic_write_text(file.path(h$dir, sub("\\.sumstats\\.gz$", "_preparation.json", f)), jsonlite::toJSON(c(list(trait = tr, dataset_sha256 = ds[[k]]$recipe$original$sha256), attr(m, "ledger")), auto_unbox = TRUE, pretty = TRUE))
      manifest_register(h$dir, c(f, sub("\\.sumstats\\.gz$", "_preparation.json", f)), role = "ldsc_input") } }
  sp <- vapply(ds, function(x) { f <- x$recipe$facts; if (!identical(f$kind, "binary")) return(NA_real_); a <- suppressWarnings(as.numeric(f$ncase)); b <- suppressWarnings(as.numeric(f$ncontrol)); if (length(a) && length(b) && is.finite(a) && is.finite(b)) a / (a + b) else NA_real_ }, 0)
  pprev <- as.numeric(pprev); pprev[is.na(sp)] <- NA_real_; sp[is.na(pprev)] <- NA_real_   # liability scale only with both proportions
  what <- if (identical(source, "official")) sprintf("LDSC authors' European LD scores (Zenodo %s, MD5 of archive %s)", LDSC_REF$w_ld$doi, LDSC_REF$w_ld$md5) else sprintf("computed by MR Studio from the 1000 Genomes Phase 3 %s panel", source)
  rd <- file.path("04_across_pairs", "ldsc"); dir.create(file.path(h$dir, rd), recursive = TRUE, showWarnings = FALSE); sf <- file.path(h$dir, rd, paste0(pair$folder, ".R"))
  writeLines(ldsc_script(pair$folder, files, c(pair$exName, pair$outName), normalizePath(ld, winslash = "/"), what, sp, pprev), sf)
  old <- setwd(h$dir); res <- tryCatch({ utils::capture.output(sys.source(sf, envir = new.env(parent = globalenv()), keep.source = FALSE))
    utils::read.csv(file.path(rd, paste0(pair$folder, "_results.csv")), stringsAsFactors = FALSE) }, finally = setwd(old))
  manifest_register(h$dir, file.path(rd, paste0(pair$folder, c(".R", "_results.csv", "_ldsc.log"))), role = "ldsc")
  attr(res, "labels") <- c(x = pair$exName, y = pair$outName); attr(res, "source") <- what; res
}

#' Results already written for a pair (NULL when LDSC has not been run).
ldsc_pair_read <- function(h, pair) { f <- file.path(h$dir, "04_across_pairs", "ldsc", paste0(pair$folder, "_results.csv"))
  if (!file.exists(f)) return(NULL); r <- utils::read.csv(f, stringsAsFactors = FALSE); attr(r, "labels") <- c(x = pair$exName, y = pair$outName); r }

LDSC_METHODS_TEXT <- "SNP heritability, genetic correlation and the LD score regression intercepts were estimated by LD score regression [Bulik-Sullivan 2015] with GenomicSEM [Grotzinger 2019], on HapMap3 variants with alleles matched to the LDSC reference list; the cross-trait intercept reflects sample overlap and shared confounding between the two GWAS."
ldsc_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); g <- function(q, t = NULL) { k <- res$quantity == q; if (!is.null(t)) k <- k & res$trait == t; res[k, ][1, ] }
  f <- function(r, d = 3) sprintf("%.*f (SE %.*f)", d, r$estimate, d, r$se)
  rg <- g("rg"); ci <- g("cross_trait_intercept")
  sprintf("LD score regression (%s variants): SNP heritability %s for %s and %s for %s; genetic correlation %s (95%% CI %.3f to %.3f; P = %s). Intercepts %s and %s; cross-trait intercept %s%s.",
          format(res$snps[1], big.mark = ","), f(g("h2", l[["x"]])), l[["x"]], f(g("h2", l[["y"]])), l[["y"]], f(rg), rg$ci_lower, rg$ci_upper, format(signif(rg$p, 2)),
          f(g("intercept", l[["x"]])), f(g("intercept", l[["y"]])), f(ci), if (ci$p < 0.05) ", which suggests overlapping samples or shared confounding" else "")
}

# ---- Background computation (separate R process, progress in a file), so the app stays usable for the 10 to 30 minutes it takes ----
.ldsc_job_file <- function(dir) file.path(dir, "compute_status.json")
ldsc_job_read <- function(dir = ldsc_ref_dir()) { f <- .ldsc_job_file(dir); if (!file.exists(f)) return(NULL); tryCatch(jsonlite::fromJSON(f), error = function(e) NULL) }
ldsc_status <- function(dir = ldsc_ref_dir()) c(ldsc_ref_status(dir), list(job = ldsc_job_read(dir), genomicsem = requireNamespace("GenomicSEM", quietly = TRUE)))
#' Map positions (cM) of the official files, used for 1 cM windows when they are installed.
ldsc_official_cm <- function(dir = ldsc_ref_dir()) { d <- ldsc_ld_folder("official", dir); if (!.ldsc_complete(d)) return(NULL)
  unlist(lapply(1:22, function(c) { o <- utils::read.table(gzfile(file.path(d, paste0(c, ".l2.ldscore.gz"))), header = TRUE, stringsAsFactors = FALSE); stats::setNames(o$CM, o$SNP) })) }
ldsc_compute_start <- function(pop = "EUR", dir = ldsc_ref_dir(), engine = getOption("mrgui.engine_dir", "")) {
  j <- ldsc_job_read(dir); if (!is.null(j) && identical(j$state, "running")) stop("LD scores are already being computed.")
  if (exists("mem_check")) { mc <- mem_check(LDSCORE_MEMORY_GB, "Computing LD scores"); if (!mc$ok) stop(mc$msg, call. = FALSE) }   # Round 66
  if (!file.exists(file.path(dir, "w_hm3.snplist"))) ldsc_ref_fetch("hm3", dir)
  job <- file.path(dir, "compute_job.R")
  writeLines(c("a <- commandArgs(TRUE); engine <- a[1]; dir <- a[2]; pop <- a[3]",
    "st <- function(state, pct = 0, chr = NA, msg = '') writeLines(jsonlite::toJSON(list(state = state, pop = pop, pct = pct, chr = chr, msg = msg, time = format(Sys.time(), '%Y-%m-%d %H:%M:%S')), auto_unbox = TRUE), file.path(dir, 'compute_status.json'))",
    "tryCatch({ `%||%` <- function(a, b) if (is.null(a)) b else a; source(file.path(engine, 'ld_core.R')); source(file.path(engine, 'ldsc_core.R')); options(mrgui.engine_dir = engine)",
    "  out <- ldsc_ld_folder(pop, dir); tmp <- paste0(out, '.part'); unlink(tmp, recursive = TRUE)",
    "  ldsc_ldscores_compute(pop, out = tmp, hm3 = file.path(dir, 'w_hm3.snplist'), cm = ldsc_official_cm(dir), progress = function(ch, f) st('running', round(100 * f), ch))",
    "  unlink(out, recursive = TRUE); file.rename(tmp, out); st('done', 100)",
    "}, error = function(e) st('error', msg = conditionMessage(e)))"), job)
  writeLines(jsonlite::toJSON(list(state = "running", pop = pop, pct = 0, chr = NA, msg = "", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), .ldsc_job_file(dir))
  rscript_background(c("--vanilla", shQuote(job), shQuote(engine), shQuote(dir), pop), log = file.path(dir, "compute.log"), vars = list(MRGUI_PLINK = plink_bin(), MRGUI_LD_DIR = ld_panel_dir()))
  invisible(ldsc_status(dir))
}
