# Data import (FR-71 to FR-78): read a GWAS summary file, suggest the column roles, check it, save a standard copy, remember the layout.
# Plain base R on purpose, so the bundled runtime needs nothing extra. Individual-level data is refused (S-06).

IMPORT_ROLES <- c(snp = "Variant ID (rsID)", chr = "Chromosome", pos = "Position", ea = "Effect allele", oa = "Other allele", eaf = "Effect allele frequency",
                  beta = "Effect size", se = "Standard error", p = "P value", n = "Sample size")
IMPORT_REQUIRED <- c("snp", "ea", "oa", "beta", "se")

.norm_col <- function(x) gsub("[^a-z0-9]", "", tolower(x))

IMPORT_SYNONYMS <- list(
  snp = c("snp", "rsid", "rsids", "rs", "markername", "marker", "variantid", "variant", "snpid", "id", "rsnumber", "rsnum"),
  chr = c("chr", "chrom", "chromosome", "chrcode", "chromosomenumber", "hgchr"),
  pos = c("pos", "position", "bp", "basepair", "basepairposition", "bppos", "bpposition", "genpos", "chrpos", "hgpos"),
  ea = c("ea", "effectallele", "a1", "allele1", "alt", "testedallele", "eallele", "effallele", "codedallele", "tested"),
  oa = c("oa", "otherallele", "a2", "allele2", "ref", "nea", "noneffectallele", "reference", "refallele", "nonefallele", "othalle"),
  eaf = c("eaf", "effectallelefreq", "effectallelefrequency", "freq", "frq", "af", "a1freq", "freq1", "afreq", "af1", "effectaf", "a1frq", "maf"),
  beta = c("beta", "b", "effect", "effectsize", "es", "logor", "logodds", "estimate", "betahat", "or", "oddsratio"),
  se = c("se", "stderr", "stderror", "standarderror", "sebeta", "sebetahat", "selogor"),
  p = c("p", "pval", "pvalue", "pvalu", "pvalues", "pvalueassoc", "plrt", "pwald", "pbolt", "pboltlmminf", "pboltlmm", "negativelog10p"),
  n = c("n", "samplesize", "neff", "totaln", "nsamples", "ntotal", "nsample", "weight", "nobs"))

#' Layers 1 and 2 of FR-78: normalised synonyms (high) and fuzzy matches (medium, with a score). Layer 3 looks at the values (low).
import_suggest <- function(cols, rows) {
  nc <- .norm_col(cols); out <- list(); used <- character()
  pick <- function(role, col, conf, why) { if (!(col %in% used)) { out[[role]] <<- list(column = col, confidence = conf, why = why); used <<- c(used, col) } }
  for (role in names(IMPORT_SYNONYMS)) { h <- which(nc %in% IMPORT_SYNONYMS[[role]]); if (length(h)) pick(role, cols[h[1]], "high", "Matches a known heading") }
  for (role in setdiff(names(IMPORT_SYNONYMS), names(out))) {
    free <- setdiff(seq_along(cols), match(used, cols)); if (!length(free)) break
    syn <- IMPORT_SYNONYMS[[role]]; best <- NULL; bs <- 1
    for (i in free) { d <- min(utils::adist(nc[i], syn) / pmax(nchar(nc[i]), nchar(syn))); if (d < bs) { bs <- d; best <- i } }
    if (!is.null(best) && bs <= 0.3 && nchar(nc[best]) >= 3) pick(role, cols[best], "medium", sprintf("Looks similar to a known heading (score %.2f)", 1 - bs))
  }
  vals <- function(col) { v <- vapply(rows, function(r) as.character(r[[col]] %||% ""), ""); v[nzchar(v) & v != "NA"] }
  free <- setdiff(cols, used)
  for (role in setdiff(c("snp", "ea", "oa", "pos", "p", "eaf"), names(out))) {
    for (col in free) {
      v <- vals(col); if (length(v) < 3) next; num <- suppressWarnings(as.numeric(v))
      ok <- switch(role, snp = all(grepl("^rs[0-9]+$", v)), ea = , oa = all(grepl("^[ACGTacgt]$", v)),
                   pos = all(!is.na(num)) && all(num == round(num)) && all(num > 1000), p = all(!is.na(num)) && all(num >= 0 & num <= 1) && any(num < 0.01) && stats::median(num) > 0.1,
                   eaf = all(!is.na(num)) && all(num > 0 & num < 1), FALSE)
      if (ok) { pick(role, col, "low", "Guessed from the values, please confirm"); free <- setdiff(free, col); break }
    }
  }
  out
}

.detect_sep <- function(line) { cands <- c("\t", ",", ";", " "); cnt <- vapply(cands, function(s) lengths(regmatches(line, gregexpr(s, line, fixed = TRUE))), 0L); cands[which.max(cnt)] }
.read_tbl <- function(path, ...) { con <- .open_any(path); on.exit(close(con)); utils::read.table(con, ...) }
.open_any <- function(path) if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")

#' The start of a file: separator, headings and the first rows. Refuses formats and layouts that are not supported.
import_read_head <- function(path, name = basename(path), nrows = 8L) {
  if (grepl("\\.(zip|bgz|vcf)(\\.gz)?$", name, ignore.case = TRUE)) stop("This file type is not supported yet. Supported: text tables (TSV, CSV, space separated), optionally gzip-compressed (.gz).")
  con <- .open_any(path); on.exit(close(con))
  lines <- readLines(con, n = nrows + 1L, warn = FALSE)
  if (!length(lines)) stop("The file is empty.")
  sep <- .detect_sep(lines[1])
  split <- function(l) if (sep == " ") strsplit(trimws(l), "[ ]+")[[1]] else strsplit(l, sep, fixed = TRUE)[[1]]
  cols <- gsub('^"|"$', "", split(lines[1])); cols[!nzchar(cols)] <- paste0("column_", seq_along(cols))[!nzchar(cols)]
  if (anyDuplicated(cols)) stop("Two columns have the same heading: ", paste(unique(cols[duplicated(cols)]), collapse = ", "), ".")
  if (length(cols) < 4) stop("The file has fewer than 4 columns. A GWAS summary file needs at least a variant, alleles, an effect and its standard error.")
  nm <- .norm_col(cols)
  if (any(nm %in% c("iid", "fid", "individualid", "sampleid")) || sum(grepl("^(gt|genotype|geno)", nm)) >= 2)
    stop("This looks like individual-level data (one row per person). MR Studio only accepts summary statistics, so nothing was read. Individual-level data stays out of the app by design.")
  rows <- lapply(lines[-1], function(l) { v <- gsub('^"|"$', "", split(l)); length(v) <- length(cols); as.list(setNames(v, cols)) })
  list(name = name, size = file.size(path), sep = sep, columns = cols, rows = rows, fingerprint = sort(unique(nm)))
}

#' Read the mapped columns completely and report what is wrong. errors block the import; warnings are passed on to the user.
import_validate <- function(path, sep, columns, mapping, facts, max_rows = 5e6) {
  err <- character(); warn <- character(); mp <- mapping[nzchar(unlist(mapping))]
  miss <- setdiff(IMPORT_REQUIRED, names(mp))
  if (!("snp" %in% names(mp)) && all(c("chr", "pos") %in% names(mp))) miss <- setdiff(miss, "snp")
  for (m in miss) err <- c(err, sprintf("No column chosen for \"%s\". It is required.", IMPORT_ROLES[[m]]))
  if (!all(unlist(mp) %in% columns)) err <- c(err, "A chosen column is not in the file.")
  if (anyDuplicated(unlist(mp))) err <- c(err, "The same column is used for two different roles.")
  if (length(err)) return(list(ok = FALSE, errors = err, warnings = warn, stats = list()))
  cc <- ifelse(columns %in% unlist(mp), "character", "NULL"); names(cc) <- columns
  d <- .read_tbl(path, header = TRUE, sep = if (sep == " ") "" else sep, quote = "\"", comment.char = "", colClasses = cc, nrows = max_rows, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  nr <- nrow(d); rn <- setNames(names(mp), unlist(mp)); names(d) <- unname(rn[names(d)])
  num <- function(x) suppressWarnings(as.numeric(x))
  for (r in intersect(c("eaf", "beta", "se", "p", "n", "pos"), names(d))) { bad <- sum(is.na(num(d[[r]])) & !is.na(d[[r]]) & nzchar(d[[r]]) & !(toupper(d[[r]]) %in% c("NA", "NAN", "."))); if (bad) warn <- c(warn, sprintf("%s values in \"%s\" are not numbers and will be dropped.", format(bad, big.mark = ","), IMPORT_ROLES[[r]])); d[[r]] <- num(d[[r]]) }
  if (nr >= max_rows) warn <- c(warn, sprintf("Only the first %s rows were read.", format(max_rows, big.mark = ",")))
  bad_se <- sum(d$se <= 0 | is.na(d$se)); if (bad_se) warn <- c(warn, sprintf("%s variants have a missing or non-positive standard error and will be dropped.", format(bad_se, big.mark = ",")))
  if ("p" %in% names(d)) { bp <- sum(d$p < 0 | d$p > 1, na.rm = TRUE); if (bp) warn <- c(warn, sprintf("%s P values are outside 0 to 1.", format(bp, big.mark = ","))) }
  ok_al <- grepl("^[ACGT]$", toupper(d$ea)) & grepl("^[ACGT]$", toupper(d$oa)); nin <- sum(!ok_al)
  if (nin) warn <- c(warn, sprintf("%s variants are not single-letter SNPs (for example insertions or deletions) and will be dropped.", format(nin, big.mark = ",")))
  if (mean(!ok_al) > 0.5) err <- c(err, "More than half of the rows do not have single-letter alleles in the chosen allele columns. Check that the effect allele and other allele columns are not swapped with another column.")
  pal <- sum(ok_al & ((toupper(d$ea) == "A" & toupper(d$oa) == "T") | (toupper(d$ea) == "T" & toupper(d$oa) == "A") | (toupper(d$ea) == "C" & toupper(d$oa) == "G") | (toupper(d$ea) == "G" & toupper(d$oa) == "C")))
  key <- if ("snp" %in% names(d)) d$snp else paste(d$chr, d$pos, d$ea, d$oa); dup <- sum(duplicated(key)); if (dup) warn <- c(warn, sprintf("%s rows repeat a variant. Only the first of each is kept.", format(dup, big.mark = ",")))
  if ("eaf" %in% names(d)) {
    f <- d$eaf[!is.na(d$eaf)]
    if (length(f) && any(f < 0 | f > 1)) err <- c(err, "Frequencies outside 0 to 1 were found in the frequency column.")
    if (length(f) > 100 && max(f) <= 0.5 + 1e-9 && identical(facts$frq, "eaf")) warn <- c(warn, "All frequencies are 0.5 or less, which suggests a minor allele frequency, but you said it is the effect allele frequency. A minor allele frequency cannot be used as effect allele frequency. Please check.")
    if (identical(facts$frq, "maf")) warn <- c(warn, "You said the frequency column is a minor allele frequency. It will not be used as an effect allele frequency.")
  }
  if (pal > 0 && !("eaf" %in% names(d) && identical(facts$frq, "eaf"))) warn <- c(warn, sprintf("%s palindromic variants (A/T or C/G) and no usable effect allele frequency. %s", format(pal, big.mark = ","), if (identical(facts$strand, "forward")) "You declared the forward strand: they will be kept when the other file of a pair is declared forward strand too, and left out otherwise." else "Palindromic variants cannot be aligned to another dataset by alleles alone, so they will be left out when the datasets are harmonised, unless both files are declared to be on the forward strand."))
  one <- function(x) { v <- if (is.null(x)) NA_real_ else num(x); if (length(v) == 1L && is.finite(v)) v else NA_real_ }
  cases_ok <- isTRUE(one(facts$ncase) >= 1 && one(facts$ncontrol) >= 1)
  given <- nzchar(paste0(if (is.null(facts$ncase)) "" else facts$ncase, if (is.null(facts$ncontrol)) "" else facts$ncontrol))
  if (identical(facts$kind, "binary") && given && !cases_ok) err <- c(err, "The numbers of cases and controls must both be whole numbers of 1 or more.")
  if (!("n" %in% names(d)) && !(is.finite(num(facts$n)) && num(facts$n) > 0) && !(identical(facts$kind, "binary") && cases_ok)) warn <- c(warn, "No sample size column and no total sample size given. Methods that need N (for example the Steiger test) will not be offered.")
  if (identical(facts$kind, "binary") && !cases_ok) warn <- c(warn, "Binary trait without the numbers of cases and controls. The Steiger direction test cannot use this dataset until they are known (import it again with them).")
  if (identical(facts$scope, "genomewide") && nr < 100000) warn <- c(warn, sprintf("You said the file is genome-wide but it has only %s rows. If it holds only selected variants, change that answer.", format(nr, big.mark = ",")))
  if (identical(facts$scale, "z")) err <- c(err, "Z-scores are not supported yet. Please supply the effect size and its standard error.")
  if (identical(facts$independent, "clump") && !identical(facts$role, "outcome")) {
    if (!("snp" %in% names(d)) || !any(grepl("^rs[0-9]+$", d$snp))) err <- c(err, "Clumping on this computer needs rsIDs (for example rs12345), and this file has none in the chosen column.")
    r2 <- suppressWarnings(as.numeric(facts$clump_r2 %||% 0.001)); kb <- suppressWarnings(as.numeric(facts$clump_kb %||% 10000))
    if (length(r2) != 1L || !is.finite(r2) || r2 <= 0 || r2 >= 1) err <- c(err, "The LD threshold r2 must be between 0 and 1.")
    if (length(kb) != 1L || !is.finite(kb) || kb < 1 || kb > 100000) err <- c(err, "The clumping window must be between 1 and 100,000 kb.")
    if (!(toupper(facts$clump_pop %||% facts$ancestry %||% "") %in% c("EUR", "EAS", "AFR", "AMR", "SAS"))) err <- c(err, "Choose the 1000 Genomes ancestry to clump with (European, East Asian, African, Admixed American or South Asian).")
  }
  if (identical(facts$scale, "or") && any(d$beta <= 0, na.rm = TRUE)) err <- c(err, "Odds ratios must be positive, but zero or negative values were found. If the column holds log odds ratios, choose \"Beta or log odds\" instead.")
  keep <- ok_al & !is.na(d$beta) & !is.na(d$se) & d$se > 0 & !duplicated(key)
  list(ok = !length(err), errors = err, warnings = warn, stats = list(rows = nr, kept = sum(keep), palindromic = pal, snps = if ("snp" %in% names(d)) sum(grepl("^rs[0-9]+$", d$snp)) else 0L, has_eaf = "eaf" %in% names(d) && identical(facts$frq, "eaf")))
}

#' Where imported datasets, with their recipe, live.
library_dir <- function(root = studies_root()) file.path(root, "Data library")

#' Write the standard copy and the recipe (file name, size, sha256, mapping, facts: FR-75). The original file is never changed.
import_save <- function(path, orig_name, sep, columns, mapping, facts, root = studies_root()) {
  v <- import_validate(path, sep, columns, mapping, facts); if (!v$ok) stop(paste(v$errors, collapse = " "))
  mp <- mapping[nzchar(unlist(mapping))]
  cc <- ifelse(columns %in% unlist(mp), "character", "NULL"); names(cc) <- columns
  d <- .read_tbl(path, header = TRUE, sep = if (sep == " ") "" else sep, quote = "\"", comment.char = "", colClasses = cc, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  rn <- setNames(names(mp), unlist(mp)); names(d) <- unname(rn[names(d)])
  for (r in intersect(c("eaf", "beta", "se", "p", "n", "pos"), names(d))) d[[r]] <- suppressWarnings(as.numeric(d[[r]]))
  d$ea <- toupper(d$ea); d$oa <- toupper(d$oa)
  # "No, the effect size is for the other allele": the labels swap, so the effect allele column always holds the allele the effect belongs to.
  if (identical(facts$ea_ok, "swap")) { tmp <- d$ea; d$ea <- d$oa; d$oa <- tmp }
  key <- if ("snp" %in% names(d)) d$snp else paste(d$chr, d$pos, d$ea, d$oa)
  keep <- grepl("^[ACGT]$", d$ea) & grepl("^[ACGT]$", d$oa) & !is.na(d$beta) & !is.na(d$se) & d$se > 0 & !duplicated(key)
  d <- d[keep, , drop = FALSE]
  if (identical(facts$scale, "or")) { if (any(d$beta <= 0)) stop("Odds ratios must be positive."); d$beta <- log(d$beta) }
  if (!identical(facts$frq, "eaf")) d$eaf <- NULL
  if (!("n" %in% names(d)) && is.finite(suppressWarnings(as.numeric(facts$n)))) d$n <- as.numeric(facts$n)
  nc <- suppressWarnings(as.numeric(if (is.null(facts$ncase)) NA else facts$ncase)); nn <- suppressWarnings(as.numeric(if (is.null(facts$ncontrol)) NA else facts$ncontrol))
  if (!("n" %in% names(d)) && identical(facts$kind, "binary") && length(nc) == 1L && length(nn) == 1L && is.finite(nc) && is.finite(nn)) d$n <- nc + nn   # total N of a case-control study
  nm <- clean_study_name(facts$name); if (!nzchar(nm)) stop("Give the trait a name.")
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) stop("A dataset named \"", nm, "\" already exists in the Data library. Choose another name.")
  p <- study_path_problem(library_dir(root)); if (!is.null(p)) stop(p)
  dir.create(dest, recursive = TRUE)
  out <- file.path(dest, "standardised.csv"); utils::write.csv(d, out, row.names = FALSE, na = "")
  rec <- list(trait = nm, role = facts$role, facts = facts, mapping = mp, original = list(file = orig_name, bytes = file.size(path), sha256 = sha256_file(path)),
              rows_in_file = v$stats$rows, rows_kept = nrow(d), warnings = v$warnings, saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              usable_as_exposure = isTRUE(facts$role %in% c("exposure", "either")) && (identical(facts$independent, "yes") || identical(facts$independent, "clump")),
              usable_as_outcome = facts$role %in% c("outcome", "either"), usable_in_runs = TRUE)
  atomic_write_text(file.path(dest, "recipe.json"), jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  rec
}

library_list <- function(root = studies_root()) {
  ld <- library_dir(root); if (!dir.exists(ld)) return(list())
  d <- list.dirs(ld, recursive = FALSE, full.names = TRUE); d <- d[file.exists(file.path(d, "recipe.json"))]
  recs <- lapply(d, function(x) tryCatch(jsonlite::fromJSON(file.path(x, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL))
  Filter(Negate(is.null), recs)
}

# ---- column-mapping templates (FR-77) ------------------------------------------------------------------------------------------
templates_dir <- function() file.path(settings_dir(), "templates")
template_list <- function() {
  d <- templates_dir(); if (!dir.exists(d)) return(list())
  f <- list.files(d, pattern = "\\.json$", full.names = TRUE)
  Filter(Negate(is.null), lapply(f, function(x) tryCatch(jsonlite::fromJSON(x, simplifyVector = FALSE), error = function(e) NULL)))
}
template_save <- function(name, fingerprint, mapping, facts) {
  nm <- clean_study_name(name); if (!nzchar(nm)) stop("Give the template a name.")
  dir.create(templates_dir(), recursive = TRUE, showWarnings = FALSE)
  facts$name <- NULL
  atomic_write_text(file.path(templates_dir(), paste0(nm, ".json")), jsonlite::toJSON(list(name = nm, fingerprint = as.list(fingerprint), mapping = mapping[nzchar(unlist(mapping))], facts = facts), auto_unbox = TRUE, pretty = TRUE))
  nm
}
template_delete <- function(name) { f <- file.path(templates_dir(), paste0(clean_study_name(name), ".json")); if (file.exists(f)) unlink(f); invisible(TRUE) }
#' The saved template whose heading set is the same or nearly the same (Jaccard 0.8 or more). NULL when none fits.
template_match <- function(fingerprint) {
  best <- NULL; bs <- 0
  for (t in template_list()) { fp <- unlist(t$fingerprint); j <- length(intersect(fp, fingerprint)) / length(union(fp, fingerprint)); if (j > bs) { bs <- j; best <- t } }
  if (!is.null(best) && bs >= 0.8) c(best, list(similarity = bs)) else NULL
}

#' One imported dataset, read back from the Data library: its recipe and its standard copy.
library_dataset <- function(trait, root = studies_root()) {
  dir <- file.path(library_dir(root), basename(trait))
  rf <- file.path(dir, "recipe.json"); cf <- file.path(dir, "standardised.csv")
  if (!file.exists(rf) || !file.exists(cf)) stop("The dataset \"", trait, "\" is not in the Data library any more.")
  rec <- jsonlite::fromJSON(rf, simplifyVector = FALSE)
  hd <- names(utils::read.csv(cf, nrows = 1, check.names = FALSE)); cc <- c(snp = "character", chr = "character", ea = "character", oa = "character")
  d <- utils::read.csv(cf, stringsAsFactors = FALSE, colClasses = cc[intersect(names(cc), hd)])
  list(recipe = rec, data = d, dir = dir)
}
