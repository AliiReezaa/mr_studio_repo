# Data import (FR-71 to FR-78): read a GWAS summary file, suggest the column roles, check it, save a standard copy, remember the layout.
# Plain base R on purpose, so the bundled runtime needs nothing extra. Individual-level data is refused (S-06).

IMPORT_ROLES <- c(snp = "Variant ID (rsID)", chr = "Chromosome", pos = "Position", ea = "Effect allele", oa = "Other allele", eaf = "Effect allele frequency",
                  beta = "Effect size", se = "Standard error", p = "P value", n = "Sample size", z = "Z score")
IMPORT_REQUIRED <- c("snp", "ea", "oa", "beta", "se")   # Round 77: the effect or the SE may be derived (see import_missing)

.norm_col <- function(x) gsub("[^a-z0-9]", "", tolower(x))

IMPORT_SYNONYMS <- list(
  snp = c("snp", "rsid", "rsids", "rs", "markername", "marker", "variantid", "variant", "snpid", "id", "rsnumber", "rsnum"),
  chr = c("chr", "chrom", "chromosome", "chrcode", "chromosomenumber", "hgchr"),
  pos = c("pos", "position", "bp", "basepair", "basepairposition", "bppos", "bpposition", "genpos", "chrpos", "hgpos"),
  ea = c("ea", "effectallele", "a1", "allele1", "alt", "testedallele", "eallele", "effallele", "codedallele", "tested"),
  oa = c("oa", "otherallele", "a2", "allele2", "ref", "nea", "noneffectallele", "reference", "refallele", "nonefallele", "othalle"),
  eaf = c("eaf", "effectallelefreq", "effectallelefrequency", "freq", "frq", "af", "a1freq", "freq1", "afreq", "af1", "effectaf", "a1frq",
          "maf", "minorallelefreq", "minorallelefrequency"),   # Round 77: a MAF heading is matched here but marked as a minor allele frequency (import_hints)
  beta = c("beta", "b", "effect", "effectsize", "es", "logor", "logodds", "estimate", "betahat", "or", "oddsratio"),
  se = c("se", "stderr", "stderror", "standarderror", "sebeta", "sebetahat", "selogor"),
  p = c("p", "pval", "pvalue", "pvalu", "pvalues", "pvalueassoc", "plrt", "pwald", "pbolt", "pboltlmminf", "pboltlmm",
        "negativelog10p", "neglog10p", "mlog10p", "minuslog10p", "log10p", "neglogp", "mlogp", "logp", "lp", "neglog10pval", "mlog10pval", "log10pvalue", "neglog10pvalue"),   # Round 77: -log10(P) headings are converted (import_hints)
  n = c("n", "samplesize", "neff", "totaln", "nsamples", "ntotal", "nsample", "weight", "nobs"),
  z = c("z", "zscore", "zstat", "zstatistic", "zvalue", "zval"))

#' Layers 1 and 2 of FR-78: normalised synonyms (high) and fuzzy matches (medium, with a score). Layer 3 looks at the values (low).
import_suggest <- function(cols, rows) {
  nc <- sub("(hg19|hg38|grch37|grch38|b37|b38|build37|build38)$", "", .norm_col(cols)); out <- list(); used <- character()   # Round 77: POS_hg38 is a position
  pick <- function(role, col, conf, why) { if (!(col %in% used)) { out[[role]] <<- list(column = col, confidence = conf, why = why); used <<- c(used, col) } }
  # Round 77: when several headings match, the earliest synonym wins (a plain P before a -log10 P, an EAF before a MAF)
  for (role in names(IMPORT_SYNONYMS)) { h <- which(nc %in% IMPORT_SYNONYMS[[role]]); h <- h[!(cols[h] %in% used)]; if (length(h)) pick(role, cols[h[order(match(nc[h], IMPORT_SYNONYMS[[role]]))[1]]], "high", "Matches a known heading") }
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

IMPORT_MAF_HEADS <- c("maf", "minorallelefreq", "minorallelefrequency")
IMPORT_MLOG10_HEADS <- c("negativelog10p", "neglog10p", "mlog10p", "minuslog10p", "log10p", "neglogp", "mlogp", "logp", "lp", "neglog10pval", "mlog10pval", "log10pvalue", "neglog10pvalue")

#' Round 77: what the headings and first rows say about the answers on the facts step (each one is shown and can be changed).
import_hints <- function(cols, rows, sug) {
  hc <- function(r) if (is.null(sug[[r]])) "" else .norm_col(sug[[r]]$column)
  vals <- function(col) suppressWarnings(as.numeric(vapply(rows, function(r) as.character(r[[col]] %||% NA), "")))
  out <- list()
  if (!is.null(sug$eaf)) out$frq <- if (hc("eaf") %in% IMPORT_MAF_HEADS) "maf" else "eaf"
  if (!is.null(sug$p)) { v <- vals(sug$p$column); out$p_scale <- if (hc("p") %in% IMPORT_MLOG10_HEADS || any(v > 1, na.rm = TRUE)) "mlog10" else "p" }
  if (!is.null(sug$beta) && hc("beta") %in% c("or", "oddsratio")) out$scale <- "or"
  if (is.null(sug$beta) && !is.null(sug$z)) out$scale <- "z"
  all_h <- tolower(paste(cols, collapse = " "))
  if (grepl("hg38|grch38|b38|build38", all_h)) out$build <- "38" else if (grepl("hg19|grch37|b37|build37", all_h)) out$build <- "37"
  out
}

#' Round 77: files that need a cleaned copy before reading: Excel workbooks (.xlsx, .xls), text files with comment or note lines before the
#' heading, a byte-order mark, or decimal commas. The copy is a tab-separated file in the session's temporary folder; the user's file is
#' not changed. Returns the path to read and what was done.
import_prepare <- function(path, name = basename(path), sheet = NULL) {
  if (grepl("\\.xlsx?$", name, ignore.case = TRUE)) return(.prepare_excel(path, name, sheet))
  con <- .open_any(path); on.exit(close(con)); top <- readLines(con, n = 60L, warn = FALSE)
  if (!length(top)) stop("The file is empty.")
  notes <- character(); bom <- startsWith(top[1], "\ufeff") || grepl("^\xef\xbb\xbf", top[1], useBytes = TRUE)
  if (bom) top[1] <- sub("^\ufeff", "", sub("^\xef\xbb\xbf", "", top[1], useBytes = TRUE))
  skip <- 0L; while (skip < length(top) && (grepl("^##", top[skip + 1L]) || !nzchar(trimws(top[skip + 1L])) || (grepl("^#", top[skip + 1L]) && skip + 2L <= length(top) && grepl("^#", top[skip + 2L])))) skip <- skip + 1L
  nf <- function(l) { sp <- .detect_sep(l); length(if (sp == " ") strsplit(trimws(l), "[ ]+")[[1]] else strsplit(l, sp, fixed = TRUE)[[1]]) }
  hash_head <- FALSE
  if (skip < length(top) && grepl("^#", top[skip + 1L])) {   # "#CHROM POS ..." is a heading; "# a note" before an ordinary heading is skipped
    h1 <- sub("^#", "", top[skip + 1L]); nxt <- if (skip + 2L <= length(top)) top[skip + 2L] else ""
    if (nf(h1) >= 4 && nf(h1) == nf(nxt)) hash_head <- TRUE else { skip <- skip + 1L; while (skip < length(top) && !nzchar(trimws(top[skip + 1L]))) skip <- skip + 1L }
  }
  rest <- top[-seq_len(skip + 1L)]; sep <- .detect_sep(sub("^#", "", top[skip + 1L]))
  dec_comma <- sep != "," && length(rest) && mean(grepl("(^|[\t;]) *-?[0-9]+,[0-9]+([eE][-+]?[0-9]+)? *($|[\t;])", rest)) > 0.5
  if (!bom && !skip && !hash_head && !dec_comma) return(list(path = path, notes = notes))
  if (bom) notes <- c(notes, "A byte-order mark at the start of the file was removed.")
  if (skip) notes <- c(notes, sprintf("%d comment or empty line%s before the heading %s skipped.", skip, if (skip > 1) "s" else "", if (skip > 1) "were" else "was"))
  if (hash_head) notes <- c(notes, "The heading line starts with #; the # was removed.")
  if (dec_comma) notes <- c(notes, "Numbers use a decimal comma (for example 0,05); they were read as decimal points.")
  out <- tempfile("import-", fileext = ".tsv"); oc <- file(out, "w"); on.exit(close(oc), add = TRUE)
  close(con); con <- .open_any(path); first <- TRUE; n_skip <- skip
  repeat {
    L <- readLines(con, n = 200000L, warn = FALSE); if (!length(L)) break
    if (n_skip) { k <- min(n_skip, length(L)); L <- L[-seq_len(k)]; n_skip <- n_skip - k; if (!length(L)) next }
    if (first) { L[1] <- sub("^\ufeff", "", sub("^\xef\xbb\xbf", "", L[1], useBytes = TRUE)); if (hash_head) L[1] <- sub("^#", "", L[1]); first <- FALSE }
    if (sep != "\t") L <- if (sep == " ") gsub("[ ]+", "\t", trimws(L)) else gsub(sep, "\t", L, fixed = TRUE)
    if (dec_comma) L <- gsub("(^|\t) *(-?[0-9]+),([0-9]+([eE][-+]?[0-9]+)?) *(?=\t|$)", "\\1\\2.\\3", L, perl = TRUE)
    writeLines(L, oc)
  }
  list(path = out, notes = notes)
}

.prepare_excel <- function(path, name, sheet = NULL) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Reading Excel files needs the readxl package, which is not in this installation.")
  sheets <- readxl::excel_sheets(path); if (!length(sheets)) stop("The workbook has no sheets.")
  if (is.null(sheet) || !(sheet %in% sheets)) {   # the first sheet with a table of at least 4 columns
    sheet <- sheets[1]; for (sh in sheets) { x <- tryCatch(readxl::read_excel(path, sh, col_names = FALSE, col_types = "text", n_max = 30, .name_repair = "minimal"), error = function(e) NULL); if (!is.null(x) && ncol(x) >= 4) { sheet <- sh; break } }
  }
  x <- as.data.frame(readxl::read_excel(path, sheet, col_names = FALSE, col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
  if (!nrow(x)) stop(sprintf("The sheet \"%s\" is empty.", sheet))
  filled <- rowSums(!is.na(x) & nzchar(trimws(as.matrix(x))))
  numish <- function(r) { v <- trimws(unlist(r)); v <- v[!is.na(v) & nzchar(v)]; if (!length(v)) 0 else mean(!is.na(suppressWarnings(as.numeric(v)))) }
  hr <- which(filled >= 4 & vapply(seq_len(nrow(x)), function(i) numish(x[i, ]) < 0.5, TRUE))[1]   # the heading: 4 or more cells, mostly words
  if (is.na(hr)) stop(sprintf("No heading row with at least 4 columns was found on the sheet \"%s\".", sheet))
  keep_c <- which(!is.na(unlist(x[hr, ])) & nzchar(trimws(unlist(x[hr, ]))))
  body <- x[seq_len(nrow(x)) > hr, keep_c, drop = FALSE]; head <- trimws(unlist(x[hr, keep_c]))
  full <- rowSums(!is.na(body) & nzchar(trimws(as.matrix(body)))) >= max(2, ceiling(length(keep_c) / 2))   # notes and footnotes are dropped
  notes <- c(if (length(sheets) > 1) sprintf("The workbook has %d sheets; read: \"%s\".", length(sheets), sheet),
             if (hr > 1) sprintf("%d line%s above the heading (titles or notes) %s skipped.", hr - 1, if (hr > 2) "s" else "", if (hr > 2) "were" else "was"),
             if (sum(!full)) sprintf("%d row%s with mostly empty cells (notes or footnotes) %s left out.", sum(!full), if (sum(!full) > 1) "s" else "", if (sum(!full) > 1) "were" else "was"))
  body <- body[full, , drop = FALSE]; body[is.na(body)] <- ""
  for (j in seq_along(body)) body[[j]] <- gsub("[\t\r\n]+", " ", body[[j]])
  out <- tempfile("import-", fileext = ".tsv")
  utils::write.table(rbind(head, as.matrix(body)), out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  list(path = out, notes = notes, sheets = as.list(sheets), sheet = sheet)
}

.detect_sep <- function(line) { cands <- c("\t", ",", ";", " "); cnt <- vapply(cands, function(s) lengths(regmatches(line, gregexpr(s, line, fixed = TRUE))), 0L); cands[which.max(cnt)] }
.read_tbl <- function(path, ...) { con <- .open_any(path); on.exit(close(con)); utils::read.table(con, ...) }
.open_any <- function(path) if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")

#' The start of a file: separator, headings and the first rows. Refuses formats and layouts that are not supported.
import_read_head <- function(path, name = basename(path), nrows = 8L, sheet = NULL) {
  if (grepl("\\.(zip|bgz|vcf)(\\.gz)?$", name, ignore.case = TRUE)) stop("This file type is not supported here. Supported: text tables (TSV, CSV, space or semicolon separated, optionally gzip-compressed .gz) and Excel workbooks (.xlsx, .xls); GWAS-VCF files have their own import.")
  orig_size <- file.size(path); pr <- import_prepare(path, name, sheet); path <- pr$path   # Round 77: Excel, comment lines, byte-order mark, decimal commas
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
  list(name = name, size = orig_size, sep = sep, columns = cols, rows = rows, fingerprint = sort(unique(nm)), path = path, notes = as.list(pr$notes), sheets = pr$sheets, sheet = pr$sheet)
}

#' Round 77: what the file must have, given the answers. The effect and its SE can be given, or derived: the SE from the effect and its P
#' value, or both from a Z score with the sample size and an allele frequency. Positions without rsIDs need a known genome build.
import_missing <- function(mp, facts) {
  has <- function(r) r %in% names(mp); err <- character()
  z <- identical(facts$scale, "z")
  if (!has("snp") && !(has("chr") && has("pos"))) err <- c(err, "Choose the variant ID (rsID) column, or the chromosome and position columns.")
  if (!has("snp") && has("chr") && has("pos") && !(facts$build %||% "unknown") %in% c("37", "38")) err <- c(err, "The file has positions but no rsIDs, so the genome build must be known: choose GRCh37 or GRCh38.")
  if (identical(facts$build, "38") && !(has("chr") && has("pos"))) err <- c(err, "To convert GRCh38 positions to GRCh37, choose the chromosome and position columns.")
  for (r in c("ea", "oa")) if (!has(r)) err <- c(err, sprintf("No column chosen for \"%s\". It is required.", IMPORT_ROLES[[r]]))
  if (z) {
    if (!has("z")) err <- c(err, "The effect size is a Z score: choose the Z score column.")
    if (identical(facts$kind, "binary")) err <- c(err, "Effect sizes cannot be derived from Z scores for a binary trait (they would not be log odds ratios). Please supply the log odds ratios or odds ratios and their standard errors.")
    nn <- suppressWarnings(as.numeric(facts$n %||% NA))
    if (!has("n") && !(length(nn) == 1L && is.finite(nn) && nn > 0)) err <- c(err, "Deriving the effect size from a Z score needs the sample size: choose the sample size column or enter the total sample size.")
    if (!has("eaf") && !identical(facts$eaf_fill, "panel")) err <- c(err, "Deriving the effect size from a Z score needs an allele frequency: choose the frequency column, or take the frequencies from the 1000 Genomes panel.")
  } else {
    if (!has("beta")) err <- c(err, sprintf("No column chosen for \"%s\". It is required (or choose Z score as the effect size).", IMPORT_ROLES[["beta"]]))
    if (!has("se") && !has("p")) err <- c(err, "No standard error column. Choose it, or choose the P value column so that the standard error can be derived from the effect and its P value.")
  }
  if (identical(facts$eaf_fill, "panel")) {
    pop <- toupper(facts$ancestry %||% "")
    if (!pop %in% c("EUR", "EAS", "AFR", "AMR", "SAS")) err <- c(err, "Frequencies from the 1000 Genomes panel need one ancestry (European, East Asian, African, Admixed American or South Asian).")
    else if (!has("snp")) err <- c(err, "Frequencies from the 1000 Genomes panel are matched by rsID, and the file has no rsID column.")
    else if (!.panel_ready(pop)) err <- c(err, sprintf("The 1000 Genomes %s panel is not installed. Install it in Settings (LD reference panel) first.", pop))
  }
  err
}
.panel_ready <- function(pop) isTRUE(tryCatch(all(file.exists(file.path(ld_panel_dir(), paste0(pop, c(".bed", ".bim", ".fam"))))) && { b <- plink_bin(); length(b) == 1L && !is.na(b) && nzchar(b) }, error = function(e) FALSE))

IMPORT_NA_TOKENS <- c("", "NA", "NAN", ".", "-", "NULL", "NONE", "N/A", "INF", "-INF", "+INF")

.import_read <- function(path, sep, columns, mp, max_rows = Inf) {
  cc <- ifelse(columns %in% unlist(mp), "character", "NULL"); names(cc) <- columns
  d <- .read_tbl(path, header = TRUE, sep = if (sep == " ") "" else sep, quote = "\"", comment.char = "", colClasses = cc, nrows = if (is.finite(max_rows)) max_rows else -1L, fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
  rn <- setNames(names(mp), unlist(mp)); names(d) <- unname(rn[names(d)]); d
}

#' Round 77: every change made to the values, in order, each counted (shown on the check step and kept in the recipe). The same function
#' runs for the check and for the save, so what the user approves is what is saved.
.import_transform <- function(d, facts, final = FALSE, ds_dir = NULL) {
  ch <- character(); warn <- character(); err <- character(); fm <- function(x) format(x, big.mark = ","); vv <- function(x) paste(fm(x), if (x == 1) "variant" else "variants")
  for (r in intersect(c("eaf", "beta", "se", "p", "n", "pos", "z"), names(d))) {
    raw <- trimws(d[[r]]); tok <- toupper(raw) %in% IMPORT_NA_TOKENS; v <- suppressWarnings(as.numeric(raw))
    bad <- sum(is.na(v) & !tok); if (bad) warn <- c(warn, sprintf("%s values in \"%s\" are not numbers and will be dropped.", fm(bad), IMPORT_ROLES[[r]]))
    v[!is.finite(v)] <- NA_real_; d[[r]] <- v
  }
  if ("chr" %in% names(d)) { c0 <- d$chr; d$chr <- .chr_norm(d$chr); nch <- sum(c0 != d$chr, na.rm = TRUE); if (nch) ch <- c(ch, sprintf("Chromosome names written in one form (1 to 22, X, Y, MT) in %s rows (for example chr1 or 23).", fm(nch))) }
  # P values: -log10(P) is converted; the log P is kept so that very small P values (which round to 0) still give a standard error
  lp <- NULL
  if ("p" %in% names(d)) {
    if (identical(facts$p_scale, "mlog10")) {
      neg <- sum(d$p < 0, na.rm = TRUE); if (neg) err <- c(err, sprintf("%s values in the -log10(P) column are negative, so it is not -log10(P). Choose \"P value\" for that column instead.", fm(neg)))
      lp <- -d$p * log(10); d$p <- 10^(-d$p); ch <- c(ch, "P values converted from -log10(P).")
    } else {
      bp <- sum(d$p < 0 | d$p > 1, na.rm = TRUE)
      if (bp) { if (mean(d$p > 1, na.rm = TRUE) > 0.01) err <- c(err, sprintf("%s P values are above 1. If this column holds -log10(P), choose \"-log10(P)\" for the P value column.", fm(bp))) else warn <- c(warn, sprintf("%s P values are outside 0 to 1.", fm(bp))) }
      lp <- log(d$p)
    }
  }
  if (identical(facts$scale, "or") && "beta" %in% names(d)) {
    if (any(d$beta <= 0, na.rm = TRUE)) err <- c(err, "Odds ratios must be positive, but zero or negative values were found. If the column holds log odds ratios, choose \"Beta or log odds\" instead.")
    else { d$beta <- log(d$beta); ch <- c(ch, "Odds ratios converted to log odds ratios.") }
  }
  # GRCh38 to GRCh37 (MR Studio works on GRCh37)
  if (identical(facts$build, "38") && all(c("chr", "pos") %in% names(d))) {
    m <- liftover_38_to_37(d$chr, d$pos); why <- table(m$why[is.na(m$pos)])
    rev <- !is.na(m$strand) & m$strand == "-"
    d$chr <- m$chr; d$pos <- m$pos
    if (any(rev)) { d$ea[rev] <- .complement(d$ea[rev]); d$oa[rev] <- .complement(d$oa[rev]) }
    ch <- c(ch, sprintf("Positions converted from GRCh38 to GRCh37 for %s (Ensembl chain file)%s.", vv(sum(!is.na(m$pos))), if (any(rev)) sprintf("; %s on reverse-strand blocks got complemented alleles", fm(sum(rev))) else ""))
    if (length(why)) warn <- c(warn, sprintf("%s could not be converted to GRCh37 and will be dropped (%s).", vv(sum(why)), paste(sprintf("%s: %s", names(why), fm(as.integer(why))), collapse = "; ")))
  }
  d$ea <- toupper(trimws(d$ea)); d$oa <- toupper(trimws(d$oa))
  if (identical(facts$ea_ok, "swap")) { tmp <- d$ea; d$ea <- d$oa; d$oa <- tmp }
  # Z scores to effect and SE: b = z / sqrt(2 f (1 - f) (n + z^2)), se = 1 / sqrt(2 f (1 - f) (n + z^2)) (Zhu et al. 2016), per SD of the trait.
  # 2 f (1 - f) is the same for either allele, so a minor allele frequency serves as well.
  if (identical(facts$scale, "z") && "z" %in% names(d)) {
    nv <- if ("n" %in% names(d)) d$n else rep(suppressWarnings(as.numeric(facts$n)), nrow(d))
    f <- if ("eaf" %in% names(d) && !identical(facts$frq, "none")) d$eaf else rep(NA_real_, nrow(d))
    if (all(is.na(f)) && !final) ch <- c(ch, "Effect sizes and standard errors will be derived from the Z scores once the frequencies are taken from the 1000 Genomes panel (when saving).")
    else {
      den <- sqrt(2 * f * (1 - f) * (nv + d$z^2)); okz <- is.finite(den) & den > 0
      d$beta <- ifelse(okz, d$z / den, NA_real_); d$se <- ifelse(okz, 1 / den, NA_real_)
      ch <- c(ch, sprintf("Effect sizes and standard errors derived from the Z scores, the sample size and the allele frequency for %s (Zhu et al. 2016; per standard deviation of the trait).", vv(sum(okz))))
      if (sum(!okz)) warn <- c(warn, sprintf("%s without a usable Z score, sample size or frequency will be dropped.", vv(sum(!okz))))
    }
  }
  # SE from the effect and its P value: se = |b| / z, with z the normal quantile of P / 2 (Altman and Bland 2011)
  if (!("se" %in% names(d)) && "beta" %in% names(d) && !is.null(lp) && !identical(facts$scale, "z")) {
    zz <- stats::qnorm(lp - log(2), lower.tail = FALSE, log.p = TRUE); se <- abs(d$beta) / zz
    oks <- is.finite(se) & se > 0; d$se <- ifelse(oks, se, NA_real_)
    ch <- c(ch, sprintf("Standard errors derived from the effect and its P value for %s (Altman and Bland 2011). The P values must be two-sided and from the same test as the effect; rounded P values give approximate standard errors.", vv(sum(oks))))
    if (sum(!oks)) warn <- c(warn, sprintf("%s with an effect of 0, a P value of 1 or no P value will be dropped (no standard error can be derived).", vv(sum(!oks))))
  }
  list(d = d, changes = ch, warnings = warn, errors = err)
}

#' Read the mapped columns completely and report what is wrong. errors block the import; warnings are passed on to the user.
import_validate <- function(path, sep, columns, mapping, facts, max_rows = 5e6) {
  err <- character(); warn <- character(); mp <- mapping[nzchar(unlist(mapping))]
  err <- c(err, import_missing(mp, facts))
  if (!all(unlist(mp) %in% columns)) err <- c(err, "A chosen column is not in the file.")
  if (anyDuplicated(unlist(mp))) err <- c(err, "The same column is used for two different roles.")
  if (length(err)) return(list(ok = FALSE, errors = err, warnings = warn, changes = list(), stats = list()))
  d <- .import_read(path, sep, columns, mp, max_rows); nr <- nrow(d)
  num <- function(x) suppressWarnings(as.numeric(x))
  tr <- .import_transform(d, facts); d <- tr$d; err <- c(err, tr$errors); warn <- c(warn, tr$warnings)
  if (nr >= max_rows) warn <- c(warn, sprintf("Only the first %s rows were read.", format(max_rows, big.mark = ",")))
  if ("se" %in% names(d)) { bad_se <- sum(d$se <= 0 | is.na(d$se)); if (bad_se) warn <- c(warn, sprintf("%s variants have a missing or non-positive standard error and will be dropped.", format(bad_se, big.mark = ","))) }
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
  if (identical(facts$independent, "clump") && !identical(facts$role, "outcome")) {
    if (!("snp" %in% names(d)) || !any(grepl("^rs[0-9]+$", d$snp))) err <- c(err, "Clumping on this computer needs rsIDs (for example rs12345), and this file has none in the chosen column.")
    r2 <- suppressWarnings(as.numeric(facts$clump_r2 %||% 0.001)); kb <- suppressWarnings(as.numeric(facts$clump_kb %||% 10000))
    if (length(r2) != 1L || !is.finite(r2) || r2 <= 0 || r2 >= 1) err <- c(err, "The LD threshold r2 must be between 0 and 1.")
    if (length(kb) != 1L || !is.finite(kb) || kb < 1 || kb > 100000) err <- c(err, "The clumping window must be between 1 and 100,000 kb.")
    if (!(toupper(facts$clump_pop %||% facts$ancestry %||% "") %in% c("EUR", "EAS", "AFR", "AMR", "SAS"))) err <- c(err, "Choose the 1000 Genomes ancestry to clump with (European, East Asian, African, Admixed American or South Asian).")
  }
  bse <- if (all(c("beta", "se") %in% names(d))) !is.na(d$beta) & !is.na(d$se) & d$se > 0 else rep(TRUE, nr)   # Z scores waiting for panel frequencies
  keep <- ok_al & bse & !duplicated(key) & (if ("pos" %in% names(d) && identical(facts$build, "38")) !is.na(d$pos) else TRUE)
  list(ok = !length(err), errors = err, warnings = warn, changes = as.list(tr$changes), stats = list(rows = nr, kept = sum(keep), palindromic = pal, snps = if ("snp" %in% names(d)) sum(grepl("^rs[0-9]+$", d$snp)) else 0L, has_eaf = "eaf" %in% names(d) && identical(facts$frq, "eaf")))
}

#' Where imported datasets, with their recipe, live.
library_dir <- function(root = studies_root()) file.path(root, "Data library")

# ---- Data library tables (Round 35) ---------------------------------------------------------------------------------
# Every table the Data library keeps for a dataset (the standard copy "standardised", and for OpenGWAS downloads
# "candidates", "instrument_selection" and "instrument_selection_ties") is one Parquet file (zstd), written with nanoparquet.
# Parquet keeps every number exactly and stays small and fast for full genome-wide files. Tables saved before Round 35 as
# CSV are converted the first time they are read: the CSV is removed only after the Parquet copy reads back identical.
# Files the analyst opens (study folders, analysis tables, exports) stay CSV.
LIB_CHAR_COLS <- c("snp", "chr", "ea", "oa", "variant_id", "status", "index_variant")
lib_table_path <- function(dir, name) file.path(dir, paste0(name, ".parquet"))
lib_table_exists <- function(dir, name) file.exists(lib_table_path(dir, name)) || file.exists(file.path(dir, paste0(name, ".csv")))
lib_table_write <- function(d, dir, name) {
  d <- as.data.frame(d, stringsAsFactors = FALSE); rownames(d) <- NULL
  tmp <- tempfile(paste0(".", name, "-"), dir, fileext = ".parquet")
  nanoparquet::write_parquet(d, tmp, compression = "zstd")
  out <- lib_table_path(dir, name)
  if (!file.rename(tmp, out)) { unlink(tmp); stop("Could not save ", basename(out), " in ", dir, call. = FALSE) }
  invisible(out)
}
# Character columns read back as "" where empty, as the CSV copies were read before Round 35.
.lib_norm <- function(d) { d <- as.data.frame(d, stringsAsFactors = FALSE); rownames(d) <- NULL
  for (n in names(d)) if (is.character(d[[n]])) d[[n]][is.na(d[[n]])] <- ""
  d }
lib_table_read <- function(dir, name) {
  pq <- lib_table_path(dir, name)
  if (file.exists(pq)) return(.lib_norm(nanoparquet::read_parquet(pq)))
  cv <- file.path(dir, paste0(name, ".csv"))
  if (!file.exists(cv)) stop("The table \"", name, "\" of \"", basename(dir), "\" is not in the Data library any more.", call. = FALSE)
  hd <- names(utils::read.csv(cv, nrows = 1, check.names = FALSE)); ch <- intersect(LIB_CHAR_COLS, hd)
  d <- .lib_norm(utils::read.csv(cv, stringsAsFactors = FALSE, check.names = FALSE, colClasses = stats::setNames(rep("character", length(ch)), ch)))
  ok <- tryCatch({ lib_table_write(d, dir, name); isTRUE(all.equal(d, .lib_norm(nanoparquet::read_parquet(pq)), check.attributes = FALSE)) }, error = function(e) FALSE)
  if (ok) unlink(cv) else unlink(pq)   # a copy that does not read back identical is dropped and the CSV kept
  d
}

#' Write the standard copy and the recipe (file name, size, sha256, mapping, facts, every change made: FR-75). The original file is never
#' changed. path is the file read (a cleaned copy for Excel and similar files); orig_path is the user's file, whose checksum is recorded.
import_save <- function(path, orig_name, sep, columns, mapping, facts, root = studies_root(), orig_path = path, prep_notes = character()) {
  v <- import_validate(path, sep, columns, mapping, facts); if (!v$ok) stop(paste(v$errors, collapse = " "))
  mp <- mapping[nzchar(unlist(mapping))]
  nm <- clean_study_name(facts$name); if (!nzchar(nm)) stop("Give the trait a name.")
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) stop("A dataset named \"", nm, "\" already exists in the Data library. Choose another name.")
  p <- study_path_problem(library_dir(root)); if (!is.null(p)) stop(p)
  d <- .import_read(path, sep, columns, mp)
  # Round 77: frequencies from the 1000 Genomes panel, matched by rsID and allele; palindromic variants only on the forward strand
  pan <- character()
  if (identical(facts$eaf_fill, "panel")) {
    dir.create(dest, recursive = TRUE); pop <- toupper(facts$ancestry)
    fq <- ld_panel_freq_cached(dest, d$snp, pop); i <- match(d$snp, fq$SNP); ea <- toupper(trimws(d$ea)); oa <- toupper(trimws(d$oa))
    if (identical(facts$ea_ok, "swap")) { tmp <- ea; ea <- oa; oa <- tmp }
    f <- ifelse(!is.na(i) & ea == fq$A1[i] & oa == fq$A2[i], fq$MAF[i], ifelse(!is.na(i) & ea == fq$A2[i] & oa == fq$A1[i], 1 - fq$MAF[i], NA_real_))
    pal <- (ea == "A" & oa == "T") | (ea == "T" & oa == "A") | (ea == "C" & oa == "G") | (ea == "G" & oa == "C")
    if (!identical(facts$strand, "forward")) f[pal] <- NA_real_
    d$eaf <- as.character(f); mp$eaf <- "eaf"; facts$frq <- "eaf"
    pan <- sprintf("Effect allele frequencies taken from the 1000 Genomes %s panel for %s of %s variants (matched by rsID and alleles%s).", pop, format(sum(!is.na(f)), big.mark = ","), format(nrow(d), big.mark = ","), if (identical(facts$strand, "forward")) "" else "; palindromic variants left without a frequency because the strand is not known")
  }
  tr <- .import_transform(d, facts, final = TRUE); d <- tr$d
  key <- if ("snp" %in% names(d)) d$snp else paste(d$chr, d$pos, d$ea, d$oa)
  keep <- grepl("^[ACGT]$", d$ea) & grepl("^[ACGT]$", d$oa) & !is.na(d$beta) & !is.na(d$se) & d$se > 0 & !duplicated(key) & (if (identical(facts$build, "38")) !is.na(d$pos) else TRUE)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) { unlink(dest, recursive = TRUE); stop("No variant is left after the checks, so nothing was saved.") }
  if (!identical(facts$frq, "eaf")) d$eaf <- NULL
  d$z <- NULL; d$p_mlog10 <- NULL
  if (!("n" %in% names(d)) && is.finite(suppressWarnings(as.numeric(facts$n)))) d$n <- as.numeric(facts$n)
  nc <- suppressWarnings(as.numeric(if (is.null(facts$ncase)) NA else facts$ncase)); nn <- suppressWarnings(as.numeric(if (is.null(facts$ncontrol)) NA else facts$ncontrol))
  if (!("n" %in% names(d)) && identical(facts$kind, "binary") && length(nc) == 1L && length(nn) == 1L && is.finite(nc) && is.finite(nn)) d$n <- nc + nn   # total N of a case-control study
  if (identical(facts$scale, "z")) facts$scale_saved <- "beta"   # per SD of the trait
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  lib_table_write(d, dest, "standardised")
  rec <- list(trait = nm, role = facts$role, facts = facts, mapping = mp, original = list(file = orig_name, bytes = file.size(orig_path), sha256 = sha256_file(orig_path)),
              rows_in_file = v$stats$rows, rows_kept = nrow(d), warnings = v$warnings, changes = as.list(c(prep_notes, pan, tr$changes)), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              build_saved = if ((facts$build %||% "unknown") %in% c("37", "38")) "GRCh37" else "as in the file (build not known)", liftover = if (identical(facts$build, "38")) LIFTOVER_CHAIN[c("file", "sha256", "source")] else NULL,
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
  rf <- file.path(dir, "recipe.json")
  if (!file.exists(rf) || !lib_table_exists(dir, "standardised")) stop("The dataset \"", trait, "\" is not in the Data library any more.")
  rec <- jsonlite::fromJSON(rf, simplifyVector = FALSE)
  list(recipe = rec, data = lib_table_read(dir, "standardised"), dir = dir)
}
