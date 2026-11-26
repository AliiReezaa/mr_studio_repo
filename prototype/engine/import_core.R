# Data import (FR-71 to FR-78): read a GWAS summary file, suggest the column roles, check it, save a standard copy, remember the layout.
# Plain base R on purpose, so the bundled runtime needs nothing extra. Individual-level data is refused (S-06).

IMPORT_ROLES <- c(snp = "Variant ID (rsID)", chr = "Chromosome", pos = "Position", ea = "Effect allele", oa = "Other allele", eaf = "Effect allele frequency",
                  beta = "Effect size", se = "Standard error", p = "P value", n = "Sample size", z = "Z score",
                  ci_lower = "Lower confidence limit", ci_upper = "Upper confidence limit", test = "Test (PLINK 2, REGENIE)", oa_alt = "Second other-allele column (PLINK 2 REF and ALT)")
IMPORT_REQUIRED <- c("snp", "ea", "oa", "beta", "se")   # Round 77: the effect or the SE may be derived (see import_missing)

.norm_col <- function(x) gsub("[^a-z0-9]", "", tolower(x))

IMPORT_SYNONYMS <- list(
  snp = c("snp", "rsid", "rsids", "rs", "markername", "marker", "variantid", "variant", "snpid", "id", "rsnumber", "rsnum", "markerid"),
  chr = c("chr", "chrom", "chromosome", "chrcode", "chromosomenumber", "hgchr"),
  pos = c("pos", "position", "bp", "basepair", "basepairposition", "bppos", "bpposition", "genpos", "chrpos", "hgpos", "basepairlocation", "bplocation"),
  ea = c("ea", "effectallele", "a1", "allele1", "alt", "testedallele", "eallele", "effallele", "codedallele", "tested"),
  oa = c("oa", "otherallele", "a2", "allele2", "ref", "nea", "noneffectallele", "reference", "refallele", "nonefallele", "othalle", "allele0", "a0", "nonea"),
  eaf = c("eaf", "effectallelefreq", "effectallelefrequency", "freq", "frq", "af", "a1freq", "freq1", "afreq", "af1", "effectaf", "a1frq", "afalt", "altfreq", "afallele2", "a2freq", "freqtestedallele", "testedallelefreq", "frqa1", "freqa1", "eafa1",
          "maf", "minorallelefreq", "minorallelefrequency"),   # Round 77: a MAF heading is matched here but marked as a minor allele frequency (import_hints)
  beta = c("beta", "b", "effect", "effectsize", "es", "logor", "logodds", "estimate", "betahat", "or", "oddsratio"),
  se = c("se", "stderr", "stderror", "standarderror", "sebeta", "sebetahat", "selogor", "logorse", "orse", "selnor", "selogodds"),
  p = c("p", "pval", "pvalue", "pvalu", "pvalues", "pvalueassoc", "plrt", "pwald", "pbolt", "pboltlmminf", "pboltlmm",
        "negativelog10p", "neglog10p", "mlog10p", "minuslog10p", "log10p", "neglogp", "mlogp", "logp", "lp", "neglog10pval", "mlog10pval", "log10pvalue", "neglog10pvalue"),   # Round 77: -log10(P) headings are converted (import_hints)
  n = c("n", "samplesize", "neff", "totaln", "nsamples", "ntotal", "nsample", "weight", "nobs", "ncompletesamples", "obsct"),
  z = c("z", "zscore", "zstat", "zstatistic", "zvalue", "zval"),
  # Round 79: a confidence interval instead of an SE (papers, GWAS Catalog), and the PLINK 2 / REGENIE test column
  ci_lower = c("cilower", "lci", "lowerci", "lower95", "l95", "lower95ci", "ci95lower", "cilow", "lowerlimit", "orlower95", "orlci", "lowerbound", "orl95", "betalower", "betal95", "lcl"),
  ci_upper = c("ciupper", "uci", "upperci", "upper95", "u95", "upper95ci", "ci95upper", "cihigh", "upperlimit", "orupper95", "oruci", "upperbound", "oru95", "betaupper", "betau95", "ucl"),
  test = c("test"))

#' Layers 1 and 2 of FR-78: normalised synonyms (high) and fuzzy matches (medium, with a score). Layer 3 looks at the values (low).
import_suggest <- function(cols, rows) {
  nc <- sub("(hg19|hg38|grch37|grch38|b37|b38|build37|build38)$", "", .norm_col(cols)); out <- list(); used <- character()   # Round 77: POS_hg38 is a position
  nc <- ifelse(grepl("^hm.", nc), sub("^hm", "", nc), nc)   # Round 79: GWAS Catalog harmonised columns (hm_beta, hm_rsid, ...)
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
                   eaf = all(!is.na(num)) && all(num > 0 & num < 1) && grepl("fr|af", .norm_col(col)), FALSE)   # Round 79: an interval column is not a frequency
      if (ok) { pick(role, col, "low", "Guessed from the values, please confirm"); free <- setdiff(free, col); break }
    }
  }
  # Round 79: a variant column without rsIDs (chr:pos:ref:alt, chr_pos_a1_a2) is not an rsID column when chromosome and position are there
  if (!is.null(out$snp) && !is.null(out$chr) && !is.null(out$pos)) { v <- vals(out$snp$column); if (length(v) && !any(grepl("^rs[0-9]", v))) out$snp <- NULL }
  # Round 79: SAIGE reports the effect of Allele2 (AF_Allele2), not Allele1
  hn <- function(r) if (is.null(out[[r]])) "" else .norm_col(out[[r]]$column)
  if (grepl("allele2|a2freq", hn("eaf")) && hn("ea") %in% c("allele1", "a1") && hn("oa") %in% c("allele2", "a2")) {
    e <- out$ea; out$ea <- out$oa; out$oa <- e
    out$ea$why <- "The frequency column is for Allele2, so Allele2 is the effect allele (as in SAIGE output)"; out$oa$why <- "Allele1 is the reference allele (as in SAIGE output)"
  }
  # Round 79: PLINK 2 --glm: A1 is the effect allele and is either REF or ALT; the other allele is the one of REF and ALT that is not A1
  if (hn("ea") == "a1" && all(c("ref", "alt") %in% nc)) {
    out$oa <- list(column = cols[match("ref", nc)], confidence = "high", why = "PLINK 2: the other allele is REF, or ALT where A1 is REF")
    out$oa_alt <- list(column = cols[match("alt", nc)], confidence = "high", why = "PLINK 2: used where A1 is the REF allele")
  }
  out
}

IMPORT_MAF_HEADS <- c("maf", "minorallelefreq", "minorallelefrequency")
IMPORT_MLOG10_HEADS <- c("negativelog10p", "neglog10p", "mlog10p", "minuslog10p", "log10p", "neglogp", "mlogp", "logp", "lp", "neglog10pval", "mlog10pval", "log10pvalue", "neglog10pvalue")

#' Round 77: what the headings and first rows say about the answers on the facts step (each one is shown and can be changed).
import_hints <- function(cols, rows, sug, build_check = TRUE) {
  hc <- function(r) if (is.null(sug[[r]])) "" else { x <- .norm_col(sug[[r]]$column); if (grepl("^hm.", x)) sub("^hm", "", x) else x }
  vals <- function(col) suppressWarnings(as.numeric(vapply(rows, function(r) as.character(r[[col]] %||% NA), "")))
  out <- list()
  if (!is.null(sug$eaf)) out$frq <- if (hc("eaf") %in% IMPORT_MAF_HEADS) "maf" else "eaf"
  if (!is.null(sug$p)) { v <- vals(sug$p$column); out$p_scale <- if (hc("p") %in% IMPORT_MLOG10_HEADS || any(v > 1, na.rm = TRUE)) "mlog10" else "p" }
  if (!is.null(sug$beta) && hc("beta") %in% c("or", "oddsratio")) out$scale <- "or"
  if (!is.null(sug$beta) && is.null(out$scale)) { v <- vals(sug$beta$column); v <- v[is.finite(v)]; if (length(v) >= 3 && all(v > 0) && stats::median(v) > 0.7 && stats::median(v) < 1.4 && grepl("or|odds|hr|hazard", hc("beta"))) out$scale <- "or" }
  if (is.null(sug$beta) && !is.null(sug$z)) out$scale <- "z"
  if (identical(out$scale, "z") && is.null(sug$eaf) && !is.null(sug$snp) && !is.na(tryCatch(.panel_pop_any(), error = function(e) NA))) out$eaf_fill <- "panel"   # Round 79: Z scores need a frequency
  all_h <- tolower(paste(cols, collapse = " "))
  if (grepl("hg38|grch38|b38|build38", all_h)) out$build <- "38" else if (grepl("hg19|grch37|b37|build37", all_h)) out$build <- "37"
  if (build_check) { bg <- tryCatch(import_build_guess(rows, sug), error = function(e) NULL); if (!is.null(bg)) { out$build <- bg$build; out$build_why <- bg$why } }   # Round 79
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
  # Round 79: a chr:pos:ref:alt variant column with no position columns, and combined columns of a published table (a sample of rows)
  spl <- function(l) gsub('^"|"$', "", if (sep == " ") strsplit(trimws(l), "[ ]+")[[1]] else strsplit(l, sep, fixed = TRUE)[[1]])
  hdr <- spl(sub("^#", "", top[skip + 1L])); smp <- lapply(utils::head(rest[nzchar(trimws(rest))], 40), spl)
  sdf <- as.data.frame(lapply(seq_along(hdr), function(j) vapply(smp, function(r) if (j <= length(r)) r[j] else "", "")), stringsAsFactors = FALSE, col.names = make.unique(hdr))
  vj <- if (nrow(sdf)) .split_varid(hdr, sdf) else NULL
  comb <- nrow(sdf) > 0 && !identical(.split_combined(hdr, sdf)$head, hdr)
  if (comb && file.size(path) > 50e6) { comb <- FALSE; notes <- c(notes, "Some columns seem to hold two values (for example an estimate with its interval). The file is too large to split them here; split them first or choose the columns by hand.") }
  if (!bom && !skip && !hash_head && !dec_comma && is.null(vj) && !comb) return(list(path = path, notes = notes))
  if (bom) notes <- c(notes, "A byte-order mark at the start of the file was removed.")
  if (skip) notes <- c(notes, sprintf("%d comment or empty line%s before the heading %s skipped.", skip, if (skip > 1) "s" else "", if (skip > 1) "were" else "was"))
  if (hash_head) notes <- c(notes, "The heading line starts with #; the # was removed.")
  if (dec_comma) notes <- c(notes, "Numbers use a decimal comma (for example 0,05); they were read as decimal points.")
  out <- tempfile("import-", fileext = ".tsv"); oc <- file(out, "w")
  close(con); con <- .open_any(path); first <- TRUE; n_skip <- skip; n_head_done <- FALSE
  repeat {
    L <- readLines(con, n = 200000L, warn = FALSE); if (!length(L)) break
    if (n_skip) { k <- min(n_skip, length(L)); L <- L[-seq_len(k)]; n_skip <- n_skip - k; if (!length(L)) next }
    if (first) { L[1] <- sub("^\ufeff", "", sub("^\xef\xbb\xbf", "", L[1], useBytes = TRUE)); if (hash_head) L[1] <- sub("^#", "", L[1]); first <- FALSE }
    if (sep != "\t") L <- if (sep == " ") gsub("[ ]+", "\t", trimws(L)) else gsub(sep, "\t", L, fixed = TRUE)
    if (dec_comma) L <- gsub("(^|\t) *(-?[0-9]+),([0-9]+([eE][-+]?[0-9]+)?) *(?=\t|$)", "\\1\\2.\\3", L, perl = TRUE)
    if (!is.null(vj)) {   # Round 79: chromosome, position and alleles appended from the variant column
      hd1 <- !n_head_done; v <- vapply(strsplit(L, "\t", fixed = TRUE), function(r) if (length(r) >= vj) gsub('^"|"$', "", r[vj]) else "", "")
      vc <- .varid_cols(v); L <- paste(L, vc$CHR, vc$POS, vc$REF, vc$ALT, sep = "\t"); if (hd1) L[1] <- paste(sub("(\t[^\t]*){4}$", "", L[1]), "CHR", "POS", "REF", "ALT", sep = "\t")
    }
    n_head_done <- TRUE
    writeLines(L, oc)
  }
  close(oc); on.exit(close(con)); oc <- NULL
  if (!is.null(vj)) notes <- c(notes, sprintf("Chromosome, position and alleles were read from the column \"%s\" (chr:pos:ref:alt); the effect is taken to be for ALT, the last allele, as in Neale lab and most pipelines.", hdr[vj]))
  if (comb) {   # a published table: read whole, split the combined columns, write again
    d <- utils::read.table(out, header = TRUE, sep = "\t", quote = "\"", colClasses = "character", check.names = FALSE, comment.char = "", fill = TRUE, na.strings = character())
    sc <- .split_combined(names(d), d); notes <- c(notes, sc$notes); names(sc$body) <- make.unique(sc$head)
    .write_tsv_bytes(rbind(names(sc$body), as.matrix(sc$body)), out)
  }
  list(path = out, notes = notes)
}

# Round 79: a table written byte for byte, so typeset or accented characters are not turned into <U+00D7> text in a non-UTF-8 locale
.write_tsv_bytes <- function(m, path) {
  m <- as.matrix(m); m[is.na(m)] <- ""; m[] <- .ascii_typeset(m)
  con <- file(path, "wb"); on.exit(close(con)); writeLines(apply(m, 1, paste, collapse = "\t"), con, useBytes = TRUE)
}
# Round 79: numbers as papers and spreadsheets print them: "1.2 × 10−8", "1.2x10^-8", "1.2E-08", a Unicode minus sign, "< 1e-300"
# Typeset characters become ASCII first, byte by byte, so the patterns below are plain ASCII and work in any locale (Windows included)
.ascii_typeset <- function(x) {
  x <- as.character(x)   # no enc2utf8: in a non-UTF-8 locale it turns these characters into text such as <c3><97>
  # built with intToUtf8: a "\\u2212" written in the code becomes the text <U+2212> when R loads it in a non-UTF-8 locale
  rp <- list(c(0x2212, "-"), c(0x2013, " to "), c(0x2014, " to "), c(0x00d7, "x"), c(0x00b7, "x"), c(0x2264, "<"), c(0x2265, ">"), c(0x00a0, " "))
  for (k in rp) x <- gsub(rawToChar(charToRaw(intToUtf8(as.integer(k[1])))), k[2], x, fixed = TRUE, useBytes = TRUE)   # compared as bytes
  x
}
.num_parse <- function(x) {
  x <- trimws(.ascii_typeset(x)); x <- sub("^[<>]\\s*", "", x)
  x <- sub("^([-+]?[0-9]*\\.?[0-9]+)\\s*[x*]\\s*10\\s*\\^?\\s*([-+]?[0-9]+)$", "\\1e\\2", x, perl = TRUE)
  x <- sub("^10\\s*\\^\\s*([-+]?[0-9]+)$", "1e\\1", x)
  suppressWarnings(as.numeric(x))
}

# Round 79: combined columns of published tables, split into one value per column:
#   "T/C" (effect/other allele), "0.15 (0.01)" (estimate and SE), "1.12 (1.05–1.19)" (estimate and 95% CI), "1:109814880" (chromosome and position)
.NUM <- "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
.split_combined <- function(head, body) {
  notes <- character(); out_h <- character(); out <- list()
  share <- function(v, re) { v <- v[nzchar(v)]; length(v) >= 1 && mean(grepl(re, v, perl = TRUE)) >= 0.8 }
  for (j in seq_along(head)) {
    h <- head[j]; v <- trimws(.ascii_typeset(body[[j]])); hl <- tolower(.ascii_typeset(h))
    re_al <- "^[ACGTacgt]+\\s*/\\s*[ACGTacgt]+$"
    re_se <- sprintf("^(%s)\\s*[\\(\\[]\\s*(%s)\\s*[\\)\\]]$", .NUM, .NUM)
    re_ci <- sprintf("^(%s)\\s*[\\(\\[]\\s*(%s)\\s*(?:to|,|;|\\s-\\s|-)\\s*(%s)\\s*[\\)\\]]$", .NUM, .NUM, .NUM)
    re_cp <- "^(chr)?([0-9]{1,2}|X|Y|MT):[0-9,]+$"
    if (share(v, re_al)) {
      second <- grepl("^(other|non|ref|a2|allele 2)", hl); a <- toupper(sub("\\s*/.*$", "", v)); b <- toupper(sub("^.*/\\s*", "", v))
      out_h <- c(out_h, "Effect allele", "Other allele"); out <- c(out, list(if (second) b else a, if (second) a else b))
      notes <- c(notes, sprintf("The column \"%s\" holds two alleles (for example %s); it was split into the effect allele (%s) and the other allele.", h, v[nzchar(v)][1], if (second) "the second" else "the first"))
    } else if (share(v, re_ci)) {
      e <- if (grepl("\\bor\\b|odds", hl)) "OR" else if (grepl("\\bhr\\b|hazard", hl)) "HR" else "Beta"
      out_h <- c(out_h, e, paste(e, "CI lower"), paste(e, "CI upper")); out <- c(out, list(sub(re_ci, "\\1", v, perl = TRUE), sub(re_ci, "\\2", v, perl = TRUE), sub(re_ci, "\\3", v, perl = TRUE)))
      notes <- c(notes, sprintf("The column \"%s\" holds an estimate with its confidence interval (for example %s); it was split into the estimate and the two limits, from which the standard error is derived.", h, v[nzchar(v)][1]))
    } else if (share(v, re_se) && grepl("se|standard error|s\\.e\\.", hl)) {
      e <- if (grepl("\\bor\\b|odds", hl)) "OR" else "Beta"
      out_h <- c(out_h, e, "SE"); out <- c(out, list(sub(re_se, "\\1", v, perl = TRUE), sub(re_se, "\\2", v, perl = TRUE)))
      notes <- c(notes, sprintf("The column \"%s\" holds an estimate with its standard error (for example %s); it was split into the two.", h, v[nzchar(v)][1]))
    } else if (share(v, re_cp)) {
      bs <- regmatches(hl, regexpr("hg19|hg38|grch37|grch38|b37|b38", hl)); sfx <- if (length(bs)) paste0(" (", bs, ")") else ""
      out_h <- c(out_h, "Chromosome", paste0("Position", sfx)); out <- c(out, list(sub("^(chr)?([^:]+):.*$", "\\2", v), gsub(",", "", sub("^[^:]+:", "", v))))
      notes <- c(notes, sprintf("The column \"%s\" holds chromosome and position together (for example %s); it was split into the two.", h, v[nzchar(v)][1]))
    } else { out_h <- c(out_h, h); out <- c(out, list(body[[j]])) }
  }
  list(head = out_h, body = as.data.frame(stats::setNames(out, make.unique(out_h)), stringsAsFactors = FALSE, check.names = FALSE), notes = notes)
}

# Round 79: a variant column of the form chr:pos:ref:alt (Neale lab UK Biobank, many pipelines) with no separate chromosome, position or allele
# columns: these are read from it. The effect is for ALT (the last allele), as in those sources; the allele question on the facts step checks it.
.RE_VARID <- "^(chr)?([0-9]{1,2}|X|Y|MT)[:_]([0-9]+)[:_]([ACGT]+)[:_]([ACGT]+)$"
.split_varid <- function(head, body) {
  nh <- sub("^hm(.)", "\\1", .norm_col(head)); if (any(nh %in% c(IMPORT_SYNONYMS$chr, IMPORT_SYNONYMS$pos))) return(NULL)   # positions already there
  for (j in seq_along(head)) { v <- trimws(body[[j]]); v <- v[nzchar(v)]
    if (length(v) && mean(grepl(.RE_VARID, v, perl = TRUE)) >= 0.8 && !any(grepl("^rs[0-9]", v))) return(j) }
  NULL
}
.varid_cols <- function(v) { v <- trimws(v); ok <- grepl(.RE_VARID, v, perl = TRUE)
  f <- function(k) ifelse(ok, sub(.RE_VARID, paste0("\\", k), v, perl = TRUE), "")
  data.frame(CHR = f(2), POS = f(3), REF = f(4), ALT = f(5), stringsAsFactors = FALSE) }

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
  sc <- .split_combined(head, body); head <- sc$head; body <- sc$body; notes <- c(notes, sc$notes)   # Round 79
  vj <- .split_varid(head, body); if (!is.null(vj)) { body <- cbind(body, .varid_cols(body[[vj]])); head <- c(head, "CHR", "POS", "REF", "ALT")
    notes <- c(notes, sprintf("Chromosome, position and alleles were read from the column \"%s\" (chr:pos:ref:alt); the effect is taken to be for ALT, the last allele.", head[vj])) }
  out <- tempfile("import-", fileext = ".tsv")
  .write_tsv_bytes(rbind(head, as.matrix(body)), out)
  list(path = out, notes = notes, sheets = as.list(sheets), sheet = sheet)
}

.detect_sep <- function(line) { cands <- c("\t", ",", ";", " "); cnt <- vapply(cands, function(s) lengths(regmatches(line, gregexpr(s, line, fixed = TRUE))), 0L); cands[which.max(cnt)] }
.read_tbl <- function(path, ...) { con <- .open_any(path); on.exit(close(con)); utils::read.table(con, ...) }
.open_any <- function(path) if (grepl("\\.b?gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")   # Round 79: .bgz (block gzip) reads as gzip

#' The start of a file: separator, headings and the first rows. Refuses formats and layouts that are not supported.
import_read_head <- function(path, name = basename(path), nrows = 8L, sheet = NULL) {
  if (grepl("\\.(zip|vcf)(\\.b?gz)?$", name, ignore.case = TRUE)) stop("This file type is not supported here. Supported: text tables (TSV, CSV, space or semicolon separated, optionally compressed as .gz or .bgz) and Excel workbooks (.xlsx, .xls); GWAS-VCF files have their own import.")
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

# Round 79: the 1000 Genomes panel (GRCh37) as a reference for rsIDs and positions, read by PLINK (the variant list is about 250 MB per
# ancestry, so it is never loaded into R). Any installed ancestry serves: positions and rsIDs are the same in all of them.
.panel_pop_any <- function(pref = "EUR") { for (p in unique(c(toupper(pref), "EUR", "EAS", "AFR", "AMR", "SAS"))) if (.panel_ready(p)) return(p); NA_character_ }
.panel_bim <- function(snps = NULL, ranges = NULL, pop = .panel_pop_any()) {
  if (is.na(pop)) return(NULL)
  tmp <- tempfile("bim-"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  ex <- if (!is.null(snps)) { writeLines(unique(snps[!is.na(snps)]), file.path(tmp, "x.txt")); c("--extract", shQuote(file.path(tmp, "x.txt"))) }
        else { utils::write.table(ranges, file.path(tmp, "r.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); c("--extract", "range", shQuote(file.path(tmp, "r.txt"))) }
  suppressWarnings(system2(plink_bin(), c("--bfile", shQuote(file.path(ld_panel_dir(), pop)), ex, "--make-just-bim", "--threads", "1", "--out", shQuote(file.path(tmp, "o"))), stdout = FALSE, stderr = FALSE))
  f <- file.path(tmp, "o.bim"); if (!file.exists(f)) return(data.frame(chr = character(), snp = character(), pos = numeric(), a1 = character(), a2 = character()))
  b <- utils::read.table(f, colClasses = c("character", "character", "NULL", "numeric", "character", "character"), col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
  b$chr <- .chr_norm(b$chr); b
}
#' Round 79: the genome build of a file with rsIDs and positions: where its first rsIDs sit in the panel (GRCh37), directly or after
#' converting the file's positions from GRCh38. NULL when it cannot be told (no panel, no rsIDs or positions, too few matches).
import_build_guess <- function(rows, sug) {
  if (is.null(sug$snp) || is.null(sug$chr) || is.null(sug$pos)) return(NULL)
  g <- function(col) vapply(rows, function(r) as.character(r[[col]] %||% NA), "")
  rs <- g(sug$snp$column); chr <- .chr_norm(g(sug$chr$column)); pos <- suppressWarnings(as.numeric(gsub(",", "", g(sug$pos$column))))
  ok <- grepl("^rs[0-9]+$", rs) & is.finite(pos); if (sum(ok) < 3) return(NULL)
  b <- tryCatch(.panel_bim(snps = rs[ok]), error = function(e) NULL); if (is.null(b) || nrow(b) < 3) return(NULL)
  i <- match(rs[ok], b$snp); k <- !is.na(i); if (sum(k) < 3) return(NULL)
  same37 <- mean(chr[ok][k] == b$chr[i[k]] & pos[ok][k] == b$pos[i[k]])
  l <- liftover_38_to_37(chr[ok][k], pos[ok][k]); same38 <- mean(!is.na(l$pos) & l$chr == b$chr[i[k]] & l$pos == b$pos[i[k]])
  why <- sprintf("Detected: %d of the first rsIDs in the file sit at the same %s positions in the 1000 Genomes panel.", sum(k), "%s")
  if (same37 >= 0.8) return(list(build = "37", why = sprintf(why, "GRCh37"))) else if (same38 >= 0.8) return(list(build = "38", why = sprintf(why, "GRCh38 (converted)")))
  NULL
}
#' Round 79: rsIDs for a file with positions only (GRCh37 after any conversion), from the panel by chromosome, position and alleles
.rsid_from_panel <- function(d) {
  ok <- !is.na(d$chr) & is.finite(d$pos); if (!any(ok)) return(list(snp = rep(NA_character_, nrow(d)), n = 0L))
  ps <- format(d$pos[ok], scientific = FALSE, trim = TRUE)   # PLINK reads 1e+06 as no position
  b <- .panel_bim(ranges = data.frame(chr = d$chr[ok], s = ps, e = ps, l = seq_len(sum(ok))))
  if (is.null(b)) return(NULL)
  kd <- paste(d$chr, format(d$pos, scientific = FALSE, trim = TRUE)); kb <- paste(b$chr, format(b$pos, scientific = FALSE, trim = TRUE))
  al <- function(x, y) paste(pmin(x, y), pmax(x, y)); m <- match(paste(kd, al(toupper(d$ea), toupper(d$oa))), paste(kb, al(b$a1, b$a2)))
  list(snp = b$snp[m], n = sum(!is.na(m)))
}

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
  # Round 79: rows of other tests (PLINK 2 and REGENIE write covariate rows such as AGE or SEX next to ADD)
  if ("test" %in% names(d)) { tt <- toupper(trimws(d$test)); add <- tt %in% c("ADD", "ADDITIVE")
    if (any(add) && !all(add)) { ch <- c(ch, sprintf("%s rows of other tests (%s) left out; the additive test (ADD) is kept.", fm(sum(!add)), paste(utils::head(unique(tt[!add]), 4), collapse = ", "))); d <- d[add, , drop = FALSE] }
    d$test <- NULL }
  for (r in intersect(c("eaf", "beta", "se", "p", "n", "pos", "z", "ci_lower", "ci_upper"), names(d))) {
    raw <- trimws(d[[r]]); tok <- toupper(raw) %in% IMPORT_NA_TOKENS; v <- .num_parse(raw)   # Round 79: typeset numbers (1.2 × 10−8)
    bad <- sum(is.na(v) & !tok); if (bad) warn <- c(warn, sprintf("%s values in \"%s\" are not numbers and will be dropped.", fm(bad), IMPORT_ROLES[[r]]))
    v[!is.finite(v)] <- NA_real_; d[[r]] <- v
  }
  # Round 79: a cell with several rsIDs (FinnGen rsids "rs1,rs2") keeps the first one
  if ("snp" %in% names(d)) { d$snp <- trimws(d$snp); ml <- grepl("[,;]", d$snp); if (any(ml)) { d$snp[ml] <- sub("\\s*[,;].*$", "", d$snp[ml]); ch <- c(ch, sprintf("%s cells listed several variant IDs; the first one was kept.", fm(sum(ml)))) }
    d$snp[d$snp %in% c(".", "", "NA")] <- NA_character_ }
  if ("oa_alt" %in% names(d)) {   # Round 79: PLINK 2: where A1 is REF, the other allele is ALT
    sw <- toupper(trimws(d$ea)) == toupper(trimws(d$oa)); d$oa[sw] <- d$oa_alt[sw]; d$oa_alt <- NULL
    if (any(sw)) ch <- c(ch, sprintf("For %s where the effect allele A1 is the REF allele, the other allele is ALT (PLINK 2).", vv(sum(sw)))) }
  # Round 79: a "beta" column that holds odds ratios (all positive, centred on 1) would give wrong results without the right answer
  if ("beta" %in% names(d) && identical(facts$scale %||% "beta", "beta")) { b <- d$beta[is.finite(d$beta)]
    if (length(b) >= 5 && all(b > 0) && stats::median(b) > 0.7 && stats::median(b) < 1.4) warn <- c(warn, "Every effect size is positive and close to 1, as odds ratios are. If the column holds odds ratios, choose Odds ratio under \"Effect size is\"; otherwise the results will be wrong.") }
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
  # Round 79: SE from a 95% confidence interval: (upper - lower) / (2 × 1.96), on the log scale for odds ratios
  if (all(c("ci_lower", "ci_upper") %in% names(d)) && "beta" %in% names(d)) {
    if (!("se" %in% names(d))) d$se <- NA_real_
    zc <- stats::qnorm(0.975); orl <- identical(facts$scale, "or")
    sec <- if (orl) suppressWarnings((log(d$ci_upper) - log(d$ci_lower)) / (2 * zc)) else (d$ci_upper - d$ci_lower) / (2 * zc)
    fill <- is.na(d$se) & is.finite(sec) & sec > 0; d$se[fill] <- sec[fill]
    if (any(fill)) ch <- c(ch, sprintf("Standard errors derived from the 95%% confidence interval for %s (%s).", vv(sum(fill)), if (orl) "log(upper / lower) / (2 × 1.96)" else "(upper - lower) / (2 × 1.96)"))
    d$ci_lower <- NULL; d$ci_upper <- NULL }
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
  sem <- if ("se" %in% names(d)) is.na(d$se) else rep(TRUE, nrow(d))   # Round 79: also rows whose SE is missing, not only a missing column
  if (any(sem) && "beta" %in% names(d) && !is.null(lp) && !identical(facts$scale, "z")) {
    zz <- stats::qnorm(lp - log(2), lower.tail = FALSE, log.p = TRUE); se <- abs(d$beta) / zz
    if (!("se" %in% names(d))) d$se <- NA_real_
    oks <- sem & is.finite(se) & se > 0; d$se[oks] <- se[oks]; left <- sem & !oks
    if (any(oks)) ch <- c(ch, sprintf("Standard errors derived from the effect and its P value for %s (Altman and Bland 2011). The P values must be two-sided and from the same test as the effect; rounded P values give approximate standard errors.", vv(sum(oks))))
    if (sum(left)) warn <- c(warn, sprintf("%s with an effect of 0, a P value of 1 or no P value will be dropped (no standard error can be derived).", vv(sum(left))))
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
  if (nr >= max_rows) warn <- c(warn, sprintf("This check looked at the first %s rows only; saving reads the whole file.", format(max_rows, big.mark = ",", scientific = FALSE)))   # Round 99: the whole file is saved
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
  d <- .import_read(path, sep, columns, mp); n_all <- nrow(d)   # Round 99: the row count of the whole file, not of the check
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
  if (!("snp" %in% names(d)) && all(c("chr", "pos") %in% names(d))) {   # Round 79: position-only files get rsIDs where the panel has the variant
    rr <- tryCatch(.rsid_from_panel(d), error = function(e) NULL)
    if (!is.null(rr) && rr$n > 0) { n0 <- nrow(d); d$snp <- rr$snp; d <- d[!is.na(d$snp), , drop = FALSE]; mp$snp <- "rsID from the 1000 Genomes panel"
      tr$changes <- c(tr$changes, sprintf("rsIDs added from the 1000 Genomes panel for %s of %s variants, matched by chromosome, position (GRCh37) and alleles%s.", format(rr$n, big.mark = ","), format(n0, big.mark = ","), if (n0 > rr$n) sprintf("; the %s not in the panel are left out, because datasets are matched by rsID", format(n0 - rr$n, big.mark = ",")) else "")) }
    else tr$changes <- c(tr$changes, "No rsIDs: this dataset can be paired only with datasets that also identify variants by position. Install the 1000 Genomes panel (Settings) to add rsIDs when importing.")
  }
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
  if (length(facts$pmid) && nzchar(as.character(facts$pmid))) { facts$pmid <- sub("^PMID:? *", "", trimws(as.character(facts$pmid)), ignore.case = TRUE)   # Round 109: optional PubMed ID of the GWAS (item 10)
    if (!grepl("^[0-9]{1,9}$", facts$pmid)) stop("The PubMed ID must be digits only (for example 30054458), or empty.", call. = FALSE) } else facts$pmid <- NULL
  rec <- list(trait = nm, role = facts$role, facts = facts, mapping = mp, original = list(file = orig_name, bytes = file.size(orig_path), sha256 = sha256_file(orig_path)),
              rows_in_file = n_all, rows_kept = nrow(d), warnings = v$warnings[!startsWith(v$warnings, "This check looked at the first")], changes = as.list(c(prep_notes, pan, tr$changes)), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
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
#' Round 84: the recipe of a Data library dataset, without reading its table (a complete GWAS has millions of rows: reading one to know
#' its variant count held up the Results page for 10 to 14 s).
library_recipe <- function(trait, root = studies_root()) {
  rf <- file.path(library_dir(root), basename(trait), "recipe.json")
  if (!file.exists(rf)) stop("The dataset \"", trait, "\" is not in the Data library any more.")
  jsonlite::fromJSON(rf, simplifyVector = FALSE)
}
