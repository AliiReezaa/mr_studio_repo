# Public GWAS from OpenGWAS (Round 12): search the catalogue, download independent instruments for an exposure, look the same variants up in an outcome.
# Uses the public packages ieugwasr and TwoSampleMR, exactly as their documentation describes. The access token is passed to them as an
# argument and is held in memory only: it is never written to a file, to the study folders, to scripts, to logs or to the Data library.
# All network calls go through a "backend" list, so the tests can replace the internet with a simulated OpenGWAS.

OGW_POPS <- c(EUR = "European", EAS = "East Asian", SAS = "South Asian", AFR = "African", AMR = "Admixed American")
OGW_POP_LABELS <- c("European" = "EUR", "East Asian" = "EAS", "South Asian" = "SAS", "African American or Afro-Caribbean" = "AFR", "African" = "AFR",
                    "Hispanic or Latin American" = "AMR")

#' OpenGWAS describes ancestry in words. The clumping reference panels use five codes. Anything else (Mixed, NR, Middle Eastern) has no panel.
ogw_pop_code <- function(label) {
  l <- as.character(label %||% NA_character_)[1]
  if (is.na(l) || !nzchar(l)) return(NA_character_)
  if (l %in% names(OGW_POPS)) return(l)
  if (l %in% names(OGW_POP_LABELS)) return(unname(OGW_POP_LABELS[[l]]))
  NA_character_
}

# ---- the access token ----------------------------------------------------------------------------------------------------------------
#' Tidy what the user pasted (spaces, quotes, "Bearer ") and check that it has the shape of a token. Returns NULL when it does not.
ogw_clean_token <- function(x) {
  x <- trimws(as.character(x %||% "")[1], whitespace = "[[:space:]\"']")
  x <- sub("^bearer\\s+", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]\"']", "", x)
  if (grepl("^[A-Za-z0-9_-]{10,}\\.[A-Za-z0-9_-]{10,}\\.[A-Za-z0-9_-]{10,}$", x)) x else NULL
}
#' What a token says about itself without contacting anyone: its expiry date. Nothing else is read or returned (the account is not shown).
ogw_token_info <- function(tok, now = Sys.time()) {
  out <- list(set = FALSE, valid_shape = FALSE, expires = NA_character_, days_left = NA_real_, expired = NA)
  if (is.null(tok) || !nzchar(tok)) return(out)
  out$set <- TRUE
  parts <- strsplit(tok, ".", fixed = TRUE)[[1]]; if (length(parts) != 3L) return(out)
  out$valid_shape <- TRUE
  pl <- tryCatch({
    x <- chartr("-_", "+/", parts[2]); x <- paste0(x, strrep("=", (4 - nchar(x) %% 4) %% 4))
    jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(x)))
  }, error = function(e) NULL)
  ex <- suppressWarnings(as.numeric(pl$exp))
  if (length(ex) == 1L && is.finite(ex)) {
    t <- as.POSIXct(ex, origin = "1970-01-01", tz = "UTC")
    out$expires <- format(t, "%Y-%m-%d"); out$days_left <- round(as.numeric(difftime(t, now, units = "days")), 1); out$expired <- out$days_left < 0
  }
  out
}

#' The token lives here, in memory, for as long as the engine runs. There is deliberately no function that writes it anywhere.
.ogw_mem <- new.env(parent = emptyenv())
ogw_token_set <- function(tok) { if (is.null(tok)) rm(list = intersect("jwt", ls(.ogw_mem)), envir = .ogw_mem) else assign("jwt", tok, envir = .ogw_mem); invisible(TRUE) }
ogw_token_get <- function() if (exists("jwt", envir = .ogw_mem, inherits = FALSE)) get("jwt", envir = .ogw_mem) else NULL

# ---- certificate file ------------------------------------------------------------------------------------------------------------------
# The bundled R is built on a different computer, and its libcurl remembers that computer's certificate path
# (for example /Users/runner/work/_temp/build_runtime/ssl/cacert.pem). On the user's machine that file does not exist and every
# https request fails with "Problem with the SSL CA cert". So we choose a certificate file that really exists here, in this order:
# the user's own setting, the file inside the bundled runtime, the operating system's file, and finally a copy shipped with the
# app (engine/certs/cacert.pem). If a request still fails with a certificate error, the next file in the list is tried.
.ogw_ca_state <- new.env(parent = emptyenv())
.ogw_ca_valid <- function(f) {
  if (!length(f) || is.na(f) || !nzchar(f) || !file.exists(f) || dir.exists(f)) return(FALSE)
  if (isTRUE(file.info(f)$size < 1000)) return(FALSE)
  l <- tryCatch(readLines(f, n = 60L, warn = FALSE), error = function(e) character())
  any(grepl("BEGIN CERTIFICATE", l, fixed = TRUE))
}
ogw_ca_candidates <- function() {
  e <- Sys.getenv(c("MRGUI_CA_BUNDLE", "CURL_CA_BUNDLE", "SSL_CERT_FILE", "REQUESTS_CA_BUNDLE"), unset = "")
  rh <- tryCatch(R.home(), error = function(e) "")
  ed <- getOption("mrgui.engine_dir", "")
  unique(c(e[nzchar(e)],
    file.path(rh, "..", "..", "ssl", "cacert.pem"), file.path(rh, "..", "..", "ssl", "cert.pem"),
    file.path(rh, "..", "..", "Library", "ssl", "cacert.pem"), file.path(rh, "..", "..", "etc", "ssl", "cacert.pem"),
    file.path(rh, "..", "..", "..", "ssl", "cacert.pem"), file.path(rh, "..", "..", "..", "Library", "ssl", "cacert.pem"),
    "/etc/ssl/cert.pem", "/private/etc/ssl/cert.pem", "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/tls/certs/ca-bundle.crt", "/etc/ssl/ca-bundle.pem", "/etc/pki/tls/cacert.pem",
    if (nzchar(ed)) file.path(ed, "certs", "cacert.pem")))
}
.ogw_ca_apply <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  Sys.setenv(CURL_CA_BUNDLE = path, SSL_CERT_FILE = path)
  if (requireNamespace("httr", quietly = TRUE)) try(httr::set_config(httr::config(cainfo = path)), silent = TRUE)
  .ogw_ca_state$path <- path
  path
}
ogw_ca_setup <- function(candidates = ogw_ca_candidates(), force = FALSE) {
  if (!force && !is.null(.ogw_ca_state$path)) return(invisible(.ogw_ca_state$path))
  ok <- Filter(.ogw_ca_valid, candidates)
  .ogw_ca_state$list <- vapply(ok, normalizePath, "", winslash = "/", mustWork = FALSE, USE.NAMES = FALSE)
  if (!length(ok)) { .ogw_ca_state$path <- ""; return(invisible("")) }
  invisible(.ogw_ca_apply(ok[1]))
}
#' Move to the next usable certificate file (used after a certificate error). FALSE when there is none left.
ogw_ca_next <- function() {
  ogw_ca_setup(); l <- .ogw_ca_state$list %||% character(); cur <- .ogw_ca_state$path %||% ""
  i <- match(cur, l); nx <- if (is.na(i)) l[1] else l[i + 1L]
  if (is.na(nx) || !length(nx)) return(FALSE)
  .ogw_ca_apply(nx); TRUE
}
.ogw_is_cert_error <- function(e) grepl("SSL|certificate|CA cert|trust anchor|peer verification", conditionMessage(e), ignore.case = TRUE)

# ---- talking to OpenGWAS (replaceable for tests) --------------------------------------------------------------------------------------
.ogw_real_backend <- function() { ogw_ca_setup(); list(
  info = function(jwt, id = NULL) as.data.frame(ieugwasr::gwasinfo(id, opengwas_jwt = jwt), stringsAsFactors = FALSE),
  tophits = function(id, p, r2, kb, pop, jwt) as.data.frame(ieugwasr::tophits(id, pval = p, clump = 1, r2 = r2, kb = kb, pop = pop, opengwas_jwt = jwt), stringsAsFactors = FALSE),
  outcome = function(snps, id, proxies, rsq, maf, jwt) as.data.frame(suppressMessages(TwoSampleMR::extract_outcome_data(snps, id, proxies = proxies, rsq = rsq, maf_threshold = maf, opengwas_jwt = jwt)), stringsAsFactors = FALSE)) }
ogw_backend <- function() getOption("mrgui.ogw_backend", .ogw_real_backend())

#' Turn the packages' technical errors into sentences that say what to do. Tokens are removed from any text that is passed on.
ogw_friendly <- function(e) {
  m <- paste(conditionMessage(e), collapse = " "); m <- gsub("\\s+", " ", redact(m, character()))
  detail <- paste0(" Technical detail: ", substr(m, 1, 220))
  if (grepl("Status code from OpenGWAS API: 401|\\b401\\b|must provide a token|token.*(expired|invalid)", m, ignore.case = TRUE))
    return("OpenGWAS did not accept the access token. It may be missing, mistyped or expired (tokens last for a limited time). Create a new token at api.opengwas.io and paste it in Settings.")
  if (grepl("\\b429\\b|too many|allowance|rate limit", m, ignore.case = TRUE))
    return("OpenGWAS has limited how many requests this token can make for now. Nothing was lost. Wait a while and try again.")
  if (grepl("\\b403\\b|forbidden", m, ignore.case = TRUE)) return("OpenGWAS refused access to this dataset with your token.")
  if (grepl("\\b404\\b|not found", m, ignore.case = TRUE)) return("OpenGWAS does not know this dataset ID. Check the ID, or update the catalogue.")
  if (grepl("\\b5[0-9][0-9]\\b|bad gateway|unavailable|server error", m, ignore.case = TRUE))
    return("OpenGWAS is not answering properly at the moment (a server error on their side). Try again later. Imported files and saved studies still work.")
  if (grepl("resolve host|could not resolve|timed out|timeout|failed to connect|connection (refused|reset)|ssl|certificate|curl|network|unreachable|cannot open the connection", m, ignore.case = TRUE))
    return(paste0("Cannot reach OpenGWAS. Check the internet connection (or a firewall or proxy) and try again. Imported files and saved studies still work.", detail))
  paste0("OpenGWAS returned an unexpected answer: ", substr(m, 1, 240))
}
.ogw_do <- function(expr) {
  ex <- substitute(expr); env <- parent.frame()
  repeat {
    r <- tryCatch(list(v = eval(ex, env)), error = function(e) e)
    if (!inherits(r, "error")) return(r$v)
    if (.ogw_is_cert_error(r) && ogw_ca_next()) next
    stop(ogw_friendly(r), call. = FALSE)
  }
}
.ogw_need <- function(jwt, offline = FALSE) {
  if (isTRUE(offline)) stop("Offline mode is on in Settings, so nothing is downloaded. Turn it off to use OpenGWAS.", call. = FALSE)
  if (is.null(jwt) || !nzchar(jwt)) stop("Paste your OpenGWAS access token in Settings first. It is free: create it at api.opengwas.io.", call. = FALSE)
  ti <- ogw_token_info(jwt); if (isTRUE(ti$expired)) stop("The OpenGWAS access token expired on ", ti$expires, ". Create a new one at api.opengwas.io and paste it in Settings.", call. = FALSE)
  invisible(TRUE)
}

# ---- the catalogue ---------------------------------------------------------------------------------------------------------------------
OGW_COLS_CHR <- c("id", "trait", "population", "sex", "build", "author", "consortium", "category", "subcategory", "unit", "note", "pmid")
OGW_COLS_NUM <- c("ncase", "ncontrol", "sample_size", "nsnp", "year", "mr")
ogw_cache_file <- function() file.path(settings_dir(), "opengwas_catalogue.rds")
#' Keep only what the search shows (no personal data is in the catalogue), with every column present and typed.
ogw_tidy_catalogue <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  for (k in OGW_COLS_CHR) d[[k]] <- if (k %in% names(d)) as.character(d[[k]]) else NA_character_
  for (k in OGW_COLS_NUM) d[[k]] <- if (k %in% names(d)) suppressWarnings(as.numeric(d[[k]])) else NA_real_
  d <- d[!is.na(d$id) & nzchar(d$id), c(OGW_COLS_CHR, OGW_COLS_NUM), drop = FALSE]
  d <- d[!duplicated(d$id), , drop = FALSE]; rownames(d) <- NULL; d
}
ogw_catalogue_status <- function() {
  f <- ogw_cache_file(); if (!file.exists(f)) return(list(have = FALSE))
  x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x)) return(list(have = FALSE))
  list(have = TRUE, n = nrow(x$data), fetched = x$fetched, age_days = round(as.numeric(difftime(Sys.time(), as.POSIXct(x$fetched, tz = "UTC"), units = "days")), 1))
}
ogw_catalogue_read <- function() { f <- ogw_cache_file(); if (!file.exists(f)) return(NULL); x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x)) NULL else x$data }
#' Download the list of all public datasets (about 50 000 rows) and keep it on this computer, so searching is instant and works offline.
ogw_catalogue_update <- function(jwt, offline = FALSE, backend = ogw_backend()) {
  .ogw_need(jwt, offline)
  d <- .ogw_do(backend$info(jwt))
  if (!is.data.frame(d) || !nrow(d)) stop("OpenGWAS returned an empty catalogue. Try again later.", call. = FALSE)
  d <- ogw_tidy_catalogue(d)
  dir.create(settings_dir(), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(ogw_cache_file(), ".tmp"); saveRDS(list(fetched = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), data = d), tmp); file.rename(tmp, ogw_cache_file())
  ogw_catalogue_status()
}
ogw_kind <- function(ncase) ifelse(is.finite(ncase) & ncase > 0, "binary", "continuous")

#' Search the saved catalogue. Every word typed must appear in the trait, ID, author, consortium or note. Best matches first, then larger studies.
ogw_search <- function(cat, query, population = "", min_n = 0, mr_only = TRUE, kind = "any", limit = 40L) {
  q <- tolower(trimws(query %||% "")); words <- strsplit(q, "[[:space:]]+")[[1]]; words <- words[nzchar(words)]
  if (is.null(cat)) return(list(rows = data.frame(), total = 0L))
  if (!length(words)) return(list(rows = cat[0, , drop = FALSE], total = 0L))
  hay <- tolower(paste(cat$id, cat$trait, cat$author, cat$consortium, cat$note, cat$pmid))
  ok <- rep(TRUE, nrow(cat)); for (w in words) ok <- ok & grepl(w, hay, fixed = TRUE)
  if (nzchar(population %||% "")) { map <- c(OGW_POP_LABELS, stats::setNames(names(OGW_POPS), names(OGW_POPS))); ok <- ok & (unname(map[cat$population]) %in% population) }
  if (isTRUE(mr_only)) ok <- ok & (cat$mr %in% 1)
  if (is.finite(min_n) && min_n > 0) ok <- ok & (is.finite(cat$sample_size) & cat$sample_size >= min_n)
  if (kind %in% c("binary", "continuous")) ok <- ok & (ogw_kind(cat$ncase) == kind)
  r <- cat[ok, , drop = FALSE]; total <- nrow(r)
  if (total) {
    tl <- tolower(r$trait %||% ""); score <- 4 * (r$id == query | tolower(r$id) == q) + 3 * (tl == q) + 2 * startsWith(tl, q)
    r <- r[order(-score, -ifelse(is.finite(r$sample_size), r$sample_size, 0)), , drop = FALSE]
  }
  list(rows = utils::head(r, limit), total = total)
}
ogw_search_payload <- function(res) {
  r <- res$rows
  list(total = res$total, rows = unname(lapply(seq_len(nrow(r)), function(i) {
    x <- as.list(r[i, , drop = FALSE]); x$kind <- ogw_kind(x$ncase); x$pop <- ogw_pop_code(x$population); x[!vapply(x, function(v) length(v) == 1L && is.na(v), TRUE)] })))
}

# ---- downloads -------------------------------------------------------------------------------------------------------------------------
.ogw_label <- function(trait, id, lim = 60L) { suffix <- paste0(" (", id, ")"); clean_study_name(paste0(substr(trait, 1L, max(8L, lim - nchar(suffix))), suffix)) }
.ogw_build <- function(b) if (grepl("38", b %||% "")) "38" else "37"
.ogw_n <- function(x) if (length(x) == 1L && is.finite(suppressWarnings(as.numeric(x)))) as.numeric(x) else ""
.ogw_is_ukb <- function(meta) grepl("^ukb-", meta$id %||% "") || grepl("uk biobank|ukbiobank", paste(meta$consortium, meta$author, meta$note), ignore.case = TRUE)
.ogw_versions <- function() list(ieugwasr = as.character(utils::packageVersion("ieugwasr")), TwoSampleMR = as.character(utils::packageVersion("TwoSampleMR")))

#' One row of catalogue information for a dataset ID (from the saved catalogue when it is there, else asked from OpenGWAS).
ogw_meta <- function(id, jwt, backend = ogw_backend()) {
  cat <- ogw_catalogue_read(); m <- if (!is.null(cat)) cat[cat$id == id, , drop = FALSE] else NULL
  if (is.null(m) || !nrow(m)) { d <- .ogw_do(backend$info(jwt, id)); m <- if (is.data.frame(d) && nrow(d)) ogw_tidy_catalogue(d) else NULL }
  if (is.null(m) || !nrow(m)) stop("OpenGWAS does not know a dataset with the ID \"", id, "\".", call. = FALSE)
  as.list(m[1, , drop = FALSE])
}

.ogw_clean_rows <- function(d, what) {
  n0 <- nrow(d); d$ea <- toupper(as.character(d$ea)); d$oa <- toupper(as.character(d$oa))
  bad_allele <- !(grepl("^[ACGT]$", d$ea) & grepl("^[ACGT]$", d$oa)); bad_num <- !(is.finite(d$beta) & is.finite(d$se) & d$se > 0)
  no_id <- is.na(d$snp) | !nzchar(d$snp); dup <- duplicated(d$snp)
  keep <- !bad_allele & !bad_num & !no_id & !dup
  w <- character()
  if (any(bad_allele)) w <- c(w, sprintf("%d %s dropped: not single-letter SNPs (insertions or deletions).", sum(bad_allele), if (sum(bad_allele) == 1) "variant was" else "variants were"))
  if (any(bad_num & !bad_allele)) w <- c(w, sprintf("%d %s dropped: missing or non-positive standard error or effect.", sum(bad_num & !bad_allele), if (sum(bad_num & !bad_allele) == 1) "variant was" else "variants were"))
  if (any(dup & !bad_allele & !bad_num)) w <- c(w, sprintf("%d repeated %s dropped.", sum(dup & !bad_allele & !bad_num), if (sum(dup & !bad_allele & !bad_num) == 1) "variant was" else "variants were"))
  list(d = d[keep, , drop = FALSE], warnings = w, dropped = n0 - sum(keep))
}
.ogw_write <- function(name, d, rec, root) {
  dest <- file.path(library_dir(root), name)
  p <- study_path_problem(library_dir(root)); if (!is.null(p)) stop(p, call. = FALSE)
  dir.create(dest, recursive = TRUE)
  out <- file.path(dest, "standardised.csv"); utils::write.csv(d, out, row.names = FALSE, na = "")
  rec$original$bytes <- file.size(out); rec$original$sha256 <- sha256_file(out)
  atomic_write_text(file.path(dest, "recipe.json"), jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  rec
}
.ogw_same <- function(a, b) isTRUE(all.equal(unclass(a), unclass(b)))

#' Independent instruments for an exposure, clumped on the OpenGWAS server against the reference panel of the chosen ancestry.
ogw_get_exposure <- function(id, jwt, p = 5e-8, r2 = 0.001, kb = 10000, pop = NULL, name = NULL, offline = FALSE, root = studies_root(), backend = ogw_backend()) {
  .ogw_need(jwt, offline)
  if (!is.finite(p) || p < 1e-30 || p > 1e-3) stop("The instrument P-value threshold must be between 1e-30 and 0.001.", call. = FALSE)
  if (!is.finite(r2) || r2 <= 0 || r2 > 0.5) stop("The clumping r² must be above 0 and at most 0.5.", call. = FALSE)
  if (!is.finite(kb) || kb < 10 || kb > 100000) stop("The clumping window must be between 10 and 100,000 kb.", call. = FALSE)
  meta <- ogw_meta(id, jwt, backend)
  pop <- pop %||% ogw_pop_code(meta$population)
  if (is.null(pop) || is.na(pop) || !pop %in% names(OGW_POPS)) stop("The ancestry of this dataset (\"", meta$population %||% "not reported", "\") has no clumping reference panel. Choose the reference population yourself.", call. = FALSE)
  nm <- clean_study_name(if (is.null(name) || !nzchar(name)) .ogw_label(meta$trait, id) else name); if (!nzchar(nm)) stop("Give the dataset a name.", call. = FALSE)
  params <- list(pval = p, clump = 1, r2 = r2, kb = kb, pop = pop)
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) {
    old <- tryCatch(jsonlite::fromJSON(file.path(dest, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(old) && identical(old$source$id, id) && identical(old$source$role, "exposure") && .ogw_same(old$source$params, params)) { old$reused <- TRUE; return(old) }
    stop("A dataset named \"", nm, "\" is already in the Data library, with other settings or another source. Choose another name.", call. = FALSE)
  }
  d <- .ogw_do(backend$tophits(id, p, r2, kb, pop, jwt))
  if (!is.data.frame(d) || !nrow(d)) stop("No variants reach P < ", format(p, scientific = TRUE), " in this dataset after clumping. A less strict threshold (for example 5e-6) would find some, but they are weaker instruments.", call. = FALSE)
  std <- data.frame(snp = as.character(d$rsid), chr = as.character(d$chr), pos = suppressWarnings(as.numeric(d$position)), ea = d$ea, oa = d$nea, eaf = suppressWarnings(as.numeric(d$eaf)),
                    beta = suppressWarnings(as.numeric(d$beta)), se = suppressWarnings(as.numeric(d$se)), p = suppressWarnings(as.numeric(d$p)), n = suppressWarnings(as.numeric(d$n)), stringsAsFactors = FALSE)
  std$n[!is.finite(std$n)] <- meta$sample_size
  cl <- .ogw_clean_rows(std, "exposure"); std <- cl$d; std <- std[order(std$p), , drop = FALSE]; rownames(std) <- NULL
  w <- cl$warnings
  if (nrow(std) < 2L) stop("Fewer than 2 usable instruments remained (", nrow(std), "). A less strict threshold would find more, but weaker, instruments.", call. = FALSE)
  if (nrow(std) < 10L) w <- c(w, sprintf("Only %d instruments. MR-Egger, the median methods and the outlier detectors need more variants to be informative.", nrow(std)))
  if (!is.finite(meta$sample_size)) w <- c(w, "The catalogue does not give a sample size for this dataset.")
  rec <- list(trait = nm, role = "exposure",
              facts = list(name = nm, role = "exposure", kind = ogw_kind(meta$ncase), scale = "beta", ea_ok = "yes", frq = "eaf", build = .ogw_build(meta$build), ancestry = pop, n = .ogw_n(meta$sample_size),
                           scope = "selected", independent = "yes", ukb = .ogw_is_ukb(meta)),
              mapping = list(snp = "rsid", chr = "chr", pos = "position", ea = "ea", oa = "nea", eaf = "eaf", beta = "beta", se = "se", p = "p", n = "n"),
              original = list(file = paste0("OpenGWAS:", id)), rows_in_file = nrow(d), rows_kept = nrow(std), warnings = as.list(w), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              usable_as_exposure = TRUE, usable_as_outcome = FALSE, usable_in_runs = TRUE,
              source = list(type = "opengwas", role = "exposure", id = id, trait = meta$trait, population = meta$population, sample_size = meta$sample_size, ncase = meta$ncase, ncontrol = meta$ncontrol,
                            consortium = meta$consortium, author = meta$author, year = meta$year, pmid = meta$pmid, unit = meta$unit, params = params, downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                            packages = .ogw_versions(), note = "Instruments clumped on the OpenGWAS server (ieugwasr::tophits) against the 1000 Genomes reference panel of the chosen population."))
  .ogw_write(nm, std, rec, root)
}

#' The same variants looked up in an outcome. Variants the outcome lacks are replaced by proxies (r² at least `rsq`) when asked.
ogw_get_outcome <- function(id, exposure_name, jwt, proxies = TRUE, rsq = 0.8, maf = 0.3, name = NULL, offline = FALSE, root = studies_root(), backend = ogw_backend()) {
  .ogw_need(jwt, offline)
  if (!is.finite(rsq) || rsq < 0.5 || rsq > 1) stop("The proxy r² must be between 0.5 and 1.", call. = FALSE)
  if (!is.finite(maf) || maf <= 0 || maf > 0.5) stop("The allele frequency limit for palindromic proxies must be above 0 and at most 0.5.", call. = FALSE)
  ex <- library_dataset(exposure_name, root)
  if (!isTRUE(ex$recipe$usable_as_exposure)) stop("\"", exposure_name, "\" is not a set of independent instruments, so there is nothing to look up.", call. = FALSE)
  snps <- unique(ex$data$snp[!is.na(ex$data$snp) & grepl("^rs[0-9]+$", ex$data$snp)])
  if (length(snps) < 2L) stop("The exposure \"", exposure_name, "\" has no rsID variants to look up. Looking up by position is not built.", call. = FALSE)
  meta <- ogw_meta(id, jwt, backend)
  ex_short <- clean_study_name(substr(if (!is.null(ex$recipe$source$id)) ex$recipe$source$id else exposure_name, 1L, 18L))
  nm <- clean_study_name(if (is.null(name) || !nzchar(name)) paste0(.ogw_label(meta$trait, id, 38L), " for ", ex_short) else name); if (!nzchar(nm)) stop("Give the dataset a name.", call. = FALSE)
  params <- list(proxies = isTRUE(proxies), rsq = rsq, maf_threshold = maf, exposure = exposure_name, exposure_sha256 = ex$recipe$original$sha256 %||% "")
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) {
    old <- tryCatch(jsonlite::fromJSON(file.path(dest, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(old) && identical(old$source$id, id) && identical(old$source$role, "outcome") && .ogw_same(old$source$params, params)) { old$reused <- TRUE; return(old) }
    stop("A dataset named \"", nm, "\" is already in the Data library, with other settings or another source. Choose another name.", call. = FALSE)
  }
  o <- .ogw_do(backend$outcome(snps, id, isTRUE(proxies), rsq, maf, jwt))
  if (!is.data.frame(o) || !nrow(o)) stop("None of the ", length(snps), " instruments was found in this outcome", if (isTRUE(proxies)) ", not even through proxies" else " (proxies were off)", ". Check that the outcome covers the same genome build and ancestry.", call. = FALSE)
  col <- function(k) if (k %in% names(o)) o[[k]] else rep(NA, nrow(o))
  tgt <- as.character(col("target_snp.outcome")); key <- ifelse(!is.na(tgt) & nzchar(tgt), tgt, as.character(col("SNP")))
  isprox <- as.logical(col("proxy.outcome")) %in% TRUE
  std <- data.frame(snp = key, chr = as.character(col("chr")), pos = suppressWarnings(as.numeric(col("pos"))), ea = col("effect_allele.outcome"), oa = col("other_allele.outcome"),
                    eaf = suppressWarnings(as.numeric(col("eaf.outcome"))), beta = suppressWarnings(as.numeric(col("beta.outcome"))), se = suppressWarnings(as.numeric(col("se.outcome"))),
                    p = suppressWarnings(as.numeric(col("pval.outcome"))), n = suppressWarnings(as.numeric(col("samplesize.outcome"))),
                    proxy_snp = ifelse(isprox, as.character(col("proxy_snp.outcome")), ""), stringsAsFactors = FALSE)
  std$n[!is.finite(std$n)] <- meta$sample_size
  # a direct hit is preferred over a proxy for the same instrument
  std <- std[order(std$proxy_snp != ""), , drop = FALSE]
  cl <- .ogw_clean_rows(std, "outcome"); std <- cl$d; rownames(std) <- NULL
  w <- cl$warnings
  found <- sum(std$snp %in% snps); n_prox <- sum(std$proxy_snp != ""); missing <- setdiff(snps, std$snp)
  if (found < 2L) stop("Only ", found, " of ", length(snps), " instruments were found in this outcome. At least 2 are needed.", call. = FALSE)
  if (length(missing)) w <- c(w, sprintf("%d of %d instruments are not in this outcome%s and are lost: %s.", length(missing), length(snps), if (isTRUE(proxies)) " and had no proxy" else " (proxies were off)", paste(utils::head(missing, 8), collapse = ", ")))
  if (n_prox) w <- c(w, sprintf("%d instruments were replaced by a proxy variant (r² at least %s). The proxy's own ID is in the column proxy_snp.", n_prox, format(rsq)))
  if (!is.finite(meta$sample_size)) w <- c(w, "The catalogue does not give a sample size for this dataset.")
  rec <- list(trait = nm, role = "outcome",
              facts = list(name = nm, role = "outcome", kind = ogw_kind(meta$ncase), scale = "beta", ea_ok = "yes", frq = "eaf", build = .ogw_build(meta$build), ancestry = ogw_pop_code(meta$population) %||% (meta$population %||% "?"), n = .ogw_n(meta$sample_size),
                           scope = "selected", independent = "no", ukb = .ogw_is_ukb(meta)),
              mapping = list(snp = "rsid", chr = "chr", pos = "position", ea = "ea", oa = "nea", eaf = "eaf", beta = "beta", se = "se", p = "p", n = "n"),
              original = list(file = paste0("OpenGWAS:", id)), rows_in_file = nrow(o), rows_kept = nrow(std), warnings = as.list(w), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE,
              lookup = list(requested = length(snps), found = found, direct = found - n_prox, by_proxy = n_prox, missing = length(missing)),
              source = list(type = "opengwas", role = "outcome", id = id, trait = meta$trait, population = meta$population, sample_size = meta$sample_size, ncase = meta$ncase, ncontrol = meta$ncontrol,
                            consortium = meta$consortium, author = meta$author, year = meta$year, pmid = meta$pmid, unit = meta$unit, params = params, for_exposure = exposure_name, downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                            packages = .ogw_versions(), note = "Only the instruments of the exposure named in for_exposure were looked up (TwoSampleMR::extract_outcome_data, with proxies when proxies is true)."))
  .ogw_write(nm, std, rec, root)
}

#' Ask OpenGWAS about one small, well-known dataset to find out whether it accepts the token.
ogw_token_test <- function(jwt, offline = FALSE, backend = ogw_backend()) {
  .ogw_need(jwt, offline); d <- .ogw_do(backend$info(jwt, "ieu-a-2"))
  if (!is.data.frame(d) || !nrow(d)) stop("OpenGWAS answered, but with nothing. Try again later.", call. = FALSE)
  TRUE
}

#' Is OpenGWAS reachable at all from this engine? Uses its public status page, which needs no token. Any answer from the server counts as reachable.
.ogw_ping_real <- function() {
  ca <- ogw_ca_setup()
  repeat {
    r <- tryCatch(httr::GET("https://api.opengwas.io/api/status", httr::timeout(20)), error = function(e) e)
    if (inherits(r, "error") && .ogw_is_cert_error(r) && ogw_ca_next()) { ca <- .ogw_ca_state$path; next }
    break
  }
  if (inherits(r, "error")) return(list(ok = FALSE, detail = substr(paste0(gsub("\\s+", " ", redact(conditionMessage(r), character())), " [certificates: ", if (nzchar(ca)) ca else "none found", "]"), 1, 400)))
  list(ok = TRUE, detail = paste0("api.opengwas.io answered (HTTP ", r$status_code, "). Certificates: ", if (nzchar(ca)) ca else "none found"))
}
ogw_ping <- function(backend = ogw_backend()) { f <- backend$ping %||% .ogw_ping_real; tryCatch(f(), error = function(e) list(ok = FALSE, detail = substr(conditionMessage(e), 1, 300))) }
