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

#' Ancestry as stored for a download: the panel code where one exists, otherwise OpenGWAS's own word (Round 72c: "Mixed" was shown as "?").
.ogw_ancestry <- function(label) { pc <- ogw_pop_code(label); if (is.na(pc)) { l <- as.character(label %||% NA_character_)[1]; if (is.na(l) || !nzchar(l)) "not reported" else l } else pc }

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
  tophits = function(id, p, r2, kb, pop, jwt, clump = 1, force = FALSE) as.data.frame(ieugwasr::tophits(id, pval = p, clump = clump, r2 = r2, kb = kb, pop = pop, force_server = force, opengwas_jwt = jwt), stringsAsFactors = FALSE),
  outcome = function(snps, id, proxies, rsq, maf, jwt) as.data.frame(suppressMessages(TwoSampleMR::extract_outcome_data(snps, id, proxies = proxies, rsq = rsq, maf_threshold = maf, opengwas_jwt = jwt)), stringsAsFactors = FALSE),
  # Round 72c: joint clumping on the OpenGWAS server (multivariable sets when no 1000 Genomes panel is installed); returns the kept rsIDs
  clump = function(snps, p, r2, kb, pop, jwt) as.character(ieugwasr::ld_clump(data.frame(rsid = snps, pval = p, stringsAsFactors = FALSE), clump_kb = kb, clump_r2 = r2, pop = pop, opengwas_jwt = jwt)$rsid),
  # Round 49 (FR-106, FR-112): other-trait associations of flagged variants, and their nearest gene from Ensembl VEP (GRCh37, the build of OpenGWAS).
  files = function(id, jwt) ieugwasr::gwasinfo_files(id, opengwas_jwt = jwt),
  phewas = function(snps, p, jwt) as.data.frame(ieugwasr::phewas(snps, pval = p, opengwas_jwt = jwt), stringsAsFactors = FALSE),
  genes = function(rsids) { r <- httr::POST("https://grch37.rest.ensembl.org/vep/human/id", body = list(ids = as.list(rsids)), encode = "json", query = list(nearest = "symbol", pick = 1),
      httr::add_headers(Accept = "application/json"), httr::timeout(60)); httr::stop_for_status(r)
    x <- httr::content(r, as = "parsed", simplifyVector = FALSE)
    data.frame(rsid = vapply(x, function(v) as.character(v$id %||% NA_character_), ""), chr = vapply(x, function(v) as.character(v$seq_region_name %||% NA_character_), ""),
               pos = vapply(x, function(v) as.numeric(v$start %||% NA_real_), 0),
               gene = vapply(x, function(v) { tc <- v$transcript_consequences; ig <- v$intergenic_consequences
                 g <- if (length(tc)) tc[[1]]$gene_symbol %||% unlist(tc[[1]]$nearest)[1] else if (length(ig)) unlist(ig[[1]]$nearest)[1] else NULL; if (is.null(g)) NA_character_ else as.character(g) }, ""),
               stringsAsFactors = FALSE) |> .ogw_intergenic_genes() }) }
#' Round 49: VEP names no gene for an intergenic variant; the closest protein-coding gene within 500 kb (then 2 Mb) is taken from Ensembl's
#' overlap endpoint (distance to the nearer end of the gene), and marked with its distance.
.ogw_intergenic_genes <- function(d) {
  for (i in which(is.na(d$gene) & !is.na(d$chr) & is.finite(d$pos))) for (w in c(5e5, 2e6)) {
    r <- tryCatch(httr::GET(sprintf("https://grch37.rest.ensembl.org/overlap/region/human/%s:%d-%d", d$chr[i], max(1, round(d$pos[i] - w)), round(d$pos[i] + w)),
                            query = list(feature = "gene", biotype = "protein_coding"), httr::add_headers(Accept = "application/json"), httr::timeout(30)), error = function(e) NULL)
    if (is.null(r) || httr::status_code(r) != 200) break
    g <- httr::content(r, as = "parsed")
    if (!length(g)) { if (w == 2e6) d$gene[i] <- "no protein-coding gene within 2 Mb"; next }
    dist <- vapply(g, function(x) max(0, x$start - d$pos[i], d$pos[i] - x$end), 0); k <- which.min(dist)
    nm <- g[[k]]$external_name %||% g[[k]]$id; d$gene[i] <- sprintf("%s (%s kb away)", nm, format(round(dist[k] / 1000))); break
  }
  d
}
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
  out <- lib_table_write(d, dest, "standardised")
  rec$original$bytes <- file.size(out); rec$original$sha256 <- sha256_file(out)
  atomic_write_text(file.path(dest, "recipe.json"), jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null"))
  rec
}
.ogw_same <- function(a, b) isTRUE(all.equal(unclass(a), unclass(b)))

#' Independent instruments for an exposure, clumped against the 1000 Genomes reference panel of the chosen ancestry, either on the
#' OpenGWAS server (where = "server", ieugwasr::tophits with clump = 1) or on this computer (where = "local", Round 22): every variant
#' below the threshold is downloaded unclumped (tophits with clump = 0) and clumped at once with the bundled PLINK (ld_clump_here, the
#' same call as ieugwasr::ld_clump_local). Clumping at download keeps outcome lookups as small as with server clumping. The full
#' candidate list and the status of every candidate are kept in the dataset folder (Parquet tables candidates and instrument_selection, see lib_table_write).
#' Round 25: where = "server_live" asks the server to clump now (tophits force_server = TRUE). With the default settings (5e-8, 0.001,
#' 10,000 kb) a plain server request ("server") returns OpenGWAS's stored, pre-computed clumps instead (tophits sets preclumped = 1);
#' on BMI (ieu-b-40) those kept a weaker lead variant in 183 regions, while live and local clumping gave the same 510 instruments.
ogw_stored_clumps <- function(params) is.null(params$where) && isTRUE(params$clump == 1) && isTRUE(params$pval == 5e-8) && isTRUE(params$r2 == 0.001) && isTRUE(params$kb == 10000)
ogw_get_exposure <- function(id, jwt, p = 5e-8, r2 = 0.001, kb = 10000, pop = NULL, name = NULL, offline = FALSE, root = studies_root(), backend = ogw_backend(), where = "server", replace = "") {
  where <- if (identical(where, "local")) "local" else if (identical(where, "server_live")) "server_live" else "server"
  .ogw_need(jwt, offline)
  if (!is.finite(p) || p < 1e-30 || p > 1e-3) stop("The instrument P-value threshold must be between 1e-30 and 0.001.", call. = FALSE)
  if (!is.finite(r2) || r2 <= 0 || r2 > 0.5) stop("The clumping r² must be above 0 and at most 0.5.", call. = FALSE)
  if (!is.finite(kb) || kb < 10 || kb > 100000) stop("The clumping window must be between 10 and 100,000 kb.", call. = FALSE)
  meta <- ogw_meta(id, jwt, backend)
  pop <- pop %||% ogw_pop_code(meta$population)
  if (is.null(pop) || is.na(pop) || !pop %in% names(OGW_POPS)) stop("The ancestry of this dataset (\"", meta$population %||% "not reported", "\") has no clumping reference panel. Choose the reference population yourself.", call. = FALSE)
  nm <- clean_study_name(if (is.null(name) || !nzchar(name)) paste0(.ogw_label(meta$trait, id, if (where == "server") 60L else 48L), switch(where, local = " clumped here", server_live = " live clump", "")) else name); if (!nzchar(nm)) stop("Give the dataset a name.", call. = FALSE)
  if (where == "local") {   # check before asking the server, so nothing is downloaded that cannot be clumped
    if (is.na(plink_bin())) stop("PLINK was not found in this installation, so the variants cannot be clumped on this computer. Clump on the OpenGWAS server, or run diagnostics in Settings.", call. = FALSE)
    if (!isTRUE(ld_panel_status()$installed[[pop]])) stop(sprintf("The %s (%s) LD reference panel is not installed, so the variants cannot be clumped on this computer. Download it in Settings > Data and access, or clump on the OpenGWAS server.", LD_POP_NAMES[[pop]], pop), call. = FALSE)
  }
  params <- if (where == "local") list(pval = p, clump = 0, r2 = r2, kb = kb, pop = pop, where = "local") else if (where == "server_live") list(pval = p, clump = 1, r2 = r2, kb = kb, pop = pop, where = "server_live") else list(pval = p, clump = 1, r2 = r2, kb = kb, pop = pop)
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) {
    old <- tryCatch(jsonlite::fromJSON(file.path(dest, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(old) && identical(old$source$id, id) && identical(old$source$role, "exposure") && .ogw_same(old$source$params, params)) { old$reused <- TRUE; return(old) }
    # Round 109: the same window searched again with another cut point replaces the download it just made (same dataset, same name)
    if (!is.null(old) && identical(old$source$id, id) && identical(old$source$role, "exposure") && identical(replace, nm)) unlink(dest, recursive = TRUE) else
    stop("A dataset named \"", nm, "\" is already in the Data library, with other settings or another source. Choose another name.", call. = FALSE)
  }
  .say(sprintf("Asking OpenGWAS for the variants of %s (%s) with P < %s%s", meta$trait %||% id, id, format(p), if (where == "local") ", unclumped, to clump on this computer" else if (where == "server_live") sprintf(", clumped live on the server (r\u00b2 < %s within %s kb, %s)", r2, format(kb, big.mark = ","), pop) else ", using the stored clumps"), "These become the instruments: genome-wide significant variants, one per region, so that they are independent.")
  d <- if (where == "local") .ogw_do(backend$tophits(id, p, r2, kb, pop, jwt, clump = 0)) else if (where == "server_live") .ogw_do(backend$tophits(id, p, r2, kb, pop, jwt, force = TRUE)) else .ogw_do(backend$tophits(id, p, r2, kb, pop, jwt))
  if (!is.data.frame(d) || !nrow(d)) stop("No variants reach P < ", format(p, scientific = TRUE), " in this dataset", if (where != "local") " after clumping", ". A less strict threshold (for example 5e-6) would find some, but they are weaker instruments.", call. = FALSE)
  std <- data.frame(snp = as.character(d$rsid), chr = as.character(d$chr), pos = suppressWarnings(as.numeric(d$position)), ea = d$ea, oa = d$nea, eaf = suppressWarnings(as.numeric(d$eaf)),
                    beta = suppressWarnings(as.numeric(d$beta)), se = suppressWarnings(as.numeric(d$se)), p = suppressWarnings(as.numeric(d$p)), n = suppressWarnings(as.numeric(d$n)), stringsAsFactors = FALSE)
  std$n[!is.finite(std$n)] <- meta$sample_size
  cl <- .ogw_clean_rows(std, "exposure"); std <- cl$d; std <- std[order(std$p), , drop = FALSE]; rownames(std) <- NULL
  w <- cl$warnings
  cand <- NULL; sel <- NULL; crec <- NULL
  if (where == "local") {
    if (!nrow(std)) stop("None of the variants below the threshold could be used (see the reasons above), so there is nothing to clump.", call. = FALSE)
    cand <- std; sel <- ld_clump_here(std$snp, std$p, r2, kb, pop)   # Round 28: as ieugwasr and the server (PLINK keeps the tied variant that comes first in the panel)
    # Tie sensitivity analysis (Round 28): when candidates share a P value, they are clumped a second time with the larger |beta/SE| first.
    alt <- if (anyDuplicated(std$p)) ld_clump_here(std$snp, std$p, r2, kb, pop, z = std$beta / std$se, tie_order = "strongest") else NULL
    std <- std[std$snp %in% sel$variant_id[sel$status == "index"], , drop = FALSE]; rownames(std) <- NULL
    crec <- c(attr(sel, "settings"), list(candidates = nrow(cand), clumped = sum(sel$status == "clumped"), not_in_panel = sum(sel$status == "not_in_panel"), no_rsid = sum(sel$status == "no_rsid"), kept = nrow(std)))
    tsum <- ld_tie_summary(sel, alt); if (!is.null(tsum)) crec$ties <- tsum
    if (crec$not_in_panel > 0) w <- c(w, sprintf("%d of %d candidate %s not in the reference panel, so %s LD could not be judged and %s removed (as the OpenGWAS server does).", crec$not_in_panel, crec$candidates, if (crec$not_in_panel == 1) "variant is" else "variants are", if (crec$not_in_panel == 1) "its" else "their", if (crec$not_in_panel == 1) "it was" else "they were"))
  }
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
                            packages = .ogw_versions(), note = if (where == "local") "Every variant below the threshold was downloaded without clumping (ieugwasr::tophits, clump = 0) and clumped on this computer with PLINK against the 1000 Genomes reference panel of the chosen population. The Data library table candidates (candidates.parquet) holds the downloaded variants and instrument_selection (instrument_selection.parquet) the status of each." else if (where == "server_live") "Instruments clumped live on the OpenGWAS server (ieugwasr::tophits with force_server = TRUE) against the 1000 Genomes reference panel of the chosen population." else if (ogw_stored_clumps(params)) "Instruments taken from the stored, pre-computed clumps of OpenGWAS (ieugwasr::tophits with the default settings uses preclumped = 1), made against the 1000 Genomes reference panel of the chosen population." else "Instruments clumped on the OpenGWAS server (ieugwasr::tophits) against the 1000 Genomes reference panel of the chosen population."))
  if (!is.null(crec)) rec$clump <- crec
  rec <- .ogw_write(nm, std, rec, root)
  if (!is.null(sel)) {
    dest <- file.path(library_dir(root), nm)
    lib_table_write(cand, dest, "candidates"); lib_table_write(sel, dest, "instrument_selection")
    if (!is.null(alt)) lib_table_write(alt, dest, "instrument_selection_ties")
  }
  rec
}

#' The same variants looked up in an outcome. Variants the outcome lacks are replaced by proxies (r² at least `rsq`) when asked.
ogw_get_outcome <- function(id, exposure_name, jwt, proxies = TRUE, rsq = 0.8, maf = 0.3, name = NULL, offline = FALSE, root = studies_root(), backend = ogw_backend()) {
  .ogw_need(jwt, offline)
  if (!is.finite(rsq) || rsq < 0.5 || rsq > 1) stop("The proxy r² must be between 0.5 and 1.", call. = FALSE)
  if (!is.finite(maf) || maf <= 0 || maf > 0.5) stop("The allele frequency limit for palindromic proxies must be above 0 and at most 0.5.", call. = FALSE)
  ex <- library_dataset(exposure_name, root)
  if (!isTRUE(ex$recipe$usable_as_exposure)) stop("\"", exposure_name, "\" is not a set of independent instruments, so there is nothing to look up.", call. = FALSE)
  snps <- unique(ex$data$snp[!is.na(ex$data$snp) & grepl("^rs[0-9]+$", ex$data$snp)])
  # Round 28: the lead variants that only the tie sensitivity analysis keeps are looked up too (they are not counted as instruments).
  extra <- setdiff(unlist(ex$recipe$clump$ties[["in"]]), snps)
  if (length(snps) < 2L) stop("The exposure \"", exposure_name, "\" has no rsID variants to look up. Looking up by position is not built.", call. = FALSE)
  meta <- ogw_meta(id, jwt, backend)
  # Short exposure tag (Round 25): the dataset ID, plus " live clump" or " clumped here" for those copies, or the first characters of
  # a name the user gave, so two copies of one exposure get different outcome names; names are cut at 60 characters.
  sid <- ex$recipe$source$id; suf <- NULL
  if (!is.null(sid)) for (lim in c(60L, 48L)) { base <- .ogw_label(ex$recipe$source$trait %||% "", sid, lim); x <- if (startsWith(exposure_name, base)) substring(exposure_name, nchar(base) + 1L) else NULL; if (!is.null(x) && x %in% c("", " live clump", " clumped here")) { suf <- x; break } }
  ex_short <- clean_study_name(if (!is.null(suf) && suf %in% c("", " live clump", " clumped here")) paste0(substr(sid, 1L, 18L), suf) else substr(exposure_name, 1L, 22L))   # suf: NULL (user's own name) or a known suffix
  nm <- clean_study_name(if (is.null(name) || !nzchar(name)) paste0(.ogw_label(meta$trait, id, if (identical(suf, "")) 38L else max(24L, 55L - nchar(ex_short))), " for ", ex_short) else name); if (!nzchar(nm)) stop("Give the dataset a name.", call. = FALSE)
  params <- list(proxies = isTRUE(proxies), rsq = rsq, maf_threshold = maf, exposure = exposure_name, exposure_sha256 = ex$recipe$original$sha256 %||% "")
  dest <- file.path(library_dir(root), nm)
  if (dir.exists(dest)) {
    old <- tryCatch(jsonlite::fromJSON(file.path(dest, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(old) && identical(old$source$id, id) && identical(old$source$role, "outcome") && .ogw_same(old$source$params, params)) { old$reused <- TRUE; return(old) }
    stop("A dataset named \"", nm, "\" is already in the Data library, with other settings or another source. Choose another name.", call. = FALSE)
  }
  .say(sprintf("Looking up the %d instruments of %s in %s (%s)%s", length(snps), exposure_name, meta$trait %||% id, id, if (isTRUE(proxies)) sprintf(", with proxies (r\u00b2 of at least %s) where one is missing", rsq) else ""), "OpenGWAS returns each instrument's effect on the outcome; a proxy is a nearby variant in strong linkage, aligned to the instrument's alleles.")
  o <- .ogw_do(backend$outcome(c(snps, extra), id, isTRUE(proxies), rsq, maf, jwt))
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
  found <- sum(std$snp %in% snps); n_prox <- sum(std$proxy_snp != "" & std$snp %in% snps); missing <- setdiff(snps, std$snp)
  if (found < 2L) stop("Only ", found, " of ", length(snps), " instruments were found in this outcome. At least 2 are needed.", call. = FALSE)
  if (length(missing)) w <- c(w, sprintf("%d of %d instruments are not in this outcome%s and are lost: %s.", length(missing), length(snps), if (isTRUE(proxies)) " and had no proxy" else " (proxies were off)", paste(utils::head(missing, 8), collapse = ", ")))
  # Round 109: why each lost instrument had no proxy, from the 1000 Genomes panel of the outcome's ancestry when it is installed
  why <- if (length(missing) && isTRUE(proxies)) tryCatch(proxy_reasons(missing, rsq, maf, ld_pop_for(OGW_POP_LABELS[[meta$population %||% ""]] %||% "EUR")), error = function(e) NULL) else NULL
  if (length(missing)) w <- c(w, if (!isTRUE(proxies)) "Proxies were off, so no proxy was searched." else if (is.null(why)) "Why no proxy was found cannot be told without the 1000 Genomes reference panel of this ancestry (Settings > Data and access)." else paste0("Why no proxy: ", paste(sprintf("%s: %s", why$variant_id, why$why), collapse = "; "), "."))
  if (n_prox) w <- c(w, sprintf("%d instruments were replaced by a proxy variant (r² at least %s). The proxy's own ID is in the column proxy_snp.", n_prox, format(rsq)))
  if (!is.finite(meta$sample_size)) w <- c(w, "The catalogue does not give a sample size for this dataset.")
  rec <- list(trait = nm, role = "outcome",
              facts = list(name = nm, role = "outcome", kind = ogw_kind(meta$ncase), scale = "beta", ea_ok = "yes", frq = "eaf", build = .ogw_build(meta$build), ancestry = .ogw_ancestry(meta$population), n = .ogw_n(meta$sample_size),
                           scope = "selected", independent = "no", ukb = .ogw_is_ukb(meta)),
              mapping = list(snp = "rsid", chr = "chr", pos = "position", ea = "ea", oa = "nea", eaf = "eaf", beta = "beta", se = "se", p = "p", n = "n"),
              original = list(file = paste0("OpenGWAS:", id)), rows_in_file = nrow(o), rows_kept = nrow(std), warnings = as.list(w), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              usable_as_exposure = FALSE, usable_as_outcome = TRUE, usable_in_runs = TRUE,
              lookup = list(requested = length(snps), found = found, direct = found - n_prox, by_proxy = n_prox, missing = length(missing), missing_why = if (!is.null(why)) stats::setNames(as.list(why$why), why$variant_id) else NULL, tie_alternatives = length(extra), tie_alternatives_found = sum(extra %in% std$snp)),
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
