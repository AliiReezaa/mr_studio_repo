# Complete GWAS inside the app (Round 75): the genome-wide summary statistics of an OpenGWAS dataset, downloaded as its GWAS-VCF file,
# and GWAS-VCF files the user already has. LD score regression, MRlap, MR-APSS, CAUSE and LHC-MR need every variant of both traits;
# an OpenGWAS instrument download holds only the instruments and their lookups.
#
# Before a download starts, the file's size is asked from the server and compared with the free space on the disk that holds the
# Data library: the download itself, plus the standard copy MR Studio writes from it, plus a margin. Too little space stops it before
# anything is written; little space left afterwards gives a warning. The download runs in a separate R process (the app stays
# usable), resumes where it stopped (also after the app was closed, or when the link has expired: a new link is asked for), and is
# then read into the Data library like an imported genome-wide file. The token never leaves this process: the background job
# receives only the download link.
#
# GWAS-VCF (Lyon et al. 2021, Genome Biology 22:32): one variant per row, the effect of the ALT allele in the sample column, whose
# fields are named by FORMAT: ES effect, SE standard error, LP -log10 P, AF ALT allele frequency, SS sample size, NC cases, ID rsID.

`%||%` <- function(a, b) if (is.null(a)) b else a

FULLGWAS_MARGIN <- 2 * 1024^3          # space that should stay free after the download and the standard copy
FULLGWAS_STD_PER_ROW <- 60             # bytes per variant of the standard copy (Parquet, zstd): 53 measured on 1 million variants
FULLGWAS_STD_RATIO <- 2.2              # or 2.2 times the .vcf.gz (1.9 measured), whichever is larger
FULLGWAS_MEM <- function(variants) 0.8 + 0.6 * max(variants, 0, na.rm = TRUE) / 1e6    # GB at the peak of reading (0.56 per million measured)

#' Free bytes on the disk that holds `path` (NA when the system does not tell). macOS and Linux: df; Windows: .NET DriveInfo.
disk_free_bytes <- function(path) {
  p <- path; while (!dir.exists(p) && dirname(p) != p) p <- dirname(p)
  tryCatch({
    if (.Platform$OS.type == "windows") {
      drv <- substr(normalizePath(p, winslash = "\\"), 1, 2)
      if (!grepl("^[A-Za-z]:$", drv)) return(NA_real_)   # a network path: not measured
      o <- system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf("[System.IO.DriveInfo]::new('%s').AvailableFreeSpace", drv))), stdout = TRUE, stderr = FALSE)
      v <- suppressWarnings(as.numeric(trimws(o[length(o)])))
    } else {
      o <- system2("df", c("-Pk", shQuote(p)), stdout = TRUE, stderr = FALSE)
      v <- suppressWarnings(as.numeric(strsplit(trimws(o[length(o)]), "\\s+")[[1]][4])) * 1024
    }
    if (length(v) == 1L && is.finite(v)) v else NA_real_
  }, error = function(e) NA_real_)
}

fmt_gb <- function(b) if (!is.finite(b)) "unknown" else if (b >= 1024^3) sprintf("%.1f GB", b / 1024^3) else sprintf("%.0f MB", max(b, 0) / 1024^2)

#' Size of the file behind a link, in bytes (NA when the server does not say). Asked without downloading: a HEAD request, and a
#' one-byte range request when the server refuses HEAD (signed download links often allow GET only).
http_size <- function(url) {
  hdr <- function(r) { if (is.null(r) || r$status_code >= 400) return(NA_real_); h <- curl::parse_headers_list(r$headers)
    cr <- h[["content-range"]]; if (!is.null(cr) && grepl("/[0-9]+$", cr)) return(as.numeric(sub(".*/", "", cr)))
    cl <- h[["content-length"]]; if (!is.null(cl) && r$status_code < 300 && is.null(h[["content-encoding"]])) return(as.numeric(cl)); NA_real_ }
  get <- function(...) tryCatch(curl::curl_fetch_memory(url, curl::new_handle(followlocation = TRUE, connecttimeout = 30, timeout = 60, ...)), error = function(e) NULL)
  s <- hdr(get(nobody = TRUE)); if (is.finite(s) && s > 0) return(s)
  hdr(get(range = "0-0"))
}

#' The download link of an OpenGWAS dataset's complete summary statistics (its GWAS-VCF), from ieugwasr::gwasinfo_files. The link
#' is signed and expires after a while, so it is asked for again for every start or resume, and never stored.
.ogw_real_files <- function(id, jwt) ieugwasr::gwasinfo_files(id, opengwas_jwt = jwt)
ogw_vcf_link <- function(id, jwt, backend = ogw_backend()) {
  f <- backend$files %||% .ogw_real_files
  u <- unique(grep("^https?://", unlist(.ogw_do(f(id, jwt)), use.names = FALSE), value = TRUE))
  v <- u[grepl("\\.vcf\\.(b)?gz($|[?#])", u)]
  if (!length(v)) stop("OpenGWAS lists no complete summary statistics file (GWAS-VCF) for ", id, ". Some datasets are available only as top hits; use it as an exposure or outcome download instead.", call. = FALSE)
  v[1]
}

.fg_lib <- function(root) library_dir(root)
.fg_dl_dir <- function(key, root) file.path(.fg_lib(root), ".downloads", gsub("[^A-Za-z0-9._-]", "_", key))
.fg_status_file <- function(dir) file.path(dir, "status.json")
.lib_recipe <- function(name, root = studies_root()) { f <- file.path(library_dir(root), basename(name), "recipe.json")
  if (file.exists(f)) tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL) else NULL }

#' Name of the complete GWAS of OpenGWAS dataset `id` in the Data library (NULL when there is none).
fullgwas_have <- function(id, root = studies_root()) {
  for (r in tryCatch(library_list(root), error = function(e) list())) if (identical(r$source$type, "opengwas_vcf") && identical(r$source$id, id)) return(r$trait)
  NULL
}

#' The datasets the complete-GWAS methods use for a pair: an OpenGWAS instrument download (exposure or outcome) is replaced by the
#' complete GWAS of the same OpenGWAS dataset when that is in the Data library, so a pair made from OpenGWAS downloads can run them.
full_gwas_ids <- function(ids, root = studies_root()) unname(vapply(ids, function(i) {
  if (!startsWith(i, "lib:")) return(i)
  r <- .lib_recipe(sub("^lib:", "", i), root)
  if (!is.null(r) && identical(r$source$type, "opengwas") && !is.null(r$source$id)) { c <- fullgwas_have(r$source$id, root); if (!is.null(c)) return(paste0("lib:", c)) }
  i }, ""))

#' Before a download: the file, its size, the space it needs and the space that is free, and whether it may start.
#' level: "ok", "warn" (it fits, but little space is left), "block" (it does not fit) or "unknown" (free space could not be measured).
fullgwas_plan <- function(id, jwt, offline = FALSE, root = studies_root(), backend = ogw_backend(), size_fun = http_size, free_fun = disk_free_bytes) {
  .ogw_need(jwt, offline)
  meta <- ogw_meta(id, jwt, backend); url <- ogw_vcf_link(id, jwt, backend)
  bytes <- size_fun(url); nsnp <- suppressWarnings(as.numeric(meta$nsnp %||% NA))
  part <- file.path(.fg_dl_dir(id, root), "gwas.vcf.gz.part"); have <- if (file.exists(part)) file.size(part) else 0
  std <- max(FULLGWAS_STD_RATIO * bytes, FULLGWAS_STD_PER_ROW * nsnp, na.rm = TRUE); if (!is.finite(std)) std <- NA_real_
  need <- max(bytes - have, 0) + std
  dir.create(.fg_lib(root), recursive = TRUE, showWarnings = FALSE); free <- free_fun(.fg_lib(root))
  level <- if (!is.finite(bytes) || !is.finite(free)) "unknown" else if (free < need) "block" else if (free - need < FULLGWAS_MARGIN) "warn" else "ok"
  msg <- switch(level,
    block = sprintf("Not enough free space: this needs about %s (the %s download and the %s standard copy MR Studio makes from it), and %s is free on the disk that holds the Data library. Free some space, or choose another folder for MR Studio's studies and Data library in Settings, then try again.", fmt_gb(need), fmt_gb(bytes - have), fmt_gb(std), fmt_gb(free)),
    warn = sprintf("It fits, but only about %s would stay free afterwards (%s needed, %s free). The computer can slow down when its disk is nearly full. Consider freeing space first, or not keeping the original VCF file.", fmt_gb(free - need), fmt_gb(need), fmt_gb(free)),
    unknown = if (!is.finite(bytes)) "OpenGWAS did not say how large the file is, so the free space cannot be checked in advance. Complete GWAS files are usually 0.2 to 1.5 GB; make sure several GB are free." else "The free space on this disk could not be measured (for example on a network drive). Make sure enough space is free before you start.",
    ok = sprintf("%s needed, %s free on this disk.", fmt_gb(need), fmt_gb(free)))
  mem <- FULLGWAS_MEM(nsnp); mi <- mem_info(); low <- is.finite(mi$available_gb) && mem > RES_SAFE_SHARE * mi$available_gb
  mmsg <- if (low) sprintf("Reading the file needs about %.1f GB of memory, and %.1f GB is free now. Close other programs before the download finishes, so that reading it does not slow the computer down.", mem, mi$available_gb) else ""
  list(id = id, trait = meta$trait, population = meta$population, nsnp = if (is.finite(nsnp)) nsnp else NULL, ncase = meta$ncase, ncontrol = meta$ncontrol, sample_size = meta$sample_size,
       bytes = if (is.finite(bytes)) bytes else NULL, have_bytes = have, std_bytes = if (is.finite(std)) std else NULL, need_bytes = if (is.finite(need)) need else NULL,
       free_bytes = if (is.finite(free)) free else NULL, level = level, msg = msg, mem_gb = mem, mem_msg = mmsg,
       folder = .fg_lib(root), existing = fullgwas_have(id, root), name = .ogw_label(paste(meta$trait, "complete GWAS"), id),
       note = "Downloads of complete GWAS count against your OpenGWAS allowance. The download continues in the background; you can keep working.")
}

# ---- the background job ---------------------------------------------------------------------------------------------------------
.fg_write_status <- function(dir, ...) { f <- .fg_status_file(dir); s <- if (file.exists(f)) tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) list()) else list()
  x <- list(...); for (k in names(x)) s[[k]] <- x[[k]]; s$time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  tmp <- paste0(f, ".tmp"); writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null"), tmp); file.rename(tmp, f); invisible(s) }

#' Start (or resume) the background job. `spec`: key, mode ("opengwas" with url and bytes, or "local" with path), name, facts and
#' source (written into the recipe), keep_vcf.
fullgwas_job_start <- function(spec, root = studies_root(), engine = getOption("mrgui.engine_dir", "")) {
  dir <- .fg_dl_dir(spec$key, root); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  j <- fullgwas_job_read(spec$key, root); if (!is.null(j) && identical(j$state %||% "", "running") && isTRUE(j$alive)) stop("This download is already running.", call. = FALSE)
  nm <- clean_study_name(spec$name %||% ""); if (!nzchar(nm)) stop("Give the dataset a name.", call. = FALSE)
  if (dir.exists(file.path(.fg_lib(root), nm))) stop("A dataset named \"", nm, "\" already exists in the Data library. Choose another name.", call. = FALSE)
  spec$name <- nm; spec$library <- normalizePath(.fg_lib(root), winslash = "/", mustWork = FALSE)
  writeLines(jsonlite::toJSON(spec, auto_unbox = TRUE, null = "null", digits = NA), file.path(dir, "spec.json"))
  unlink(file.path(dir, "cancel"))
  .fg_write_status(dir, key = spec$key, name = nm, mode = spec$mode, state = "running", step = if (identical(spec$mode, "local")) "reading" else "downloading",
                   msg = "", bytes = spec$bytes %||% NULL, rows = 0, trait = spec$source$trait %||% nm, id = spec$source$id %||% NULL)
  code <- sprintf("source(file.path(%s, 'fullgwas_core.R')); fullgwas_job_run(%s)", deparse(engine), deparse(normalizePath(dir, winslash = "/")))
  rscript_background(c("--vanilla", "-e", shQuote(code)), log = file.path(dir, "job.log"))
  invisible(fullgwas_job_read(spec$key, root))
}

#' Download the rest of a file and add it to the part already downloaded. Streamed, so what arrives before a dropped connection is
#' kept. (curl::curl_download with mode "ab" writes to a new temporary file and renames it over the part, losing what was there.)
fetch_append <- function(url, part, handle) {
  con <- curl::curl(url, "rb", handle = handle); on.exit(close(con))
  out <- file(part, "ab"); on.exit(close(out), add = TRUE)
  repeat { b <- readBin(con, "raw", 1048576L); if (!length(b)) break; writeBin(b, out) }
  invisible(TRUE)
}

#' The job itself (runs in the background R process).
fullgwas_job_run <- function(dir) {
  spec <- jsonlite::fromJSON(file.path(dir, "spec.json"), simplifyVector = FALSE)
  st <- function(...) .fg_write_status(dir, ...)
  st(pid = Sys.getpid())
  tryCatch({
    vcf <- if (identical(spec$mode, "local")) spec$path else file.path(dir, "gwas.vcf.gz.part")
    if (!identical(spec$mode, "local")) {
      bytes <- as.numeric(spec$bytes %||% NA)
      for (attempt in 1:30) {
        have <- if (file.exists(vcf)) file.size(vcf) else 0
        if (is.finite(bytes) && have >= bytes) break
        h <- curl::new_handle(resume_from = have, low_speed_limit = 100, low_speed_time = 120, connecttimeout = 60, followlocation = TRUE, failonerror = TRUE)
        err <- tryCatch({ fetch_append(spec$url, vcf, h); NULL }, error = function(e) conditionMessage(e))
        if (is.null(err) && !is.finite(bytes)) break
        if (!is.null(err) && grepl("40[13]|expired|Forbidden|Unauthorized", err)) stop("The download link has expired. Choose Resume: MR Studio asks OpenGWAS for a new link and continues where it stopped.")
        if (!is.null(err) && grepl("range|resume", err, ignore.case = TRUE)) { unlink(vcf); next }   # the server cannot continue a part: start the file again
        if (!is.null(err)) Sys.sleep(10)
      }
      got <- file.size(vcf)
      if (is.finite(bytes) && got != bytes) { if (got > bytes) unlink(vcf); stop(sprintf("The download stopped at %s of %s. Choose Resume to continue where it stopped.", fmt_gb(got), fmt_gb(bytes))) }
    }
    st(step = "reading", msg = "")
    sha <- digest::digest(file = vcf, algo = "sha256")
    res <- vcf_read_gwas(vcf, progress = function(n) { if (file.exists(file.path(dir, "cancel"))) stop("cancelled"); st(rows = n) })
    st(step = "saving", rows = res$stats$rows_read)
    rec <- fullgwas_save(res, spec, sha, vcf_bytes = file.size(vcf))
    if (!identical(spec$mode, "local")) {
      dest <- file.path(spec$library, spec$name)
      if (isTRUE(spec$keep_vcf)) file.rename(vcf, file.path(dest, paste0(spec$source$id %||% "gwas", ".vcf.gz"))) else unlink(vcf)
    }
    st(state = "done", step = "done", rows = rec$rows_kept, msg = sprintf("%s variants saved as \"%s\".", format(rec$rows_kept, big.mark = ","), spec$name))
  }, error = function(e) {
    m <- conditionMessage(e)
    if (identical(m, "cancelled") || file.exists(file.path(dir, "cancel"))) st(state = "cancelled", msg = "Stopped. What was downloaded is kept, so Resume continues where it stopped.")
    else st(state = "error", msg = m)
  })
}

#' Status of one job: state running / done / error / cancelled, with the bytes downloaded so far. A job whose process has gone
#' (the app or the computer was closed) is reported as stopped, ready to resume.
fullgwas_job_read <- function(key, root = studies_root()) {
  dir <- .fg_dl_dir(key, root); f <- .fg_status_file(dir); if (!file.exists(f)) return(NULL)
  s <- tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL); if (is.null(s)) return(NULL)
  part <- file.path(dir, "gwas.vcf.gz.part"); s$done_bytes <- if (file.exists(part)) file.size(part) else if (identical(s$state, "done")) s$bytes else 0
  s$alive <- identical(s$state, "running") && .fg_alive(s$pid)
  if (identical(s$state, "running") && !is.null(s$pid) && !s$alive) { s$state <- "cancelled"; s$msg <- "Stopped when MR Studio or the computer was closed. Resume continues where it stopped." }
  s$key <- key; s
}
.fg_alive <- function(pid) { if (is.null(pid)) return(TRUE)   # not started yet
  pid <- as.integer(pid)
  if (.Platform$OS.type == "windows") any(grepl(paste0(" ", pid, " "), suppressWarnings(system2("tasklist", c("/FI", shQuote(paste("PID eq", pid))), stdout = TRUE, stderr = FALSE))))
  else identical(suppressWarnings(system2("kill", c("-0", pid), stdout = FALSE, stderr = FALSE)), 0L) }
fullgwas_jobs <- function(root = studies_root()) {
  d <- file.path(.fg_lib(root), ".downloads"); if (!dir.exists(d)) return(list())
  Filter(Negate(is.null), lapply(basename(list.dirs(d, recursive = FALSE)), function(k) fullgwas_job_read(k, root)))
}
fullgwas_cancel <- function(key, root = studies_root()) {
  dir <- .fg_dl_dir(key, root); writeLines("cancel", file.path(dir, "cancel")); s <- fullgwas_job_read(key, root)
  if (!is.null(s$pid) && isTRUE(s$alive)) tools::pskill(as.integer(s$pid))
  .fg_write_status(dir, state = "cancelled", msg = "Stopped. What was downloaded is kept, so Resume continues where it stopped."); invisible(TRUE)
}
#' Remove a job and its partial download (never a finished dataset in the Data library).
fullgwas_discard <- function(key, root = studies_root()) {
  s <- fullgwas_job_read(key, root); if (!is.null(s) && isTRUE(s$alive)) fullgwas_cancel(key, root)
  unlink(.fg_dl_dir(key, root), recursive = TRUE); invisible(TRUE)
}
#' Resume a stopped OpenGWAS download: a new link is asked for (the old one may have expired), the rest of the file is downloaded.
fullgwas_resume <- function(key, jwt, offline = FALSE, root = studies_root(), backend = ogw_backend()) {
  dir <- .fg_dl_dir(key, root); spec <- jsonlite::fromJSON(file.path(dir, "spec.json"), simplifyVector = FALSE)
  if (identical(spec$mode, "opengwas")) { .ogw_need(jwt, offline); spec$url <- ogw_vcf_link(spec$source$id, jwt, backend) }
  spec$library <- NULL; fullgwas_job_start(spec, root)
}

#' Start an OpenGWAS complete-GWAS download after fullgwas_plan. Refused when the plan says the disk is too full.
fullgwas_download_start <- function(id, name, jwt, keep_vcf = TRUE, offline = FALSE, root = studies_root(), backend = ogw_backend(), plan = NULL) {
  p <- plan %||% fullgwas_plan(id, jwt, offline, root, backend)
  if (identical(p$level, "block")) stop(p$msg, call. = FALSE)
  if (!is.null(p$existing)) stop("The complete GWAS of ", id, " is already in the Data library as \"", p$existing, "\".", call. = FALSE)
  meta <- ogw_meta(id, jwt, backend); url <- ogw_vcf_link(id, jwt, backend)
  kind <- ogw_kind(suppressWarnings(as.numeric(meta$ncase %||% NA)))
  anc <- ogw_pop_code(meta$population); facts <- list(name = name, role = "either", kind = kind, scale = "beta", ea_ok = "yes", frq = "eaf", strand = "forward", build = .ogw_build(meta$build),
    ancestry = if (is.na(anc)) "" else anc, n = .ogw_n(meta$sample_size), ncase = .ogw_n(meta$ncase), ncontrol = .ogw_n(meta$ncontrol), scope = "genomewide", independent = "clump",
    clump_r2 = 0.001, clump_kb = 10000, clump_pop = if (is.na(anc)) "EUR" else anc, ukb = .ogw_is_ukb(meta))
  src <- list(type = "opengwas_vcf", id = id, trait = meta$trait, population = meta$population, sample_size = meta$sample_size, ncase = meta$ncase, ncontrol = meta$ncontrol,
              author = meta$author, year = meta$year, consortium = meta$consortium, pmid = meta$pmid, build = meta$build,
              file = paste0(id, ".vcf.gz"), link = sub("[?#].*$", "", url), downloaded = format(Sys.Date()), versions = .ogw_versions())
  fullgwas_job_start(list(key = id, mode = "opengwas", url = url, bytes = p$bytes, name = name, keep_vcf = isTRUE(keep_vcf), facts = facts, source = src), root)
}

# ---- GWAS-VCF -------------------------------------------------------------------------------------------------------------------
.vcf_sample_meta <- function(hdr) {
  s <- grep("^##SAMPLE=<", hdr, value = TRUE); out <- list()
  if (length(s)) { kv <- regmatches(s[1], gregexpr("[A-Za-z]+=(\"[^\"]*\"|[^,>]*)", s[1]))[[1]]
    for (x in kv) { k <- sub("=.*", "", x); out[[k]] <- gsub('^"|"$', "", sub("^[^=]*=", "", x)) } }
  a <- regmatches(hdr, regexpr("assembly=[^,>]+", hdr)); if (length(a)) out$assembly <- sub("assembly=", "", a[1])
  out
}

#' Header and first rows of a VCF, to check that it is GWAS-VCF and to fill the import form. Reads only the start of the file.
vcf_peek <- function(path, rows = 5L) {
  con <- gzfile(path, "rt"); on.exit(close(con)); hdr <- character(); cols <- NULL
  repeat { l <- readLines(con, n = 200L, warn = FALSE); if (!length(l)) break
    i <- which(!startsWith(l, "##")); hdr <- c(hdr, l[startsWith(l, "##")]); if (length(i)) { cols <- l[i[1]]; first <- l[-seq_len(i[1])]; break } }
  if (is.null(cols) || !startsWith(cols, "#CHROM")) stop("This is not a VCF file: the #CHROM header line was not found.", call. = FALSE)
  cols <- strsplit(sub("^#", "", cols), "\t", fixed = TRUE)[[1]]
  if (length(cols) < 10L) stop("This VCF has no sample column, so it holds no summary statistics (GWAS-VCF keeps them in the sample column).", call. = FALSE)
  first <- c(first, readLines(con, n = max(0L, rows - length(first)), warn = FALSE))[seq_len(rows)]; first <- first[!is.na(first)]
  fm <- if (length(first)) strsplit(first[1], "\t", fixed = TRUE)[[1]][9] else ""
  keys <- strsplit(fm %||% "", ":", fixed = TRUE)[[1]]
  if (!all(c("ES", "SE") %in% keys)) stop("This VCF has no ES and SE fields, so it is not a GWAS summary statistics file (GWAS-VCF). A genotype VCF cannot be used.", call. = FALSE)
  m <- .vcf_sample_meta(hdr)
  list(columns = cols, samples = cols[-(1:9)], format = keys, meta = m, rows = first, header_lines = length(hdr),
       kind = if (grepl("case", m$StudyType %||% "", ignore.case = TRUE)) "binary" else if (nzchar(m$StudyType %||% "")) "continuous" else "",
       ncase = m$TotalCases %||% "", ncontrol = m$TotalControls %||% "", variants = m$TotalVariants %||% "", build = m$assembly %||% "")
}

#' Read a GWAS-VCF (plain or gzipped) into MR Studio's standard columns, in chunks. The effect belongs to the ALT allele, so ea = ALT,
#' oa = REF and eaf = AF. P = 10^-LP. Only single-nucleotide variants with one ALT allele are kept (as for imported files); every
#' row left out is counted with its reason. With several samples (studies) in one file, the first is used.
vcf_read_gwas <- function(path, chunk = 250000L, progress = NULL) {
  pk <- vcf_peek(path, rows = 1L)
  con <- gzfile(path, "rt"); on.exit(close(con))
  repeat { l <- readLines(con, n = 1L, warn = FALSE); if (!length(l) || startsWith(l, "#CHROM")) break }
  num <- function(x) { x[x %in% c(".", "", "NA")] <- NA; suppressWarnings(as.numeric(x)) }
  parts <- list(); n <- 0; drop <- c(not_snv = 0, multiallelic = 0, missing = 0, duplicate = 0)
  repeat {
    x <- readLines(con, n = chunk, warn = FALSE); if (!length(x)) break
    x <- x[nzchar(x)]; n <- n + length(x)
    t <- data.table::fread(text = x, sep = "\t", header = FALSE, colClasses = "character", select = c(1:5, 9, 10), quote = "", na.strings = character(0), showProgress = FALSE)
    data.table::setnames(t, c("chr", "pos", "id", "ref", "alt", "fmt", "smp"))
    v <- stats::setNames(rep(list(rep(NA_character_, nrow(t))), 6L), c("ES", "SE", "LP", "AF", "SS", "ID"))
    for (f in unique(t$fmt)) { i <- which(t$fmt == f); keys <- strsplit(f, ":", fixed = TRUE)[[1]]
      sp <- data.table::tstrsplit(t$smp[i], ":", fixed = TRUE, fill = NA_character_)
      for (k in intersect(names(v), keys)) { j <- match(k, keys); if (j <= length(sp)) v[[k]][i] <- sp[[j]] } }
    multi <- grepl(",", t$alt, fixed = TRUE); snv <- grepl("^[ACGTacgt]$", t$ref) & grepl("^[ACGTacgt]$", t$alt)
    d <- data.frame(snp = ifelse(grepl("^rs[0-9]+$", t$id), t$id, ifelse(grepl("^rs[0-9]+$", v$ID %||% ""), v$ID, paste0(t$chr, ":", t$pos))),
                    chr = sub("^chr", "", t$chr), pos = num(t$pos), ea = toupper(t$alt), oa = toupper(t$ref), beta = num(v$ES), se = num(v$SE),
                    p = 10^(-num(v$LP)), eaf = num(v$AF), n = num(v$SS), stringsAsFactors = FALSE)
    ok <- is.finite(d$beta) & is.finite(d$se) & d$se > 0
    drop["multiallelic"] <- drop["multiallelic"] + sum(multi); drop["not_snv"] <- drop["not_snv"] + sum(!multi & !snv)
    drop["missing"] <- drop["missing"] + sum(!multi & snv & !ok)
    parts[[length(parts) + 1L]] <- d[!multi & snv & ok, , drop = FALSE]
    if (is.function(progress)) progress(n)
  }
  d <- as.data.frame(data.table::rbindlist(parts)); rm(parts)
  if (!nrow(d)) stop("No usable variant was read from this VCF (each row needs ES and SE for one single-nucleotide ALT allele).", call. = FALSE)
  dup <- duplicated(d$snp); drop["duplicate"] <- sum(dup); d <- d[!dup, , drop = FALSE]
  d$p[!is.finite(d$p)] <- 2 * stats::pnorm(-abs(d$beta / d$se))[!is.finite(d$p)]
  for (k in c("eaf", "n")) if (all(is.na(d[[k]]))) d[[k]] <- NULL
  rownames(d) <- NULL
  list(data = d, peek = pk, stats = list(rows_read = n, rows_kept = nrow(d), dropped = as.list(drop), rsids = sum(grepl("^rs[0-9]+$", d$snp)),
       header_variants = suppressWarnings(as.numeric(pk$meta$TotalVariants %||% NA)), sample = pk$samples[1], samples = length(pk$samples)))
}

#' Save what vcf_read_gwas read as a Data library dataset (standard copy + recipe), written to a temporary folder and renamed, so an
#' unfinished dataset never appears in the library.
fullgwas_save <- function(res, spec, sha, vcf_bytes) {
  lib <- spec$library; dest <- file.path(lib, spec$name); tmp <- file.path(lib, paste0(".", spec$name, ".partial"))
  unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
  nanoparquet::write_parquet(res$data, file.path(tmp, "standardised.parquet"), compression = "zstd")
  f <- spec$facts; f$name <- spec$name; m <- res$peek$meta
  if (!nzchar(f$kind %||% "")) f$kind <- if (nzchar(res$peek$kind)) res$peek$kind else "continuous"
  if (!nzchar(as.character(f$ncase %||% "")) && nzchar(res$peek$ncase)) f$ncase <- res$peek$ncase
  if (!nzchar(as.character(f$ncontrol %||% "")) && nzchar(res$peek$ncontrol)) f$ncontrol <- res$peek$ncontrol
  if (!nzchar(as.character(f$n %||% "")) && !("n" %in% names(res$data))) { nc <- suppressWarnings(as.numeric(f$ncase)); nn <- suppressWarnings(as.numeric(f$ncontrol)); if (is.finite(nc) && is.finite(nn)) f$n <- nc + nn }
  warn <- character()
  if (res$stats$rsids < 0.5 * res$stats$rows_kept) warn <- c(warn, "Fewer than half of the variants have an rsID, so this dataset can be an outcome but cannot be clumped as an exposure.")
  if (res$stats$samples > 1) warn <- c(warn, sprintf("The file holds %d studies; the first (%s) was read.", res$stats$samples, res$stats$sample))
  hv <- res$stats$header_variants; if (is.finite(hv) && hv > 0 && abs(res$stats$rows_read - hv) > 0.001 * hv) warn <- c(warn, sprintf("The file's header lists %s variants, but %s rows were read. The file may be incomplete.", format(hv, big.mark = ","), format(res$stats$rows_read, big.mark = ",")))
  src <- spec$source %||% list(type = "vcf_file"); src$sha256 <- sha; src$bytes <- vcf_bytes; src$vcf_sample <- res$stats$sample
  rec <- list(trait = spec$name, role = "either", facts = f, mapping = list(format = "GWAS-VCF", ea = "ALT", oa = "REF", beta = "ES", se = "SE", p = "10^-LP", eaf = "AF", n = "SS", snp = "ID"),
              original = list(file = src$file %||% basename(spec$path %||% "gwas.vcf.gz"), bytes = vcf_bytes, sha256 = sha), source = src,
              rows_in_file = res$stats$rows_read, rows_kept = res$stats$rows_kept, dropped = res$stats$dropped, warnings = as.list(warn), saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              usable_as_exposure = res$stats$rsids >= 0.5 * res$stats$rows_kept, usable_as_outcome = TRUE, usable_in_runs = TRUE, complete_gwas = TRUE)
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), file.path(tmp, "recipe.json"))
  if (dir.exists(dest)) { unlink(tmp, recursive = TRUE); stop("A dataset named \"", spec$name, "\" appeared in the Data library meanwhile.") }
  if (!file.rename(tmp, dest)) stop("Could not save the dataset in ", lib)
  rec
}

#' Import a GWAS-VCF the user has (from "Import a GWAS file"). Read in the background like a download; the user's file is not changed.
vcf_import_start <- function(path, name, facts = list(), root = studies_root()) {
  pk <- vcf_peek(path); nm <- clean_study_name(name %||% "")
  need <- FULLGWAS_STD_RATIO * file.size(path); free <- disk_free_bytes(.fg_lib(root))
  if (is.finite(free) && free < need) stop(sprintf("Not enough free space: the standard copy of this file needs about %s, and %s is free on the disk that holds the Data library.", fmt_gb(need), fmt_gb(free)), call. = FALSE)
  anc <- toupper(facts$ancestry %||% ""); kind <- facts$kind %||% pk$kind
  f <- list(name = nm, role = "either", kind = kind, scale = "beta", ea_ok = "yes", frq = "eaf", strand = "forward", build = facts$build %||% pk$build, ancestry = anc,
            n = facts$n %||% "", ncase = facts$ncase %||% pk$ncase, ncontrol = facts$ncontrol %||% pk$ncontrol, scope = "genomewide", independent = "clump",
            clump_r2 = 0.001, clump_kb = 10000, clump_pop = if (anc %in% c("EUR", "EAS", "AFR", "AMR", "SAS")) anc else "EUR")
  src <- list(type = "vcf_file", file = facts$file_name %||% basename(path), trait = nm, imported = format(Sys.Date()))
  fullgwas_job_start(list(key = paste0("vcf_", nm), mode = "local", path = normalizePath(path, winslash = "/"), name = nm, facts = f, source = src), root)
}
