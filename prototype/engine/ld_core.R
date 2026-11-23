# Clumping on this computer (Round 19): PLINK 1.9 plus the 1000 Genomes Phase 3 LD reference panel of the MRC IEU,
# the same tool, panel and command that ieugwasr::ld_clump_local() and the OpenGWAS server use.
#   PLINK:  Purcell et al. 2007 (PMID 17701901); PLINK 1.9: Chang et al. 2015 (PMID 25722852). GPL-3; licence shipped next to the binary.
#   Panel:  1000 Genomes Project Consortium 2015 (PMID 26432245), MRC IEU build "1kg.v3" (one .bed/.bim/.fam per super-population).
# PLINK is bundled with the app (engine/tools/plink, put there by the build pipeline from a pinned, checksummed official build).
# The panel is an optional one-time download (1.57 GB) kept once per computer, outside every study folder.

PLINK_PIN <- list(version = "1.9.0 (27 Sep 2026)", build = "20260927",
                  url = "https://s3.amazonaws.com/plink1-assets/plink_%s_20260927.zip",
                  sha256 = c(mac = "927ecadbefb3e0ceefbe4978d71bf63566e9ea422d6595c2ee23a83ad76eeb62",
                             win64 = "625538e0062eb9fad6c26a47c9f36b19d3ec0afe5053ec68ff5ff877a56676c6",
                             linux_x86_64 = "ee4d8616dc64c04ed395e81a7bb04f9b3d92e4bf23d12a8b6233b19f7227eccd"))

LD_PANEL <- list(name = "1000 Genomes Phase 3 LD reference panel (MRC IEU build 1kg.v3)",
                 url = "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz", bytes = 1567805147,
                 sha256 = "06125160770fc65adf0296fe9a461011a910338bd106b4e319add13e81ab3670",
                 pops = c("EUR", "EAS", "AFR", "AMR", "SAS"))
# Pinned 2026-11-06 (Round 29). MRC IEU publishes no checksum, so the value was taken from two separate copies of 1kg.v3.tgz on the
# analyst's computer (file dates 2026-10-15), each 1,567,805,147 bytes, the Content-Length the server reports (its Last-Modified is
# 2020-05-07). Both copies have this SHA-256 and pass gzip -t, and their EUR .bed/.bim/.fam equal byte for byte the EUR files whose
# clumping matched ieugwasr and the OpenGWAS server (Round 28).
# SHA-256 of each unpacked file in that archive, so that panel files copied in from elsewhere (Round 24) can be compared too.
LD_FILE_SHA256 <- c(
  "AFR.bed" = "39d6f1fcc15caed279cd1aa07f3c92751f28588bcb277370c111d046a8ad54d2",
  "AFR.bim" = "94dacd3a98ec2c14c619086e6f72159e9847c3b06b4c06422b544727127a56fb",
  "AFR.fam" = "e1d9efed2d6ffa125ae0a302f00e3c271b42627d7da5d370ec9c32738045ed48",
  "AMR.bed" = "a3e45a95552a8439f773968d5af40e68ec0ea24c986aea1dc17c61443070b8bf",
  "AMR.bim" = "aa6983bc70a5708a4ec2beb16743f6c80994087fb5b1638c8bf1a4f927c1f679",
  "AMR.fam" = "8b55c3e54f16a86521e6a705e04f3d3b015fca23ad951b4aac27d92d80dd8439",
  "EAS.bed" = "8d00e7e094028499ab319def7cd5d7ca56541b07e09146ad50a0d7f739d1e8c9",
  "EAS.bim" = "e54e818b47f080caede23a1f48321c85ce3948be5db04e56bfc7e915235f5804",
  "EAS.fam" = "98ea977eb4a687df48c40a15f5c084da20a7b32467cc92b2620b75176ec3255f",
  "EUR.bed" = "ffe4adc0cd22d53397050ff8584323a31958bc4cbce281b1aeb935d5555e8363",
  "EUR.bim" = "fa1084f8104ea52c24fa71bedca04b28721e0e5519b38d723b65fa68c4b689b5",
  "EUR.fam" = "a6971a3b9b7fd344009619d7d2803d7c7d77aae2b1974cf6b1b29a0fbae293cd",
  "SAS.bed" = "8ca65d089cd71942bb4c4774d6d1870c0664198a6ab01a233529b798ff25df1a",
  "SAS.bim" = "cd3a0e8f7d9823c90914557655e50fa2dd2fe8f2aed6db9e170cef2bba6f79e0",
  "SAS.fam" = "2d307143700115cf0d299b7edcf4e21a3bba2de2a8706e32f293adc0d49f875c")
# People per super-population in 1000 Genomes Phase 3 (2,504 in all; Auton 2015, PMID 26432245). Used to recognise a full panel
# file among derived subsets (Round 24). All five counts match the .fam files of an MRC IEU 1kg.v3 copy on the analyst's computer (checked 2026-11-05).
LD_PEOPLE <- c(EUR = 503, EAS = 504, AFR = 661, AMR = 347, SAS = 489)
LD_POP_NAMES <- c(EUR = "European", EAS = "East Asian", AFR = "African", AMR = "Admixed American", SAS = "South Asian")

#' The bundled PLINK program (NA when it is missing). MRGUI_PLINK overrides it (tests, support).
plink_bin <- function() {
  e <- Sys.getenv("MRGUI_PLINK", "")
  if (nzchar(e)) return(if (file.exists(e)) normalizePath(e) else NA_character_)
  ed <- getOption("mrgui.engine_dir", "")
  f <- file.path(ed, "tools", "plink", if (.Platform$OS.type == "windows") "plink.exe" else "plink")
  if (nzchar(ed) && file.exists(f)) normalizePath(f) else NA_character_
}

#' Version line of PLINK, or NA when it does not run (for diagnostics and the report).
plink_version <- function(bin = plink_bin()) {
  if (is.na(bin)) return(NA_character_)
  v <- tryCatch(suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE)), error = function(e) character())
  v <- grep("PLINK", v, value = TRUE)
  if (length(v)) trimws(v[1]) else NA_character_
}

#' Where the panel is kept: once per computer, in the app's data folder (never inside a study or a synced folder).
ld_panel_dir <- function() {
  d <- Sys.getenv("MRGUI_LD_DIR", "")
  if (nzchar(d)) d else file.path(tools::R_user_dir("MRStudio", "data"), "ld_ref", "1kg.v3")
}

#' Which ancestries are installed (all three files present) and how much space they take.
ld_panel_status <- function(dir = ld_panel_dir()) {
  pops <- LD_PANEL$pops
  have <- vapply(pops, function(p) all(file.exists(file.path(dir, paste0(p, c(".bed", ".bim", ".fam"))))), TRUE)
  size <- vapply(pops, function(p) sum(file.size(file.path(dir, paste0(p, c(".bed", ".bim", ".fam")))), na.rm = TRUE), 0)
  st <- .ld_job_read(dir)
  if (!is.null(st) && st$state %in% c("downloading", "verifying", "unpacking", "copying") && !ld_job_active(dir)) { st$state <- "stopped"; st$msg <- if (identical(st$mode, "folder")) "Copying the panel files stopped (for example because MR Studio was closed). Choose the files again to restart." else if (identical(st$mode, "local")) "Checking and unpacking your file stopped (for example because MR Studio was closed). Choose the file again to restart." else "The download stopped (for example because MR Studio was closed). Start it again to continue where it stopped." }
  part <- file.path(dir, "1kg.v3.tgz.part")
  ver <- vapply(pops, function(p) have[[p]] && ld_panel_matches_pin(p, dir), TRUE)
  # Round 81: the job's "pinned" flag was written when the job ran (with the pins of that version); the note follows the current pins
  if (!is.null(st) && identical(st$state, "done") && identical(st$mode, "folder")) { jp <- intersect(unlist(st$pops), pops); if (length(jp)) st$pinned <- all(ver[jp]) }
  list(dir = dir, verified = as.list(setNames(ver, pops)), partial_bytes = if (file.exists(part)) file.size(part) else 0, installed = as.list(setNames(have, pops)), bytes = as.list(setNames(ifelse(have, size, 0), pops)),
       total_bytes = sum(size[have]), plink = plink_bin(), plink_version = plink_version(),
       download_bytes = LD_PANEL$bytes, job = st)
}

#' Round 29: does an installed ancestry match the pinned panel? Judged from the hashes recorded when it was installed: the archive's
#' SHA-256 in SOURCE.txt (download or own archive), or the three file hashes in SOURCE_<POP>.txt (files copied in). Nothing is re-hashed.
ld_panel_matches_pin <- function(pop, dir = ld_panel_dir(), pin = LD_PANEL$sha256, files = LD_FILE_SHA256) {
  arch <- function(x) { l <- grep("^sha256 of the archive: ", x, value = TRUE); length(l) > 0 && identical(sub(".*: ", "", l[1]), pin) }
  f1 <- file.path(dir, paste0("SOURCE_", pop, ".txt")); f0 <- file.path(dir, "SOURCE.txt")
  if (file.exists(f1)) { x <- readLines(f1, warn = FALSE); ext <- c(".bed", ".bim", ".fam")
    got <- vapply(ext, function(e) { l <- grep(paste0("^sha256 ", pop, "\\", e, ": "), x, value = TRUE); if (length(l)) sub(".*: ", "", l[1]) else "" }, "")
    if (any(nzchar(got))) return(identical(unname(got), unname(files[paste0(pop, ext)])))
    return(arch(x)) }
  if (file.exists(f0)) { x <- readLines(f0, warn = FALSE); a <- grep("^ancestries: ", x, value = TRUE)   # installs before Round 29
    return(arch(x) && (!length(a) || pop %in% strsplit(sub("^ancestries: ", "", a[1]), " ")[[1]])) }
  FALSE
}

.ld_job_file <- function(dir) file.path(dir, "download_status.json")
.ld_job_read <- function(dir) {
  f <- .ld_job_file(dir); if (!file.exists(f)) return(NULL)
  st <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL); if (is.null(st)) return(NULL)
  part <- file.path(dir, "1kg.v3.tgz.part")
  if (identical(st$state, "downloading") && file.exists(part)) st$bytes <- file.size(part)
  if (!is.null(st$pops)) st$pops <- as.list(st$pops)   # always a list for the screen, also for one ancestry
  st
}

# The .bed check shared by the background jobs: PLINK signature, and size = 3 + variants x ceiling(people / 4) from .bim and .fam.
LD_BEDCHECK_R <- c(
    "bedcheck <- function(b) {   # PLINK .bed: 3 signature bytes, then ceiling(people / 4) bytes per variant",
    "  nl <- function(f) { con <- file(f, 'rb'); on.exit(close(con)); n <- 0; repeat { x <- readBin(con, 'raw', 1e7); if (!length(x)) break; n <- n + sum(x == as.raw(10)) }; n }",
    "  sig <- readBin(paste0(b, '.bed'), 'raw', 3); nv <- nl(paste0(b, '.bim')); ni <- nl(paste0(b, '.fam'))",
    "  if (!identical(sig, as.raw(c(0x6c, 0x1b, 0x01)))) stop(basename(b), '.bed is not a PLINK genotype file (wrong signature). The archive may be damaged.')",
    "  if (file.size(paste0(b, '.bed')) != 3 + nv * ceiling(ni / 4)) stop(basename(b), '.bed does not match its .bim and .fam (', nv, ' variants, ', ni, ' people). The archive may be damaged.')",
    "  c(variants = nv, people = ni) }")

#' Start the one-time panel download in a separate R process, so the app stays usable. The download resumes where it stopped
#' (also after the app was closed), retries when the connection stalls, checks the size (and the SHA-256 once it is pinned), and
#' unpacks only the ancestries asked for (each is about 2 GB unpacked; all five need about 10 GB). Progress is read from files.
#' With `local_file`, nothing is downloaded: the archive the user downloaded in a browser is checked and unpacked in place
#' (Round 23). The user's file is never moved or deleted.
ld_panel_download_start <- function(pops = "EUR", dir = ld_panel_dir(), url = LD_PANEL$url, bytes = LD_PANEL$bytes, sha256 = LD_PANEL$sha256, local_file = NULL) {
  pops <- intersect(toupper(pops), LD_PANEL$pops); if (!length(pops)) stop("Choose at least one ancestry to install.")
  if (!is.null(local_file)) { p <- ld_archive_problem(local_file, bytes); if (!is.null(p)) stop(p, call. = FALSE); local_file <- normalizePath(local_file, winslash = "/") }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (ld_job_active(dir)) return(invisible(.ld_job_read(dir)))
  job <- file.path(dir, "download_job.R")
  writeLines(c(
    "args <- commandArgs(TRUE); dir <- args[1]; url <- args[2]; bytes <- as.numeric(args[3]); sha <- args[4]; pops <- strsplit(args[5], ',')[[1]]; loc <- if (length(args) >= 6) args[6] else ''",
    "mode <- if (nzchar(loc)) 'local' else 'download'",
    "st <- function(state, msg = '', ...) writeLines(jsonlite::toJSON(c(list(state = state, msg = msg, pops = pops, mode = mode, time = format(Sys.time(), '%Y-%m-%d %H:%M:%S')), list(...)), auto_unbox = TRUE), file.path(dir, 'download_status.json'))",
    "part <- if (mode == 'local') loc else file.path(dir, '1kg.v3.tgz.part')",
    LD_BEDCHECK_R,
    "fetch_append <- function(url, part, h) { con <- curl::curl(url, 'rb', handle = h); on.exit(close(con)); out <- file(part, 'ab'); on.exit(close(out), add = TRUE); repeat { b <- readBin(con, 'raw', 1048576L); if (!length(b)) break; writeBin(b, out) } }   # Round 75: curl_download(mode = 'ab') replaced the part instead of adding to it",
    "tryCatch({",
    "  st(if (mode == 'local') 'verifying' else 'downloading')",
    "  if (mode == 'download') for (attempt in 1:50) {   # resume after a stall or a dropped connection",
    "    have <- if (file.exists(part)) file.size(part) else 0",
    "    if (have >= bytes) break",
    "    h <- curl::new_handle(resume_from = have, low_speed_limit = 100, low_speed_time = 120, connecttimeout = 60)",
    "    ok <- tryCatch({ fetch_append(url, part, h); TRUE }, error = function(e) FALSE)",
    "    if (!ok) Sys.sleep(10)",
    "  }",
    "  st('verifying')",
    "  if (file.size(part) != bytes) stop(if (mode == 'local') sprintf('The file has %s bytes, not %s: it is incomplete or a different file. Let the browser finish, or download it again.', file.size(part), bytes) else sprintf('The download stopped at %s of %s bytes. Start it again to continue where it stopped.', file.size(part), bytes))",
    "  h <- digest::digest(file = part, algo = 'sha256')",
    "  if (!identical(sha, 'PANEL_SHA256_PENDING') && !identical(h, sha)) { if (mode == 'download') unlink(part); stop(if (mode == 'local') 'The file does not match the published panel (SHA-256 differs), so it was not used. Your file was left where it is.' else 'The downloaded file does not match the published panel (SHA-256 differs). It was deleted and not used. Please try again.') }",
    "  st('unpacking'); tmp <- file.path(dir, 'unpack'); unlink(tmp, recursive = TRUE); dir.create(tmp)",
    "  want <- as.vector(outer(pops, c('.bed', '.bim', '.fam'), paste0))",
    "  utils::untar(part, files = want, exdir = tmp, tar = 'internal')",
    "  got <- file.exists(file.path(tmp, want)); if (!all(got)) stop('The archive did not hold ', paste(want[!got], collapse = ', '), '.')",
    "  dims <- lapply(pops, function(p) bedcheck(file.path(tmp, p)))",
    "  for (f in want) file.rename(file.path(tmp, f), file.path(dir, f))",
    "  unlink(tmp, recursive = TRUE); if (mode == 'download') unlink(part)",
    "  chk <- if (identical(sha, 'PANEL_SHA256_PENDING')) 'size and file structure checked; no official SHA-256 pinned yet, so the hash is recorded only' else 'size, SHA-256 and file structure checked'",
    "  writeLines(c(paste('source:', if (mode == 'local') paste0('archive downloaded by the user in a browser (', loc, '), from ', url) else url), paste('sha256 of the archive:', h), paste('checked:', chk), paste('installed:', format(Sys.time(), '%Y-%m-%d')), paste('ancestries:', paste(pops, collapse = ' ')),",
    "    vapply(seq_along(pops), function(i) sprintf('%s: %s variants, %s people', pops[i], format(dims[[i]][['variants']], big.mark = ','), dims[[i]][['people']]), '')), file.path(dir, 'SOURCE.txt'))",
    "  for (i in seq_along(pops)) writeLines(c(paste('source: MRC IEU archive 1kg.v3.tgz,', if (mode == 'local') loc else url), paste('sha256 of the archive:', h), paste('checked:', chk),",
    "    sprintf('%s: %s variants, %s people', pops[i], format(dims[[i]][['variants']], big.mark = ','), dims[[i]][['people']]), paste('installed:', format(Sys.time(), '%Y-%m-%d'))), file.path(dir, paste0('SOURCE_', pops[i], '.txt')))   # one record per ancestry (Round 29)",
    "  st('done', sha256 = h, pinned = !identical(sha, 'PANEL_SHA256_PENDING'))",
    "}, error = function(e) st('error', conditionMessage(e)))"), job)
  writeLines(jsonlite::toJSON(list(state = if (is.null(local_file)) "downloading" else "verifying", msg = "", pops = pops, mode = if (is.null(local_file)) "download" else "local", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), .ld_job_file(dir))
  old <- Sys.getenv("RHOME", unset = NA); Sys.setenv(RHOME = R.home())
  on.exit(if (is.na(old)) Sys.unsetenv("RHOME") else Sys.setenv(RHOME = old), add = TRUE)
  system2(default_rscript(), c("--vanilla", shQuote(job), shQuote(dir), shQuote(url), format(bytes, scientific = FALSE), sha256, paste(pops, collapse = ","), if (!is.null(local_file)) shQuote(local_file)),
          stdout = file.path(dir, "download.log"), stderr = file.path(dir, "download.log"), wait = FALSE)
  invisible(.ld_job_read(dir))
}

#' Why a file the user points to cannot be the panel archive (NULL when it can be used).
ld_archive_problem <- function(path, bytes = LD_PANEL$bytes) {
  if (is.null(path) || !nzchar(trimws(path))) return("Give the path of the downloaded 1kg.v3.tgz file.")
  path <- path.expand(trimws(gsub('^["\']|["\']$', "", trimws(path))))
  if (!file.exists(path)) return("No file was found at that path. On a Mac, select the file in Finder and press Option+Cmd+C to copy its path.")
  if (dir.exists(path)) return(if (grepl("\\.download$", path)) "This is an unfinished Safari download. Wait until Safari has finished, then use the 1kg.v3.tgz file." else "This is a folder, not the archive. Your browser may have unpacked it: point to the 1kg.v3.tgz file itself (in Safari, turn off \"Open safe files after downloading\" and download again).")
  if (grepl("\\.(crdownload|part|download)$", path)) return("This is an unfinished browser download. Wait until the browser has finished, then use the 1kg.v3.tgz file.")
  if (grepl("\\.tar$", path)) return("This is the unpacked .tar, not the .tgz archive (Safari unpacks archives by itself). Turn off \"Open safe files after downloading\" in Safari's settings and download again, or use another browser.")
  b <- file.size(path)
  if (!isTRUE(b == bytes)) return(sprintf("The file has %s bytes, but the panel archive has %s. It is incomplete or a different file.", format(b, big.mark = ","), format(bytes, big.mark = ",")))
  NULL
}

#' Panel archives already on this computer, in the usual places a browser saves to (newest first), each with a verdict.
ld_archive_find <- function(places = path.expand(c("~/Downloads", "~/Desktop", "~")), bytes = LD_PANEL$bytes) {
  fs <- unique(unlist(lapply(places[dir.exists(places)], function(d) list.files(d, pattern = "^1kg\\.v3", full.names = TRUE, include.dirs = TRUE))))
  if (!length(fs)) return(list())
  fs <- fs[order(file.mtime(fs), decreasing = TRUE)]
  lapply(fs, function(f) { p <- ld_archive_problem(f, bytes); list(path = normalizePath(f, winslash = "/"), name = basename(f), bytes = if (dir.exists(f)) 0 else file.size(f), ok = is.null(p), why = p %||% "") })
}

#' Round 24: an ancestry the user already has unpacked (PLINK .bed/.bim/.fam, for example from an earlier project). `path` is
#' the .bed file or the folder holding it. The ancestry is read from the file name (EUR.bed, 1000G_EUR.bed, ...) or given.
#' Returns list(bed prefix, pop, problem). Only quick checks here; the job checks the structure and copies.
ld_unpacked_problem <- function(path, pop = NULL, people = LD_PEOPLE) {
  if (is.null(path) || !nzchar(trimws(path))) return(list(problem = "Give the path of the EUR.bed file (or the folder that holds it)."))
  path <- path.expand(trimws(gsub('^["\']|["\']$', "", trimws(path))))
  if (!file.exists(path)) return(list(problem = "No file or folder was found at that path. On a Mac, select it in Finder and press Option+Cmd+C to copy its path."))
  if (dir.exists(path)) {
    beds <- list.files(path, pattern = "\\.bed$", full.names = TRUE)
    code <- function(f) { t <- toupper(strsplit(sub("\\.bed$", "", basename(f)), "[^A-Za-z]+")[[1]]); t[t %in% names(people)] }
    beds <- beds[vapply(beds, function(f) length(code(f)) == 1L, TRUE)]
    if (!is.null(pop)) beds <- beds[vapply(beds, function(f) identical(code(f), toupper(pop)), TRUE)]
    if (!length(beds)) return(list(problem = "This folder has no .bed file named after a 1000 Genomes ancestry (for example EUR.bed)."))
    if (length(beds) > 1L) return(list(problem = paste0("This folder has several panel files (", paste(basename(beds), collapse = ", "), "). Give the path of the one you want.")))
    path <- beds[1]
  }
  if (!grepl("\\.bed$", path)) return(list(problem = "Give the .bed file of the panel (it needs its .bim and .fam next to it)."))
  b <- sub("\\.bed$", "", path); miss <- c(".bim", ".fam")[!file.exists(paste0(b, c(".bim", ".fam")))]
  if (length(miss)) return(list(problem = paste0("The ", paste(miss, collapse = " and "), " file next to ", basename(path), " is missing. All three files (.bed, .bim, .fam) are needed.")))
  t <- toupper(strsplit(basename(b), "[^A-Za-z]+")[[1]]); cd <- t[t %in% names(people)]
  p <- if (length(cd) == 1L) cd else toupper(pop %||% "")
  if (!(p %in% names(people))) return(list(problem = "The ancestry cannot be read from the file name. Rename it after its 1000 Genomes ancestry (for example EUR.bed) or choose the ancestry."))
  if (!is.null(pop) && nzchar(pop) && length(cd) == 1L && !identical(cd, toupper(pop))) return(list(problem = paste0("The file name says ", cd, " but ", toupper(pop), " was chosen.")))
  ni <- length(readLines(paste0(b, ".fam"), warn = FALSE))
  if (ni != people[[p]]) return(list(problem = sprintf("%s.fam lists %d people, but the full 1000 Genomes %s panel has %d. This looks like a subset made from the panel (for example by an earlier analysis), not the panel itself, so it is not used.", basename(b), ni, LD_POP_NAMES[[p]], people[[p]])))
  list(bed = normalizePath(b, winslash = "/", mustWork = FALSE), pop = p, people = ni, problem = NULL)
}

#' Unpacked panel files already on this computer. On macOS Spotlight finds them anywhere; elsewhere the usual folders are listed.
ld_unpacked_find <- function(people = LD_PEOPLE) {
  pat <- paste0("^(", paste(names(people), collapse = "|"), ")\\.bed$")
  fs <- character()
  if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("mdfind")))
    fs <- tryCatch(suppressWarnings(system2("mdfind", c("-name", ".bed"), stdout = TRUE, stderr = FALSE)), error = function(e) character())
  places <- path.expand(c("~/Downloads", "~/Desktop", "~/Documents"))
  fs <- c(fs, unlist(lapply(places[dir.exists(places)], function(d) list.files(d, pattern = pat, full.names = TRUE))))
  fs <- unique(fs[grepl(pat, basename(fs)) & file.exists(fs)])
  fs <- fs[!grepl("test_fixtures|/Library/|/\\.Trash/", fs)]
  dir_panel <- normalizePath(ld_panel_dir(), winslash = "/", mustWork = FALSE); fs <- fs[!startsWith(normalizePath(fs, winslash = "/"), dir_panel)]
  if (!length(fs)) return(list())
  fs <- fs[order(grepl("/Documents/|Mobile Documents", fs), file.mtime(fs), decreasing = c(FALSE, TRUE), method = "radix")]
  lapply(fs, function(f) { r <- ld_unpacked_problem(f, people = people); list(path = normalizePath(f, winslash = "/"), pop = r$pop %||% "", ok = is.null(r$problem), why = r$problem %||% "",
    icloud = grepl("/Documents/|/Desktop/|Mobile Documents", f)) })
}

#' Copy an unpacked ancestry into the panel folder in a background job: structure checked before and after copying, SHA-256 of
#' each copied file recorded in SOURCE_<POP>.txt. The user's files are never moved or changed.
#' Round 29: the three hashes are compared with the files of the pinned archive (`pins`). A difference is recorded and shown, not refused:
#' the structure and people checks already passed, and a copy from another build may be what the analyst intends.
ld_panel_copy_start <- function(path, pop = NULL, dir = ld_panel_dir(), people = LD_PEOPLE, pins = LD_FILE_SHA256) {
  r <- ld_unpacked_problem(path, pop, people); if (!is.null(r$problem)) stop(r$problem, call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (ld_job_active(dir)) stop("Another panel job is still working. Wait until it has finished.", call. = FALSE)
  job <- file.path(dir, "copy_job.R")
  writeLines(c(
    "args <- commandArgs(TRUE); dir <- args[1]; src <- args[2]; pop <- args[3]; pin <- strsplit(args[4], ',')[[1]]; pops <- pop; mode <- 'folder'",
    "st <- function(state, msg = '', ...) writeLines(jsonlite::toJSON(c(list(state = state, msg = msg, pops = pops, mode = mode, time = format(Sys.time(), '%Y-%m-%d %H:%M:%S')), list(...)), auto_unbox = TRUE), file.path(dir, 'download_status.json'))",
    LD_BEDCHECK_R,
    "tryCatch({",
    "  st('verifying'); d0 <- bedcheck(src)",
    "  st('copying'); tmp <- file.path(dir, 'unpack'); unlink(tmp, recursive = TRUE); dir.create(tmp)",
    "  ext <- c('.bed', '.bim', '.fam')",
    "  for (e in ext) if (!file.copy(paste0(src, e), file.path(tmp, paste0(pop, e)), overwrite = TRUE)) stop('Could not copy ', basename(src), e, '. Is there enough free disk space?')",
    "  if (any(file.size(file.path(tmp, paste0(pop, ext))) != file.size(paste0(src, ext)))) stop('A copied file has a different size from the original. Check the free disk space and try again.')",
    "  d1 <- bedcheck(file.path(tmp, pop))",
    "  h <- vapply(ext, function(e) digest::digest(file = file.path(tmp, paste0(pop, e)), algo = 'sha256'), '')",
    "  for (e in ext) file.rename(file.path(tmp, paste0(pop, e)), file.path(dir, paste0(pop, e)))",
    "  unlink(tmp, recursive = TRUE)",
    "  same <- length(pin) == 3L && identical(unname(h), pin)",
    "  writeLines(c(paste('source: unpacked panel files copied from', src, '(.bed, .bim, .fam)'), paste('checked: PLINK file structure and number of people (full 1000 Genomes Phase 3 super-population);', if (same) 'SHA-256 of all three files identical to the files in the pinned MRC IEU 1kg.v3 archive' else 'SHA-256 differs from the files in the pinned MRC IEU 1kg.v3 archive (another build or an edited copy), so clumping may differ from OpenGWAS'),",
    "    paste0('sha256 ', pop, ext, ': ', h), sprintf('%s: %s variants, %s people', pop, format(d1[['variants']], big.mark = ','), d1[['people']]), paste('installed:', format(Sys.time(), '%Y-%m-%d'))), file.path(dir, paste0('SOURCE_', pop, '.txt')))",
    "  st('done', sha256 = unname(h[1]), pinned = same, from = src)",
    "}, error = function(e) { unlink(file.path(dir, 'unpack'), recursive = TRUE); st('error', conditionMessage(e)) })"), job)
  writeLines(jsonlite::toJSON(list(state = "verifying", msg = "", pops = r$pop, mode = "folder", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), .ld_job_file(dir))
  old <- Sys.getenv("RHOME", unset = NA); Sys.setenv(RHOME = R.home())
  on.exit(if (is.na(old)) Sys.unsetenv("RHOME") else Sys.setenv(RHOME = old), add = TRUE)
  system2(default_rscript(), c("--vanilla", shQuote(job), shQuote(dir), shQuote(r$bed), r$pop, paste(c(pins[paste0(r$pop, c(".bed", ".bim", ".fam"))], recursive = TRUE), collapse = ",")), stdout = file.path(dir, "copy.log"), stderr = file.path(dir, "copy.log"), wait = FALSE)
  invisible(.ld_job_read(dir))
}

#' Is a download job still working? Judged by its files changing recently (the .part file grows while it downloads).
ld_job_active <- function(dir = ld_panel_dir()) {
  st <- .ld_job_read(dir); if (is.null(st) || !(st$state %in% c("downloading", "verifying", "unpacking", "copying"))) return(FALSE)
  fs <- c(.ld_job_file(dir), file.path(dir, "1kg.v3.tgz.part"), file.path(dir, "unpack")); fs <- fs[file.exists(fs)]
  any(difftime(Sys.time(), file.mtime(fs), units = "mins") < 10)
}

#' Remove one ancestry (or all) to free disk space.
ld_panel_remove <- function(pops = LD_PANEL$pops, dir = ld_panel_dir()) {
  pops <- intersect(pops, LD_PANEL$pops)
  unlink(file.path(dir, c(as.vector(outer(pops, c(".bed", ".bim", ".fam"), paste0)), paste0("SOURCE_", pops, ".txt"))))
  invisible(ld_panel_status(dir))
}

#' The super-population to clump with, from a declared ancestry (NA when it has no single 1000 Genomes match, e.g. "Mixed").
ld_pop_for <- function(ancestry) { a <- toupper(as.character(ancestry %||% "")); if (a %in% LD_PANEL$pops) a else NA_character_ }

#' Clump candidate instruments on this computer. Same PLINK call as ieugwasr::ld_clump_local() (and the OpenGWAS server):
#'   plink --bfile <panel/POP> --clump <file of SNP, P> --clump-p1 1 --clump-p2 1 --clump-r2 r2 --clump-kb kb
#' The candidates are already below the instrument P threshold, so p1 = p2 = 1 keeps every one of them in play (as in
#' TwoSampleMR::clump_data). Returns every candidate with its status: "index" (kept), "clumped" (with the index variant it is in
#' LD with) or "not_in_panel" (absent from the reference, so its LD cannot be judged; dropped, as ieugwasr does).
#' Ties (Round 28): P values below what OpenGWAS stores (1e-200) or below about 1e-308 (stored as 0) are equal. PLINK then keeps the
#' tied variant that comes first in the reference panel (by position), whatever the order of the input. tie_order = "server" (the
#' default) gives PLINK the P values as they are, exactly as ieugwasr::ld_clump_local() and the OpenGWAS server do, so the instruments
#' equal those of the reference tools (validation study 2: 79 of 79 on ieu-a-300).
#' tie_order = "strongest" (the tie sensitivity analysis) breaks equal P by the larger |beta/SE| (`z`): PLINK is given the rank of each
#' variant by (P, then larger |z| first) scaled into (0, 1). PLINK uses P only to order the variants (p1 = p2 = 1 keeps them all), so
#' without ties both orders give exactly the result of clumping on P.
ld_clump_here <- function(snp, p, r2 = 0.001, kb = 10000, pop = "EUR", dir = ld_panel_dir(), bin = plink_bin(), z = NULL, tie_order = c("server", "strongest")) {
  tie_order <- match.arg(tie_order)
  if (identical(tie_order, "strongest") && is.null(z)) stop("The tie sensitivity analysis needs beta / SE for every candidate.", call. = FALSE)
  if (is.na(bin)) stop("PLINK is not available in this installation, so variants cannot be clumped on this computer.", call. = FALSE)
  pop <- ld_pop_for(pop)
  if (is.na(pop)) stop("Clumping needs one of the 1000 Genomes ancestries (European, East Asian, African, Admixed American, South Asian).", call. = FALSE)
  bfile <- file.path(dir, pop)
  if (!all(file.exists(paste0(bfile, c(".bed", ".bim", ".fam")))))
    stop(sprintf("The %s (%s) LD reference panel is not installed. Download it in Settings > Data and access.", LD_POP_NAMES[[pop]], pop), call. = FALSE)
  ok <- !is.na(snp) & grepl("^rs[0-9]+$", snp) & is.finite(p)
  if (!any(ok)) stop("Clumping needs rsIDs (for example rs12345). This dataset identifies variants by position only.", call. = FALSE)
  if (anyDuplicated(snp[ok])) stop("The candidate variants contain duplicated rsIDs. Each variant may appear only once.", call. = FALSE)
  tmp <- tempfile("clump-"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  inp <- file.path(tmp, "candidates.txt"); out <- file.path(tmp, "clumped")
  ties <- identical(tie_order, "strongest") && anyDuplicated(p[ok]) > 0
  pp <- if (ties) { zz <- abs(z[ok]); zz[!is.finite(zz)] <- 0; o <- order(p[ok], -zz, method = "radix"); rk <- integer(sum(ok)); rk[o] <- seq_along(o); rk / (sum(ok) + 1) } else p[ok]
  utils::write.table(data.frame(SNP = snp[ok], P = formatC(pp, digits = 15, format = "g")), inp, row.names = FALSE, col.names = TRUE, quote = FALSE)
  args <- c("--bfile", shQuote(bfile), "--clump", shQuote(inp), "--clump-p1", "1", "--clump-p2", "1",
            "--clump-r2", format(r2, scientific = FALSE), "--clump-kb", format(kb, scientific = FALSE), "--threads", "1", "--out", shQuote(out))
  log <- suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE))
  st <- attr(log, "status"); if (is.null(st)) st <- 0L
  res <- data.frame(variant_id = snp, p = p, status = ifelse(ok, "not_in_panel", "no_rsid"), index_variant = NA_character_, stringsAsFactors = FALSE)
  cf <- paste0(out, ".clumped")
  if (file.exists(cf)) {
    cl <- utils::read.table(cf, header = TRUE, stringsAsFactors = FALSE)
    res$status[res$variant_id %in% cl$SNP] <- "index"
    for (i in seq_len(nrow(cl))) {
      mem <- sub("\\(.*$", "", trimws(strsplit(cl$SP2[i], ",")[[1]])); mem <- mem[mem != "NONE" & nzchar(mem)]
      hit <- res$variant_id %in% mem & res$status == "not_in_panel"
      res$status[hit] <- "clumped"; res$index_variant[hit] <- cl$SNP[i]
    }
  } else if (!any(grepl("No significant --clump results", log))) {
    stop("PLINK could not clump the variants (exit status ", st, "). End of its log: ", paste(utils::tail(log, 5), collapse = " "), call. = FALSE)
  }
  res$index_variant[res$status == "index"] <- res$variant_id[res$status == "index"]
  attr(res, "settings") <- list(r2 = r2, kb = kb, pop = pop, panel = LD_PANEL$name, panel_matches_pin = ld_panel_matches_pin(pop, dir), plink = plink_version(bin), tie_break = if (!anyDuplicated(p[ok])) "none needed" else if (ties) "equal P ordered by |beta/SE|" else "equal P kept in the order of the reference panel (PLINK), as ieugwasr and the OpenGWAS server",
                                tied = .n_tied(p[ok]), tie_order = tie_order,
                                command = paste("plink --bfile", pop, "--clump candidates.txt --clump-p1 1 --clump-p2 1 --clump-r2", format(r2, scientific = FALSE), "--clump-kb", format(kb, scientific = FALSE)))
  res
}

#' Clumping result for a dataset, cached in the dataset's library folder (the same inputs never run PLINK twice).
ld_clump_cached <- function(ds_dir, snp, p, r2, kb, pop, z = NULL, tie_order = "server") {
  key <- substr(digest::digest(list(snp = snp, p = p, r2 = r2, kb = kb, pop = pop, panel = LD_PANEL$sha256, plink = PLINK_PIN$build, tie_order = tie_order, z = if (identical(tie_order, "strongest") && anyDuplicated(p)) signif(z, 10)), algo = "sha256"), 1, 16)
  f <- file.path(ds_dir, "clumping", paste0(key, ".csv")); fs <- file.path(ds_dir, "clumping", paste0(key, ".json"))
  if (file.exists(f) && file.exists(fs)) {
    res <- utils::read.csv(f, stringsAsFactors = FALSE, colClasses = c(variant_id = "character", status = "character", index_variant = "character"))
    attr(res, "settings") <- jsonlite::fromJSON(fs); return(res)
  }
  res <- ld_clump_here(snp, p, r2, kb, pop, z = z, tie_order = tie_order)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, f, row.names = FALSE); writeLines(jsonlite::toJSON(attr(res, "settings"), auto_unbox = TRUE), fs)
  res
}

#' Number of candidates whose P value equals that of another candidate.
.n_tied <- function(p) { p <- p[is.finite(p)]; sum(duplicated(p) | duplicated(p, fromLast = TRUE)) }

#' Which lead variants a tie decided (Round 28): the instruments of the primary selection (`main`, tie_order = "server") compared with
#' those of a selection that breaks equal P by the larger |beta/SE| (`alt`). `out` are leads only the primary keeps, `in` leads only
#' the alternative keeps; `regions` is the larger of the two counts. NULL when no candidate shares its P value.
ld_tie_summary <- function(main, alt) {
  if (is.null(alt)) return(NULL)
  a <- main$variant_id[main$status == "index"]; b <- alt$variant_id[alt$status == "index"]
  list(tied = .n_tied(main$p[main$status %in% c("index", "clumped")]), regions = max(length(setdiff(a, b)), length(setdiff(b, a))),
       out = as.list(setdiff(a, b)), `in` = as.list(setdiff(b, a)), alt_kept = length(b))
}

# ---- Round 33: allele-frequency check against the reference panel (data_limits_solutions.md item 5, FR-91) ----------------
# The effect-allele frequency (EAF) a dataset reports is compared with the frequency of the same allele in the 1000 Genomes
# panel of the dataset's ancestry. A difference above EAF_FLIP_THRESHOLD marks a possible allele or strand error (the EAF may
# belong to the other allele). Calibrated in Round 32 (validation/calibration/README.md): no false flag among 2,320 European
# variants, every flip found below a minor allele frequency of 0.3; near 0.5 a flip cannot be seen, so variants with a panel
# minor allele frequency of EAF_FLIP_MAF_LIMIT or more are not judged. Flagged variants are reported, never removed.
EAF_FLIP_THRESHOLD <- 0.2
EAF_FLIP_MAF_LIMIT <- 0.4

#' Panel allele frequencies of the given rsIDs (PLINK --freq), cached in the dataset's folder. Columns SNP, A1, A2, MAF, NCHROBS.
ld_panel_freq_cached <- function(ds_dir, snp, pop, dir = ld_panel_dir(), bin = plink_bin()) {
  snp <- sort(unique(snp[!is.na(snp) & grepl("^rs[0-9]+$", snp)]))
  empty <- data.frame(SNP = character(), A1 = character(), A2 = character(), MAF = numeric(), NCHROBS = integer(), stringsAsFactors = FALSE)
  if (!length(snp)) return(empty)
  key <- substr(digest::digest(list(snp = snp, pop = pop, panel = LD_PANEL$sha256, plink = PLINK_PIN$build)), 1, 16)
  f <- file.path(ds_dir, "panel_freq", paste0(pop, "_", key, ".csv"))
  cc <- c(SNP = "character", A1 = "character", A2 = "character")
  if (file.exists(f)) return(utils::read.csv(f, stringsAsFactors = FALSE, colClasses = cc))
  tmp <- tempfile("freq-"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  writeLines(snp, file.path(tmp, "snps.txt"))
  log <- suppressWarnings(system2(bin, c("--bfile", shQuote(file.path(dir, pop)), "--extract", shQuote(file.path(tmp, "snps.txt")), "--freq", "--threads", "1", "--out", shQuote(file.path(tmp, "frq"))), stdout = TRUE, stderr = TRUE))
  ff <- file.path(tmp, "frq.frq")
  res <- if (file.exists(ff)) { r <- utils::read.table(ff, header = TRUE, stringsAsFactors = FALSE, colClasses = cc); r[is.finite(r$MAF), c("SNP", "A1", "A2", "MAF", "NCHROBS")] }
         else if (any(grepl("No variants remaining", log))) empty
         else stop("PLINK could not read the panel frequencies. End of its log: ", paste(utils::tail(log, 5), collapse = " "), call. = FALSE)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); utils::write.csv(res, f, row.names = FALSE)
  res
}

#' Compare the EAF of one side of a harmonised pair with the panel. `d` is the harmonised table (ea, oa aligned to the exposure),
#' `eaf` the EAF column of that side (already aligned to `ea`). Returns the summary for the ledger and one row per judged variant.
eaf_panel_check <- function(d, eaf, ds_dir, ancestry, trait, dir = ld_panel_dir(), bin = plink_bin()) {
  base <- list(trait = trait, threshold = EAF_FLIP_THRESHOLD, maf_limit = EAF_FLIP_MAF_LIMIT)
  no <- function(why) c(base, list(run = FALSE, why = why))
  kept <- d$mr_keep %in% TRUE; e <- eaf; e[!kept] <- NA_real_
  if (!any(is.finite(e))) return(no(sprintf("\"%s\" reports no effect-allele frequency.", trait)))
  pop <- ld_pop_for(ancestry)
  if (is.na(pop)) return(no(sprintf("Its ancestry (%s) is not one of the five 1000 Genomes ancestries.", ancestry %||% "unknown")))
  if (is.na(bin)) return(no("PLINK is not available in this installation."))
  if (!all(file.exists(file.path(dir, paste0(pop, c(".bed", ".bim", ".fam"))))))
    return(no(sprintf("The %s (%s) reference panel is not installed (Settings > Data and access).", LD_POP_NAMES[[pop]], pop)))
  fq <- ld_panel_freq_cached(ds_dir, d$variant_id[kept], pop, dir, bin)
  pal <- (d$ea %in% c("A", "T") & d$oa %in% c("A", "T")) | (d$ea %in% c("C", "G") & d$oa %in% c("C", "G"))
  m <- match(d$variant_id, fq$SNP); a1 <- fq$A1[m]; a2 <- fq$A2[m]
  pan <- ifelse(toupper(a1) == d$ea & toupper(a2) == d$oa, fq$MAF[m], ifelse(toupper(a2) == d$ea & toupper(a1) == d$oa, 1 - fq$MAF[m], NA_real_))
  maf <- pmin(pan, 1 - pan); judged <- kept & is.finite(e) & !pal & is.finite(pan) & maf < EAF_FLIP_MAF_LIMIT
  rows <- data.frame(variant_id = d$variant_id, side = trait, ea = d$ea, oa = d$oa, eaf = e, panel_eaf = pan, difference = e - pan,
                     judged = judged, flagged = judged & abs(e - pan) > EAF_FLIP_THRESHOLD, stringsAsFactors = FALSE)[kept & is.finite(e), , drop = FALSE]
  c(base, list(run = TRUE, pop = pop, panel = LD_PANEL$name, with_eaf = sum(kept & is.finite(e)), judged = sum(judged), flagged = sum(rows$flagged),
               flagged_ids = rows$variant_id[rows$flagged], palindromic = sum(kept & is.finite(e) & pal), not_in_panel = sum(kept & is.finite(e) & !pal & !is.finite(pan)),
               near_half = sum(kept & is.finite(e) & !pal & is.finite(pan) & maf >= EAF_FLIP_MAF_LIMIT),
               max_difference = if (any(judged)) max(abs(e - pan)[judged]) else NA_real_), list(rows = rows))
}
