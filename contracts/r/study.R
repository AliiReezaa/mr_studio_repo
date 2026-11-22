# MR GUI contracts v1.1 (DRAFT) - visible study folder, manifest and the script task (C6.5, C8, C9)
#
# Engine-side code. Requires canonical.R, adapter.R, runner.R and script.R to be sourced first.
# The study folder is what the user sees and keeps; the hidden `.engine/` folder inside it holds
# project.sqlite and the blobs (unchanged from contracts v1.0.0).

# Round 96: 05_report_and_supplement and 06_reproducibility were never written (the package goes to 05_publication), so they
# are no longer created; folders a study uses only sometimes are created when first written.
# Round 98: 02_datasets, 04_across_pairs and archive are used only by some studies, so they too are created when first written.
STUDY_FOLDERS <- c("01_study_plan", "03_pairs")
# Round 96: 05_figures and 06_method_reports were never written and 08_outlier_analysis is created by the outlier analysis
# itself, so a pair no longer gets them empty. Folder names stay as they were, so earlier studies and scripts keep working.
PAIR_SUBFOLDERS <- c("01_data", "02_scripts", "03_console", "04_tables", "07_replay_check")
MAX_PATH_CHARS <- 240L          # Windows limit with a safety margin (C9)
MANIFEST_COLS <- c("relative_path", "sha256", "bytes", "role", "task_id", "created_by", "version")

# ---- names and paths (FR-54) -------------------------------------------------
#' ASCII slug, at most `max` characters, never empty.
slugify <- function(x, max = 40L) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("[^A-Za-z0-9]+", "-", x)
  x <- gsub("^-+|-+$", "", x)
  x <- substr(x, 1L, max)
  x <- gsub("-+$", "", x)
  if (is.na(x) || !nzchar(x)) "x" else x
}

#' Refuse folders the contracts forbid (D2): cloud-sync folders and network drives.
study_path_problem <- function(path) {
  p <- gsub("\\\\", "/", path)
  if (grepl("^//", p)) return("network drive (UNC path)")
  cloud <- c("Library/Mobile Documents", "iCloud", "OneDrive", "Dropbox", "Google Drive", "GoogleDrive",
             "CloudStorage", "Box Sync", "pCloud", "Sync.com")
  hit <- cloud[vapply(cloud, function(k) grepl(k, p, fixed = TRUE), logical(1))]
  if (length(hit)) return(paste("cloud-sync folder:", hit[[1]]))
  NULL
}

check_path_length <- function(path) {
  if (nchar(path, type = "chars") >= MAX_PATH_CHARS)
    stop("SYS_PATH_TOO_LONG: the path is ", nchar(path), " characters; choose a shorter study folder or shorter names: ", path)
  invisible(path)
}

#' Create the study folder tree (C9). Returns the normalised study folder.
study_create <- function(study_dir) {
  prob <- study_path_problem(study_dir)
  if (!is.null(prob)) stop("SYS_STUDY_FOLDER_REFUSED: ", prob, ". Choose a folder on this computer that is not synchronised.")
  check_path_length(study_dir)
  dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)
  for (d in STUDY_FOLDERS) dir.create(file.path(study_dir, d), showWarnings = FALSE)
  dir.create(file.path(study_dir, ".engine"), showWarnings = FALSE)
  normalizePath(study_dir, winslash = "/", mustWork = TRUE)
}

pair_folder_name <- function(index, exposure_label, outcome_label) {
  sprintf("P%02d_%s__to__%s", as.integer(index), slugify(exposure_label, 30L), slugify(outcome_label, 30L))
}

# ---- atomic writes and the manifest (C9, FR-56) --------------------------------
#' Write bytes to a temporary file in the same folder, then rename (atomic on one volume).
atomic_write_raw <- function(path, bytes) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  con <- file(tmp, open = "wb"); writeBin(bytes, con); flush(con); close(con)
  if (!file.rename(tmp, path)) { unlink(tmp); stop("SYS_DISK_FULL: could not move ", basename(path), " into place") }
  invisible(path)
}
atomic_write_text <- function(path, text) atomic_write_raw(path, charToRaw(enc2utf8(text)))

manifest_path <- function(study_dir) file.path(study_dir, "MANIFEST.csv")

manifest_read <- function(study_dir) {
  p <- manifest_path(study_dir)
  if (!file.exists(p)) {
    m <- data.frame(relative_path = character(), sha256 = character(), bytes = numeric(), role = character(),
                    task_id = character(), created_by = character(), version = integer(), stringsAsFactors = FALSE)
    return(m)
  }
  m <- utils::read.csv(p, stringsAsFactors = FALSE, na.strings = character())
  m$version <- as.integer(m$version); m$bytes <- as.numeric(m$bytes)
  m
}

manifest_write <- function(study_dir, m) {
  m <- m[order(m$relative_path), MANIFEST_COLS, drop = FALSE]
  out <- utils::capture.output(utils::write.csv(m, row.names = FALSE, quote = TRUE, eol = "\n"))
  atomic_write_text(manifest_path(study_dir), paste0(paste(out, collapse = "\n"), "\n"))
}

#' Register files (paths relative to the study folder) in the manifest. A path that already exists with a
#' different hash gets the next version number; the manifest always shows the newest.
manifest_register <- function(study_dir, rel_paths, role, task_id = NA_character_, created_by = "app") {
  m <- manifest_read(study_dir)
  for (i in seq_along(rel_paths)) {
    rel <- rel_paths[[i]]
    abs <- file.path(study_dir, rel)
    if (nzchar(Sys.readlink(abs))) stop("symbolic links are not allowed in the study folder: ", rel)
    h <- sha256_file(abs); b <- file.size(abs)
    role_i <- if (length(role) == length(rel_paths)) role[[i]] else role[[1]]
    old <- m[m$relative_path == rel, , drop = FALSE]
    ver <- if (nrow(old)) max(old$version) else 0L
    if (nrow(old) && old$sha256[old$version == ver] == h) next          # unchanged
    m <- m[m$relative_path != rel, , drop = FALSE]                       # keep only the newest row per path
    m <- rbind(m, data.frame(relative_path = rel, sha256 = h, bytes = b, role = role_i, task_id = as.character(task_id),
                             created_by = created_by, version = ver + 1L, stringsAsFactors = FALSE))
  }
  manifest_write(study_dir, m)
  invisible(m)
}

#' Compare the manifest with the disk. Returns a data frame (relative_path, status) where status is
#' ok | modified | missing | untracked. Temporary files of interrupted writes are ignored.
manifest_check <- function(study_dir, tracked_folders = c("01_study_plan", "02_datasets", "03_pairs",
                                                          "04_across_pairs", "05_report_and_supplement", "06_reproducibility")) {
  m <- manifest_read(study_dir)
  st <- character(nrow(m))
  for (i in seq_len(nrow(m))) {
    p <- file.path(study_dir, m$relative_path[i])
    st[i] <- if (!file.exists(p)) "missing" else if (sha256_file(p) != m$sha256[i]) "modified" else "ok"
  }
  out <- data.frame(relative_path = m$relative_path, status = st, stringsAsFactors = FALSE)
  on_disk <- unlist(lapply(tracked_folders, function(d) {
    base <- file.path(study_dir, d)
    if (!dir.exists(base)) return(character())
    f <- list.files(base, recursive = TRUE, all.files = TRUE)
    f <- f[!grepl("\\.tmp-[0-9]+$", f)]
    file.path(d, f)
  }))
  extra <- setdiff(on_disk, m$relative_path)
  if (length(extra)) out <- rbind(out, data.frame(relative_path = extra, status = "untracked", stringsAsFactors = FALSE))
  out
}

# ---- the script task ---------------------------------------------------------
.default_flags <- function(adapter, params) {
  nm <- names(adapter$params)
  lapply(setNames(nm, nm), function(p) list(value = params[[p]], is_default = identical(params[[p]], adapter$params[[p]]$default)))
}

.rel <- function(study_dir, path) {
  s <- normalizePath(study_dir, winslash = "/", mustWork = TRUE)
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  substring(p, nchar(s) + 2L)
}

#' Move earlier outputs of this pair to archive/run_<datetime>/ (C9: re-running with changed settings).
archive_pair_outputs <- function(study_dir, pair_dir) {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
  dest <- file.path(study_dir, "archive", paste0("run_", stamp), basename(pair_dir))
  moved <- FALSE
  for (d in c("02_scripts", "03_console", "04_tables", "07_replay_check")) {
    src <- file.path(pair_dir, d)
    if (dir.exists(src) && length(list.files(src))) {
      dir.create(file.path(dest, d), recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(src, full.names = TRUE), file.path(dest, d), recursive = TRUE)
      unlink(list.files(src, full.names = TRUE), recursive = TRUE)
      moved <- TRUE
    }
  }
  # Round 41: the moved files are no longer in the pair folder, so their manifest rows go too. Otherwise a later run
  # sees a renumbered script (for example 05_loo.R that became 07_loo.R) as missing and stops with SCRIPT_USER_MODIFIED.
  if (moved) {
    m <- manifest_read(study_dir); rel <- .rel(study_dir, pair_dir)
    gone <- startsWith(m$relative_path, paste0(rel, "/")) & !file.exists(file.path(study_dir, m$relative_path))
    if (any(gone)) manifest_write(study_dir, m[!gone, , drop = FALSE])
  }
  if (moved) dest else NA_character_
}

#' Round 41: studies saved before this fix can list files the app itself moved to archive/ (renumbered or dropped
#' scripts). A listed file that is missing, but whose identical copy (same SHA-256) is in archive/ under this pair,
#' was archived by the app, so its row is dropped. A file the analyst deleted has no such copy and is still reported.
.manifest_heal_archived <- function(study_dir, pair_dir) {
  m <- manifest_read(study_dir); rel <- .rel(study_dir, pair_dir)
  miss <- which(startsWith(m$relative_path, paste0(rel, "/")) & !file.exists(file.path(study_dir, m$relative_path)))
  if (!length(miss)) return(invisible(FALSE))
  arch <- Sys.glob(file.path(study_dir, "archive", "run_*", basename(pair_dir), "*", "*"))
  if (!length(arch)) return(invisible(FALSE))
  ah <- vapply(arch, sha256_file, ""); drop <- miss[m$sha256[miss] %in% ah]
  if (length(drop)) manifest_write(study_dir, m[-drop, , drop = FALSE])
  invisible(length(drop) > 0)
}

#' Run the methods of one pair as generated scripts (C6.2, C8 task type `script`), read the results back,
#' commit them to the project database, run the replay check (C6.5) and update MANIFEST.csv.
#'
#' methods: list of list(adapter, input, params); input is a pair_iv input (C6.3).
#' on_modified: what to do when the user changed or deleted a tracked file of this pair:
#'   "stop" (default) refuses and names the files (code SCRIPT_USER_MODIFIED); "archive" moves them to archive/ first.
run_script_task <- function(con, study_dir, pair_index, exposure_label, outcome_label, methods,
                            replay = TRUE, on_modified = c("stop", "archive"), seed = NULL) {
  on_modified <- match.arg(on_modified)
  engine <- file.path(study_dir, ".engine")
  pair_dir <- file.path(study_dir, "03_pairs", pair_folder_name(pair_index, exposure_label, outcome_label))
  check_path_length(file.path(pair_dir, "02_scripts", "00_run_all.R"))
  for (d in PAIR_SUBFOLDERS) dir.create(file.path(pair_dir, d), recursive = TRUE, showWarnings = FALSE)
  pair_id <- basename(pair_dir)
  if (is.null(seed)) seed <- task_seed(pair_id)

  # 1. fingerprints: the same composition as v1.0.0 (inputs are the hashes of the analysis tables)
  tmp_in <- file.path(pair_dir, "01_data", ".fp-check"); dir.create(tmp_in, showWarnings = FALSE)
  plans <- lapply(methods, function(m) {
    a <- m$adapter; params <- resolve_params(a, m$params)
    csv <- file.path(tmp_in, "d.csv"); write_analysis_csv(m$input$data, csv)
    data_sha <- sha256_file(csv)
    fp <- task_fingerprint("method", list(analysis_data = data_sha, exposure = m$input$exposure_id,
                                          outcome = m$input$outcome_id, scale = m$input$outcome_scale), params, a)
    list(adapter = a, input = m$input, params = params, data_sha = data_sha, fp = fp,
         task_id = sprintf("%s:%s", pair_id, a$id))
  })
  unlink(tmp_in, recursive = TRUE)

  # 2. reuse (C8.3): every method done with the same fingerprint, files on disk untouched
  done_same <- vapply(plans, function(p) {
    t <- db_task(con, p$task_id); nrow(t) == 1L && t$state == "done" && identical(t$fingerprint, p$fp)
  }, logical(1))
  .manifest_heal_archived(study_dir, pair_dir)   # Round 41: rows left behind by archiving in older versions
  chk <- manifest_check(study_dir)
  pair_rel <- paste0("03_pairs/", pair_id, "/")
  touched <- chk[startsWith(chk$relative_path, pair_rel) & chk$status %in% c("modified", "missing"), , drop = FALSE]
  if (all(done_same) && !nrow(touched)) return(invisible(list(reused = TRUE, pair_dir = pair_dir, task_ids = vapply(plans, `[[`, "", "task_id"))))

  # 3. never overwrite the user's files silently
  if (nrow(touched)) {
    if (on_modified == "stop")
      stop("SCRIPT_USER_MODIFIED: these files were changed or deleted outside the app: ",
           paste(touched$relative_path, collapse = ", "), ". Choose 'archive' to keep your versions in archive/ and regenerate.")
    archive_pair_outputs(study_dir, pair_dir)
  } else if (!all(done_same)) {
    archive_pair_outputs(study_dir, pair_dir)       # settings changed: earlier outputs move to archive/
  }

  # 4. write the scripts and data, run in a fresh R, read the results back
  for (p in plans) { db_upsert_task(con, p$task_id, "script", p$adapter$id, p$adapter$adapter_version)
                     db_mark_running(con, p$task_id, p$fp, seed) }
  mats <- lapply(plans, function(p) list(adapter = p$adapter, input = p$input, params = p$params))
  materialise_pair(pair_dir, mats, seed)
  started <- Sys.time()
  run <- run_in_dir(pair_dir, run_pair_script(pair_dir))
  finished <- Sys.time()
  console_rel <- file.path("03_console", "run_all_console.txt")

  results <- list()
  # Round 70c: each method is judged on its own. The master script lists the scripts that stopped with an error; a method
  # whose script finished and wrote its tables is "ok" even when another method of the pair failed.
  ff <- file.path(pair_dir, "03_console", "failed_scripts.txt"); failed_files <- if (file.exists(ff)) trimws(readLines(ff, warn = FALSE)) else character()
  for (p in plans) {
    a <- p$adapter; codes <- character(); messages <- character()
    script_name <- sprintf("%02d_%s.R", match(a$id, vapply(plans, function(q) q$adapter$id, "")), a$id)
    out <- if (!(script_name %in% failed_files)) tryCatch(read_method_results(pair_dir, a$id), error = function(e) NULL) else NULL
    if (is.null(out)) { codes <- "SCRIPT_RUN_FAILED"
      own <- grep(paste0("FAILED ", script_name, ": "), run$console, fixed = TRUE, value = TRUE)
      messages <- if (length(own)) paste(own, collapse = "\n") else paste(utils::tail(run$console, 8), collapse = "\n") }
    status <- if (is.null(out)) "failed" else "ok"
    script_file <- file.path(pair_dir, "02_scripts", script_name)
    res <- list(contract_version = CONTRACT_VERSION, task_id = p$task_id, fingerprint = p$fp, method_id = a$id,
                adapter_version = a$adapter_version, status = status, codes = unique(codes), messages = messages,
                estimates = if (is.null(out)) NULL else out$estimates,
                diagnostics = if (is.null(out)) NULL else out$diagnostics, variant_level_blob = NULL,
                script_sha256 = sha256_file(script_file), analysis_data_sha256 = p$data_sha,
                params_effective = .default_flags(a, p$params),
                customised = !all(vapply(.default_flags(a, p$params), function(z) z$is_default, logical(1))),
                user_modified = FALSE,
                provenance = list(engine_version = "0.0.0-contracts", r_version = paste(R.version$major, R.version$minor, sep = "."),
                  os = os_family_id(), packages = a$packages, seed = seed, rng_kind = "Mersenne-Twister/Inversion/Rejection",
                  blas = basename(extSoftVersion()[["BLAS"]]), threads = 1L,
                  started_at = format(started, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
                  finished_at = format(finished, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
                  elapsed_s = as.double(difftime(finished, started, units = "secs"))))
    if (isTRUE(a$diagnostic_only)) res$diagnostic_only <- TRUE   # contracts v1.1 (additive): only set for diagnostic-only methods, so other results are unchanged
    json <- result_to_json(res)
    blob <- blob_put_text(engine, json)
    db_commit_result(con, p$task_id, p$fp, status, blob, json)
    DBI::dbExecute(con, "UPDATE results SET script_sha256=?, analysis_data_sha256=?, params_effective_json=?, customised=?, user_modified=0 WHERE task_id=?",
                   params = list(res$script_sha256, res$analysis_data_sha256, canonical_json(res$params_effective), as.integer(res$customised), p$task_id))
    DBI::dbExecute(con, "INSERT OR IGNORE INTO script_versions(task_id, version, script_sha256, created_by, created_at) VALUES (?, 1, ?, 'app', ?)",
                   params = list(p$task_id, res$script_sha256, .now()))
    results[[a$id]] <- res
  }

  # 5. replay check (C6.5): automatic for final runs, optional for exploratory ones (D3)
  rp <- NULL
  ok_ids <- names(Filter(function(r) identical(r$status, "ok"), results))
  if (replay && length(ok_ids)) {   # Round 70c: the methods that finished are replayed and compared
    stored <- lapply(ok_ids, function(id) read_method_results(pair_dir, id)); names(stored) <- ok_ids
    rp <- replay_check(pair_dir, names(stored), stored)
    atomic_write_text(file.path(pair_dir, "07_replay_check", "replay_check.json"), paste0(canonical_json(rp), "\n"))
    for (p in plans) if (p$adapter$id %in% ok_ids) DBI::dbExecute(con, "UPDATE results SET replay_json=? WHERE task_id=?", params = list(canonical_json(rp), p$task_id))
  } else if (!replay) {
    rp <- list(status = "not_run")
    for (p in plans) DBI::dbExecute(con, "UPDATE results SET replay_json=? WHERE task_id=?", params = list(canonical_json(rp), p$task_id))
  }

  # 6. manifest: every delivered file is listed with its hash
  files <- list.files(pair_dir, recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("\\.tmp-[0-9]+$", files)]
  roles <- vapply(files, function(f) {
    d <- basename(dirname(f))
    switch(d, `01_data` = "analysis_data", `02_scripts` = "script", `03_console` = "console", `04_tables` = "result_table",
           `07_replay_check` = "replay_check", "file")
  }, "")
  manifest_register(study_dir, vapply(files, function(f) .rel(study_dir, f), ""), unname(roles))
  DBI::dbWriteTable(con, "manifest_tmp", transform(manifest_read(study_dir), created_at = .now()), overwrite = TRUE, temporary = TRUE)
  DBI::dbExecute(con, "INSERT OR REPLACE INTO manifest(relative_path, version, sha256, bytes, role, task_id, created_by, created_at)
                       SELECT relative_path, version, sha256, CAST(bytes AS INTEGER), role, NULLIF(task_id, ''), created_by, created_at FROM manifest_tmp")
  DBI::dbExecute(con, "DROP TABLE IF EXISTS manifest_tmp")
  invisible(list(reused = FALSE, pair_dir = pair_dir, results = results, replay = rp, run_status = run$status,
                 task_ids = vapply(plans, `[[`, "", "task_id")))
}
