# MR GUI contracts v1 — blob store, project DB, task runner, resume (C8, C9)
# Requires DBI + RSQLite. Single writer: only the main process calls db_* functions.

# ---- blob store (C8.7 step 1) ---------------------------------------------
blob_path <- function(root, hash) file.path(root, "blobs", substr(hash, 1, 2), hash)

blob_put_text <- function(root, text) {
  tmp_dir <- file.path(root, "blobs", "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = tmp_dir)
  con <- file(tmp, open = "wb"); writeBin(charToRaw(enc2utf8(text)), con); close(con)
  hash <- sha256_file(tmp)
  dest <- blob_path(root, hash)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest)) {          # identical content already stored
    unlink(tmp)
  } else if (!file.rename(tmp, dest)) {
    unlink(tmp); stop("SYS_DISK_FULL: could not move blob into place")
  }
  hash
}
blob_get_text <- function(root, hash) {
  p <- blob_path(root, hash)
  if (sha256_file(p) != hash) stop("blob ", hash, " failed its integrity check")
  rawToChar(readBin(p, "raw", file.size(p)))
}

# ---- project database (C9) --------------------------------------------------
db_open <- function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(root, "project.sqlite"))
  DBI::dbExecute(con, "PRAGMA journal_mode=WAL")
  DBI::dbExecute(con, "PRAGMA synchronous=FULL")
  DBI::dbExecute(con, "PRAGMA foreign_keys=ON")
  if (!DBI::dbExistsTable(con, "meta")) {
    ddl <- paste(readLines(.schema_path(), warn = FALSE), collapse = "\n")
    for (stmt in strsplit(ddl, ";\\s*\n")[[1]]) if (nzchar(trimws(gsub("--[^\n]*", "", stmt)))) DBI::dbExecute(con, stmt)
  }
  con
}
.schema_path <- function() getOption("mrgui.schema", "contracts/sql/schema_v1.sql")
.now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

db_upsert_task <- function(con, task_id, task_type, method_id = NA, adapter_version = NA) {
  DBI::dbExecute(con,
    "INSERT INTO tasks(task_id, task_type, state, method_id, adapter_version, updated_at)
     VALUES (?, ?, 'pending', ?, ?, ?) ON CONFLICT(task_id) DO NOTHING",
    params = list(task_id, task_type, method_id, adapter_version, .now()))
}
db_task <- function(con, task_id) DBI::dbGetQuery(con, "SELECT * FROM tasks WHERE task_id = ?", params = list(task_id))

db_mark_running <- function(con, task_id, fingerprint, seed) {
  DBI::dbExecute(con,
    "UPDATE tasks SET state='running', fingerprint=?, seed=?, attempt=attempt+1, started_at=?, updated_at=? WHERE task_id=?",
    params = list(fingerprint, seed, .now(), .now(), task_id))
}

#' C8.7 step 3: one transaction makes the result visible.
db_commit_result <- function(con, task_id, fingerprint, status, result_blob, result_json) {
  state <- if (status %in% c("ok", "warn")) "done" else status
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "INSERT OR IGNORE INTO blobs(sha256, size, media_type, created_at) VALUES (?, ?, 'application/json', ?)",
                   params = list(result_blob, nchar(result_json, type = "bytes"), .now()))
    DBI::dbExecute(con, "INSERT OR REPLACE INTO results(task_id, fingerprint, status, result_blob, created_at) VALUES (?, ?, ?, ?, ?)",
                   params = list(task_id, fingerprint, status, result_blob, .now()))
    DBI::dbExecute(con, "UPDATE tasks SET state=?, finished_at=?, updated_at=? WHERE task_id=?",
                   params = list(state, .now(), .now(), task_id))
  })
}

#' C8.2/C8.8: anything left 'running' by a dead process goes back to 'pending'.
db_recover <- function(con) {
  DBI::dbExecute(con, "UPDATE tasks SET state='pending', updated_at=? WHERE state='running'", params = list(.now()))
}

# ---- running one method task (C6.3, C6.4, C7, C8.4) -------------------------
run_method_task <- function(con, root, adapter, input, params, task_id, input_blobs, time_budget_s = 3600) {
  params <- resolve_params(adapter, params)
  fp <- task_fingerprint("method", input_blobs, params, adapter)
  db_upsert_task(con, task_id, "method", adapter$id, adapter$adapter_version)
  prev <- db_task(con, task_id)
  if (nrow(prev) && prev$state == "done" && identical(prev$fingerprint, fp)) {
    return(invisible(list(reused = TRUE, fingerprint = fp)))     # C8.3 reuse
  }
  seed <- task_seed(task_id)
  db_mark_running(con, task_id, fp, seed)

  scratch <- file.path(root, "scratch", substr(fp, 1, 16))
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  ctx <- list(seed = seed, scratch_dir = scratch, time_budget_s = time_budget_s,
              progress = function(fraction, message = "") invisible(NULL),
              is_cancelled = function() FALSE,
              checkpoint_save = function(obj) invisible(NULL), checkpoint_load = function() NULL,
              log = function(level, message) invisible(NULL))

  started <- Sys.time()
  codes <- character(); messages <- character()
  v <- adapter$check(input, params)
  out <- NULL
  if (v$status == "block") {
    status <- "blocked"; codes <- v$codes; messages <- v$messages
  } else {
    codes <- v$codes; messages <- v$messages
    old_kind <- RNGkind(); on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
    out <- tryCatch(
      withCallingHandlers(adapter$run(input, params, ctx), warning = function(w) {
        codes <<- c(codes, "RUN_PACKAGE_WARNING"); messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) { codes <<- c(codes, "RUN_PACKAGE_ERROR"); messages <<- c(messages, conditionMessage(e)); NULL })
    status <- if (is.null(out)) "failed" else if (length(codes)) "warn" else "ok"
  }
  finished <- Sys.time()

  res <- list(
    contract_version = CONTRACT_VERSION, task_id = task_id, fingerprint = fp,
    method_id = adapter$id, adapter_version = adapter$adapter_version, status = status,
    codes = unique(codes), messages = messages,
    estimates = if (is.null(out)) NULL else out$estimates,
    diagnostics = if (is.null(out)) NULL else out$diagnostics,
    variant_level_blob = NULL,
    provenance = list(
      engine_version = "0.0.0-contracts", r_version = paste(R.version$major, R.version$minor, sep = "."),
      os = os_family_id(), packages = adapter$packages, seed = seed,
      rng_kind = paste(RNGkind(), collapse = "/"), blas = basename(extSoftVersion()[["BLAS"]]),
      threads = 1L, started_at = format(started, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
      finished_at = format(finished, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
      elapsed_s = as.double(difftime(finished, started, units = "secs"))))
  json <- result_to_json(res)
  blob <- blob_put_text(root, json)                         # C8.7 step 1
  db_commit_result(con, task_id, fp, status, blob, json)   # C8.7 step 3
  unlink(scratch, recursive = TRUE)
  invisible(list(reused = FALSE, fingerprint = fp, result = res, blob = blob))
}
