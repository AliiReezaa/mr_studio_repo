# Round 83: long analyses run in a separate R process, so the app stays usable while they run.
#
# Before, a study's run (every pair, every kind of study), the sensitivity, tie and outlier reanalyses, mediation, bidirectional
# MR-cML, LD score regression and the publication package ran inside the engine. The engine handles one message at a time, so for
# the minutes a pair took it answered nothing: the Studies screen did not list or open other studies, and the running study could not
# be opened again (reported by the user on 2026-11-21). Now each such action is started as a job: a background R process loads the
# engine, opens the study itself and runs the same code; the messages it would have sent to the screen are written, one JSON line
# each, to the job's events file, which the engine reads every half second and passes on, tagged with the study they belong to. The
# engine is free meanwhile. A job belongs to its study folder, so opening another study during a run no longer moves the run's
# remaining pairs into the study that was opened (the queue used whichever study was open).
#
# Two processes now share a study's database: each connection waits up to 15 s for the other's write (busy timeout), and opening a
# study whose job is running does not reset that job's running task to pending (the crash recovery of contracts C8 is for tasks
# left running by a process that is gone).

`%||%` <- function(a, b) if (is.null(a)) b else a

BG_ACTIONS <- list()          # name -> function(h, a, ev); filled in app.R next to the screen messages they answer
BG_POLL_S <- 0.5

.bg_jobs_dir <- function(study_dir) file.path(study_dir, ".engine", "jobs")
.bg_alive <- function(pid) {
  if (is.null(pid)) return(TRUE)   # starting
  pid <- as.integer(pid)
  if (.Platform$OS.type == "windows") any(grepl(paste0(" ", pid, " "), suppressWarnings(system2("tasklist", c("/FI", shQuote(paste("PID eq", pid))), stdout = TRUE, stderr = FALSE))))
  else identical(suppressWarnings(system2("kill", c("-0", pid), stdout = FALSE, stderr = FALSE)), 0L)
}
.bg_status <- function(job) { f <- file.path(job, "status.json"); if (!file.exists(f)) return(list()); tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) list()) }
.bg_write_status <- function(job, ...) { s <- .bg_status(job); s[names(list(...))] <- list(...); f <- file.path(job, "status.json"); tmp <- paste0(f, ".tmp")
  writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null"), tmp); file.rename(tmp, f); invisible(s) }

#' The jobs of a study that are still running (their process is alive).
bg_running <- function(study_dir) {
  d <- .bg_jobs_dir(study_dir); if (!dir.exists(d)) return(character())
  j <- list.dirs(d, recursive = FALSE)
  j[vapply(j, function(x) { s <- .bg_status(x); identical(s$state, "running") && .bg_alive(s$pid) }, TRUE)]
}

#' Open a study for a job, or for the screen while a job of it runs: like studio_open, without the crash recovery.
studio_attach <- function(dir) {
  study <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  list(dir = study, con = db_open(file.path(study, ".engine")))
}

#' Start action `t` with arguments `a` for the study in `study_dir`. Returns the job folder.
bg_start <- function(study_dir, t, a, act = "", engine = getOption("mrgui.engine_dir", "")) {
  if (is.null(BG_ACTIONS[[t]])) stop("Not a background action: ", t)
  job <- file.path(.bg_jobs_dir(study_dir), sprintf("%s_%s_%04d", t, format(Sys.time(), "%Y%m%d-%H%M%S"), sample.int(9999, 1)))
  dir.create(job, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(t = t, a = a, act = act, study = normalizePath(study_dir, winslash = "/"), studies_root = studies_root()), file.path(job, "args.rds"))
  file.create(file.path(job, "events.jsonl"))
  .bg_write_status(job, state = "running", t = t, started = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  code <- sprintf("source(file.path(%s, 'bg_core.R')); bg_job_run(%s, %s)", deparse(engine), deparse(engine), deparse(normalizePath(job, winslash = "/")))
  rscript_background(c("--vanilla", "-e", shQuote(code)), log = file.path(job, "job.log"))
  job
}

#' The job itself, in the background process: loads the engine as main.R does, runs the action, writes its messages.
bg_job_run <- function(engine, job) {
  .bg_write_status(job, pid = Sys.getpid())
  x <- readRDS(file.path(job, "args.rds")); Sys.setenv(MRGUI_STUDIES_DIR = x$studies_root)
  ev <- function(t, ...) { line <- as.character(shiny:::toJSON(c(list(t = t), list(...))))   # as Shiny would send it
    cat(gsub("\n", " ", line), "\n", sep = "", file = file.path(job, "events.jsonl"), append = TRUE) }
  ok <- tryCatch({
    args_dir <- engine
    mf <- readLines(file.path(engine, "main.R"), warn = FALSE)
    files <- regmatches(mf, regexpr("(?<=source\\(file\\.path\\(args_dir, \")[^\"]+", mf, perl = TRUE))
    for (f in files) { source(file.path(engine, f), local = globalenv()); if (identical(f, "studio_core.R")) load_contracts(find_contracts_dir(engine)) }
    options(mrgui.engine_dir = engine, mrstudio.act = x$act %||% "",
            mrstudio.say = function(msg, why = "") ev("live", msg = msg, why = why, act = getOption("mrstudio.act", "")),
            mrstudio.script_progress = function(dir, console_name) live_scripts(dir, console_name))
    h <- studio_attach(x$study); on.exit(try(DBI::dbDisconnect(h$con), silent = TRUE), add = TRUE)
    get("BG_ACTIONS", globalenv())[[x$t]](h, x$a, ev)
    TRUE
  }, error = function(e) { ev("error", msg = conditionMessage(e)); FALSE })
  if (nzchar(x$act %||% "")) ev("live_end", act = x$act)
  .bg_write_status(job, state = if (isTRUE(ok)) "done" else "error", finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  invisible(ok)
}

#' New messages of a job since byte `from`: list(lines, at, done). Only complete lines are taken; a job whose process is gone without
#' saying so is done (with an error message for the screen).
bg_poll <- function(job, from = 0) {
  f <- file.path(job, "events.jsonl"); n <- if (file.exists(f)) file.size(f) else 0
  lines <- character(); at <- from
  if (n > from) { con <- file(f, "rb"); on.exit(close(con)); seek(con, from); raw <- readBin(con, "raw", n - from)
    cut <- max(c(0L, which(raw == as.raw(10L)))); if (cut > 0) { lines <- strsplit(rawToChar(raw[seq_len(cut)]), "\n", fixed = TRUE)[[1]]; at <- from + cut } }
  s <- .bg_status(job); gone <- identical(s$state, "running") && !is.null(s$pid) && !.bg_alive(s$pid)
  if (gone) lines <- c(lines, as.character(jsonlite::toJSON(list(t = "error", msg = "The analysis stopped unexpectedly (its process ended). What was saved is kept; run it again to finish it."), auto_unbox = TRUE)))
  list(lines = lines[nzchar(lines)], at = at, done = s$state %in% c("done", "error") || gone, ok = identical(s$state, "done"))
}
