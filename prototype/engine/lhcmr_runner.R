# ---- Round 67: adaptive, resumable LHC-MR runner --------------------------------------------------------------------------------
# lhc_mr() hands its 50 starting-point fits and its 200 jackknife fits to parallel::mclapply with a fixed number of workers, so a
# long run cannot adapt to the memory free on the computer, and an interrupted run starts again from nothing. This runner evaluates
# the statements of the installed, pinned lhc_mr() itself, in order, and replaces only its two mclapply calls by a scheduler that
#   - runs the package's own unit function (the same optim() call) in forked workers, a batch at a time;
#   - saves every finished unit (keyed by a hash of its input), so a run that was stopped resumes where it was;
#   - chooses the number of workers before each batch from the memory free now (at most `max_workers`, growing by one per batch),
#     stops the batch and retries with fewer workers when free memory falls below its floor, and pauses (checking every minute)
#     when even one worker does not fit, instead of failing.
# The units are deterministic (Nelder-Mead from fixed starting points), so the results do not depend on how they were scheduled;
# the bimodality check of the jackknife summary (normal mixtures with random starts) is seeded for the same reason. Validated: on the
# same input, the fit tables and SummarisedResults equal lhc_mr()'s exactly, and a resumed run reads its checkpoints (Round 67 tests).
# The source below is evaluated by the app and written verbatim into the study script, so the script replays without the app.
LHCMR_RUNNER_SRC <- r"--(lhcmr_schedule <- function(X, FUN, stage = 1L, ck = "units", max_workers = 2L, per_worker_gb = 1, safe_share = 0.5, floor_share = 0.15,
                           status = NULL, poll = 5, pause = 60, mem = function() list(total_gb = NA_real_, available_gb = NA_real_), state = new.env()) {
  st <- function(s, msg) if (!is.null(status)) writeLines(jsonlite::toJSON(list(state = s, msg = msg, pid = Sys.getpid(), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), status)
  if (is.null(state$k)) { state$k <- 1L; state$cool <- 0L }
  lab <- c("starting points", "jackknife blocks")[min(stage, 2L)]; n <- length(X); dir.create(ck, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(ck, sprintf("stage%d_%03d.rds", stage, seq_len(n))); key <- vapply(X, function(x) digest::digest(x), "")
  have <- function(i) file.exists(f[i]) && identical(tryCatch(readRDS(f[i])$key, error = function(e) ""), key[i])
  fork <- .Platform$OS.type != "windows"
  repeat {
    todo <- which(!vapply(seq_len(n), have, TRUE)); if (!length(todo)) break
    mi <- mem(); fits <- if (is.finite(mi$available_gb)) floor(safe_share * mi$available_gb / per_worker_gb) else 1L
    target <- max(0L, min(as.integer(max_workers), as.integer(fits), if (fork) length(todo) else 1L))
    state$k <- if (target > state$k && state$cool == 0L) state$k + 1L else min(state$k, target); state$cool <- max(0L, state$cool - 1L); k <- state$k
    if (k < 1L) { st("paused", sprintf("Paused: not enough free memory for one worker (%.1f GB free; about %.1f GB is needed). It continues by itself when memory is free again.", mi$available_gb, per_worker_gb / safe_share)); Sys.sleep(pause); next }
    st("running", sprintf("LHC-MR %s: %d of %d done, %d worker%s.", lab, n - length(todo), n, k, if (k == 1L) "" else "s"))
    batch <- utils::head(todo, k)
    if (!fork) { for (i in batch) saveRDS(list(key = key[i], value = FUN(X[[i]])), f[i]); next }
    left <- lapply(batch, function(i) parallel::mcparallel({ saveRDS(list(key = key[i], value = FUN(X[[i]])), f[i]); TRUE }, silent = TRUE))
    while (length(left)) {
      done <- parallel::mccollect(left, wait = FALSE, timeout = poll)
      if (length(done)) left <- Filter(function(j) !(as.character(j$pid) %in% names(done)), left)
      if (!length(left)) break
      mi <- mem()
      if (is.finite(mi$available_gb) && mi$available_gb < max(1, floor_share * mi$total_gb)) {   # memory low: stop this batch, fewer workers next
        for (j in left) tools::pskill(j$pid); invisible(suppressWarnings(parallel::mccollect(left, wait = TRUE)))
        state$k <- max(0L, k %/% 2L); state$cool <- 3L
        st("running", sprintf("Free memory fell to %.1f GB: continuing with %d worker%s.", mi$available_gb, state$k, if (state$k == 1L) "" else "s")); break
      }
    }
  }
  lapply(seq_len(n), function(i) readRDS(f[i])$value)
}
lhcmr_adaptive <- function(SP_list, trait.names = c("X", "Y"), nBlock = 200L, ck = "units", max_workers = 2L, per_worker_gb = 1,
                           safe_share = 0.5, floor_share = 0.15, status = NULL, poll = 5, pause = 60, mem = function() list(total_gb = NA_real_, available_gb = NA_real_), seed = NULL) {
  ns <- asNamespace("lhcMR"); stmts <- as.list(body(base::get("lhc_mr", envir = ns)))[-1]   # base::get: a function named get in the session must not shadow it
  is_if <- function(s, cond) is.call(s) && identical(s[[1]], as.name("if")) && identical(paste(deparse(s[[2]]), collapse = ""), cond)
  is_mcl <- function(s) is.call(s) && (identical(s[[1]], as.name("<-")) || identical(s[[1]], as.name("="))) && is.call(s[[3]]) && identical(paste(deparse(s[[3]][[1]]), collapse = ""), "parallel::mclapply")
  sched <- new.env(); stage <- 0L
  adaptive_lapply <- function(X, FUN) { stage <<- stage + 1L
    lhcmr_schedule(X, FUN, stage, ck, max_workers, per_worker_gb, safe_share, floor_share, status, poll, pause, mem, sched) }
  e <- new.env(parent = ns)
  e$SP_list <- SP_list; e$trait.names <- trait.names; e$partition <- NA; e$account <- NA; e$param <- "comp"; e$paral_method <- "lapply"; e$nCores <- 1L; e$nBlock <- nBlock; e$M <- 1e7
  e$adaptive_lapply <- adaptive_lapply
  run <- function(s) if (is_mcl(s)) { s[[3]] <- as.call(list(as.name("adaptive_lapply"), s[[3]][[2]], s[[3]][[3]])); eval(s, e) } else eval(s, e)
  k2 <- which(vapply(stmts, is_if, TRUE, cond = "nStep == 2"))[1]; k1 <- which(vapply(stmts, is_if, TRUE, cond = "nStep == 1"))[1]
  for (s in stmts[seq_len(k2 - 1L)]) run(s)
  if (!identical(e$nStep, 2)) stop("The adaptive runner covers the two-step LHC-MR (nStep = 2), the setting MR Studio uses.")
  blk <- as.list(stmts[[k2]][[3]])[-1]; kl <- which(vapply(blk, is_if, TRUE, cond = "paral_method == \"lapply\""))[1]
  for (s in blk[seq_len(kl - 1L)]) if (!is_if(s, "paral_method == \"rslurm\"")) run(s)
  for (s in as.list(blk[[kl]][[3]])[-1]) run(s)
  if (!is.null(seed)) set.seed(seed)   # the jackknife summary fits normal mixtures with random starts (mixtools): seeded, so it does not depend on the schedule
  for (s in stmts[-seq_len(k1)]) run(s)
  if (!is.null(status)) writeLines(jsonlite::toJSON(list(state = "done", msg = "", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), status)
  invisible(TRUE)
})--"
