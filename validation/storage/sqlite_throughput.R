# Round 36: is SQLite with synchronous=FULL fast enough for MR Studio's study database (contracts_v1.md, open item 2)?
# Uses the runner's own functions (db_open, db_upsert_task, db_mark_running, db_commit_result) to write the same
# sequence of commits a method task makes, for TASKS tasks, once with synchronous=FULL (what the app uses) and once
# with synchronous=NORMAL for comparison. Run from the contracts folder:
#   Rscript --vanilla ../validation/storage/sqlite_throughput.R [tasks, default 500]
suppressMessages({ library(DBI); library(RSQLite) })
for (f in c("r/canonical.R", "r/runner.R")) source(f); options(mrgui.schema = "sql/schema_v1.sql")
tasks <- suppressWarnings(as.integer(commandArgs(TRUE)[1])); if (!is.finite(tasks)) tasks <- 500L
bench <- function(sync) {
  root <- tempfile("mrgui-sqlite-"); con <- db_open(root); on.exit({ DBI::dbDisconnect(con); unlink(root, recursive = TRUE) })
  DBI::dbExecute(con, paste0("PRAGMA synchronous=", sync))
  res <- paste(rep("{\"estimate\":0.341,\"se\":0.0517}", 40), collapse = "")   # about 1.3 kB, the size of a result envelope
  t <- system.time(for (i in seq_len(tasks)) {
    id <- sprintf("task_%05d", i); db_upsert_task(con, id, "method", "ivw_mre", "1.1.0")
    db_mark_running(con, id, sprintf("%064d", i), i)
    db_commit_result(con, id, sprintf("%064d", i), "ok", blob_put_text(root, paste0(res, i)), paste0(res, i)) })[["elapsed"]]
  data.frame(synchronous = sync, tasks = tasks, seconds = round(t, 2), ms_per_task = round(1000 * t / tasks, 2), tasks_per_second = round(tasks / t))
}
r <- rbind(bench("FULL"), bench("NORMAL"))
r$system <- paste(Sys.info()[["sysname"]], Sys.info()[["machine"]]); r$sqlite <- RSQLite::rsqliteVersion()[["library"]]
print(r, row.names = FALSE)
