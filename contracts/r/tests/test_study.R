# Tests for contracts v1.1 (draft): study folder, manifest, script task, archive, replay, schema v2.
# Run from the contracts/ folder:   Rscript r/tests/test_study.R <output_dir>
# Needs DBI, RSQLite, digest and TwoSampleMR at the pinned commit.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- tempfile("run12-", tmpdir = if (length(args)) args[1] else tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c("r/canonical.R", "r/adapter.R", "r/security.R", "r/runner.R", "r/adapters/reference_adapters.R",
            "r/script.R", "r/adapters/reference_adapters_v11.R", "r/study.R")) source(f)
options(mrgui.schema = "sql/schema_v2.sql")

n_pass <- 0L; n_fail <- 0L
expect <- function(desc, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("PASS ", desc, "\n") } else { n_fail <<- n_fail + 1L; cat("FAIL ", desc, "\n") }
}
near <- function(a, b, tol = 1e-12) isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))
err_of <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

set.seed(20261029); k <- 30
bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE); sx <- runif(k, 0.005, 0.012); sy <- runif(k, 0.008, 0.02)
by <- 0.3 * bx + rnorm(k, 0.002, 0.01) + rnorm(k, 0, sy)
dat <- data.frame(variant_id = sprintf("1:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G", beta_x = bx, se_x = sx, eaf_x = 0.3,
                  beta_y = by, se_y = sy, eaf_y = 0.3, p_x = 2 * pnorm(-abs(bx / sx)), p_y = 2 * pnorm(-abs(by / sy)),
                  n_x = 1e5, n_y = 1e5, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
dat$mr_keep[c(4, 17)] <- FALSE
input <- list(kind = "pair_iv", data = dat, exposure_id = "src:opengwas:ieu-a-300", outcome_id = "src:opengwas:ieu-a-7", outcome_scale = "log_or")
mk <- function(ci = 0.95) list(list(adapter = ADAPTER_IVW_MRE_V11, input = input, params = list(ci_level = ci)),
                               list(adapter = ADAPTER_EGGER_V11, input = input, params = list()))

# ---- names, paths, refusals --------------------------------------------------------------------
expect("slug: ASCII, at most 40 characters, never empty",
       identical(slugify("Omega-3 fatty acids / \u03a9"), "Omega-3-fatty-acids") && nchar(slugify(strrep("ab ", 40))) <= 40 && identical(slugify("\u4e2d\u6587"), "x"))
expect("pair folder name follows P<nn>_<exposure>__to__<outcome>", identical(pair_folder_name(1, "Omega-3", "Triglycerides"), "P01_Omega-3__to__Triglycerides"))
expect("cloud-sync folders are refused", !is.null(study_path_problem("/Users/a/Library/Mobile Documents/com~apple~CloudDocs/x")) &&
       !is.null(study_path_problem("C:/Users/a/OneDrive/Study")) && !is.null(study_path_problem("/Users/a/Dropbox/Study")))
expect("network drives are refused", !is.null(study_path_problem("\\\\server\\share\\Study")) && is.null(study_path_problem("/Users/a/Documents/MR Studio")))
expect("a path of 240 characters or more is refused with SYS_PATH_TOO_LONG",
       grepl("SYS_PATH_TOO_LONG", err_of(check_path_length(strrep("a", 240)))) && is.na(err_of(check_path_length(strrep("a", 239)))))
expect("creating a study in a cloud folder fails with SYS_STUDY_FOLDER_REFUSED",
       grepl("SYS_STUDY_FOLDER_REFUSED", err_of(study_create(file.path(out_dir, "OneDrive", "S")))))

# ---- schema v2 ---------------------------------------------------------------------------------
study <- study_create(file.path(out_dir, "My Study"))        # a folder name with a space on purpose
expect("study tree has the numbered folders and a hidden .engine", all(dir.exists(file.path(study, c(STUDY_FOLDERS, ".engine")))))
con <- db_open(file.path(study, ".engine"))
expect("schema v2 is applied", identical(DBI::dbGetQuery(con, "SELECT value FROM meta WHERE key='schema_version'")$value, "2"))
expect("tasks accept the new task types",
       DBI::dbExecute(con, "INSERT INTO tasks(task_id, task_type, state, updated_at) VALUES ('t_u','user_script','pending','x')") == 1L)
invisible(DBI::dbExecute(con, "DELETE FROM tasks WHERE task_id = 't_u'"))

# ---- first run of a pair -----------------------------------------------------------------------
r1 <- run_script_task(con, study, 1, "Omega-3", "Triglycerides", mk(0.9))
pair <- r1$pair_dir
expect("first run is not a reuse and the scripts ran (exit status 0)", !r1$reused && r1$run_status == 0L)
expect("both methods finished with status ok", all(vapply(r1$results, function(x) x$status, "") == "ok"))
t <- DBI::dbGetQuery(con, "SELECT task_id, task_type, state FROM tasks ORDER BY task_id")
expect("two tasks of type script are done", nrow(t) == 2L && all(t$task_type == "script") && all(t$state == "done"))
expect("all pair subfolders exist", all(dir.exists(file.path(pair, PAIR_SUBFOLDERS))))
expect("delivered files: data CSV, scripts, console, result tables, replay record",
       all(file.exists(file.path(pair, c("01_data/analysis_data_ivw_mre.csv", "02_scripts/00_run_all.R", "02_scripts/01_ivw_mre.R",
         "02_scripts/02_egger.R", "03_console/run_all_console.txt", "04_tables/ivw_mre_estimates.csv", "04_tables/egger_estimates.csv",
         "07_replay_check/replay_check.json")))))
ctx0 <- list(seed = 1L)
o_ivw <- ADAPTER_IVW_MRE$run(input, resolve_params(ADAPTER_IVW_MRE, list(ci_level = 0.9)), ctx0)
o_eg <- ADAPTER_EGGER$run(input, resolve_params(ADAPTER_EGGER, list()), ctx0)
cmp <- function(a, b) { num <- names(a)[vapply(a, is.numeric, logical(1))]; all(vapply(num, function(n) near(a[[n]], b[[n]]), logical(1))) }
expect("results equal the in-process oracle within 1e-12",
       cmp(r1$results$ivw_mre$estimates, o_ivw$estimates) && cmp(r1$results$egger$estimates, o_eg$estimates))
expect("replay check ran automatically and is identical", identical(r1$replay$status, "identical") && r1$replay$n_values_compared == 27L)
rr <- DBI::dbGetQuery(con, "SELECT * FROM results ORDER BY task_id")
expect("C7 fields are stored: script hash, data hash, effective params, customised, replay",
       nrow(rr) == 2L && all(nchar(rr$script_sha256) == 64) && all(nchar(rr$analysis_data_sha256) == 64) &&
       grepl("is_default", rr$params_effective_json[1]) && all(grepl("identical", rr$replay_json)) && any(rr$customised == 1L))
expect("the stored script hash equals the file on disk", identical(rr$script_sha256[rr$task_id == basename(pair) | grepl("ivw", rr$task_id)][1], sha256_file(file.path(pair, "02_scripts", "01_ivw_mre.R"))))

# ---- manifest ----------------------------------------------------------------------------------
chk <- manifest_check(study)
expect("MANIFEST.csv exists and lists every delivered file (none missing, modified or untracked)", file.exists(file.path(study, "MANIFEST.csv")) && nrow(chk) == 11L && all(chk$status == "ok"))
mf <- manifest_read(study)
expect("manifest has the agreed columns and 64-character hashes", identical(names(mf), MANIFEST_COLS) && all(nchar(mf$sha256) == 64))
expect("the project database mirrors the manifest", DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM manifest")$n == nrow(mf))
files_all <- list.files(study, recursive = TRUE, full.names = TRUE, all.files = TRUE)
expect("no symbolic links and no temporary files are left in the study folder",
       !any(nzchar(Sys.readlink(files_all))) && !any(grepl("\\.tmp-[0-9]+$", files_all)))
stale <- file.path(study, "03_pairs", "stale.tmp-999"); writeLines("partial", stale)
expect("a temporary file from an interrupted write is ignored by the check", !("03_pairs/stale.tmp-999" %in% manifest_check(study)$relative_path))
unlink(stale)

# ---- reuse -------------------------------------------------------------------------------------
hash_before <- sha256_file(file.path(pair, "02_scripts", "01_ivw_mre.R"))
r2 <- run_script_task(con, study, 1, "Omega-3", "Triglycerides", mk(0.9))
expect("same inputs and settings are reused and nothing is rewritten", isTRUE(r2$reused) && identical(sha256_file(file.path(pair, "02_scripts", "01_ivw_mre.R")), hash_before))

# ---- changed settings archive the earlier outputs -----------------------------------------------
r3 <- run_script_task(con, study, 1, "Omega-3", "Triglycerides", mk(0.8))
arch <- list.files(file.path(study, "archive"), pattern = "^run_", full.names = TRUE)
expect("changed settings re-run and move earlier outputs to archive/run_<datetime>/", !r3$reused && length(arch) == 1L &&
       file.exists(file.path(arch, basename(pair), "02_scripts", "01_ivw_mre.R")) &&
       identical(sha256_file(file.path(arch, basename(pair), "02_scripts", "01_ivw_mre.R")), hash_before))
expect("the new script carries the new setting", any(grepl("ci_level <- 0.8", readLines(file.path(pair, "02_scripts", "01_ivw_mre.R")))))
expect("manifest is consistent after the re-run", all(manifest_check(study)$status[!grepl("^archive/", manifest_check(study)$relative_path)] == "ok"))

# ---- the user's edits are never overwritten silently ------------------------------------------
sp <- file.path(pair, "02_scripts", "01_ivw_mre.R")
mine <- readBin(sp, "raw", file.size(sp))
writeBin(c(mine, charToRaw("# my own note\n")), sp); edited <- readBin(sp, "raw", file.size(sp))
expect("a script edited by the user is reported as modified", "modified" %in% manifest_check(study)$status[manifest_check(study)$relative_path == sprintf("03_pairs/%s/02_scripts/01_ivw_mre.R", basename(pair))])
e <- err_of(run_script_task(con, study, 1, "Omega-3", "Triglycerides", mk(0.8)))
expect("re-running refuses with SCRIPT_USER_MODIFIED and leaves the user's file untouched",
       grepl("SCRIPT_USER_MODIFIED", e) && identical(readBin(sp, "raw", file.size(sp)), edited))
r4 <- run_script_task(con, study, 1, "Omega-3", "Triglycerides", mk(0.8), on_modified = "archive")
arch2 <- list.files(file.path(study, "archive"), pattern = "^run_", full.names = TRUE)
expect("with 'archive' the user's version is kept in archive/ and the script is regenerated",
       !r4$reused && length(arch2) >= 1L && any(vapply(arch2, function(a) { f <- file.path(a, basename(pair), "02_scripts", "01_ivw_mre.R")
         file.exists(f) && identical(readBin(f, "raw", file.size(f)), edited) }, logical(1))) &&
       !any(grepl("my own note", readLines(sp))))

# ---- exploratory run without the replay check (D3) ----------------------------------------------
r5 <- run_script_task(con, study, 2, "Iron", "PCOS", mk(0.95), replay = FALSE)
expect("with the switch off the replay is recorded as not_run, never silently skipped",
       identical(r5$replay$status, "not_run") && any(grepl("not_run", DBI::dbGetQuery(con, "SELECT replay_json FROM results")$replay_json)))
expect("two pairs live side by side", dir.exists(file.path(study, "03_pairs", "P01_Omega-3__to__Triglycerides")) && dir.exists(file.path(study, "03_pairs", "P02_Iron__to__PCOS")))

# ---- a script that fails -----------------------------------------------------------------------
broken <- ADAPTER_IVW_MRE_V11
broken$emit_script <- function(spec, params, seed) {
  b <- ADAPTER_IVW_MRE_V11$emit_script(spec, params, seed)
  b$method_script <- paste0(b$method_script, "stop(\"broken on purpose\")\n"); b
}
r6 <- run_script_task(con, study, 3, "Folate", "CHD", list(list(adapter = broken, input = input, params = list())))
expect("a failing script gives status failed with SCRIPT_RUN_FAILED and no replay",
       r6$results$ivw_mre$status == "failed" && "SCRIPT_RUN_FAILED" %in% r6$results$ivw_mre$codes && r6$run_status != 0L &&
       DBI::dbGetQuery(con, "SELECT state FROM tasks WHERE task_id = 'P03_Folate__to__CHD:ivw_mre'")$state == "failed")

# ---- the database survives a restart -----------------------------------------------------------
DBI::dbDisconnect(con)
con2 <- db_open(file.path(study, ".engine"))
expect("after reopening, tasks and manifest rows are still there", DBI::dbGetQuery(con2, "SELECT COUNT(*) n FROM tasks")$n >= 5L &&
       DBI::dbGetQuery(con2, "SELECT COUNT(*) n FROM manifest")$n > 12L)
DBI::dbDisconnect(con2)

# keep the first pair for the cross-OS comparison
dir.create(file.path(dirname(out_dir), "keep"), showWarnings = FALSE)
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
