# Headless tests for the screen logic and the Shiny server (no browser, no desktop shell).
# Run:  Rscript prototype/engine/test_studio.R        (from the repository root)
args_dir <- normalizePath(file.path("prototype", "engine"), winslash = "/")
source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir)); source(file.path(args_dir, "app.R"))
root <- file.path(tempfile("studies-"), "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root)

n_pass <- 0L; n_fail <- 0L
expect <- function(desc, ok) { if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("PASS ", desc, "\n") } else { n_fail <<- n_fail + 1L; cat("FAIL ", desc, "\n") } }
err_of <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

expect("studies root honours MRGUI_STUDIES_DIR", identical(studies_root(), root))
expect("study names are cleaned for the file system", identical(clean_study_name("  Omega-3: TG/risk? "), "Omega-3 TG risk") && identical(clean_study_name("..."), ""))
expect("empty study name is refused", grepl("Give the study a name", err_of(studio_new("   "))))
expect("a study in a cloud folder is refused", grepl("SYS_STUDY_FOLDER_REFUSED", err_of(studio_new("S", file.path(tempdir(), "OneDrive")))))
d <- demo_pair(args_dir)
expect("demo data: 40 simulated instruments with the contract columns",
       nrow(d$input$data) == 40L && all(c("beta_x", "se_x", "beta_y", "se_y", "mr_keep") %in% names(d$input$data)) && is.logical(d$input$data$mr_keep))

h <- studio_new("Test study")
expect("new study has the folder tree and appears in the list", all(dir.exists(file.path(h$dir, STUDY_FOLDERS))) && "Test study" %in% studies_list())
expect("creating the same name twice is refused", grepl("already exists", err_of(studio_new("Test study"))))
r <- studio_run(h, "Demo exposure", "Demo outcome", d$input)
expect("first run finished with both methods ok and an identical replay",
       !r$reused && all(vapply(r$results, function(x) x$status, "") == "ok") && identical(r$replay$status, "identical"))
res <- pair_results(h$dir, r$pair_folder)
expect("results table has IVW and Egger rows with estimates, intervals and P values",
       !is.null(res) && all(c("Method", "estimate", "ci_lower", "ci_upper", "p") %in% names(res)) && any(grepl("IVW", res$Method)) && any(grepl("Egger", res$Method)))
fr <- format_results(res)
expect("formatted table: readable headers, interval with level, P never printed as 0",
       identical(names(fr), c("Method", "Estimate", "Std. error", "Interval", "P value", "Variants", "Estimand")) &&
       all(grepl("^-?[0-9.]+ to -?[0-9.]+ \\(95%\\)$", fr$Interval[!is.na(res$ci_lower)])) && !any(fr$`P value` %in% c("0.000", "0")) &&
       "Pleiotropy intercept" %in% fr$Estimand)
expect("size formatting", identical(fmt_bytes(c(500, 2048, 3 * 1024^2)), c("500 B", "2.0 KB", "3.0 MB")))
expect("the IVW estimate recovers the simulated effect (0.25) within 0.1", abs(res$estimate[grepl("IVW", res$Method)][1] - 0.25) < 0.1)
expect("same labels give the same pair folder and the second run is a reuse",
       identical(studio_run(h, "Demo exposure", "Demo outcome", d$input)$pair_folder, r$pair_folder) && isTRUE(studio_run(h, "Demo exposure", "Demo outcome", d$input)$reused))
r2 <- studio_run(h, "Demo exposure", "Demo outcome", d$input, ci_level = 0.9)
expect("a changed confidence level re-runs (not a reuse) in the same pair folder", !isTRUE(r2$reused) && identical(r2$pair_folder, r$pair_folder))
r3 <- studio_run(h, "Other exposure", "Demo outcome", d$input, methods = "ivw_mre")
expect("a new pair gets the next index and one method can be run alone", grepl("^P02_", r3$pair_folder) && all(vapply(r3$results, function(x) x$status, "") == "ok") && length(r3$results) == 1L)
expect("only the chosen method was run for the new pair", identical(unique(pair_results(h$dir, r3$pair_folder)$Method), "IVW (multiplicative random effects)"))
expect("no method chosen is refused", grepl("at least one method", err_of(studio_run(h, "A", "B", d$input, character()))))
expect("file list comes from the manifest and is not empty", nrow(study_files(h$dir)) >= 11L)
expect("opening only works inside the studies root", grepl("Only study folders", err_of(open_in_file_manager(tempdir(), root))))
DBI::dbDisconnect(h$con)
h2 <- studio_open("Test study")
expect("reopening shows the same pairs and results", identical(studio_pair_folders(h2$dir), c(r$pair_folder, r3$pair_folder)) && !is.null(pair_results(h2$dir, r$pair_folder)))
DBI::dbDisconnect(h2$con)

# ---- the Shiny server, driven without a browser -------------------------------------------------------------------------
shiny::testServer(mr_server(args_dir), {
  session$setInputs(study_name = "Server study"); session$setInputs(study_create = 1)
  expect("server: creating a study sets the current study", !is.null(rv$h) && grepl("Server study$", rv$h$dir))
  session$setInputs(exp_label = "Demo exposure", out_label = "Demo outcome", methods = c("ivw_mre", "egger"), ci_level = 0.95)
  session$setInputs(run = 1)
  expect("server: run finishes with a success message", identical(rv$run_class, "ok") && grepl("finished", rv$run_msg))
  session$setInputs(pair_pick = studio_pair_folders(rv$h$dir)[1])
  expect("server: results table renders with IVW and Egger rows", grepl("<table", output$results_tbl) && grepl("IVW", output$results_tbl) && grepl("Egger", output$results_tbl))
  session$setInputs(run = 2)
  expect("server: running again without changes reports a reuse", grepl("reused", rv$run_msg))
  session$setInputs(exp_label = "", run = 3)
  expect("server: missing labels give an error message, not a crash", identical(rv$run_class, "err"))
  session$setInputs(study_name = "Server study", study_create = 2)
  expect("server: duplicate study name gives an error message", identical(rv$msg_class, "err") && grepl("already exists", rv$msg))
})
shiny::testServer(mr_server(args_dir), {
  session$setInputs(exp_label = "A", out_label = "B", methods = "ivw_mre", ci_level = 0.95, run = 1)
  expect("server: run without a study asks to create or open one", grepl("first", rv$run_msg) && identical(rv$run_class, "err"))
})
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail)); if (n_fail) quit(status = 1)
