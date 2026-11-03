# MR GUI contracts v1.1 (DRAFT) - script emission, execution and replay check (C6.2, C6.5)
#
# Engine-side code. The delivered scripts never source this file: they call public
# package functions directly and read/write only files inside the pair folder.

# ---- number and string formatting for generated code --------------------------
.q <- function(x) encodeString(as.character(x), quote = "\"")      # R string literal
.num <- function(x) {                                              # R number literal, must round-trip
  s <- as.character(x)
  if (!identical(as.numeric(s), as.numeric(x))) stop("parameter value does not round-trip as text: ", x)
  s
}

# ---- analysis-data CSV: doubles written with 17 significant digits ------------
ANALYSIS_QUOTED_COLS <- c("variant_id", "ea", "oa", "harmonise_action")
write_analysis_csv <- function(data, path) {
  d <- data
  for (nm in names(d)) if (is.double(d[[nm]]))
    d[[nm]] <- ifelse(is.na(d[[nm]]), NA_character_, sprintf("%.17g", d[[nm]]))
  q <- which(names(d) %in% ANALYSIS_QUOTED_COLS)
  utils::write.table(d, path, sep = ",", row.names = FALSE, quote = q, na = "NA", qmethod = "double")
  invisible(path)
}

# ---- pieces shared by every generated method script ----------------------------
.script_lines_head <- function(adapter, spec) c(
  "# ==============================================================================",
  "# MR Studio generated script (contracts v1.1)",
  sprintf("# Method   : %s  [%s %s]", adapter$label, adapter$id, adapter$adapter_version),
  sprintf("# Package  : %s", paste(sprintf("TwoSampleMR commit %s", adapter$packages$TwoSampleMR), collapse = "; ")),
  sprintf("# Exposure : %s", spec$exposure_id),
  sprintf("# Outcome  : %s   (effects on the %s scale)", spec$outcome_id, spec$outcome_scale),
  "# Run from the pair folder (the folder that contains 01_data). In RStudio, open the .Rproj file.",
  "# Every argument below is written out, including defaults, so results do not depend on package defaults.",
  "# ==============================================================================",
  "")

.script_lines_guard <- function(spec) c(
  "# ---- 1. Working directory check ----------------------------------------------",
  sprintf("DATA_FILE <- file.path(\"01_data\", %s)", .q(spec$data_file)),
  "if (!file.exists(DATA_FILE))",
  "  stop(\"Wrong working directory. Set it to the pair folder (the folder that contains 01_data). In RStudio, open the .Rproj file.\")",
  "if (!requireNamespace(\"TwoSampleMR\", quietly = TRUE))",
  "  stop(\"Package TwoSampleMR is not installed. See README_RERUN.md (renv::restore()).\")",
  "")

.script_lines_repro <- function(seed) c(
  "# ---- 3. Reproducibility ------------------------------------------------------",
  "RNGkind(\"Mersenne-Twister\", \"Inversion\", \"Rejection\")",
  sprintf("set.seed(%d)", as.integer(seed)),
  "")

.script_lines_data <- c(
  "# ---- 4. Data (the exact table the calculation uses) ----------------------------",
  "dat <- utils::read.csv(DATA_FILE, stringsAsFactors = FALSE)",
  "use <- dat$mr_keep & is.finite(dat$beta_x) & is.finite(dat$beta_y) &",
  "       is.finite(dat$se_x) & is.finite(dat$se_y) & dat$se_y > 0",
  "dat <- dat[use, , drop = FALSE]   # variants kept after harmonisation",
  "")

.script_lines_write <- function(spec) c(
  "# ---- 6. Write results (read back by MR Studio; also your final result tables) --",
  "fmt <- function(x) if (is.double(x)) ifelse(is.na(x), NA_character_, sprintf(\"%.17g\", x)) else x",
  "dir.create(\"04_tables\", showWarnings = FALSE)",
  sprintf("utils::write.csv(as.data.frame(lapply(estimates, fmt)), file.path(\"04_tables\", %s), row.names = FALSE, na = \"NA\")", .q(spec$estimates_file)),
  sprintf("utils::write.csv(as.data.frame(lapply(diagnostics, fmt)), file.path(\"04_tables\", %s), row.names = FALSE, na = \"NA\")", .q(spec$diagnostics_file)),
  "print(estimates)",
  "print(diagnostics)",
  "",
  "# ---- 7. Session information --------------------------------------------------",
  "print(sessionInfo())",
  "")

.param_comment_lines <- function(adapter, params) {
  unlist(lapply(names(adapter$params), function(p) {
    sp <- adapter$params[[p]]; v <- params[[p]]
    c(sprintf("# %s: %s", if (!is.null(sp$label)) sp$label else p, if (!is.null(sp$hint)) sp$hint else ""),
      if (identical(sp$type, "enum")) sprintf("#   choices: %s", paste(sp$values, collapse = ", ")),
      sprintf("%s <- %s", p, if (is.character(v)) encodeString(v, quote = "\"") else .num(v)))
  }))
}

#' Bundle returned by emit_script (C6.2).
script_bundle <- function(method_script, packages, outputs, seed, data_file) {
  list(method_script = method_script, packages = packages, outputs = outputs, seed = seed, data_file = data_file)
}

#' Master script that runs every method script in order (00_run_all.R).
emit_master_script <- function(method_files) {
  paste0(paste(c(
    "# MR Studio: run all methods for this pair, in order. Run from the pair folder.",
    "if (!dir.exists(\"01_data\") || !dir.exists(\"02_scripts\"))",
    "  stop(\"Wrong working directory. To check these results in R, open VERIFY_IN_R.R in the pair folder and click Source.\")",
    unlist(lapply(method_files, function(f) c(
      sprintf("cat(\"\\n=== %s ===\\n\")", f),
      sprintf("source(file.path(\"02_scripts\", %s), echo = FALSE, print.eval = FALSE)", encodeString(f, quote = "\""))))),
    ""), collapse = "\n"))
}


# ---- files that let a user check the results in their own R (VERIFY_IN_R.R, .Rproj, README) -----
# The text is the same for every pair: the script finds its own folder and reads the packages from the method scripts.
VERIFY_SCRIPT_TEXT <- r"---(# ==============================================================================
# MR Studio: check these results in your own copy of R
#
# HOW TO USE
#   1. Open this file in RStudio.
#   2. Click the "Source" button at the top right of this editor. Nothing else is needed.
#
# WHAT IT DOES
#   - Copies the data and the scripts of this pair into a NEW folder: 09_verification/run_<date and time>
#   - Runs the scripts there, in a fresh R process.
#   - Compares every number with the results MR Studio saved in 04_tables.
#   Your original results are never changed or overwritten.
# ==============================================================================

local({
  lines <- character()
  say <- function(...) { x <- paste0(..., collapse = "\n"); lines <<- c(lines, x); cat(x, "\n", sep = "") }

  # ---- 1. Find the pair folder (the folder that holds this file) -----------------
  find_here <- function() {
    for (i in seq_len(sys.nframe())) {                       # sourced (the Source button in RStudio)
      of <- sys.frame(i)$ofile
      if (is.character(of) && length(of) == 1L) return(dirname(normalizePath(of, mustWork = FALSE)))
    }
    a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
    f <- gsub("~+~", " ", f, fixed = TRUE)                  # Rscript writes spaces in the path as ~+~
    if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {   # lines run one by one
      p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
      if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
    }
    getwd()
  }
  here <- find_here()
  if (!dir.exists(file.path(here, "01_data")) || !dir.exists(file.path(here, "02_scripts")))
    stop("Cannot find the pair folder. Open VERIFY_IN_R.R from inside the pair folder and click Source.", call. = FALSE)
  here <- normalizePath(here, winslash = "/")

  # ---- 2. Check the R packages the scripts need ----------------------------------
  scripts <- sort(list.files(file.path(here, "02_scripts"), pattern = "^[0-9]+_.*\\.R$"))
  scripts <- scripts[!grepl("^00_", scripts)]
  if (!length(scripts)) stop("There are no method scripts in 02_scripts.", call. = FALSE)
  pins <- list()
  for (s in scripts) {
    h <- grep("^# Package  :", readLines(file.path(here, "02_scripts", s), n = 12, warn = FALSE), value = TRUE)
    if (length(h)) for (part in trimws(strsplit(sub("^# Package  :\\s*", "", h[1]), ";")[[1]])) {
      m <- regmatches(part, regexec("^(\\S+) commit (\\S+)$", part))[[1]]
      if (length(m) == 3) pins[[m[2]]] <- m[3]
    }
  }
  repo <- c(TwoSampleMR = "MRCIEU/TwoSampleMR", MRPRESSO = "rondolab/MR-PRESSO", RadialMR = "WSpiller/RadialMR")
  absent <- names(pins)[!vapply(names(pins), requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    say("These R packages are not installed in your R yet: ", paste(absent, collapse = ", "))
    say("Install them once by running these lines in the Console, then click Source again:")
    say("  install.packages(\"remotes\")")
    for (p in absent) say("  remotes::install_github(\"", if (p %in% names(repo)) repo[[p]] else p, "@", pins[[p]], "\")")
    stop("Missing packages. See the lines above.", call. = FALSE)
  }

  # ---- 3. Copy data and scripts into a new folder --------------------------------
  run <- file.path(here, "09_verification", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  k <- 1L; base_run <- run; while (dir.exists(run)) { k <- k + 1L; run <- paste0(base_run, "_", k) }
  dir.create(run, recursive = TRUE)
  for (d in c("01_data", "02_scripts")) file.copy(file.path(here, d), run, recursive = TRUE)
  dir.create(file.path(run, "03_console")); dir.create(file.path(run, "04_tables"))

  # ---- 4. Run the scripts in a fresh R process (same way MR Studio runs them) ------
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  say("Running ", length(scripts), " method script(s) in a fresh R process. Slow methods can take a few minutes.")
  old <- setwd(run); on.exit(setwd(old), add = TRUE)
  env_new <- c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), RHOME = R.home(), OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")   # RHOME: Rscript finds this R even when it was moved
  env_old <- Sys.getenv(names(env_new), unset = NA)    # the same on every operating system; restored afterwards
  do.call(Sys.setenv, as.list(env_new))
  on.exit(for (n in names(env_old)) if (is.na(env_old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(env_old[[n]]), n)), add = TRUE)
  out <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(file.path("02_scripts", "00_run_all.R"))), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  writeLines(out, file.path(run, "03_console", "verify_console.txt"))
  if (status != 0L) {
    say("The scripts stopped with an error. Last lines of the output:")
    say(paste0("  ", utils::tail(out, 12)))
    stop("The run failed. The full output is in ", file.path(run, "03_console", "verify_console.txt"), call. = FALSE)
  }

  # ---- 5. Compare every number with the results MR Studio saved ---------------------
  tol <- 1e-8
  rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)
  files <- list.files(file.path(here, "04_tables"), pattern = "_(estimates|diagnostics)\\.csv$")
  worst <- 0; n <- 0L; bad <- character(); headline <- list()
  for (f in files) {
    pa <- file.path(here, "04_tables", f); pb <- file.path(run, "04_tables", f)
    if (!file.exists(pb)) { bad <- c(bad, paste0(f, " (not produced by this run)")); next }
    a <- rd(pa); b <- rd(pb)
    if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) { bad <- c(bad, paste0(f, " (different layout)")); next }
    for (nm in names(a)) {
      x <- a[[nm]]; y <- b[[nm]]
      if (is.numeric(x) && is.numeric(y)) {
        if (any(is.na(x) != is.na(y))) { bad <- c(bad, paste0(f, ", column ", nm, " (missing values differ)")); next }
        k <- !is.na(x)
        if (any(k)) {
          dif <- abs(x[k] - y[k])
          d <- dif / pmax(abs(x[k]), abs(y[k]), 1e-300); d[dif <= 1e-14] <- 0
          worst <- max(worst, d); n <- n + sum(k)
          if (any(d > tol)) bad <- c(bad, sprintf("%s, column %s (relative difference %.2e)", f, nm, max(d)))
        }
      } else if (!identical(as.character(x), as.character(y))) {
        bad <- c(bad, paste0(f, ", column ", nm, " (text differs)"))
      }
    }
    if (grepl("_estimates\\.csv$", f) && "estimand" %in% names(a)) {
      r <- which(a$estimand == "causal_effect")
      if (length(r)) headline[[f]] <- data.frame(method = sub("_estimates\\.csv$", "", f),
                                                 MR_Studio = a$estimate[r], This_run = b$estimate[r])
    }
  }

  # ---- 6. Report ----------------------------------------------------------------------
  say("")
  say("==================== RESULT ====================")
  say(R.version.string)
  for (p in names(pins)) say(p, " version in your R: ", as.character(utils::packageVersion(p)))
  if (length(headline)) {
    say("")
    say("Causal effect estimates (MR Studio against this run):")
    say(paste0("  ", utils::capture.output(print(do.call(rbind, headline), row.names = FALSE, digits = 10))))
  }
  say("")
  if (!length(files)) {
    say("RESULT: no saved result tables were found in 04_tables, so there is nothing to compare.")
  } else if (!length(bad)) {
    say("RESULT: ALL ", n, " NUMBERS MATCH THE RESULTS SAVED BY MR STUDIO.")
    say("Largest relative difference: ", format(worst, digits = 3), "  (0 means bit for bit identical; the allowed limit is ", format(tol), ")")
  } else {
    say("RESULT: DIFFERENCES FOUND. These do not match the results saved by MR Studio:")
    say(paste0("  - ", bad))
  }
  say("")
  say("This run is saved in: ", run)
  writeLines(lines, file.path(run, "verification_result.txt"))
  invisible(list(match = length(files) > 0 && !length(bad), folder = run))
})
)---"

RPROJ_TEXT <- "Version: 1.0\n\nRestoreWorkspace: No\nSaveWorkspace: No\nAlwaysSaveHistory: No\n\nEncoding: UTF-8\n"

RERUN_README_TEXT <- paste(c(
  "How to check this pair in your own R",
  "====================================",
  "",
  "1. Open VERIFY_IN_R.R in RStudio (or double-click the .Rproj file, then open VERIFY_IN_R.R).",
  "2. Click Source. If R packages are missing, the script prints the exact lines to install them once.",
  "3. The script runs the saved scripts on a copy, in 09_verification/run_<date and time>,",
  "   and compares every number with 04_tables. Your original results are never changed.",
  "",
  "You do not need to set a working folder or copy anything.",
  ""), collapse = "\n")

#' Write the verification files at the top of a pair folder when they are missing. Existing files are never rewritten
#' (they are listed in the manifest, and a pair keeps the check it was delivered with). Older pair folders get them when opened.
ensure_verify_files <- function(pair_dir) {
  w <- function(name, text) if (!file.exists(file.path(pair_dir, name))) writeBin(charToRaw(enc2utf8(text)), file.path(pair_dir, name))
  w("VERIFY_IN_R.R", paste0(VERIFY_SCRIPT_TEXT, "\n"))
  w("README_RERUN.md", RERUN_README_TEXT)
  w("MR_Studio_pair.Rproj", RPROJ_TEXT)
  invisible(file.path(pair_dir, "VERIFY_IN_R.R"))
}

# ---- materialise a pair folder, run it in a fresh R, read results back --------
PAIR_DIRS <- c("01_data", "02_scripts", "03_console", "04_tables")

materialise_pair <- function(pair_dir, methods, seed) {
  # methods: list of list(adapter, input, params)
  for (d in PAIR_DIRS) dir.create(file.path(pair_dir, d), recursive = TRUE, showWarnings = FALSE)
  files <- character(); info <- list()
  for (i in seq_along(methods)) {
    m <- methods[[i]]; a <- m$adapter
    spec <- list(exposure_id = m$input$exposure_id, outcome_id = m$input$outcome_id,
                 outcome_scale = m$input$outcome_scale,
                 data_file = sprintf("analysis_data_%s.csv", a$id),
                 estimates_file = sprintf("%s_estimates.csv", a$id),
                 diagnostics_file = sprintf("%s_diagnostics.csv", a$id))
    write_analysis_csv(m$input$data, file.path(pair_dir, "01_data", spec$data_file))
    params <- resolve_params(a, m$params)
    b <- a$emit_script(spec, params, seed)
    fname <- sprintf("%02d_%s.R", i, a$id)
    writeBin(charToRaw(enc2utf8(b$method_script)), file.path(pair_dir, "02_scripts", fname))
    files <- c(files, fname); info[[a$id]] <- c(spec, list(script_file = fname, bundle = b))
  }
  writeBin(charToRaw(emit_master_script(files)), file.path(pair_dir, "02_scripts", "00_run_all.R"))
  ensure_verify_files(pair_dir)
  invisible(info)
}

#' Run 00_run_all.R in a fresh R process (--vanilla) with the pair folder as working directory.
default_rscript <- function() {
  p <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows" && file.exists(paste0(p, ".exe"))) paste0(p, ".exe") else p
}

run_pair_script <- function(pair_dir, rscript = default_rscript(), console_name = "run_all_console.txt") {
  pair_dir <- normalizePath(pair_dir)
  # one thread for the numerical libraries (same on every operating system); restored afterwards
  # RHOME tells Rscript where this R lives (a moved runtime has no valid built-in location)
  old <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RHOME"), unset = NA)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", RHOME = R.home())
  on.exit(for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, setNames(list(old[[n]]), n)), add = TRUE)
  out <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(file.path("02_scripts", "00_run_all.R"))),
                                  stdout = TRUE, stderr = TRUE, wait = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  writeLines(out, file.path(pair_dir, "03_console", console_name))
  list(status = status, console = out)
}

# run_pair_script must execute with the pair folder as the working directory
run_in_dir <- function(dir, expr) { old <- setwd(dir); on.exit(setwd(old)); force(expr) }

read_method_results <- function(pair_dir, method_id) {
  rd <- function(f) utils::read.csv(file.path(pair_dir, "04_tables", f), stringsAsFactors = FALSE)
  list(estimates = rd(sprintf("%s_estimates.csv", method_id)),
       diagnostics = rd(sprintf("%s_diagnostics.csv", method_id)))
}

# ---- C6.5 replay check -------------------------------------------------------
.max_rel_diff <- function(a, b) {
  stopifnot(identical(dim(a), dim(b)))
  num <- vapply(a, is.numeric, logical(1))
  d <- 0
  for (nm in names(a)[num]) {
    x <- a[[nm]]; y <- b[[nm]]
    both_na <- is.na(x) & is.na(y)
    if (any(is.na(x) != is.na(y))) return(Inf)
    ok <- !both_na
    if (any(ok)) d <- max(d, abs(x[ok] - y[ok]) / pmax(abs(y[ok]), .Machine$double.xmin))
  }
  d
}

#' Copy the delivered files to a temporary folder, run the master script in a fresh R,
#' and compare every numeric output with the stored values (FR-59). Returns the C7 `replay` object.
replay_check <- function(pair_dir, method_ids, stored, rscript = default_rscript(),
                         same_machine_digits = 15L) {
  tmp <- tempfile("replay-"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  for (d in c("01_data", "02_scripts")) file.copy(file.path(pair_dir, d), tmp, recursive = TRUE)
  dir.create(file.path(tmp, "03_console")); dir.create(file.path(tmp, "04_tables"))   # outputs start empty
  run <- run_in_dir(tmp, run_pair_script(tmp, rscript, "replay_console.txt"))
  if (run$status != 0L) return(list(status = "failed", checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                                    r_version = paste(R.version$major, R.version$minor, sep = "."), os = os_family_id(),
                                    max_rel_diff = NA_real_, tolerance = 10^-(same_machine_digits - 1L),
                                    n_values_compared = 0L, console_tail = utils::tail(run$console, 15)))
  worst <- 0; n <- 0L
  for (id in method_ids) {
    now <- read_method_results(tmp, id)
    for (part in c("estimates", "diagnostics")) {
      worst <- max(worst, .max_rel_diff(now[[part]], stored[[id]][[part]]))
      n <- n + sum(vapply(now[[part]], is.numeric, logical(1))) * nrow(now[[part]])
    }
  }
  bitwise <- identical(worst, 0)
  tol <- 10^-(same_machine_digits - 1L)           # 15 significant digits
  status <- if (bitwise) "identical" else if (worst <= tol) "within_tolerance" else "differs"
  list(status = status, checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
       r_version = paste(R.version$major, R.version$minor, sep = "."), os = os_family_id(),
       max_rel_diff = worst, tolerance = tol, n_values_compared = as.integer(n))
}
