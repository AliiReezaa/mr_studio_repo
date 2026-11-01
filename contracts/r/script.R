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
    sp <- adapter$params[[p]]
    c(sprintf("# %s: %s", if (!is.null(sp$label)) sp$label else p, if (!is.null(sp$hint)) sp$hint else ""),
      sprintf("%s <- %s", p, .num(params[[p]])))
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
    "  stop(\"Wrong working directory. Set it to the pair folder (the folder that contains 01_data).\")",
    unlist(lapply(method_files, function(f) c(
      sprintf("cat(\"\\n=== %s ===\\n\")", f),
      sprintf("source(file.path(\"02_scripts\", %s), echo = FALSE, print.eval = FALSE)", encodeString(f, quote = "\""))))),
    ""), collapse = "\n"))
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
  old <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"), unset = NA)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
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
