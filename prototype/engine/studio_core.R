# MR Studio engine: study logic used by the screen. No Shiny code in this file, so it can be tested without a browser.
# Depends on the contracts v1.1 R code (loaded by load_contracts()).

CONTRACT_FILES <- c("canonical.R", "adapter.R", "security.R", "runner.R", file.path("adapters", "reference_adapters.R"),
                    "script.R", file.path("adapters", "reference_adapters_v11.R"), "study.R")

#' Find the contracts folder: MRGUI_CONTRACTS_DIR, the installed app (next to the engine), or the repository.
find_contracts_dir <- function(engine_dir) {
  cands <- c(Sys.getenv("MRGUI_CONTRACTS_DIR", ""), file.path(engine_dir, "..", "contracts"),
             file.path(engine_dir, "..", "..", "contracts"))
  for (d in cands[nzchar(cands)]) if (file.exists(file.path(d, "r", "study.R"))) return(normalizePath(d, winslash = "/"))
  stop("The contracts folder was not found. Looked in: ", paste(cands[nzchar(cands)], collapse = "; "))
}

load_contracts <- function(contracts_dir) {
  for (f in CONTRACT_FILES) source(file.path(contracts_dir, "r", f), local = .GlobalEnv)
  options(mrgui.schema = file.path(contracts_dir, "sql", "schema_v2.sql"))
  invisible(contracts_dir)
}

#' Default place for studies (decision D2). MRGUI_STUDIES_DIR overrides it (tests, support).
studies_root <- function() {
  r <- Sys.getenv("MRGUI_STUDIES_DIR", "")
  if (nzchar(r)) r else file.path(path.expand("~"), "MR Studio")
}

#' Study folder names are what the user sees in the file manager: keep them readable, remove only what a file system refuses.
clean_study_name <- function(x) {
  x <- trimws(gsub("[\\\\/:*?\"<>|[:cntrl:]]", " ", x %||% ""))
  x <- gsub("\\s+", " ", x)
  x <- sub("^\\.+", "", x); x <- sub("[ .]+$", "", x)
  substr(x, 1L, 60L)
}
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

studies_list <- function(root = studies_root()) {
  if (!dir.exists(root)) return(character())
  d <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  d[file.exists(file.path(root, d, ".engine", "project.sqlite"))]
}

#' Create a study folder under `root` and open its database.
studio_new <- function(name, root = studies_root()) {
  nm <- clean_study_name(name)
  if (!nzchar(nm)) stop("Give the study a name.")
  dir <- file.path(root, nm)
  if (dir.exists(dir)) stop("A folder named \"", nm, "\" already exists in ", root, ". Choose another name, or open that study.")
  study <- study_create(dir)
  list(dir = study, con = db_open(file.path(study, ".engine")))
}

studio_open <- function(name, root = studies_root()) {
  dir <- file.path(root, name)
  if (!file.exists(file.path(dir, ".engine", "project.sqlite"))) stop("This folder is not a study: ", dir)
  study <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  con <- db_open(file.path(study, ".engine"))
  db_recover(con)    # a task left "running" by a crash goes back to pending (contracts C8)
  list(dir = study, con = con)
}

studio_pair_folders <- function(study) {
  d <- list.dirs(file.path(study, "03_pairs"), recursive = FALSE, full.names = FALSE)
  sort(d[grepl("^P[0-9]+_", d)])
}

#' The same exposure and outcome labels always map to the same pair folder (so a second run is a reuse, not a duplicate).
studio_pair_index <- function(study, exposure_label, outcome_label) {
  existing <- studio_pair_folders(study)
  key <- sub("^P[0-9]+_", "", pair_folder_name(1L, exposure_label, outcome_label))
  hit <- which(sub("^P[0-9]+_", "", existing) == key)
  if (length(hit)) return(as.integer(sub("^P([0-9]+)_.*", "\\1", existing[hit[1]])))
  length(existing) + 1L
}

# ---- demo data (SIMULATED; the real public GWAS download is a later milestone) ------------------------------------------
demo_pair <- function(engine_dir) {
  f <- file.path(engine_dir, "demo", "demo_pair_iv.csv")
  if (!file.exists(f)) stop("Demo data file is missing: ", f)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$mr_keep <- as.logical(d$mr_keep)
  list(exposure_label = "Demo exposure", outcome_label = "Demo outcome",
       input = list(kind = "pair_iv", data = d, exposure_id = "demo:simulated-exposure",
                    outcome_id = "demo:simulated-outcome", outcome_scale = "beta"))
}

METHOD_CHOICES <- c("Inverse-variance weighted, multiplicative random effects" = "ivw_mre", "MR-Egger regression" = "egger")

studio_run <- function(h, exposure_label, outcome_label, input, methods = c("ivw_mre", "egger"), ci_level = 0.95) {
  if (!length(methods)) stop("Choose at least one method.")
  exposure_label <- trimws(exposure_label); outcome_label <- trimws(outcome_label)
  if (!nzchar(exposure_label) || !nzchar(outcome_label)) stop("Give the exposure and the outcome a name.")
  ad <- list(ivw_mre = ADAPTER_IVW_MRE_V11, egger = ADAPTER_EGGER_V11)
  spec <- lapply(methods, function(m) list(adapter = ad[[m]], input = input,
                                           params = if (m == "ivw_mre") list(ci_level = ci_level) else list()))
  idx <- studio_pair_index(h$dir, exposure_label, outcome_label)
  r <- run_script_task(h$con, h$dir, idx, exposure_label, outcome_label, spec)
  r$pair_folder <- basename(r$pair_dir)
  r
}

METHOD_NAMES <- c(ivw_mre = "IVW (multiplicative random effects)", egger = "MR-Egger")

#' Estimates of a pair as one table, read from the delivered CSV files (works for a reopened study too).
pair_results <- function(study, pair_folder) {
  tdir <- file.path(study, "03_pairs", pair_folder, "04_tables")
  out <- list()
  for (m in names(METHOD_NAMES)) {
    f <- file.path(tdir, paste0(m, "_estimates.csv"))
    if (!file.exists(f)) next
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    keep <- intersect(c("estimand", "estimate", "se", "ci_lower", "ci_upper", "ci_level", "p", "n_ivs_used"), names(d))
    out[[m]] <- cbind(Method = METHOD_NAMES[[m]], d[, keep, drop = FALSE])
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}

fmt_num <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "fg", flag = "#"))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, formatC(p, digits = 2, format = "e"), formatC(p, digits = 3, format = "f")))
fmt_bytes <- function(b) ifelse(b < 1024, paste0(b, " B"), ifelse(b < 1024^2, sprintf("%.1f KB", b / 1024), sprintf("%.1f MB", b / 1024^2)))

#' The results table as the screen shows it: readable headers, interval in one column, P values that never print as 0.
format_results <- function(res) {
  est <- c(causal_effect = "Causal effect", pleiotropy_intercept = "Pleiotropy intercept")
  lvl <- if ("ci_level" %in% names(res)) res$ci_level else rep(0.95, nrow(res))
  data.frame(Method = res$Method, Estimate = fmt_num(res$estimate), `Std. error` = fmt_num(res$se),
             Interval = ifelse(is.na(res$ci_lower), "", sprintf("%s to %s (%g%%)", fmt_num(res$ci_lower), fmt_num(res$ci_upper), 100 * lvl)),
             `P value` = fmt_p(res$p), Variants = as.character(res$n_ivs_used),
             Estimand = unname(ifelse(res$estimand %in% names(est), est[res$estimand], res$estimand)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

pair_replay <- function(study, pair_folder) {
  f <- file.path(study, "03_pairs", pair_folder, "07_replay_check", "replay_check.json")
  if (!file.exists(f)) return(NULL)
  jsonlite::fromJSON(f)
}

study_files <- function(study) {
  m <- manifest_read(study)
  if (!nrow(m)) return(m)
  m <- m[!duplicated(m$relative_path, fromLast = TRUE), c("relative_path", "role", "bytes"), drop = FALSE]
  rownames(m) <- NULL
  m
}

#' Open a folder in Finder or Explorer. Only paths inside the studies root are opened.
open_in_file_manager <- function(path, root = studies_root()) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  r <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!startsWith(p, paste0(r, "/")) && !identical(p, r)) stop("Only study folders can be opened from here.")
  if (.Platform$OS.type == "windows") utils::shell.exec(normalizePath(p))
  else system2(if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open", shQuote(p), wait = FALSE)
  invisible(TRUE)
}
