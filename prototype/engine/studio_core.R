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
  if (nzchar(r)) return(r)
  s <- settings_read()$studies_root
  if (is.character(s) && length(s) == 1L && nzchar(s)) s else file.path(path.expand("~"), "MR Studio")
}

# ---- user settings (theme, toggles, studies folder): one small JSON file, remembered per user (UI-1) ----------------
settings_dir <- function() {
  d <- Sys.getenv("MRGUI_CONFIG_DIR", "")
  if (nzchar(d)) d else tools::R_user_dir("MRStudio", "config")
}
SETTINGS_DEFAULT <- list(theme = "aurora", setup_done = FALSE, offline = FALSE, update = TRUE, awake = TRUE, notify = TRUE,
                         studies_root = "")
SETTINGS_THEMES <- c("aurora", "paper", "graphite", "ocean", "rose", "lilac")
settings_read <- function() {
  f <- file.path(settings_dir(), "settings.json")
  s <- SETTINGS_DEFAULT
  if (file.exists(f)) {
    u <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) list())
    for (k in names(s)) if (!is.null(u[[k]]) && length(u[[k]]) == 1L && identical(class(u[[k]]), class(s[[k]]))) s[[k]] <- u[[k]]
  }
  if (!s$theme %in% SETTINGS_THEMES) s$theme <- "aurora"
  s
}
settings_set <- function(key, value) {
  if (!key %in% names(SETTINGS_DEFAULT)) stop("Unknown setting: ", key)
  if (!identical(class(value), class(SETTINGS_DEFAULT[[key]]))) stop("Wrong type for setting ", key)
  if (identical(key, "theme") && !value %in% SETTINGS_THEMES) stop("Unknown theme.")
  if (identical(key, "studies_root") && nzchar(value)) {
    prob <- study_path_problem(value); if (!is.null(prob)) stop("SYS_STUDY_FOLDER_REFUSED: ", prob, ". Choose a folder on this computer that is not synchronised.")
    check_path_length(value)
  }
  s <- settings_read(); s[[key]] <- value
  dir.create(settings_dir(), recursive = TRUE, showWarnings = FALSE)
  atomic_write_text(file.path(settings_dir(), "settings.json"), jsonlite::toJSON(s, auto_unbox = TRUE, pretty = TRUE))
  invisible(s)
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

# ---- demo catalogue: SIMULATED traits so every screen can be tested before the real catalogue exists --------------------
DEMO_EXPOSURES <- data.frame(id = c("expA", "expB", "expC", "expD"), name = paste("Simulated exposure", LETTERS[1:4]),
                             n = c(173082, 681275, 757601, 204402), stringsAsFactors = FALSE)
DEMO_OUTCOMES <- data.frame(id = c("outA", "outB", "outC", "outD"), name = paste("Simulated outcome", LETTERS[1:4]),
                            n = c(184305, 392423, 440328, 487511), stringsAsFactors = FALSE)
DEMO_EFFECT <- matrix(c(0.30, 0.00, -0.15, 0.10,
                        0.00, 0.25,  0.00, 0.05,
                        0.20, 0.00,  0.00, 0.35,
                        0.00, 0.10, -0.25, 0.00), nrow = 4, byrow = TRUE, dimnames = list(DEMO_EXPOSURES$id, DEMO_OUTCOMES$id))

demo_catalogue <- function() {
  mk <- function(d) lapply(seq_len(nrow(d)), function(i) list(id = d$id[i], name = d$name[i], src = "Demo, simulated", n = format(d$n[i], big.mark = ","), anc = "EUR"))
  list(exposures = mk(DEMO_EXPOSURES), outcomes = mk(DEMO_OUTCOMES))
}
.seed_of <- function(x) sum(utf8ToInt(x) * seq_along(utf8ToInt(x))) %% 2147483647

#' One simulated pair of 40 instruments. The instruments belong to the exposure; the outcome noise belongs to the pair.
demo_input <- function(exp_id, out_id) {
  ei <- match(exp_id, DEMO_EXPOSURES$id); oi <- match(out_id, DEMO_OUTCOMES$id)
  if (is.na(ei) || is.na(oi)) stop("Unknown demo trait.")
  k <- 40L
  set.seed(.seed_of(exp_id))
  bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE); sx <- runif(k, 0.005, 0.012); eaf <- round(runif(k, 0.1, 0.9), 3)
  set.seed(.seed_of(paste(exp_id, out_id)))
  sy <- runif(k, 0.008, 0.02); by <- DEMO_EFFECT[ei, oi] * bx + rnorm(k, 0.001, 0.006) + rnorm(k, 0, sy)
  d <- data.frame(variant_id = sprintf("demo:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G", beta_x = round(bx, 6), se_x = round(sx, 6),
                  eaf_x = eaf, beta_y = round(by, 6), se_y = round(sy, 6), eaf_y = eaf, p_x = signif(2 * pnorm(-abs(bx / sx)), 6),
                  p_y = signif(2 * pnorm(-abs(by / sy)), 6), n_x = DEMO_EXPOSURES$n[ei], n_y = DEMO_OUTCOMES$n[oi],
                  harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  list(exposure_label = DEMO_EXPOSURES$name[ei], outcome_label = DEMO_OUTCOMES$name[oi],
       input = list(kind = "pair_iv", data = d, exposure_id = paste0("demo:", exp_id), outcome_id = paste0("demo:", out_id), outcome_scale = "beta"))
}

#' What the engine can say about a pair before running it (Methods and Check screens). Real numbers from the data table.
pair_facts <- function(exp_id, out_id) {
  p <- demo_input(exp_id, out_id); d <- p$input$data; kept <- d[d$mr_keep, ]
  Fs <- (kept$beta_x / kept$se_x)^2
  pal <- sum((kept$ea == "A" & kept$oa == "T") | (kept$ea == "T" & kept$oa == "A") | (kept$ea == "C" & kept$oa == "G") | (kept$ea == "G" & kept$oa == "C"))
  list(exposure = p$exposure_label, outcome = p$outcome_label, instruments = nrow(kept), mean_f = round(mean(Fs), 1), min_f = round(min(Fs), 1),
       n_weak = sum(Fs < 10), palindromic = pal, scale = p$input$outcome_scale)
}

METHOD_CHOICES <- c("Inverse-variance weighted, multiplicative random effects" = "ivw_mre", "MR-Egger regression" = "egger")

studio_run <- function(h, exposure_label, outcome_label, input, methods = c("ivw_mre", "egger"), ci_level = 0.95, idx = NULL) {
  if (!length(methods)) stop("Choose at least one method.")
  exposure_label <- trimws(exposure_label); outcome_label <- trimws(outcome_label)
  if (!nzchar(exposure_label) || !nzchar(outcome_label)) stop("Give the exposure and the outcome a name.")
  ad <- list(ivw_mre = ADAPTER_IVW_MRE_V11, egger = ADAPTER_EGGER_V11)
  spec <- lapply(methods, function(m) list(adapter = ad[[m]], input = input,
                                           params = if (m == "ivw_mre") list(ci_level = ci_level) else list()))
  if (is.null(idx)) idx <- studio_pair_index(h$dir, exposure_label, outcome_label)
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


# ---- UI state kept in the study database (plan, pairs) --------------------------------------------------------------------
ui_state_get <- function(con, key) {
  r <- DBI::dbGetQuery(con, "SELECT value_json FROM ui_state WHERE key = ?", params = list(key))
  if (!nrow(r)) NULL else jsonlite::fromJSON(r$value_json, simplifyVector = FALSE)
}
ui_state_set <- function(con, key, value) {
  DBI::dbExecute(con, "INSERT OR REPLACE INTO ui_state(key, value_json, updated_at) VALUES (?, ?, ?)",
                 params = list(key, as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null")), format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  invisible(TRUE)
}

#' Save the question page of a new study (what the user typed) as a readable file and in the database.
study_plan_save <- function(h, plan) {
  plan$saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  rel <- "01_study_plan/study_plan.json"
  atomic_write_text(file.path(h$dir, rel), jsonlite::toJSON(plan, auto_unbox = TRUE, pretty = TRUE))
  manifest_register(h$dir, rel, role = "study_plan")
  ui_state_set(h$con, "plan", plan)
  invisible(plan)
}

age_label <- function(t) {
  d <- as.numeric(difftime(Sys.time(), t, units = "days"))
  if (d < 1) "today" else if (d < 2) "yesterday" else if (d < 14) paste(floor(d), "days ago") else if (d < 60) paste(floor(d / 7), "weeks ago") else format(t, "%b %Y")
}

#' Recent studies for Home and Studies: status comes from the files (results present means complete).
studies_overview <- function(root = studies_root()) {
  nm <- studies_list(root)
  out <- lapply(nm, function(n) {
    dir <- file.path(root, n)
    pf <- list.dirs(file.path(dir, "03_pairs"), recursive = FALSE, full.names = FALSE); pf <- pf[grepl("^P[0-9]+_", pf)]
    has_res <- any(vapply(pf, function(p) length(Sys.glob(file.path(dir, "03_pairs", p, "04_tables", "*_estimates.csv"))) > 0, logical(1)))
    nmeth <- length(unique(unlist(lapply(pf, function(p) sub("_estimates\\.csv$", "", basename(Sys.glob(file.path(dir, "03_pairs", p, "04_tables", "*_estimates.csv"))))))))
    mt <- file.mtime(file.path(dir, ".engine", "project.sqlite"))
    list(name = n, status = if (has_res) "COMPLETE" else "DRAFT", pairs = length(pf), methods = nmeth, when = age_label(mt), mtime = as.numeric(mt))
  })
  out[order(-vapply(out, function(x) x$mtime, 0))]
}

#' Everything the Results screen needs for one pair, read back from the delivered files.
pair_payload <- function(study, pair_folder) {
  pd <- file.path(study, "03_pairs", pair_folder)
  rd <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  est <- list(); dg <- list()
  for (m in names(METHOD_NAMES)) {
    e <- rd(file.path(pd, "04_tables", paste0(m, "_estimates.csv")))
    if (!is.null(e)) est[[m]] <- lapply(seq_len(nrow(e)), function(i) c(list(method = m, method_name = METHOD_NAMES[[m]]), as.list(e[i, , drop = FALSE])))
    g <- rd(file.path(pd, "04_tables", paste0(m, "_diagnostics.csv")))
    if (!is.null(g)) dg[[m]] <- lapply(seq_len(nrow(g)), function(i) c(list(method = m), as.list(g[i, , drop = FALSE])))
  }
  v <- rd(Sys.glob(file.path(pd, "01_data", "analysis_data_*.csv"))[1])
  vars <- if (!is.null(v)) lapply(seq_len(nrow(v)), function(i) list(bx = v$beta_x[i], by = v$beta_y[i], sx = v$se_x[i], sy = v$se_y[i], keep = isTRUE(as.logical(v$mr_keep[i])))) else list()
  list(folder = pair_folder, estimates = unname(unlist(est, recursive = FALSE)), diagnostics = unname(unlist(dg, recursive = FALSE)), variants = vars,
       replay = pair_replay(study, pair_folder))
}

folder_size_bytes <- function(path) {
  f <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  if (!length(f)) 0 else sum(file.size(f), na.rm = TRUE)
}

#' Self-check for Settings > Run diagnostics. Each item: name, ok, detail.
run_diagnostics <- function(root = studies_root()) {
  chk <- function(name, expr) { r <- tryCatch(expr, error = function(e) list(ok = FALSE, detail = conditionMessage(e))); c(list(name = name), r) }
  list(
    chk("R engine", list(ok = TRUE, detail = R.version.string)),
    chk("TwoSampleMR package", list(ok = requireNamespace("TwoSampleMR", quietly = TRUE), detail = if (requireNamespace("TwoSampleMR", quietly = TRUE)) paste("version", as.character(utils::packageVersion("TwoSampleMR"))) else "not installed")),
    chk("Script runner (Rscript)", { r <- default_rscript(); list(ok = file.exists(r), detail = if (file.exists(r)) "found" else paste("missing:", r)) }),
    chk("Studies folder is writable", { dir.create(root, recursive = TRUE, showWarnings = FALSE); f <- tempfile(".mrstudio-check-", root); ok <- tryCatch({ writeLines("x", f); TRUE }, error = function(e) FALSE, warning = function(w) FALSE); unlink(f); list(ok = ok, detail = root) }),
    chk("Studies folder is not synchronised", { p <- study_path_problem(root); list(ok = is.null(p), detail = if (is.null(p)) "local folder" else p) })
  )
}
