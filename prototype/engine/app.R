# MR Studio engine: page and bridge. The screens are drawn by ui/studio.js; R does the work and answers with messages.
# Everything is served inline (no resource paths), so the gate rule "user data is never a resource path" holds.
library(shiny)
options(shiny.maxRequestSize = 2 * 1024^3)   # imports up to 2 GB; the uploaded copy is only a temporary file

read_ui_file <- function(engine_dir, f) paste(readLines(file.path(engine_dir, "ui", f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

mr_ui_for <- function(engine_dir) {
  css <- read_ui_file(engine_dir, "studio.css"); js <- read_ui_file(engine_dir, "studio.js")
  function(req) tags$html(lang = "en",
    tags$head(tags$meta(charset = "utf-8"), tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"), tags$title("MR Studio"),
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap"),
      tags$style(HTML(css))),
    tags$body(
      div(id = "root"),
      div(style = "position:absolute;left:-9999px;width:1px;height:1px;overflow:hidden", fileInput("impfile", "Choose a GWAS file", accept = NULL)),
      div(style = "position:fixed;right:10px;bottom:6px;font-size:10px;color:#8a97a6;opacity:.7;pointer-events:none", "MR GUI engine prototype. Engine: ", textOutput("hello", inline = TRUE)),
      tags$script(HTML(js))))
}

# The Customise panel is generated from the adapters' own parameter schemas (FR-62), so a label or hint is written once.
#' Settings sent by the screen, checked against the adapters before anything runs (unknown names are dropped, bad values are refused later by the adapter).
run_params <- function(x) {
  out <- list(ci_level = 0.95)
  x <- as.list(x)
  if (!is.null(x$ci_level)) out$ci_level <- as.numeric(x$ci_level)
  if (!is.null(x$min_f)) out$min_f <- as.numeric(x$min_f)
  if (!is.null(x$ivw_model)) out$ivw_model <- as.character(x$ivw_model)
  out
}

param_schemas <- function() {
  one <- function(ad) lapply(ad$params, function(p) list(label = p$label, hint = p$hint, min = p$min, max = p$max, default = p$default, type = p$type,
                                                         values = if (!is.null(p$values)) as.list(p$values), value_labels = if (!is.null(p$value_labels)) as.list(p$value_labels), value_hints = if (!is.null(p$value_hints)) as.list(p$value_hints)))
  list(ivw_mre = one(ADAPTER_IVW_MRE_V11), egger = one(ADAPTER_EGGER_V11))
}

mr_server <- function(engine_dir) {
  function(input, output, session) {
    output$hello <- renderText("connected")
    rv <- reactiveValues(h = NULL, imp = NULL)
    ev <- function(t, ...) session$sendCustomMessage("ev", c(list(t = t), list(...)))
    close_study <- function() { if (!is.null(rv$h)) try(DBI::dbDisconnect(rv$h$con), silent = TRUE); rv$h <- NULL }
    onSessionEnded(function() isolate(close_study()))
    fail <- function(e) ev("error", msg = conditionMessage(e))

    study_msg <- function(t, h, to = "auto") {
      root <- studies_root(); name <- basename(h$dir)
      plan <- ui_state_get(h$con, "plan") %||% list()
      pl <- ui_state_get(h$con, "pairs") %||% list()
      payloads <- list()
      pairs <- lapply(pl, function(p) {
        if (!is.null(p$folder) && dir.exists(file.path(h$dir, "03_pairs", p$folder)) && !is.null(pair_results(h$dir, p$folder))) payloads[[p$folder]] <<- pair_payload(h$dir, p$folder)
        p
      })
      pairs <- Filter(function(p) !is.null(p$folder), pairs)
      sens <- ui_state_get(h$con, "sens") %||% list()
      for (sn in sens) if (!is.null(sn$folder) && dir.exists(file.path(h$dir, "03_pairs", sn$folder)) && !is.null(pair_results(h$dir, sn$folder))) payloads[[sn$folder]] <- pair_payload(h$dir, sn$folder)
      ev(t, name = name, dir = h$dir, root = root, plan = plan, pairs = pairs, sens = sens, payloads = payloads, to = to, studies = studies_overview(root))
    }

    observeEvent(input$impfile, {
      f <- input$impfile
      tryCatch({
        hd <- import_read_head(f$datapath, f$name)
        rv$imp <- list(path = f$datapath, name = f$name, sep = hd$sep, columns = hd$columns, fingerprint = hd$fingerprint)
        tm <- template_match(hd$fingerprint)
        ev("import_head", name = hd$name, size = hd$size, sep = hd$sep, columns = as.list(hd$columns), rows = unname(hd$rows),
           suggest = import_suggest(hd$columns, hd$rows), template = tm)
      }, error = function(e) { rv$imp <- NULL; ev("import_error", msg = conditionMessage(e)) })
    })

    observeEvent(input$act, {
      a <- input$act; t <- a$t
      tryCatch((function() switch(t,
        init = {
          s <- settings_read(); root <- studies_root()
          ev("init", settings = s, library = library_list(), templates = template_list(), studies = studies_overview(root), catalogue = demo_catalogue(), schemas = param_schemas(), root = root, diag = run_diagnostics(root))
        },
        study_folder = open_in_file_manager(file.path(studies_root(), basename(a$name))),
        library = ev("library", items = library_list(), templates = template_list()),
        library_folder = { d <- library_dir(); dir.create(d, recursive = TRUE, showWarnings = FALSE); open_in_file_manager(d) },
        import_check = {
          req(rv$imp); v <- import_validate(rv$imp$path, rv$imp$sep, rv$imp$columns, a$mapping, a$facts)
          ev("import_report", ok = v$ok, errors = as.list(v$errors), warnings = as.list(v$warnings), stats = v$stats)
        },
        import_save = {
          req(rv$imp)
          rec <- import_save(rv$imp$path, rv$imp$name, rv$imp$sep, rv$imp$columns, a$mapping, a$facts)
          tn <- if (nzchar(a$template %||% "")) template_save(a$template, rv$imp$fingerprint, a$mapping, a$facts) else ""
          ev("import_saved", rec = rec, template = tn, items = library_list(), templates = template_list())
        },
        template_delete = { template_delete(a$name); ev("library", items = library_list(), templates = template_list()) },
        studies = ev("studies", studies = studies_overview(), root = studies_root()),
        set = {
          s <- settings_set(a$key, a$value)
          if (a$key == "studies_root") ev("set_ok", settings = s, root = studies_root())
        },
        settings_info = {
          root <- studies_root(); b <- if (dir.exists(root)) folder_size_bytes(root) else 0
          ev("settings_info", info = list(size = fmt_bytes(b)))
        },
        diagnostics = ev("diagnostics", items = run_diagnostics()),
        study_create = {
          nm <- clean_study_name(a$name)
          if (!nzchar(nm)) { ev("name_error", msg = "Give the study a name."); return() }
          if (dir.exists(file.path(studies_root(), nm))) { ev("name_error", msg = paste0("A study named \"", nm, "\" already exists. Choose another name, or open it from Studies.")); return() }
          close_study(); h <- studio_new(nm); rv$h <- h
          study_plan_save(h, list(kind = a$kind, aim = a$aim %||% "", prereg = isTRUE(a$prereg), question_kind = a$kind))
          study_msg("study_created", h)
        },
        plan_update = { req(rv$h); study_plan_save(rv$h, list(kind = a$kind, aim = a$aim %||% "", prereg = isTRUE(a$prereg))) },
        study_open = {
          close_study()
          h <- tryCatch(studio_open(a$name), error = function(e) { ev("study_open_failed", name = a$name, msg = conditionMessage(e)); NULL })
          if (is.null(h)) return(invisible())
          rv$h <- h
          study_msg("study_opened", h, to = a$to %||% "auto")
        },
        facts = {
          f <- lapply(a$pairs, function(p) pair_facts(p[[1]], p[[2]]))
          ev("facts", facts = f)
        },
        run_pair = {
          req(rv$h)
          inp <- demo_input(a$ex, a$out)
          meths <- unlist(a$methods); prm <- run_params(a$params); ci <- prm$ci_level
          # The contract slugs labels to 30 characters, so two different traits can share a folder name. The pair number keeps them apart.
          old <- Filter(function(p) identical(p$ex, a$ex) && identical(p$out, a$out), ui_state_get(rv$h$con, "pairs") %||% list())
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          r <- tryCatch(studio_run(rv$h, inp$exposure_label, inp$outcome_label, inp$input, meths, ci, idx = idx, params = prm), error = function(e) e)
          if (inherits(r, "error")) { ev("pair_done", i = a$i, ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("pair_done", i = a$i, ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "), paste(unlist(lapply(r$results[bad], `[[`, "messages")), collapse = " "))); return() }
          pl <- ui_state_get(rv$h$con, "pairs") %||% list()
          pl <- Filter(function(p) !(identical(p$ex, a$ex) && identical(p$out, a$out)), pl)
          pl[[length(pl) + 1L]] <- list(ex = a$ex, out = a$out, exName = inp$exposure_label, outName = inp$outcome_label, folder = r$pair_folder, methods = as.list(meths), ci = ci, params = as.list(prm))
          ui_state_set(rv$h$con, "pairs", pl)
          ev("pair_done", i = a$i, ok = TRUE, folder = r$pair_folder, reused = isTRUE(r$reused), payload = pair_payload(rv$h$dir, r$pair_folder))
        },
        instruments = ev("instruments", folder = a$folder, rows = { it <- instrument_table(a$ex, a$out); unname(lapply(seq_len(nrow(it)), function(i) as.list(it[i, , drop = FALSE]))) }),
        run_sens = {
          req(rv$h)
          rem <- sort(unlist(a$remove)); parent <- a$folder
          if (!length(rem)) stop("Choose at least one variant to leave out.")
          sens <- ui_state_get(rv$h$con, "sens") %||% list()
          old <- Filter(function(x) identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem), sens)
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          r <- tryCatch(studio_run_without(rv$h, a$ex, a$out, rem, unlist(a$methods), run_params(a$params)$ci_level, idx, params = run_params(a$params)), error = function(e) e)
          if (inherits(r, "error")) { ev("sens_done", ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("sens_done", ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "))); return() }
          sens <- Filter(function(x) !(identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem)), sens)
          rec <- list(parent = parent, remove = as.list(rem), folder = r$pair_folder, label = r$label, kept = r$kept)
          sens[[length(sens) + 1L]] <- rec; ui_state_set(rv$h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(rv$h$dir, r$pair_folder), reused = isTRUE(r$reused))
        },
        open_folder = { req(rv$h); open_in_file_manager(rv$h$dir) },
        open_rstudio = { req(rv$h); r <- open_in_rstudio(rv$h$dir, a$folder); ev("rstudio", opened = isTRUE(r$opened), script = r$script) },
        open_script = { req(rv$h); sp <- file.path(rv$h$dir, "03_pairs", a$folder, "02_scripts", "00_run_all.R"); open_in_file_manager(sp) },
        compare_default = {
          req(rv$h); pay <- pair_payload(rv$h$dir, a$folder)
          ev("compare", folder = a$folder, rows = ci_default_rows(pay$estimates))
        },
        NULL))(), error = function(e) {
          if (identical(t, "set") && identical(a$key, "studies_root")) ev("set_err", key = "studies_root", msg = conditionMessage(e)) else fail(e)
        })
    })
  }
}
