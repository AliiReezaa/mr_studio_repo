# MR Studio engine: page and bridge. The screens are drawn by ui/studio.js; R does the work and answers with messages.
# Everything is served inline (no resource paths), so the gate rule "user data is never a resource path" holds.
library(shiny)

read_ui_file <- function(engine_dir, f) paste(readLines(file.path(engine_dir, "ui", f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

mr_ui_for <- function(engine_dir) {
  css <- read_ui_file(engine_dir, "studio.css"); js <- read_ui_file(engine_dir, "studio.js")
  function(req) tags$html(lang = "en",
    tags$head(tags$meta(charset = "utf-8"), tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"), tags$title("MR Studio"),
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap"),
      tags$style(HTML(css))),
    tags$body(
      div(id = "root"),
      div(style = "position:fixed;right:10px;bottom:6px;font-size:10px;color:#8a97a6;opacity:.7;pointer-events:none", "MR GUI engine prototype. Engine: ", textOutput("hello", inline = TRUE)),
      tags$script(HTML(js))))
}

PARAM_SCHEMAS <- list(
  ivw_mre = list(ci_level = list(label = "Confidence level", hint = "Width of the interval reported for the estimate.", min = 0.5, max = 0.999, default = 0.95)),
  egger = list())

mr_server <- function(engine_dir) {
  function(input, output, session) {
    output$hello <- renderText("connected")
    rv <- reactiveValues(h = NULL)
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
      ev(t, name = name, dir = h$dir, root = root, plan = plan, pairs = pairs, payloads = payloads, to = to, studies = studies_overview(root))
    }

    observeEvent(input$act, {
      a <- input$act; t <- a$t
      tryCatch((function() switch(t,
        init = {
          s <- settings_read(); root <- studies_root()
          ev("init", settings = s, studies = studies_overview(root), catalogue = demo_catalogue(), schemas = PARAM_SCHEMAS, root = root, diag = run_diagnostics(root))
        },
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
          close_study(); h <- studio_open(a$name); rv$h <- h
          study_msg("study_opened", h, to = a$to %||% "auto")
        },
        facts = {
          f <- lapply(a$pairs, function(p) pair_facts(p[[1]], p[[2]]))
          ev("facts", facts = f)
        },
        run_pair = {
          req(rv$h)
          inp <- demo_input(a$ex, a$out)
          meths <- unlist(a$methods); ci <- as.numeric(a$ci %||% 0.95)
          # The contract slugs labels to 18 characters, so two different traits can share a folder name. The pair number keeps them apart.
          old <- Filter(function(p) identical(p$ex, a$ex) && identical(p$out, a$out), ui_state_get(rv$h$con, "pairs") %||% list())
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          r <- tryCatch(studio_run(rv$h, inp$exposure_label, inp$outcome_label, inp$input, meths, ci, idx = idx), error = function(e) e)
          if (inherits(r, "error")) { ev("pair_done", i = a$i, ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("pair_done", i = a$i, ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "), paste(unlist(lapply(r$results[bad], `[[`, "messages")), collapse = " "))); return() }
          pl <- ui_state_get(rv$h$con, "pairs") %||% list()
          pl <- Filter(function(p) !(identical(p$ex, a$ex) && identical(p$out, a$out)), pl)
          pl[[length(pl) + 1L]] <- list(ex = a$ex, out = a$out, exName = inp$exposure_label, outName = inp$outcome_label, folder = r$pair_folder, methods = as.list(meths), ci = ci)
          ui_state_set(rv$h$con, "pairs", pl)
          ev("pair_done", i = a$i, ok = TRUE, folder = r$pair_folder, reused = isTRUE(r$reused), payload = pair_payload(rv$h$dir, r$pair_folder))
        },
        open_folder = { req(rv$h); open_in_file_manager(rv$h$dir) },
        NULL))(), error = function(e) {
          if (identical(t, "set") && identical(a$key, "studies_root")) ev("set_err", key = "studies_root", msg = conditionMessage(e)) else fail(e)
        })
    })
  }
}
