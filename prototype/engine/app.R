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
  if (!is.null(x$nboot)) out$nboot <- as.integer(round(as.numeric(x$nboot)))
  if (!is.null(x$phi)) out$phi <- as.numeric(x$phi)
  if (!is.null(x$nbdist)) out$nbdist <- as.integer(round(as.numeric(x$nbdist)))
  if (!is.null(x$signif)) out$signif <- as.numeric(x$signif)
  if (!is.null(x$alpha)) out$alpha <- as.numeric(x$alpha)
  if (!is.null(x$weights)) out$weights <- as.character(x$weights)
  if (!is.null(x$tol)) out$tol <- as.numeric(x$tol)
  if (!is.null(x$r_xxo)) out$r_xxo <- as.numeric(x$r_xxo)
  if (!is.null(x$r_yyo)) out$r_yyo <- as.numeric(x$r_yyo)
  out
}

param_schemas <- function() {
  one <- function(ad) lapply(ad$params, function(p) list(label = p$label, hint = p$hint, min = p$min, max = p$max, default = p$default, type = p$type,
                                                         values = if (!is.null(p$values)) as.list(p$values), value_labels = if (!is.null(p$value_labels)) as.list(p$value_labels), value_hints = if (!is.null(p$value_hints)) as.list(p$value_hints)))
  list(ivw_mre = one(ADAPTER_IVW_MRE_V11), egger = one(ADAPTER_EGGER_V11), wm = one(ADAPTER_WM_V11), smed = one(ADAPTER_SMED_V11), mode = one(ADAPTER_MODE_V11), smode = one(ADAPTER_SMODE_V11), presso = one(ADAPTER_PRESSO_V11), radial = one(ADAPTER_RADIAL_V11), radegger = one(ADAPTER_RADEGGER_V11), steiger = one(ADAPTER_STEIGER_V11))
}

mr_server <- function(engine_dir) {
  function(input, output, session) {
    output$hello <- renderText("connected")
    rv <- reactiveValues(h = NULL, imp = NULL)
    ev <- function(t, ...) session$sendCustomMessage("ev", c(list(t = t), list(...)))
    ogw_ev <- function() ev("ogw_state", token = ogw_token_info(ogw_token_get()), catalogue = ogw_catalogue_status(), offline = isTRUE(settings_read()$offline))
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
      ev(t, name = name, dir = h$dir, root = root, plan = plan, pairs = pairs, sens = sens, payloads = payloads, overlap = overlap_get(h), to = to, studies = studies_overview(root))
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
          ev("init", settings = s, library = library_list(), templates = template_list(), studies = studies_overview(root), catalogue = catalogue_all(), schemas = param_schemas(), root = root, diag = run_diagnostics(root))
          ogw_ev()
        },
        study_folder = open_in_file_manager(file.path(studies_root(), basename(a$name))),
        library = { ev("library", items = library_list(), templates = template_list(), catalogue = catalogue_all()); ogw_ev() },
        library_folder = { d <- library_dir(); dir.create(d, recursive = TRUE, showWarnings = FALSE); open_in_file_manager(d) },
        import_check = {
          req(rv$imp); v <- import_validate(rv$imp$path, rv$imp$sep, rv$imp$columns, a$mapping, a$facts)
          ev("import_report", ok = v$ok, errors = as.list(v$errors), warnings = as.list(v$warnings), stats = v$stats)
        },
        import_save = {
          req(rv$imp)
          rec <- import_save(rv$imp$path, rv$imp$name, rv$imp$sep, rv$imp$columns, a$mapping, a$facts)
          tn <- if (nzchar(a$template %||% "")) template_save(a$template, rv$imp$fingerprint, a$mapping, a$facts) else ""
          ev("import_saved", rec = rec, template = tn, items = library_list(), templates = template_list(), catalogue = catalogue_all())
        },
        template_delete = { template_delete(a$name); ev("library", items = library_list(), templates = template_list(), catalogue = catalogue_all()) },
        studies = ev("studies", studies = studies_overview(), root = studies_root()),
        set = {
          s <- settings_set(a$key, a$value)
          if (a$key == "studies_root") ev("set_ok", settings = s, root = studies_root())
        },
        settings_info = {
          root <- studies_root(); b <- if (dir.exists(root)) folder_size_bytes(root) else 0
          ev("settings_info", info = list(size = fmt_bytes(b))); ogw_ev()
        },
        diagnostics = ev("diagnostics", items = run_diagnostics(net = !isTRUE(settings_read()$offline))),
        study_create = {
          nm <- clean_study_name(a$name)
          if (!nzchar(nm)) { ev("name_error", msg = "Give the study a name."); return() }
          if (dir.exists(file.path(studies_root(), nm))) { ev("name_error", msg = paste0("A study named \"", nm, "\" already exists. Choose another name, or open it from Studies.")); return() }
          close_study(); h <- studio_new(nm); rv$h <- h
          study_plan_save(h, list(kind = a$kind, aim = a$aim %||% "", prereg = isTRUE(a$prereg), question_kind = a$kind))
          study_msg("study_created", h)
        },
        overlap_save = {
          req(rv$h)
          for (it in a$items) overlap_set(rv$h, it$ex, it$out, it$answer, it$pct, it$note %||% "")
          ev("overlap", overlap = overlap_get(rv$h))
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
          pt <- as.numeric(a$pthr %||% 5e-8); f <- lapply(a$pairs, function(p) tryCatch(pair_facts(p[[1]], p[[2]], pt), error = function(e) list(error = conditionMessage(e))))
          ev("facts", facts = f)
        },
        run_pair = {
          req(rv$h)
          pt <- as.numeric(a$pthr %||% 5e-8)
          inp <- tryCatch(pair_input(a$ex, a$out, pt), error = function(e) e)
          if (inherits(inp, "error")) { ev("pair_done", i = a$i, ok = FALSE, msg = conditionMessage(inp)); return() }
          meths <- unlist(a$methods); prm <- run_params(a$params); ci <- prm$ci_level
          # The contract slugs labels to 30 characters, so two different traits can share a folder name. The pair number keeps them apart.
          old <- Filter(function(p) identical(p$ex, a$ex) && identical(p$out, a$out), ui_state_get(rv$h$con, "pairs") %||% list())
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          sti <- if ("steiger" %in% meths) steiger_input(inp, list(x = a$prev$x, y = a$prev$y)) else NULL   # prevalence the analyst entered for a binary trait
          r <- tryCatch(studio_run(rv$h, inp$exposure_label, inp$outcome_label, inp$input, meths, ci, idx = idx, params = prm, steiger_in = sti), error = function(e) e)
          if (inherits(r, "error")) { ev("pair_done", i = a$i, ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("pair_done", i = a$i, ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "), paste(unlist(lapply(r$results[bad], `[[`, "messages")), collapse = " "))); return() }
          pl <- ui_state_get(rv$h$con, "pairs") %||% list()
          pl <- Filter(function(p) !(identical(p$ex, a$ex) && identical(p$out, a$out)), pl)
          pl[[length(pl) + 1L]] <- list(ex = a$ex, out = a$out, exName = inp$exposure_label, outName = inp$outcome_label, folder = r$pair_folder, methods = as.list(meths), ci = ci, params = as.list(prm), prev = list(x = a$prev$x, y = a$prev$y), pthr = pt, simulated = isTRUE(inp$simulated), ledger = inp$ledger, anc_x = inp$anc_x, anc_y = inp$anc_y)
          ui_state_set(rv$h$con, "pairs", pl)
          ev("pair_done", i = a$i, ok = TRUE, folder = r$pair_folder, reused = isTRUE(r$reused), payload = pair_payload(rv$h$dir, r$pair_folder))
        },
        instruments = ev("instruments", folder = a$folder, rows = { it <- instrument_table(a$ex, a$out, as.numeric(a$pthr %||% 5e-8)); unname(lapply(seq_len(nrow(it)), function(i) as.list(it[i, , drop = FALSE]))) }),
        run_sens = {
          req(rv$h)
          rem <- sort(unlist(a$remove)); parent <- a$folder
          if (!length(rem)) stop("Choose at least one variant to leave out.")
          sens <- ui_state_get(rv$h$con, "sens") %||% list()
          old <- Filter(function(x) identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem), sens)
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          r <- tryCatch(studio_run_without(rv$h, a$ex, a$out, rem, unlist(a$methods), run_params(a$params)$ci_level, idx, params = run_params(a$params), pthr = as.numeric(a$pthr %||% 5e-8)), error = function(e) e)
          if (inherits(r, "error")) { ev("sens_done", ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("sens_done", ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "))); return() }
          sens <- Filter(function(x) !(identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem)), sens)
          rec <- list(parent = parent, remove = as.list(rem), folder = r$pair_folder, label = r$label, kept = r$kept)
          sens[[length(sens) + 1L]] <- rec; ui_state_set(rv$h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(rv$h$dir, r$pair_folder), reused = isTRUE(r$reused))
        },
        build_package = {
          req(rv$h)
          r <- tryCatch(build_package(rv$h, unlist(a$items), a$width_mm, a$pub), error = function(e) e)
          if (inherits(r, "error")) { ev("package_done", ok = FALSE, msg = conditionMessage(r)); return() }
          ev("package_done", ok = TRUE, folder = r$folder, files = as.list(r$files), notes = as.list(r$notes), simulated = r$simulated)
        },
        open_package = { req(rv$h); open_in_file_manager(file.path(rv$h$dir, "05_publication", a$folder)) },
        open_folder = { req(rv$h); open_in_file_manager(rv$h$dir) },
        open_rstudio = { req(rv$h); r <- open_in_rstudio(rv$h$dir, a$folder); ev("rstudio", opened = isTRUE(r$opened), script = r$script, pair = r$pair) },
        open_script = { req(rv$h); sp <- verify_script_path(rv$h$dir, a$folder)$script; open_in_text_editor(sp) },
        compare_default = {
          req(rv$h)
          r <- tryCatch(compare_defaults(rv$h$dir, a$folder, a$ex, a$out, unlist(a$methods), as.numeric(a$pthr %||% 5e-8)), error = function(e) e)
          if (inherits(r, "error")) { ev("compare", folder = a$folder, ok = FALSE, msg = conditionMessage(r)); return() }
          ev("compare", folder = a$folder, ok = TRUE, rows = r$rows, rerun = r$rerun)
        },
        ogw_state = ogw_ev(),
        ogw_token_save = {
          tok <- ogw_clean_token(a$token)
          if (is.null(tok)) { ev("ogw_token", ok = FALSE, msg = "That does not look like an OpenGWAS token. A token is one long line made of three parts separated by dots, and it starts with eyJ. Copy all of it from api.opengwas.io."); return() }
          ogw_token_set(tok); ev("ogw_token", ok = TRUE, token = ogw_token_info(tok), msg = "Token kept in memory until you close MR Studio.")
        },
        ogw_token_clear = { ogw_token_set(NULL); ev("ogw_token", ok = TRUE, token = ogw_token_info(NULL), msg = "Token removed.") },
        ogw_token_test = {
          r <- tryCatch(ogw_token_test(ogw_token_get(), isTRUE(settings_read()$offline)), error = function(e) e)
          ev("ogw_token_tested", ok = !inherits(r, "error"), msg = if (inherits(r, "error")) conditionMessage(r) else "OpenGWAS accepted the token.")
        },
        ogw_catalogue_update = {
          r <- tryCatch(ogw_catalogue_update(ogw_token_get(), isTRUE(settings_read()$offline)), error = function(e) e)
          ev("ogw_cat", ok = !inherits(r, "error"), status = if (!inherits(r, "error")) r, msg = if (inherits(r, "error")) conditionMessage(r) else "")
        },
        ogw_search = {
          cat <- ogw_catalogue_read()
          if (is.null(cat)) { ev("ogw_results", q = a$q, total = 0L, rows = list(), nocat = TRUE); return() }
          res <- ogw_search(cat, a$q %||% "", a$pop %||% "", as.numeric(a$min_n %||% 0), !identical(a$mr_only, FALSE), a$kind %||% "any")
          pl <- ogw_search_payload(res); ev("ogw_results", q = a$q, total = pl$total, rows = pl$rows, nocat = FALSE)
        },
        ogw_dl_exposure = {
          r <- tryCatch(ogw_get_exposure(a$id, ogw_token_get(), as.numeric(a$p), as.numeric(a$r2), as.numeric(a$kb), a$pop, a$name, isTRUE(settings_read()$offline)), error = function(e) e)
          if (inherits(r, "error")) { ev("ogw_dl", ok = FALSE, msg = conditionMessage(r)); return() }
          ev("ogw_dl", ok = TRUE, rec = r, reused = isTRUE(r$reused), items = library_list(), templates = template_list(), catalogue = catalogue_all())
        },
        ogw_dl_outcome = {
          r <- tryCatch(ogw_get_outcome(a$id, a$exposure, ogw_token_get(), !identical(a$proxies, FALSE), as.numeric(a$rsq), as.numeric(a$maf), a$name, isTRUE(settings_read()$offline)), error = function(e) e)
          if (inherits(r, "error")) { ev("ogw_dl", ok = FALSE, msg = conditionMessage(r)); return() }
          ev("ogw_dl", ok = TRUE, rec = r, reused = isTRUE(r$reused), items = library_list(), templates = template_list(), catalogue = catalogue_all())
        },
        NULL))(), error = function(e) {
          if (identical(t, "set") && identical(a$key, "studies_root")) ev("set_err", key = "studies_root", msg = conditionMessage(e)) else fail(e)
        })
    })
  }
}
