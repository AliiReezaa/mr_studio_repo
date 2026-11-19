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
      div(style = "position:absolute;left:-9999px;width:1px;height:1px;overflow:hidden", fileInput("impfile", "Choose a GWAS file", accept = NULL), fileInput("mvcovfile", "Choose an individual-level file", accept = NULL)),
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
  if (!is.null(x$psi)) out$psi <- as.numeric(x$psi)   # Round 40
  if (!is.null(x$ci_step)) out$ci_step <- as.numeric(x$ci_step)
  if (!is.null(x$simex_b)) out$simex_b <- as.integer(round(as.numeric(x$simex_b)))   # Round 41
  if (!is.null(x$mv_orient)) out$mv_orient <- as.integer(round(as.numeric(x$mv_orient)))   # Round 50
  if (!is.null(x$mv_iter)) out$mv_iter <- as.integer(round(as.numeric(x$mv_iter)))
  for (k in intersect(c("ml_psi", "mbe_phi", "piv_lambda"), names(x))) out[[k]] <- as.numeric(x[[k]])   # Round 58
  for (k in intersect(c("mbe_iter", "cml_pert", "mv_cml_pert"), names(x))) out[[k]] <- as.integer(round(as.numeric(x[[k]])))
  if (!is.null(x$mv_qhet_iter)) out$mv_qhet_iter <- as.integer(round(as.numeric(x$mv_qhet_iter)))   # Round 54
  for (k in grep("^mv_r_[1-5]_[2-6]$", names(x), value = TRUE)) out[[k]] <- as.numeric(x[[k]])   # Round 53: phenotypic correlations of the exposures
  if (!is.null(x$lso_fraction)) out$lso_fraction <- as.numeric(x$lso_fraction)   # Round 43
  if (!is.null(x$lso_repeats)) out$lso_repeats <- as.integer(round(as.numeric(x$lso_repeats)))
  if (!is.null(x$r_xxo)) out$r_xxo <- as.numeric(x$r_xxo)
  if (!is.null(x$r_yyo)) out$r_yyo <- as.numeric(x$r_yyo)
  # Round 38: every enum setting a method's schema declares is passed on, so a new method's choices are never dropped.
  for (a in param_adapters()) for (k in names(a$params))
    if (is.null(out[[k]]) && !is.null(x[[k]]) && identical(a$params[[k]]$type, "enum") && as.character(x[[k]]) %in% a$params[[k]]$values) out[[k]] <- as.character(x[[k]])
  out
}

param_schemas <- function() {
  one <- function(ad) lapply(ad$params, function(p) list(label = p$label, hint = p$hint, min = p$min, max = p$max, default = p$default, type = p$type,
                                                         values = if (!is.null(p$values)) as.list(p$values), value_labels = if (!is.null(p$value_labels)) as.list(p$value_labels), value_hints = if (!is.null(p$value_hints)) as.list(p$value_hints)))
  list(ivw_mre = one(ADAPTER_IVW_MRE_V11), egger = one(ADAPTER_EGGER_V11), divw = one(ADAPTER_DIVW_V11), raps = one(ADAPTER_RAPS_V11), conmix = one(ADAPTER_CONMIX_V11), lasso = one(ADAPTER_LASSO_V11), egger_simex = one(ADAPTER_EGGER_SIMEX_V11), wm = one(ADAPTER_WM_V11), smed = one(ADAPTER_SMED_V11), mode = one(ADAPTER_MODE_V11), smode = one(ADAPTER_SMODE_V11), presso = one(ADAPTER_PRESSO_V11), radial = one(ADAPTER_RADIAL_V11), radegger = one(ADAPTER_RADEGGER_V11), steiger = one(ADAPTER_STEIGER_V11), loo = one(ADAPTER_LOO_V11))
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
      ev(t, name = name, dir = h$dir, root = root, plan = plan, pairs = pairs, sens = sens, mvcov = ui_state_get(h$con, "mv_cov") %||% list(), outliers = ui_state_get(h$con, "outliers") %||% list(), payloads = payloads, overlap = overlap_get(h), to = to, studies = studies_overview(root))
    }

    # Round 55: individual-level file for the MVMR covariances. Only its header is read here; the covariances are computed on request.
    observeEvent(input$mvcovfile, {
      f <- input$mvcovfile
      tryCatch({
        hd <- data.table::fread(f$datapath, nrows = 50L, data.table = FALSE, showProgress = FALSE)
        rv$mvcov <- list(path = f$datapath, name = f$name)
        ev("mvcov_head", name = f$name, columns = as.list(names(hd)[vapply(hd, is.numeric, TRUE)]))
      }, error = function(e) { rv$mvcov <- NULL; ev("mvcov_error", msg = conditionMessage(e)) })
    })

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
          ev("init", settings = s, library = library_list(), templates = template_list(), studies = studies_overview(root), catalogue = catalogue_all(), schemas = param_schemas(), root = root, diag = run_diagnostics(root), ld = ld_panel_status())
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
          pt <- as.numeric(a$pthr %||% 5e-8); f <- lapply(a$pairs, function(p) tryCatch({ mi <- if (length(p) > 2 && !is.null(p[[3]])) mediation_pair_input(p[[3]], pt) else NULL   # Round 57: mediation steps
            if (length(unlist(p[[1]])) > 1) mv_pair_facts(p[[1]], p[[2]], pt, mi) else pair_facts(p[[1]], p[[2]], pt, mi) }, error = function(e) list(error = conditionMessage(e))))   # Round 50: MVMR pairs
          rv$facts <- Filter(function(x) is.null(x$error), f)
          ev("facts", facts = f, feas = feasibility(Filter(function(x) !isTRUE(x$mv), rv$facts), run_params(a$params %||% list()), if (!is.null(rv$h)) overlap_get(rv$h) else list()))
        },
        mediation = {   # Round 57: write and run the mediation script over the study's mediation pairs
          req(rv$h); res <- mediation_run(rv$h, ui_state_get(rv$h$con, "pairs") %||% list(), as.numeric(a$ci %||% 0.95))
          ev("mediation", rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), spec = attr(res, "spec"))   # Shiny sends a data frame by columns
        },
        mvcov_compute = {   # Round 55: covariances from the chosen individual-level file, for these exposures in model order
          req(rv$h); if (is.null(rv$mvcov)) stop("Choose the individual-level file first.", call. = FALSE)
          cov <- mv_cov_from_file(rv$mvcov$path, unlist(a$cols))
          st <- mv_cov_save(rv$h, a$ex, cov, rv$mvcov$name); rv$mvcov <- NULL
          ev("mvcov_saved", mvcov = st)
        },
        mvcov_clear = { req(rv$h); ev("mvcov_saved", mvcov = mv_cov_clear(rv$h, a$ex)) },
        run_pair = {
          req(rv$h)
          pt <- as.numeric(a$pthr %||% 5e-8)
          mvp <- length(unlist(a$ex)) > 1   # Round 50: several exposures in one model
          inp <- tryCatch(if (!is.null(a$med)) mediation_pair_input(a$med, pt) else if (mvp) mv_pair_input(a$ex, a$out, pt) else pair_input(a$ex, a$out, pt), error = function(e) e)   # Round 57: mediation steps
          if (mvp && !inherits(inp, "error")) inp <- tryCatch({ inp$input <- mv_with_cov(rv$h, a$ex, inp$input); inp }, error = function(e) e)   # Round 55
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
          tryCatch(save_instrument_selection(rv$h$dir, r$pair_folder, inp$selection), error = function(e) NULL)
          tryCatch(save_eaf_check(rv$h$dir, r$pair_folder, inp$eaf_rows), error = function(e) NULL)   # Round 33
          pl <- ui_state_get(rv$h$con, "pairs") %||% list()
          pl <- Filter(function(p) !(identical(p$ex, a$ex) && identical(p$out, a$out)), pl)
          pl[[length(pl) + 1L]] <- list(ex = a$ex, out = a$out, rev = isTRUE(a$rev), mv = mvp, med = a$med, exNames = if (mvp) as.list(inp$exposure_labels) else NULL, exName = inp$exposure_label, outName = inp$outcome_label, folder = r$pair_folder, methods = as.list(meths), ci = ci, params = as.list(prm), prev = list(x = a$prev$x, y = a$prev$y), pthr = pt, simulated = isTRUE(inp$simulated), ledger = inp$ledger, anc_x = inp$anc_x, anc_y = inp$anc_y)
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
          if (identical(a$kind, "eaf")) rec$kind <- "eaf"   # Round 34: without the variants flagged by the allele-frequency check
          sens[[length(sens) + 1L]] <- rec; ui_state_set(rv$h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(rv$h$dir, r$pair_folder), reused = isTRUE(r$reused))
        },
        run_outliers = {   # Round 42: outlier module (FR-90 to FR-102). Every eligible set is analysed as its own labelled pair; the primary stays as it is.
          req(rv$h)
          parent <- a$folder; pay <- pair_payload(rv$h$dir, parent)
          sig <- outlier_signals(pay); trig <- outlier_triggered(sig); sets <- outlier_sets(pay)
          if (!length(attr(sets, "detectors_run"))) { ev("outliers_done", ok = FALSE, parent = parent, msg = "No outlier detector was run for this pair. Add MR-PRESSO, radial IVW, radial MR-Egger or MR-Lasso on the Methods screen and run again."); return() }
          prm <- run_params(a$params); sens <- ui_state_get(rv$h$con, "sens") %||% list(); pays <- list(); done <- list()
          for (i in which(sets$eligible)) {
            rem <- sort(sets$ids[[i]]); key <- paste(rem, collapse = ";")
            if (!is.null(done[[key]])) { pays[[sets$set[i]]] <- done[[key]]$pay; next }
            old <- Filter(function(x) identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem), sens)
            idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
            r <- tryCatch(studio_run_without(rv$h, a$ex, a$out, rem, OUTLIER_REANALYSIS, prm$ci_level, idx, params = prm, pthr = as.numeric(a$pthr %||% 5e-8), keep_presso = "presso" %in% trig, keep_loo = any(vapply(pay$diagnostics, function(d) identical(d$method, "loo"), TRUE))), error = function(e) e)
            if (inherits(r, "error")) { ev("outliers_done", ok = FALSE, parent = parent, msg = paste0(sets$label[i], ": ", conditionMessage(r))); return() }
            bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
            if (length(bad)) { ev("outliers_done", ok = FALSE, parent = parent, msg = paste0(sets$label[i], ": some methods did not finish: ", paste(bad, collapse = ", "))); return() }
            sens <- Filter(function(x) !(identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem)), sens)
            rec <- list(parent = parent, remove = as.list(rem), folder = r$pair_folder, label = r$label, kept = r$kept, kind = "outlier", sets = as.list(sets$label[vapply(sets$ids, function(x) identical(sort(x), rem), TRUE) & sets$eligible]))
            sens[[length(sens) + 1L]] <- rec; done[[key]] <- list(pay = pair_payload(rv$h$dir, r$pair_folder), rec = rec); pays[[sets$set[i]]] <- done[[key]]$pay
          }
          ui_state_set(rv$h$con, "sens", sens)
          scen <- outlier_scenarios(pay, sets, pays); sel <- outlier_select(scen, trig)
          res <- outlier_record(parent, sig, trig, scen, sel, NULL)
          flds <- list(); for (d in done) for (k in sets$set[vapply(sets$ids, function(x) identical(sort(x), sort(unlist(d$rec$remove))), TRUE) & sets$eligible]) flds[[k]] <- d$rec$folder
          tryCatch(outlier_write(rv$h$dir, parent, pay, sets, scen, sel, folders = flds), error = function(e) NULL)
          ol <- Filter(function(x) !identical(x$parent, parent), ui_state_get(rv$h$con, "outliers") %||% list()); ol[[length(ol) + 1L]] <- res; ui_state_set(rv$h$con, "outliers", ol)
          ev("outliers_done", ok = TRUE, parent = parent, res = res, recs = unname(lapply(done, `[[`, "rec")), payloads = stats::setNames(lapply(done, `[[`, "pay"), vapply(done, function(d) d$rec$folder, "")))
        },
        outlier_lookup = {   # Round 49 (FR-106, FR-112): nearest gene (Ensembl) and other-trait associations (OpenGWAS) of the flagged variants
          req(rv$h)
          ol <- ui_state_get(rv$h$con, "outliers") %||% list(); j <- which(vapply(ol, function(x) identical(x$parent, a$folder), TRUE))
          if (!length(j)) return(invisible())
          own <- unlist(lapply(c(a$ex, a$out), function(id) if (startsWith(id %||% "", "lib:")) tryCatch(library_dataset(sub("^lib:", "", id))$recipe$source$id, error = function(e) NULL)))   # the pair's own GWAS are not "other" traits
          lk <- tryCatch(outlier_lookup(ol[[j[1]]], ogw_token_get(), exclude = own), error = function(e) e)
          if (inherits(lk, "error")) { ev("outlier_lookup_done", ok = FALSE, parent = a$folder, msg = conditionMessage(lk)); return() }
          ol[[j[1]]]$lookup <- list(rows = unname(lapply(seq_len(nrow(lk$rows)), function(i) as.list(lk$rows[i, , drop = FALSE]))), notes = as.list(lk$notes), date = lk$date, skipped = as.list(lk$skipped), sources = as.list(lk$sources))
          ui_state_set(rv$h$con, "outliers", ol); tryCatch(outlier_write_lookup(rv$h$dir, a$folder, lk), error = function(e) NULL)
          ev("outlier_lookup_done", ok = TRUE, parent = a$folder, res = ol[[j[1]]])
        },
        outlier_choose = {   # FR-102: the analyst adopts another eligible set, or keeps all instruments; the choice is recorded beside the rule's
          req(rv$h)
          ol <- ui_state_get(rv$h$con, "outliers") %||% list(); j <- which(vapply(ol, function(x) identical(x$parent, a$folder), TRUE))
          if (!length(j)) return(invisible())
          ol[[j[1]]]$chosen <- a$set; ol[[j[1]]]$review <- "reviewed"; ui_state_set(rv$h$con, "outliers", ol)
          tryCatch(outlier_write_choice(rv$h$dir, ol[[j[1]]]), error = function(e) NULL)
          ev("outliers_chosen", parent = a$folder, res = ol[[j[1]]])
        },
        run_ties = {   # Round 28: tie sensitivity analysis (the strongest |beta/SE| as lead where equal P decided it)
          req(rv$h); parent <- a$folder
          sens <- ui_state_get(rv$h$con, "sens") %||% list()
          old <- Filter(function(x) identical(x$parent, parent) && identical(x$kind, "ties"), sens)
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(rv$h$dir)) + 1L
          r <- tryCatch(studio_run_ties(rv$h, a$ex, a$out, unlist(a$methods), run_params(a$params)$ci_level, idx, params = run_params(a$params), pthr = as.numeric(a$pthr %||% 5e-8)), error = function(e) e)
          if (inherits(r, "error")) { ev("sens_done", ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("sens_done", ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "))); return() }
          tryCatch(save_instrument_selection(rv$h$dir, r$pair_folder, r$selection), error = function(e) NULL)
          sens <- Filter(function(x) !(identical(x$parent, parent) && identical(x$kind, "ties")), sens)
          rec <- list(parent = parent, kind = "ties", remove = list(), folder = r$pair_folder, label = r$label, kept = r$kept, regions = r$ties$regions, swaps = list(out = r$ties$out, `in` = r$ties[["in"]]))
          sens[[length(sens) + 1L]] <- rec; ui_state_set(rv$h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(rv$h$dir, r$pair_folder), reused = isTRUE(r$reused))
        },
        build_package = {
          req(rv$h)
          r <- tryCatch(build_package(rv$h, unlist(a$items), a$width_mm, a$pub), error = function(e) e)
          if (inherits(r, "error")) { ev("package_done", ok = FALSE, msg = conditionMessage(r)); return() }
          ev("package_done", ok = TRUE, folder = r$folder, files = as.list(r$files), notes = as.list(r$notes), simulated = r$simulated, strobe = r$strobe)
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
        feas = { req(rv$facts); ev("feas", feas = feasibility(Filter(function(x) !isTRUE(x$mv), rv$facts), run_params(a$params %||% list()), if (!is.null(rv$h)) overlap_get(rv$h) else list())) },
        ld_status = ev("ld", status = ld_panel_status()),
        ld_download = { r <- tryCatch(ld_panel_download_start(unlist(a$pops) %||% "EUR"), error = function(e) e); if (inherits(r, "error")) ev("error", msg = conditionMessage(r)); ev("ld", status = ld_panel_status()) },
        ld_find = ev("ld_found", items = ld_archive_find()),
        ld_use_file = { r <- tryCatch(ld_panel_download_start(unlist(a$pops) %||% "EUR", local_file = a$path), error = function(e) e)
          if (inherits(r, "error")) ev("ld_file_err", msg = conditionMessage(r)); ev("ld", status = ld_panel_status()) },
        ld_find_unpacked = ev("ld_found_unpacked", items = ld_unpacked_find()),
        ld_copy = { r <- tryCatch(ld_panel_copy_start(a$path, if (nzchar(a$pop %||% "")) a$pop else NULL), error = function(e) e)
          if (inherits(r, "error")) ev("ld_copy_err", msg = conditionMessage(r)); ev("ld", status = ld_panel_status()) },
        ld_remove = { ld_panel_remove(as.character(a$pop)); ev("ld", status = ld_panel_status()) },
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
          r <- tryCatch(ogw_get_exposure(a$id, ogw_token_get(), as.numeric(a$p), as.numeric(a$r2), as.numeric(a$kb), a$pop, a$name, isTRUE(settings_read()$offline), where = if (isTRUE(a$where %in% c("local", "server_live"))) a$where else "server"), error = function(e) e)
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
