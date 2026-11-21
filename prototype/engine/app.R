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
      div(style = "position:fixed;right:10px;bottom:6px;font-size:10px;color:#8a97a6;opacity:0;pointer-events:none", "MR GUI engine prototype. Engine: ", textOutput("hello", inline = TRUE)),   # Round 91: kept for the engine checks, no longer drawn over the content
      tags$script(HTML(js))))
}

# The Customise panel is generated from the adapters' own parameter schemas (FR-62), so a label or hint is written once.
#' Settings sent by the screen, checked against the adapters before anything runs (unknown names are dropped, bad values are refused later by the adapter).
#' Round 81: what an LD-based method needs for a pair (panel, LD scores, the analyst's acceptance for Mixed data); NULL without full GWAS.
.pair_ldref <- function(p, method) tryCatch({ ids <- full_gwas_ids(c(p$ex, p$out)); if (!all(startsWith(ids, "lib:"))) return(NULL)
  recs <- lapply(ids, function(i) .lib_recipe(sub("^lib:", "", i))); if (any(vapply(recs, is.null, TRUE))) NULL else ldref_info(recs, method) }, error = function(e) NULL)
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
  for (k in intersect(c("pwm_penk", "smn_phi", "wmn_phi"), names(x))) out[[k]] <- as.numeric(x[[k]])   # Round 59
  for (k in intersect(c("eb_nboot", "pwm_nboot", "smn_nboot", "wmn_nboot", "rk_nboot"), names(x))) out[[k]] <- as.integer(round(as.numeric(x[[k]])))
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
  # Round 72e: every method the engine knows (the six MVMR methods were missing, so their Customise panels were empty)
  ad <- param_adapters(); s <- lapply(ad, one); names(s) <- names(ad); s
}

# Round 83: the long analyses, as functions of the study (h), the screen's arguments (a) and the message sender (ev). The engine runs
# them in a background process (bg_core.R), so it stays free for the screen; the code is what the handlers ran before.
pair_inp_h <- function(h, a) {   # Round 75g: the input of a saved mediation step, rebuilt from its own record
  rec <- if (!is.null(h) && !is.null(a$folder)) Filter(function(p) identical(p$folder, a$folder), ui_state_get(h$con, "pairs") %||% list()) else list()
  med <- if (length(rec)) rec[[1]]$med else NULL
  if (is.null(med)) return(NULL)
  mediation_pair_input(med, as.numeric(a$pthr %||% rec[[1]]$pthr %||% 5e-8))
}
BG_ACTIONS$run_sens <- function(h, a, ev) { t <- "run_sens"
          req(h)
          rem <- sort(unlist(a$remove)); parent <- a$folder
          if (!length(rem)) stop("Choose at least one variant to leave out.")
          sens <- ui_state_get(h$con, "sens") %||% list()
          old <- Filter(function(x) identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem), sens)
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(h$dir)) + 1L
          r <- tryCatch(studio_run_without(h, a$ex, a$out, rem, unlist(a$methods), run_params(a$params)$ci_level, idx, params = run_params(a$params), pthr = as.numeric(a$pthr %||% 5e-8), inp = pair_inp_h(h, a)), error = function(e) e)
          if (inherits(r, "error")) { ev("sens_done", ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("sens_done", ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "))); return() }
          sens <- Filter(function(x) !(identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem)), sens)
          rec <- list(parent = parent, remove = as.list(rem), folder = r$pair_folder, label = r$label, kept = r$kept)
          if (identical(a$kind, "eaf")) rec$kind <- "eaf"
          if (!is.null(a$why)) rec$why <- as.character(a$why)[1]   # Round 72e: why these variants were left out (MVMR review)   # Round 34: without the variants flagged by the allele-frequency check
          sens[[length(sens) + 1L]] <- rec; ui_state_set(h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(h$dir, r$pair_folder), reused = isTRUE(r$reused))
        }

BG_ACTIONS$run_outliers <- function(h, a, ev) { t <- "run_outliers"   # Round 42: outlier module (FR-90 to FR-102). Every eligible set is analysed as its own labelled pair; the primary stays as it is.
          req(h)
          parent <- a$folder; pay <- pair_payload(h$dir, parent)
          sig <- outlier_signals(pay); trig <- outlier_triggered(sig); sets <- outlier_sets(pay)
          if (!length(attr(sets, "detectors_run"))) { ev("outliers_done", ok = FALSE, parent = parent, msg = "No outlier detector was run for this pair. Add MR-PRESSO, radial IVW, radial MR-Egger or MR-Lasso on the Methods screen and run again."); return() }
          prm <- run_params(a$params); sens <- ui_state_get(h$con, "sens") %||% list(); pays <- list(); done <- list()
          pin <- pair_inp_h(h, a)   # Round 75g: a mediation step's own datasets
          for (i in which(sets$eligible)) {
            rem <- sort(sets$ids[[i]]); key <- paste(rem, collapse = ";")
            if (!is.null(done[[key]])) { pays[[sets$set[i]]] <- done[[key]]$pay; next }
            old <- Filter(function(x) identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem), sens)
            idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(h$dir)) + 1L
            r <- tryCatch(studio_run_without(h, a$ex, a$out, rem, OUTLIER_REANALYSIS, prm$ci_level, idx, params = prm, pthr = as.numeric(a$pthr %||% 5e-8), keep_presso = "presso" %in% trig, keep_loo = any(vapply(pay$diagnostics, function(d) identical(d$method, "loo"), TRUE)), inp = pin), error = function(e) e)
            if (inherits(r, "error")) { ev("outliers_done", ok = FALSE, parent = parent, msg = paste0(sets$label[i], ": ", conditionMessage(r))); return() }
            bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
            if (length(bad)) { ev("outliers_done", ok = FALSE, parent = parent, msg = paste0(sets$label[i], ": some methods did not finish: ", paste(bad, collapse = ", "))); return() }
            sens <- Filter(function(x) !(identical(x$parent, parent) && identical(sort(unlist(x$remove)), rem)), sens)
            rec <- list(parent = parent, remove = as.list(rem), folder = r$pair_folder, label = r$label, kept = r$kept, kind = "outlier", sets = as.list(sets$label[vapply(sets$ids, function(x) identical(sort(x), rem), TRUE) & sets$eligible]))
            sens[[length(sens) + 1L]] <- rec; done[[key]] <- list(pay = pair_payload(h$dir, r$pair_folder), rec = rec); pays[[sets$set[i]]] <- done[[key]]$pay
          }
          ui_state_set(h$con, "sens", sens)
          scen <- outlier_scenarios(pay, sets, pays); sel <- outlier_select(scen, trig)
          res <- outlier_record(parent, sig, trig, scen, sel, NULL)
          flds <- list(); for (d in done) for (k in sets$set[vapply(sets$ids, function(x) identical(sort(x), sort(unlist(d$rec$remove))), TRUE) & sets$eligible]) flds[[k]] <- d$rec$folder
          tryCatch(outlier_write(h$dir, parent, pay, sets, scen, sel, folders = flds), error = function(e) NULL)
          ol <- Filter(function(x) !identical(x$parent, parent), ui_state_get(h$con, "outliers") %||% list()); ol[[length(ol) + 1L]] <- res; ui_state_set(h$con, "outliers", ol)
          ev("outliers_done", ok = TRUE, parent = parent, res = res, recs = unname(lapply(done, `[[`, "rec")), payloads = stats::setNames(lapply(done, `[[`, "pay"), vapply(done, function(d) d$rec$folder, "")))
        }

BG_ACTIONS$run_ties <- function(h, a, ev) { t <- "run_ties"   # Round 28: tie sensitivity analysis (the strongest |beta/SE| as lead where equal P decided it)
          req(h); parent <- a$folder
          sens <- ui_state_get(h$con, "sens") %||% list()
          old <- Filter(function(x) identical(x$parent, parent) && identical(x$kind, "ties"), sens)
          idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(h$dir)) + 1L
          r <- tryCatch(studio_run_ties(h, a$ex, a$out, unlist(a$methods), run_params(a$params)$ci_level, idx, params = run_params(a$params), pthr = as.numeric(a$pthr %||% 5e-8)), error = function(e) e)
          if (inherits(r, "error")) { ev("sens_done", ok = FALSE, msg = conditionMessage(r)); return() }
          bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
          if (length(bad)) { ev("sens_done", ok = FALSE, msg = paste("Some methods did not finish:", paste(bad, collapse = ", "))); return() }
          tryCatch(save_instrument_selection(h$dir, r$pair_folder, r$selection), error = function(e) NULL)
          sens <- Filter(function(x) !(identical(x$parent, parent) && identical(x$kind, "ties")), sens)
          rec <- list(parent = parent, kind = "ties", remove = list(), folder = r$pair_folder, label = r$label, kept = r$kept, regions = r$ties$regions, swaps = list(out = r$ties$out, `in` = r$ties[["in"]]))
          sens[[length(sens) + 1L]] <- rec; ui_state_set(h$con, "sens", sens)
          ev("sens_done", ok = TRUE, rec = rec, payload = pair_payload(h$dir, r$pair_folder), reused = isTRUE(r$reused))
        }

BG_ACTIONS$mediation <- function(h, a, ev) { t <- "mediation"   # Round 57: write and run the mediation script over the study's mediation pairs
          req(h); .say("Mediation: combining the steps into the total, direct and indirect effects and the proportion mediated, with delta-method intervals", "The total effect comes from univariable MR, the direct effect from multivariable MR; the indirect effect is their product (exposure on mediator times mediator on outcome) and also total minus direct.")
          res <- mediation_run(h, ui_state_get(h$con, "pairs") %||% list(), as.numeric(a$ci %||% 0.95))
          ev("mediation", rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), spec = attr(res, "spec"))   # Shiny sends a data frame by columns
        }

BG_ACTIONS$bidir_cml <- function(h, a, ev) { t <- "bidir_cml"   # Round 63: bidirectional constrained maximum likelihood over each forward/reverse couple
          req(h); if (!requireNamespace("BiDirectCausal", quietly = TRUE)) { ev("bidir_cml", rows = list(), msg = "BiDirectCausal is not installed in this runtime."); return() }
          .say("Bidirectional MR-cML: estimating both directions together for each pair of traits", "Fitting both directions jointly allows for the reverse effect when estimating each one.")
          res <- bidir_cml_run(h, ui_state_get(h$con, "pairs") %||% list())
          ev("bidir_cml", rows = unlist(lapply(names(res), function(k) { r <- res[[k]]; if (!is.data.frame(r)) return(list())
            lapply(seq_len(nrow(r)), function(i) c(list(forward = k), as.list(r[i, , drop = FALSE]))) }), recursive = FALSE), errors = Filter(Negate(is.data.frame), res))
        }

BG_ACTIONS$build_package <- function(h, a, ev) { t <- "build_package"
          req(h)
          r <- tryCatch(build_package(h, unlist(a$items), a$width_mm, a$pub), error = function(e) e)
          if (inherits(r, "error")) { ev("package_done", ok = FALSE, msg = conditionMessage(r)); return() }
          ev("package_done", ok = TRUE, folder = r$folder, files = as.list(r$files), notes = as.list(r$notes), simulated = r$simulated, strobe = r$strobe, fonts = .fig_cairo())   # Round 91: whether the figure fonts are embedded
        }

BG_ACTIONS$ldsc <- function(h, a, ev) { t <- "ldsc"
          req(h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(h$con, "pairs") %||% list()); if (is.null(p)) return()
          res <- if (identical(t, "ldsc_get")) ldsc_pair_read(h, p) else {
            if (!requireNamespace("GenomicSEM", quietly = TRUE)) { ev("ldsc", folder = a$folder, rows = list(), msg = "GenomicSEM is not installed in this runtime."); return() }
            prev <- vapply(1:2, function(i) suppressWarnings(as.numeric(a$prev[[i]] %||% NA)), 0)
            tryCatch(ldsc_pair_run(h, p, a$source, prev, accept = a$accept), error = function(e) e) }
          if (inherits(res, "error")) { ev("ldsc", folder = a$folder, rows = list(), msg = conditionMessage(res), anc = .pair_ldref(p, "ldsc"), ldref = ldref_read(h, p, "ldsc")); return() }
          ev("ldsc", folder = a$folder, rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else ldsc_text(res), anc = .pair_ldref(p, "ldsc"), ldref = ldref_read(h, p, "ldsc"))
        }

BG_ACTIONS$run_pair <- function(h, a, ev) { t0 <- Sys.time()
  pair_done <- function(...) { t1 <- Sys.time(); ev("pair_done", ..., started = format(t0, "%H:%M:%S"), finished = format(t1, "%H:%M:%S"), secs = round(as.numeric(difftime(t1, t0, units = "secs")))) }
      req(h)
      pt <- as.numeric(a$pthr %||% 5e-8)
      mvp <- length(unlist(a$ex)) > 1   # Round 50: several exposures in one model
      inp <- tryCatch(if (!is.null(a$med)) mediation_pair_input(a$med, pt) else if (mvp) mv_pair_input(a$ex, a$out, pt) else pair_input(a$ex, a$out, pt), error = function(e) e)   # Round 57: mediation steps
      if (mvp && !inherits(inp, "error")) inp <- tryCatch({ inp$input <- mv_with_cov(h, a$ex, inp$input); inp }, error = function(e) e)   # Round 55
      if (inherits(inp, "error")) { pair_done(i = a$i, ok = FALSE, msg = conditionMessage(inp)); return() }
      meths <- unlist(a$methods); prm <- run_params(a$params); ci <- prm$ci_level
      live_reset(); .say(sprintf("Pair %s: %s → %s, %d instruments; preparing the analysis data and writing one R script per method", a$i + 1L, inp$exposure_label, inp$outcome_label, sum(inp$input$data$mr_keep %in% TRUE)), "The scripts and the data are saved in the pair folder, so every result can be rerun in R and checked.")
      # The contract slugs labels to 30 characters, so two different traits can share a folder name. The pair number keeps them apart.
      old <- Filter(function(p) identical(p$ex, a$ex) && identical(p$out, a$out), ui_state_get(h$con, "pairs") %||% list())
      idx <- if (length(old)) as.integer(sub("^P([0-9]+)_.*", "\\1", old[[1]]$folder)) else length(studio_pair_folders(h$dir)) + 1L
      sti <- if (any(c("steiger", "ivw_steiger") %in% meths)) steiger_input(inp, list(x = a$prev$x, y = a$prev$y)) else NULL   # prevalence the analyst entered for a binary trait
      r <- tryCatch(studio_run(h, inp$exposure_label, inp$outcome_label, inp$input, meths, ci, idx = idx, params = prm, steiger_in = sti), error = function(e) e)
      if (inherits(r, "error")) { pair_done(i = a$i, ok = FALSE, msg = conditionMessage(r)); return() }
      bad <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
      # Round 70c: each method runs on its own; the pair fails only when no method finished, otherwise it is saved and the failures are named
      if (length(bad) && length(bad) == length(r$results)) { pair_done(i = a$i, ok = FALSE, msg = paste("No method finished:", paste(bad, collapse = ", "), paste(unlist(lapply(r$results[bad], `[[`, "messages")), collapse = " "))); return() }
      if (length(bad)) r$left_out <- c(r$left_out, stats::setNames(vapply(r$results[bad], function(x) paste0("stopped with an error (", gsub("^FAILED [^:]+: ", "", paste(x$messages, collapse = " ")), ")"), ""), bad))
      tryCatch(save_instrument_selection(h$dir, r$pair_folder, inp$selection), error = function(e) NULL)
      tryCatch(save_eaf_check(h$dir, r$pair_folder, inp$eaf_rows), error = function(e) NULL)   # Round 33
      pl <- ui_state_get(h$con, "pairs") %||% list()
      pl <- Filter(function(p) !(identical(p$ex, a$ex) && identical(p$out, a$out)), pl)
      pl[[length(pl) + 1L]] <- list(ex = a$ex, out = a$out, rev = isTRUE(a$rev), mv = mvp, med = a$med, exNames = if (mvp) as.list(inp$exposure_labels) else NULL, exName = inp$exposure_label, outName = inp$outcome_label, folder = r$pair_folder, methods = as.list(meths), ci = ci, params = as.list(prm), prev = list(x = a$prev$x, y = a$prev$y), pthr = pt, simulated = isTRUE(inp$simulated), ledger = inp$ledger, anc_x = inp$anc_x, anc_y = inp$anc_y)
      ui_state_set(h$con, "pairs", pl)
      rp <- if (isTRUE(r$reused)) NULL else tryCatch(pair_replay(h$dir, r$pair_folder), error = function(e) NULL)   # Round 91: the run's last line says how the replay check ended
      .say(if (isTRUE(r$reused)) "These results were already saved with the same data and settings, so they are reused" else if (is.null(rp)) "Saved. The replay check did not run" else
           switch(as.character(rp$status %||% ""), identical = sprintf("Replay check finished: all %s numbers identical to the first run. Saved.", format(rp$n_values_compared %||% NA, big.mark = ",")),
                  sprintf("Replay check finished: %s (largest relative difference %s; see 07_replay_check). Saved.", gsub("_", " ", rp$status %||% "unknown"), format(signif(rp$max_rel_diff %||% NA, 2)))))
      pair_done(i = a$i, ok = TRUE, folder = r$pair_folder, reused = isTRUE(r$reused), payload = pair_payload(h$dir, r$pair_folder), left_out = if (length(r$left_out)) paste(sprintf("%s: %s.", names(r$left_out), r$left_out), collapse = " ") else "")
}

mr_server <- function(engine_dir) {
  function(input, output, session) {
    output$hello <- renderText("connected")
    rv <- reactiveValues(h = NULL, imp = NULL)
    ev <- function(t, ...) session$sendCustomMessage("ev", c(list(t = t), list(...)))
    options(mrstudio.say = function(msg, why = "") ev("live", msg = msg, why = why, act = getOption("mrstudio.act", "")))   # Round 75d: live narration of the background work; Round 80: each message names its action (its own box on the screen)
    options(mrstudio.script_progress = function(dir, console_name) live_scripts(dir, console_name))   # Round 75e: which method is running, on the Run screen
    # Round 75f: one pair (what the run_pair event did), and the queue the engine works through by itself. Between pairs the queue returns
    # to Shiny's event loop (later::later), so Pause and Stop from the window are handled; the run no longer needs the window to ask for
    # each pair, so it continues when the screen is locked or the window hidden. At the end the computer may sleep again and the
    # notification is shown by the engine.
    # Round 75g: the input of a saved pair, rebuilt from its own record. A mediation step uses its own datasets (for exposure -> mediator,
    # the mediator looked up at the exposure's instruments), which its plain exposure and outcome labels do not say; instrument review,
    # sensitivity reruns and the outlier analysis failed on such pairs ("TG was imported as an exposure only").
    pair_inp <- function(a) {
      rec <- if (!is.null(rv$h) && !is.null(a$folder)) Filter(function(p) identical(p$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()) else list()
      med <- if (length(rec)) rec[[1]]$med else NULL
      if (is.null(med)) return(NULL)
      mediation_pair_input(med, as.numeric(a$pthr %||% rec[[1]]$pthr %||% 5e-8))
    }
    # Round 80: a finished pair carries its start, finish and duration as the engine measured them, so the screen shows the right times
    # even when the window was asleep (a locked screen) and receives the news late
    pair_done <- function(...) { t0 <- rq$t_start %||% Sys.time(); t1 <- Sys.time()
      ev("pair_done", ..., started = format(t0, "%H:%M:%S"), finished = format(t1, "%H:%M:%S"), secs = round(as.numeric(difftime(t1, t0, units = "secs")))) }
    # Round 83: one pair is run by BG_ACTIONS$run_pair in a background process (bg_core.R); the queue below starts it and passes on its messages
    rq <- new.env()
    run_queue_step <- function() {
      if (is.null(rq$pairs)) return(invisible())
      if (isTRUE(rq$pause) && !isTRUE(rq$stop)) return(invisible())   # continued by run_pause(on = FALSE)
      if (isTRUE(rq$stop) || rq$k > length(rq$pairs)) {
        stopped <- isTRUE(rq$stop) && rq$k <= length(rq$pairs); rq$pairs <- NULL
        keep_awake_stop()
        if (isTRUE(settings_read()$notify)) notify_done("MR Studio", sprintf("%s %s: %d of %d pair%s saved.", if (nzchar(rq$study)) rq$study else "The analysis", if (stopped) "stopped" else "finished", rq$ok + rq$done_before, rq$total, if (rq$total == 1) "" else "s"))
        evs(rq$name, "run_queue_done", stopped = stopped); return(invisible()) }
      a <- rq$pairs[[rq$k]]; rq$k <- rq$k + 1L; rq$t_start <- Sys.time()
      # Round 83: the pair runs in a background process; its messages reach the screen as it sends them, and the next pair starts when it ends
      job <- tryCatch(bg_start(rq$dir, "run_pair", a, act = "run_all"), error = function(e) e)
      if (inherits(job, "error")) { evs(rq$name, "pair_done", i = a$i, ok = FALSE, msg = conditionMessage(job)); later::later(run_queue_step, 0.05); return(invisible()) }
      bg_watch(job, rq$name, on_line = function(ln) if (startsWith(ln, '{"t":"pair_done"') && isTRUE(jsonlite::parse_json(ln)$ok)) rq$ok <- rq$ok + 1L,
               on_done = function(r) later::later(run_queue_step, 0.05))
      invisible()
    }
    # Round 83: background jobs. Their messages are passed on every half second, tagged with their study, so the screen files them under
    # that study even when another study is open.
    bg <- new.env(); bg$jobs <- list()
    evs <- function(study, t, ...) tryCatch(session$sendCustomMessage("ev", c(list(t = t), list(...), list(`_study` = study))), error = function(e) NULL)
    bg_tick <- function() {
      for (job in names(bg$jobs)) { j <- bg$jobs[[job]]
        r <- tryCatch(bg_poll(job, j$at), error = function(e) list(lines = character(), at = j$at, done = TRUE))
        for (ln in r$lines) { tryCatch(session$sendCustomMessage("evraw", list(s = j$study, j = ln)), error = function(e) NULL); if (is.function(j$on_line)) tryCatch(j$on_line(ln), error = function(e) NULL) }
        if (isTRUE(r$done)) { bg$jobs[[job]] <- NULL; bg_studies(); if (is.function(j$on_done)) j$on_done(r) } else bg$jobs[[job]]$at <- r$at }
      if (length(bg$jobs)) later::later(bg_tick, BG_POLL_S); invisible()
    }
    bg_watch <- function(job, study, on_line = NULL, on_done = NULL) {
      first <- !length(bg$jobs); bg$jobs[[job]] <- list(at = 0, study = study, on_line = on_line, on_done = on_done)
      if (first) later::later(bg_tick, BG_POLL_S); bg_studies(); invisible(job)
    }
    bg_studies <- function() tryCatch(ev("studies", studies = studies_overview(), root = studies_root()), error = function(e) NULL)   # the list shows RUNNING while a job of a study runs
    bg_go <- function(t, a) { req(rv$h); bg_watch(bg_start(rv$h$dir, t, a, act = t), basename(rv$h$dir)) }
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
        if (!is.null(p$folder) && dir.exists(file.path(h$dir, "03_pairs", p$folder)) && !is.null(pair_results(h$dir, p$folder))) payloads[[p$folder]] <<- pair_payload_json(h$dir, p$folder)
        p
      })
      pairs <- Filter(function(p) !is.null(p$folder), pairs)
      sens <- ui_state_get(h$con, "sens") %||% list()
      for (sn in sens) if (!is.null(sn$folder) && dir.exists(file.path(h$dir, "03_pairs", sn$folder)) && !is.null(pair_results(h$dir, sn$folder))) payloads[[sn$folder]] <- pair_payload_json(h$dir, sn$folder)
      ev(t, name = name, dir = h$dir, root = root, plan = plan, pairs = pairs, sens = sens, mvcov = ui_state_get(h$con, "mv_cov") %||% list(), outliers = ui_state_get(h$con, "outliers") %||% list(), payloads = payloads, overlap = overlap_get(h), to = to, jobs = length(bg_running(h$dir)), queue = !is.null(rq$pairs) && identical(rq$dir, h$dir), studies = studies_overview(root))
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
        if (grepl("\\.vcf(\\.b?gz)?$", f$name, ignore.case = TRUE)) {   # Round 75: a GWAS-VCF is read by its own reader, no column mapping
          pk <- vcf_peek(f$datapath); rv$vcf <- list(path = f$datapath, name = f$name)
          ev("vcf_head", name = f$name, size = file.size(f$datapath), meta = pk$meta, format = as.list(pk$format), samples = as.list(pk$samples), kind = pk$kind, ncase = pk$ncase, ncontrol = pk$ncontrol, variants = pk$variants, build = pk$build)
          return()
        }
        hd <- import_read_head(f$datapath, f$name)
        rv$imp <- list(path = hd$path, orig = f$datapath, name = f$name, sep = hd$sep, columns = hd$columns, fingerprint = hd$fingerprint, notes = unlist(hd$notes))   # Round 77: hd$path may be a cleaned copy (Excel and similar)
        tm <- template_match(hd$fingerprint); sg <- import_suggest(hd$columns, hd$rows)
        ev("import_head", name = hd$name, size = hd$size, sep = hd$sep, columns = as.list(hd$columns), rows = unname(hd$rows),
           suggest = sg, hints = import_hints(hd$columns, hd$rows, sg), template = tm, notes = hd$notes, sheets = hd$sheets, sheet = hd$sheet)
      }, error = function(e) { rv$imp <- NULL; ev("import_error", msg = conditionMessage(e)) })
    })

    observeEvent(input$act, {
      a <- input$act; t <- a$t
      old_act <- options(mrstudio.act = t); on.exit({ options(old_act); if (!identical(t, "run_all") && !(t %in% names(BG_ACTIONS))) ev("live_end", act = t) }, add = TRUE)   # Round 80: the screen closes this action's box
      tryCatch((function() switch(t,
        init = {
          s <- settings_read(); root <- studies_root()
          ev("init", settings = s, library = library_list(), templates = template_list(), studies = studies_overview(root), catalogue = catalogue_all(), schemas = param_schemas(), root = root, diag = run_diagnostics(root), ld = ld_panel_status(), fgjobs = fullgwas_jobs())
          ogw_ev()
        },
        study_folder = open_in_file_manager(file.path(studies_root(), basename(a$name))),
        library = { ev("library", items = library_list(), templates = template_list(), catalogue = catalogue_all()); ogw_ev() },
        library_folder = { d <- library_dir(); dir.create(d, recursive = TRUE, showWarnings = FALSE); open_in_file_manager(d) },
        import_sheet = {   # Round 77: another sheet of the same Excel workbook
          req(rv$imp); hd <- import_read_head(rv$imp$orig, rv$imp$name, sheet = a$sheet)
          rv$imp <- modifyList(rv$imp, list(path = hd$path, sep = hd$sep, columns = hd$columns, fingerprint = hd$fingerprint, notes = unlist(hd$notes)))
          sg <- import_suggest(hd$columns, hd$rows)
          ev("import_head", name = hd$name, size = hd$size, sep = hd$sep, columns = as.list(hd$columns), rows = unname(hd$rows), suggest = sg, hints = import_hints(hd$columns, hd$rows, sg),
             template = template_match(hd$fingerprint), notes = hd$notes, sheets = hd$sheets, sheet = hd$sheet)
        },
        import_check = {
          req(rv$imp); v <- import_validate(rv$imp$path, rv$imp$sep, rv$imp$columns, a$mapping, a$facts)
          ev("import_report", ok = v$ok, errors = as.list(v$errors), warnings = as.list(v$warnings), changes = as.list(c(rv$imp$notes, unlist(v$changes))), stats = v$stats)
        },
        import_save = {
          req(rv$imp)
          rec <- import_save(rv$imp$path, rv$imp$name, rv$imp$sep, rv$imp$columns, a$mapping, a$facts, orig_path = rv$imp$orig %||% rv$imp$path, prep_notes = rv$imp$notes %||% character())
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
          pt <- as.numeric(a$pthr %||% 5e-8); np <- length(a$pairs); ip <- 0L; f <- lapply(a$pairs, function(p) tryCatch({ ip <<- ip + 1L
            .say(sprintf("Pair %d of %d, %s → %s: reading both datasets, keeping the variants below P %s and matching each one in the outcome", ip, np, paste(sub("^lib:", "", unlist(p[[1]])), collapse = " + "), sub("^lib:", "", p[[2]]), format(pt)), "Each instrument's effect on the exposure and on the outcome must refer to the same allele, so the alleles are aligned and ambiguous variants are handled before anything is estimated.")
            mi <- if (length(p) > 2 && !is.null(p[[3]])) mediation_pair_input(p[[3]], pt) else NULL   # Round 57: mediation steps
            if (length(unlist(p[[1]])) > 1) mv_pair_facts(p[[1]], p[[2]], pt, mi) else pair_facts(p[[1]], p[[2]], pt, mi) }, error = function(e) list(error = conditionMessage(e))))   # Round 50: MVMR pairs
          rv$facts <- Filter(function(x) is.null(x$error), f)
          .say("Running the checks: instrument strength (F statistics), ancestry, sample overlap, units, palindromic variants, allele frequencies and the direction of effect", "Each check is shown with what it means; items marked REVIEW do not stop the study and are written into the limitations.")
          ev("facts", facts = f, feas = feasibility(Filter(function(x) !isTRUE(x$mv), rv$facts), run_params(a$params %||% list()), if (!is.null(rv$h)) overlap_get(rv$h) else list()))
        },
        mediation = bg_go("mediation", a),
        bidir_cml = bg_go("bidir_cml", a),
        ldsc_status = ev("ldsc_ref", status = ldsc_status()),   # Round 64: LD scores and LD score regression
        ldsc_fetch = { r <- tryCatch(ldsc_fetch_start(), error = function(e) e);   # Round 82: in the background
          if (inherits(r, "error")) ev("error", msg = conditionMessage(r)); ev("ldsc_ref", status = ldsc_status()) },
        ldsc_compute = { pop <- a$pop %||% tryCatch({ p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); x <- ld_pop_for(library_recipe(sub("^lib:", "", p$ex))$facts$ancestry); if (is.na(x)) "EUR" else x }, error = function(e) "EUR")
          r <- tryCatch(ldsc_compute_start(pop), error = function(e) e); if (inherits(r, "error")) ev("error", msg = conditionMessage(r)); ev("ldsc_ref", status = ldsc_status()) },
        ldsc_get = , ldsc = { if (identical(t, "ldsc")) return(bg_go("ldsc", a))
          req(rv$h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); if (is.null(p)) return()
          res <- if (identical(t, "ldsc_get")) ldsc_pair_read(rv$h, p) else {
            if (!requireNamespace("GenomicSEM", quietly = TRUE)) { ev("ldsc", folder = a$folder, rows = list(), msg = "GenomicSEM is not installed in this runtime."); return() }
            prev <- vapply(1:2, function(i) suppressWarnings(as.numeric(a$prev[[i]] %||% NA)), 0)
            tryCatch(ldsc_pair_run(rv$h, p, a$source, prev, accept = a$accept), error = function(e) e) }
          if (inherits(res, "error")) { ev("ldsc", folder = a$folder, rows = list(), msg = conditionMessage(res), anc = .pair_ldref(p, "ldsc"), ldref = ldref_read(rv$h, p, "ldsc")); return() }
          ev("ldsc", folder = a$folder, rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else ldsc_text(res), anc = .pair_ldref(p, "ldsc"), ldref = ldref_read(rv$h, p, "ldsc"))
        },
        mrapss_get = , mrapss_start = {   # Round 69: MR-APSS, run in a separate process (minutes)
          req(rv$h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); if (is.null(p)) return()
          r <- if (identical(t, "mrapss_start")) tryCatch({ if (!requireNamespace("MRAPSS", quietly = TRUE)) stop("MRAPSS is not installed in this runtime.", call. = FALSE); mrapss_pair_start(rv$h, p, accept = a$accept) }, error = function(e) e) else NULL
          res <- mrapss_pair_read(rv$h, p)
          ev("mrapss", folder = a$folder, job = mrapss_job_read(rv$h, p), anc = .pair_ldref(p, "mrapss"), ldref = ldref_read(rv$h, p, "mrapss"), ldsc_official = ldsc_ref_status()$official, panel = !is.na(plink_bin()) && file.exists(file.path(ld_panel_dir(), "EUR.bed")),
             msg = if (inherits(r, "error")) conditionMessage(r) else "", rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else mrapss_text(res))
        },
        cause_get = , cause_start = {   # Round 70: CAUSE, run in a separate process (minutes)
          req(rv$h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); if (is.null(p)) return()
          r <- if (identical(t, "cause_start")) tryCatch({ if (!requireNamespace("cause", quietly = TRUE)) stop("cause is not installed in this runtime.", call. = FALSE); cause_pair_start(rv$h, p, accept = a$accept) }, error = function(e) e) else NULL
          res <- cause_pair_read(rv$h, p)
          ev("cause", folder = a$folder, job = cause_job_read(rv$h, p), anc = .pair_ldref(p, "cause"), ldref = ldref_read(rv$h, p, "cause"), panel = !is.na(plink_bin()) && file.exists(file.path(ld_panel_dir(), "EUR.bed")),
             msg = if (inherits(r, "error")) conditionMessage(r) else "", rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else cause_text(res), ties = if (is.null(res)) "" else paste(cause_ties_text(res), collapse = " "))   # Round 87
        },
        mrlap_get = , mrlap_start = {   # Round 68: MRlap, run in a separate process (minutes)
          req(rv$h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); if (is.null(p)) return()
          r <- if (identical(t, "mrlap_start")) tryCatch({ if (!requireNamespace("MRlap", quietly = TRUE)) stop("MRlap is not installed in this runtime.", call. = FALSE); mrlap_pair_start(rv$h, p, accept = a$accept) }, error = function(e) e) else NULL
          res <- mrlap_pair_read(rv$h, p)
          ev("mrlap", folder = a$folder, job = mrlap_job_read(rv$h, p), anc = .pair_ldref(p, "mrlap"), ldref = ldref_read(rv$h, p, "mrlap"), ldsc_official = ldsc_ref_status()$official, msg = if (inherits(r, "error")) conditionMessage(r) else "",
             rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else mrlap_text(res))
        },
        lhcmr_get = , lhcmr_fetch = , lhcmr_start = {   # Round 65: LHC-MR, run in a separate process (one to two hours)
          req(rv$h); p <- Find(function(q) identical(q$folder, a$folder), ui_state_get(rv$h$con, "pairs") %||% list()); if (is.null(p)) return()
          r <- tryCatch(switch(t, lhcmr_fetch = lhcmr_fetch_start(), lhcmr_start = { if (!requireNamespace("lhcMR", quietly = TRUE)) stop("lhcMR is not installed in this runtime.", call. = FALSE); w <- suppressWarnings(as.integer(a$workers %||% NA)); lhcmr_pair_start(rv$h, p, cores = if (is.na(w)) NULL else w) }, NULL), error = function(e) e)
          res <- lhcmr_pair_read(rv$h, p)
          op <- tryCatch(lhcmr_pair_options(p), error = function(e) NULL); mi <- if (is.null(op)) mem_info() else attr(op, "memory")   # Round 66: only safe settings
          ev("lhcmr", folder = a$folder, status = lhcmr_status(rv$h, p), anc = .pair_ldref(p, "lhcmr"), memory = mi, options = if (is.null(op) || !nrow(op)) list() else lapply(seq_len(nrow(op)), function(i) as.list(op[i, c("workers", "label", "memory_gb", "hours")])), msg = if (inherits(r, "error")) conditionMessage(r) else "",
             rows = if (is.null(res)) list() else lapply(seq_len(nrow(res)), function(i) as.list(res[i, , drop = FALSE])), text = if (is.null(res)) "" else lhcmr_text(res))
        },
        mvcov_compute = {   # Round 55: covariances from the chosen individual-level file, for these exposures in model order
          req(rv$h); if (is.null(rv$mvcov)) stop("Choose the individual-level file first.", call. = FALSE)
          cov <- mv_cov_from_file(rv$mvcov$path, unlist(a$cols))
          st <- mv_cov_save(rv$h, a$ex, cov, rv$mvcov$name); rv$mvcov <- NULL
          ev("mvcov_saved", mvcov = st)
        },
        mvcov_clear = { req(rv$h); ev("mvcov_saved", mvcov = mv_cov_clear(rv$h, a$ex)) },
        run_pair = bg_go("run_pair", a),   # Round 83: in the background
        run_all = {   # Round 75f: the engine runs the queue of pairs itself (see run_queue_step)
          req(rv$h); if (!is.null(rq$pairs)) stop(sprintf("A run is already going on in \"%s\". Wait until it ends, or stop it there, before starting another.", rq$name))
          rq$dir <- rv$h$dir; rq$name <- basename(rv$h$dir); rq$pairs <- a$pairs; rq$k <- 1L; rq$pause <- FALSE; rq$stop <- FALSE; rq$ok <- 0L; rq$total <- as.integer(a$total %||% length(a$pairs)); rq$study <- a$study %||% ""
          rq$done_before <- rq$total - length(a$pairs); later::later(run_queue_step, 0.05) },
        run_pause = { rq$pause <- isTRUE(a$on); if (!rq$pause && !is.null(rq$pairs)) later::later(run_queue_step, 0.05) },
        run_stop = { rq$stop <- TRUE; if (isTRUE(rq$pause) || is.null(rq$pairs)) { rq$pause <- FALSE; later::later(run_queue_step, 0.05) } },
        instruments = ev("instruments", folder = a$folder, rows = { it <- instrument_table(a$ex, a$out, as.numeric(a$pthr %||% 5e-8), pair_inp(a)); unname(lapply(seq_len(nrow(it)), function(i) as.list(it[i, , drop = FALSE]))) }),
        full_status = { r <- tryCatch(full_gwas_status(a$ex, a$out), error = function(e) list(ok = FALSE, why = "instruments")); ev("full_status", folder = a$folder, ok = isTRUE(r$ok), why = r$why) },   # Round 72e
        mv_review = {   # Round 72e: per-variant contributions to the MVMR-IVW heterogeneity, and MVMR-Lasso's flags
          req(rv$h); r <- tryCatch(mv_review(rv$h, a$folder), error = function(e) e)
          if (inherits(r, "error")) ev("mv_review", folder = a$folder, error = conditionMessage(r))
          else ev("mv_review", folder = a$folder, q_total = r$q_total, df = r$df, n = r$n, lasso = as.list(r$lasso), rows = lapply(seq_len(nrow(r$rows)), function(i) as.list(r$rows[i, ])))   # Round 75c: every variant, to tick any of them
 },
        run_sens = bg_go("run_sens", a),
        run_outliers = bg_go("run_outliers", a),
        outlier_lookup = {   # Round 49 (FR-106, FR-112): nearest gene (Ensembl) and other-trait associations (OpenGWAS) of the flagged variants
          req(rv$h)
          ol <- ui_state_get(rv$h$con, "outliers") %||% list(); j <- which(vapply(ol, function(x) identical(x$parent, a$folder), TRUE))
          if (!length(j)) return(invisible())
          own <- unlist(lapply(c(a$ex, a$out), function(id) if (startsWith(id %||% "", "lib:")) tryCatch(library_recipe(sub("^lib:", "", id))$source$id, error = function(e) NULL)))   # the pair's own GWAS are not "other" traits
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
        run_ties = bg_go("run_ties", a),
        build_package = bg_go("build_package", a),
        open_package = { req(rv$h); open_in_file_manager(file.path(rv$h$dir, "05_publication", a$folder)) },
        open_folder = { req(rv$h); open_in_file_manager(rv$h$dir) },
        med_prepare = {   # Round 72c: the mediators of a mediation study in both forms, and its multivariable sets, from OpenGWAS where needed
          r <- tryCatch(c(mediation_prepare(a$x, unlist(a$ms), a$y, ogw_token_get(), isTRUE(settings_read()$offline)), list(ok = TRUE)), error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
          ev("med_ready", r = r) },
        mv_prepare = {   # Round 72c: multivariable sets from OpenGWAS downloads (no-op for imported full summary statistics)
          rows <- lapply(a$sets, function(q) tryCatch({ ex <- unlist(q[[1]]); out <- q[[2]]
              rs <- lapply(ex, function(e) tryCatch(library_recipe(sub("^lib:", "", e)), error = function(z) NULL))
              if (all(startsWith(ex, "lib:")) && all(vapply(rs, function(r) !is.null(r) && .is_ogw_instruments(r), TRUE))) {
                r <- ogw_mv_build(ex, out, ogw_token_get(), isTRUE(settings_read()$offline)); list(ok = TRUE, note = if (isTRUE(r$reused)) "" else sprintf("prepared \"%s\": %d joint instruments (clumped on %s); every exposure and the outcome were looked up on OpenGWAS at all of them, so the exposure an outcome download was made for does not matter", r$name, r$ledger$final, r$ledger$clumped_where)) }
              else list(ok = TRUE, note = "") }, error = function(e) list(ok = FALSE, msg = conditionMessage(e))))
          ev("mv_ready", rows = rows) },
        out_prepare = {   # Round 88: a two-sample study's outcome looked up again for each exposure it was not downloaded for (in the engine: the token stays in its memory)
          pref <- unique(unlist(lapply(a$pairs, `[[`, 2)))   # the study's own outcomes are preferred to older lookups in the library
          rows <- lapply(a$pairs, function(q) tryCatch(c(pair_outcome_for(q[[1]], q[[2]], ogw_token_get(), isTRUE(settings_read()$offline), prefer = pref), list(ok = TRUE)),
                                                        error = function(e) list(ok = FALSE, msg = conditionMessage(e))))
          ev("out_ready", rows = rows) },
        bidir_prepare = {   # Round 72c: the reverse pairs of a bidirectional study, prepared from OpenGWAS where needed
          rows <- lapply(a$pairs, function(q) tryCatch(c(bidir_reverse(q[[1]], q[[2]], ogw_token_get(), isTRUE(settings_read()$offline)), list(ok = TRUE)),
                                                        error = function(e) list(ok = FALSE, msg = conditionMessage(e))))
          ev("bidir_ready", rows = rows) },
        run_begin = { if (isTRUE(settings_read()$awake)) keep_awake_start() },   # Round 72c
        run_end = { keep_awake_stop()
          if (isTRUE(settings_read()$notify)) { d <- as.integer(a$done %||% 0); n <- as.integer(a$total %||% 0)
            notify_done("MR Studio", sprintf("%s finished: %d of %d pair%s saved.", if (nzchar(a$study %||% "")) a$study else "The analysis", d, n, if (n == 1) "" else "s")) } },
        open_guide = {   # Round 72: the user guide and tutorial, in the default browser (a fixed address, never one from the page)
          r <- tryCatch({ .with_user_env(utils::browseURL(GUIDE_URL)); TRUE }, error = function(e) e)
          if (inherits(r, "error")) ev("error", msg = paste("The browser could not be opened. The guide is at", GUIDE_URL)) },
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
        fullgwas_plan = {
          r <- tryCatch(fullgwas_plan(a$id, ogw_token_get(), isTRUE(settings_read()$offline)), error = function(e) e)
          if (inherits(r, "error")) ev("fullgwas_plan", ok = FALSE, id = a$id, msg = conditionMessage(r)) else ev("fullgwas_plan", ok = TRUE, id = a$id, plan = r)
        },
        fullgwas_start = {
          r <- tryCatch(fullgwas_download_start(a$id, a$name, ogw_token_get(), !identical(a$keep, FALSE), isTRUE(settings_read()$offline)), error = function(e) e)
          ev("fullgwas_jobs", jobs = fullgwas_jobs(), msg = if (inherits(r, "error")) conditionMessage(r) else "", started = !inherits(r, "error"), id = a$id)
        },
        fullgwas_status = ev("fullgwas_jobs", jobs = fullgwas_jobs(), items = library_list(), catalogue = catalogue_all()),
        fullgwas_cancel = { fullgwas_cancel(a$key); ev("fullgwas_jobs", jobs = fullgwas_jobs()) },
        fullgwas_discard = { fullgwas_discard(a$key); ev("fullgwas_jobs", jobs = fullgwas_jobs()) },
        fullgwas_resume = {
          r <- tryCatch(fullgwas_resume(a$key, ogw_token_get(), isTRUE(settings_read()$offline)), error = function(e) e)
          ev("fullgwas_jobs", jobs = fullgwas_jobs(), msg = if (inherits(r, "error")) conditionMessage(r) else "")
        },
        vcf_import = {
          if (is.null(rv$vcf)) { ev("fullgwas_jobs", jobs = fullgwas_jobs(), msg = "Choose the VCF file again."); return() }
          r <- tryCatch(vcf_import_start(rv$vcf$path, a$name, list(ancestry = a$ancestry, kind = a$kind, ncase = a$ncase, ncontrol = a$ncontrol, n = a$n, file_name = rv$vcf$name)), error = function(e) e)
          ev("fullgwas_jobs", jobs = fullgwas_jobs(), msg = if (inherits(r, "error")) conditionMessage(r) else "", vcf_started = !inherits(r, "error"))
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
