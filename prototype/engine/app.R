# MR Studio engine: the first real screen (study, pair, IVW and Egger, results).
# Plain on purpose: the visual theme comes later. Defines `mr_ui` and `mr_server`; main.R starts them behind the gate.
library(shiny)

mr_css <- "
body { font-family: -apple-system, 'Segoe UI', system-ui, sans-serif; background:#f6f7f9; color:#1d2530; margin:0; }
.wrap { max-width:1280px; margin:0 auto; padding:18px 24px 40px; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:15px; margin:0 0 10px; text-transform:uppercase; letter-spacing:.04em; color:#4a5a6c; }
.sub { color:#5b6b7d; margin-bottom:18px; font-size:13px; }
.card { background:#fff; border:1px solid #dfe4ea; border-radius:10px; padding:14px 16px; margin-bottom:14px; }
.msg { margin:8px 0 0; font-size:13px; white-space:pre-wrap; } .msg.err { color:#a32121; } .msg.ok { color:#1c6b3c; }
.path { font-family: ui-monospace, Menlo, Consolas, monospace; font-size:12px; color:#4a5a6c; word-break:break-all; }
td { word-break:break-word; } .note { font-size:12px; color:#6a7a8c; } table { font-size:13px; } .btn-primary { background:#1f6feb; border-color:#1f6feb; }
.foot { margin-top:24px; font-size:11px; color:#8a97a6; }
"

mr_ui <- fluidPage(
  tags$head(tags$title("MR Studio"), tags$style(HTML(mr_css))),
  div(class = "wrap",
    h1("MR Studio"),
    div(class = "sub", "Two-sample Mendelian randomization. Every result comes with its script, its data table and a replay check."),
    fluidRow(
      column(5,
        div(class = "card",
          h2("1. Study"),
          selectInput("study_pick", "Open a study", choices = NULL),
          actionButton("study_open", "Open"),
          tags$hr(),
          textInput("study_name", "Or create a new study", placeholder = "For example: Omega-3 and triglycerides"),
          actionButton("study_create", "Create study", class = "btn-primary"),
          uiOutput("study_status")
        ),
        div(class = "card",
          h2("2. Pair and methods"),
          div(class = "note", "Demo data: 40 simulated instruments, not real GWAS data. Downloading public GWAS data comes in a later step."),
          textInput("exp_label", "Exposure", value = "Demo exposure"),
          textInput("out_label", "Outcome", value = "Demo outcome"),
          checkboxGroupInput("methods", "Methods", choices = METHOD_CHOICES, selected = unname(METHOD_CHOICES)),
          numericInput("ci_level", "Confidence level", value = 0.95, min = 0.5, max = 0.999, step = 0.01),
          actionButton("run", "Run analysis", class = "btn-primary"),
          uiOutput("run_status")
        )
      ),
      column(7,
        div(class = "card",
          h2("3. Results"),
          selectInput("pair_pick", "Pair", choices = NULL),
          tableOutput("results_tbl"),
          uiOutput("replay_line")
        ),
        div(class = "card",
          h2("4. Files of this study"),
          div(class = "note", "Everything below is a normal file in the study folder. The scripts rerun outside the app."),
          actionButton("show_folder", "Show study folder"),
          tableOutput("files_tbl")
        )
      )
    ),
    div(class = "foot", "MR GUI engine prototype. Engine: ", textOutput("hello", inline = TRUE))
  )
)

mr_server <- function(engine_dir) {
  function(input, output, session) {
    output$hello <- renderText("connected")
    rv <- reactiveValues(h = NULL, msg = NULL, msg_class = "ok", run_msg = NULL, run_class = "ok", tick = 0L)
    demo <- demo_pair(engine_dir)

    refresh_picker <- function(selected = NULL) updateSelectInput(session, "study_pick", choices = studies_list(), selected = selected)
    observe(refresh_picker())

    set_study <- function(h) {
      if (!is.null(rv$h)) try(DBI::dbDisconnect(rv$h$con), silent = TRUE)
      rv$h <- h; rv$tick <- rv$tick + 1L
      pf <- studio_pair_folders(h$dir)
      updateSelectInput(session, "pair_pick", choices = pf, selected = if (length(pf)) pf[length(pf)] else NULL)
    }
    say <- function(text, cls = "ok") { rv$msg <- text; rv$msg_class <- cls }
    onSessionEnded(function() isolate(if (!is.null(rv$h)) try(DBI::dbDisconnect(rv$h$con), silent = TRUE)))

    observeEvent(input$study_create, {
      tryCatch({
        h <- studio_new(input$study_name); set_study(h)
        refresh_picker(basename(h$dir)); say("Study created.")
      }, error = function(e) say(conditionMessage(e), "err"))
    })
    observeEvent(input$study_open, {
      req(nzchar(input$study_pick))
      tryCatch({ set_study(studio_open(input$study_pick)); say("Study opened.") },
               error = function(e) say(conditionMessage(e), "err"))
    })
    output$study_status <- renderUI({
      if (is.null(rv$h) && is.null(rv$msg)) return(div(class = "note", "Studies are kept in: ", span(class = "path", studies_root())))
      tagList(if (!is.null(rv$h)) div(class = "msg", "Current study: ", span(class = "path", rv$h$dir)),
              if (!is.null(rv$msg)) div(class = paste("msg", rv$msg_class), rv$msg))
    })

    observeEvent(input$run, {
      if (is.null(rv$h)) { rv$run_msg <- "Create or open a study first."; rv$run_class <- "err"; return() }
      tryCatch({
        withProgress(message = "Running the analysis", detail = "Generating scripts, running them and checking the replay", value = 0.3, {
          r <- studio_run(rv$h, input$exp_label, input$out_label, demo$input, input$methods, input$ci_level)
          setProgress(1)
        })
        pf <- studio_pair_folders(rv$h$dir); rv$tick <- rv$tick + 1L
        updateSelectInput(session, "pair_pick", choices = pf, selected = r$pair_folder)
        failed <- if (isTRUE(r$reused)) character() else names(Filter(function(x) !identical(x$status, "ok"), r$results))
        if (isTRUE(r$reused)) { rv$run_msg <- "Nothing changed since the last run, so the saved results were reused."; rv$run_class <- "ok" }
        else if (length(failed)) { rv$run_msg <- paste0("Some methods did not finish: ", paste(failed, collapse = ", "), ". ",
                                                          paste(unlist(lapply(r$results[failed], `[[`, "messages")), collapse = " ")); rv$run_class <- "err" }
        else { rv$run_msg <- "Analysis finished. The scripts and result files are in the study folder."; rv$run_class <- "ok" }
      }, error = function(e) { rv$run_msg <- conditionMessage(e); rv$run_class <- "err" })
    })
    output$run_status <- renderUI(if (!is.null(rv$run_msg)) div(class = paste("msg", rv$run_class), rv$run_msg))

    output$results_tbl <- renderTable({
      rv$tick; req(rv$h, nzchar(input$pair_pick))
      res <- pair_results(rv$h$dir, input$pair_pick)
      validate(need(!is.null(res), "No results for this pair yet."))
      format_results(res)
    }, na = "")
    output$replay_line <- renderUI({
      rv$tick; req(rv$h, nzchar(input$pair_pick))
      rp <- pair_replay(rv$h$dir, input$pair_pick)
      if (is.null(rp)) return(NULL)
      ok <- identical(rp$status, "identical")
      div(class = paste("msg", if (ok) "ok" else "err"),
          if (ok) sprintf("Replay check: the saved script, run again in a fresh R session, gave identical results (%s values compared).", rp$n_values_compared)
          else paste("Replay check:", rp$status))
    })
    output$files_tbl <- renderTable({
      rv$tick; req(rv$h)
      f <- study_files(rv$h$dir)
      validate(need(nrow(f) > 0, "No files yet."))
      data.frame(File = f$relative_path, Kind = f$role, Size = fmt_bytes(f$bytes), stringsAsFactors = FALSE)
    }, na = "")
    observeEvent(input$show_folder, {
      req(rv$h)
      tryCatch(open_in_file_manager(rv$h$dir), error = function(e) say(conditionMessage(e), "err"))
    })
  }
}
