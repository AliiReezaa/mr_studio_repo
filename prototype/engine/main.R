# MR GUI engine entry point (contracts C10 handshake).
# The shell starts:  Rscript main.R   and writes ONE JSON line to stdin:
#   {"secret": "<64 hex>", "login_token": "<64 hex>", "parent_pid": <n, optional>}
# The engine replies on stdout with ONE JSON line once it is listening:
#   {"ready": true, "port": <n>, "pid": <n>}

args_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(args_dir, "gate.R"))

con <- file("stdin", open = "r")
init <- jsonlite::fromJSON(readLines(con, n = 1L, warn = FALSE))
close(con)
stopifnot(grepl("^[0-9a-f]{64}$", init$secret), grepl("^[0-9a-f]{64}$", init$login_token))

port <- httpuv::randomPort(host = "127.0.0.1")
install_gate(init$secret, init$login_token, port)
parent_pid <- init$parent_pid
rm(init)

# Exit when the desktop shell is gone (crash or force-quit), so no orphan engine
# keeps running. macOS/Linux: poll the parent with signal 0. Windows: the shell
# places the engine in a job object that the OS kills together with the shell.
watch_parent <- function() {
  if (!is.null(parent_pid) && .Platform$OS.type == "unix" && !isTRUE(tools::pskill(parent_pid, 0L))) {
    message("desktop shell is gone; engine exiting")
    shiny::stopApp()
    return(invisible())
  }
  later::later(watch_parent, 2)
}

# The screen (see studio_core.R and app.R). Contracts v1.1 code is loaded once at start-up.
source(file.path(args_dir, "studio_core.R"))
load_contracts(find_contracts_dir(args_dir))
source(file.path(args_dir, "import_core.R"))
options(mrgui.engine_dir = args_dir)
source(file.path(args_dir, "opengwas_core.R"))
source(file.path(args_dir, "ld_core.R"))
source(file.path(args_dir, "report_core.R"))
source(file.path(args_dir, "app.R"))
ui <- mr_ui_for(args_dir)
server <- mr_server(args_dir)

app <- shiny::shinyApp(ui, server, onStart = function() {
  later::later(watch_parent, 2)
  later::later(function() {
    cat(jsonlite::toJSON(list(ready = TRUE, port = port, pid = Sys.getpid()), auto_unbox = TRUE), "\n", sep = "")
    flush(stdout())
  }, 0)
})
shiny::runApp(app, host = "127.0.0.1", port = port, launch.browser = FALSE, quiet = TRUE)
