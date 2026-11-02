# Real https request with the moved runtime (run by the build on every operating system).
# Fails the build if the bundled R cannot make a verified https call to OpenGWAS, which is what broke on the first real test.
args_dir <- normalizePath(file.path("prototype", "engine"), winslash = "/")
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a
redact <- function(x, ...) x
options(mrgui.engine_dir = args_dir)
source(file.path(args_dir, "opengwas_core.R"))
ogw_ca_setup(force = TRUE)
cat("usable certificate files, in the order they are tried:\n"); print(.ogw_ca_state$list)
r <- ogw_ping()
cat("ping:", r$ok, "-", r$detail, "\n")
cat("certificate file in use:", .ogw_ca_state$path, "\n")
if (!isTRUE(r$ok)) quit(status = 1)
