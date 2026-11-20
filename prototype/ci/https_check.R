# Real https request with the moved runtime (run by the build on every operating system).
# What this guards: the bundled R must make verified https calls after the runtime was moved (its libcurl first remembered the
# build folder's certificate file, which broke the first real test). A certificate or SSL failure fails the build.
# What it must not do: fail the build because OpenGWAS is slow or briefly down (run #23 on macOS arm64: a 20 s timeout with
# 0 bytes received). So OpenGWAS is tried four times; if it still does not answer for a network reason, the same certificates are
# checked against a second host (CRAN). If that verified call works, https is fine and the build goes on with a warning.
args_dir <- normalizePath(file.path("prototype", "engine"), winslash = "/")
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a
redact <- function(x, ...) x
options(mrgui.engine_dir = args_dir)
source(file.path(args_dir, "opengwas_core.R"))
ogw_ca_setup(force = TRUE)
cat("usable certificate files, in the order they are tried:\n"); print(.ogw_ca_state$list)
cert_like <- function(d) grepl("certificate|SSL|CA file|CAfile|self.signed|unable to get local issuer|verify", d, ignore.case = TRUE)
r <- NULL
for (i in 1:4) {
  r <- ogw_ping(); cat(sprintf("ping %d: %s - %s\n", i, r$ok, r$detail))
  if (isTRUE(r$ok) || cert_like(r$detail)) break
  if (i < 4) Sys.sleep(10 * 2^(i - 1))   # 10, 20, 40 s
}
cat("certificate file in use:", .ogw_ca_state$path, "\n")
if (isTRUE(r$ok)) quit(status = 0)
if (cert_like(r$detail)) { cat("FAIL: a certificate or SSL error, which is what this check exists to catch.\n"); quit(status = 1) }
# OpenGWAS did not answer for a network reason: verify https with the same certificates against another host
alt <- tryCatch(httr::GET("https://cloud.r-project.org/", httr::timeout(30)), error = function(e) e)
if (inherits(alt, "error")) {
  cat("FAIL: OpenGWAS did not answer and the second host failed too:", conditionMessage(alt), "\n"); quit(status = 1)
}
cat(sprintf("Verified https to cloud.r-project.org works with the moved runtime (HTTP %d, certificates %s).\n", alt$status_code, .ogw_ca_state$path))
cat("::warning title=OpenGWAS unreachable::OpenGWAS did not answer from this runner (", substr(r$detail, 1, 150), "); https itself is verified, so the build goes on.\n", sep = "")
quit(status = 0)
