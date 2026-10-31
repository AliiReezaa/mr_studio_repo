# MR GUI engine — connection gate in front of Shiny (contracts C10)
#
# Shiny serves through httpuv. Every request that reaches R passes through
# shiny:::handlerManager (HTTP) or its WebSocket handler list. We prepend one
# gate to each list, so the gate runs before any Shiny code, for every HTTP
# request handled in R and for every WebSocket connection.
#
# Requests for Shiny's bundled static files (/shared/..., resource paths) are
# served by httpuv in C++ without entering R, so they are NOT gated. Rule: only
# public library files may be static; user data (results, figures, reports) is
# never registered as a resource path and is served only through gated routes.
#
# Uses shiny internals (handlerManager); pinned and covered by the gate test.

COOKIE_NAME <- "mrgui_s"

.const_eq <- function(a, b) {
  a <- charToRaw(enc2utf8(a %||% "")); b <- charToRaw(enc2utf8(b %||% ""))
  if (length(a) != length(b)) return(FALSE)
  sum(bitwXor(as.integer(a), as.integer(b))) == 0   # bitwise; base xor() is logical and would accept any token
}
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.cookie_value <- function(header, name) {
  if (is.null(header) || !nzchar(header)) return(NULL)
  for (part in strsplit(header, ";", fixed = TRUE)[[1]]) {
    kv <- strsplit(trimws(part), "=", fixed = TRUE)[[1]]
    if (length(kv) >= 2 && kv[1] == name) return(paste(kv[-1], collapse = "="))
  }
  NULL
}

.forbidden <- function() shiny::httpResponse(403L, "text/plain", "Forbidden")

install_gate <- function(secret, login_token, port) {
  force(secret); force(login_token); force(port)   # the caller deletes its copies
  expected_host <- sprintf("127.0.0.1:%d", port)
  expected_origin <- sprintf("http://127.0.0.1:%d", port)
  token_state <- new.env(); token_state$unused <- TRUE

  host_ok <- function(req) identical(req$HTTP_HOST, expected_host)
  cookie_ok <- function(req) .const_eq(.cookie_value(req$HTTP_COOKIE, COOKIE_NAME), secret)

  # SameSite=Lax, not Strict: the shell window starts on tauri://localhost and navigates to
  # http://127.0.0.1:<port>/auth. WebKit treats that chain as cross-site initiated and withholds
  # Strict cookies on the redirect to "/", which produced "Forbidden" in the real app (found
  # 2026-11-01). Lax is sent on top-level GET navigations; the Host check, the Origin check on
  # WebSockets and the secret in the cookie remain the real protection.
  http_gate <- function(req) {
    if (!host_ok(req)) return(.forbidden())                    # DNS rebinding
    if (identical(req$PATH_INFO, "/auth")) {                    # one-time login
      t <- shiny::parseQueryString(req$QUERY_STRING %||% "")$t
      if (!isTRUE(token_state$unused) || !.const_eq(t, login_token)) return(.forbidden())
      token_state$unused <- FALSE
      return(shiny::httpResponse(302L, "text/plain", "", headers = list(
        Location = "/",
        `Set-Cookie` = sprintf("%s=%s; Path=/; HttpOnly; SameSite=Lax", COOKIE_NAME, secret),
        `Cache-Control` = "no-store")))
    }
    if (!cookie_ok(req)) return(.forbidden())
    NULL                                                         # pass to Shiny
  }

  ws_gate <- function(ws) {
    req <- ws$request
    if (!host_ok(req) || !identical(req$HTTP_ORIGIN, expected_origin) || !cookie_ok(req)) {
      ws$close()
      return(TRUE)                                               # handled: stop here
    }
    NULL                                                         # continue to Shiny (HandlerList stops at the first non-NULL result)
  }

  hm <- shiny:::handlerManager
  hm$addHandler(http_gate, "mrgui_gate", tail = FALSE)
  hm$addWSHandler(ws_gate, "mrgui_gate", tail = FALSE)
  invisible(TRUE)
}
