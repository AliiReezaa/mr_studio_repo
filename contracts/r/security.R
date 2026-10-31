# MR GUI contracts v1 — security helpers (C12)

#' Formula-injection guard for text cells in CSV/XLSX exports.
csv_safe_text <- function(x) {
  x <- as.character(x)
  bad <- !is.na(x) & grepl("^[=+@\t\r-]", x)
  x[bad] <- paste0("'", x[bad])
  x
}

#' Escape text from metadata before it reaches HTML (UI or reports).
html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE); gsub("'", "&#39;", x, fixed = TRUE)
}

#' Join a user- or archive-supplied relative path under `root`, refusing escapes
#' (zip-slip, `..`, absolute paths, drive letters).
safe_path_join <- function(root, rel) {
  if (grepl("^([A-Za-z]:|/|\\\\)", rel) || grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", rel))
    stop("unsafe path refused: ", rel)
  root_n <- normalizePath(root, winslash = "/", mustWork = TRUE)
  file.path(root_n, rel)
}

#' Remove secrets from log or diagnostic text.
redact <- function(text, secrets) {
  for (s in secrets[nzchar(secrets)]) text <- gsub(s, "[REDACTED]", text, fixed = TRUE)
  # JWT-shaped strings (OpenGWAS tokens are JWTs)
  gsub("eyJ[A-Za-z0-9_-]{10,}\\.[A-Za-z0-9_-]{10,}\\.[A-Za-z0-9_-]{10,}", "[REDACTED]", text)
}
