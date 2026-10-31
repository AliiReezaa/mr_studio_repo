# MR GUI contracts v1 — canonical JSON, hashing, seeds (C8.3–C8.5)
# Base R + digest only. This is the ONLY place in the engine that may produce
# JSON for hashing.

CONTRACT_VERSION <- "1.0.0"

# Works on raw UTF-8 bytes so the result never depends on the session locale
# (a C locale would otherwise turn non-ASCII characters into "<c3><a9>").
.json_escape <- function(s) {
  b <- as.integer(.utf8_bytes(s))
  out <- vector("list", length(b))
  for (i in seq_along(b)) {
    v <- b[i]
    out[[i]] <- if (v == 0x22L) charToRaw("\\\"")
                else if (v == 0x5CL) charToRaw("\\\\")
                else if (v < 0x20L) charToRaw(sprintf("\\u%04x", v))
                else as.raw(v)
  }
  .utf8(c(charToRaw("\""), unlist(out), charToRaw("\"")))
}
.utf8 <- function(r) { s <- rawToChar(r); Encoding(s) <- "UTF-8"; s }
# Bytes of a string as UTF-8 without locale-dependent re-encoding: strings that
# are already valid UTF-8 are taken as they are (in a C locale, enc2utf8() would
# otherwise mangle them); anything else is converted from its declared encoding.
.utf8_bytes <- function(s) {
  if (Encoding(s) %in% c("UTF-8", "unknown") && validUTF8(s)) charToRaw(s) else charToRaw(enc2utf8(s))
}
.cat_utf8 <- function(...) .utf8(unlist(lapply(list(...), function(p) .utf8_bytes(p))))
.join_utf8 <- function(open, parts, close) {
  if (!length(parts)) return(.cat_utf8(open, close))
  sep <- rep(list(","), length(parts)); sep[[length(parts)]] <- close
  do.call(.cat_utf8, c(list(open), as.vector(rbind(parts, sep))))
}

.json_scalar <- function(x) {
  if (length(x) != 1L) stop("internal: scalar expected")
  if (is.na(x)) return("null")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.integer(x)) return(as.character(x))
  if (is.double(x)) {
    if (!is.finite(x)) stop("canonical_json: Inf/NaN are not allowed")
    if (x == 0) return("0")  # also maps -0 to 0
    return(sprintf("%.17g", x))
  }
  if (is.character(x)) return(.json_escape(x))
  stop("canonical_json: unsupported scalar type ", typeof(x))
}

#' Canonical JSON (C8.5). Length-1 atomic vectors are scalars; wrap in I() to
#' force an array. Named lists are objects with keys sorted by code point;
#' unnamed lists are arrays. Factors, Inf and NaN are refused.
canonical_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.factor(x)) stop("canonical_json: factors are not allowed; use character")
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) {
      return(.join_utf8("[", lapply(x, canonical_json), "]"))
    }
    if (any(nm == "") || anyDuplicated(nm)) stop("canonical_json: object keys must be unique and non-empty")
    nm8 <- nm
    ord <- order(nm8, method = "radix")  # radix = byte order = code-point order for UTF-8
    parts <- lapply(ord, function(i) .cat_utf8(.json_escape(nm8[i]), ":", canonical_json(x[[i]])))
    return(.join_utf8("{", parts, "}"))
  }
  if (is.atomic(x)) {
    if (!is.null(names(x))) stop("canonical_json: named atomic vectors are ambiguous; use a named list")
    if (length(x) == 1L && !inherits(x, "AsIs")) return(.json_scalar(unclass(x)))
    v <- unclass(x)
    return(.join_utf8("[", lapply(seq_along(v), function(i) .json_scalar(v[i])), "]"))
  }
  stop("canonical_json: unsupported type ", typeof(x))
}

sha256_text <- function(s) digest::digest(.utf8_bytes(s), algo = "sha256", serialize = FALSE)
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

#' Fingerprint of a task (C8.3). `input_blobs` is a named list of blob hashes;
#' secrets must never be passed here.
task_fingerprint <- function(task_type, input_blobs, params, adapter = NULL, os_family = os_family_id()) {
  payload <- list(
    contract_version = CONTRACT_VERSION,
    task_type = task_type,
    inputs = input_blobs,
    params = params,
    adapter = if (is.null(adapter)) NULL else list(
      id = adapter$id, version = adapter$adapter_version, packages = adapter$packages),
    os = os_family
  )
  sha256_text(canonical_json(payload))
}

#' Seed = first 31 bits of SHA-256(task_id) (C8.4).
task_seed <- function(task_id) {
  h <- sha256_text(task_id)
  first32 <- sum(strtoi(strsplit(substr(h, 1, 8), "")[[1]], 16L) * 16^(7:0))
  as.integer(floor(first32 / 2))
}

os_family_id <- function() {
  s <- Sys.info()[["sysname"]]
  m <- R.version$arch
  if (s == "Darwin") return(if (grepl("aarch64|arm64", m)) "macos-arm64" else "macos-x64")
  if (s == "Windows") return("windows-x64")
  "linux-x64"
}
