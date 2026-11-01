# MR GUI contracts v1 — adapter interface, verdicts, results (C6, C7, C11)

INPUT_KINDS <- c("pair_iv", "pair_iv_ld", "multi_iv", "region", "genome_wide_pair")
ESTIMANDS <- c("causal_effect", "direct_effect", "indirect_effect", "total_effect",
               "pleiotropy_intercept", "proportion_mediated", "colocalisation_pp",
               "genetic_correlation", "reverse_effect")
SCALES <- c("beta", "log_or", "log_hr", "probability", "correlation", "proportion")
CI_METHODS <- c("wald", "wald_t", "profile", "bootstrap", "bayesian", "robust", "as_reported")
RUNTIME_CLASSES <- c("fast", "moderate", "slow", "very_slow")

.adapter_registry <- new.env(parent = emptyenv())

#' Build and validate an adapter (C6.1). `params` is a named list of
#' list(type = "double"|"integer"|"logical"|"enum", default, min, max, values).
mrgui_adapter <- function(id, adapter_version, label, input_kind, estimands, requires,
                          packages, stochastic, supports_checkpoint, runtime_class,
                          params = list(), check, run,
                          emit_script = NULL, column_requirements = NULL) {
  stopifnot(is.character(id), length(id) == 1L, nzchar(id))
  stopifnot(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", adapter_version))
  stopifnot(input_kind %in% INPUT_KINDS, all(estimands %in% ESTIMANDS))
  stopifnot(runtime_class %in% RUNTIME_CLASSES, is.logical(stochastic), is.logical(supports_checkpoint))
  stopifnot(is.list(requires), is.numeric(requires$min_ivs))
  stopifnot(is.list(packages), !is.null(names(packages)))
  stopifnot(is.function(check), is.function(run))
  stopifnot(is.null(emit_script) || is.function(emit_script))   # v1.1 (additive)
  for (p in names(params)) {
    spec <- params[[p]]
    stopifnot(spec$type %in% c("double", "integer", "logical", "enum"), !is.null(spec$default))
  }
  structure(list(id = id, adapter_version = adapter_version, label = label,
                 input_kind = input_kind, estimands = estimands, requires = requires,
                 packages = packages, stochastic = stochastic,
                 supports_checkpoint = supports_checkpoint, runtime_class = runtime_class,
                 params = params, check = check, run = run,
                 emit_script = emit_script, column_requirements = column_requirements),
            class = "mrgui_adapter")
}

register_adapter <- function(adapter, inventory_ids = NULL) {
  stopifnot(inherits(adapter, "mrgui_adapter"))
  if (!is.null(inventory_ids) && !(adapter$id %in% inventory_ids))
    stop("adapter id '", adapter$id, "' is not a method id in mr_inventory.json")
  assign(adapter$id, adapter, envir = .adapter_registry)
  invisible(adapter)
}
get_adapter <- function(id) get(id, envir = .adapter_registry, inherits = FALSE)

#' Merge user params with defaults and validate (only declared params allowed).
resolve_params <- function(adapter, params = list()) {
  unknown <- setdiff(names(params), names(adapter$params))
  if (length(unknown)) stop("unknown parameter(s) for ", adapter$id, ": ", paste(unknown, collapse = ", "))
  out <- lapply(adapter$params, `[[`, "default")
  for (p in names(params)) {
    spec <- adapter$params[[p]]; v <- params[[p]]
    v <- switch(spec$type,
      double = as.double(v), integer = as.integer(v), logical = as.logical(v),
      enum = { if (!(v %in% spec$values)) stop(p, " must be one of ", paste(spec$values, collapse = "/")); as.character(v) })
    if (length(v) != 1L || is.na(v)) stop(p, " must be a single non-missing value")
    if (!is.null(spec$min) && v < spec$min) stop(p, " below minimum ", spec$min)
    if (!is.null(spec$max) && v > spec$max) stop(p, " above maximum ", spec$max)
    out[[p]] <- v
  }
  out
}

verdict <- function(status = c("ok", "warn", "block"), codes = character(), messages = character()) {
  status <- match.arg(status)
  list(status = status, codes = as.character(codes), messages = as.character(messages))
}

#' One row of the estimates table (C7.2).
estimate_row <- function(estimand, exposure_id, outcome_id, scale, estimate, se,
                         ci_lower, ci_upper, ci_method, p, n_ivs_used, ci_level = 0.95) {
  stopifnot(estimand %in% ESTIMANDS, scale %in% SCALES, ci_method %in% CI_METHODS)
  data.frame(estimand = estimand, exposure_id = exposure_id, outcome_id = outcome_id,
             scale = scale, estimate = as.double(estimate), se = as.double(se),
             ci_lower = as.double(ci_lower), ci_upper = as.double(ci_upper),
             ci_level = as.double(ci_level), ci_method = ci_method,
             p = as.double(p), n_ivs_used = as.integer(n_ivs_used), stringsAsFactors = FALSE)
}

diagnostic_row <- function(name, value, df = NA_real_, p = NA_real_, note = NA_character_) {
  data.frame(name = name, value = as.double(value), df = as.double(df), p = as.double(p),
             note = as.character(note), stringsAsFactors = FALSE)
}

.df_rows <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(list())
  lapply(seq_len(nrow(df)), function(i) {
    r <- as.list(df[i, , drop = FALSE])
    lapply(r, function(v) { v <- unname(v); if (is.double(v) && !is.na(v) && !is.finite(v)) NA else v })
  })
}

#' Serialise a MethodResult (C7) to canonical JSON for storage.
result_to_json <- function(res) {
  x <- res
  x$estimates <- .df_rows(res$estimates)
  x$diagnostics <- .df_rows(res$diagnostics)
  x$codes <- I(as.character(res$codes))
  x$messages <- I(as.character(res$messages))
  x$provenance$packages <- as.list(res$provenance$packages)
  canonical_json(x)
}

#' Manifest as canonical JSON (functions excluded), for the UI and feasibility engine.
manifest_json <- function(adapter) {
  m <- unclass(adapter)[c("id", "adapter_version", "label", "input_kind", "estimands", "requires",
                          "packages", "stochastic", "supports_checkpoint", "runtime_class", "params")]
  m$estimands <- I(m$estimands)
  m$requires$min_ivs <- as.integer(m$requires$min_ivs)
  m$requires$min_exposures <- as.integer(m$requires$min_exposures)
  m$requires$max_exposures <- as.integer(m$requires$max_exposures)
  canonical_json(m)
}
