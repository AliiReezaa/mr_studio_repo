# Reference adapters v1.1.0 (DRAFT): ivw_mre and egger gain a parameter schema with hints,
# column requirements and emit_script() (contracts v1.1.0, C6.1 and C6.2).
# The in-process run() of v1.0.0 is kept unchanged and is used only as a test oracle:
# the engine runs the generated script (no second code path).
# Requires reference_adapters.R, adapter.R and script.R to be sourced first.

.ci_param_v11 <- list(ci_level = list(
  type = "double", default = 0.95, min = 0.5, max = 0.999,
  label = "Confidence level",
  hint = "Width of the confidence interval for the causal effect. 95% is the usual reporting convention. It changes only the interval, not the estimate or the P value.",
  package_arg = NA_character_, kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "Reporting convention (95% interval); changes only the interval limits.",
  affects = c("ci")))

.minf_param_v11 <- list(min_f = list(
  type = "double", default = 0, min = 0, max = 100,
  label = "Minimum instrument strength (F)",
  hint = "Variants whose single-variant F statistic, (beta_x / se_x)^2, is below this value are left out before the analysis. 10 is the usual rule of thumb for weak instruments. 0 keeps every variant.",
  package_arg = NA_character_, kind = "float", conventional = FALSE, advanced = FALSE,
  evidence = NA_character_, affects = c("estimate", "ci", "p")))

.ivw_model_param_v11 <- list(ivw_model = list(
  type = "enum", default = "mre", values = c("mre", "mre_floor", "fixed"),
  value_labels = c(mre = "Multiplicative random effects", mre_floor = "Random effects when the variants are over-dispersed", fixed = "Fixed effect"),
  value_hints = c(mre = "The standard error grows with heterogeneity and may also shrink when the variants agree better than chance. TwoSampleMR function mr_ivw_mre. The usual default.",
                  mre_floor = "Random-effects standard error when heterogeneity is present, the fixed-effect standard error otherwise. TwoSampleMR function mr_ivw.",
                  fixed = "Assumes every variant estimates the same effect, so heterogeneity does not widen the interval. TwoSampleMR function mr_ivw_fe."),
  label = "Effects model",
  hint = "How between-variant heterogeneity enters the standard error. The estimate is the same in all three; only its standard error, interval and P value change.",
  package_arg = NA_character_, kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Multiplicative random effects is the default of TwoSampleMR::mr_ivw_mre and the usual choice when variants may differ in effect.",
  affects = c("se", "ci", "p")))

# Single-variant F filter, written into the script as plain base R (no estimator).
.script_lines_minf <- function(min_ivs) c(
  "# ---- 4b. Instrument strength filter (setting: min_f) ---------------------------",
  "f_stat <- (dat$beta_x / dat$se_x)^2   # single-variant F statistic",
  "n_before <- nrow(dat)",
  "dat <- dat[f_stat >= min_f, , drop = FALSE]",
  "if (n_before > nrow(dat)) cat(sprintf(\"Left out %d variant(s) with F below %g.\\n\", n_before - nrow(dat), min_f))",
  sprintf("if (nrow(dat) < %dL) stop(\"Too few variants remain after the F filter (%d needed). Lower the minimum F.\")", as.integer(min_ivs), as.integer(min_ivs)),
  "")

.colreq_basic <- list(
  list(column = "beta_x", level = "required", why = "Effect of each variant on the exposure.", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
  list(column = "se_x", level = "required", why = "Standard error of the exposure effect (needed for instrument strength).", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
  list(column = "beta_y", level = "required", why = "Effect of each variant on the outcome.", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
  list(column = "se_y", level = "required", why = "Standard error of the outcome effect (gives each variant its weight).", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
  list(column = "eaf_x", level = "optional", why = "Not used by this method; helps allele checks upstream.", if_missing = list(action = "ignore", code = NA_character_)))

.upgrade_adapter <- function(a, emit, extra = list()) {
  mrgui_adapter(id = a$id, adapter_version = "1.1.0", label = a$label, input_kind = a$input_kind,
                estimands = a$estimands, requires = a$requires, packages = a$packages,
                stochastic = a$stochastic, supports_checkpoint = a$supports_checkpoint,
                runtime_class = a$runtime_class, params = c(.ci_param_v11, extra, .minf_param_v11), check = a$check, run = a$run,
                emit_script = emit, column_requirements = .colreq_basic)
}

.est_cols <- function(spec) c(
  sprintf("exposure_id = %s, outcome_id = %s, scale = %s,", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale)))

.emit_ivw_mre <- function(spec, params, seed) {
  a <- ADAPTER_IVW_MRE_V11
  txt <- c(
    .script_lines_head(a, spec), .script_lines_guard(spec),
    "# ---- 2. Settings (every adjustable argument of this method) --------------------",
    .param_comment_lines(a, params),
    sprintf("seed <- %d   # fixed seed (this method is deterministic; kept for the record)", as.integer(seed)),
    "",
    .script_lines_repro(seed), .script_lines_data, .script_lines_minf(a$requires$min_ivs),
    "# ---- 5. Analysis: inverse-variance weighted estimate (TwoSampleMR) ----------------",
    "fit_fun <- switch(ivw_model,",
    "  mre = TwoSampleMR::mr_ivw_mre,         # multiplicative random effects",
    "  mre_floor = TwoSampleMR::mr_ivw,       # random effects when over-dispersed, fixed effect otherwise",
    "  fixed = TwoSampleMR::mr_ivw_fe)        # fixed effect",
    "fit <- fit_fun(",
    "  b_exp = dat$beta_x, b_out = dat$beta_y,",
    "  se_exp = dat$se_x, se_out = dat$se_y,",
    "  parameters = TwoSampleMR::default_parameters())",
    "z <- stats::qnorm(1 - (1 - ci_level) / 2)   # normal quantile for the interval",
    "estimates <- data.frame(",
    "  estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit$b, se = fit$se,",
    "  ci_lower = fit$b - z * fit$se, ci_upper = fit$b + z * fit$se, ci_level = ci_level,",
    "  ci_method = \"wald\", p = fit$pval, n_ivs_used = as.integer(fit$nsnp),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- data.frame(name = \"cochran_q\", value = fit$Q, df = fit$Q_df, p = fit$Q_pval,",
    "                          note = NA_character_, stringsAsFactors = FALSE)",
    "",
    .script_lines_write(spec))
  script_bundle(paste0(paste(txt, collapse = "\n"), "\n"), a$packages,
                list(list(file = file.path("04_tables", spec$estimates_file), role = "estimates"),
                     list(file = file.path("04_tables", spec$diagnostics_file), role = "diagnostics")),
                seed, spec$data_file)
}

.emit_egger <- function(spec, params, seed) {
  a <- ADAPTER_EGGER_V11
  txt <- c(
    .script_lines_head(a, spec), .script_lines_guard(spec),
    "# ---- 2. Settings (every adjustable argument of this method) --------------------",
    .param_comment_lines(a, params),
    sprintf("seed <- %d   # fixed seed (this method is deterministic; kept for the record)", as.integer(seed)),
    "",
    .script_lines_repro(seed), .script_lines_data, .script_lines_minf(a$requires$min_ivs),
    "# ---- 5. Analysis: TwoSampleMR::mr_egger_regression (orients variants internally) -",
    "fit <- TwoSampleMR::mr_egger_regression(",
    "  b_exp = dat$beta_x, b_out = dat$beta_y,",
    "  se_exp = dat$se_x, se_out = dat$se_y,",
    "  parameters = TwoSampleMR::default_parameters())",
    "tq <- stats::qt(1 - (1 - ci_level) / 2, fit$Q_df)   # t quantile, df = number of variants - 2",
    "estimates <- rbind(",
    "  data.frame(estimand = \"causal_effect\",",
    paste0("    ", .est_cols(spec)),
    "    estimate = fit$b, se = fit$se, ci_lower = fit$b - tq * fit$se, ci_upper = fit$b + tq * fit$se,",
    "    ci_level = ci_level, ci_method = \"wald_t\", p = fit$pval, n_ivs_used = as.integer(fit$nsnp),",
    "    stringsAsFactors = FALSE),",
    "  data.frame(estimand = \"pleiotropy_intercept\",",
    paste0("    ", .est_cols(spec)),
    "    estimate = fit$b_i, se = fit$se_i, ci_lower = fit$b_i - tq * fit$se_i, ci_upper = fit$b_i + tq * fit$se_i,",
    "    ci_level = ci_level, ci_method = \"wald_t\", p = fit$pval_i, n_ivs_used = as.integer(fit$nsnp),",
    "    stringsAsFactors = FALSE))",
    "diagnostics <- data.frame(name = \"rucker_q\", value = fit$Q, df = fit$Q_df, p = fit$Q_pval,",
    "                          note = NA_character_, stringsAsFactors = FALSE)",
    "",
    .script_lines_write(spec))
  script_bundle(paste0(paste(txt, collapse = "\n"), "\n"), a$packages,
                list(list(file = file.path("04_tables", spec$estimates_file), role = "estimates"),
                     list(file = file.path("04_tables", spec$diagnostics_file), role = "diagnostics")),
                seed, spec$data_file)
}

# ---- Round 8: robust estimators (TwoSampleMR). All three bootstrap their standard error, so they are stochastic:
# the script sets a fixed seed and the replay check repeats the same draws. ----
.nboot_param_v11 <- list(nboot = list(
  type = "integer", default = 1000L, min = 100L, max = 10000L,
  label = "Bootstrap replications",
  hint = "How many times the data are resampled to get the standard error. More replications give a steadier standard error but take longer. The seed is fixed, so the same number always gives the same answer.",
  package_arg = "parameters$nboot", kind = "int", conventional = TRUE, advanced = FALSE,
  evidence = "1000 replications is the TwoSampleMR default and the usual choice in published analyses.",
  affects = c("se", "ci", "p")))

.phi_param_v11 <- list(phi = list(
  type = "double", default = 1, min = 0.05, max = 5,
  label = "Kernel bandwidth multiplier (phi)",
  hint = "Scales the smoothing of the density whose peak is the mode estimate. Smaller values follow the data more closely but are noisier. 1 is the default; 0.5 and 0.25 are common sensitivity checks.",
  package_arg = "parameters$phi", kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "phi = 1 is the default of Hartwig, Davey Smith and Bowden (2017); 0.5 and 0.25 are the sensitivity values they suggest.",
  affects = c("estimate", "se", "ci", "p")))

.emit_robust <- function(a, fun, title, spec, params, seed, pars) {
  txt <- c(
    .script_lines_head(a, spec), .script_lines_guard(spec),
    "# ---- 2. Settings (every adjustable argument of this method) --------------------",
    .param_comment_lines(a, params),
    sprintf("seed <- %d   # fixed seed: the standard error is bootstrapped, so the same seed gives the same numbers", as.integer(seed)),
    "",
    .script_lines_repro(seed), .script_lines_data, .script_lines_minf(a$requires$min_ivs),
    sprintf("# ---- 5. Analysis: %s (TwoSampleMR) ----------------", title),
    sprintf("fit <- TwoSampleMR::%s(", fun),
    "  b_exp = dat$beta_x, b_out = dat$beta_y,",
    "  se_exp = dat$se_x, se_out = dat$se_y,",
    sprintf("  parameters = utils::modifyList(TwoSampleMR::default_parameters(), list(%s)))", pars),
    "z <- stats::qnorm(1 - (1 - ci_level) / 2)   # normal interval around the estimate, using the bootstrapped standard error",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit$b, se = fit$se, ci_lower = fit$b - z * fit$se, ci_upper = fit$b + z * fit$se,",
    "  ci_level = ci_level, ci_method = \"wald\", p = fit$pval, n_ivs_used = as.integer(fit$nsnp),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- data.frame(name = \"bootstrap_replications\", value = as.double(nboot), df = NA_real_, p = NA_real_,",
    "                          note = NA_character_, stringsAsFactors = FALSE)",
    "",
    .script_lines_write(spec))
  script_bundle(paste0(paste(txt, collapse = "\n"), "\n"), a$packages,
                list(list(file = file.path("04_tables", spec$estimates_file), role = "estimates"),
                     list(file = file.path("04_tables", spec$diagnostics_file), role = "diagnostics")),
                seed, spec$data_file)
}

.robust_adapter <- function(id, label, fun, title, extra, pars_fn) {
  ad <- NULL
  ad <- mrgui_adapter(
    id = id, adapter_version = "1.1.0", label = label, input_kind = "pair_iv", estimands = "causal_effect",
    requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                    needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
    packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"),
    stochastic = TRUE, supports_checkpoint = FALSE, runtime_class = "fast",
    params = c(.ci_param_v11, extra, .minf_param_v11),
    check = function(input, params) {
      d <- .usable(input$data)
      if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("%s needs at least 3 usable instruments; %d available.", label, nrow(d))))
      w <- .common_warnings(d)
      verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
    },
    run = function(input, params, ctx) {
      d <- .usable(input$data)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(ctx$seed))
      fit <- getExportedValue("TwoSampleMR", fun)(d$beta_x, d$beta_y, d$se_x, d$se_y, do.call(function(nboot = 1000L, phi = 1) utils::modifyList(TwoSampleMR::default_parameters(), list(nboot = nboot, phi = phi)), params[intersect(names(params), c("nboot", "phi"))]))
      z <- stats::qnorm(1 - (1 - params$ci_level) / 2)
      list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
             fit$b, fit$se, fit$b - z * fit$se, fit$b + z * fit$se, "wald", fit$pval, as.integer(fit$nsnp), params$ci_level),
           diagnostics = diagnostic_row("bootstrap_replications", as.double(params$nboot)))
    },
    emit_script = function(spec, params, seed) .emit_robust(ad, fun, title, spec, params, seed, pars_fn),
    column_requirements = .colreq_basic)
  ad
}
ADAPTER_WM_V11 <- .robust_adapter("wm", "Weighted median", "mr_weighted_median", "weighted median estimator", .nboot_param_v11, "nboot = nboot")
ADAPTER_SMED_V11 <- .robust_adapter("smed", "Simple median", "mr_simple_median", "simple median estimator", .nboot_param_v11, "nboot = nboot")
ADAPTER_MODE_V11 <- .robust_adapter("mode", "Weighted mode", "mr_weighted_mode", "weighted mode estimator (MBE)", c(.nboot_param_v11, .phi_param_v11), "nboot = nboot, phi = phi")

ADAPTER_IVW_MRE_V11 <- .upgrade_adapter(ADAPTER_IVW_MRE, function(spec, params, seed) .emit_ivw_mre(spec, params, seed), .ivw_model_param_v11)
ADAPTER_EGGER_V11 <- .upgrade_adapter(ADAPTER_EGGER, function(spec, params, seed) .emit_egger(spec, params, seed))

#' C13 conformance 1: parameter schema is complete and every default is valid.
validate_param_schema <- function(adapter) {
  problems <- character()
  for (p in names(adapter$params)) {
    s <- adapter$params[[p]]
    for (f in c("label", "hint", "kind", "conventional", "advanced"))
      if (is.null(s[[f]])) problems <- c(problems, sprintf("%s: missing %s", p, f))
    if (!is.null(s$hint) && nchar(s$hint) > 280) problems <- c(problems, sprintf("%s: hint longer than 280 characters", p))
    if (isTRUE(s$conventional) && (is.null(s$evidence) || !nzchar(s$evidence))) problems <- c(problems, sprintf("%s: conventional without evidence", p))
    ok <- !inherits(try(resolve_params(adapter, setNames(list(s$default), p)), silent = TRUE), "try-error")
    if (!ok) problems <- c(problems, sprintf("%s: default is not valid", p))
  }
  problems
}
