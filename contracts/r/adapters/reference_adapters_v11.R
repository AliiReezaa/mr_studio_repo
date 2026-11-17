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
  type = "enum", default = "mre_floor", values = c("mre_floor", "mre", "fixed"),
  value_labels = c(mre_floor = "Multiplicative random effects (standard)", mre = "Multiplicative random effects without the fixed-effect floor", fixed = "Fixed effect"),
  value_hints = c(mre_floor = "The standard error grows with heterogeneity but never falls below the fixed-effect standard error. TwoSampleMR function mr_ivw, the IVW method TwoSampleMR::mr() runs by default; MendelianRandomization::mr_ivw gives the same result with 4 or more variants. Recommended by the Burgess et al. guidelines.",
                  mre = "The standard error grows with heterogeneity and also shrinks below the fixed-effect standard error when the variants agree better than chance, which makes the P value too small with few variants. TwoSampleMR function mr_ivw_mre.",
                  fixed = "Assumes every variant estimates the same effect, so heterogeneity does not widen the interval. TwoSampleMR function mr_ivw_fe."),
  label = "Effects model",
  hint = "How between-variant heterogeneity enters the standard error. The estimate is the same in all three; only its standard error, interval and P value change.",
  package_arg = NA_character_, kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Burgess et al. guidelines (PMID 32760811): multiplicative random effects as the default IVW model; the residual standard error is not allowed below 1, as in TwoSampleMR::mr_ivw and MendelianRandomization::mr_ivw.",
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
    "  mre_floor = TwoSampleMR::mr_ivw,       # multiplicative random effects, SE not below the fixed-effect SE (TwoSampleMR default)",
    "  mre = TwoSampleMR::mr_ivw_mre,         # multiplicative random effects without the fixed-effect floor",
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
ADAPTER_SMODE_V11 <- .robust_adapter("smode", "Simple mode", "mr_simple_mode", "simple mode estimator (unweighted MBE)", c(.nboot_param_v11, .phi_param_v11), "nboot = nboot, phi = phi")

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

# ---- Round 10: outlier detectors (MR-PRESSO, radial IVW, radial Egger). They report an estimate AND the list of variants
# they flag, so the app can offer the flagged variants as a labelled sensitivity analysis. They never remove anything themselves. ----
PIN_MRPRESSO <- "3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87"
PIN_RADIALMR <- "7ca35c052489001890347871e38998ca53a02b03"

.nbdist_param_v11 <- list(nbdist = list(
  type = "integer", default = 1000L, min = 500L, max = 20000L,
  label = "Simulated datasets (NbDistribution)",
  hint = "How many random datasets MR-PRESSO simulates for its null distribution. Must exceed the number of variants; the smallest possible P is 1 divided by it. More gives finer P values but takes longer.",
  package_arg = "NbDistribution", kind = "int", conventional = TRUE, advanced = FALSE,
  evidence = "1000 is the default of MR-PRESSO; Verbanck et al. (2018) used 10000 for finer P values.",
  affects = c("p", "outliers")))

.signif_param_v11 <- list(signif = list(
  type = "double", default = 0.05, min = 0.001, max = 0.5,
  label = "Significance threshold for outliers",
  hint = "A variant is flagged when its outlier test P value is below this threshold, and outlier detection only runs when the global test is below it. MR-PRESSO applies a Bonferroni correction over all variants inside its outlier test.",
  package_arg = "SignifThreshold", kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "0.05 is the default of MR-PRESSO (Verbanck et al. 2018).",
  affects = c("outliers", "estimate", "se", "ci", "p")))

.radial_alpha_param_v11 <- list(alpha = list(
  type = "double", default = 0.05, min = 0.0001, max = 0.5,
  label = "Outlier significance level (alpha)",
  hint = "A variant is flagged when its contribution to the Q statistic is significant at this level. The package applies it per variant with no multiplicity correction; dividing 0.05 by the number of variants gives a Bonferroni level.",
  package_arg = "alpha", kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "0.05 is the example value in the RadialMR documentation; Bowden et al. (2018) discuss a Bonferroni-corrected level.",
  affects = c("outliers")))

.radial_weights_param_v11 <- list(weights = list(
  type = "enum", default = "1", values = c("1", "2", "3"),
  value_labels = c("1" = "First order", "2" = "Second order", "3" = "Modified second order"),
  value_hints = c("1" = "Weights each variant by the inverse variance of its ratio estimate using first-order terms. Same estimate as IVW.",
                  "2" = "Adds the second-order term that accounts for uncertainty in the variant-exposure effect.",
                  "3" = "Modified second-order weights (Bowden et al. 2018), which behave better when instruments are weak."),
  label = "Weights",
  hint = "How the variants are weighted in the radial regression and in the outlier Q statistics.",
  package_arg = "weights", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Options of RadialMR; first order is the package default.",
  affects = c("estimate", "se", "ci", "p", "outliers")))

.radial_tol_param_v11 <- list(tol = list(
  type = "double", default = 0.0001, min = 1e-8, max = 0.01,
  label = "Convergence tolerance",
  hint = "Stopping rule of the iterative radial IVW estimate. Smaller values iterate longer. Rarely worth changing.",
  package_arg = "tol", kind = "float", conventional = FALSE, advanced = TRUE,
  evidence = NA_character_, affects = c("estimate")))

.script_lines_head_pk <- function(a, spec) {
  h <- .script_lines_head(a, spec)
  pins <- unlist(a$packages); kind <- ifelse(grepl("^[0-9]+(\\.[0-9]+)+$", pins), "version", "commit")   # Round 37: CRAN packages are pinned by version
  h[grepl("^# Package  :", h)] <- sprintf("# Package  : %s", paste(sprintf("%s %s %s", names(a$packages), kind, pins), collapse = "; "))
  h
}
.script_lines_guard_pk <- function(spec, pkg) c(
  "# ---- 1. Working directory check ----------------------------------------------",
  sprintf("DATA_FILE <- file.path(\"01_data\", %s)", .q(spec$data_file)),
  "if (!file.exists(DATA_FILE))",
  "  stop(\"Wrong working directory. Set it to the pair folder (the folder that contains 01_data). In RStudio, open the .Rproj file.\")",
  sprintf("if (!requireNamespace(\"%s\", quietly = TRUE))", pkg),
  sprintf("  stop(\"Package %s is not installed. See README_RERUN.md (renv::restore()).\")", pkg),
  "")

.detector_emit <- function(a, spec, params, seed, pkg, body, seed_note) {
  txt <- c(.script_lines_head_pk(a, spec), .script_lines_guard_pk(spec, pkg),
           "# ---- 2. Settings (every adjustable argument of this method) --------------------",
           .param_comment_lines(a, params),
           sprintf("seed <- %d   # %s", as.integer(seed), seed_note),
           "", .script_lines_repro(seed), .script_lines_data, .script_lines_minf(a$requires$min_ivs),
           body, "", .script_lines_write(spec))
  script_bundle(paste0(paste(txt, collapse = "\n"), "\n"), a$packages,
                list(list(file = file.path("04_tables", spec$estimates_file), role = "estimates"),
                     list(file = file.path("04_tables", spec$diagnostics_file), role = "diagnostics")),
                seed, spec$data_file)
}

.ids_note <- function(ids) if (length(ids)) paste(ids, collapse = ";") else NA_character_

# MR-PRESSO ----------------------------------------------------------------------------------------
ADAPTER_PRESSO_V11 <- NULL
ADAPTER_PRESSO_V11 <- mrgui_adapter(
  id = "presso", adapter_version = "1.1.0", label = "MR-PRESSO", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 4, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(MRPRESSO = PIN_MRPRESSO), stochastic = TRUE, supports_checkpoint = FALSE, runtime_class = "slow",
  params = c(.ci_param_v11, .nbdist_param_v11, .signif_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data); nb <- if (is.null(params$nbdist)) 1000L else params$nbdist
    if (nrow(d) < 4) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("MR-PRESSO needs at least 4 usable instruments; %d available.", nrow(d))))
    if (nrow(d) >= nb) return(verdict("block", "FEAS_PRESSO_NBDIST", sprintf("MR-PRESSO needs more simulated datasets (%d) than variants (%d). Raise the number of simulated datasets.", nb, nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(ctx$seed))
    pd <- data.frame(by = d$beta_y, bx = d$beta_x, sy = d$se_y, sx = d$se_x)
    fit <- suppressWarnings(MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", data = pd,
             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, SignifThreshold = params$signif, NbDistribution = params$nbdist, seed = as.integer(ctx$seed)))
    pr <- .presso_pick(fit, d$variant_id, nrow(d), params$ci_level, params$nbdist)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           pr$b, pr$se, pr$lo, pr$hi, "wald_t", pr$p, pr$n, params$ci_level),
         diagnostics = pr$diag)
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_PRESSO_V11, spec, params, seed, "MRPRESSO", c(
    "# ---- 5. Analysis: MR-PRESSO (global test, outlier test, distortion test) -------",
    "pd <- data.frame(by = dat$beta_y, bx = dat$beta_x, sy = dat$se_y, sx = dat$se_x)",
    "fit <- suppressWarnings(MRPRESSO::mr_presso(",
    "  BetaOutcome = \"by\", BetaExposure = \"bx\", SdOutcome = \"sy\", SdExposure = \"sx\", data = pd,",
    "  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,",
    "  SignifThreshold = signif, NbDistribution = nbdist, seed = seed))",
    "main <- fit[[\"Main MR results\"]]",
    "gt <- fit[[\"MR-PRESSO results\"]][[\"Global Test\"]]",
    "dt <- fit[[\"MR-PRESSO results\"]][[\"Distortion Test\"]]",
    "idx <- if (!is.null(dt) && is.numeric(dt[[\"Outliers Indices\"]])) dt[[\"Outliers Indices\"]] else integer(0)",
    "n_all <- nrow(dat); n_out <- length(idx)",
    "use_corrected <- n_out > 0 && is.finite(main[2, \"Causal Estimate\"])",
    "row <- if (use_corrected) 2 else 1   # the outlier-corrected row exists only when outliers were found",
    "df_t <- if (use_corrected) n_all - n_out - 1 else n_all - 1",
    "tq <- stats::qt(1 - (1 - ci_level) / 2, df_t)   # t quantile, as in the lm() fit MR-PRESSO uses",
    "b <- main[row, \"Causal Estimate\"]; se <- main[row, \"Sd\"]",
    "p_num <- function(x) suppressWarnings(as.numeric(sub(\"^<\", \"\", as.character(x))))",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = b, se = se, ci_lower = b - tq * se, ci_upper = b + tq * se,",
    "  ci_level = ci_level, ci_method = \"wald_t\", p = main[row, \"P-value\"], n_ivs_used = as.integer(n_all - n_out),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- rbind(",
    "  data.frame(name = \"presso_global_test\", value = gt$RSSobs, df = NA_real_, p = p_num(gt$Pvalue),",
    "             note = if (grepl(\"^<\", as.character(gt$Pvalue))) \"P is below the smallest value the simulation can give (1 / NbDistribution)\" else NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"presso_outliers\", value = n_out, df = NA_real_, p = NA_real_,",
    "             note = if (n_out) paste(dat$variant_id[idx], collapse = \";\") else NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"presso_distortion\", value = if (is.null(dt)) NA_real_ else as.double(dt[[\"Distortion Coefficient\"]]), df = NA_real_,",
    "             p = if (is.null(dt)) NA_real_ else p_num(dt$Pvalue), note = \"percent change of the estimate caused by removing the outliers\", stringsAsFactors = FALSE),",
    "  data.frame(name = \"presso_raw_estimate\", value = main[1, \"Causal Estimate\"], df = NA_real_, p = NA_real_,",
    "             note = if (use_corrected) \"the reported estimate is the outlier-corrected one\" else \"no outliers found: the reported estimate is the raw one\", stringsAsFactors = FALSE),",
    "  data.frame(name = \"presso_simulations\", value = as.double(nbdist), df = NA_real_, p = NA_real_, note = NA_character_, stringsAsFactors = FALSE))"),
    "fixed seed: the null distribution is simulated, so the same seed gives the same numbers"),
  column_requirements = .colreq_basic)

.presso_pick <- function(fit, ids, n_all, ci_level, nbdist) {
  main <- fit[["Main MR results"]]; gt <- fit[["MR-PRESSO results"]][["Global Test"]]; dt <- fit[["MR-PRESSO results"]][["Distortion Test"]]
  idx <- if (!is.null(dt) && is.numeric(dt[["Outliers Indices"]])) dt[["Outliers Indices"]] else integer(0)
  n_out <- length(idx); corr <- n_out > 0 && is.finite(main[2, "Causal Estimate"]); row <- if (corr) 2 else 1
  df_t <- if (corr) n_all - n_out - 1 else n_all - 1; tq <- stats::qt(1 - (1 - ci_level) / 2, df_t)
  b <- main[row, "Causal Estimate"]; se <- main[row, "Sd"]
  pn <- function(x) suppressWarnings(as.numeric(sub("^<", "", as.character(x))))
  list(b = b, se = se, lo = b - tq * se, hi = b + tq * se, p = main[row, "P-value"], n = n_all - n_out,
       diag = rbind(diagnostic_row("presso_global_test", gt$RSSobs, NA_real_, pn(gt$Pvalue), if (grepl("^<", as.character(gt$Pvalue))) "P is below the smallest value the simulation can give (1 / NbDistribution)" else NA_character_),
                    diagnostic_row("presso_outliers", n_out, NA_real_, NA_real_, .ids_note(ids[idx])),
                    diagnostic_row("presso_distortion", if (is.null(dt)) NA_real_ else as.double(dt[["Distortion Coefficient"]]), NA_real_, if (is.null(dt)) NA_real_ else pn(dt$Pvalue), "percent change of the estimate caused by removing the outliers"),
                    diagnostic_row("presso_raw_estimate", main[1, "Causal Estimate"], NA_real_, NA_real_, if (corr) "the reported estimate is the outlier-corrected one" else "no outliers found: the reported estimate is the raw one"),
                    diagnostic_row("presso_simulations", as.double(nbdist))))
}

# Radial IVW and radial Egger (RadialMR) --------------------------------------------------------------
.radial_q_p <- function(q, df) stats::pchisq(q, df, lower.tail = FALSE)
.radial_ids <- function(o) if (is.data.frame(o)) as.character(o$SNP) else character(0)

ADAPTER_RADIAL_V11 <- NULL
ADAPTER_RADIAL_V11 <- mrgui_adapter(
  id = "radial", adapter_version = "1.1.0", label = "Radial IVW", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(RadialMR = PIN_RADIALMR), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .radial_alpha_param_v11, .radial_weights_param_v11, .radial_tol_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Radial IVW needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    r <- RadialMR::format_radial(d$beta_x, d$beta_y, d$se_x, d$se_y, d$variant_id)
    fit <- suppressWarnings(RadialMR::ivw_radial(r, params$alpha, as.integer(params$weights), params$tol, FALSE))
    cf <- fit$coef[1, ]; tq <- stats::qt(1 - (1 - params$ci_level) / 2, fit$df); ids <- .radial_ids(fit$outliers)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           cf[[1]], cf[[2]], cf[[1]] - tq * cf[[2]], cf[[1]] + tq * cf[[2]], "wald_t", cf[[4]], nrow(d), params$ci_level),
         diagnostics = rbind(diagnostic_row("radial_q", fit$qstatistic, fit$df, .radial_q_p(fit$qstatistic, fit$df)),
                             diagnostic_row("radial_outliers", length(ids), NA_real_, NA_real_, .ids_note(ids)),
                             diagnostic_row("radial_alpha", params$alpha)))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_RADIAL_V11, spec, params, seed, "RadialMR", c(
    "# ---- 5. Analysis: radial IVW (RadialMR) ------------------------------------------",
    "r_in <- RadialMR::format_radial(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, dat$variant_id)",
    "fit <- suppressWarnings(RadialMR::ivw_radial(r_in, alpha, as.integer(weights), tol, FALSE))",
    "cf <- fit$coef[1, ]   # first row is the radial IVW effect for the chosen weights",
    "tq <- stats::qt(1 - (1 - ci_level) / 2, fit$df)   # t quantile, df = number of variants - 1",
    "ids <- if (is.data.frame(fit$outliers)) as.character(fit$outliers$SNP) else character(0)",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = cf[[1]], se = cf[[2]], ci_lower = cf[[1]] - tq * cf[[2]], ci_upper = cf[[1]] + tq * cf[[2]],",
    "  ci_level = ci_level, ci_method = \"wald_t\", p = cf[[4]], n_ivs_used = nrow(dat),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- rbind(",
    "  data.frame(name = \"radial_q\", value = fit$qstatistic, df = fit$df, p = stats::pchisq(fit$qstatistic, fit$df, lower.tail = FALSE), note = NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"radial_outliers\", value = length(ids), df = NA_real_, p = NA_real_, note = if (length(ids)) paste(ids, collapse = \";\") else NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"radial_alpha\", value = alpha, df = NA_real_, p = NA_real_, note = NA_character_, stringsAsFactors = FALSE))"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

ADAPTER_RADEGGER_V11 <- NULL
ADAPTER_RADEGGER_V11 <- mrgui_adapter(
  id = "radegger", adapter_version = "1.1.0", label = "Radial MR-Egger", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 4, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(RadialMR = PIN_RADIALMR), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .radial_alpha_param_v11, .radial_weights_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 4) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Radial MR-Egger needs at least 4 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    r <- RadialMR::format_radial(d$beta_x, d$beta_y, d$se_x, d$se_y, d$variant_id)
    fit <- suppressWarnings(RadialMR::egger_radial(r, params$alpha, as.integer(params$weights), FALSE))
    sl <- fit$coef[2, ]; ic <- fit$coef[1, ]; tq <- stats::qt(1 - (1 - params$ci_level) / 2, fit$df); ids <- .radial_ids(fit$outliers)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           sl[[1]], sl[[2]], sl[[1]] - tq * sl[[2]], sl[[1]] + tq * sl[[2]], "wald_t", sl[[4]], nrow(d), params$ci_level),
         diagnostics = rbind(diagnostic_row("radial_q", fit$qstatistic, fit$df, .radial_q_p(fit$qstatistic, fit$df)),
                             diagnostic_row("radial_outliers", length(ids), NA_real_, NA_real_, .ids_note(ids)),
                             diagnostic_row("radial_intercept", ic[[1]], NA_real_, ic[[4]], "intercept of the radial Egger regression, on the weighted scale (not the same unit as the MR-Egger intercept)"),
                             diagnostic_row("radial_alpha", params$alpha)))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_RADEGGER_V11, spec, params, seed, "RadialMR", c(
    "# ---- 5. Analysis: radial MR-Egger (RadialMR) -------------------------------------",
    "r_in <- RadialMR::format_radial(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, dat$variant_id)",
    "fit <- suppressWarnings(RadialMR::egger_radial(r_in, alpha, as.integer(weights), FALSE))",
    "sl <- fit$coef[2, ]; ic <- fit$coef[1, ]   # slope = causal effect; intercept = directional pleiotropy",
    "tq <- stats::qt(1 - (1 - ci_level) / 2, fit$df)   # t quantile, df = number of variants - 2",
    "ids <- if (is.data.frame(fit$outliers)) as.character(fit$outliers$SNP) else character(0)",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = sl[[1]], se = sl[[2]], ci_lower = sl[[1]] - tq * sl[[2]], ci_upper = sl[[1]] + tq * sl[[2]],",
    "  ci_level = ci_level, ci_method = \"wald_t\", p = sl[[4]], n_ivs_used = nrow(dat),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- rbind(",
    "  data.frame(name = \"radial_q\", value = fit$qstatistic, df = fit$df, p = stats::pchisq(fit$qstatistic, fit$df, lower.tail = FALSE), note = NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"radial_outliers\", value = length(ids), df = NA_real_, p = NA_real_, note = if (length(ids)) paste(ids, collapse = \";\") else NA_character_, stringsAsFactors = FALSE),",
    "  data.frame(name = \"radial_intercept\", value = ic[[1]], df = NA_real_, p = ic[[4]], note = \"intercept of the radial Egger regression, on the weighted scale (not the same unit as the MR-Egger intercept)\", stringsAsFactors = FALSE),",
    "  data.frame(name = \"radial_alpha\", value = alpha, df = NA_real_, p = NA_real_, note = NA_character_, stringsAsFactors = FALSE))"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

# ---- MR Steiger directionality test (Hemani, Tilling, Davey Smith 2017, PLoS Genet 13:e1007081). ----
# A diagnostic-only method (contracts v1.1, additive): it gives no causal estimate, only a check of the assumed direction
# (does the exposure explain more variance in the instruments than the outcome does?). In the app it is shown as part of IVW
# (Round 16), but it stays its own adapter and script so saved studies and their replay checks are unchanged.
# Continuous trait: variance explained from effect, SE and sample size (TwoSampleMR::get_r_from_bsen).
# Binary trait (Round 16): variance explained on the liability scale from the log odds ratio, the allele frequency, the numbers of
# cases and controls and the population prevalence (TwoSampleMR::get_r_from_lor, logit model; Lee et al. 2012, Genet Epidemiol
# 36:214, PMID 22714935). The case and control counts and the prevalence are written into the analysis data file as the columns
# ncase_*, ncontrol_*, prevalence_*; they are present only for a binary trait, and that presence is what the script reads.
.steiger_r_param <- function(id, label, what, arg) {
  setNames(list(list(
    type = "double", default = 1, min = 0.1, max = 1,
    label = label,
    hint = sprintf("Correlation between the measured and the true %s (its measurement reliability). 1 assumes no measurement error. Values below 1 re-run the test after correcting the variance explained for that error.", what),
    package_arg = arg, kind = "float", conventional = TRUE, advanced = FALSE,
    evidence = "1 is the default of TwoSampleMR::mr_steiger; Hemani et al. (2017) describe the sensitivity analysis over values below 1.",
    affects = c("direction"))), id)
}

.steiger_empty_estimates <- function() data.frame(
  estimand = character(0), exposure_id = character(0), outcome_id = character(0), scale = character(0),
  estimate = numeric(0), se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0), ci_level = numeric(0),
  ci_method = character(0), p = numeric(0), n_ivs_used = integer(0), stringsAsFactors = FALSE)

# Is this side of the pair binary in the analysis data (case and control counts present)?
.steiger_is_binary <- function(d, s) paste0("ncase_", s) %in% names(d)

# Signed correlation of each variant with one trait, and the sample size the test uses for that trait.
.steiger_side <- function(d, s) {
  b <- d[[paste0("beta_", s)]]; se <- d[[paste0("se_", s)]]
  if (.steiger_is_binary(d, s)) {
    nc <- d[[paste0("ncase_", s)]]; nn <- d[[paste0("ncontrol_", s)]]
    list(r = TwoSampleMR::get_r_from_lor(lor = b, af = d[[paste0("eaf_", s)]], ncase = nc, ncontrol = nn,
                                         prevalence = d[[paste0("prevalence_", s)]], model = "logit", correction = FALSE),
         n = nc + nn)
  } else list(r = TwoSampleMR::get_r_from_bsen(b, se, d[[paste0("n_", s)]]), n = d[[paste0("n_", s)]])
}

.steiger_scale_note <- function(d, s) {
  if (.steiger_is_binary(d, s)) sprintf("binary: liability scale (logit model), %s cases, %s controls, prevalence %s",
                                        format(d[[paste0("ncase_", s)]][1], scientific = FALSE), format(d[[paste0("ncontrol_", s)]][1], scientific = FALSE),
                                        format(d[[paste0("prevalence_", s)]][1], digits = 15))
  else "continuous: from effect, standard error and sample size"
}

.steiger_core <- function(d, r_xxo, r_yyo) {
  x <- .steiger_side(d, "x"); y <- .steiger_side(d, "y")
  st <- TwoSampleMR::mr_steiger(p_exp = rep(NA_real_, nrow(d)), p_out = rep(NA_real_, nrow(d)), n_exp = x$n, n_out = y$n,
                                r_exp = x$r, r_out = y$r, r_xxo = r_xxo, r_yyo = r_yyo)
  wrong <- abs(y$r) > abs(x$r)
  rbind(
    diagnostic_row("steiger_r2_exposure", st$r2_exp, NA_real_, NA_real_, "variance in the exposure explained by the variants (sum over variants)"),
    diagnostic_row("steiger_r2_outcome", st$r2_out, NA_real_, NA_real_, "variance in the outcome explained by the same variants (sum over variants)"),
    diagnostic_row("steiger_direction", as.double(st$correct_causal_direction), NA_real_, st$steiger_test,
                   "1 = the exposure explains more variance than the outcome (the assumed direction); 0 = the reverse. P tests whether the two differ"),
    diagnostic_row("steiger_direction_adjusted", as.double(st$correct_causal_direction_adj), NA_real_, st$steiger_test_adj,
                   "same test after correcting for measurement reliability (r_xxo, r_yyo); equals the unadjusted test when both are 1"),
    diagnostic_row("steiger_variants", nrow(d), NA_real_, NA_real_, "variants used in the test"),
    diagnostic_row("steiger_variants_wrong_way", sum(wrong), NA_real_, NA_real_, .ids_note(d$variant_id[wrong])),
    diagnostic_row("steiger_exposure_binary", as.double(.steiger_is_binary(d, "x")), NA_real_, NA_real_, .steiger_scale_note(d, "x")),
    diagnostic_row("steiger_outcome_binary", as.double(.steiger_is_binary(d, "y")), NA_real_, NA_real_, .steiger_scale_note(d, "y")))
}

# What a binary side needs, checked before running. Returns NULL when fine, else list(code, message).
.steiger_side_problem <- function(d, s, scale, what) {
  if (!identical(scale, "log_or")) {
    nn <- d[[paste0("n_", s)]]
    if (is.null(nn) || any(!is.finite(nn) | nn <= 2))
      return(list(code = "FEAS_STEIGER_NEEDS_N", msg = sprintf("The Steiger test needs the sample size of the %s for every variant, and some are missing.", what)))
    return(NULL)
  }
  col <- function(k) d[[paste0(k, "_", s)]]
  nc <- col("ncase"); nn <- col("ncontrol")
  if (is.null(nc) || is.null(nn) || any(!is.finite(nc) | !is.finite(nn) | nc < 1 | nn < 1))
    return(list(code = "FEAS_STEIGER_NEEDS_CASES", msg = sprintf("The %s is binary. The Steiger test then needs its numbers of cases and controls, and they are not known for this dataset.", what)))
  pv <- col("prevalence")
  if (is.null(pv) || any(!is.finite(pv) | pv <= 0 | pv >= 1))
    return(list(code = "FEAS_STEIGER_NEEDS_PREVALENCE", msg = sprintf("The %s is binary. The Steiger test then needs its prevalence in the population (a proportion between 0 and 1), and it has not been entered.", what)))
  af <- col("eaf")
  if (is.null(af) || any(!is.finite(af) | af <= 0 | af >= 1))
    return(list(code = "FEAS_STEIGER_NEEDS_EAF", msg = sprintf("The %s is binary. The Steiger test then needs the effect allele frequency of every variant in that dataset, and some are missing.", what)))
  NULL
}

.steiger_script_side <- function(s, what) c(
  sprintf("if (\"ncase_%s\" %%in%% names(dat)) {   # binary %s: case and control counts and prevalence are columns of the data file", s, what),
  sprintf("  if (any(!is.finite(dat$eaf_%s) | !is.finite(dat$ncase_%s) | !is.finite(dat$ncontrol_%s) | !is.finite(dat$prevalence_%s) | dat$prevalence_%s <= 0 | dat$prevalence_%s >= 1))", s, s, s, s, s, s),
  sprintf("    stop(\"The binary %s needs eaf_%s, ncase_%s, ncontrol_%s and prevalence_%s (between 0 and 1) for every variant.\")", what, s, s, s, s),
  sprintf("  r_%s <- TwoSampleMR::get_r_from_lor(lor = dat$beta_%s, af = dat$eaf_%s, ncase = dat$ncase_%s, ncontrol = dat$ncontrol_%s,", s, s, s, s, s),
  sprintf("    prevalence = dat$prevalence_%s, model = \"logit\", correction = FALSE)   # liability scale (Lee et al. 2012)", s),
  sprintf("  n_%s_used <- dat$ncase_%s + dat$ncontrol_%s", s, s, s),
  sprintf("  scale_note_%s <- sprintf(\"binary: liability scale (logit model), %%s cases, %%s controls, prevalence %%s\",", s),
  sprintf("    format(dat$ncase_%s[1], scientific = FALSE), format(dat$ncontrol_%s[1], scientific = FALSE), format(dat$prevalence_%s[1], digits = 15))", s, s, s),
  "} else {   # continuous: variance explained from effect, standard error and sample size",
  sprintf("  if (any(!is.finite(dat$n_%s) | dat$n_%s <= 2)) stop(\"The Steiger test needs the sample size of the %s for every variant (column n_%s).\")", s, s, what, s),
  sprintf("  r_%s <- TwoSampleMR::get_r_from_bsen(dat$beta_%s, dat$se_%s, dat$n_%s)", s, s, s, s),
  sprintf("  n_%s_used <- dat$n_%s", s, s),
  sprintf("  scale_note_%s <- \"continuous: from effect, standard error and sample size\"", s),
  "}")

ADAPTER_STEIGER_V11 <- NULL
ADAPTER_STEIGER_V11 <- mrgui_adapter(
  id = "steiger", adapter_version = "1.2.0", label = "MR Steiger directionality", input_kind = "pair_iv", estimands = character(0),
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = TRUE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.steiger_r_param("r_xxo", "Exposure measurement reliability (r_xxo)", "exposure", "r_xxo"),
             .steiger_r_param("r_yyo", "Outcome measurement reliability (r_yyo)", "outcome", "r_yyo"), .minf_param_v11),
  diagnostic_only = TRUE,
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("The Steiger test needs at least 3 usable instruments; %d available.", nrow(d))))
    for (s in list(list("x", if (is.null(input$exposure_scale)) "beta" else input$exposure_scale, "exposure"), list("y", input$outcome_scale, "outcome"))) {
      pr <- .steiger_side_problem(d, s[[1]], s[[2]], s[[3]])
      if (!is.null(pr)) return(verdict("block", pr$code, pr$msg))
    }
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    list(estimates = .steiger_empty_estimates(), diagnostics = .steiger_core(d, params$r_xxo, params$r_yyo))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_STEIGER_V11, spec, params, seed, "TwoSampleMR", c(
    "# ---- 5. Analysis: MR Steiger directionality test (TwoSampleMR::mr_steiger) --------",
    "# Does the exposure explain more variance in the instruments than the outcome does? If so the assumed direction is supported.",
    "# This method gives no causal estimate; it only checks the direction.",
    "# Variance explained: continuous trait from effect, SE and N; binary trait on the liability scale from the log odds ratio,",
    "# allele frequency, cases, controls and population prevalence (TwoSampleMR::get_r_from_lor; Lee et al. 2012, PMID 22714935).",
    .steiger_script_side("x", "exposure"),
    .steiger_script_side("y", "outcome"),
    "st <- TwoSampleMR::mr_steiger(p_exp = rep(NA_real_, nrow(dat)), p_out = rep(NA_real_, nrow(dat)),",
    "  n_exp = n_x_used, n_out = n_y_used, r_exp = r_x, r_out = r_y, r_xxo = r_xxo, r_yyo = r_yyo)",
    "wrong <- abs(r_y) > abs(r_x)   # variants that explain more variance in the outcome than in the exposure",
    "estimates <- data.frame(estimand = character(0), exposure_id = character(0), outcome_id = character(0), scale = character(0),",
    "  estimate = numeric(0), se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0), ci_level = numeric(0),",
    "  ci_method = character(0), p = numeric(0), n_ivs_used = integer(0), stringsAsFactors = FALSE)   # no causal estimate by design",
    "dg <- function(name, value, p, note) data.frame(name = name, value = as.double(value), df = NA_real_, p = as.double(p), note = as.character(note), stringsAsFactors = FALSE)",
    "diagnostics <- rbind(",
    "  dg(\"steiger_r2_exposure\", st$r2_exp, NA_real_, \"variance in the exposure explained by the variants (sum over variants)\"),",
    "  dg(\"steiger_r2_outcome\", st$r2_out, NA_real_, \"variance in the outcome explained by the same variants (sum over variants)\"),",
    "  dg(\"steiger_direction\", st$correct_causal_direction, st$steiger_test,",
    "     \"1 = the exposure explains more variance than the outcome (the assumed direction); 0 = the reverse. P tests whether the two differ\"),",
    "  dg(\"steiger_direction_adjusted\", st$correct_causal_direction_adj, st$steiger_test_adj,",
    "     \"same test after correcting for measurement reliability (r_xxo, r_yyo); equals the unadjusted test when both are 1\"),",
    "  dg(\"steiger_variants\", nrow(dat), NA_real_, \"variants used in the test\"),",
    "  dg(\"steiger_variants_wrong_way\", sum(wrong), NA_real_, if (any(wrong)) paste(dat$variant_id[wrong], collapse = \";\") else NA_character_),",
    "  dg(\"steiger_exposure_binary\", \"ncase_x\" %in% names(dat), NA_real_, scale_note_x),",
    "  dg(\"steiger_outcome_binary\", \"ncase_y\" %in% names(dat), NA_real_, scale_note_y))"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = list(
    list(column = "beta_x", level = "required", why = "Effect of each variant on the exposure (log odds ratio for a binary exposure).", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
    list(column = "se_x", level = "required", why = "Standard error of the exposure effect.", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
    list(column = "beta_y", level = "required", why = "Effect of each variant on the outcome (log odds ratio for a binary outcome).", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
    list(column = "se_y", level = "required", why = "Standard error of the outcome effect.", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
    list(column = "n_x", level = "required", why = "Sample size of the exposure study (continuous exposure).", if_missing = list(action = "block", code = "FEAS_STEIGER_NEEDS_N")),
    list(column = "n_y", level = "required", why = "Sample size of the outcome study (continuous outcome).", if_missing = list(action = "block", code = "FEAS_STEIGER_NEEDS_N")),
    list(column = "ncase_y", level = "recommended", why = "Number of cases, for a binary outcome (also ncontrol_y). Same for a binary exposure (ncase_x, ncontrol_x).", if_missing = list(action = "block", code = "FEAS_STEIGER_NEEDS_CASES")),
    list(column = "prevalence_y", level = "recommended", why = "Population prevalence, for a binary outcome (prevalence_x for a binary exposure). Entered by the analyst.", if_missing = list(action = "block", code = "FEAS_STEIGER_NEEDS_PREVALENCE")),
    list(column = "eaf_y", level = "recommended", why = "Effect allele frequency, for a binary outcome (eaf_x for a binary exposure).", if_missing = list(action = "block", code = "FEAS_STEIGER_NEEDS_EAF"))))

# ---- Round 37: debiased IVW (Ye, Shao and Kang 2021, Ann Stat 49:2079-2100, doi:10.1214/20-AOS2027) ------------------
# MendelianRandomization::mr_divw from CRAN (conda-forge r-mendelianrandomization), pinned by its released version.
PIN_MENDELIANRANDOMIZATION <- "0.10.0"
.divw_disp_param_v11 <- list(dispersion = list(
  type = "enum", default = "balanced", values = c("balanced", "none"),
  value_labels = c(balanced = "Allow for balanced pleiotropy (standard)", none = "No pleiotropy"),
  value_hints = c(balanced = "The standard error includes a term for balanced horizontal pleiotropy estimated from the data (over.dispersion = TRUE, the MendelianRandomization default).",
                  none = "Assumes every variant is a valid instrument, so the standard error is smaller. Use it only with a strong reason."),
  label = "Pleiotropy in the standard error",
  hint = "Whether the standard error allows for balanced horizontal pleiotropy. The estimate is the same either way; only the standard error, interval and P value change.",
  package_arg = "over.dispersion", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Ye et al. 2021 (Ann Stat, doi:10.1214/20-AOS2027); over.dispersion = TRUE is the default of MendelianRandomization::mr_divw.",
  affects = c("se", "ci", "p")))
DIVW_CONDITION_NOTE <- "(mean F - 1) x sqrt(number of variants); the package recommends more than 20 for a reliable interval and P value"

ADAPTER_DIVW_V11 <- NULL
ADAPTER_DIVW_V11 <- mrgui_adapter(
  id = "divw", adapter_version = "1.1.0", label = "Debiased IVW", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 2, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(MendelianRandomization = PIN_MENDELIANRANDOMIZATION), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .divw_disp_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 2) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Debiased IVW needs at least 2 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    fit <- MendelianRandomization::mr_divw(MendelianRandomization::mr_input(bx = d$beta_x, bxse = d$se_x, by = d$beta_y, byse = d$se_y),
                                           over.dispersion = identical(params$dispersion, "balanced"), alpha = 1 - params$ci_level)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, "wald", fit@Pvalue, as.integer(fit@SNPs), params$ci_level),
         diagnostics = diagnostic_row("divw_condition", fit@Condition, NA_real_, NA_real_, DIVW_CONDITION_NOTE))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_DIVW_V11, spec, params, seed, "MendelianRandomization", c(
    "# ---- 5. Analysis: debiased inverse-variance weighted estimate (MendelianRandomization) ----",
    "fit <- MendelianRandomization::mr_divw(",
    "  MendelianRandomization::mr_input(bx = dat$beta_x, bxse = dat$se_x, by = dat$beta_y, byse = dat$se_y),",
    "  over.dispersion = identical(dispersion, \"balanced\"),   # TRUE allows for balanced pleiotropy (package default)",
    "  alpha = 1 - ci_level)   # normal interval: estimate +/- z * standard error",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit@Estimate, se = fit@StdError, ci_lower = fit@CILower, ci_upper = fit@CIUpper,",
    "  ci_level = ci_level, ci_method = \"wald\", p = fit@Pvalue, n_ivs_used = as.integer(fit@SNPs),",
    "  stringsAsFactors = FALSE)",
    sprintf("diagnostics <- data.frame(name = \"divw_condition\", value = fit@Condition, df = NA_real_, p = NA_real_, note = %s,", .q(DIVW_CONDITION_NOTE)),
    "                          stringsAsFactors = FALSE)"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

# ---- Round 38: MR-RAPS (Zhao, Wang, Hemani, Bowden and Small 2020, Ann Stat 48:1742-1769, doi:10.1214/19-AOS1866) -----
# mr.raps 0.4.3 from GitHub (qingyuanzhao/mr.raps, GPL-3), pinned by commit like MR-PRESSO; its imports come from conda-forge.
# shrinkage = FALSE gives the published robust adjusted profile score estimator, as TwoSampleMR::mr_raps does by default.
PIN_MRRAPS <- "dd79b5bd74d10b699503cdd226d70c726f11c796"
.raps_disp_param_v11 <- list(raps_dispersion = list(
  type = "enum", default = "balanced", values = c("balanced", "none"),
  value_labels = c(balanced = "Allow for balanced pleiotropy (standard)", none = "No pleiotropy"),
  value_hints = c(balanced = "Each variant may have its own pleiotropic effect, normally distributed around zero; its variance is estimated (over.dispersion = TRUE, the default of mr.raps and TwoSampleMR).",
                  none = "Assumes every variant is a valid instrument (over.dispersion = FALSE). Use it only with a strong reason."),
  label = "Pleiotropy model",
  hint = "Whether the model allows for balanced horizontal pleiotropy. It changes the estimate a little and its standard error more.",
  package_arg = "over.dispersion", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Zhao et al. 2020 (Ann Stat, doi:10.1214/19-AOS1866); over.dispersion = TRUE is the default of mr.raps and TwoSampleMR.",
  affects = c("estimate", "se", "ci", "p")))
.raps_loss_param_v11 <- list(loss = list(
  type = "enum", default = "huber", values = c("huber", "tukey", "l2"),
  value_labels = c(huber = "Huber (standard)", tukey = "Tukey biweight", l2 = "Squared error (not robust)"),
  value_hints = c(huber = "Limits the influence of variants with large residuals (k = 1.345). Default of mr.raps and TwoSampleMR.",
                  tukey = "Gives variants with very large residuals no influence at all (k = 4.685); more robust to a few strong outliers.",
                  l2 = "Every variant keeps full influence: no protection against outlying variants."),
  label = "Robust loss function",
  hint = "How strongly variants with large residuals (possible outliers from pleiotropy) are down-weighted.",
  package_arg = "loss.function", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "Zhao et al. 2020; Huber loss is the default of mr.raps and TwoSampleMR::default_parameters().",
  affects = c("estimate", "se", "ci", "p")))

ADAPTER_RAPS_V11 <- NULL
ADAPTER_RAPS_V11 <- mrgui_adapter(
  id = "raps", adapter_version = "1.1.0", label = "MR-RAPS", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(mr.raps = PIN_MRRAPS), stochastic = TRUE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .raps_disp_param_v11, .raps_loss_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("MR-RAPS needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(ctx$seed))
    fit <- suppressWarnings(mr.raps::mr.raps(data.frame(beta.exposure = d$beta_x, beta.outcome = d$beta_y, se.exposure = d$se_x, se.outcome = d$se_y),
             diagnostics = FALSE, over.dispersion = identical(params$raps_dispersion, "balanced"), loss.function = params$loss, shrinkage = FALSE))
    z <- stats::qnorm(1 - (1 - params$ci_level) / 2)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           fit$beta.hat, fit$beta.se, fit$beta.hat - z * fit$beta.se, fit$beta.hat + z * fit$beta.se, "wald",
           2 * stats::pnorm(-abs(fit$beta.hat / fit$beta.se)), nrow(d), params$ci_level),
         diagnostics = diagnostic_row("raps_overdispersion", fit$tau2.hat, NA_real_, NA_real_, sprintf("estimated variance of the pleiotropic effects (SE %s)", format(signif(fit$tau2.se, 4)))))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_RAPS_V11, spec, params, seed, "mr.raps", c(
    "# ---- 5. Analysis: MR-RAPS, robust adjusted profile score (mr.raps) -------------",
    "fit <- suppressWarnings(mr.raps::mr.raps(",
    "  data.frame(beta.exposure = dat$beta_x, beta.outcome = dat$beta_y, se.exposure = dat$se_x, se.outcome = dat$se_y),",
    "  diagnostics = FALSE,",
    "  over.dispersion = identical(raps_dispersion, \"balanced\"),   # TRUE allows for balanced pleiotropy (package default)",
    "  loss.function = loss,                                   # robust loss: huber (default), tukey or l2",
    "  shrinkage = FALSE))                                     # the published RAPS estimator, as TwoSampleMR::mr_raps",
    "z <- stats::qnorm(1 - (1 - ci_level) / 2)   # normal interval, as in the paper",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit$beta.hat, se = fit$beta.se, ci_lower = fit$beta.hat - z * fit$beta.se, ci_upper = fit$beta.hat + z * fit$beta.se,",
    "  ci_level = ci_level, ci_method = \"wald\", p = 2 * stats::pnorm(-abs(fit$beta.hat / fit$beta.se)), n_ivs_used = nrow(dat),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- data.frame(name = \"raps_overdispersion\", value = fit$tau2.hat, df = NA_real_, p = NA_real_,",
    "  note = sprintf(\"estimated variance of the pleiotropic effects (SE %s)\", format(signif(fit$tau2.se, 4))), stringsAsFactors = FALSE)"),
    "fixed seed: mr.raps draws random starting values, so the same seed gives the same numbers"),
  column_requirements = .colreq_basic)

# ---- Round 39: leave-one-out and single-variant analyses (part of IVW; STROBE-MR item 13e) ----------------------------
# Diagnostic-only. Leave-one-out: the IVW estimate with each variant left out in turn, with the same IVW model and
# TwoSampleMR function as the primary analysis (setting ivw_model, shared with the IVW card). Single-variant: the Wald
# ratio of each variant (TwoSampleMR::mr_wald_ratio, first-order standard error), which the funnel plot shows.
# Every number is a diagnostics row keyed by the variant in `note`, so the replay check compares all of them.
.loo_core <- function(d, ivw_model) {
  fit_fun <- switch(ivw_model %||% "mre_floor", mre_floor = TwoSampleMR::mr_ivw, mre = TwoSampleMR::mr_ivw_mre, fixed = TwoSampleMR::mr_ivw_fe)
  k <- nrow(d); row <- function(n, v, p, id) data.frame(name = n, value = as.double(v), df = NA_real_, p = as.double(p), note = id, stringsAsFactors = FALSE)
  out <- list(row("loo_variants", k, NA_real_, NA_character_))
  for (i in seq_len(k)) {
    r <- fit_fun(d$beta_x[-i], d$beta_y[-i], d$se_x[-i], d$se_y[-i])
    w <- TwoSampleMR::mr_wald_ratio(d$beta_x[i], d$beta_y[i], d$se_x[i], d$se_y[i])
    out[[length(out) + 1L]] <- rbind(row("loo_estimate", r$b, r$pval, d$variant_id[i]), row("loo_se", r$se, NA_real_, d$variant_id[i]),
                                     row("wald_ratio", w$b, w$pval, d$variant_id[i]), row("wald_ratio_se", w$se, NA_real_, d$variant_id[i]))
  }
  do.call(rbind, out)
}

ADAPTER_LOO_V11 <- NULL
ADAPTER_LOO_V11 <- mrgui_adapter(
  id = "loo", adapter_version = "1.1.0", label = "Leave-one-out and single-variant analyses", input_kind = "pair_iv", estimands = character(0),
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ivw_model_param_v11, .minf_param_v11),
  diagnostic_only = TRUE,
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Leave-one-out analysis needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) list(estimates = .steiger_empty_estimates(), diagnostics = .loo_core(.usable(input$data), params$ivw_model)),
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_LOO_V11, spec, params, seed, "TwoSampleMR", c(
    "# ---- 5. Analysis: leave-one-out IVW and single-variant Wald ratios (TwoSampleMR) ----",
    "fit_fun <- switch(ivw_model,                # the same IVW model as the primary analysis",
    "  mre_floor = TwoSampleMR::mr_ivw, mre = TwoSampleMR::mr_ivw_mre, fixed = TwoSampleMR::mr_ivw_fe)",
    "row <- function(n, v, p, id) data.frame(name = n, value = as.double(v), df = NA_real_, p = as.double(p), note = id, stringsAsFactors = FALSE)",
    "out <- list(row(\"loo_variants\", nrow(dat), NA_real_, NA_character_))",
    "for (i in seq_len(nrow(dat))) {",
    "  r <- fit_fun(dat$beta_x[-i], dat$beta_y[-i], dat$se_x[-i], dat$se_y[-i])                     # IVW without variant i",
    "  w <- TwoSampleMR::mr_wald_ratio(dat$beta_x[i], dat$beta_y[i], dat$se_x[i], dat$se_y[i])        # variant i on its own",
    "  out[[length(out) + 1L]] <- rbind(row(\"loo_estimate\", r$b, r$pval, dat$variant_id[i]), row(\"loo_se\", r$se, NA_real_, dat$variant_id[i]),",
    "                                   row(\"wald_ratio\", w$b, w$pval, dat$variant_id[i]), row(\"wald_ratio_se\", w$se, NA_real_, dat$variant_id[i]))",
    "}",
    "diagnostics <- do.call(rbind, out)   # one row per variant and quantity; the variant is in 'note'",
    "estimates <- data.frame(estimand = character(0), exposure_id = character(0), outcome_id = character(0), scale = character(0),",
    "  estimate = numeric(0), se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0), ci_level = numeric(0),",
    "  ci_method = character(0), p = numeric(0), n_ivs_used = integer(0), stringsAsFactors = FALSE)   # diagnostic only"),
    "fixed seed (this analysis is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

# ---- Round 40: contamination mixture (Burgess et al. 2020, Nat Commun 11:376, PMID 31953392) and MR-Lasso ------------
# (Rees, Wood, Dudbridge and Burgess 2019, PLoS One 14:e0222362, PMID 31545794), both from MendelianRandomization 0.10.0.
.conmix_psi_param_v11 <- list(psi = list(
  type = "double", default = 0, min = 0, max = 1000,
  label = "Spread of the invalid estimates (psi)",
  hint = "Standard deviation assumed for the ratio estimates of invalid variants. 0 lets the package choose 1.5 times the standard deviation of all ratio estimates, which is its default.",
  package_arg = "psi", kind = "float", conventional = TRUE, advanced = TRUE,
  evidence = "MendelianRandomization::mr_conmix default (psi = 0); Burgess et al. 2020 (PMID 31953392).",
  affects = c("estimate", "ci", "p")))
.conmix_step_param_v11 <- list(ci_step = list(
  type = "double", default = 0.01, min = 1e-5, max = 1,
  label = "Grid step",
  hint = "The estimate and its confidence interval are found on a grid of candidate values with this spacing, so they are precise only to this step. 0.01 is the package default; choose a smaller step when the effect is small (for example below 0.1).",
  package_arg = "CIStep", kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "MendelianRandomization::mr_conmix default (CIStep = 0.01).",
  affects = c("estimate", "ci")))
.lasso_dist_param_v11 <- list(lasso_dist = list(
  type = "enum", default = "normal", values = c("normal", "t-dist"),
  value_labels = c(normal = "Normal distribution (standard)", "t-dist" = "t distribution"),
  value_hints = c(normal = "Confidence interval and P value from the normal distribution (MendelianRandomization default).",
                  "t-dist" = "Confidence interval and P value from a t distribution, a little wider with few variants."),
  label = "Distribution for the interval",
  hint = "Distribution used for the confidence interval and P value of the post-lasso estimate. The estimate itself does not change.",
  package_arg = "distribution", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "MendelianRandomization::mr_lasso default (distribution = \"normal\").",
  affects = c("ci", "p")))

.mr_input_line <- "MendelianRandomization::mr_input(bx = dat$beta_x, bxse = dat$se_x, by = dat$beta_y, byse = dat$se_y, snps = dat$variant_id)"
.mrin <- function(d) MendelianRandomization::mr_input(bx = d$beta_x, bxse = d$se_x, by = d$beta_y, byse = d$se_y, snps = d$variant_id)
.mr_pkg_check <- function(label, min) function(input, params) {
  d <- .usable(input$data)
  if (nrow(d) < min) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("%s needs at least %d usable instruments; %d available.", label, min, nrow(d))))
  w <- .common_warnings(d); verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
}
.conmix_diag <- function(fit, ids) {
  nr <- length(fit@CILower); inv <- setdiff(ids, ids[fit@Valid])
  rbind(diagnostic_row("conmix_psi", fit@Psi, NA_real_, NA_real_, "standard deviation of the invalid estimates used"),
        diagnostic_row("conmix_invalid", length(inv), NA_real_, NA_real_, if (length(inv)) paste(inv, collapse = ";") else NA_character_),
        diagnostic_row("conmix_ci_ranges", nr, NA_real_, NA_real_, if (nr > 1) paste(sprintf("%.6g to %.6g", fit@CILower, fit@CIUpper), collapse = "; ") else "a single range of values"),
        diagnostic_row("conmix_grid_step", fit@CIStep, NA_real_, NA_real_, "the estimate and interval are precise to this step"))
}
ADAPTER_CONMIX_V11 <- NULL
ADAPTER_CONMIX_V11 <- mrgui_adapter(
  id = "conmix", adapter_version = "1.1.0", label = "Contamination mixture", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(MendelianRandomization = PIN_MENDELIANRANDOMIZATION), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .conmix_psi_param_v11, .conmix_step_param_v11, .minf_param_v11),
  check = .mr_pkg_check("The contamination mixture method", 3L),
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    fit <- MendelianRandomization::mr_conmix(.mrin(d), psi = params$psi, CIStep = params$ci_step, alpha = 1 - params$ci_level)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           fit@Estimate, NA_real_, min(fit@CILower), max(fit@CIUpper), "profile", fit@Pvalue, as.integer(fit@SNPs), params$ci_level),
         diagnostics = .conmix_diag(fit, d$variant_id))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_CONMIX_V11, spec, params, seed, "MendelianRandomization", c(
    "# ---- 5. Analysis: contamination mixture method (MendelianRandomization) --------",
    sprintf("fit <- MendelianRandomization::mr_conmix(%s,", .mr_input_line),
    "  psi = psi,            # 0: package default, 1.5 x SD of the ratio estimates",
    "  CIStep = ci_step,     # grid step: estimate and interval are precise to this step",
    "  alpha = 1 - ci_level)",
    "# The method gives no standard error. Its profile-likelihood interval can be several separate ranges;",
    "# the outer limits are reported here and every range is listed in the diagnostics.",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit@Estimate, se = NA_real_, ci_lower = min(fit@CILower), ci_upper = max(fit@CIUpper),",
    "  ci_level = ci_level, ci_method = \"profile\", p = fit@Pvalue, n_ivs_used = as.integer(fit@SNPs),",
    "  stringsAsFactors = FALSE)",
    "nr <- length(fit@CILower); inv <- setdiff(dat$variant_id, dat$variant_id[fit@Valid])",
    "diagnostics <- data.frame(name = c(\"conmix_psi\", \"conmix_invalid\", \"conmix_ci_ranges\", \"conmix_grid_step\"),",
    "  value = c(fit@Psi, length(inv), nr, fit@CIStep), df = NA_real_, p = NA_real_,",
    "  note = c(\"standard deviation of the invalid estimates used\", if (length(inv)) paste(inv, collapse = \";\") else NA_character_,",
    "           if (nr > 1) paste(sprintf(\"%.6g to %.6g\", fit@CILower, fit@CIUpper), collapse = \"; \") else \"a single range of values\",",
    "           \"the estimate and interval are precise to this step\"), stringsAsFactors = FALSE)"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

ADAPTER_LASSO_V11 <- NULL
ADAPTER_LASSO_V11 <- mrgui_adapter(
  id = "lasso", adapter_version = "1.1.0", label = "MR-Lasso", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(MendelianRandomization = PIN_MENDELIANRANDOMIZATION), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .lasso_dist_param_v11, .minf_param_v11),
  check = .mr_pkg_check("MR-Lasso", 3L),
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    fit <- MendelianRandomization::mr_lasso(.mrin(d), distribution = params$lasso_dist, alpha = 1 - params$ci_level)
    inv <- setdiff(d$variant_id, fit@ValidSNPs)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, if (identical(params$lasso_dist, "t-dist")) "wald_t" else "wald", fit@Pvalue, as.integer(fit@Valid), params$ci_level),
         diagnostics = rbind(diagnostic_row("lasso_excluded", length(inv), NA_real_, NA_real_, if (length(inv)) paste(inv, collapse = ";") else NA_character_),
                             diagnostic_row("lasso_lambda", fit@Lambda, NA_real_, NA_real_, "tuning parameter chosen by the package (heterogeneity stopping rule)"),
                             diagnostic_row("lasso_regression_estimate", fit@RegEstimate, NA_real_, NA_real_, "penalised regression estimate before the post-lasso IVW step")))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_LASSO_V11, spec, params, seed, "MendelianRandomization", c(
    "# ---- 5. Analysis: MR-Lasso (MendelianRandomization) ---------------------------",
    sprintf("fit <- MendelianRandomization::mr_lasso(%s,", .mr_input_line),
    "  distribution = lasso_dist,   # interval and P value: normal (default) or t-dist",
    "  alpha = 1 - ci_level)",
    "# The reported estimate is the post-lasso IVW estimate on the variants the lasso kept as valid.",
    "inv <- setdiff(dat$variant_id, fit@ValidSNPs)",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = fit@Estimate, se = fit@StdError, ci_lower = fit@CILower, ci_upper = fit@CIUpper,",
    "  ci_level = ci_level, ci_method = if (identical(lasso_dist, \"t-dist\")) \"wald_t\" else \"wald\", p = fit@Pvalue, n_ivs_used = as.integer(fit@Valid),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- data.frame(name = c(\"lasso_excluded\", \"lasso_lambda\", \"lasso_regression_estimate\"),",
    "  value = c(length(inv), fit@Lambda, fit@RegEstimate), df = NA_real_, p = NA_real_,",
    "  note = c(if (length(inv)) paste(inv, collapse = \";\") else NA_character_, \"tuning parameter chosen by the package (heterogeneity stopping rule)\",",
    "           \"penalised regression estimate before the post-lasso IVW step\"), stringsAsFactors = FALSE)"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)
