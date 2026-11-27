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

.radial_level <- function(params, n) if (identical(params$alpha_correction, "none")) params$alpha else params$alpha / n
.radial_level_note <- function(params, n) if (identical(params$alpha_correction, "none")) "no correction" else sprintf("%s divided by %d variants", format(params$alpha), n)

.radial_alpha_param_v11 <- list(alpha = list(
  type = "double", default = 0.05, min = 0.0001, max = 0.5,
  label = "Outlier significance level (alpha)",
  hint = "A variant is flagged when its contribution to the Q statistic is significant at this level. The package applies it per variant with no multiplicity correction; the setting below divides it by the number of variants (Bonferroni) by default.",
  package_arg = "alpha", kind = "float", conventional = TRUE, advanced = FALSE,
  evidence = "0.05 is the example value in the RadialMR documentation; Bowden et al. (2018) discuss a Bonferroni-corrected level.",
  affects = c("outliers")))

.radial_correction_param_v11 <- list(alpha_correction = list(
  type = "enum", default = "bonferroni", values = c("bonferroni", "none"),
  value_labels = c(bonferroni = "Divide by the number of variants", none = "None"),
  value_hints = c(bonferroni = "Each variant is tested at alpha divided by the number of variants (Bonferroni), so about 5% of analyses flag a valid variant by chance.",
                  none = "Each variant is tested at alpha itself, as the package does. With many variants, valid ones are flagged by chance (about 5% of them)."),
  label = "Correction for the number of variants",
  hint = "The level at which each variant is tested. MR Studio passes alpha divided by the number of variants to the package unless you choose None.",
  package_arg = "alpha", kind = "enum", conventional = FALSE, advanced = FALSE,
  evidence = "Bowden et al. (2018) discuss a Bonferroni-corrected level. MR Studio calibration (Round 32, 48,000 simulated data sets): at 0.05 radial detectors flagged about 5% of valid variants and lowered IVW coverage after removal to 0.90 to 0.92; at 0.05 divided by the number of variants under 0.5%, with coverage 0.94 to 0.95.",
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
  if (!length(a$packages)) { h[grepl("^# Package  :", h)] <- "# Package  : none beyond base R (MR Studio's own implementation of the published algorithm, below)"; return(h) }
  pins <- unlist(a$packages); kind <- ifelse(grepl("^[0-9]+([.-][0-9]+)+$", pins), "version", "commit")   # CRAN versions may contain a hyphen (boot 1.3-32)   # Round 37: CRAN packages are pinned by version
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

# Round 110a: mv = TRUE (MVMR GRAPPLE) writes the multivariable data head: the table as saved, no single-exposure F filter (its min_f
# setting does not exist for a multivariable method, so the script stopped with "object 'min_f' not found"); the body keeps the variants
# with every exposure and the outcome, as .mv_core.
.detector_emit <- function(a, spec, params, seed, pkg, body, seed_note, mv = FALSE) {
  data_lines <- if (mv) c("# ---- 4. Data (the exact table the calculation uses) ----------------------------", "dat <- utils::read.csv(DATA_FILE, stringsAsFactors = FALSE)", "") else c(.script_lines_data, .script_lines_minf(a$requires$min_ivs))
  txt <- c(.script_lines_head_pk(a, spec), .script_lines_guard_pk(spec, pkg),
           "# ---- 2. Settings (every adjustable argument of this method) --------------------",
           .param_comment_lines(a, params),
           sprintf("seed <- %d   # %s", as.integer(seed), seed_note),
           "", .script_lines_repro(seed), data_lines,
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
  params = c(.ci_param_v11, .radial_alpha_param_v11, .radial_correction_param_v11, .radial_weights_param_v11, .radial_tol_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Radial IVW needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    r <- RadialMR::format_radial(d$beta_x, d$beta_y, d$se_x, d$se_y, d$variant_id)
    au <- .radial_level(params, nrow(d))
    fit <- suppressWarnings(RadialMR::ivw_radial(r, au, as.integer(params$weights), params$tol, FALSE))
    cf <- fit$coef[1, ]; tq <- stats::qt(1 - (1 - params$ci_level) / 2, fit$df); ids <- .radial_ids(fit$outliers)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           cf[[1]], cf[[2]], cf[[1]] - tq * cf[[2]], cf[[1]] + tq * cf[[2]], "wald_t", cf[[4]], nrow(d), params$ci_level),
         diagnostics = rbind(diagnostic_row("radial_q", fit$qstatistic, fit$df, .radial_q_p(fit$qstatistic, fit$df)),
                             diagnostic_row("radial_outliers", length(ids), NA_real_, NA_real_, .ids_note(ids)),
                             diagnostic_row("radial_alpha", au, NA_real_, NA_real_, .radial_level_note(params, nrow(d)))))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_RADIAL_V11, spec, params, seed, "RadialMR", c(
    "# ---- 5. Analysis: radial IVW (RadialMR) ------------------------------------------",
    "r_in <- RadialMR::format_radial(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, dat$variant_id)",
    "alpha_used <- if (identical(alpha_correction, \"bonferroni\")) alpha / nrow(dat) else alpha   # level per variant",
    "fit <- suppressWarnings(RadialMR::ivw_radial(r_in, alpha_used, as.integer(weights), tol, FALSE))",
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
    "  data.frame(name = \"radial_alpha\", value = alpha_used, df = NA_real_, p = NA_real_, note = if (identical(alpha_correction, \"bonferroni\")) sprintf(\"%s divided by %d variants\", format(alpha), nrow(dat)) else \"no correction\", stringsAsFactors = FALSE))"),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)

ADAPTER_RADEGGER_V11 <- NULL
ADAPTER_RADEGGER_V11 <- mrgui_adapter(
  id = "radegger", adapter_version = "1.1.0", label = "Radial MR-Egger", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 4, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(RadialMR = PIN_RADIALMR), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .radial_alpha_param_v11, .radial_correction_param_v11, .radial_weights_param_v11, .minf_param_v11),
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 4) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Radial MR-Egger needs at least 4 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    r <- RadialMR::format_radial(d$beta_x, d$beta_y, d$se_x, d$se_y, d$variant_id)
    au <- .radial_level(params, nrow(d))
    fit <- suppressWarnings(RadialMR::egger_radial(r, au, as.integer(params$weights), FALSE))
    sl <- fit$coef[2, ]; ic <- fit$coef[1, ]; tq <- stats::qt(1 - (1 - params$ci_level) / 2, fit$df); ids <- .radial_ids(fit$outliers)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
           sl[[1]], sl[[2]], sl[[1]] - tq * sl[[2]], sl[[1]] + tq * sl[[2]], "wald_t", sl[[4]], nrow(d), params$ci_level),
         diagnostics = rbind(diagnostic_row("radial_q", fit$qstatistic, fit$df, .radial_q_p(fit$qstatistic, fit$df)),
                             diagnostic_row("radial_outliers", length(ids), NA_real_, NA_real_, .ids_note(ids)),
                             diagnostic_row("radial_intercept", ic[[1]], NA_real_, ic[[4]], "intercept of the radial Egger regression, on the weighted scale (not the same unit as the MR-Egger intercept)"),
                             diagnostic_row("radial_alpha", au, NA_real_, NA_real_, .radial_level_note(params, nrow(d)))))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_RADEGGER_V11, spec, params, seed, "RadialMR", c(
    "# ---- 5. Analysis: radial MR-Egger (RadialMR) -------------------------------------",
    "r_in <- RadialMR::format_radial(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, dat$variant_id)",
    "alpha_used <- if (identical(alpha_correction, \"bonferroni\")) alpha / nrow(dat) else alpha   # level per variant",
    "fit <- suppressWarnings(RadialMR::egger_radial(r_in, alpha_used, as.integer(weights), FALSE))",
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
    "  data.frame(name = \"radial_alpha\", value = alpha_used, df = NA_real_, p = NA_real_, note = if (identical(alpha_correction, \"bonferroni\")) sprintf(\"%s divided by %d variants\", format(alpha), nrow(dat)) else \"no correction\", stringsAsFactors = FALSE))"),
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
# Round 43 (FR-105, Burgess et al. guidelines): the share of the total IVW weight carried by the strongest variant (first-order
# weights beta_x^2 / se_y^2; its leave-one-out estimate is the estimate without it), and, with 10 or more variants, repeated
# leave-out of random subsets (a fraction at a time, a fixed number of repeats, drawn with the script's seed) summarised by
# quantiles of the IVW estimate. The 10-variant floor, 30% and 200 repeats are [design] defaults.
.lso_frac_param_v11 <- list(lso_fraction = list(
  type = "double", default = 0.3, min = 0.05, max = 0.5,
  label = "Leave-subset-out: share left out each time",
  hint = "Share of the variants left out at random in each repeat of the leave-subset-out check (runs with 10 or more variants).",
  package_arg = NA_character_, kind = "float", conventional = FALSE, advanced = TRUE,
  evidence = "Burgess et al. guidelines (PMID 32760811) suggest leaving out subsets when there are many variants; 30% is an MR Studio design default.",
  affects = c("diagnostics")))
.lso_rep_param_v11 <- list(lso_repeats = list(
  type = "integer", default = 200L, min = 50L, max = 5000L,
  label = "Leave-subset-out: repeats",
  hint = "Number of random subsets left out. More repeats give steadier quantiles and take longer.",
  package_arg = NA_character_, kind = "integer", conventional = FALSE, advanced = TRUE,
  evidence = "MR Studio design default (200), drawn with the fixed seed written in the script.",
  affects = c("diagnostics")))

.loo_core <- function(d, ivw_model, lso_fraction = 0.3, lso_repeats = 200L) {
  fit_fun <- switch(ivw_model %||% "mre_floor", mre_floor = TwoSampleMR::mr_ivw, mre = TwoSampleMR::mr_ivw_mre, fixed = TwoSampleMR::mr_ivw_fe)
  k <- nrow(d); row <- function(n, v, p, id) data.frame(name = n, value = as.double(v), df = NA_real_, p = as.double(p), note = id, stringsAsFactors = FALSE)
  out <- list(row("loo_variants", k, NA_real_, NA_character_))
  for (i in seq_len(k)) {
    r <- fit_fun(d$beta_x[-i], d$beta_y[-i], d$se_x[-i], d$se_y[-i])
    w <- TwoSampleMR::mr_wald_ratio(d$beta_x[i], d$beta_y[i], d$se_x[i], d$se_y[i])
    out[[length(out) + 1L]] <- rbind(row("loo_estimate", r$b, r$pval, d$variant_id[i]), row("loo_se", r$se, NA_real_, d$variant_id[i]),
                                     row("wald_ratio", w$b, w$pval, d$variant_id[i]), row("wald_ratio_se", w$se, NA_real_, d$variant_id[i]))
  }
  wt <- d$beta_x^2 / d$se_y^2; top <- which.max(wt)
  out[[length(out) + 1L]] <- row("dominant_share", wt[top] / sum(wt), NA_real_, d$variant_id[top])
  if (k >= 10) {
    m <- max(1L, floor(lso_fraction * k)); b <- vapply(seq_len(lso_repeats), function(j) { o <- sample.int(k, m); fit_fun(d$beta_x[-o], d$beta_y[-o], d$se_x[-o], d$se_y[-o])$b }, 0)
    q <- stats::quantile(b, c(0.025, 0.5, 0.975), names = FALSE)
    out[[length(out) + 1L]] <- rbind(row("lso_repeats", lso_repeats, NA_real_, sprintf("%d of %d variants left out each time", m, k)),
      row("lso_q025", q[1], NA_real_, NA_character_), row("lso_median", q[2], NA_real_, NA_character_), row("lso_q975", q[3], NA_real_, NA_character_),
      row("lso_min", min(b), NA_real_, NA_character_), row("lso_max", max(b), NA_real_, NA_character_))
  }
  do.call(rbind, out)
}

ADAPTER_LOO_V11 <- NULL
ADAPTER_LOO_V11 <- mrgui_adapter(
  id = "loo", adapter_version = "1.1.0", label = "Leave-one-out and single-variant analyses", input_kind = "pair_iv", estimands = character(0),
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"), stochastic = TRUE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ivw_model_param_v11, .lso_frac_param_v11, .lso_rep_param_v11, .minf_param_v11),
  diagnostic_only = TRUE,
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("Leave-one-out analysis needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(ctx$seed))
    list(estimates = .steiger_empty_estimates(), diagnostics = .loo_core(.usable(input$data), params$ivw_model, params$lso_fraction, as.integer(params$lso_repeats))) },
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
    "k <- nrow(dat); wt <- dat$beta_x^2 / dat$se_y^2; top <- which.max(wt)   # first-order IVW weights",
    "out[[length(out) + 1L]] <- row(\"dominant_share\", wt[top] / sum(wt), NA_real_, dat$variant_id[top])   # strongest variant's share of the weight",
    "if (k >= 10) {   # leave-subset-out: random subsets left out, drawn with the seed above",
    "  m <- max(1L, floor(lso_fraction * k))",
    "  b <- vapply(seq_len(lso_repeats), function(j) { o <- sample.int(k, m); fit_fun(dat$beta_x[-o], dat$beta_y[-o], dat$se_x[-o], dat$se_y[-o])$b }, 0)",
    "  q <- stats::quantile(b, c(0.025, 0.5, 0.975), names = FALSE)",
    "  out[[length(out) + 1L]] <- rbind(row(\"lso_repeats\", lso_repeats, NA_real_, sprintf(\"%d of %d variants left out each time\", m, k)),",
    "    row(\"lso_q025\", q[1], NA_real_, NA_character_), row(\"lso_median\", q[2], NA_real_, NA_character_), row(\"lso_q975\", q[3], NA_real_, NA_character_),",
    "    row(\"lso_min\", min(b), NA_real_, NA_character_), row(\"lso_max\", max(b), NA_real_, NA_character_))",
    "}",
    "diagnostics <- do.call(rbind, out)   # one row per variant and quantity; the variant is in 'note'",
    "estimates <- data.frame(estimand = character(0), exposure_id = character(0), outcome_id = character(0), scale = character(0),",
    "  estimate = numeric(0), se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0), ci_level = numeric(0),",
    "  ci_method = character(0), p = numeric(0), n_ivs_used = integer(0), stringsAsFactors = FALSE)   # diagnostic only"),
    "fixed seed: the leave-subset-out check draws its random subsets with this seed"),
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

# ---- Round 41: MR-Egger with SIMEX correction (Bowden et al. 2016, Int J Epidemiol 45:1961-1974, PMID 27616674) ------
# Variants oriented to a positive exposure effect, MR-Egger as a regression weighted by 1/se_y^2, and simulation
# extrapolation (simex package, quadratic extrapolation) on the exposure effects with their standard errors as the
# measurement error, as in the paper. simex 1.8 (CRAN, GPL-3, base R only) is installed from its CRAN source file, checked by SHA-256.
PIN_SIMEX <- "1.8"
.simex_b_param_v11 <- list(simex_b = list(
  type = "integer", default = 100L, min = 50L, max = 10000L,
  label = "Simulations per step",
  hint = "Number of simulated data sets for each amount of added measurement error. 100 is the default of the simex package; more gives a smoother extrapolation and takes longer.",
  package_arg = "B", kind = "int", conventional = TRUE, advanced = TRUE,
  evidence = "simex::simex default (B = 100); Bowden et al. 2016 (PMID 27616674) used the simex package with quadratic extrapolation.",
  affects = c("estimate", "se", "ci", "p")))
.simex_core <- function(d, simex_b, ci_level) {
  eg <- data.frame(BXG = abs(d$beta_x), BYG = d$beta_y * sign(d$beta_x)); k <- nrow(d)
  # simex refits the model with update(), which looks up the weights in the global environment only (as in a script
  # run at top level); the weights are placed there for the call and removed afterwards.
  assign(".mrgui_simex_w", 1 / d$se_y^2, envir = globalenv()); on.exit(rm(".mrgui_simex_w", envir = globalenv()), add = TRUE)
  fit <- eval(quote(stats::lm(BYG ~ BXG, data = eg, weights = .mrgui_simex_w, x = TRUE, y = TRUE)))
  sim <- simex::simex(fit, SIMEXvariable = "BXG", measurement.error = d$se_x, B = simex_b, fitting.method = "quadratic", jackknife.estimation = "quadratic", asymptotic = FALSE)
  cf <- summary(sim)$coefficients$jackknife; tq <- stats::qt(1 - (1 - ci_level) / 2, k - 2)
  list(b = cf["BXG", 1], se = cf["BXG", 2], p = cf["BXG", 4], lo = cf["BXG", 1] - tq * cf["BXG", 2], hi = cf["BXG", 1] + tq * cf["BXG", 2], k = k,
       b0 = unname(stats::coef(fit)["BXG"]), int = cf["(Intercept)", 1], int_p = cf["(Intercept)", 4])
}
ADAPTER_EGGER_SIMEX_V11 <- NULL
ADAPTER_EGGER_SIMEX_V11 <- mrgui_adapter(
  id = "egger_simex", adapter_version = "1.1.0", label = "MR-Egger with SIMEX correction", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(simex = PIN_SIMEX), stochastic = TRUE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .simex_b_param_v11, .minf_param_v11),
  check = .mr_pkg_check("MR-Egger with SIMEX correction", 3L),
  run = function(input, params, ctx) {
    d <- .usable(input$data); RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(as.integer(ctx$seed))
    r <- .simex_core(d, params$simex_b, params$ci_level)
    list(estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale, r$b, r$se, r$lo, r$hi, "wald_t", r$p, as.integer(r$k), params$ci_level),
         diagnostics = rbind(diagnostic_row("simex_uncorrected_slope", r$b0, NA_real_, NA_real_, "weighted MR-Egger slope before the correction"),
                             diagnostic_row("simex_intercept", r$int, NA_real_, r$int_p, "SIMEX-corrected MR-Egger intercept (directional pleiotropy)")))
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_EGGER_SIMEX_V11, spec, params, seed, "simex", c(
    "# ---- 5. Analysis: MR-Egger with SIMEX correction (Bowden et al. 2016) ---------",
    "eg <- data.frame(BXG = abs(dat$beta_x), BYG = dat$beta_y * sign(dat$beta_x))   # variants oriented to a positive exposure effect",
    "w <- 1 / dat$se_y^2   # weights; kept outside 'eg' because simex refits the model with update() at top level",
    "fit <- stats::lm(BYG ~ BXG, data = eg, weights = w, x = TRUE, y = TRUE)   # weighted MR-Egger regression",
    "sim <- simex::simex(fit, SIMEXvariable = \"BXG\", measurement.error = dat$se_x, B = simex_b,",
    "  fitting.method = \"quadratic\", jackknife.estimation = \"quadratic\", asymptotic = FALSE)",
    "cf <- summary(sim)$coefficients$jackknife  # jackknife standard errors",
    "tq <- stats::qt(1 - (1 - ci_level) / 2, nrow(dat) - 2)   # t interval on k - 2 degrees of freedom, as MR-Egger",
    "estimates <- data.frame(estimand = \"causal_effect\",",
    paste0("  ", .est_cols(spec)),
    "  estimate = cf[\"BXG\", 1], se = cf[\"BXG\", 2], ci_lower = cf[\"BXG\", 1] - tq * cf[\"BXG\", 2], ci_upper = cf[\"BXG\", 1] + tq * cf[\"BXG\", 2],",
    "  ci_level = ci_level, ci_method = \"wald_t\", p = cf[\"BXG\", 4], n_ivs_used = nrow(dat),",
    "  stringsAsFactors = FALSE)",
    "diagnostics <- data.frame(name = c(\"simex_uncorrected_slope\", \"simex_intercept\"),",
    "  value = c(unname(stats::coef(fit)[\"BXG\"]), cf[\"(Intercept)\", 1]), df = NA_real_, p = c(NA_real_, cf[\"(Intercept)\", 4]),",
    "  note = c(\"weighted MR-Egger slope before the correction\", \"SIMEX-corrected MR-Egger intercept (directional pleiotropy)\"), stringsAsFactors = FALSE)"),
    "fixed seed: SIMEX adds simulated measurement error, so the same seed gives the same numbers"),
  column_requirements = .colreq_basic)

# ---- Round 50: multivariable MR (MVMR) from MendelianRandomization 0.10.0 (contracts input kind multi_iv) --------------------
# Input data: one beta_x_<k> / se_x_<k> (and optional n_x_<k>) column pair per exposure, plus beta_y, se_y and mr_keep; exposure_id
# holds the exposure IDs joined by "|". Each method reports one direct_effect row per exposure. The analysis code below is written
# once and used both by run() (evaluated in a clean environment) and by the generated script, so the replay compares the same code.
.mv_core <- c(
  "K <- sum(grepl(\"^beta_x_[0-9]+$\", names(dat)))   # number of exposures",
  "bxn <- sprintf(\"beta_x_%d\", seq_len(K)); sxn <- sprintf(\"se_x_%d\", seq_len(K))",
  "use <- dat$mr_keep & is.finite(dat$beta_y) & is.finite(dat$se_y) & dat$se_y > 0 & stats::complete.cases(dat[, c(bxn, sxn)])",
  "dat <- dat[use, , drop = FALSE]   # variants with every exposure and the outcome",
  "EXPOSURE_IDS <- strsplit(EXPOSURE_ID, \"|\", fixed = TRUE)[[1]]",
  "inp <- MendelianRandomization::mr_mvinput(bx = as.matrix(dat[, bxn]), bxse = as.matrix(dat[, sxn]), by = dat$beta_y, byse = dat$se_y,",
  "  snps = dat$variant_id, exposure = EXPOSURE_IDS, outcome = OUTCOME_ID)",
  "est_rows <- function(est, se, lo, hi, p, n) data.frame(estimand = \"direct_effect\", exposure_id = EXPOSURE_IDS, outcome_id = OUTCOME_ID, scale = SCALE,",
  "  estimate = as.numeric(est), se = as.numeric(se), ci_lower = as.numeric(lo), ci_upper = as.numeric(hi), ci_level = ci_level, ci_method = \"wald\",",
  "  p = as.numeric(p), n_ivs_used = as.integer(n), stringsAsFactors = FALSE)",
  "drow <- function(name, value, df = NA_real_, p = NA_real_, note = NA_character_) data.frame(name = name, value = as.numeric(value), df = as.numeric(df), p = as.numeric(p), note = as.character(note), stringsAsFactors = FALSE)")
# Round 53: the correlation matrix of the exposures, shared by MVMR-IVW (conditional F, Q_A) and the Q-minimisation estimator
.mv_cov_lines <- c("# Covariance of the variant-exposure associations of each variant (C[l, i, j]). Either computed from individual-level data",
  "# (columns cov_x_i_j, Round 55, as MVMR::snpcov_mvmr) or the phenotypic correlation of the exposures (mv_r_i_j) times the product of",
  "# their standard errors (as MVMR::phenocov_mvmr). All correlations 0 means non-overlapping exposure samples.",
  "bx <- as.matrix(dat[, bxn]); sx <- as.matrix(dat[, sxn]); J <- nrow(bx)",
  "cvn <- outer(seq_len(K), seq_len(K), function(i, j) sprintf(\"cov_x_%d_%d\", pmin(i, j), pmax(i, j)))",
  "R <- diag(K); for (i in seq_len(K - 1)) for (j in (i + 1):K) R[i, j] <- R[j, i] <- get(sprintf(\"mv_r_%d_%d\", i, j))",
  "if (all(cvn %in% names(dat))) { C <- array(as.numeric(as.matrix(dat[, as.vector(cvn)])), c(J, K, K)); COV_SOURCE <- \"individual-level data\"",
  "} else { C <- array(0, c(J, K, K)); for (i in seq_len(K)) for (j in seq_len(K)) C[, i, j] <- R[i, j] * sx[, i] * sx[, j]; COV_SOURCE <- \"phenotypic correlation\" }",
  "vq <- function(w, rows = seq_len(J)) as.vector(matrix(C[rows, , , drop = FALSE], length(rows)) %*% as.vector(outer(w, w)))   # t(w) %*% C[l, , ] %*% w for every variant",
  "cov_note <- if (COV_SOURCE == \"phenotypic correlation\") paste(sprintf(\"%d-%d %s\", row(R)[upper.tri(R)], col(R)[upper.tri(R)], format(R[upper.tri(R)])), collapse = \"; \") else COV_SOURCE")
.mv_body <- list(
  mvivw = c("# MVMR-IVW",
    "fit <- MendelianRandomization::mr_mvivw(inp, model = mv_model, alpha = 1 - ci_level)",
    "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
    .mv_cov_lines,
    "# Conditional F statistic of Sanderson, Spiller and Bowden 2021 for each exposure (as MVMR::strength_mvmr); values below 10 indicate weak instruments",
    "cond_f <- vapply(seq_len(K), function(k) { d <- stats::lm.fit(bx[, -k, drop = FALSE], bx[, k])$coefficients; w <- numeric(K); w[k] <- -1; w[-k] <- d",
    "  sum((bx[, k] - as.vector(bx[, -k, drop = FALSE] %*% d))^2 / vq(w)) / (J - K + 1) }, 0)",
    "# Q_A of the same paper: heterogeneity at the MVMR-IVW estimate allowing for the uncertainty of the exposure associations (as MVMR::pleiotropy_mvmr)",
    "b <- as.vector(fit@Estimate); q_a <- sum((dat$beta_y - as.vector(bx %*% b))^2 / (dat$se_y^2 + vq(b)))",
    "diagnostics <- rbind(drow(\"mv_cond_f\", cond_f, note = EXPOSURE_IDS),",
    "  drow(\"mv_q\", fit@Heter.Stat[1], fit@SNPs - K, fit@Heter.Stat[2], \"Cochran's Q of the MVMR-IVW fit\"),",
    "  drow(\"mv_q_a\", q_a, J - K, stats::pchisq(q_a, J - K, lower.tail = FALSE), \"Q_A (Sanderson et al. 2021), allowing for exposure uncertainty\"), drow(\"mv_model\", NA, note = fit@Model),",
    "  drow(\"mv_pheno_cor\", as.numeric(COV_SOURCE == \"individual-level data\" || any(R[upper.tri(R)] != 0)), note = cov_note))"),
  mvegger = c("# MVMR-Egger, variants oriented to a positive effect on the exposure given by mv_orient",
    "fit <- MendelianRandomization::mr_mvegger(inp, orientate = mv_orient, alpha = 1 - ci_level)",
    "estimates <- rbind(est_rows(fit@Estimate, fit@StdError.Est, fit@CILower.Est, fit@CIUpper.Est, fit@Pvalue.Est, fit@SNPs),",
    "  data.frame(estimand = \"pleiotropy_intercept\", exposure_id = EXPOSURE_ID, outcome_id = OUTCOME_ID, scale = SCALE, estimate = fit@Intercept, se = fit@StdError.Int,",
    "    ci_lower = fit@CILower.Int, ci_upper = fit@CIUpper.Int, ci_level = ci_level, ci_method = \"wald\", p = fit@Pvalue.Int, n_ivs_used = as.integer(fit@SNPs), stringsAsFactors = FALSE))",
    "diagnostics <- rbind(drow(\"mv_q\", fit@Heter.Stat[1], fit@SNPs - K - 1, fit@Heter.Stat[2], \"residual heterogeneity of the MVMR-Egger fit\"), drow(\"mv_orient\", mv_orient, note = EXPOSURE_IDS[mv_orient]))"),
  mvqhet = c("# MVMR Q-minimisation (Sanderson, Spiller and Bowden 2021), as MVMR::qhet_mvmr 0.4.9: the direct effects minimise Q_A, which allows",
    "# for the uncertainty (and correlation) of the variant-exposure associations, so they are robust to weak instruments. If the minimum Q_A",
    "# exceeds its degrees of freedom (J - K), an additive heterogeneity variance tau2 is added to every weight until it equals them.",
    .mv_cov_lines,
    "qfit <- function(rows) {",
    "  bx <- bx[rows, , drop = FALSE]; by <- dat$beta_y[rows]; sy <- dat$se_y[rows]",
    "  Qs <- function(b, tau2) sum((by - as.vector(bx %*% b))^2 / (sy^2 + vq(b, rows) + tau2))",
    "  b0 <- stats::lm.wfit(bx, by, 1 / sy^2)$coefficients   # start at MVMR-IVW",
    "  qm <- function(tau2) stats::optim(b0, Qs, tau2 = tau2, method = \"BFGS\", control = list(reltol = 1e-12))",
    "  target <- nrow(bx) - ncol(bx); tau2 <- 0",
    "  if (qm(0)$value > target) { up <- stats::median(sy^2); while (qm(up)$value > target) up <- up * 10",
    "    tau2 <- stats::uniroot(function(t) qm(t)$value - target, interval = c(0, up), tol = up * 1e-8)$root }",
    "  f <- qm(tau2); list(b = unname(f$par), tau2 = tau2, q = f$value) }",
    "est <- qfit(seq_len(J)); z <- stats::qnorm(1 - (1 - ci_level) / 2)",
    "if (mv_qhet_ci == \"jackknife\") {   # leave one variant out at a time; deterministic",
    "  jk <- matrix(vapply(seq_len(J), function(l) qfit(seq_len(J)[-l])$b, numeric(K)), nrow = K)",
    "  se <- sqrt((J - 1) / J * rowSums((jk - rowMeans(jk))^2)); lo <- est$b - z * se; hi <- est$b + z * se; cim <- \"wald\"",
    "} else {   # nonparametric bootstrap over variants with bias-corrected and accelerated (BCa) intervals, as qhet_mvmr(CI = TRUE)",
    "  bt <- boot::boot(data = data.frame(variant = seq_len(J)), statistic = function(d, ix) qfit(ix)$b, R = mv_qhet_iter)",
    "  se <- apply(bt$t, 2, stats::sd); ci <- vapply(seq_len(K), function(k) boot::boot.ci(bt, conf = ci_level, type = \"bca\", index = k)$bca[4:5], numeric(2))",
    "  lo <- ci[1, ]; hi <- ci[2, ]; cim <- \"bootstrap\" }",
    "estimates <- est_rows(est$b, se, lo, hi, 2 * stats::pnorm(-abs(est$b / se)), J); estimates$ci_method <- cim   # P from the jackknife or bootstrap standard error",
    "diagnostics <- rbind(drow(\"mv_qhet_tau2\", est$tau2, note = \"additive heterogeneity variance (0 when the minimum Q_A is within its degrees of freedom)\"),",
    "  drow(\"mv_qhet_q\", est$q, J - K, stats::pchisq(est$q, J - K, lower.tail = FALSE), \"Q_A at the estimate\"),",
    "  drow(\"mv_qhet_ci\", if (mv_qhet_ci == \"bootstrap\") mv_qhet_iter else J, note = mv_qhet_ci),",
    "  drow(\"mv_pheno_cor\", as.numeric(COV_SOURCE == \"individual-level data\" || any(R[upper.tri(R)] != 0)), note = cov_note))"),
  mvcml = c("# MVMR-cML (Lin, Xue and Pan 2023): constrained maximum likelihood with BIC selection of invalid variants, optionally data perturbation",
    "n_cml <- min(c(vapply(seq_len(K), function(k) stats::median(dat[[sprintf(\"n_x_%d\", k)]]), 0), stats::median(dat$n_y)))   # smallest GWAS sample size",
    "fit <- MendelianRandomization::mr_mvcML(inp, n = n_cml, DP = identical(mv_cml_dp, \"dp\"), num_pert = mv_cml_pert, seed = seed, alpha = 1 - ci_level)",
    "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
    "inv <- if (length(fit@BIC_invalid)) dat$variant_id[fit@BIC_invalid] else character()",
    "diagnostics <- rbind(drow(\"mv_cml_invalid\", length(inv), note = if (length(inv)) paste(inv, collapse = \";\") else NA_character_), drow(\"mv_cml_n\", n_cml, note = \"sample size used\"))"),
  mvmedian = c("# MVMR median (weighted), bootstrap standard errors with the script's seed",
    "fit <- MendelianRandomization::mr_mvmedian(inp, alpha = 1 - ci_level, iterations = mv_iter, seed = seed)",
    "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
    "diagnostics <- drow(\"mv_bootstrap\", mv_iter, note = \"bootstrap iterations\")"),
  mvlasso = c("# MVMR-Lasso: penalised intercepts identify invalid variants; the estimate is the post-lasso MVMR-IVW on the valid ones",
    "fit <- MendelianRandomization::mr_mvlasso(inp, orientate = mv_orient, alpha = 1 - ci_level)",
    "inv <- setdiff(dat$variant_id, fit@ValidSNPs)",
    "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@Valid)",
    "diagnostics <- rbind(drow(\"mv_lasso_excluded\", length(inv), note = if (length(inv)) paste(inv, collapse = \";\") else NA_character_), drow(\"mv_lasso_lambda\", fit@Lambda))"))
.mv_model_param <- list(mv_model = list(type = "enum", default = "default", values = c("default", "random", "fixed"),
  value_labels = c(default = "Package default (random effects with more than 3 variants)", random = "Multiplicative random effects", fixed = "Fixed effect"),
  label = "Effects model", hint = "How heterogeneity between variants enters the standard errors. The estimates do not change.",
  package_arg = "model", kind = "enum", conventional = TRUE, advanced = FALSE, evidence = "MendelianRandomization::mr_mvivw default (model = \"default\").", affects = c("se", "ci", "p")))
.mv_orient_param <- list(mv_orient = list(type = "integer", default = 1L, min = 1L, max = 6L, label = "Orient variants to exposure",
  hint = "MVMR-Egger and MVMR-Lasso orient the variants so that their effect on this exposure (1 = the first listed) is positive; the intercept depends on this choice.",
  package_arg = "orientate", kind = "integer", conventional = TRUE, advanced = TRUE, evidence = "MendelianRandomization default (orientate = 1); Rees et al. 2017.", affects = c("estimate", "se", "ci", "p")))
# Round 53: phenotypic correlation of each pair of exposures (up to 6 exposures, so 15 pairs), for the conditional F and Q_A
.mv_cor_params <- local({ out <- list(); for (i in 1:5) for (j in (i + 1):6) out[[sprintf("mv_r_%d_%d", i, j)]] <- list(type = "double", default = 0, min = -0.99, max = 0.99,
  label = sprintf("Phenotypic correlation, exposures %d and %d", i, j),
  hint = sprintf("Correlation between exposures %d and %d when their GWAS share participants (for example lipids from one GWAS); 0 for non-overlapping samples. Used only for the conditional F and Q_A, not the estimates.", i, j),
  package_arg = "pcor (MVMR::phenocov_mvmr)", kind = "double", conventional = TRUE, advanced = TRUE,
  evidence = "Sanderson, Spiller and Bowden 2021; MVMR package (phenocov_mvmr). Default 0 (zero covariance), as strength_mvmr(gencov = 0).", affects = c("diagnostics"))
  out })
.mv_cor_ok <- function(params, K) { R <- diag(K); if (K > 1) for (i in seq_len(K - 1)) for (j in (i + 1):K) R[i, j] <- R[j, i] <- as.numeric(params[[sprintf("mv_r_%d_%d", i, j)]] %||% 0)
  min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 1e-8 }
PIN_BOOT <- "1.3-32"   # Round 54: boot, for the bootstrap BCa interval of the Q-minimisation estimate (as MVMR::qhet_mvmr); the conda-forge build for R 4.5 (noarch), written as on CRAN
.mv_qhet_params <- list(mv_qhet_ci = list(type = "enum", default = "jackknife", values = c("jackknife", "bootstrap"),
    value_labels = c(jackknife = "Jackknife standard error (leave one variant out; deterministic)", bootstrap = "Bootstrap BCa interval (slow: one fit per iteration)"),
    label = "Uncertainty", hint = "How the standard error and interval of the Q-minimisation estimate are obtained. MVMR::qhet_mvmr offers both; its default when an interval is asked for is the bootstrap, which takes many minutes, so MR Studio defaults to the jackknife.",
    package_arg = "CI_method", kind = "enum", conventional = FALSE, advanced = FALSE, evidence = "MVMR::qhet_mvmr 0.4.9 (CI_method = bootstrap or jackknife).", affects = c("se", "ci", "p")),
  mv_qhet_iter = list(type = "integer", default = 1000L, min = 50L, max = 10000L, label = "Bootstrap iterations", hint = "Used only with the bootstrap.",
    package_arg = "iterations", kind = "integer", conventional = TRUE, advanced = TRUE, evidence = "MVMR::qhet_mvmr default (iterations = 1000).", affects = c("se", "ci", "p")))
.mv_iter_param <- list(mv_iter = list(type = "integer", default = 10000L, min = 500L, max = 100000L, label = "Bootstrap iterations",
  hint = "Iterations for the bootstrap standard error of the MVMR median estimate.", package_arg = "iterations", kind = "integer", conventional = TRUE, advanced = TRUE,
  evidence = "MendelianRandomization::mr_mvmedian default (iterations = 10000).", affects = c("se", "ci", "p")))
.mv_colreq <- list(
  list(column = "beta_x_<k>", level = "required", why = "Effect of each variant on exposure k (one column per exposure).", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
  list(column = "se_x_<k>", level = "required", why = "Standard error of each exposure effect.", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
  list(column = "beta_y", level = "required", why = "Effect of each variant on the outcome.", if_missing = list(action = "block", code = "IMPORT_NO_BETA")),
  list(column = "se_y", level = "required", why = "Standard error of the outcome effect.", if_missing = list(action = "block", code = "IMPORT_NO_SE")),
  list(column = "cov_x_<i>_<j>", level = "optional", why = "Covariance of the effects of each variant on exposures i and j (i <= j), from individual-level data; replaces the phenotypic correlations for the conditional F, Q_A and Q-minimisation.", if_missing = list(action = "ignore", code = NA_character_)),
  list(column = "n_x_<k>", level = "optional", why = "Sample size of each exposure, for the record.", if_missing = list(action = "ignore", code = NA_character_)))
.mv_usable <- function(d) { K <- sum(grepl("^beta_x_[0-9]+$", names(d))); bxn <- sprintf("beta_x_%d", seq_len(K)); sxn <- sprintf("se_x_%d", seq_len(K))
  d[d$mr_keep & is.finite(d$beta_y) & is.finite(d$se_y) & d$se_y > 0 & stats::complete.cases(d[, c(bxn, sxn)]), , drop = FALSE] }
.mv_adapter <- function(id, label, extra_k, params, stochastic, body_id) {
  a <- NULL
  a <- mrgui_adapter(id = id, adapter_version = "1.1.0", label = label, input_kind = "multi_iv", estimands = c("direct_effect", if (id == "mvegger") "pleiotropy_intercept"),
    requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 2, max_exposures = 6),
    packages = c(list(MendelianRandomization = PIN_MENDELIANRANDOMIZATION), if (id == "mvqhet") list(boot = PIN_BOOT)), stochastic = stochastic, supports_checkpoint = FALSE,
    runtime_class = if (id == "mvqhet") "moderate" else "fast",
    params = c(.ci_param_v11, params),
    check = function(input, params) {
      K <- sum(grepl("^beta_x_[0-9]+$", names(input$data))); n <- nrow(.mv_usable(input$data))
      if (K < 2) return(verdict("block", "FEAS_TOO_FEW_EXPOSURES", sprintf("%s needs at least 2 exposures; %d given.", label, K)))
      if (n < K + extra_k) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("%s needs at least %d variants with all %d exposures (exposures + %d); %d available.", label, K + extra_k, K, extra_k, n)))
      if (id %in% c("mvcml", "mvgmm") && !(all(sprintf("n_x_%d", seq_len(K)) %in% names(input$data)) && "n_y" %in% names(input$data) && all(is.finite(unlist(.mv_usable(input$data)[, c(sprintf("n_x_%d", seq_len(K)), "n_y")])))))
        return(verdict("block", "FEAS_NEEDS_N", sprintf("%s needs the sample size of every exposure GWAS and of the outcome GWAS.", label)))
      if (id %in% c("mvivw", "mvqhet", "mvivwme", "mvgmm") && !any(grepl("^cov_x_", names(input$data))) && !.mv_cor_ok(params, K)) return(verdict("block", "FEAS_INVALID_CORRELATION", "The phenotypic correlations of the exposures do not form a valid correlation matrix (it must be positive definite). Check the values."))
      verdict("ok", character(), character())
    },
    run = function(input, params, ctx) {
      e <- new.env(parent = globalenv()); e$dat <- input$data; for (k in names(params)) assign(k, params[[k]], envir = e)
      e$seed <- as.integer(ctx$seed); e$EXPOSURE_ID <- input$exposure_id; e$OUTCOME_ID <- input$outcome_id; e$SCALE <- input$outcome_scale
      RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(e$seed)
      eval(parse(text = c(.mv_core, .mv_body[[body_id]])), envir = e)
      list(estimates = e$estimates, diagnostics = e$diagnostics)
    },
    emit_script = function(spec, params, seed) {
      txt <- c(.script_lines_head_pk(a, spec), .script_lines_guard_pk(spec, "MendelianRandomization"),
        "# ---- 2. Settings (every adjustable argument of this method) --------------------", .param_comment_lines(a, params),
        sprintf("seed <- %d   # %s", as.integer(seed), if (stochastic) "fixed seed: the bootstrap uses this seed" else "fixed seed (this method is deterministic; kept for the record)"),
        sprintf("EXPOSURE_ID <- %s   # exposures, joined by |", .q(spec$exposure_id)), sprintf("OUTCOME_ID <- %s", .q(spec$outcome_id)), sprintf("SCALE <- %s", .q(spec$outcome_scale)),
        "", .script_lines_repro(seed),
        "# ---- 4. Data (the exact table the calculation uses) ----------------------------",
        "dat <- utils::read.csv(DATA_FILE, stringsAsFactors = FALSE)", .mv_core, "",
        sprintf("# ---- 5. Analysis: %s (MendelianRandomization) ----", label), .mv_body[[body_id]], "", .script_lines_write(spec))
      script_bundle(paste0(paste(txt, collapse = "\n"), "\n"), a$packages,
                    list(list(file = file.path("04_tables", spec$estimates_file), role = "estimates"), list(file = file.path("04_tables", spec$diagnostics_file), role = "diagnostics")),
                    seed, spec$data_file)
    },
    column_requirements = .mv_colreq)
  a
}
ADAPTER_MVIVW_V11 <- .mv_adapter("mvivw", "MVMR-IVW", 1L, c(.mv_model_param, .mv_cor_params), FALSE, "mvivw")
ADAPTER_MVEGGER_V11 <- .mv_adapter("mvegger", "MVMR-Egger", 2L, .mv_orient_param, FALSE, "mvegger")
ADAPTER_MVQHET_V11 <- .mv_adapter("mvqhet", "MVMR Q-minimisation", 2L, c(.mv_qhet_params, .mv_cor_params), TRUE, "mvqhet")
.mv_cml_params <- list(mv_cml_dp = list(type = "enum", default = "dp", values = c("dp", "ma"), value_labels = c(dp = "With data perturbation (MVMR-cML-DP, recommended)", ma = "Without data perturbation (MVMR-cML-BIC)"),
    label = "Variant", hint = "Data perturbation allows for uncertainty in which variants are selected as invalid.", package_arg = "DP", kind = "enum", conventional = TRUE, advanced = FALSE,
    evidence = "Lin et al. 2023; MendelianRandomization::mr_mvcML default (DP = TRUE).", affects = c("estimate", "se", "ci", "p")),
  mv_cml_pert = list(type = "integer", default = 200L, min = 50L, max = 2000L, label = "Perturbations", hint = "Number of data perturbations (DP only).", package_arg = "num_pert", kind = "integer",
    conventional = TRUE, advanced = TRUE, evidence = "MendelianRandomization::mr_mvcML default (num_pert = 200).", affects = c("se", "ci", "p")))
ADAPTER_MVCML_V11 <- .mv_adapter("mvcml", "MVMR-cML", 2L, .mv_cml_params, TRUE, "mvcml")
ADAPTER_MVMEDIAN_V11 <- .mv_adapter("mvmedian", "MVMR median", 2L, .mv_iter_param, TRUE, "mvmedian")
ADAPTER_MVLASSO_V11 <- .mv_adapter("mvlasso", "MVMR-Lasso", 2L, .mv_orient_param, FALSE, "mvlasso")

# ---- Round 58: more MendelianRandomization estimators ---------------------------------------------------------------------
# One body per method, run by the app (evaluated in an environment) and written into the script, so the two cannot drift apart.
# The body sees dat (usable variants after the F filter), the settings, seed, ci_level and the .est_cols values, and sets
# estimates and diagnostics.
.mr_body_common <- c(
  "mri <- if (requireNamespace(\"MendelianRandomization\", quietly = TRUE)) MendelianRandomization::mr_input(bx = dat$beta_x, bxse = dat$se_x, by = dat$beta_y, byse = dat$se_y) else NULL",
  "est1 <- function(b, se, lo, hi, p, n, method = \"wald\") data.frame(estimand = \"causal_effect\", exposure_id = EXPOSURE_ID, outcome_id = OUTCOME_ID, scale = SCALE,",
  "  estimate = as.numeric(b), se = as.numeric(se), ci_lower = as.numeric(lo), ci_upper = as.numeric(hi), ci_level = ci_level, ci_method = method, p = as.numeric(p),",
  "  n_ivs_used = as.integer(n), stringsAsFactors = FALSE)",
  "drow <- function(name, value, df = NA_real_, p = NA_real_, note = NA_character_) data.frame(name = name, value = as.numeric(value), df = as.numeric(df), p = as.numeric(p),",
  "  note = as.character(note), stringsAsFactors = FALSE)")
.mrpkg_adapter <- function(id, label, min_ivs, params, stochastic, body, needs_n = FALSE, check_extra = NULL, pkg = "MendelianRandomization", estimands = "causal_effect") {
  a <- NULL
  a <- mrgui_adapter(id = id, adapter_version = "1.1.0", label = label, input_kind = "pair_iv", estimands = estimands,
    requires = list(min_ivs = min_ivs, needs_ld = "none", needs_eaf = FALSE, needs_n = needs_n, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
    packages = switch(pkg, stats = stats::setNames(list(), character(0)), mrclust = list(mrclust = PIN_MRCLUST_V11), TwoSampleMR = list(TwoSampleMR = PIN_TWOSAMPLEMR_V11), MRMix = list(MRMix = PIN_MRMIX_V11), list(MendelianRandomization = PIN_MENDELIANRANDOMIZATION)), stochastic = stochastic, supports_checkpoint = FALSE, runtime_class = "fast",
    params = c(.ci_param_v11, params, .minf_param_v11),
    check = function(input, params) {
      d <- .usable(input$data)
      if (nrow(d) < min_ivs) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("%s needs at least %d usable instruments; %d available.", label, min_ivs, nrow(d))))
      if (!is.null(check_extra)) { v <- if (length(formals(check_extra)) >= 3L) check_extra(d, params, input) else check_extra(d, params); if (!is.null(v)) return(v) }
      w <- .common_warnings(d); verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
    },
    run = function(input, params, ctx) {
      e <- new.env(parent = globalenv()); for (k in names(params)) assign(k, params[[k]], envir = e)
      d <- .usable(input$data); d <- d[(d$beta_x / d$se_x)^2 >= params$min_f, , drop = FALSE]
      e$dat <- d; e$seed <- as.integer(ctx$seed); e$EXPOSURE_ID <- input$exposure_id; e$OUTCOME_ID <- input$outcome_id; e$SCALE <- input$outcome_scale; e$XSCALE <- input$exposure_scale %||% "beta"
      RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(e$seed)
      eval(parse(text = c(.mr_body_common, body)), envir = e)
      list(estimates = e$estimates, diagnostics = e$diagnostics)
    },
    emit_script = function(spec, params, seed) .detector_emit(a, spec, params, seed, pkg, c(
      sprintf("EXPOSURE_ID <- %s; OUTCOME_ID <- %s; SCALE <- %s; XSCALE <- %s", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale), .q(spec$exposure_scale %||% "beta")),
      sprintf("# ---- 5. Analysis: %s (%s) ----", label, pkg), "RNGkind(\"Mersenne-Twister\", \"Inversion\", \"Rejection\"); set.seed(seed)",
      .mr_body_common, body),
      if (stochastic) "fixed seed: the method draws random numbers with this seed" else "fixed seed (this method is deterministic; kept for the record)"),
    column_requirements = .colreq_basic)
  a
}
.enum_p <- function(name, default, values, labels, label, hint, arg, evidence, affects = c("estimate", "se", "ci", "p"), advanced = FALSE)
  setNames(list(list(type = "enum", default = default, values = values, value_labels = labels, label = label, hint = hint, package_arg = arg, kind = "enum",
                     conventional = TRUE, advanced = advanced, evidence = evidence, affects = affects)), name)
.num_p <- function(name, type, default, min, max, label, hint, arg, evidence, affects = c("estimate", "se", "ci", "p"), advanced = TRUE)
  setNames(list(list(type = type, default = default, min = min, max = max, label = label, hint = hint, package_arg = arg, kind = type, conventional = TRUE,
                     advanced = advanced, evidence = evidence, affects = affects)), name)

# Two-sample maximum likelihood (Burgess, Butterworth and Thompson 2013): allows for uncertainty in the variant-exposure associations.
ADAPTER_MAXLIK_V11 <- .mrpkg_adapter("maxlik", "Maximum likelihood", 3L, c(
  .enum_p("ml_model", "default", c("default", "random", "fixed"), c(default = "Package default (random effects with more than 3 variants)", random = "Random effects", fixed = "Fixed effect"),
          "Effects model", "How heterogeneity enters the standard error.", "model", "MendelianRandomization::mr_maxlik default (model = \"default\").", c("se", "ci", "p")),
  .num_p("ml_psi", "double", 0, -0.99, 0.99, "Correlation between the exposure and outcome associations (psi)",
         "Set it when the exposure and outcome GWAS share participants; 0 for separate samples.", "psi", "MendelianRandomization::mr_maxlik default (psi = 0).")), FALSE, c(
  "fit <- MendelianRandomization::mr_maxlik(mri, model = ml_model, psi = ml_psi, alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
  "diagnostics <- rbind(drow(\"ml_heterogeneity\", fit@Heter.Stat[1], fit@SNPs - 1, fit@Heter.Stat[2], \"likelihood-ratio heterogeneity test\"), drow(\"ml_model\", NA, note = fit@Model))"))

# Mode-based estimate (Hartwig, Davey Smith and Bowden 2017), MendelianRandomization implementation.
ADAPTER_MBE_V11 <- .mrpkg_adapter("mbe", "Mode-based estimate (MendelianRandomization)", 3L, c(
  .enum_p("mbe_weighting", "weighted", c("weighted", "simple"), c(weighted = "Weighted mode", simple = "Simple mode"), "Weighting",
          "Weighted by the inverse variance of each variant's ratio estimate, or equal weights.", "weighting", "MendelianRandomization::mr_mbe default (weighting = \"weighted\")."),
  .enum_p("mbe_stderror", "delta", c("delta", "simple"), c(delta = "Second-order (delta) weights (NOME not assumed)", simple = "First-order weights (NOME assumed)"), "Ratio standard errors",
          "Whether the ratio standard errors allow for uncertainty in the variant-exposure associations.", "stderror", "MendelianRandomization::mr_mbe default (stderror = \"delta\")."),
  .num_p("mbe_phi", "double", 1, 0.1, 5, "Bandwidth multiplier (phi)", "Smaller values give a narrower kernel and a less biased but more variable estimate.", "phi",
         "MendelianRandomization::mr_mbe default (phi = 1)."),
  .num_p("mbe_iter", "integer", 10000L, 500L, 100000L, "Bootstrap iterations", "Iterations for the bootstrap standard error.", "iterations", "MendelianRandomization::mr_mbe default (iterations = 10000).", c("se", "ci", "p"))), TRUE, c(
  "fit <- MendelianRandomization::mr_mbe(mri, weighting = if (identical(mbe_weighting, \"simple\")) \"unweighted\" else mbe_weighting, stderror = mbe_stderror, phi = mbe_phi, seed = seed, iterations = mbe_iter, alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
  "diagnostics <- drow(\"mbe_bootstrap\", mbe_iter, note = \"bootstrap iterations\")"))

# MR-cML (Xue, Shen and Pan 2021): constrained maximum likelihood, model-averaged by BIC, optionally with data perturbation (DP).
.cml_n_check <- function(d, params) if (!all(c("n_x", "n_y") %in% names(d)) || !all(is.finite(d$n_x) & is.finite(d$n_y))) verdict("block", "FEAS_NEEDS_N", "MR-cML needs the sample sizes of the exposure and outcome GWAS (n_x, n_y).") else NULL
ADAPTER_CML_V11 <- .mrpkg_adapter("cml", "MR-cML", 3L, c(
  .enum_p("cml_dp", "dp", c("dp", "ma"), c(dp = "cML-MA-BIC-DP: with data perturbation (recommended)", ma = "cML-MA-BIC: without data perturbation"), "Variant",
          "Data perturbation allows for uncertainty in which variants are selected as invalid, giving wider and better-calibrated intervals.", "DP", "Xue et al. 2021; MendelianRandomization::mr_cML default (DP = TRUE)."),
  .num_p("cml_pert", "integer", 200L, 50L, 2000L, "Perturbations", "Number of data perturbations (cML-DP only).", "num_pert", "MendelianRandomization::mr_cML default (num_pert = 200).", c("se", "ci", "p"))), TRUE, c(
  "# the sample size is the smaller of the exposure and outcome GWAS medians, as the package recommends",
  "n_cml <- min(stats::median(dat$n_x), stats::median(dat$n_y))",
  "fit <- MendelianRandomization::mr_cML(mri, MA = TRUE, DP = identical(cml_dp, \"dp\"), num_pert = cml_pert, random_seed = seed, n = n_cml, Alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
  "inv <- if (length(fit@BIC_invalid)) dat$variant_id[fit@BIC_invalid] else character()",
  "diagnostics <- rbind(drow(\"cml_invalid\", length(inv), note = if (length(inv)) paste(inv, collapse = \";\") else NA_character_),",
  "  drow(\"cml_gof1\", NA, p = fit@GOF1_p, note = \"goodness of fit, test 1\"), drow(\"cml_gof2\", NA, p = fit@GOF2_p, note = \"goodness of fit, test 2\"), drow(\"cml_n\", n_cml, note = \"sample size used\"))"),
  needs_n = TRUE, check_extra = .cml_n_check)

# MR-Robust (Burgess et al. 2016; Rees et al. 2019): IVW with robust regression (MM-estimation) and/or penalised weights.
ADAPTER_ROBUST_V11 <- .mrpkg_adapter("robust_ivw", "MR-Robust", 3L,
  .enum_p("rob_variant", "both", c("both", "robust", "penalized"), c(both = "Robust regression with penalised weights (MR-Robust)", robust = "Robust regression (MM-estimation) only", penalized = "Penalised weights only"),
          "Variant", "Robust regression down-weights outlying variants; penalised weights shrink the weight of variants whose ratio estimate is heterogeneous.", "robust, penalized",
          "Rees et al. 2019 (PLoS One): the robust and penalised IVW performed best together."), FALSE, c(
  "fit <- MendelianRandomization::mr_ivw(mri, robust = rob_variant %in% c(\"both\", \"robust\"), penalized = rob_variant %in% c(\"both\", \"penalized\"), alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
  "diagnostics <- drow(\"robust_variant\", NA, note = rob_variant)"))

# Penalised IVW (Xu, Wang, Fung and Liu 2023): penalised and debiased IVW with a bootstrap Fieller interval.
ADAPTER_PIVW_V11 <- .mrpkg_adapter("pivw", "Penalised IVW", 3L, c(
  .num_p("piv_lambda", "double", 1, 0, 10, "Penalty (lambda)", "Penalty on weak variants; 1 is the package default.", "lambda", "MendelianRandomization::mr_pivw default (lambda = 1)."),
  .enum_p("piv_disp", "balanced", c("balanced", "none"), c(balanced = "Allow for balanced pleiotropy (standard)", none = "No pleiotropy"), "Pleiotropy in the standard error",
          "Whether the standard error allows for balanced horizontal pleiotropy.", "over.dispersion", "MendelianRandomization::mr_pivw default (over.dispersion = TRUE).", c("se", "ci", "p")),
  .enum_p("piv_ci", "fieller", c("fieller", "wald"), c(fieller = "Bootstrap Fieller interval (standard)", wald = "Normal (Wald) interval"), "Interval",
          "The bootstrap Fieller interval is valid with weak instruments.", "Boot.Fieller", "MendelianRandomization::mr_pivw default (Boot.Fieller = TRUE).", c("ci"))), TRUE, c(
  "fit <- MendelianRandomization::mr_pivw(mri, lambda = piv_lambda, over.dispersion = identical(piv_disp, \"balanced\"), Boot.Fieller = identical(piv_ci, \"fieller\"), alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError, min(fit@CILower), max(fit@CIUpper), fit@Pvalue, fit@SNPs, if (identical(piv_ci, \"fieller\")) \"bootstrap\" else \"wald\")",
  "diagnostics <- rbind(drow(\"pivw_condition\", fit@Condition, note = \"(mean F - 1) x sqrt(number of variants)\"), drow(\"pivw_tau2\", fit@Tau2, note = \"balanced pleiotropy variance\"))"))

# ---- Round 59: TwoSampleMR estimators and diagnostics ---------------------------------------------------------------------
PIN_TWOSAMPLEMR_V11 <- "5e93ee7d8e1f67581acb8134299759389d73efba"
.tsmr_params <- c("pars <- TwoSampleMR::default_parameters()   # the package defaults, with this method's settings below",
  "z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  "tsd <- data.frame(id.exposure = EXPOSURE_ID, id.outcome = OUTCOME_ID, exposure = EXPOSURE_ID, outcome = OUTCOME_ID, SNP = dat$variant_id, beta.exposure = dat$beta_x,",
  "  beta.outcome = dat$beta_y, se.exposure = dat$se_x, se.outcome = dat$se_y, mr_keep = TRUE, stringsAsFactors = FALSE)")
.tsmr_wald <- "est_w <- function(f) est1(f$b, f$se, f$b - z * f$se, f$b + z * f$se, f$pval, f$nsnp)"
.nboot_p <- function(name, def = 1000L) .num_p(name, "integer", def, 100L, 100000L, "Bootstrap iterations", "Iterations for the bootstrap standard error.", "nboot", "TwoSampleMR default (nboot = 1000).", c("se", "ci", "p"))

# MR-Egger with bootstrap standard errors (Bowden et al. 2015), TwoSampleMR::mr_egger_regression_bootstrap
ADAPTER_EGGER_BOOT_V11 <- .mrpkg_adapter("egger_boot", "MR-Egger (bootstrap)", 3L, .nboot_p("eb_nboot"), TRUE, c(.tsmr_params, .tsmr_wald,
  "pars$nboot <- eb_nboot",
  "f <- TwoSampleMR::mr_egger_regression_bootstrap(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "estimates <- rbind(est_w(f), transform(est1(f$b_i, f$se_i, f$b_i - z * f$se_i, f$b_i + z * f$se_i, f$pval_i, f$nsnp), estimand = \"pleiotropy_intercept\"))",
  "diagnostics <- drow(\"egger_boot_nboot\", eb_nboot, note = \"bootstrap iterations\")"), pkg = "TwoSampleMR", estimands = c("causal_effect", "pleiotropy_intercept"))

# Penalised weighted median (Burgess et al. 2016), TwoSampleMR::mr_penalised_weighted_median
ADAPTER_PWM_V11 <- .mrpkg_adapter("pwm", "Penalised weighted median", 3L, c(.nboot_p("pwm_nboot"),
  .num_p("pwm_penk", "double", 20, 1, 1000, "Penalty constant (penk)", "How strongly variants with heterogeneous ratio estimates are down-weighted.", "penk", "TwoSampleMR default (penk = 20).")), TRUE, c(.tsmr_params, .tsmr_wald,
  "pars$nboot <- pwm_nboot; pars$penk <- pwm_penk",
  "f <- TwoSampleMR::mr_penalised_weighted_median(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "estimates <- est_w(f)", "diagnostics <- drow(\"pwm_penk\", pwm_penk, note = \"penalty constant\")"), pkg = "TwoSampleMR")

# Simple and weighted mode under the NOME assumption (Hartwig et al. 2017), TwoSampleMR::mr_simple_mode_nome / mr_weighted_mode_nome
.nome_params <- function(pref) c(.nboot_p(paste0(pref, "_nboot")),
  .num_p(paste0(pref, "_phi"), "double", 1, 0.1, 5, "Bandwidth multiplier (phi)", "Smaller values give a narrower kernel.", "phi", "TwoSampleMR default (phi = 1)."))
ADAPTER_SMODE_NOME_V11 <- .mrpkg_adapter("smode_nome", "Simple mode (NOME)", 3L, .nome_params("smn"), TRUE, c(.tsmr_params, .tsmr_wald,
  "pars$nboot <- smn_nboot; pars$phi <- smn_phi",
  "f <- TwoSampleMR::mr_simple_mode_nome(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "estimates <- est_w(f)", "diagnostics <- drow(\"mode_nome\", NA, note = \"ratio standard errors assume no measurement error in the exposure associations (NOME)\")"), pkg = "TwoSampleMR")
ADAPTER_WMODE_NOME_V11 <- .mrpkg_adapter("wmode_nome", "Weighted mode (NOME)", 3L, .nome_params("wmn"), TRUE, c(.tsmr_params, .tsmr_wald,
  "pars$nboot <- wmn_nboot; pars$phi <- wmn_phi",
  "f <- TwoSampleMR::mr_weighted_mode_nome(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "estimates <- est_w(f)", "diagnostics <- drow(\"mode_nome\", NA, note = \"ratio standard errors assume no measurement error in the exposure associations (NOME)\")"), pkg = "TwoSampleMR")

# Rucker's framework (Bowden et al. 2018): Q (IVW) against Rucker's Q' (MR-Egger), model selection between IVW and MR-Egger, and its
# bootstrap and jackknife versions. TwoSampleMR::mr_rucker, mr_rucker_bootstrap, mr_rucker_jackknife.
ADAPTER_RUCKER_V11 <- .mrpkg_adapter("rucker", "Rucker model selection", 4L, c(
  .enum_p("rk_variant", "point", c("point", "bootstrap", "jackknife"), c(point = "Selected model (point estimate)", bootstrap = "Bootstrap over variants", jackknife = "Jackknife"),
          "Version", "The point estimate picks IVW or MR-Egger by comparing Q and Rucker's Q'; the bootstrap and jackknife repeat the selection on resampled variants.", "mr_rucker / _bootstrap / _jackknife",
          "Bowden et al. 2018 (Int J Epidemiol); TwoSampleMR."),
  .enum_p("rk_stat", "median", c("median", "mean"), c(median = "Median of the resampled estimates", mean = "Mean of the resampled estimates"), "Summary (bootstrap or jackknife)",
          "How the resampled estimates are summarised.", "Rucker median / mean", "TwoSampleMR reports both."),
  .nboot_p("rk_nboot")), TRUE, c(.tsmr_params,
  "pars$nboot <- rk_nboot",
  "r1 <- suppressMessages(TwoSampleMR::mr_rucker(tsd, pars))[[1]]",
  "res <- switch(rk_variant, point = r1$selected, bootstrap = suppressMessages(TwoSampleMR::mr_rucker_bootstrap(tsd, pars))$res,",
  "  jackknife = suppressMessages(TwoSampleMR::mr_rucker_jackknife(tsd, pars))[[1]]$res)",
  "row <- if (rk_variant == \"point\") res[1, ] else res[grepl(sprintf(\"Rucker %s\", rk_stat), res$Method), ][1, ]",
  "# the interval limits are put in order: the jackknife mean row of TwoSampleMR lists them the other way round",
  "estimates <- est1(row$Estimate, row$SE, min(row$CI_low, row$CI_upp), max(row$CI_low, row$CI_upp), row$P, row$nsnp, if (rk_variant == \"point\") \"wald\" else \"bootstrap\")",
  "q <- r1$Q",
  "diagnostics <- rbind(drow(\"q_ivw\", q$Q[1], q$df[1], q$P[1], \"Cochran's Q (IVW)\"), drow(\"q_rucker\", q$Q[2], q$df[2], q$P[2], \"Rucker's Q' (MR-Egger)\"),",
  "  drow(\"q_diff\", q$Q[3], q$df[3], q$P[3], \"Q - Q'\"), drow(\"rucker_model\", NA, note = r1$res))"), pkg = "TwoSampleMR")

# Sign concordance test (TwoSampleMR::mr_sign): are the variant-exposure and variant-outcome associations in the same direction more often
# than chance? A diagnostic, not an effect estimate.
ADAPTER_SIGN_V11 <- .mrpkg_adapter("sign", "Sign concordance test", 6L, list(), FALSE, c(.tsmr_params,
  "f <- TwoSampleMR::mr_sign(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "estimates <- est1(0, 0, 0, 0, 0, 0L)[0, ]   # no effect estimate: an empty table with the standard columns",
  "diagnostics <- drow(\"sign_concordance\", f$b, f$nsnp, f$pval, \"(proportion of variants with the same sign - 0.5) x 2; binomial test\")"), pkg = "TwoSampleMR", estimands = character(0))

# ---- Round 60: MRMix (Qi and Chatterjee 2019, Nat Commun 10:1941) --------------------------------------------------------------
# MRMix 0.1.0 from GitHub (gqi/MRMix), pinned to the commit the runtime installs. It works on standardised associations (per SD of the
# genotype; continuous traits also per SD of the trait, as MRMix::standardize). MR Studio reports the estimate per unit of the exposure
# like the other methods, by multiplying by the median implied SD of the outcome (continuous outcome only) and dividing by the median
# implied SD of the exposure, SD = se x sqrt(n x 2p(1 - p)); the standardised estimate is kept as a diagnostic.
PIN_MRMIX_V11 <- "56afdb2bc96760842405396f5d3f02e60e305039"
.mrmix_check <- function(d, params) {
  if (!all(c("n_x", "n_y") %in% names(d)) || !all(is.finite(d$n_x) & is.finite(d$n_y))) return(verdict("block", "FEAS_NEEDS_N", "MRMix needs the sample sizes of the exposure and outcome GWAS (n_x, n_y)."))
  f <- if ("eaf_x" %in% names(d)) d$eaf_x else rep(NA_real_, nrow(d)); if ("eaf_y" %in% names(d)) f[!is.finite(f)] <- d$eaf_y[!is.finite(f)]
  if (sum(is.finite(f) & f > 0 & f < 1) < 3) return(verdict("block", "FEAS_NEEDS_EAF", "MRMix needs effect allele frequencies (eaf_x, or eaf_y in its place) to standardise the associations."))
  NULL }
ADAPTER_MRMIX_V11 <- .mrpkg_adapter("mrmix", "MRMix", 3L, list(), FALSE, c(
  "z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  "maf <- if (\"eaf_x\" %in% names(dat)) dat$eaf_x else rep(NA_real_, nrow(dat)); if (\"eaf_y\" %in% names(dat)) maf[!is.finite(maf)] <- dat$eaf_y[!is.finite(maf)]   # eaf_y in place of a missing eaf_x, as TwoSampleMR::run_mrmix",
  "dat <- dat[is.finite(maf) & maf > 0 & maf < 1, , drop = FALSE]; maf <- maf[is.finite(maf) & maf > 0 & maf < 1]",
  "xt <- if (XSCALE == \"log_or\") \"binary\" else \"continuous\"; yt <- if (SCALE == \"log_or\") \"binary\" else \"continuous\"",
  "l <- MRMix::standardize(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, xt, yt, dat$n_x, dat$n_y, maf)",
  "f <- MRMix::MRMix(l$betahat_x_std, l$betahat_y_std, l$sx_std, l$sy_std)",
  "g <- sqrt(2 * maf * (1 - maf)); sd_x <- if (xt == \"continuous\") stats::median(dat$se_x * sqrt(dat$n_x) * g) else 1; sd_y <- if (yt == \"continuous\") stats::median(dat$se_y * sqrt(dat$n_y) * g) else 1",
  "k <- sd_y / sd_x; b <- f$theta * k; se <- f$SE_theta * k   # back to the scale of the other methods (per unit of the exposure)",
  "estimates <- est1(b, se, b - z * se, b + z * se, f$pvalue_theta, nrow(dat))",
  "diagnostics <- rbind(drow(\"mrmix_theta_std\", f$theta, p = f$pvalue_theta, note = sprintf(\"standardised estimate (SE %.4g), per SD of the exposure\", f$SE_theta)),",
  "  drow(\"mrmix_pi0\", f$pi0, note = \"estimated proportion of variants with no pleiotropic effect\"), drow(\"mrmix_sigma2\", f$sigma2, note = \"variance of the pleiotropic effects (standardised)\"),",
  "  drow(\"mrmix_scale\", k, note = sprintf(\"implied SD of exposure %.4g, of outcome %.4g\", sd_x, sd_y)))"), needs_n = TRUE, check_extra = .mrmix_check, pkg = "MRMix")

# ---- Round 61: MR-GRIP (Dudbridge et al. 2025, PLoS Genet 21:e1011967), TwoSampleMR::mr_grip ------------------------------------
# MR-Egger made independent of allele coding: the regression of beta_y x beta_x on beta_x^2 with an intercept, weighted by
# 1 / (beta_x^2 se_y^2); with the intercept fixed at zero it equals IVW. The weak-instrument adjusted version is also offered.
ADAPTER_GRIP_V11 <- .mrpkg_adapter("grip", "MR-GRIP", 3L,
  .enum_p("grip_variant", "standard", c("standard", "weak"), c(standard = "MR-GRIP", weak = "MR-GRIP adjusted for weak instruments"), "Version",
          "The adjusted version corrects for weak instruments; it is less precise.", "b / b.adj", "Dudbridge et al. 2025 report both (Table 6)."), FALSE, c(.tsmr_params,
  "f <- TwoSampleMR::mr_grip(dat$beta_x, dat$beta_y, dat$se_x, dat$se_y, pars)",
  "b <- if (grip_variant == \"weak\") f$b.adj else f$b; se <- if (grip_variant == \"weak\") f$se.adj else f$se; p <- if (grip_variant == \"weak\") f$pval.adj else f$pval",
  "estimates <- rbind(est1(b, se, b - z * se, b + z * se, p, f$nsnp), transform(est1(f$b_i, f$se_i, f$b_i - z * f$se_i, f$b_i + z * f$se_i, f$pval_i, f$nsnp), estimand = \"pleiotropy_intercept\"))",
  "diagnostics <- drow(\"grip_variant\", NA, note = grip_variant)"), pkg = "TwoSampleMR", estimands = c("causal_effect", "pleiotropy_intercept"))

# ---- Round 62, rewritten in Round 70e: BWMR (Zhao et al. 2020, Bioinformatics 36:1501-1508, PMID 31593215) ------------------------
# Bayesian weighted MR: variational EM with weights that down-weight outlying variants and a linear-response correction of the
# posterior variance. The authors' package (GitHub jiazhao97/BWMR) carries no licence, so it cannot be distributed with MR Studio;
# this is MR Studio's own implementation of the published algorithm, written from the equations of the paper and its supplement
# (arXiv 1811.10223v2), and validated against the package (validation/zhao2020_bwmr). It needs no package beyond base R, and every
# study script contains the whole function. It iterates to convergence of the evidence lower bound (relative change below 1e-12),
# where the package stops earlier (about 1e-6); see the validation README for the size of that difference.
BWMR_MRS_SRC <- strsplit(r"---(# BWMR (Zhao et al. 2020, Bioinformatics 36:1501-1508), written for MR Studio from the equations of the paper and its supplement
# (arXiv 1811.10223v2: Methods 2.1-2.3; Supplement S1.1 E-step and M-step, S1.2 linear response variance). MIT licence (MR Studio).
# Model (bx_j, by_j the observed variant-exposure and variant-outcome effects): bx_j ~ N(gamma_j, se_x_j^2); by_j ~ N(beta gamma_j, se_y_j^2 + tau^2)^w_j; gamma_j ~ N(0, sigma^2); beta ~ N(0, sigma0^2);
# w_j ~ Bernoulli(pi1), pi1 ~ Beta(a0, 1). Variational EM (mean field) for beta, gamma, w, pi1; M-step for sigma^2 and tau^2; the posterior
# variance of beta corrected by linear response variational Bayes, computed as -[(V H1 - I)^-1 V]_11 (supplement, reparameterised form).
# Starting values and the stopping rule are not given in the paper; they only affect how fast the iterations reach the same solution.
bwmr_mrs <- function(bx, by, se_x, se_y, a0 = 100, sigma0 = 1e6, tol = 1e-12, max_iter = 100000L, init = NULL) {
  N <- length(bx); sx2 <- se_x^2; sy2 <- se_y^2; G <- by; g <- bx   # g, G: variant-exposure and variant-outcome effects
  xlx <- function(p) ifelse(p > 0, p * log(p), 0)
  # starting values: IVW for beta, the observed effects for gamma, every variant weighted 1
  mu_b <- sum(g * G / sy2) / sum(g^2 / sy2); s2_b <- 1 / sum(g^2 / sy2)
  mu_g <- g; s2_g <- sx2; pw <- rep(1, N); sigma2 <- mean(g^2)
  tau2 <- max(mean((G - mu_b * g)^2) - mean(sy2), mean(sy2) / 10)
  if (!is.null(init)) { mu_b <- init$estimate; pw <- init$weights; tau2 <- init$tau2 }
  elbo <- -Inf; it <- 0L
  repeat {
    it <- it + 1L; v <- sy2 + tau2
    # E-step (supplement S1.1), in the order of the paper: beta, gamma_j, pi1, w_j
    s2_b <- 1 / (1 / sigma0^2 + sum(pw * (mu_g^2 + s2_g) / v))
    mu_b <- sum(pw * mu_g * G / v) * s2_b
    s2_g <- 1 / (1 / sx2 + pw * (mu_b^2 + s2_b) / v + 1 / sigma2)
    mu_g <- (g / sx2 + pw * mu_b * G / v) * s2_g
    a <- a0 + sum(pw); b <- N + 1 - sum(pw)
    q <- (mu_b^2 + s2_b) * (mu_g^2 + s2_g) - 2 * mu_b * mu_g * G + G^2   # E[(by_j - beta gamma_j)^2]
    l1 <- -0.5 * log(2 * pi) - 0.5 * log(v) - q / (2 * v) + digamma(a) - digamma(a + b)
    l0 <- digamma(b) - digamma(a + b)
    pw <- 1 / (1 + exp(l0 - l1))                                          # q_j1 / (q_j0 + q_j1), on the log scale
    # M-step (supplement S1.1): sigma^2 in closed form, tau^2 by the bound of the paper
    sigma2 <- mean(mu_g^2 + s2_g)
    tau2 <- sqrt(sum(pw * q * tau2^2 / v^2) / sum(pw / v))
    # evidence lower bound, to decide convergence
    v <- sy2 + tau2; da <- digamma(a) - digamma(a + b); db <- digamma(b) - digamma(a + b)
    e_new <- sum(-0.5 * log(sx2) - ((g - mu_g)^2 + s2_g) / (2 * sx2)) + sum(pw * (-0.5 * log(2 * pi) - 0.5 * log(v) - q / (2 * v))) -
      (mu_b^2 + s2_b) / (2 * sigma0^2) - N / 2 * log(sigma2) - sum(mu_g^2 + s2_g) / (2 * sigma2) + sum(pw * da + (1 - pw) * db) + (a0 - 1) * da +
      0.5 * log(s2_b) + 0.5 * sum(log(s2_g)) - sum(xlx(pw) + xlx(1 - pw)) - ((a - 1) * da + (b - 1) * db - lbeta(a, b))
    done <- is.finite(elbo) && abs(e_new - elbo) <= tol * abs(e_new); elbo <- e_new
    if (done || it >= max_iter) break
  }
  # linear response variance of beta (supplement S1.2): m = (E beta, E beta^2, [E gamma_j, E gamma_j^2, E w_j]_j, E log pi1, E log(1 - pi1)).
  # Var(beta) = -[(V H1 - I)^-1 V]_11. V is block diagonal and H1 links each variant's block only to the global blocks (beta, pi1), so
  # A = V H1 - I has an arrow shape and is solved exactly through the 4 x 4 Schur complement of the global block (time linear in N).
  v <- sy2 + tau2; a <- a0 + sum(pw); b <- N + 1 - sum(pw)          # q(pi1) consistent with the final weights
  VG <- matrix(0, 4, 4); VG[1:2, 1:2] <- c(s2_b, 2 * mu_b * s2_b, 2 * mu_b * s2_b, 2 * s2_b^2 + 4 * mu_b^2 * s2_b)
  VG[3:4, 3:4] <- c(trigamma(a) - trigamma(a + b), -trigamma(a + b), -trigamma(a + b), trigamma(b) - trigamma(a + b))
  S <- -diag(4)                                                       # A_GG = VG H1_GG - I, and H1_GG = 0
  ok <- TRUE
  for (j in seq_len(N)) {
    Vj <- matrix(c(s2_g[j], 2 * mu_g[j] * s2_g[j], 0, 2 * mu_g[j] * s2_g[j], 2 * s2_g[j]^2 + 4 * mu_g[j]^2 * s2_g[j], 0, 0, 0, pw[j] * (1 - pw[j])), 3, 3)
    Hjj <- matrix(0, 3, 3); Hjj[1, 3] <- Hjj[3, 1] <- mu_b * G[j] / v[j]; Hjj[2, 3] <- Hjj[3, 2] <- -(mu_b^2 + s2_b) / (2 * v[j])
    HjG <- matrix(0, 3, 4)                                            # rows gamma_j, gamma_j^2, w_j; columns beta, beta^2, log pi1, log(1 - pi1)
    HjG[1, 1] <- pw[j] * G[j] / v[j]; HjG[2, 2] <- -pw[j] / (2 * v[j])
    HjG[3, 1] <- mu_g[j] * G[j] / v[j]; HjG[3, 2] <- -(mu_g[j]^2 + s2_g[j]) / (2 * v[j]); HjG[3, 3] <- 1; HjG[3, 4] <- -1
    Ajj <- Vj %*% Hjj - diag(3); AjG <- Vj %*% HjG; AGj <- VG %*% t(HjG)
    Z <- tryCatch(solve(Ajj, AjG), error = function(e) NULL); if (is.null(Z)) { ok <- FALSE; break }
    S <- S - AGj %*% Z }
  x <- if (ok) tryCatch(solve(S, c(VG[1, 1], VG[2, 1], 0, 0)), error = function(e) NULL) else NULL
  var_b <- if (is.null(x)) NA_real_ else -x[1]
  se <- sqrt(var_b)
  list(estimate = mu_b, se = se, p = 2 * stats::pnorm(-abs(mu_b / se)), weights = pw, tau2 = tau2, sigma2 = sigma2,
       se_mfvb = sqrt(s2_b), mu_gamma = mu_g, s2_gamma = s2_g, a = a, b = b, iterations = it, converged = it < max_iter, elbo = elbo)
})---", "\n")[[1]]
ADAPTER_BWMR_V11 <- .mrpkg_adapter("bwmr", "BWMR", 3L, list(), FALSE, c(
  "z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  BWMR_MRS_SRC,
  "f <- bwmr_mrs(bx = dat$beta_x, by = dat$beta_y, se_x = dat$se_x, se_y = dat$se_y)",
  "if (!f$converged || !is.finite(f$se)) stop(\"BWMR did not converge on these data (the variance correction could not be computed).\")",
  "estimates <- est1(f$estimate, f$se, f$estimate - z * f$se, f$estimate + z * f$se, f$p, nrow(dat))",
  "diagnostics <- rbind(drow(\"bwmr_downweighted\", sum(f$weights < 0.5), note = \"variants with posterior weight below 0.5 (likely outliers)\"),",
  "  drow(\"bwmr_tau2\", f$tau2, note = \"variance of the polygenic (pleiotropic) effects\"),",
  "  drow(\"bwmr_iterations\", f$iterations, note = \"variational EM iterations to convergence of the evidence lower bound\"))"), pkg = "stats")

# ---- Round 68b: IVW after Steiger filtering (Hemani, Tilling and Davey Smith 2017, PLoS Genet 13:e1007081) ------------------
# A sensitivity analysis, never the primary result: variants that explain more variance in the outcome than in the exposure are left
# out and the IVW estimate (TwoSampleMR::mr_ivw, as the primary IVW) is computed on the rest. Variance explained exactly as in the
# Steiger directionality test (continuous: from effect, SE and N; binary: liability scale from the log odds ratio, allele frequency,
# cases, controls and prevalence). Filtering uses the outcome data to choose instruments, which can itself bias the estimate when
# pleiotropy or sample overlap is present, so it is offered off by default and labelled as filtered (Burgess et al. 2023 guidelines).
ADAPTER_IVW_STEIGER_V11 <- .mrpkg_adapter("ivw_steiger", "IVW after Steiger filtering", 3L, list(), FALSE, c(
  .steiger_script_side("x", "exposure"), .steiger_script_side("y", "outcome"),
  "keep <- abs(r_x) > abs(r_y)   # Steiger filtering: keep the variants that explain more variance in the exposure than in the outcome",
  "kept <- dat[keep, , drop = FALSE]",
  "if (nrow(kept) < 3) stop(sprintf(\"Only %d variants remain after Steiger filtering; IVW needs at least 3.\", nrow(kept)))",
  "fit <- TwoSampleMR::mr_ivw(kept$beta_x, kept$beta_y, kept$se_x, kept$se_y)   # multiplicative random effects, SE not below the fixed-effect SE",
  "z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  "estimates <- est1(fit$b, fit$se, fit$b - z * fit$se, fit$b + z * fit$se, fit$pval, fit$nsnp)",
  "diagnostics <- rbind(drow(\"steiger_kept\", sum(keep), note = \"variants kept (more variance explained in the exposure)\"),",
  "  drow(\"steiger_removed\", sum(!keep), note = if (any(!keep)) paste(dat$variant_id[!keep], collapse = \" \") else \"none\"))"),
  needs_n = TRUE, pkg = "TwoSampleMR",
  check_extra = function(d, params, input) {
    for (s in list(list("x", input$exposure_scale %||% "beta", "exposure"), list("y", input$outcome_scale, "outcome"))) {
      pr <- .steiger_side_problem(d, s[[1]], s[[2]], s[[3]]); if (!is.null(pr)) return(verdict("block", pr$code, pr$msg)) }
    NULL })


# ---- Round 109: further MendelianRandomization 0.10.0 estimators (help pages of mr_egger, mr_hetpen, mr_clr, mr_mvivwme, mr_mvgmm) ----
# MR-Egger with robust regression and/or penalised weights (Burgess and Thompson 2017; Rees et al. 2019).
ADAPTER_EGGER_ROB_V11 <- .mrpkg_adapter("egger_rob", "MR-Egger (robust / penalised)", 3L,
  .enum_p("egr_variant", "both", c("both", "robust", "penalized"), c(both = "Robust regression with penalised weights", robust = "Robust regression (MM-estimation) only", penalized = "Penalised weights only"),
          "Variant", "Robust regression down-weights outlying variants; penalised weights shrink the weight of variants whose ratio estimate is heterogeneous.", "robust, penalized",
          "MendelianRandomization::mr_egger (robust, penalized); Rees et al. 2019 (PLoS One)."), FALSE, c(
  "fit <- MendelianRandomization::mr_egger(mri, robust = egr_variant %in% c(\"both\", \"robust\"), penalized = egr_variant %in% c(\"both\", \"penalized\"), alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, fit@StdError.Est, fit@CILower.Est, fit@CIUpper.Est, fit@Pvalue.Est, fit@SNPs)",
  "diagnostics <- rbind(drow(\"egger_rob_intercept\", fit@Intercept, p = fit@Pvalue.Int, note = sprintf(\"intercept %.4g (SE %.4g, CI %.4g to %.4g)\", fit@Intercept, fit@StdError.Int, fit@CILower.Int, fit@CIUpper.Int)),",
  "  drow(\"egger_rob_variant\", NA, note = egr_variant))"))

# Heterogeneity-penalised model averaging (Burgess et al. 2018). The method weighs every subset of variants, so its time grows
# exponentially with the number of variants: MR Studio allows at most HETPEN_MAX_IVS. Its interval can be several separate ranges.
HETPEN_MAX_IVS <- 16L   # about 10 s at 16 variants on an Intel Mac; each 2 more variants take about four times as long
.hetpen_check <- function(d, params) if (nrow(d) > HETPEN_MAX_IVS) verdict("block", "FEAS_TOO_MANY_IVS", sprintf("Heterogeneity-penalised model averaging weighs every subset of variants, so it takes hours with many variants; MR Studio allows at most %d (this pair has %d).", HETPEN_MAX_IVS, nrow(d))) else NULL
ADAPTER_HETPEN_V11 <- .mrpkg_adapter("hetpen", "Heterogeneity-penalised", 3L, c(
  .num_p("hp_prior", "double", 0.5, 0.01, 0.99, "Prior probability of a valid variant", "Prior probability that each variant is a valid instrument.", "prior", "MendelianRandomization::mr_hetpen default (prior = 0.5).", advanced = FALSE),
  .num_p("hp_grid", "integer", 2000L, 200L, 20000L, "Grid points for the interval", "The interval is searched on this many points between the smallest and largest plausible values.", "CIStep", "MR Studio: the range of the ratio estimates widened by 1 on each side, split into this many steps.", c("ci"))), FALSE, c(
  "r <- dat$beta_y / dat$beta_x; lo <- min(r) - 1; hi <- max(r) + 1   # search range for the interval: the ratio estimates plus 1 on each side",
  "fit <- MendelianRandomization::mr_hetpen(mri, prior = hp_prior, CIMin = lo, CIMax = hi, CIStep = (hi - lo) / hp_grid, alpha = 1 - ci_level)",
  "estimates <- est1(fit@Estimate, NA, min(fit@CILower), max(fit@CIUpper), NA, fit@SNPs, \"profile\")",
  "ivs <- paste(sprintf(\"%.4g to %.4g\", fit@CILower, fit@CIUpper), collapse = \"; \")",
  "diagnostics <- drow(\"hetpen_intervals\", length(fit@CILower), note = if (length(fit@CILower) > 1) paste(\"the interval is several ranges:\", ivs) else paste(\"one range:\", ivs))"),
  check_extra = .hetpen_check)

# Conditional likelihood ratio method (Wang and Kang 2022, as in MendelianRandomization): weak-instrument robust test and intervals
# (Anderson-Rubin, Kleibergen, conditional likelihood ratio). It gives intervals, not a point estimate, so its results are diagnostics.
.clr_n_check <- function(d, params) if (!all(c("n_x", "n_y") %in% names(d)) || !all(is.finite(d$n_x) & is.finite(d$n_y))) verdict("block", "FEAS_NEEDS_N", "The conditional likelihood ratio method needs the sample sizes of the exposure and outcome GWAS (n_x, n_y).") else NULL
ADAPTER_CLR_V11 <- .mrpkg_adapter("clr", "Conditional likelihood ratio", 2L,
  .num_p("clr_grid", "integer", 2000L, 200L, 20000L, "Grid points for the intervals", "The intervals are searched on this many points.", "CIStep", "MR Studio: the range of the ratio estimates widened by 1 on each side.", c("ci")), FALSE, c(
  "r <- dat$beta_y / dat$beta_x; lo <- min(r) - 1; hi <- max(r) + 1",
  "nx <- stats::median(dat$n_x); ny <- stats::median(dat$n_y)   # median sample sizes of the two GWAS",
  "fit <- MendelianRandomization::mr_clr(mri, nx = nx, ny = ny, alpha = 1 - ci_level, CIMin = lo, CIMax = hi, CIStep = (hi - lo) / clr_grid)",
  "iv <- function(l, u) if (!length(l) || all(is.na(l))) \"no value is compatible with the data\" else paste(sprintf(\"%.4g to %.4g\", l, u), collapse = \"; \")",
  "estimates <- est1(NA, NA, NA, NA, NA, NA)[0, ]   # no point estimate",
  "diagnostics <- rbind(drow(\"clr_interval\", length(fit@CLRlower), note = iv(fit@CLRlower, fit@CLRupper)), drow(\"clr_ar_interval\", length(fit@ARlower), note = iv(fit@ARlower, fit@ARupper)),",
  "  drow(\"clr_k_interval\", length(fit@Klower), note = iv(fit@Klower, fit@Kupper)), drow(\"clr_n\", nx, note = sprintf(\"exposure N %g, outcome N %g\", nx, ny)))"),
  needs_n = TRUE, check_extra = .clr_n_check)

# MVMR-IVW allowing for measurement error in the variant-exposure associations (Yao et al. 2023) and multivariable GMM, robust to
# weak instruments and sample overlap (Patel et al. 2024). The exposure correlation matrix is the phenotypic correlations (mv_r_i_j).
.mv_rx <- "Rx <- diag(K); if (K > 1) for (i in seq_len(K - 1)) for (j in (i + 1):K) Rx[i, j] <- Rx[j, i] <- get(sprintf(\"mv_r_%d_%d\", i, j))   # exposure correlations"
.mv_body$mvivwme <- c("# MVMR-IVW with measurement error in the exposure associations (MendelianRandomization::mr_mvivwme)", .mv_rx,
  "fit <- MendelianRandomization::mr_mvivwme(inp, model = mv_model, correl.x = Rx, alpha = 1 - ci_level, seed = seed)",
  "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, fit@SNPs)",
  "diagnostics <- drow(\"mv_ivwme_heterogeneity\", fit@Heter.Stat[1], fit@SNPs - K, fit@Heter.Stat[2], \"heterogeneity statistic\")")
.mv_body$mvgmm <- c("# Multivariable GMM (MendelianRandomization::mr_mvgmm), robust to weak instruments and overlapping samples", .mv_rx,
  "nx <- vapply(seq_len(K), function(k) stats::median(dat[[sprintf(\"n_x_%d\", k)]]), 0); ny <- stats::median(dat$n_y)",
  "fit <- MendelianRandomization::mr_mvgmm(inp, nx = nx, ny = ny, cor.x = Rx, robust = identical(mv_gmm_robust, \"robust\"), alpha = 1 - ci_level)",
  "estimates <- est_rows(fit@Estimate, fit@StdError, fit@CILower, fit@CIUpper, fit@Pvalue, nrow(dat))",
  "diagnostics <- rbind(drow(\"mv_gmm_heterogeneity\", fit@Heter.Stat[1], nrow(dat) - K, fit@Heter.Stat[2], \"heterogeneity statistic\"),",
  "  do.call(rbind, lapply(seq_len(K), function(k) drow(sprintf(\"mv_gmm_cond_f_%d\", k), fit@CondFstat[k], note = sprintf(\"conditional F, exposure %d\", k)))))")
.mv_gmm_params <- list(mv_gmm_robust = list(type = "enum", default = "robust", values = c("robust", "plain"), value_labels = c(robust = "Allow for balanced pleiotropy (robust, package default)", plain = "No pleiotropy"),
  label = "Pleiotropy", hint = "Whether the estimator allows for balanced horizontal pleiotropy.", package_arg = "robust", kind = "enum", conventional = TRUE, advanced = FALSE,
  evidence = "MendelianRandomization::mr_mvgmm default (robust = TRUE).", affects = c("estimate", "se", "ci", "p")))
ADAPTER_MVIVWME_V11 <- .mv_adapter("mvivwme", "MVMR-IVW (measurement error)", 1L, c(.mv_model_param, .mv_cor_params), TRUE, "mvivwme")
ADAPTER_MVGMM_V11 <- .mv_adapter("mvgmm", "MVMR GMM", 2L, c(.mv_gmm_params, .mv_cor_params), FALSE, "mvgmm")

# Round 109: influence of each variant on the IVW fit (weighted regression of beta_y on beta_x through the origin, weights 1 / se_y^2):
# Cook's distance (flagged above 4 / k) and the externally studentized residual (flagged when |r| is above the cut-off, 2 by default).
# They give no causal estimate; the flagged variants feed the outlier analysis on the Instrument review screen.
.infl_lines <- c("fit <- stats::lm(beta_y ~ beta_x - 1, data = dat, weights = 1 / dat$se_y^2)   # the IVW (fixed-effect) regression",
  "estimates <- est1(NA, NA, NA, NA, NA, NA)[0, ]   # no causal estimate")
ADAPTER_COOK_V11 <- .mrpkg_adapter("cook", "Cook's distance (IVW)", 3L,
  .num_p("cook_mult", "double", 4, 1, 20, "Cut-off multiplier", "A variant is flagged when its Cook's distance is above this number divided by the number of variants (4 / k, the usual cut-off).", "", "Cook 1977; the 4 / k rule of thumb.", c("diagnostics")), FALSE, c(.infl_lines,
  "cd <- stats::cooks.distance(fit); k <- nrow(dat); cut <- cook_mult / k; fl <- dat$variant_id[cd > cut]",
  "diagnostics <- rbind(drow(\"cook_outliers\", length(fl), note = if (length(fl)) paste(fl, collapse = \";\") else NA_character_), drow(\"cook_cutoff\", cut, note = sprintf(\"%s / %d\", format(cook_mult), k)),",
  "  drow(\"cook_max\", max(cd), note = dat$variant_id[which.max(cd)]))"), pkg = "stats")
ADAPTER_STUDRES_V11 <- .mrpkg_adapter("studres", "Studentized residuals (IVW)", 3L,
  .enum_p("stud_cut", "2", c("2", "3"), c("2" = "|r| above 2 (about 5% expected by chance)", "3" = "|r| above 3 (stricter)"), "Cut-off", "Variants whose externally studentized residual is larger than this are flagged.", "", "Usual cut-offs for studentized residuals.", c("diagnostics")), FALSE, c(.infl_lines,
  "rs <- stats::rstudent(fit); fl <- dat$variant_id[abs(rs) > as.numeric(stud_cut)]",
  "diagnostics <- rbind(drow(\"studres_outliers\", length(fl), note = if (length(fl)) paste(fl, collapse = \";\") else NA_character_), drow(\"studres_cutoff\", as.numeric(stud_cut)),",
  "  drow(\"studres_max\", max(abs(rs)), note = dat$variant_id[which.max(abs(rs))]))"), pkg = "stats")

# ---- Round 109: MR-Clust (Foley et al. 2021, Bioinformatics 37:531-541, PMID 32915962) ----
# Clusters the variants by their ratio estimates (by / bx, SE se_y / |bx|) with a normal mixture fitted by EM, plus a null cluster
# (ratio 0) and a junk cluster; the number of clusters is chosen by BIC. As in the package vignette, a variant counts towards a cluster
# when its allocation probability is at least mc_prob and the cluster has at least mc_min such variants (pr_clust). The package gives
# no standard error for a cluster mean, so each cluster is a diagnostic (value = cluster mean, df = number of variants, note = members),
# not a causal estimate. mrclust 0.1.0 (GitHub cnfoley/mrclust, commit 039ed855, GPL-3, pure R), installed by SHA-256 in the build.
PIN_MRCLUST_V11 <- "039ed855bda28b76195222a8ddd2b7b5041b775b"
ADAPTER_MRCLUST_V11 <- .mrpkg_adapter("mrclust", "MR-Clust", 4L, c(
  .num_p("mc_prob", "double", 0.8, 0.5, 0.99, "Allocation probability to keep a variant", "A variant counts towards a cluster only when its probability of belonging to that cluster is at least this value.", "prob",
         "mrclust vignette and Foley et al. 2021: pr_clust(prob = 0.8).", c("diagnostics"), advanced = FALSE),
  .num_p("mc_min", "integer", 4L, 1L, 50L, "Smallest cluster reported", "Clusters with fewer kept variants than this are not reported.", "min_obs", "mrclust vignette and Foley et al. 2021: pr_clust(min_obs = 4).", c("diagnostics"))), TRUE, c(
  "theta <- dat$beta_y / dat$beta_x; theta_se <- dat$se_y / abs(dat$beta_x)   # ratio estimates and their first-order standard errors",
  "fit <- mrclust::mr_clust_em(theta = theta, theta_se = theta_se, bx = dat$beta_x, by = dat$beta_y, bxse = dat$se_x, byse = dat$se_y, obs_names = dat$variant_id, plot_results = NULL)",
  "best <- as.data.frame(fit$results$best); kept <- as.data.frame(mrclust::pr_clust(dta = fit$results$best, prob = mc_prob, min_obs = mc_min))",
  "sub <- kept[!kept$cluster_class %in% c(\"Null\", \"Junk\"), , drop = FALSE]",
  "cl <- unique(sub[, c(\"cluster_class\", \"cluster_mean\")]); cl <- cl[order(cl$cluster_mean), , drop = FALSE]",
  "estimates <- est1(NA, NA, NA, NA, NA, NA)[0, ]   # clusters are reported as diagnostics (no standard error)",
  "mem <- lapply(seq_len(nrow(cl)), function(i) { k <- cl$cluster_class[i]; drow(sprintf(\"mrclust_cluster_%d\", i), cl$cluster_mean[i], sum(sub$cluster_class == k), note = paste(sub$observation[sub$cluster_class == k], collapse = \";\")) })",
  "diagnostics <- do.call(rbind, c(list(drow(\"mrclust_n_clusters\", nrow(cl), note = sprintf(\"clusters with at least %d variants of probability >= %s; each mrclust_cluster_i row: value = cluster mean, df = variants\", mc_min, format(mc_prob))),",
  "  drow(\"mrclust_null\", sum(best$cluster_class == \"Null\"), note = paste(best$observation[best$cluster_class == \"Null\"], collapse = \";\")),",
  "  drow(\"mrclust_junk\", sum(best$cluster_class == \"Junk\"), note = paste(best$observation[best$cluster_class == \"Junk\"], collapse = \";\")),",
  "  drow(\"mrclust_unassigned\", nrow(dat) - nrow(kept), note = \"below the probability or cluster-size cut-off\")), mem))"),
  pkg = "mrclust")

# ---- Round 109: GSMR (item 6), from the methods helper (validated against the gsmr2 vignette example) ----
# ---- Item 6: GSMR, generalised summary-data-based MR (Zhu et al. 2018, Nat Commun 9:224, PMID 29335400) ------------------
# Every instrument is used together with the LD correlation between the instruments (from the 1000 Genomes panel, PLINK;
# see ld_matrix_here.R): generalised least squares estimate of the causal effect, after removing variants in high LD
# (ld_r2_thresh), and HEIDI-outlier filtering of pleiotropic variants. Optionally the reverse direction too (bi-directional GSMR).
# gsmr2 (GitHub JianYang-Lab/gsmr2, commit e7953cc, GPL (>= 2)) with gsmr2_beta = 0, the original HEIDI-outlier method of the paper.
# Needs survey (pchisqsum, multi-variant HEIDI) and markdown (Depends). Requires reference_adapters_v11.R sourced first.
PIN_GSMR2 <- "e7953cc136ad1cfae14674facde21e90caddea5f"
.gsmr_params <- c(
  .enum_p("gsmr_std", "no", c("no", "yes"), c(no = "As given (effects in the exposure GWAS units)", yes = "Standardise the exposure effects (per SD, from EAF and N)"),
          "Exposure effect scale", "GSMR's own example standardises the exposure effects with std_effect() (needs eaf_x and n_x, continuous exposure only), so the estimate is per SD of the exposure.",
          "std_effect", "gsmr2::std_effect; Zhu et al. 2018 (PMID 29335400)."),
  .enum_p("gsmr_dir", "forward", c("forward", "both"), c(forward = "Exposure on outcome", both = "Both directions (bi-directional GSMR)"),
          "Direction", "Both directions also estimates the outcome on the exposure with the same variants (reported as a diagnostic).", "bi_gsmr", "gsmr2::bi_gsmr.", c("diagnostics")),
  .enum_p("gsmr_heidi", "on", c("on", "off"), c(on = "Remove HEIDI outliers (package default)", off = "Keep every variant"),
          "HEIDI-outlier filter", "Removes variants whose effects are heterogeneous with the top variant (pleiotropic outliers).", "heidi_outlier_flag", "gsmr2::gsmr default (heidi_outlier_flag = TRUE)."),
  .num_p("gsmr_gwas_p", "double", 5e-8, 1e-300, 1, "Instrument P threshold", "Variants with an exposure P above this are left out (as weak). MR Studio has already selected the instruments; the package default is 5e-8.", "gwas_thresh", "gsmr2::gsmr default (gwas_thresh = 5e-8)."),
  .num_p("gsmr_heidi_single", "double", 0.01, 1e-10, 0.5, "Single-variant HEIDI P threshold", "Variants with a single-variant HEIDI P below this are removed.", "single_snp_heidi_thresh", "gsmr2::gsmr default (0.01)."),
  .num_p("gsmr_heidi_multi", "double", 0.01, 1e-10, 0.5, "Multi-variant HEIDI P threshold", "Threshold for the multi-variant HEIDI test.", "multi_snps_heidi_thresh", "gsmr2::gsmr default (0.01)."),
  .num_p("gsmr_nsnps", "integer", 10L, 3L, 1000L, "Minimum number of variants", "GSMR stops when fewer variants remain after the filters. The package default is 10.", "nsnps_thresh", "gsmr2::gsmr default (nsnps_thresh = 10)."),
  .num_p("gsmr_ld_r2", "double", 0.05, 0.001, 1, "LD r2 threshold", "Of two variants in LD above this, one is removed.", "ld_r2_thresh", "gsmr2::gsmr default (ld_r2_thresh = 0.05)."),
  .num_p("gsmr_ld_fdr", "double", 0.05, 1e-6, 1, "LD FDR threshold", "Small LD correlations that may be chance (given the size of the reference panel) are set to 0 at this FDR.", "ld_fdr_thresh", "gsmr2::gsmr default (ld_fdr_thresh = 0.05)."),
  .num_p("gsmr_n_ref", "integer", 503L, 50L, 1000000L, "Reference panel size", "Number of people in the LD reference. 503 is the 1000 Genomes European panel; MR Studio sets it from the panel used.", "n_ref", "1000 Genomes phase 3: EUR 503, EAS 504, AFR 661, AMR 347, SAS 489.", c("diagnostics")))
# One body, run by the app (evaluated in an environment) and written into the script. It sees dat, ld (LD matrix with the rsIDs as
# names, signed to the effect alleles), the settings, ci_level, est1() and drow().
.gsmr_body <- c(
  "ids <- as.character(dat$variant_id); ld <- ld[ids, ids, drop = FALSE]   # LD between the instruments, in the order of dat",
  "bzx <- dat$beta_x; bzx_se <- dat$se_x",
  "if (identical(gsmr_std, \"yes\")) { s <- gsmr2::std_effect(dat$eaf_x, dat$beta_x, dat$se_x, dat$n_x); bzx <- s$b; bzx_se <- s$se }   # per SD of the exposure",
  "px <- if (\"p_x\" %in% names(dat)) dat$p_x else 2 * stats::pnorm(-abs(dat$beta_x / dat$se_x)); py <- if (\"p_y\" %in% names(dat)) dat$p_y else 2 * stats::pnorm(-abs(dat$beta_y / dat$se_y))",
  "args <- list(bzx, bzx_se, px, dat$beta_y, dat$se_y, py, ld, ids, gsmr_n_ref, heidi_outlier_flag = identical(gsmr_heidi, \"on\"), gwas_thresh = gsmr_gwas_p,",
  "  single_snp_heidi_thresh = gsmr_heidi_single, multi_snps_heidi_thresh = gsmr_heidi_multi, nsnps_thresh = gsmr_nsnps, ld_r2_thresh = gsmr_ld_r2, ld_fdr_thresh = gsmr_ld_fdr, gsmr2_beta = 0)",
  "fit <- do.call(gsmr2::gsmr, args)   # forward GSMR (exposure on outcome)",
  "z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  "estimates <- est1(fit$bxy, fit$bxy_se, fit$bxy - z * fit$bxy_se, fit$bxy + z * fit$bxy_se, fit$bxy_pval, length(fit$used_index))",
  "lst <- function(x) if (length(x)) paste(x, collapse = \";\") else NA_character_",
  "diagnostics <- rbind(drow(\"gsmr_n_used\", length(fit$used_index), note = lst(ids[fit$used_index])), drow(\"gsmr_n_weak\", length(fit$weak_snps), note = lst(fit$weak_snps)),",
  "  drow(\"gsmr_n_linkage\", length(fit$linkage_snps), note = lst(fit$linkage_snps)), drow(\"gsmr_n_heidi_outliers\", length(fit$pleio_snps), note = lst(fit$pleio_snps)),",
  "  drow(\"gsmr_scale\", NA, note = if (identical(gsmr_std, \"yes\")) \"per SD of the exposure (standardised with std_effect)\" else \"per unit of the exposure GWAS\"))",
  "if (identical(gsmr_dir, \"both\")) {",
  "  rv <- tryCatch(do.call(gsmr2::bi_gsmr, args), error = function(e) e)   # reverse: outcome on exposure, same variants",
  "  diagnostics <- rbind(diagnostics, if (inherits(rv, \"error\")) drow(\"gsmr_reverse\", NA, note = conditionMessage(rv)) else",
  "    drow(\"gsmr_reverse\", rv$reverse_bxy, p = rv$reverse_bxy_pval, note = sprintf(\"outcome on exposure %.4g (SE %.4g), %d variants\", rv$reverse_bxy, rv$reverse_bxy_se, length(rv$reverse_index))))",
  "}")
.gsmr_check <- function(input, params) {
  d <- .usable(input$data); ld <- input$ld
  if (is.null(ld) || !is.matrix(ld)) return(verdict("block", "FEAS_NEEDS_LD", "GSMR needs the LD correlation between the instruments. Install the 1000 Genomes LD reference panel for this ancestry (Settings > Data and access)."))
  miss <- setdiff(d$variant_id, rownames(ld))
  if (length(miss)) return(verdict("block", "FEAS_LD_MISSING_VARIANTS", sprintf("%d instrument(s) are not in the LD reference panel (%s); leave them out or use another panel.", length(miss), paste(utils::head(miss, 5), collapse = ", "))))
  if (identical(params$gsmr_std, "yes") && !(all(c("eaf_x", "n_x") %in% names(d)) && all(is.finite(d$eaf_x) & d$eaf_x > 0 & d$eaf_x < 1 & is.finite(d$n_x))))
    return(verdict("block", "FEAS_NEEDS_EAF_N", "Standardising the exposure effects needs the effect allele frequency (eaf_x) and the sample size (n_x) of every variant. Choose 'As given' or use a dataset that reports them."))
  if (identical(params$gsmr_std, "yes") && identical(input$exposure_scale %||% "beta", "log_or")) return(verdict("block", "FEAS_STD_BINARY", "std_effect standardises continuous exposures only; this exposure is binary (log odds ratios). Choose 'As given'."))
  if (nrow(d) < params$gsmr_nsnps) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("GSMR needs at least %d variants (setting: minimum number of variants); %d available.", params$gsmr_nsnps, nrow(d))))
  w <- .common_warnings(d); verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
}
ADAPTER_GSMR_V11 <- NULL
ADAPTER_GSMR_V11 <- mrgui_adapter(
  id = "gsmr", adapter_version = "1.1.0", label = "GSMR (with HEIDI-outlier)", input_kind = "pair_iv_ld", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "ld_matrix", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(gsmr2 = PIN_GSMR2), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .gsmr_params, .minf_param_v11),
  check = .gsmr_check,
  run = function(input, params, ctx) {
    e <- new.env(parent = globalenv()); for (k in names(params)) assign(k, params[[k]], envir = e)
    d <- .usable(input$data); d <- d[(d$beta_x / d$se_x)^2 >= params$min_f, , drop = FALSE]
    e$dat <- d; e$ld <- input$ld; e$EXPOSURE_ID <- input$exposure_id; e$OUTCOME_ID <- input$outcome_id; e$SCALE <- input$outcome_scale
    eval(parse(text = c(.mr_body_common, .gsmr_body)), envir = e)
    list(estimates = e$estimates, diagnostics = e$diagnostics)
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_GSMR_V11, spec, params, seed, "gsmr2", c(
    sprintf("EXPOSURE_ID <- %s; OUTCOME_ID <- %s; SCALE <- %s", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale)),
    sprintf("ld <- as.matrix(utils::read.csv(file.path(\"01_data\", %s), row.names = 1, check.names = FALSE))   # LD matrix (1000 Genomes, PLINK), signed to the effect alleles", .q(spec$ld_file %||% "ld_matrix.csv")),
    "# ---- 5. Analysis: GSMR (gsmr2, original HEIDI-outlier method) ----", .mr_body_common, .gsmr_body),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = c(.colreq_basic[1:4], list(
    list(column = "eaf_x", level = "recommended", why = "Needed only to standardise the exposure effects (per SD).", if_missing = list(action = "warn", code = "FEAS_NEEDS_EAF_N")),
    list(column = "n_x", level = "recommended", why = "Needed only to standardise the exposure effects (per SD).", if_missing = list(action = "warn", code = "FEAS_NEEDS_EAF_N")))))

# Figure, laid out as the GSMR example figure (vignette): exposure effect against outcome effect for the variants GSMR used, with
# +/- 1 SE bars in both directions and a dashed line through the origin with slope bxy. `used` = rsIDs in gsmr_n_used.

# ---- Round 109: grapple (item 18), from the methods helper (validated by simulation with known truth) ----
# ---- Item 18: GRAPPLE, robust profile-likelihood MR, univariable and multivariable --------------------------------------
# (Wang, Zhao, Bowden, Hemani, Davey Smith, Small, Zhang 2021, PLoS Genet 17:e1009575, PMID 34157017, doi 10.1371/journal.pgen.1009575).
# Robustified adjusted profile likelihood with over-dispersion tau2 (pervasive balanced pleiotropy), robust to weak instruments;
# with K exposures it is multivariable MR; cor.mat is the correlation of the GWAS estimation errors (shared participants:
# sample overlap), order exposures 1..K then the outcome. With one exposure, no overlap and the same loss it is MR-RAPS.
# GRAPPLE (GitHub jingshuw/GRAPPLE, commit 317e837, GPL (>= 2)). Only grappleRobustEst() is used; getInput() (needs full GWAS
# files and PLINK clumping) and the marker annotation (haploR, a web query) are not. Requires reference_adapters_v11.R sourced first.
PIN_GRAPPLE <- "317e837340a29129c52c9d9b22f5f8bae72e27d0"
.grapple_loss_p <- .enum_p("gr_loss", "tukey", c("tukey", "huber", "l2"), c(tukey = "Tukey's biweight (robust, package default)", huber = "Huber (robust)", l2 = "Squared error (not robust)"),
  "Loss function", "Robust losses limit the influence of outlying (pleiotropic) variants.", "loss.function", "GRAPPLE::grappleRobustEst default (loss.function = \"tukey\", k = 4.685).")
.grapple_rxy_p <- .num_p("gr_rho_xy", "double", 0, -0.99, 0.99, "Error correlation, exposure(s) and outcome", "Correlation of the estimation errors of the exposure and outcome GWAS when they share participants (sample overlap); 0 for separate samples. GRAPPLE estimates it from null variants of the full GWAS files; enter that value here.", "cor.mat", "Wang et al. 2021 (PMID 34157017); 0 for non-overlapping samples.", advanced = FALSE)
# One body for both kinds. It sees bx, sx (matrices, one column per exposure), by, sy, ids, Rx (K x K exposure error correlations),
# gr_rho_xy, gr_loss, ci_level, est_rows() (direct or causal) and drow().
.grapple_body <- c(
  "K <- ncol(bx); g <- data.frame(SNP = ids, stringsAsFactors = FALSE)",
  "for (k in seq_len(K)) { g[[paste0(\"gamma_exp\", k)]] <- bx[, k]; g[[paste0(\"se_exp\", k)]] <- sx[, k] }",
  "g$gamma_out1 <- by; g$se_out1 <- sy",
  "C <- diag(K + 1); C[seq_len(K), seq_len(K)] <- Rx; C[K + 1, seq_len(K)] <- C[seq_len(K), K + 1] <- gr_rho_xy   # GWAS error correlations",
  "fit <- GRAPPLE::grappleRobustEst(g, cor.mat = C, loss.function = gr_loss, diagnosis = TRUE, plot.it = FALSE)",
  "se <- sqrt(diag(as.matrix(fit$beta.var))); z <- stats::qnorm(1 - (1 - ci_level) / 2)",
  "estimates <- est_rows(fit$beta.hat, se, fit$beta.hat - z * se, fit$beta.hat + z * se, fit$beta.p.value, length(ids))",
  "out <- if (is.data.frame(fit$outliers) && nrow(fit$outliers)) ids[vapply(fit$outliers$sample, function(v) which.min(abs(fit$std.resid - v)), 1L)] else character(0)   # QQ-plot outlier points back to their variants",
  "diagnostics <- rbind(drow(\"grapple_tau2\", fit$tau2.hat, note = sprintf(\"over-dispersion (pleiotropy variance), SE %.4g\", fit$tau2.se)),",
  "  drow(\"grapple_normality_p\", NA, p = if (length(fit$std.resid) %in% 3:5000) stats::shapiro.test(fit$std.resid)$p.value else NA_real_, note = \"test of the standardised residuals (Shapiro-Wilk, as GRAPPLE prints); a small P suggests outliers or several pleiotropic pathways\"),",
  "  drow(\"grapple_outliers\", length(out), note = if (length(out)) paste(out, collapse = \";\") else NA_character_))")
.grapple_uv_prep <- c(
  "bx <- matrix(dat$beta_x); sx <- matrix(dat$se_x); by <- dat$beta_y; sy <- dat$se_y; ids <- dat$variant_id; Rx <- diag(1)",
  "est_rows <- function(b, se, lo, hi, p, n) est1(b, se, lo, hi, p, n)")
.grapple_mv_prep <- c(
  "K <- sum(grepl(\"^beta_x_[0-9]+$\", names(dat))); bxn <- sprintf(\"beta_x_%d\", seq_len(K)); sxn <- sprintf(\"se_x_%d\", seq_len(K))",
  "use <- dat$mr_keep & is.finite(dat$beta_y) & is.finite(dat$se_y) & dat$se_y > 0 & stats::complete.cases(dat[, c(bxn, sxn)]); dat <- dat[use, , drop = FALSE]",
  "EXPOSURE_IDS <- strsplit(EXPOSURE_ID, \"|\", fixed = TRUE)[[1]]",
  "bx <- as.matrix(dat[, bxn]); sx <- as.matrix(dat[, sxn]); by <- dat$beta_y; sy <- dat$se_y; ids <- dat$variant_id",
  "Rx <- diag(K); if (K > 1) for (i in seq_len(K - 1)) for (j in (i + 1):K) Rx[i, j] <- Rx[j, i] <- get(sprintf(\"mv_r_%d_%d\", i, j))   # exposure error correlations",
  "est_rows <- function(b, se, lo, hi, p, n) data.frame(estimand = \"direct_effect\", exposure_id = EXPOSURE_IDS, outcome_id = OUTCOME_ID, scale = SCALE,",
  "  estimate = as.numeric(b), se = as.numeric(se), ci_lower = as.numeric(lo), ci_upper = as.numeric(hi), ci_level = ci_level, ci_method = \"wald\",",
  "  p = as.numeric(p), n_ivs_used = as.integer(n), stringsAsFactors = FALSE)",
  .mr_body_common[5:6])   # drow() only: the univariable mr_input line does not apply
.grapple_run <- function(prep, mv) function(input, params, ctx) {
  e <- new.env(parent = globalenv()); for (k in names(params)) assign(k, params[[k]], envir = e)
  d <- if (mv) input$data else { d <- .usable(input$data); d[(d$beta_x / d$se_x)^2 >= params$min_f, , drop = FALSE] }
  e$dat <- d; e$EXPOSURE_ID <- input$exposure_id; e$OUTCOME_ID <- input$outcome_id; e$SCALE <- input$outcome_scale
  eval(parse(text = c(if (!mv) .mr_body_common, prep, .grapple_body)), envir = e)
  list(estimates = e$estimates, diagnostics = e$diagnostics)
}
ADAPTER_GRAPPLE_V11 <- NULL
ADAPTER_GRAPPLE_V11 <- mrgui_adapter(
  id = "grapple", adapter_version = "1.1.0", label = "GRAPPLE", input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(GRAPPLE = PIN_GRAPPLE), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .grapple_loss_p, .grapple_rxy_p, .minf_param_v11),
  check = .mr_pkg_check("GRAPPLE", 3L),
  run = .grapple_run(.grapple_uv_prep, FALSE),
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_GRAPPLE_V11, spec, params, seed, "GRAPPLE", c(
    sprintf("EXPOSURE_ID <- %s; OUTCOME_ID <- %s; SCALE <- %s", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale)),
    "# ---- 5. Analysis: GRAPPLE (grappleRobustEst) ----", .mr_body_common, .grapple_uv_prep, .grapple_body),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = .colreq_basic)
# Multivariable: the same body on K exposures (multi_iv input, columns beta_x_k / se_x_k as the other MVMR adapters).
ADAPTER_MVGRAPPLE_V11 <- NULL
ADAPTER_MVGRAPPLE_V11 <- mrgui_adapter(
  id = "mvgrapple", adapter_version = "1.1.0", label = "MVMR GRAPPLE", input_kind = "multi_iv", estimands = "direct_effect",
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE, needs_full_sumstats = FALSE, min_exposures = 2, max_exposures = 6),
  packages = list(GRAPPLE = PIN_GRAPPLE), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .grapple_loss_p, .grapple_rxy_p, .mv_cor_params),
  check = function(input, params) {
    K <- sum(grepl("^beta_x_[0-9]+$", names(input$data))); n <- nrow(.mv_usable(input$data))
    if (K < 2) return(verdict("block", "FEAS_TOO_FEW_EXPOSURES", sprintf("MVMR GRAPPLE needs at least 2 exposures; %d given.", K)))
    if (n < K + 2) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("MVMR GRAPPLE needs at least %d variants with all %d exposures; %d available.", K + 2, K, n)))
    R <- diag(K + 1); for (i in seq_len(K - 1)) for (j in (i + 1):K) R[i, j] <- R[j, i] <- as.numeric(params[[sprintf("mv_r_%d_%d", i, j)]] %||% 0)
    R[K + 1, 1:K] <- R[1:K, K + 1] <- params$gr_rho_xy
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(verdict("block", "FEAS_INVALID_CORRELATION", "The error correlations do not form a valid correlation matrix (it must be positive definite). Check the values."))
    verdict("ok", character(), character())
  },
  run = .grapple_run(.grapple_mv_prep, TRUE),
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_MVGRAPPLE_V11, spec, params, seed, "GRAPPLE", c(
    sprintf("EXPOSURE_ID <- %s; OUTCOME_ID <- %s; SCALE <- %s", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale)),
    "# ---- 5. Analysis: multivariable GRAPPLE (grappleRobustEst) ----", .grapple_mv_prep, .grapple_body),
    "fixed seed (this method is deterministic; kept for the record)", mv = TRUE),   # Round 110a: the multivariable data head
  column_requirements = .mv_colreq)

# ---- Round 109: ppmr (item 19), from the methods helper (validated by simulation with known truth) ----
# ---- Item 19: PPMR, probabilistic MR with correlated instruments and horizontal pleiotropy ------------------------------
# (PMR-Egger: Yuan et al. 2020, Nat Commun 11:3861, PMID 32737316, doi 10.1038/s41467-020-17668-6; summary-data version).
# Likelihood model with the LD between the instruments (Sigma1 = Sigma2 = the 1000 Genomes LD matrix, PLINK, see
# ../gsmr2/ld_matrix_here.R), polygenic variant-exposure effects, a causal effect alpha and a horizontal pleiotropic effect gamma
# shared by the variants (an Egger-type direct effect). PPMR 1.0.1 (CRAN, GPL-3) exposes only the fit (PMR_summary); the two
# P values are likelihood ratio tests, as in the paper: refit with alpha = 0 and with gamma = 0, 2 x (loglik difference) on 1 df.
# PMR_summary works on the standardised scale (betax = X'y / (n - 1) for standardised genotypes and phenotypes), so the
# associations are entered as z / sqrt(n): N of both GWAS is required, and alpha is per SD of the exposure on SD of the outcome.
# The package gives no standard error, so the estimate has no confidence interval (as the paper, which reports P values).
# Round 110a: alpha is converted back to the pair's own scale (beta or log OR per unit of the exposure, as the other methods) with
# k = median over the variants of (se_y sqrt(n_y)) / (se_x sqrt(n_x)): z / sqrt(n) is beta / (se sqrt(n)), so by = alpha bx on the
# standardised scale is beta_y = alpha k beta_x when k is the same for every variant (it is, up to sampling, as se sqrt(n) is about
# SD(trait) / SD(genotype)). k is saved as the diagnostic ppmr_scale; gamma stays on the standardised scale (it is not a slope).
# Requires reference_adapters_v11.R sourced first.
PIN_PPMR <- "1.0.1"
.ppmr_params <- c(
  .num_p("ppmr_max_iter", "integer", 1000L, 50L, 100000L, "Maximum iterations", "Maximum number of EM iterations of each fit.", "max_iterin", "PPMR paper and the original PPMR package (max_iterin = 1000)."),
  .num_p("ppmr_eps", "double", 1e-5, 1e-10, 1e-2, "Convergence tolerance", "The EM stops when the log-likelihood changes by less than this.", "epsin", "Original PPMR package default (epsin = 1e-5)."))
.ppmr_body <- c(
  "ids <- as.character(dat$variant_id); R <- ld[ids, ids, drop = FALSE]   # LD between the instruments, in the order of dat",
  "nx <- round(stats::median(dat$n_x)); ny <- round(stats::median(dat$n_y))   # sample sizes of the two GWAS",
  "bx <- dat$beta_x / dat$se_x / sqrt(dat$n_x); by <- dat$beta_y / dat$se_y / sqrt(dat$n_y)   # standardised scale, z / sqrt(n)",
  "pmr <- function(g0, a0) PPMR::PMR_summary(bx, by, R, R, nx, ny, g0, a0, ppmr_max_iter, ppmr_eps)",
  "f <- pmr(0, 0); fa <- pmr(0, 1); fg <- pmr(1, 0)   # full model, alpha = 0, gamma = 0",
  "lrt <- function(l1, l0) stats::pchisq(max(0, 2 * (l1 - l0)), 1, lower.tail = FALSE)",
  "p_alpha <- lrt(f$loglik, fa$loglik); p_gamma <- lrt(f$loglik, fg$loglik)",
  "k <- stats::median((dat$se_y * sqrt(dat$n_y)) / (dat$se_x * sqrt(dat$n_x)))   # back to the pair's scale (beta or log OR per unit)",
  "estimates <- est1(f$alpha * k, NA, NA, NA, p_alpha, nrow(dat), \"profile\")",
  "diagnostics <- rbind(drow(\"ppmr_gamma\", f$gamma, df = 1, p = p_gamma, note = \"horizontal pleiotropic effect (likelihood ratio test of gamma = 0)\"),",
  "  drow(\"ppmr_iterations\", f$iteration, note = sprintf(\"EM iterations of the full model (maximum %d)\", ppmr_max_iter)),",
  "  drow(\"ppmr_loglik\", f$loglik, note = sprintf(\"log-likelihood; alpha = 0: %.6g; gamma = 0: %.6g\", fa$loglik, fg$loglik)),",
  "  drow(\"ppmr_scale\", k, note = sprintf(\"fitted on the standardised scale (z / sqrt(n)), alpha %.6g per SD on SD; the estimate is alpha x k, k = median (se_y sqrt(n_y)) / (se_x sqrt(n_x)), on the scale of the other methods; gamma is on the standardised scale; N exposure %g, outcome %g\", f$alpha, nx, ny)))")
.ppmr_check <- function(input, params) {
  d <- .usable(input$data); ld <- input$ld
  if (is.null(ld) || !is.matrix(ld)) return(verdict("block", "FEAS_NEEDS_LD", "PPMR needs the LD correlation between the instruments. Install the 1000 Genomes LD reference panel for this ancestry (Settings > Data and access)."))
  miss <- setdiff(d$variant_id, rownames(ld))
  if (length(miss)) return(verdict("block", "FEAS_LD_MISSING_VARIANTS", sprintf("%d instrument(s) are not in the LD reference panel (%s).", length(miss), paste(utils::head(miss, 5), collapse = ", "))))
  if (!all(c("n_x", "n_y") %in% names(d)) || !all(is.finite(d$n_x) & is.finite(d$n_y))) return(verdict("block", "FEAS_NEEDS_N", "PPMR needs the sample sizes of the exposure and outcome GWAS (n_x, n_y) to put the associations on the standardised scale."))
  if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS", sprintf("PPMR needs at least 3 variants; %d available.", nrow(d))))
  w <- .common_warnings(d); verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
}
ADAPTER_PPMR_V11 <- NULL
ADAPTER_PPMR_V11 <- mrgui_adapter(
  id = "ppmr", adapter_version = "1.1.0", label = "PPMR (correlated instruments, pleiotropy)", input_kind = "pair_iv_ld", estimands = "causal_effect",
  requires = list(min_ivs = 3, needs_ld = "ld_matrix", needs_eaf = FALSE, needs_n = TRUE, needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(PPMR = PIN_PPMR), stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = c(.ci_param_v11, .ppmr_params, .minf_param_v11),
  check = .ppmr_check,
  run = function(input, params, ctx) {
    e <- new.env(parent = globalenv()); for (k in names(params)) assign(k, params[[k]], envir = e)
    d <- .usable(input$data); d <- d[(d$beta_x / d$se_x)^2 >= params$min_f, , drop = FALSE]
    e$dat <- d; e$ld <- input$ld; e$EXPOSURE_ID <- input$exposure_id; e$OUTCOME_ID <- input$outcome_id; e$SCALE <- input$outcome_scale
    eval(parse(text = c(.mr_body_common, .ppmr_body)), envir = e)
    list(estimates = e$estimates, diagnostics = e$diagnostics)
  },
  emit_script = function(spec, params, seed) .detector_emit(ADAPTER_PPMR_V11, spec, params, seed, "PPMR", c(
    sprintf("EXPOSURE_ID <- %s; OUTCOME_ID <- %s; SCALE <- %s", .q(spec$exposure_id), .q(spec$outcome_id), .q(spec$outcome_scale)),
    sprintf("ld <- as.matrix(utils::read.csv(file.path(\"01_data\", %s), row.names = 1, check.names = FALSE))   # LD matrix (1000 Genomes, PLINK), signed to the effect alleles", .q(spec$ld_file %||% "ld_matrix.csv")),
    "# ---- 5. Analysis: PPMR (PMR_summary, likelihood ratio tests) ----", .mr_body_common, .ppmr_body),
    "fixed seed (this method is deterministic; kept for the record)"),
  column_requirements = c(.colreq_basic[1:4], list(
    list(column = "n_x", level = "required", why = "Sample size of the exposure GWAS (standardised scale).", if_missing = list(action = "block", code = "FEAS_NEEDS_N")),
    list(column = "n_y", level = "required", why = "Sample size of the outcome GWAS (standardised scale).", if_missing = list(action = "block", code = "FEAS_NEEDS_N")))))
