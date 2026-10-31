# Reference adapters (C6). Numerics follow TwoSampleMR commit 5e93ee7 (R/mr.R):
# mr_ivw_mre() and mr_egger_regression(). Conformance is tested against those
# functions in tests/test_contracts.R.

.usable <- function(d) {
  ok <- d$mr_keep & is.finite(d$beta_x) & is.finite(d$beta_y) &
        is.finite(d$se_x) & is.finite(d$se_y) & d$se_y > 0
  d[ok, , drop = FALSE]
}

.common_warnings <- function(d) {
  codes <- character(); msgs <- character()
  f <- (d$beta_x / d$se_x)^2
  if (nrow(d) && mean(f) < 10) {
    codes <- "DATA_WEAK_IVS"
    msgs <- sprintf("Mean approximate F-statistic is %.1f (below 10).", mean(f))
  }
  list(codes = codes, messages = msgs)
}

.ci_param <- list(ci_level = list(type = "double", default = 0.95, min = 0.5, max = 0.999))

ADAPTER_IVW_MRE <- mrgui_adapter(
  id = "ivw_mre", adapter_version = "1.0.0",
  label = "Inverse-variance weighted (multiplicative random effects)",
  input_kind = "pair_iv", estimands = "causal_effect",
  requires = list(min_ivs = 2, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"),
  stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = .ci_param,
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 2) return(verdict("block", "FEAS_TOO_FEW_IVS",
      sprintf("IVW needs at least 2 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    fit <- summary(stats::lm(beta_y ~ -1 + beta_x, data = d, weights = 1 / d$se_y^2))
    b <- fit$coefficients["beta_x", "Estimate"]
    se <- fit$coefficients["beta_x", "Std. Error"]
    z <- stats::qnorm(1 - (1 - params$ci_level) / 2)
    q_df <- nrow(d) - 1; q <- fit$sigma^2 * q_df
    list(
      estimates = estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
        b, se, b - z * se, b + z * se, "wald", 2 * stats::pnorm(abs(b / se), lower.tail = FALSE),
        nrow(d), params$ci_level),
      diagnostics = diagnostic_row("cochran_q", q, q_df, stats::pchisq(q, q_df, lower.tail = FALSE)))
  })

ADAPTER_EGGER <- mrgui_adapter(
  id = "egger", adapter_version = "1.0.0", label = "MR-Egger",
  input_kind = "pair_iv", estimands = c("causal_effect", "pleiotropy_intercept"),
  requires = list(min_ivs = 3, needs_ld = "none", needs_eaf = FALSE, needs_n = FALSE,
                  needs_full_sumstats = FALSE, min_exposures = 1, max_exposures = 1),
  packages = list(TwoSampleMR = "5e93ee7d8e1f67581acb8134299759389d73efba"),
  stochastic = FALSE, supports_checkpoint = FALSE, runtime_class = "fast",
  params = .ci_param,
  check = function(input, params) {
    d <- .usable(input$data)
    if (nrow(d) < 3) return(verdict("block", "FEAS_TOO_FEW_IVS",
      sprintf("MR-Egger needs at least 3 usable instruments; %d available.", nrow(d))))
    w <- .common_warnings(d)
    verdict(if (length(w$codes)) "warn" else "ok", w$codes, w$messages)
  },
  run = function(input, params, ctx) {
    d <- .usable(input$data)
    s <- ifelse(d$beta_x == 0, 1, sign(d$beta_x))      # C5.2: orient to positive beta_x
    d$beta_y <- d$beta_y * s; d$beta_x <- abs(d$beta_x)
    fit <- summary(stats::lm(beta_y ~ beta_x, data = d, weights = 1 / d$se_y^2))
    if (nrow(stats::coefficients(fit)) < 2) stop("Collinear instruments: MR-Egger slope not estimable.")
    k <- nrow(d); dfree <- k - 2; scale_se <- min(1, fit$sigma)
    b <- fit$coefficients[2, 1]; se <- fit$coefficients[2, 2] / scale_se
    bi <- fit$coefficients[1, 1]; sei <- fit$coefficients[1, 2] / scale_se
    tq <- stats::qt(1 - (1 - params$ci_level) / 2, dfree)
    pt2 <- function(est, s) 2 * stats::pt(abs(est / s), dfree, lower.tail = FALSE)
    q <- fit$sigma^2 * dfree
    list(
      estimates = rbind(
        estimate_row("causal_effect", input$exposure_id, input$outcome_id, input$outcome_scale,
          b, se, b - tq * se, b + tq * se, "wald_t", pt2(b, se), k, params$ci_level),
        estimate_row("pleiotropy_intercept", input$exposure_id, input$outcome_id, input$outcome_scale,
          bi, sei, bi - tq * sei, bi + tq * sei, "wald_t", pt2(bi, sei), k, params$ci_level)),
      diagnostics = diagnostic_row("rucker_q", q, dfree, stats::pchisq(q, dfree, lower.tail = FALSE)))
  })
