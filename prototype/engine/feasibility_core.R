# Feasibility engine (Round 21, FR-10): the rules of feasibility_rules.csv (ids from mr_feasibility_rules.csv) applied to the
# facts of every pair. Each method gets a state per study: "ok" (Available), "check" (allowed, with warnings) or "locked",
# with every rule that fired, its message and its source. Rules that need results (stage "result") are evaluated on a saved pair.
# The adapters' own check() stays the final guard at run time; a test keeps the two from disagreeing on a block.

FEAS_METHODS <- c("ivw_mre", "divw", "raps", "egger", "wm", "smed", "mode", "smode", "presso", "radial", "radegger", "steiger")
FEAS_CAUSAL <- setdiff(FEAS_METHODS, "steiger")

#' The rule table (read once).
feas_rules <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      f <- file.path(getOption("mrgui.engine_dir", file.path("prototype", "engine")), "feasibility_rules.csv")
      cache <<- utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
    }
    cache
  }
})

.feas_targets <- function(applies) {
  a <- trimws(strsplit(applies, ",")[[1]])
  unique(unlist(lapply(a, function(x) switch(x, all = FEAS_METHODS, causal = FEAS_CAUSAL, if (x %in% FEAS_METHODS) x else character()))))
}
.fill <- function(msg, v) { for (k in names(v)) msg <- gsub(paste0("{", k, "}"), as.character(v[[k]]), msg, fixed = TRUE); msg }

#' I2GX of a set of instruments (Bowden 2016), with TwoSampleMR::Isq on the absolute exposure effects.
feas_i2gx <- function(bx, sx) { ok <- is.finite(bx) & is.finite(sx) & sx > 0; if (sum(ok) < 2) NA_real_ else TwoSampleMR::Isq(abs(bx[ok]), sx[ok]) }

#' One test: does rule `r` fire for method `m` on pair facts `f`? Returns NULL or a list of message values.
.feas_test <- function(test, r, m, f, ctx) {
  thr <- suppressWarnings(as.numeric(r$threshold))
  switch(test,
    several_pairs = if (ctx$n_pairs > 1) list(n = ctx$n_pairs, bonf = format(signif(0.05 / ctx$n_pairs, 3))),
    ancestry_differs = { ax <- f$anc_x %||% "?"; ay <- f$anc_y %||% "?"; if (!(ax %in% c("?", "")) && !(ay %in% c("?", "")) && !identical(ax, ay)) list(ax = ax, ay = ay) },
    steiger_data = { w <- f$steiger$why %||% ""; if (nzchar(w)) list(why = w) },
    iv_below_min = { mi <- ctx$min_ivs[[m]]; if (f$instruments < mi) list(min = mi, n = f$instruments) },
    iv_below_recommended = if (f$instruments < thr && f$instruments >= ctx$min_ivs[[m]]) list(n = f$instruments, thr = thr,
      median_note = if (m %in% c("wm", "smed") && f$instruments <= 4) " With 3 or 4 instruments the weighted and simple median gave false positives about twice as often as intended in the MR Studio calibration (about 10% instead of 5%)." else ""),   # Round 35
    mean_f_below = if (!(m %in% c("divw", "raps")) && is.finite(f$mean_f) && f$mean_f < thr) list(f = format(round(f$mean_f, 1), nsmall = 1), thr = thr),   # dIVW and MR-RAPS correct this bias (Round 37, 38)
    divw_condition_below = { k <- f$instruments; cn <- if (is.finite(f$mean_f %||% NA_real_) && k > 0) (f$mean_f - 1) * sqrt(k) else NA_real_   # Round 37, Ye 2021
      if (is.finite(cn) && cn < thr) list(c = format(round(cn, 1), nsmall = 1), thr = thr) },
    i2gx_below = if (is.finite(f$i2gx %||% NA_real_) && f$i2gx < thr) {   # Round 35: with overlapping samples the bias goes toward the observational association
      o <- if (isTRUE(f$simulated)) NULL else ctx$overlap[[paste0(f$ex_id, " | ", f$out_id)]]
      ov <- !is.null(o) && (identical(o$answer, "some") || identical(o$answer, "unknown"))
      list(i2 = format(round(f$i2gx, 2), nsmall = 2), thr = thr, bias = if (ov) "the samples may overlap, so the MR-Egger estimate is pulled toward the observational association between exposure and outcome rather than diluted toward the null" else "the MR-Egger estimate is diluted toward the null") },
    i2gx_near = { lo <- suppressWarnings(as.numeric(feas_rules()$threshold[feas_rules()$rule_id == "R24"]))   # Round 32: just above the R24 threshold
      if (is.finite(f$i2gx %||% NA_real_) && length(lo) == 1L && f$i2gx >= lo && f$i2gx < thr) list(i2 = format(round(f$i2gx, 2), nsmall = 2), lo = lo, thr = thr) },
    overlap = {
      if (isTRUE(f$simulated)) NULL else {
        o <- ctx$overlap[[paste0(f$ex_id, " | ", f$out_id)]]
        why <- if (is.null(o)) (if (isTRUE(f$ukb_both)) "Both datasets include UK Biobank and the overlap question is not answered yet." else "The sample overlap question is not answered yet (Check screen).")
               else if (identical(o$answer, "some")) sprintf("You stated that the studies share participants%s.", if (!is.null(o$pct)) sprintf(" (about %s%%)", o$pct) else "")
               else if (identical(o$answer, "unknown")) "You do not know whether the studies share participants."
               else if (isTRUE(f$ukb_both)) "You stated no shared participants, but both datasets include UK Biobank." else NULL
        if (!is.null(why)) list(why = why)
      }
    },
    palindromic_dropped = { n <- f$ledger$dropped_palindromic %||% 0; if (n > 0) list(n = n) },
    palindromic_kept = { n <- if (identical(f$ledger$strand, "forward")) f$ledger$kept_palindromic %||% 0 else 0; if (n > 0) list(n = n) },   # Round 26
    binary_exposure = if (identical(f$scale_x, "log_or")) list(),
    presso_nbdist = { nb <- suppressWarnings(as.numeric(ctx$params$nbdist %||% 1000)); if (nb <= f$instruments) list(nb = nb, n = f$instruments) },
    heavy_runtime = list(),
    NULL)
}

#' Feasibility of every method for the selected pairs. facts = list of pair_facts(); params = the Customise settings;
#' overlap = overlap answers of the study (named by "exposure | outcome"). Returns list(methods = list(id = list(state, reasons))).
feasibility <- function(facts, params = list(), overlap = list(), methods = FEAS_METHODS) {
  rules <- feas_rules(); rules <- rules[rules$active_v1 == "yes" & rules$stage %in% c("selection", "data", "runtime"), , drop = FALSE]
  ad <- param_adapters(); min_ivs <- lapply(ad, function(a) a$requires$min_ivs)
  ctx <- list(n_pairs = length(facts), min_ivs = min_ivs, params = params, overlap = overlap)
  out <- list()
  for (m in methods) {
    reasons <- list()
    for (i in seq_len(nrow(rules))) {
      r <- as.list(rules[i, ]); if (!(m %in% .feas_targets(r$applies_to))) next
      for (f in facts) {
        v <- .feas_test(r$test, r, m, f, ctx); if (is.null(v)) next
        pre <- if (length(facts) > 1 && r$test != "several_pairs") sprintf("%s -> %s: ", f$exposure, f$outcome) else ""
        reasons[[length(reasons) + 1L]] <- list(rule = r$rule_id, level = r$level, msg = paste0(pre, .fill(r$message, v)), recommend = r$recommend, source = r$source)
        if (r$test == "several_pairs" || r$test == "heavy_runtime") break   # study-level rules: once
      }
    }
    lv <- vapply(reasons, `[[`, "", "level")
    out[[m]] <- list(state = if (any(lv == "BLOCK")) "locked" else if (any(lv == "WARN")) "check" else "ok", reasons = reasons)
  }
  list(methods = out, rules_active = sum(feas_rules()$active_v1 == "yes"), rules_total = nrow(feas_rules()))
}

#' Rules that need results, evaluated on a saved pair (payload from pair_payload()).
feasibility_results <- function(pay) {
  rules <- feas_rules(); rules <- rules[rules$active_v1 == "yes" & rules$stage == "result", , drop = FALSE]
  dg <- pay$diagnostics %||% list(); g <- function(m, n) { x <- Filter(function(d) identical(d$method, m) && identical(d$name, n), dg); if (length(x)) x[[1]] else NULL }
  est <- pay$estimates %||% list(); ic <- Filter(function(e) identical(e$method, "egger") && identical(e$estimand, "pleiotropy_intercept"), est)
  out <- list()
  for (i in seq_len(nrow(rules))) {
    r <- as.list(rules[i, ]); thr <- suppressWarnings(as.numeric(r$threshold)); v <- NULL
    if (r$test == "q_significant") { q <- g("ivw_mre", "cochran_q"); if (!is.null(q) && is.finite(as.numeric(q$p)) && as.numeric(q$p) < thr) v <- list(p = .feas_p(q$p)) }
    if (r$test == "egger_intercept" && length(ic)) { p <- as.numeric(ic[[1]]$p); if (is.finite(p) && p < thr) v <- list(p = .feas_p(p)) }
    if (r$test == "steiger_reversed") { d <- g("steiger", "steiger_direction")
      if (!is.null(d) && !(isTRUE(as.numeric(d$value) == 1) && is.finite(as.numeric(d$p)) && as.numeric(d$p) < 0.05))
        v <- list(detail = sprintf(" (direction %s, P = %s)", if (isTRUE(as.numeric(d$value) == 1)) "as assumed but not clearly" else "reversed", .feas_p(d$p))) }
    if (!is.null(v)) out[[length(out) + 1L]] <- list(rule = r$rule_id, level = r$level, msg = .fill(r$message, v), recommend = r$recommend, source = r$source)
  }
  out
}
.feas_p <- function(p) { p <- as.numeric(p); if (!is.finite(p)) "NA" else if (p < 1e-4) formatC(p, digits = 2, format = "e") else formatC(p, digits = 3, format = "f") }
