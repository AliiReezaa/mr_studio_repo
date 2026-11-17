# Round 42: outlier module (requirements_spec.md, Part H, FR-90 to FR-103).
# The detectors are the packaged methods (MR-PRESSO, radial IVW, radial MR-Egger, MR-Lasso) run as part of the pair; this file only
# reads what they flagged, builds the candidate sets, applies the safeguards and the pre-specified selection rule, and writes the
# flag ledger and the comparison table. Every re-analysis is an ordinary pair run with its own scripts (studio_run_without), so no
# statistical step lives here beyond arithmetic on the packages' own results (Q minus df, I squared, excess reduction).
# The primary analysis is never changed: the selected set drives a labelled, data-driven sensitivity analysis (FR-100).

OUTLIER_DETECTORS <- c(presso = "MR-PRESSO", radial = "Radial IVW", radegger = "Radial MR-Egger", lasso = "MR-Lasso")
OUTLIER_FLAG_ROW <- c(presso = "presso_outliers", radial = "radial_outliers", radegger = "radial_outliers", lasso = "lasso_excluded")
OUTLIER_CAP <- 0.20      # FR-95a; kept by the Round 32 calibration
OUTLIER_FLOOR <- 5L      # FR-95b
OUTLIER_TRIGGER <- 0.05  # FR-90; kept by the Round 32 calibration
OUTLIER_REANALYSIS <- c("ivw_mre", "egger", "wm")   # FR-96; MR-PRESSO is added when its global test triggered the module

.o_num <- function(x) { x <- suppressWarnings(as.numeric(x)); if (length(x)) x[[1]] else NA_real_ }
.o_dg <- function(pay, m, n) { for (g in pay$diagnostics) if (identical(g$method, m) && identical(g$name, n)) return(g); NULL }
.o_est <- function(pay, m, estimand = "causal_effect") { for (e in pay$estimates) if (identical(e$method, m) && identical(e$estimand, estimand)) return(e); NULL }
.o_ids <- function(note) { note <- as.character(note %||% ""); if (!length(note) || is.na(note[1]) || !nzchar(note[1]) || identical(note[1], "NA")) character() else strsplit(note[1], ";", fixed = TRUE)[[1]] }

#' The instruments a pair used (mr_keep), in the order of the analysis data.
outlier_instruments <- function(pay) vapply(Filter(function(v) isTRUE(v$keep), pay$variants), function(v) as.character(v$id), "")

#' FR-90: the signals that start the module. Each is TRUE, FALSE or NA (method not run).
outlier_signals <- function(pay) {
  q <- .o_dg(pay, "ivw_mre", "cochran_q"); it <- .o_est(pay, "egger", "pleiotropy_intercept"); g <- .o_dg(pay, "presso", "presso_global_test")
  f <- function(x) if (is.null(x) || !is.finite(.o_num(x$p))) NA else .o_num(x$p) < OUTLIER_TRIGGER
  list(q = f(q), intercept = f(it), presso = f(g))
}
outlier_triggered <- function(sig) names(Filter(isTRUE, sig))

#' FR-92, FR-94, FR-95: one candidate set per detector that ran, plus the consensus set (flagged by at least two detectors).
#' Returns a data frame with one row per set; ids is a list column. The union of all detectors is not offered (it over-removes).
outlier_sets <- function(pay) {
  ids_all <- outlier_instruments(pay); k <- length(ids_all); cap <- floor(OUTLIER_CAP * k)
  ran <- names(OUTLIER_DETECTORS)[vapply(names(OUTLIER_DETECTORS), function(m) !is.null(.o_dg(pay, m, OUTLIER_FLAG_ROW[[m]])), TRUE)]
  flags <- lapply(stats::setNames(ran, ran), function(m) intersect(.o_ids(.o_dg(pay, m, OUTLIER_FLAG_ROW[[m]])$note), ids_all))
  votes <- if (length(flags)) table(unlist(flags)) else table(character())
  rows <- lapply(ran, function(m) list(set = m, label = unname(OUTLIER_DETECTORS[[m]]), ids = flags[[m]]))
  cons <- names(votes)[votes >= 2]
  if (length(ran) >= 2 && length(cons)) rows[[length(rows) + 1L]] <- list(set = "consensus", label = "Consensus (two or more detectors)", ids = ids_all[ids_all %in% cons])
  if (!length(rows)) return(data.frame(set = character(), label = character(), n_removed = integer(), agreement = numeric(), eligible = logical(), why = character(), stringsAsFactors = FALSE))
  out <- data.frame(set = vapply(rows, `[[`, "", "set"), label = vapply(rows, `[[`, "", "label"), stringsAsFactors = FALSE)
  out$ids <- lapply(rows, `[[`, "ids"); out$n_removed <- lengths(out$ids)
  out$agreement <- vapply(out$ids, function(x) if (length(x)) mean(as.numeric(votes[x])) else NA_real_, 0)
  out$why <- ifelse(out$n_removed == 0L, "nothing flagged",
             ifelse(out$n_removed > cap, sprintf("over the cap (removes %d of %d; at most %d, 20%%)", out$n_removed, k, cap),
             ifelse(k - out$n_removed < OUTLIER_FLOOR, sprintf("fewer than %d instruments would remain", OUTLIER_FLOOR), "")))
  out$eligible <- !nzchar(out$why); out$why[out$eligible] <- "eligible"
  attr(out, "k") <- k; attr(out, "detectors_run") <- ran; attr(out, "flags") <- flags
  out
}

#' FR-96: the comparison statistics of one analysed set, read from its pair payload.
outlier_stats <- function(pay) {
  q <- .o_dg(pay, "ivw_mre", "cochran_q"); rq <- .o_dg(pay, "egger", "rucker_q"); it <- .o_est(pay, "egger", "pleiotropy_intercept"); g <- .o_dg(pay, "presso", "presso_global_test")
  e <- function(m, f) { x <- .o_est(pay, m); if (is.null(x)) NA_real_ else .o_num(x[[f]]) }
  Q <- .o_num(q$value); df <- .o_num(q$df)
  data.frame(n = e("ivw_mre", "n_ivs_used"), q = Q, q_df = df, q_p = .o_num(q$p), i2 = if (is.finite(Q) && Q > 0) max(0, (Q - df) / Q) else NA_real_,
             excess = Q - df, q_rucker = .o_num(rq$value), q_diff = Q - .o_num(rq$value),
             int = if (is.null(it)) NA_real_ else .o_num(it$estimate), int_lo = if (is.null(it)) NA_real_ else .o_num(it$ci_lower),
             int_hi = if (is.null(it)) NA_real_ else .o_num(it$ci_upper), int_p = if (is.null(it)) NA_real_ else .o_num(it$p),
             presso_p = .o_num(g$p),
             ivw = e("ivw_mre", "estimate"), ivw_lo = e("ivw_mre", "ci_lower"), ivw_hi = e("ivw_mre", "ci_upper"),
             egger = e("egger", "estimate"), egger_lo = e("egger", "ci_lower"), egger_hi = e("egger", "ci_upper"),
             wm = e("wm", "estimate"), wm_lo = e("wm", "ci_lower"), wm_hi = e("wm", "ci_upper"), stringsAsFactors = FALSE)
}

#' FR-99: the pre-specified selection rule. It never looks at the causal estimates or their P values.
#' scen: one row per scenario with set, eligible, n_removed, agreement and the outlier_stats columns; the first row is "all".
#' Step 1 keeps eligible sets after which every signal that started the module is gone (descriptively); step 2 takes the one that
#' removes the fewest variants, ties to higher detector agreement, then to larger excess reduction; step 3 otherwise takes the
#' eligible set with the largest positive excess reduction; else all instruments are kept.
outlier_select <- function(scen, triggered) {
  if (!length(triggered)) return(list(set = "all", step = 0L, reason = "No heterogeneity or pleiotropy signal (Cochran's Q, MR-Egger intercept and MR-PRESSO global test all P of 0.05 or above, or not run), so all instruments are kept and nothing is removed."))
  s0 <- scen[scen$set == "all", , drop = FALSE]
  c1 <- scen[scen$set != "all" & scen$eligible, , drop = FALSE]
  if (!nrow(c1)) return(list(set = "all", step = 3L, reason = "No candidate set is eligible (each removes nothing, more than 20% of the instruments, or leaves fewer than 5), so all instruments are kept."))
  cleared <- function(r) all(vapply(triggered, function(s) switch(s, q = isTRUE(r$q_p >= OUTLIER_TRIGGER), intercept = isTRUE(r$int_p >= OUTLIER_TRIGGER), presso = isTRUE(r$presso_p >= OUTLIER_TRIGGER), FALSE), TRUE))
  ok <- c1[vapply(seq_len(nrow(c1)), function(i) cleared(c1[i, , drop = FALSE]), TRUE), , drop = FALSE]
  sig <- paste(c(q = "Cochran's Q", intercept = "the MR-Egger intercept", presso = "the MR-PRESSO global test")[triggered], collapse = " and ")
  if (nrow(ok)) {
    o <- ok[order(ok$n_removed, -ok$agreement, -ok$excess_reduction), , drop = FALSE][1, ]
    tie <- sum(ok$n_removed == o$n_removed) > 1
    return(list(set = o$set, step = 2L, reason = sprintf("%s is the eligible set that removes the fewest variants (%d) after which %s no longer %s P below 0.05%s.",
      o$label, o$n_removed, sig, if (length(triggered) > 1) "give" else "gives", if (tie) "; equal counts were decided by detector agreement, then by excess reduction" else "")))
  }
  pos <- c1[is.finite(c1$excess_reduction) & c1$excess_reduction > 0, , drop = FALSE]
  if (nrow(pos)) { o <- pos[order(-pos$excess_reduction), , drop = FALSE][1, ]
    return(list(set = o$set, step = 3L, reason = sprintf("No eligible set removes the signal from %s. %s gives the largest excess reduction (Q falls by %s more than the %d variants removed would explain).", sig, o$label, format(round(o$excess_reduction, 1), nsmall = 1), o$n_removed))) }
  list(set = "all", step = 3L, reason = sprintf("No eligible set removes the signal from %s or lowers Q by more than the number of variants it removes, so all instruments are kept (no eligible improvement).", sig))
}

#' Builds the scenario table from the primary payload and the payloads of the analysed sets (a named list by set id; NULL if not run).
outlier_scenarios <- function(pay, sets, pays) {
  base <- cbind(data.frame(set = "all", label = "All instruments", n_removed = 0L, agreement = NA_real_, eligible = TRUE, why = "primary analysis", stringsAsFactors = FALSE), outlier_stats(pay))
  base$ids <- list(character())
  rows <- lapply(seq_len(nrow(sets)), function(i) {
    p <- pays[[sets$set[i]]]; st <- if (is.null(p)) outlier_stats(list()) else outlier_stats(p)
    r <- cbind(sets[i, c("set", "label", "n_removed", "agreement", "eligible", "why"), drop = FALSE], st); r$ids <- sets$ids[i]; r })
  sc <- do.call(rbind, c(list(base), rows))
  sc$q_change <- base$q - sc$q; sc$excess_reduction <- sc$q_change - sc$n_removed
  sc
}

#' FR-93 and FR-97: writes the flag ledger and the comparison table into the pair folder (08_outlier_analysis) and registers them.
outlier_write <- function(study_dir, pair_folder, pay, sets, scen, sel, user_set = NULL) {
  rel <- file.path("03_pairs", pair_folder, "08_outlier_analysis"); dir.create(file.path(study_dir, rel), recursive = TRUE, showWarnings = FALSE)
  ids <- outlier_instruments(pay); flags <- attr(sets, "flags") %||% list(); chosen <- user_set %||% sel$set
  rm_ids <- if (identical(chosen, "all")) character() else unlist(scen$ids[scen$set == chosen])
  vv <- Filter(function(v) isTRUE(v$keep), pay$variants)
  led <- data.frame(variant_id = ids, beta_x = vapply(vv, function(v) .o_num(v$bx), 0), se_x = vapply(vv, function(v) .o_num(v$sx), 0), stringsAsFactors = FALSE)
  led$f_stat <- (led$beta_x / led$se_x)^2
  for (m in names(flags)) led[[paste0("flagged_", m)]] <- ifelse(led$variant_id %in% flags[[m]], "yes", "no")
  led$n_detectors <- if (length(flags)) rowSums(vapply(flags, function(f) led$variant_id %in% f, logical(nrow(led))) * 1) else 0
  led$removed_in_selected_set <- ifelse(led$variant_id %in% rm_ids, "yes", "no")
  led$not_removed_because <- ifelse(led$variant_id %in% rm_ids | led$n_detectors == 0, "", if (identical(chosen, "all")) "all instruments kept" else "not in the selected set")
  utils::write.csv(led, file.path(study_dir, rel, "outliers_flagged.csv"), row.names = FALSE, na = "")
  cmp <- scen; cmp$removed_ids <- vapply(cmp$ids, function(x) paste(x, collapse = ";"), ""); cmp$ids <- NULL
  cmp$decision <- ifelse(cmp$set == chosen, if (is.null(user_set) || identical(user_set, sel$set)) "selected by the rule" else "chosen by the analyst", ifelse(cmp$set == sel$set, "selected by the rule, not adopted", ""))
  cmp$note <- "Q and MR-Egger intercept P values after removal are descriptive only: choosing which variants to remove makes them optimistic. Causal estimates are shown, not used for the selection."
  utils::write.csv(cmp, file.path(study_dir, rel, "outlier_comparison.csv"), row.names = FALSE, na = "")
  for (f in c("outliers_flagged.csv", "outlier_comparison.csv")) manifest_register(study_dir, file.path(rel, f), role = "outlier_analysis")
  invisible(rel)
}

#' The record kept in the study state and sent to the screen: signals, scenarios (one list per row), the rule's choice and its reason.
#' review is "awaiting review" until the analyst confirms or changes the choice (FR-102).
outlier_record <- function(parent, sig, trig, scen, sel, chosen) {
  cols <- setdiff(names(scen), "ids"); fin <- function(v) if (is.numeric(v) && !is.finite(v)) NULL else v
  rows <- lapply(seq_len(nrow(scen)), function(i) c(lapply(stats::setNames(cols, cols), function(k) fin(scen[[k]][i])), list(ids = as.list(scen$ids[[i]]))))
  list(parent = parent, signals = lapply(sig, function(x) if (is.na(x)) NULL else x), triggered = as.list(trig), scenarios = rows,
       selected = sel$set, step = sel$step, reason = sel$reason, chosen = chosen %||% sel$set, review = "awaiting review",
       cap = OUTLIER_CAP, floor = OUTLIER_FLOOR, trigger = OUTLIER_TRIGGER)
}

#' FR-102: records the analyst's choice in the two files written by outlier_write.
outlier_write_choice <- function(study_dir, rec) {
  rel <- file.path("03_pairs", rec$parent, "08_outlier_analysis"); fc <- file.path(study_dir, rel, "outlier_comparison.csv"); fl <- file.path(study_dir, rel, "outliers_flagged.csv")
  if (!file.exists(fc) || !file.exists(fl)) return(invisible(NULL))
  cmp <- utils::read.csv(fc, stringsAsFactors = FALSE, colClasses = c(removed_ids = "character")); led <- utils::read.csv(fl, stringsAsFactors = FALSE)
  own <- !identical(rec$chosen, rec$selected)
  cmp$decision <- ifelse(cmp$set == rec$chosen, if (own) "chosen by the analyst" else "selected by the rule, confirmed by the analyst", ifelse(cmp$set == rec$selected, "selected by the rule, not adopted", ""))
  rm_ids <- if (identical(rec$chosen, "all")) character() else strsplit(cmp$removed_ids[cmp$set == rec$chosen][1] %||% "", ";", fixed = TRUE)[[1]]
  led$removed_in_selected_set <- ifelse(led$variant_id %in% rm_ids, "yes", "no")
  led$not_removed_because <- ifelse(led$variant_id %in% rm_ids | led$n_detectors == 0, "", if (identical(rec$chosen, "all")) "all instruments kept" else "not in the chosen set")
  utils::write.csv(cmp, fc, row.names = FALSE, na = ""); utils::write.csv(led, fl, row.names = FALSE, na = "")
  for (f in c("outliers_flagged.csv", "outlier_comparison.csv")) manifest_register(study_dir, file.path(rel, f), role = "outlier_analysis")
  invisible(rel)
}

#' Report text (FR-103). Methods: the module as pre-specified; Results: per pair, the rule's choice, the adopted set and the IVW estimates.
OUTLIER_METHODS_TEXT <- paste(
  "Where Cochran's Q, the MR-Egger intercept or the MR-PRESSO global test gave P below 0.05, an outlier analysis was run as a data-driven sensitivity analysis.",
  "Each outlier detector that was run (MR-PRESSO, radial IVW and radial MR-Egger with each variant tested at 0.05 divided by the number of variants unless stated otherwise, and MR-Lasso) proposed a set of variants to remove, and variants flagged by at least two detectors formed a consensus set.",
  "A set was eligible if it removed at most 20% of the instruments and left at least 5. Each eligible set was analysed again with IVW, MR-Egger and the weighted median (and MR-PRESSO when its global test had started the analysis).",
  "A rule fixed before the results were seen chose among them without using the causal estimates: the eligible set that removed the fewest variants after which every starting signal had P of 0.05 or above (ties broken by agreement between detectors, then by the fall in Q beyond the number of variants removed); otherwise the set with the largest such fall; otherwise all instruments were kept.",
  "The analyst reviewed the choice before it was adopted. The primary analysis used all instruments [Burgess 2023].")

outlier_results_text <- function(ol, ce, sc) {
  lab <- function(set) { r <- Filter(function(x) identical(x$set, set), ol$scenarios); if (length(r)) r[[1]]$label else "All instruments" }
  trig <- unlist(ol$triggered); sig <- if (length(trig)) paste(c(q = "Cochran's Q", intercept = "the MR-Egger intercept", presso = "the MR-PRESSO global test")[trig], collapse = " and ") else "no signal"
  ne <- sum(vapply(ol$scenarios, function(x) isTRUE(x$eligible) && !identical(x$set, "all"), TRUE))
  ch <- ol$chosen_row; base <- ce[ce$method == "ivw_mre", , drop = FALSE]
  what <- if (is.null(ch) || identical(ol$chosen, "all")) "all instruments were kept" else sprintf("the %s set was adopted (%d variant%s removed: %s; excess reduction in Q %s; Q P = %s after removal, descriptive), with an IVW estimate of %s%s",
    ch$label, as.integer(ch$n_removed), if (ch$n_removed == 1) "" else "s", paste(unlist(ch$ids), collapse = ", "), format(round(ch$excess_reduction %||% NA_real_, 1), nsmall = 1), .fp(ch$q_p %||% NA_real_),
    if (NROW(ol$chosen_ivw)) .ci_phrase(ol$chosen_ivw[1, ], ol$chosen_ivw$scale[1]) else if (!is.null(ch$ivw)) sprintf("%s (%s to %s)", .fnum(ch$ivw, 3), .fnum(ch$ivw_lo, 3), .fnum(ch$ivw_hi, 3)) else "not available",
    if (nrow(base)) sprintf(", compared with %s in the primary analysis", .ci_phrase(base[1, ], sc)) else "")
  rule <- if (identical(ol$selected, "all")) "to keep all instruments" else sprintf("the %s set", lab(ol$selected))
  sprintf("Outlier analysis (started by %s): %d eligible set%s analysed. The pre-specified rule chose %s. %s%s.%s",
    sig, ne, if (ne == 1) " was" else "s were", rule, toupper(substr(what, 1, 1)), substring(what, 2),
    if (!identical(ol$chosen, ol$selected)) sprintf(" The analyst adopted %s instead of the rule's choice.", lab(ol$chosen)) else if (identical(ol$review, "awaiting review")) " The choice was not yet reviewed by the analyst." else "")
}

