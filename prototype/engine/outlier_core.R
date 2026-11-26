# Round 42: outlier module (requirements_spec.md, Part H, FR-90 to FR-103).
# The detectors are the packaged methods (MR-PRESSO, radial IVW, radial MR-Egger, MR-Lasso) run as part of the pair; this file only
# reads what they flagged, builds the candidate sets, applies the safeguards and the pre-specified selection rule, and writes the
# flag ledger and the comparison table. Every re-analysis is an ordinary pair run with its own scripts (studio_run_without), so no
# statistical step lives here beyond arithmetic on the packages' own results (Q minus df, I squared, excess reduction).
# The primary analysis is never changed: the selected set drives a labelled, data-driven sensitivity analysis (FR-100).

OUTLIER_DETECTORS <- c(presso = "MR-PRESSO", radial = "Radial IVW", radegger = "Radial MR-Egger", lasso = "MR-Lasso", cook = "Cook's distance", studres = "Studentized residuals")   # Round 109: cook, studres
OUTLIER_FLAG_ROW <- c(presso = "presso_outliers", radial = "radial_outliers", radegger = "radial_outliers", lasso = "lasso_excluded", cook = "cook_outliers", studres = "studres_outliers")
OUTLIER_CAP <- 0.20      # FR-95a; kept by the Round 32 calibration
OUTLIER_FLOOR <- 5L      # FR-95b
OUTLIER_TRIGGER <- 0.05  # FR-90; kept by the Round 32 calibration
OUTLIER_REANALYSIS <- c("ivw_mre", "egger", "wm")   # FR-96; MR-PRESSO is added when its global test triggered the module
# Round 109: the outlier detectors and the methods that give no causal estimate are not rerun on the reduced sets; every other method of the study is.
OUTLIER_NOT_RERUN <- c("presso", "radial", "radegger", "lasso", "cook", "studres", "mrclust", "steiger", "loo", "sign", "clr", "rucker", "ivw_steiger")
outlier_methods <- function(pay) { m <- unique(vapply(Filter(function(e) identical(e$estimand, "causal_effect"), pay$estimates %||% list()), function(e) as.character(e$method), ""))
  union(OUTLIER_REANALYSIS, setdiff(m, OUTLIER_NOT_RERUN)) }

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
outlier_stats <- function(pay, methods = character()) {   # Round 109: methods: every other estimator of the study, for the comparison table
  q <- .o_dg(pay, "ivw_mre", "cochran_q"); rq <- .o_dg(pay, "egger", "rucker_q"); it <- .o_est(pay, "egger", "pleiotropy_intercept"); g <- .o_dg(pay, "presso", "presso_global_test")
  e <- function(m, f) { x <- .o_est(pay, m); if (is.null(x)) NA_real_ else .o_num(x[[f]]) }
  Q <- .o_num(q$value); df <- .o_num(q$df)
  out <- data.frame(n = e("ivw_mre", "n_ivs_used"), q = Q, q_df = df, q_p = .o_num(q$p), i2 = if (is.finite(Q) && Q > 0) max(0, (Q - df) / Q) else NA_real_,
             excess = Q - df, q_rucker = .o_num(rq$value), q_diff = Q - .o_num(rq$value),
             int = if (is.null(it)) NA_real_ else .o_num(it$estimate), int_lo = if (is.null(it)) NA_real_ else .o_num(it$ci_lower),
             int_hi = if (is.null(it)) NA_real_ else .o_num(it$ci_upper), int_p = if (is.null(it)) NA_real_ else .o_num(it$p),
             presso_p = .o_num(g$p),
             ivw = e("ivw_mre", "estimate"), ivw_lo = e("ivw_mre", "ci_lower"), ivw_hi = e("ivw_mre", "ci_upper"),
             egger = e("egger", "estimate"), egger_lo = e("egger", "ci_lower"), egger_hi = e("egger", "ci_upper"),
             wm = e("wm", "estimate"), wm_lo = e("wm", "ci_lower"), wm_hi = e("wm", "ci_upper"), stringsAsFactors = FALSE)
  for (m in setdiff(methods, c("ivw_mre", "egger", "wm"))) { out[[paste0("m_", m)]] <- e(m, "estimate"); out[[paste0("m_", m, "_lo")]] <- e(m, "ci_lower"); out[[paste0("m_", m, "_hi")]] <- e(m, "ci_upper") }
  out
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
  ms <- outlier_methods(pay)
  base <- cbind(data.frame(set = "all", label = "All instruments", n_removed = 0L, agreement = NA_real_, eligible = TRUE, why = "primary analysis", stringsAsFactors = FALSE), outlier_stats(pay, ms))
  base$ids <- list(character())
  rows <- lapply(seq_len(nrow(sets)), function(i) {
    p <- pays[[sets$set[i]]]; st <- if (is.null(p)) outlier_stats(list(), ms) else outlier_stats(p, ms)
    r <- cbind(sets[i, c("set", "label", "n_removed", "agreement", "eligible", "why"), drop = FALSE], st); r$ids <- sets$ids[i]; r })
  sc <- do.call(rbind, c(list(base), rows)); attr(sc, "methods") <- ms
  sc$q_change <- base$q - sc$q; sc$excess_reduction <- sc$q_change - sc$n_removed
  outlier_strength(sc, pay)
}

#' Round 48 (FR-95c): instrument strength of every set, recomputed from the remaining variants: mean F ((beta_x / se_x)^2) and I2GX
#' (Bowden 2016, TwoSampleMR::Isq on |beta_x|). A set below the R23 (mean F) or R24 (I2GX) threshold is flagged amber; it stays eligible.
outlier_strength <- function(sc, pay) {
  v <- Filter(function(x) isTRUE(x$keep), pay$variants); id <- vapply(v, function(x) as.character(x$id), ""); bx <- vapply(v, function(x) .o_num(x$bx), 0); sx <- vapply(v, function(x) .o_num(x$sx), 0)
  thr <- function(r, d) { x <- tryCatch(suppressWarnings(as.numeric(feas_rules()$threshold[feas_rules()$rule_id == r])), error = function(e) NA_real_); if (length(x) && is.finite(x[1])) x[1] else d }
  tf <- thr("R23", 10); ti <- thr("R24", 0.9)
  st <- lapply(seq_len(nrow(sc)), function(i) { k <- !(id %in% unlist(sc$ids[[i]])); b <- bx[k]; s <- sx[k]
    c(mean_f = mean((b / s)^2), i2gx = if (sum(k) >= 2) TwoSampleMR::Isq(abs(b), s) else NA_real_) })
  sc$mean_f <- vapply(st, `[[`, 0, "mean_f"); sc$i2gx <- vapply(st, `[[`, 0, "i2gx")
  sc$strength_flag <- ifelse(is.finite(sc$mean_f) & sc$mean_f < tf, sprintf("mean F below %s", format(tf)), ifelse(is.finite(sc$i2gx) & sc$i2gx < ti, sprintf("I2GX below %s", format(ti)), ""))
  sc
}

#' FR-93 and FR-97: writes the flag ledger and the comparison table into the pair folder (08_outlier_analysis) and registers them.
outlier_write <- function(study_dir, pair_folder, pay, sets, scen, sel, user_set = NULL, folders = NULL) {
  rel <- file.path("03_pairs", pair_folder, "08_outlier_analysis"); for (d in c("01_flagged", "02_scenarios", "03_tables")) dir.create(file.path(study_dir, rel, d), recursive = TRUE, showWarnings = FALSE)   # FR-101
  ids <- outlier_instruments(pay); flags <- attr(sets, "flags") %||% list(); chosen <- user_set %||% sel$set
  rm_ids <- if (identical(chosen, "all")) character() else unlist(scen$ids[scen$set == chosen])
  vv <- Filter(function(v) isTRUE(v$keep), pay$variants)
  led <- data.frame(variant_id = ids, beta_x = vapply(vv, function(v) .o_num(v$bx), 0), se_x = vapply(vv, function(v) .o_num(v$sx), 0), stringsAsFactors = FALSE)
  led$f_stat <- (led$beta_x / led$se_x)^2
  for (m in names(flags)) led[[paste0("flagged_", m)]] <- ifelse(led$variant_id %in% flags[[m]], "yes", "no")
  led$n_detectors <- if (length(flags)) rowSums(vapply(flags, function(f) led$variant_id %in% f, logical(nrow(led))) * 1) else 0
  led$removed_in_selected_set <- ifelse(led$variant_id %in% rm_ids, "yes", "no")
  led$not_removed_because <- ifelse(led$variant_id %in% rm_ids | led$n_detectors == 0, "", if (identical(chosen, "all")) "all instruments kept" else "not in the selected set")
  utils::write.csv(led, file.path(study_dir, rel, "01_flagged", "outliers_flagged.csv"), row.names = FALSE, na = "")
  cmp <- scen; cmp$removed_ids <- vapply(cmp$ids, function(x) paste(x, collapse = ";"), ""); cmp$ids <- NULL
  cmp$decision <- ifelse(cmp$set == chosen, if (is.null(user_set) || identical(user_set, sel$set)) "selected by the rule" else "chosen by the analyst", ifelse(cmp$set == sel$set, "selected by the rule, not adopted", ""))
  cmp$note <- "Q and MR-Egger intercept P values after removal are descriptive only: choosing which variants to remove makes them optimistic. Causal estimates are shown, not used for the selection."
  utils::write.csv(cmp, file.path(study_dir, rel, "03_tables", "outlier_comparison.csv"), row.names = FALSE, na = "")
  if (exists("xlsx_write")) tryCatch(xlsx_write(file.path(study_dir, rel, "03_tables", "outlier_comparison.xlsx"), list(Comparison = cmp, `Flag ledger` = led)), error = function(e) NULL)
  sc <- data.frame(set = scen$set, label = scen$label, n_removed = scen$n_removed, eligible = scen$eligible, removed_ids = vapply(scen$ids, function(x) paste(x, collapse = ";"), ""),
                   pair_folder = vapply(scen$set, function(k) if (identical(k, "all")) pair_folder else as.character(folders[[k]] %||% ""), ""), stringsAsFactors = FALSE)
  utils::write.csv(sc, file.path(study_dir, rel, "02_scenarios", "scenarios.csv"), row.names = FALSE, na = "")
  writeLines(c("Outlier analysis of this pair (MR Studio; requirements FR-90 to FR-106).", "",
               "01_flagged/   flag ledger: one row per instrument with each detector's flag, strength and the removal decision; nearest genes and other-trait associations when looked up",
               "02_scenarios/ scenarios.csv: each candidate set and the pair folder (in 03_pairs) that holds its analysis data, scripts, console log and replay check",
               "03_tables/    comparison of all sets (CSV and Excel); P values after removal are descriptive only",
               "", "Figures (comparison, detector agreement, radial plot, leave-one-out before and after) are written into the publication package."),
             file.path(study_dir, rel, "README.txt"))
  outlier_register(study_dir, rel)
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
  rel <- file.path("03_pairs", rec$parent, "08_outlier_analysis"); fc <- file.path(study_dir, rel, "03_tables", "outlier_comparison.csv"); fl <- file.path(study_dir, rel, "01_flagged", "outliers_flagged.csv")
  if (!file.exists(fc) || !file.exists(fl)) return(invisible(NULL))
  cmp <- utils::read.csv(fc, stringsAsFactors = FALSE, colClasses = c(removed_ids = "character")); led <- utils::read.csv(fl, stringsAsFactors = FALSE)
  own <- !identical(rec$chosen, rec$selected)
  cmp$decision <- ifelse(cmp$set == rec$chosen, if (own) "chosen by the analyst" else "selected by the rule, confirmed by the analyst", ifelse(cmp$set == rec$selected, "selected by the rule, not adopted", ""))
  rm_ids <- if (identical(rec$chosen, "all")) character() else strsplit(cmp$removed_ids[cmp$set == rec$chosen][1] %||% "", ";", fixed = TRUE)[[1]]
  led$removed_in_selected_set <- ifelse(led$variant_id %in% rm_ids, "yes", "no")
  led$not_removed_because <- ifelse(led$variant_id %in% rm_ids | led$n_detectors == 0, "", if (identical(rec$chosen, "all")) "all instruments kept" else "not in the chosen set")
  utils::write.csv(cmp, fc, row.names = FALSE, na = ""); utils::write.csv(led, fl, row.names = FALSE, na = "")
  if (exists("xlsx_write")) tryCatch(xlsx_write(file.path(study_dir, rel, "03_tables", "outlier_comparison.xlsx"), list(Comparison = cmp, `Flag ledger` = led)), error = function(e) NULL)
  outlier_register(study_dir, rel)
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
  sprintf("Outlier analysis (started by %s): %d eligible set%s analysed. The pre-specified rule chose %s. %s%s.%s%s%s",
    sig, ne, if (ne == 1) " was" else "s were", rule, toupper(substr(what, 1, 1)), substring(what, 2),
    if (!is.null(ch) && !identical(ol$chosen, "all") && nzchar(ch$strength_flag %||% "")) sprintf(" After removal the instruments were weaker (mean F %s, I2GX %s: %s).", format(round(ch$mean_f, 1), nsmall = 1), format(round(ch$i2gx %||% NA_real_, 2), nsmall = 2), ch$strength_flag) else "",
    outlier_lookup_text(ol),
    if (!identical(ol$chosen, ol$selected)) (if (identical(ol$chosen, "all")) " The analyst kept all instruments instead of the rule's choice." else sprintf(" The analyst adopted the %s set instead of the rule's choice.", lab(ol$chosen))) else if (identical(ol$review, "awaiting review")) " The choice was not yet reviewed by the analyst." else "")
}


# ---- Round 44: comparison figures (FR-98), drawn from the saved record. Base graphics, written by the package's .fig_write. ----
.outl_rows <- function(ol) Filter(function(x) identical(x$set, "all") || isTRUE(x$eligible), ol$scenarios)
.outl_num <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))

#' Rows: all instruments and every analysed set. Panels: IVW estimate with CI; I squared; MR-Egger intercept with CI and null line;
#' MR-PRESSO global P when any set has it. The adopted set is drawn in the accent colour, the rule's choice with an open marker.
draw_outlier_comparison <- function(ol, unit) {
  r <- .outl_rows(ol); k <- length(r); lab <- vapply(r, function(x) sprintf("%s (%d out)", if (identical(x$set, "consensus")) "Consensus" else x$label, as.integer(x$n_removed)), "")
  g <- function(f) vapply(r, function(x) .outl_num(x[[f]]), 0); pr <- g("presso_p"); np <- if (any(is.finite(pr))) 4L else 3L
  col <- ifelse(vapply(r, function(x) identical(x$set, ol$chosen), TRUE), "#2b3fc7", "grey25"); rule <- vapply(r, function(x) identical(x$set, ol$selected), TRUE)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op)); graphics::layout(matrix(seq_len(np), 1), widths = c(1.9, rep(1, np - 1L)))
  panel <- function(b, lo, hi, xlab, ref = NA, ylabs = FALSE, xlim = NULL) {
    graphics::par(mar = c(3.4, if (ylabs) 12.5 else 0.6, 1.6, 0.6), mgp = c(1.9, 0.6, 0), cex = 0.7, las = 1)
    xl <- xlim %||% { v <- c(b, lo, hi, ref); v <- v[is.finite(v)]; if (length(v)) range(v) + c(-1, 1) * max(diff(range(v)) * 0.08, 1e-6) else c(0, 1) }
    graphics::plot(NA, xlim = xl, ylim = c(k + 0.5, 0.5), yaxt = "n", xlab = xlab, ylab = "", bty = "n")
    if (ylabs) graphics::axis(2, at = seq_len(k), labels = lab, tick = FALSE, cex.axis = 0.9)
    if (is.finite(ref)) graphics::abline(v = ref, lty = 2, col = "grey55")
    ok <- is.finite(lo) & is.finite(hi); graphics::segments(lo[ok], which(ok), hi[ok], which(ok), col = col[ok], lwd = 1.3)
    graphics::points(b, seq_len(k), pch = ifelse(rule, 1, 16), cex = 1.1, col = col)
  }
  panel(g("ivw"), g("ivw_lo"), g("ivw_hi"), paste("IVW,", unit), ref = 0, ylabs = TRUE); graphics::title("Causal estimate (shown, not used)", cex.main = 0.85, font.main = 1)
  panel(100 * g("i2"), NA, NA, "I squared (%)", xlim = c(0, 100)); graphics::title("Heterogeneity", cex.main = 0.85, font.main = 1)
  panel(g("int"), g("int_lo"), g("int_hi"), "MR-Egger intercept", ref = 0); graphics::title("Directional pleiotropy", cex.main = 0.85, font.main = 1)
  if (np == 4L) { panel(-log10(pr), NA, NA, "-log10 P", ref = -log10(0.05)); graphics::title("MR-PRESSO global test", cex.main = 0.85, font.main = 1) }
}

#' Which detector flagged which variant (variants flagged by at least one detector; at most 40, most agreed first), with the variants
#' removed in the adopted set marked.
OUTLIER_AGREE_MAX <- 40L
outlier_agreement <- function(ol) {
  det <- Filter(function(x) x$set %in% names(OUTLIER_DETECTORS), ol$scenarios); if (!length(det)) return(NULL)
  ids <- unique(unlist(lapply(det, function(x) unlist(x$ids)))); if (!length(ids)) return(NULL)
  m <- vapply(det, function(x) ids %in% unlist(x$ids), logical(length(ids))); m <- matrix(m, nrow = length(ids), dimnames = list(ids, vapply(det, `[[`, "", "label")))
  ch <- Filter(function(x) identical(x$set, ol$chosen), ol$scenarios); rm <- if (length(ch)) unlist(ch[[1]]$ids) else character()
  o <- order(-rowSums(m), rownames(m)); m <- m[o, , drop = FALSE]; shown <- min(nrow(m), OUTLIER_AGREE_MAX)
  list(m = m[seq_len(shown), , drop = FALSE], removed = rownames(m)[seq_len(shown)] %in% rm, total = nrow(m))
}
draw_outlier_agreement <- function(ag) {
  m <- ag$m; k <- nrow(m); d <- ncol(m)
  op <- graphics::par(mar = c(1, 8.6, 9, 1), cex = 0.72, las = 1); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0.5, d + 1.5), ylim = c(k + 0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  graphics::axis(2, at = seq_len(k), labels = rownames(m), tick = FALSE, cex.axis = 0.85); graphics::axis(3, at = c(seq_len(d), d + 1), labels = c(colnames(m), "Removed"), tick = FALSE, las = 2, cex.axis = 0.85)
  for (j in seq_len(d)) graphics::points(rep(j, k), seq_len(k), pch = ifelse(m[, j], 15, 0), cex = 1.3, col = ifelse(m[, j], "grey20", "grey75"))
  graphics::points(rep(d + 1, k), seq_len(k), pch = ifelse(ag$removed, 16, 1), cex = 1.2, col = ifelse(ag$removed, "#2b3fc7", "grey75"))
}

#' Round 47: radial plot (Bowden et al. 2018): x = |beta_x| / se_y, y = the variant's ratio estimate times x, so the IVW estimate (first-order
#' weights) is the slope of a line through the origin. Variants flagged by any detector are filled orange; those removed in the adopted set
#' are ringed. This is a drawing of the saved data, not a new estimate.
draw_outlier_radial <- function(pay, ol, unit) {   # Round 96: ol may be NULL (radial MR without an outlier analysis)
  v <- Filter(function(x) isTRUE(x$keep), pay$variants); bx <- vapply(v, function(x) .outl_num(x$bx), 0); by <- vapply(v, function(x) .outl_num(x$by), 0); sy <- vapply(v, function(x) .outl_num(x$sy), 0)
  id <- vapply(v, function(x) as.character(x$id), ""); xr <- abs(bx) / sy; yr <- (by / bx) * xr
  fl <- if (is.null(ol)) character() else unique(unlist(lapply(Filter(function(x) x$set %in% names(OUTLIER_DETECTORS), ol$scenarios), function(x) unlist(x$ids))))
  ch <- if (is.null(ol)) list() else Filter(function(x) identical(x$set, ol$chosen), ol$scenarios); rm <- if (length(ch)) unlist(ch[[1]]$ids) else character()
  iv <- Filter(function(e) identical(e$method, "ivw_mre") && identical(e$estimand, "causal_effect"), pay$estimates); b <- if (length(iv)) .outl_num(iv[[1]]$estimate) else sum(bx * by / sy^2) / sum(bx^2 / sy^2)
  op <- graphics::par(mar = c(4.4, 4.9, 0.8, 0.8), mgp = c(3.2, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  graphics::plot(xr, yr, type = "n", xlim = c(0, max(xr) * 1.05), xlab = "Square root of weight (|beta exposure| / SE outcome)", ylab = "Ratio estimate x square root of weight", bty = "l")
  graphics::abline(h = 0, col = "grey75"); graphics::abline(a = 0, b = b, col = "#2b3fc7", lwd = 2)
  f <- id %in% fl; graphics::points(xr[!f], yr[!f], pch = 16, cex = 0.7, col = "grey35"); graphics::points(xr[f], yr[f], pch = 16, cex = 0.95, col = "#d4874f")
  r <- id %in% rm; if (any(r)) graphics::points(xr[r], yr[r], pch = 1, cex = 1.7, lwd = 1.3, col = "#2b3fc7")
  if (is.null(ol)) graphics::legend("bottomright", legend = sprintf("IVW slope (%s)", format(round(b, 3), nsmall = 3)), col = "#2b3fc7", lty = 1, lwd = 2, bty = "n", cex = 0.85) else
  graphics::legend("bottomright", legend = c(sprintf("IVW slope (%s)", format(round(b, 3), nsmall = 3)), "flagged by a detector", "removed in the adopted set"), col = c("#2b3fc7", "#d4874f", "#2b3fc7"),
                   lty = c(1, NA, NA), pch = c(NA, 16, 1), lwd = c(2, NA, 1.3), bty = "n", cex = 0.85)
}

#' Round 47: leave-one-out before and after removal. Each row is one IVW estimate with one variant left out (from the leave-one-out task of
#' the primary pair and of the adopted set); vertical lines are the estimates with all variants of each set.
draw_outlier_loo <- function(lb, la, fb, fa, unit) {
  op <- graphics::par(mar = c(3.4, 9, 0.8, 1), mgp = c(1.9, 0.6, 0), cex = 0.75, las = 1); on.exit(graphics::par(op))
  xl <- range(c(lb$loo_b, la$loo_b, fb, fa), na.rm = TRUE); xl <- xl + c(-1, 1) * diff(xl) * 0.08
  graphics::plot(NA, xlim = xl, ylim = c(2.6, 0.4), yaxt = "n", xlab = unit, ylab = "", bty = "n")
  graphics::axis(2, at = 1:2, labels = c(sprintf("All instruments (%d)", nrow(lb)), sprintf("Adopted set (%d)", nrow(la))), tick = FALSE)
  set.seed(1); jb <- 1 + stats::runif(nrow(lb), -0.22, 0.22); ja <- 2 + stats::runif(nrow(la), -0.22, 0.22)   # fixed jitter, for display only
  graphics::points(lb$loo_b, jb, pch = 16, cex = 0.6, col = "grey35"); graphics::points(la$loo_b, ja, pch = 16, cex = 0.6, col = "grey35")
  graphics::segments(fb, 0.65, fb, 1.35, col = "#2b3fc7", lwd = 2); graphics::segments(fa, 1.65, fa, 2.35, col = "#2b3fc7", lwd = 2)
}

#' The legends (FR-98), numbered from s, named by figure: comparison, agreement, radial, loo (only those drawn).
outlier_legends <- function(ol, cx, s, ag, loo = FALSE) {
  out <- c(comparison = sprintf("Supplementary Figure S%%d. Outlier analysis of %s on %s: one row per analysed set of instruments (all instruments and every set that removed at most 20%%%% of them and left at least 5). Panels show the IVW estimate with its confidence interval, which the selection rule did not use; I squared; and the MR-Egger intercept with its confidence interval%s. The adopted set is in blue; an open marker shows the rule's choice. Heterogeneity and intercept values after removal are descriptive only.",
            cx$p$exName, cx$p$outName, if (any(vapply(ol$scenarios, function(x) is.finite(.outl_num(x$presso_p)), TRUE))) ", and the MR-PRESSO global test (-log10 P; dashed line at P = 0.05)" else ""))
  if (!is.null(ag)) out["agreement"] <- sprintf("Supplementary Figure S%%d. Variants flagged by each outlier detector for %s on %s (filled square: flagged). The last column marks the variants removed in the adopted set.%s", cx$p$exName, cx$p$outName,
            if (ag$total > nrow(ag$m)) sprintf(" The %d variants flagged most often are shown; all %d are in outliers_flagged.csv.", nrow(ag$m), ag$total) else "")
  out["radial"] <- sprintf("Supplementary Figure S%%d. Radial plot of %s on %s [Bowden 2018]: each point is one variant, with the square root of its first-order weight on the x axis and its ratio estimate multiplied by that on the y axis, so the IVW estimate is the slope of the line through the origin. Orange points were flagged by at least one outlier detector; ringed points were removed in the adopted set.", cx$p$exName, cx$p$outName)
  if (loo) out["loo"] <- sprintf("Supplementary Figure S%%d. Leave-one-out analysis of %s on %s before and after outlier removal: each point is the IVW estimate with one variant left out, for all instruments and for the adopted set; vertical lines are the estimates using every variant of each set. Points are spread vertically for legibility only.", cx$p$exName, cx$p$outName)
  stats::setNames(sprintf(out, s + seq_along(out) - 1L), names(out))
}

#' Round 47: whether both the primary pair and the adopted set have a leave-one-out analysis (for the before-and-after figure).
outlier_has_loo <- function(cx) { ol <- cx$outl; !is.null(ol) && !identical(ol$chosen, "all") && !is.null(ol$chosen_pay) && !is.null(loo_facts(cx)) && !is.null(loo_facts(list(pay = ol$chosen_pay))) }

# ---- Round 49: nearest gene and other-trait associations of the flagged variants (FR-93, FR-106, FR-112) -------------------------
# Only the rsIDs of flagged variants are sent: to Ensembl VEP (GRCh37) for the nearest gene and to OpenGWAS (ieugwasr::phewas, with the
# analyst's token) for associations at P < 5e-8. The app does not decide whether a variant is horizontally or vertically pleiotropic:
# it lists the associations so that the analyst can judge (FR-106).
OUTLIER_PHEWAS_P <- 5e-8
outlier_flagged_ids <- function(ol) unique(unlist(lapply(Filter(function(x) x$set %in% names(OUTLIER_DETECTORS), ol$scenarios), function(x) unlist(x$ids))))
outlier_lookup <- function(ol, jwt, backend = ogw_backend(), exclude = character()) {
  ids <- outlier_flagged_ids(ol); rs <- ids[grepl("^rs[0-9]+$", ids)]
  if (!length(ids)) stop("No variant was flagged, so there is nothing to look up.")
  if (!length(rs)) stop("The flagged variants have no rsIDs (for example simulated traits), so genes and other-trait associations cannot be looked up.")
  out <- data.frame(variant_id = rs, nearest_gene = NA_character_, chr = NA_character_, pos = NA_real_, n_other_traits = NA_integer_, top_other_traits = NA_character_, stringsAsFactors = FALSE)
  notes <- character()
  g <- tryCatch(backend$genes(rs), error = function(e) { notes <<- c(notes, paste("Nearest genes could not be looked up (Ensembl):", substr(conditionMessage(e), 1, 160))); NULL })
  if (!is.null(g) && nrow(g)) { m <- match(out$variant_id, g$rsid); out$nearest_gene <- g$gene[m]; out$chr <- g$chr[m]; out$pos <- g$pos[m] }
  if (is.null(jwt)) notes <- c(notes, "Other-trait associations need an OpenGWAS token (Settings); they were not looked up.") else {
    ph <- tryCatch(backend$phewas(rs, OUTLIER_PHEWAS_P, jwt), error = function(e) { notes <<- c(notes, ogw_friendly(e)); NULL })
    if (!is.null(ph)) {
      ph <- ph[is.finite(suppressWarnings(as.numeric(ph$p))) & as.numeric(ph$p) < OUTLIER_PHEWAS_P & !(ph$id %in% exclude), , drop = FALSE]
      for (i in seq_len(nrow(out))) { x <- ph[ph$rsid == out$variant_id[i], , drop = FALSE]; x <- x[order(as.numeric(x$p)), , drop = FALSE]; x <- x[!duplicated(x$trait), , drop = FALSE]
        out$n_other_traits[i] <- nrow(x); out$top_other_traits[i] <- if (nrow(x)) paste(sprintf("%s (P = %s)", x$trait[seq_len(min(5, nrow(x)))], formatC(as.numeric(x$p[seq_len(min(5, nrow(x)))]), digits = 2, format = "e")), collapse = "; ") else "" } } }
  list(rows = out, notes = notes, date = format(Sys.Date()), skipped = setdiff(ids, rs), sources = c(genes = "Ensembl VEP REST, GRCh37", traits = sprintf("OpenGWAS PheWAS (ieugwasr::phewas), P < %s", format(OUTLIER_PHEWAS_P))))
}

#' Registers every file of the module folder in the study manifest (so the integrity check covers them).
outlier_register <- function(study_dir, rel) { f <- list.files(file.path(study_dir, rel), recursive = TRUE); if (length(f)) manifest_register(study_dir, file.path(rel, f), role = "outlier_analysis"); invisible(f) }

#' Round 49: adds the lookup to the flag ledger (nearest gene, number and top other-trait associations) and saves it on its own.
outlier_write_lookup <- function(study_dir, parent, lk) {
  rel <- file.path("03_pairs", parent, "08_outlier_analysis"); fl <- file.path(study_dir, rel, "01_flagged", "outliers_flagged.csv"); if (!file.exists(fl)) return(invisible(NULL))
  led <- utils::read.csv(fl, stringsAsFactors = FALSE); m <- match(led$variant_id, lk$rows$variant_id)
  led$nearest_gene <- lk$rows$nearest_gene[m]; led$n_other_traits <- lk$rows$n_other_traits[m]; led$top_other_traits <- lk$rows$top_other_traits[m]
  utils::write.csv(led, fl, row.names = FALSE, na = "")
  utils::write.csv(cbind(lk$rows, looked_up = lk$date, genes_from = lk$sources[["genes"]], traits_from = lk$sources[["traits"]]), file.path(study_dir, rel, "01_flagged", "flagged_annotations.csv"), row.names = FALSE, na = "")
  outlier_register(study_dir, rel)
}

#' Round 49: report sentences for the lookup (Methods once; Results per pair).
OUTLIER_LOOKUP_METHODS <- "For the flagged variants, the nearest gene was taken from the Ensembl Variant Effect Predictor (GRCh37) and associations with other traits at P < 5e-8 from the OpenGWAS PheWAS (ieugwasr::phewas), excluding the exposure and outcome datasets, so that the analyst could judge whether a variant was more likely horizontally or vertically pleiotropic; this judgement was not automated."
outlier_lookup_text <- function(ol) {
  L <- ol$lookup; if (is.null(L) || !length(L$rows)) return("")
  it <- vapply(L$rows, function(r) sprintf("%s (%s; %s)", r$variant_id, if (is.null(r$nearest_gene) || is.na(r$nearest_gene)) "gene not found" else r$nearest_gene,
               if (is.null(r$n_other_traits) || is.na(r$n_other_traits)) "other traits not looked up" else sprintf("%d other trait%s", as.integer(r$n_other_traits), if (identical(as.integer(r$n_other_traits), 1L)) "" else "s")), "")
  sprintf(" Flagged variants, with nearest gene and number of other traits at P < 5e-8 (looked up on %s): %s.", L$date, paste(it, collapse = ", "))
}
