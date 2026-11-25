# Live narration of long waits (Round 75d/75e). The Check screen and the Run screen show what the engine is doing, step by step, each step
# with one sentence on why. The engine reports steps with .say(); while a pair's scripts run, the runner calls live_scripts() once a
# second, which reports each method as it starts and finishes (each method writes its tables when it finishes, in script order).

#' What each method does, in one sentence, for the narration.
METHOD_WHY <- c(
  ivw_mre = "Inverse-variance weighted: combines the ratio estimates of all instruments, each weighted by its precision; valid when every instrument is valid or pleiotropy averages out.",
  egger = "MR-Egger: lets the regression have an intercept, which estimates directional pleiotropy; valid if pleiotropic effects are independent of instrument strength (InSIDE).",
  wm = "Weighted median: the median of the ratio estimates weighted by precision; valid if instruments carrying at least half of the weight are valid.",
  mode = "Weighted mode: the most common ratio estimate (a smoothed density); valid if the largest group of instruments with the same estimate is valid.",
  smode = "Simple mode: the mode of the ratio estimates without weights; a check on the weighted mode.",
  smed = "Simple median: the unweighted median of the ratio estimates; a check on the weighted median.",
  presso = "MR-PRESSO: tests for horizontal pleiotropy, removes outlying variants and reruns IVW; it simulates the data many times, so it takes longer with many instruments.",
  divw = "Debiased IVW: corrects the bias of IVW from weak instruments while keeping all of them.",
  raps = "MR-RAPS: a profile-likelihood estimate robust to weak instruments and to some pleiotropy.",
  conmix = "Contamination mixture: models the instruments as a mixture of valid and invalid ones and finds the effect most supported by the valid group.",
  lasso = "MR-Lasso: gives each variant its own pleiotropy term with a penalty, and drops variants whose term is not zero.",
  egger_simex = "MR-Egger with SIMEX: corrects MR-Egger for measurement error in the instrument effects (weak-instrument dilution) by simulation.",
  egger_boot = "MR-Egger with bootstrap: MR-Egger with confidence intervals from resampling the variants.",
  radial = "Radial IVW: IVW on the radial scale, which shows each variant's contribution to heterogeneity and flags outliers.",
  radegger = "Radial MR-Egger: MR-Egger on the radial scale, with its outliers.",
  steiger = "Steiger direction test: checks that the instruments explain more variance in the exposure than in the outcome, as expected if the effect runs from exposure to outcome.",
  ivw_steiger = "IVW after Steiger filtering: IVW without the variants that explain more variance in the outcome than in the exposure.",
  loo = "Leave-one-out: IVW repeated without each variant in turn, to see whether one variant decides the result.",
  bwmr = "Bayesian weighted MR: down-weights outlying variants and accounts for weak-instrument uncertainty with variational Bayes.",
  grip = "MR-GRIP: a regression robust to genotype coding, which removes the dependence of MR-Egger on the choice of reference allele.",
  mrmix = "MRMix: a mixture model that estimates the effect from the variants with no pleiotropy.",
  wmode_nome = "Weighted mode (NOME): the weighted mode assuming no measurement error in the instrument effects.",
  smode_nome = "Simple mode (NOME): the simple mode assuming no measurement error.",
  pwm = "Penalised weighted median: the weighted median with outlying variants down-weighted.",
  rucker = "Rucker model selection: chooses between IVW and MR-Egger (fixed or random effects) by their heterogeneity.",
  sign = "Sign concordance test: whether more variants than expected by chance have the same direction of effect on exposure and outcome.",
  pivw = "Penalised IVW: IVW with variants that contribute most to heterogeneity down-weighted.",
  maxlik = "Maximum likelihood: models the uncertainty in both the exposure and the outcome effects, unlike IVW.",
  mbe = "Mode-based estimate: the mode of the ratio estimates, with bootstrap intervals.",
  cml = "MR-cML: constrained maximum likelihood that selects invalid instruments; it fits many models, so it takes longer.",
  robust_ivw = "Robust IVW: IVW with robust regression, which down-weights outliers.",
  mvivw = "MVMR-IVW: the direct effect of each exposure, holding the others constant, combined over all joint instruments.",
  mvegger = "MVMR-Egger: multivariable MR with an intercept for directional pleiotropy not explained by the exposures.",
  mvmedian = "MVMR median: the multivariable weighted median, valid if instruments carrying half of the weight are valid.",
  mvlasso = "MVMR-Lasso: multivariable MR that drops variants with pleiotropy not explained by the exposures.",
  mvqhet = "MVMR Q-minimisation: estimates robust to weak instruments by minimising the heterogeneity statistic Q.",
  mvcml = "MVMR-cML: multivariable constrained maximum likelihood; it fits many models, so it takes longer.")
METHOD_SLOW <- c("presso", "cml", "mvcml", "egger_boot", "egger_simex", "mbe", "conmix", "mrmix", "bwmr", "mvqhet")

.live <- new.env()
.method_label <- function(id) { ad <- tryCatch(param_adapters()[[id]], error = function(e) NULL); lb <- ad$label %||% (names(METHOD_CHOICES)[METHOD_CHOICES == id])[1]; if (is.null(lb) || is.na(lb)) id else lb }

#' Called once a second while a pair's scripts run (contracts' run_pair_script, watched). Reports the method that started, the one that
#' finished or failed, and a reminder when a slow method runs for long.
live_scripts <- function(dir, console_name = "run_all_console.txt") {
  key <- paste(dir, console_name); st <- .live[[key]]
  ids <- sub("^[0-9]+_(.*)\\.R$", "\\1", setdiff(list.files(file.path(dir, "02_scripts"), pattern = "^[0-9]+_.*\\.R$"), "00_run_all.R"))
  if (!length(ids)) return(invisible())
  replay <- grepl("replay", console_name)
  if (is.null(st)) { st <- list(cur = "", t = Sys.time(), told = FALSE)
    .say(if (replay) sprintf("Replay check: rerunning the %d method scripts in a fresh R session from the saved data", length(ids)) else sprintf("Running %d method scripts for this pair, one after another", length(ids)),
         if (replay) "Every number is compared with the first run; identical results show that the analysis can be reproduced from the files in the pair folder alone." else "Each method runs on its own: if one stops with an error, the others still finish, and the error is named.") }
  ff <- file.path(dir, "03_console", "failed_scripts.txt"); failed <- if (file.exists(ff)) readLines(ff, warn = FALSE) else character(); failed <- sub("^[0-9]+_|\\.R$", "", failed); failed <- sub("\\.R$", "", failed)
  done <- ids[file.exists(file.path(dir, "04_tables", paste0(ids, "_estimates.csv"))) | ids %in% failed]
  nxt <- setdiff(ids, done)[1]; nxt <- if (is.na(nxt)) "" else nxt
  if (!identical(nxt, st$cur)) {
    # Round 80: every method that finished since the last look is reported (fast ones finish between two looks and were skipped)
    newly <- setdiff(done, st$done %||% character())
    if (length(newly) && !replay) {
      bad <- newly[newly %in% failed]; good <- setdiff(newly, bad)
      # Round 91: the method that was running is reported with its own duration, the quick ones that finished alongside it together
      if (st$cur %in% good) .say(sprintf("%s finished in %s (%d of %d)", .method_label(st$cur), .fmt_secs(difftime(Sys.time(), st$t, units = "secs")), match(st$cur, ids), length(ids)))
      oth <- setdiff(good, st$cur); if (length(oth)) .say(sprintf("Finished: %s (%d of %d)", paste(vapply(oth, .method_label, ""), collapse = ", "), length(done), length(ids)))
      for (b in bad) .say(sprintf("%s stopped with an error; the other methods continue (%d of %d)", .method_label(b), length(done), length(ids)))
    }
    st$done <- done
    if (nzchar(nxt) && !replay) .say(sprintf("Method %d of %d: %s", match(nxt, ids), length(ids), .method_label(nxt)), METHOD_WHY[[nxt]] %||% "")
    if (nzchar(nxt) && replay && (length(ids) <= 10L || match(nxt, ids) %% 10 == 1)) .say(sprintf("Replay check: method %d of %d (%s)", match(nxt, ids), length(ids), .method_label(nxt)))
    st$cur <- nxt; st$t <- Sys.time(); st$told <- FALSE
  } else if (nzchar(nxt) && !st$told && as.numeric(difftime(Sys.time(), st$t, units = "secs")) > 45) {
    .say(sprintf("%s%s is still running (started at %s)", if (replay) "Replay check: " else "", .method_label(nxt), format(st$t, "%H:%M:%S")),   # Round 91: the replay's own
         if (nxt %in% METHOD_SLOW) "This method refits the model many times (simulation, bootstrap or model selection), so it takes longer with many instruments. The app is working." else "Some methods take longer with many instruments. The app is working.")
    st$told <- TRUE
  }
  .live[[key]] <- st; invisible()
}
.fmt_secs <- function(s) { s <- round(as.numeric(s)); if (s < 60) sprintf("%d s", s) else sprintf("%d min %02d s", s %/% 60, s %% 60) }
live_reset <- function() rm(list = ls(.live), envir = .live)
