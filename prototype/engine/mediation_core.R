# ---- Round 57: Mendelian randomization for mediation (Carter et al. 2021, Eur J Epidemiol 36:465-478, PMID 33961203) ----------------
# A mediation study is a set of ordinary pairs, each run, saved and replayed like any other:
#   total   exposure X on outcome Y (univariable)                        -> total effect
#   xm      X on mediator M_j (univariable, X's instruments)             -> a_j
#   direct  MVMR of X and M_j on Y                                       -> direct effect of X (d_j) and effect of M_j given X (b_j)
#   joint   MVMR of X and all mediators on Y (only with 2+ mediators)    -> direct effect of X not through any of them
# The mediation script (04_across_pairs/mediation/mediation.R) reads their estimate tables and computes the indirect effect by the
# product of coefficients (two-step MR: a_j x b_j) and by the difference in coefficients (total - direct), and the proportion
# mediated, with delta-method standard errors, as for summary data (Carter et al. 2019 BMJ 365:l1855; Carter et al. 2021).

#' Simulated mediation world: X and every mediator have 50 instruments each; each mediator depends on X (DEMO_MED_A), and Y depends on
#' X directly (DEMO_EFFECT[X, Y]) and on each mediator (DEMO_EFFECT[M, Y]). So the total effect is d + sum_j a b_j, and the demo shows a
#' known proportion mediated. Mediators are taken from the simulated exposures.
DEMO_MED_A <- 0.4
demo_mediation_world <- function(x, ms, y) {
  tr <- c(x, sort(unique(ms))); ei <- match(tr, DEMO_EXPOSURES$id); oi <- match(y, DEMO_OUTCOMES$id)
  if (anyNA(ei) || is.na(oi)) stop("Unknown demo trait.", call. = FALSE)
  if (anyDuplicated(tr)) stop("The mediator must differ from the exposure.", call. = FALSE)
  K <- length(tr); m <- 50L; n <- K * m; grp <- rep(seq_len(K), each = m)
  set.seed(.seed_of(paste(c("med", tr, y), collapse = "|")))
  own <- runif(n, 0.05, 0.15) * sample(c(-1, 1), n, replace = TRUE)   # strong instruments and a precise outcome, so the demo sits near its truth
  B <- matrix(0, n, K); B[, 1] <- ifelse(grp == 1L, own, 0)
  if (K > 1) for (j in 2:K) B[, j] <- ifelse(grp == j, own, 0) + DEMO_MED_A * B[, 1]
  S <- matrix(runif(n * K, 0.004, 0.009), n, K); Bo <- B + matrix(rnorm(n * K, 0, S), n, K)
  sy <- runif(n, 0.005, 0.012); eff <- DEMO_EFFECT[cbind(tr, rep(y, K))]
  by <- as.vector(B %*% eff) + rnorm(n, 0, 0.002) + rnorm(n, 0, sy)
  list(traits = tr, names = DEMO_EXPOSURES$name[ei], n_x = DEMO_EXPOSURES$n[ei], n_y = DEMO_OUTCOMES$n[oi], y_name = DEMO_OUTCOMES$name[oi],
       B = round(Bo, 6), S = round(S, 6), by = round(by, 6), sy = round(sy, 6), grp = grp, ids = sprintf("demo:%d:T:C", 300000L + seq_len(n)),
       truth = list(direct = eff[1], a = DEMO_MED_A, b = eff[-1], total = eff[1] + DEMO_MED_A * sum(eff[-1])))
}
.demo_uni <- function(w, rows, bx, sx, by, sy, nx, ny, xl, yl, xid, yid) {
  d <- data.frame(variant_id = w$ids[rows], ea = "T", oa = "C", beta_x = bx, se_x = sx, eaf_x = 0.3, beta_y = by, se_y = sy, eaf_y = 0.3,
                  p_x = signif(2 * pnorm(-abs(bx / sx)), 6), p_y = signif(2 * pnorm(-abs(by / sy)), 6), n_x = nx, n_y = ny, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  list(exposure_label = xl, outcome_label = yl, simulated = TRUE, anc_x = "EUR", anc_y = "EUR", ledger = NULL,
       input = list(kind = "pair_iv", data = d, exposure_id = paste0("demo:", xid), outcome_id = paste0("demo:", yid), outcome_scale = "beta", exposure_scale = "beta"))
}
#' The input of one pair of a mediation study. med = list(role, x, ms (all mediators), y, m (this pair's mediator)).
mediation_pair_input <- function(med, pthr = 5e-8) {
  role <- med$role; x <- med$x; y <- med$y; ms <- unlist(med$ms); m <- med$m
  if (!role %in% c("total", "xm", "direct", "joint")) stop("Unknown mediation step.", call. = FALSE)
  lib <- startsWith(c(x, ms, y), "lib:")
  if (length(unique(lib)) > 1) stop("Simulated traits cannot be mixed with imported ones. Choose only simulated or only imported datasets.", call. = FALSE)
  # Round 72c: with OpenGWAS downloads a mediator is used in two forms, prepared by mediation_prepare(): looked up at the exposure's
  # instruments (exposure -> mediator) and as its own instruments (the multivariable steps). m and ms stay the labels the steps are matched by.
  m_out <- med$m_out %||% m; m_inst <- med$m_inst %||% m; ms_inst <- unlist(med$ms_inst) %||% ms
  if (all(lib)) return(switch(role, total = pair_input(x, y, pthr), xm = pair_input(x, m_out, pthr),
                              direct = mv_pair_input(c(x, m_inst), y, pthr), joint = mv_pair_input(c(x, ms_inst), y, pthr)))
  w <- demo_mediation_world(x, ms, y); K <- length(w$traits); j <- match(m, w$traits)
  if (role == "total") return(.demo_uni(w, w$grp == 1L, w$B[w$grp == 1L, 1], w$S[w$grp == 1L, 1], w$by[w$grp == 1L], w$sy[w$grp == 1L], w$n_x[1], w$n_y, w$names[1], w$y_name, x, y))
  if (role == "xm") { r <- w$grp == 1L; return(.demo_uni(w, r, w$B[r, 1], w$S[r, 1], w$B[r, j], w$S[r, j], w$n_x[1], w$n_x[j], w$names[1], w$names[j], x, m)) }
  ks <- if (role == "direct") c(1L, j) else seq_len(K); r <- w$grp %in% ks
  d <- data.frame(variant_id = w$ids[r], ea = "T", oa = "C", stringsAsFactors = FALSE)
  for (i in seq_along(ks)) { d[[sprintf("beta_x_%d", i)]] <- w$B[r, ks[i]]; d[[sprintf("se_x_%d", i)]] <- w$S[r, ks[i]]; d[[sprintf("n_x_%d", i)]] <- w$n_x[ks[i]] }
  d$beta_y <- w$by[r]; d$se_y <- w$sy[r]; d$n_y <- w$n_y; d$mr_keep <- TRUE
  list(exposure_label = paste(w$names[ks], collapse = " + "), outcome_label = w$y_name, simulated = TRUE, anc_x = "EUR", anc_y = "EUR", exposure_labels = w$names[ks],
       ledger = list(mv = list(exposures = w$names[ks], own = rep(50L, length(ks)), union = sum(r), clumped = sum(r), final = sum(r))),
       input = list(kind = "multi_iv", data = d, exposure_id = paste(paste0("demo:", w$traits[ks]), collapse = "|"), outcome_id = paste0("demo:", y), outcome_scale = "beta", exposure_scale = "beta"))
}

#' Lines of the mediation script. spec: x_id, y_id, x_label, y_label, scale, ci_level, uni_method, total (folder), mediators = list of
#' list(id, label, xm (folder), direct (folder)), joint (folder or NULL).
mediation_script <- function(spec) {
  q <- function(s) encodeString(s, quote = "\"")
  med_rows <- vapply(spec$mediators, function(m) sprintf("  list(id = %s, label = %s, xm = %s, direct = %s)", q(m$id), q(m$label), q(m$xm), q(m$direct)), "")
  c("# MR Studio: mediation analysis by Mendelian randomization (Carter et al. 2021, Eur J Epidemiol 36:465-478; Carter et al. 2019, BMJ 365:l1855).",
    "# Reads the estimate tables of the pairs below (run from the study folder) and writes 04_across_pairs/mediation/mediation_results.csv.",
    "# Product of coefficients (two-step MR): indirect = a x b, with a the effect of the exposure on the mediator (univariable) and b the effect of",
    "# the mediator on the outcome given the exposure (MVMR-IVW). Difference in coefficients: indirect = total - direct, with the direct effect from",
    "# the same MVMR. Proportion mediated = indirect / total. Delta-method standard errors, treating the estimates as independent, as for summary data.",
    "# ---- Settings ----",
    sprintf("EXPOSURE_ID <- %s   # %s", q(spec$x_id), spec$x_label), sprintf("OUTCOME_ID <- %s   # %s", q(spec$y_id), spec$y_label), sprintf("SCALE <- %s", q(spec$scale)),
    sprintf("CI_LEVEL <- %s", format(spec$ci_level)), sprintf("UNI_METHOD <- %s   # univariable method for the total effect and the exposure-mediator effect", q(spec$uni_method)),
    sprintf("TOTAL <- %s   # pair folder: exposure on outcome", q(spec$total)), "MEDIATORS <- list(", paste(med_rows, collapse = ",\n"), ")",
    sprintf("JOINT <- %s   # pair folder: MVMR of the exposure and all mediators (NULL with one mediator)", if (is.null(spec$joint)) "NULL" else q(spec$joint)),
    "# ---- Analysis ----",
    "one <- function(folder, method, estimand, exposure = NULL) {",
    "  e <- utils::read.csv(file.path(\"03_pairs\", folder, \"04_tables\", paste0(method, \"_estimates.csv\")), stringsAsFactors = FALSE)",
    "  e <- e[e$estimand == estimand, , drop = FALSE]; if (!is.null(exposure)) e <- e[e$exposure_id == exposure, , drop = FALSE]",
    "  if (nrow(e) != 1L) stop(sprintf(\"Expected one %s estimate of %s in %s.\", estimand, method, folder), call. = FALSE)",
    "  c(b = e$estimate, se = e$se) }",
    "z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)",
    "row <- function(mediator, method, quantity, est, se) data.frame(mediator = mediator, method = method, quantity = quantity, estimate = est, se = se,",
    "  ci_lower = est - z * se, ci_upper = est + z * se, p = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)",
    "tt <- one(TOTAL, UNI_METHOD, \"causal_effect\")",
    "out <- row(\"(none)\", \"univariable\", \"total\", tt[[\"b\"]], tt[[\"se\"]])",
    "for (m in MEDIATORS) {",
    "  a <- one(m$xm, UNI_METHOD, \"causal_effect\"); dd <- one(m$direct, \"mvivw\", \"direct_effect\", EXPOSURE_ID); b <- one(m$direct, \"mvivw\", \"direct_effect\", m$id)",
    "  ip <- a[[\"b\"]] * b[[\"b\"]]; sp <- sqrt(a[[\"b\"]]^2 * b[[\"se\"]]^2 + b[[\"b\"]]^2 * a[[\"se\"]]^2)   # product of coefficients (Sobel)",
    "  id <- tt[[\"b\"]] - dd[[\"b\"]]; sd <- sqrt(tt[[\"se\"]]^2 + dd[[\"se\"]]^2)                       # difference in coefficients",
    "  pp <- ip / tt[[\"b\"]]; spp <- abs(pp) * sqrt((a[[\"se\"]] / a[[\"b\"]])^2 + (b[[\"se\"]] / b[[\"b\"]])^2 + (tt[[\"se\"]] / tt[[\"b\"]])^2)",
    "  pd <- 1 - dd[[\"b\"]] / tt[[\"b\"]]; spd <- abs(dd[[\"b\"]] / tt[[\"b\"]]) * sqrt((dd[[\"se\"]] / dd[[\"b\"]])^2 + (tt[[\"se\"]] / tt[[\"b\"]])^2)",
    "  out <- rbind(out, row(m$label, \"two-step\", \"exposure on mediator\", a[[\"b\"]], a[[\"se\"]]), row(m$label, \"two-step\", \"mediator on outcome given exposure\", b[[\"b\"]], b[[\"se\"]]),",
    "    row(m$label, \"MVMR\", \"direct\", dd[[\"b\"]], dd[[\"se\"]]), row(m$label, \"product\", \"indirect\", ip, sp), row(m$label, \"difference\", \"indirect\", id, sd),",
    "    row(m$label, \"product\", \"proportion mediated\", pp, spp), row(m$label, \"difference\", \"proportion mediated\", pd, spd))",
    "}",
    "if (!is.null(JOINT)) {   # all mediators together: difference in coefficients only (the product method does not combine mediators that affect each other)",
    "  dj <- one(JOINT, \"mvivw\", \"direct_effect\", EXPOSURE_ID); ij <- tt[[\"b\"]] - dj[[\"b\"]]; sj <- sqrt(tt[[\"se\"]]^2 + dj[[\"se\"]]^2)",
    "  pj <- 1 - dj[[\"b\"]] / tt[[\"b\"]]; spj <- abs(dj[[\"b\"]] / tt[[\"b\"]]) * sqrt((dj[[\"se\"]] / dj[[\"b\"]])^2 + (tt[[\"se\"]] / tt[[\"b\"]])^2)",
    "  out <- rbind(out, row(\"(all together)\", \"MVMR\", \"direct\", dj[[\"b\"]], dj[[\"se\"]]), row(\"(all together)\", \"difference\", \"indirect\", ij, sj),",
    "    row(\"(all together)\", \"difference\", \"proportion mediated\", pj, spj))",
    "}",
    "out$scale <- SCALE",
    "dir.create(file.path(\"04_across_pairs\", \"mediation\"), recursive = TRUE, showWarnings = FALSE)",
    "utils::write.csv(out, file.path(\"04_across_pairs\", \"mediation\", \"mediation_results.csv\"), row.names = FALSE)",
    "out")
}

#' Write and run the mediation script for a study. pairs: the study's pair records (with med and folder). Returns the results table, or
#' NULL when a step has not run yet. The table is what the script writes, so the app shows exactly what the script reproduces.
mediation_run <- function(h, pairs, ci_level = 0.95, uni_method = "ivw_mre") {
  med <- Filter(function(p) !is.null(p$med) && !is.null(p$folder), pairs); if (!length(med)) return(NULL)
  role <- function(r, m = NULL) Filter(function(p) identical(p$med$role, r) && (is.null(m) || identical(p$med$m, m)), med)
  tot <- role("total"); if (!length(tot)) return(NULL); f0 <- tot[[1]]$med; ms <- unlist(f0$ms)
  mm <- lapply(ms, function(m) { a <- role("xm", m); d <- role("direct", m); if (!length(a) || !length(d)) return(NULL)
    xs <- strsplit(.mv_ids(h, d[[1]]$folder), "|", fixed = TRUE)[[1]]
    list(id = xs[2], label = (unlist(d[[1]]$exNames) %||% xs)[2], xm = a[[1]]$folder, direct = d[[1]]$folder) })
  if (any(vapply(mm, is.null, TRUE))) return(NULL)
  jt <- role("joint"); x_id <- strsplit(.mv_ids(h, mm[[1]]$direct), "|", fixed = TRUE)[[1]][1]
  sc <- tryCatch(.est_df(pair_payload(h$dir, tot[[1]]$folder))$scale[1], error = function(e) "beta")
  spec <- list(x_id = x_id, y_id = f0$y, x_label = tot[[1]]$exName, y_label = tot[[1]]$outName, scale = sc %||% "beta", ci_level = ci_level, uni_method = uni_method,
               total = tot[[1]]$folder, mediators = mm, joint = if (length(jt) && length(ms) > 1) jt[[1]]$folder else NULL)
  rel <- file.path("04_across_pairs", "mediation"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE)
  sf <- file.path(h$dir, rel, "mediation.R"); writeLines(mediation_script(spec), sf)
  old <- setwd(h$dir); on.exit(setwd(old), add = TRUE)
  res <- sys.source(sf, envir = new.env(parent = globalenv()), keep.source = FALSE); res <- utils::read.csv(file.path(rel, "mediation_results.csv"), stringsAsFactors = FALSE)
  manifest_register(h$dir, file.path(rel, c("mediation.R", "mediation_results.csv")), role = "mediation")
  attr(res, "spec") <- spec; res
}
#' The exposure IDs of an MVMR pair, in model order, from its MVMR-IVW estimate table.
.mv_ids <- function(h, folder) { e <- utils::read.csv(file.path(h$dir, "03_pairs", folder, "04_tables", "mvivw_estimates.csv"), stringsAsFactors = FALSE)
  paste(unique(e$exposure_id[e$estimand == "direct_effect"]), collapse = "|") }

# ---- report text ----
MEDIATION_METHODS_TEXT <- "Mediation was assessed with Mendelian randomization [Carter 2021]. The total effect of the exposure on the outcome was estimated by univariable MR, the effect of the exposure on each mediator by univariable MR with the exposure's instruments, and the direct effect of the exposure and the effect of each mediator given the exposure by multivariable MR (MVMR-IVW) with the instruments of both [Sanderson 2019]. The indirect effect was estimated by the product of coefficients (two-step MR: the effect of the exposure on the mediator times the effect of the mediator on the outcome given the exposure) and by the difference in coefficients (total minus direct effect), and the proportion mediated as the indirect effect divided by the total effect; standard errors were obtained with the delta method, treating the estimates as independent, as is usual with summary data [Carter 2019]."
MEDIATION_LIMITS_TEXT <- "The mediation estimates assume no interaction between the exposure and the mediator, linear effects, valid instruments for both the exposure and the mediator, and no overlap between the samples of the exposure, mediator and outcome GWAS; the delta-method intervals ignore the correlation between the estimates. A proportion mediated is hard to interpret when the total effect is small or its confidence interval includes zero, or when the direct and indirect effects have opposite signs."
mediation_results_text <- function(res, spec) {
  if (is.null(res) || !nrow(res)) return(character())
  log_or <- identical(spec$scale, "log_or"); f <- function(r) if (log_or) sprintf("OR %.2f (%s CI %.2f to %.2f)", exp(r$estimate), "95%", exp(r$ci_lower), exp(r$ci_upper)) else sprintf("%.3f (%s CI %.3f to %.3f)", r$estimate, "95%", r$ci_lower, r$ci_upper)
  pc <- function(r) sprintf("%.0f%% (%s CI %.0f%% to %.0f%%)", 100 * r$estimate, "95%", 100 * r$ci_lower, 100 * r$ci_upper)
  g <- function(m, me, q) res[res$mediator == m & res$method == me & res$quantity == q, , drop = FALSE]
  out <- sprintf("Mediation of the effect of %s on %s. Total effect: %s.", spec$x_label, spec$y_label, f(res[res$quantity == "total", ][1, ]))
  for (m in setdiff(unique(res$mediator), c("(none)", "(all together)"))) {
    out <- c(out, sprintf("Through %s: effect of the exposure on the mediator %s; effect of the mediator on the outcome given the exposure %s; direct effect %s; indirect effect %s by the product of coefficients and %s by the difference in coefficients; proportion mediated %s (product) and %s (difference).",
      m, sprintf("%.3f (%s CI %.3f to %.3f)", g(m, "two-step", "exposure on mediator")$estimate, "95%", g(m, "two-step", "exposure on mediator")$ci_lower, g(m, "two-step", "exposure on mediator")$ci_upper),
      f(g(m, "two-step", "mediator on outcome given exposure")), f(g(m, "MVMR", "direct")), f(g(m, "product", "indirect")), f(g(m, "difference", "indirect")),
      pc(g(m, "product", "proportion mediated")), pc(g(m, "difference", "proportion mediated"))))
  }
  if (any(res$mediator == "(all together)")) out <- c(out, sprintf("All mediators together: direct effect %s; proportion mediated %s (difference in coefficients).", f(g("(all together)", "MVMR", "direct")), pc(g("(all together)", "difference", "proportion mediated"))))
  tt <- res[res$quantity == "total", ][1, ]
  if (tt$ci_lower <= 0 && tt$ci_upper >= 0) out <- c(out, "The total effect's confidence interval includes no effect, so the proportions mediated should not be interpreted.")
  out
}

# ---- Round 72c: mediation from OpenGWAS downloads ------------------------------------------------------------------------------------
#' For exposure x, mediators ms and outcome y (library ids): each mediator in both forms (outcome at x's instruments; own instruments, with
#' x's instrument settings), then the multivariable sets of x with each mediator, and with all of them, on y. Imported full summary
#' statistics serve both forms as they are. Returns list(m_out, m_inst, notes).
mediation_prepare <- function(x, ms, y, jwt = NULL, offline = FALSE, root = studies_root(), get_exposure = ogw_get_exposure, get_outcome = ogw_get_outcome, backend = ogw_backend(), panel = NULL) {
  ms <- unlist(ms); rx <- library_dataset(sub("^lib:", "", x), root)$recipe; ry <- library_dataset(sub("^lib:", "", y), root)$recipe
  if (!.is_ogw_instruments(rx)) return(list(m_out = as.list(ms), m_inst = as.list(ms), notes = character()))   # full summary statistics: nothing to prepare
  notes <- character(); m_out <- m_inst <- character(length(ms))
  for (i in seq_along(ms)) { rm <- library_dataset(sub("^lib:", "", ms[i]), root)$recipe; src <- rm$source
    if (!identical(src$type, "opengwas")) { m_out[i] <- m_inst[i] <- ms[i]; next }
    if (isTRUE(rm$usable_as_exposure)) { m_inst[i] <- ms[i]
      r <- ogw_ensure_outcome(src, rm$trait, rx$trait, ry$source$params %||% list(proxies = TRUE, rsq = 0.8, maf_threshold = 0.3), jwt, offline, root, get_outcome); m_out[i] <- paste0("lib:", r$name); notes <- c(notes, r$note)
    } else { m_out[i] <- ms[i]
      r <- ogw_ensure_exposure(src, rm$trait, rx$source$params, jwt, offline, root, get_exposure); m_inst[i] <- paste0("lib:", r$name); notes <- c(notes, r$note) } }
  for (i in seq_along(ms)) { r <- ogw_mv_build(c(x, m_inst[i]), y, jwt, offline, root, backend, panel); if (!isTRUE(r$reused)) notes <- c(notes, sprintf("prepared \"%s\"", r$name)) }
  if (length(ms) > 1) { r <- ogw_mv_build(c(x, m_inst), y, jwt, offline, root, backend, panel); if (!isTRUE(r$reused)) notes <- c(notes, sprintf("prepared \"%s\"", r$name)) }
  list(m_out = as.list(m_out), m_inst = as.list(m_inst), notes = notes[nzchar(notes)])
}
