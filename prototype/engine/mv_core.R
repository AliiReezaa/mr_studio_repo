# Round 50: data for multivariable MR (contracts input kind multi_iv). One table with a beta_x_<k> / se_x_<k> / n_x_<k> column set per
# exposure, plus the outcome. Instruments are the union of every exposure's own instruments, clumped jointly, and every exposure's effect
# is looked up at every one of them, aligned to the same effect allele (TwoSampleMR::harmonise_data, as for outcomes). Losses are counted.

#' The simulated traits: K exposures, 30 instruments each. Every variant also affects the other exposures (correlated exposures), and the
#' outcome depends on the exposures through their direct effects (DEMO_EFFECT), so MVMR can recover them where univariable MR cannot.
demo_mv_input <- function(exp_ids, out_id) {
  ei <- match(exp_ids, DEMO_EXPOSURES$id); oi <- match(out_id, DEMO_OUTCOMES$id)
  if (anyNA(ei) || is.na(oi)) stop("Unknown demo trait.")
  K <- length(ei); m <- 30L; n <- K * m
  set.seed(.seed_of(paste(c("mv", exp_ids), collapse = "|")))
  own <- runif(n, 0.03, 0.10) * sample(c(-1, 1), n, replace = TRUE); grp <- rep(seq_len(K), each = m)
  bx <- sapply(seq_len(K), function(k) ifelse(grp == k, own, 0.35 * own + rnorm(n, 0, 0.01)))
  sx <- matrix(runif(n * K, 0.004, 0.009), n, K)
  set.seed(.seed_of(paste(c("mv", exp_ids, out_id), collapse = "|")))
  sy <- runif(n, 0.008, 0.02); by <- as.vector(bx %*% DEMO_EFFECT[ei, oi]) + rnorm(n, 0, 0.002) + rnorm(n, 0, sy)
  d <- data.frame(variant_id = sprintf("demo:%d:G:A", 200000L + seq_len(n)), ea = "G", oa = "A", stringsAsFactors = FALSE)
  for (k in seq_len(K)) { d[[sprintf("beta_x_%d", k)]] <- round(bx[, k], 6); d[[sprintf("se_x_%d", k)]] <- round(sx[, k], 6); d[[sprintf("n_x_%d", k)]] <- DEMO_EXPOSURES$n[ei[k]] }
  d$beta_y <- round(by, 6); d$se_y <- round(sy, 6); d$n_y <- DEMO_OUTCOMES$n[oi]; d$mr_keep <- TRUE
  list(exposure_label = paste(DEMO_EXPOSURES$name[ei], collapse = " + "), outcome_label = DEMO_OUTCOMES$name[oi], simulated = TRUE, anc_x = "EUR", anc_y = "EUR",
       exposure_labels = DEMO_EXPOSURES$name[ei], ledger = list(mv = list(exposures = DEMO_EXPOSURES$name[ei], own = rep(m, K), union = n, clumped = n, final = n)),
       input = list(kind = "multi_iv", data = d, exposure_id = paste(paste0("demo:", exp_ids), collapse = "|"), outcome_id = paste0("demo:", out_id), outcome_scale = "beta", exposure_scale = "beta"))
}

.lib_fmt <- function(D, type) {
  D2 <- data.frame(SNP = D$key, beta = D$beta, se = D$se, effect_allele = D$ea, other_allele = D$oa, stringsAsFactors = FALSE)
  if ("eaf" %in% names(D) && any(is.finite(D$eaf))) D2$eaf <- D$eaf
  if ("p" %in% names(D)) D2$pval <- D$p
  if ("n" %in% names(D)) D2$samplesize <- D$n
  suppressMessages(TwoSampleMR::format_data(D2, type = type, snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize"))
}

#' Imported or downloaded datasets. Every exposure must hold full summary statistics (it is looked up at the other exposures' instruments),
#' like an outcome. The union is clumped on this computer (PLINK + 1000 Genomes, r2 < 0.001 within 10,000 kb) by the smallest P over the exposures.
library_mv_input <- function(exp_ids, out_id, pthr = 5e-8, root = studies_root(), clump = function(snp, p, pop) ld_clump_here(snp, p, 0.001, 10000, pop)) {
  ds <- lapply(exp_ids, function(e) library_dataset(sub("^lib:", "", e), root)); K <- length(exp_ids)
  for (k in seq_len(K)) if (!isTRUE(ds[[k]]$recipe$usable_as_outcome) || !is.null(ds[[k]]$recipe$source$for_exposure))   # checked first: it decides whether MVMR is possible at all
    stop("\"", ds[[k]]$recipe$trait, "\" holds only its own instruments, so its effects at the other exposures' instruments cannot be looked up. For MVMR every exposure must be imported or downloaded with full summary statistics (as an outcome).")
  ps <- lapply(exp_ids, function(e) library_pair_input(e, out_id, pthr, root))   # each exposure's own instruments, harmonised with the outcome
  labs <- vapply(ps, `[[`, "", "exposure_label")
  own <- lapply(ps, function(p) { d <- p$input$data; d[d$mr_keep, , drop = FALSE] })
  U <- do.call(rbind, lapply(seq_len(K), function(k) data.frame(key = own[[k]]$variant_id, ea = own[[k]]$ea, oa = own[[k]]$oa, eaf = own[[k]]$eaf_x, p = own[[k]]$p_x,
         beta_y = own[[k]]$beta_y, se_y = own[[k]]$se_y, n_y = own[[k]]$n_y %||% NA_real_, from = k, stringsAsFactors = FALSE)))
  U <- U[order(U$p), , drop = FALSE]; U <- U[!duplicated(U$key), , drop = FALSE]; led <- list(exposures = labs, own = vapply(own, nrow, 1L), union = nrow(U))
  if (!all(grepl("^rs[0-9]+$", U$key))) stop("The joint instruments must be clumped together, which needs rsIDs; at least one dataset identifies variants by position only.")
  pop <- ld_pop_for(ds[[1]]$recipe$facts$ancestry)
  if (is.na(pop)) stop("The joint instruments are clumped with a 1000 Genomes panel, which needs a European, East Asian, African, Admixed American or South Asian ancestry for the first exposure.")
  cl <- tryCatch(clump(U$key, U$p, pop), error = function(e) stop("The joint instruments of the exposures must be clumped together, and that needs the 1000 Genomes reference panel on this computer (Settings). ", conditionMessage(e), call. = FALSE))
  U <- U[U$key %in% cl$variant_id[cl$status == "index"], , drop = FALSE]; led$clumped <- nrow(U)
  d <- data.frame(variant_id = U$key, ea = U$ea, oa = U$oa, stringsAsFactors = FALSE); found <- integer(K)
  ref <- .lib_fmt(data.frame(key = U$key, beta = 1, se = 1, ea = U$ea, oa = U$oa, eaf = U$eaf), "exposure")
  # Palindromic variants are kept (action 1) only when every exposure was declared to be on the forward strand, as library_pair_input
  # does for a pair (Round 26); otherwise they are aligned by allele frequency or dropped (action 2). Found in the Rees 2017 validation.
  fwd <- all(vapply(ds, function(x) identical(x$recipe$facts$strand, "forward"), TRUE)); led$palindromic_kept <- fwd
  for (k in seq_len(K)) {
    X <- ds[[k]]$data; X$key <- .vkey(X)$key; X <- X[X$key %in% U$key, , drop = FALSE]
    H <- if (nrow(X)) suppressWarnings(suppressMessages(TwoSampleMR::harmonise_data(ref, .lib_fmt(X, "outcome"), action = if (fwd) 1 else 2))) else NULL
    H <- if (is.null(H)) data.frame(SNP = character(), beta.outcome = numeric(), se.outcome = numeric(), mr_keep = logical()) else H[H$mr_keep, , drop = FALSE]
    m <- match(d$variant_id, H$SNP); found[k] <- sum(!is.na(m))
    d[[sprintf("beta_x_%d", k)]] <- H$beta.outcome[m]; d[[sprintf("se_x_%d", k)]] <- H$se.outcome[m]
    nn <- if ("samplesize.outcome" %in% names(H)) H$samplesize.outcome[m] else rep(NA_real_, nrow(d))
    d[[sprintf("n_x_%d", k)]] <- ifelse(is.finite(nn), nn, suppressWarnings(as.numeric(ds[[k]]$recipe$facts$n %||% NA)))
  }
  d$beta_y <- U$beta_y; d$se_y <- U$se_y; d$n_y <- U$n_y
  d$mr_keep <- stats::complete.cases(d[, c(sprintf("beta_x_%d", seq_len(K)), sprintf("se_x_%d", seq_len(K)), "beta_y", "se_y")])
  led$found <- found; led$final <- sum(d$mr_keep)
  list(exposure_label = paste(labs, collapse = " + "), outcome_label = ps[[1]]$outcome_label, simulated = FALSE, anc_x = ps[[1]]$anc_x, anc_y = ps[[1]]$anc_y,
       exposure_labels = labs, ledger = list(mv = led),
       input = list(kind = "multi_iv", data = d, exposure_id = paste(exp_ids, collapse = "|"), outcome_id = ps[[1]]$input$outcome_id, outcome_scale = ps[[1]]$input$outcome_scale, exposure_scale = "beta"))
}

mv_pair_input <- function(exp_ids, out_id, pthr = 5e-8) {
  exp_ids <- unlist(exp_ids); lib <- startsWith(exp_ids, "lib:")
  if (length(exp_ids) < 2) stop("Multivariable MR needs at least two exposures.")
  if (length(exp_ids) > 6) stop("Multivariable MR in MR Studio takes at most six exposures.")
  if (length(unique(c(lib, startsWith(out_id, "lib:")))) > 1) stop("Simulated traits cannot be mixed with imported ones. Choose only simulated or only imported datasets.")
  if (all(lib)) library_mv_input(exp_ids, out_id, pthr) else demo_mv_input(exp_ids, out_id)
}

#' What the Check and Methods screens need before an MVMR run (the same fields as pair_facts where they apply).
mv_pair_facts <- function(exp_ids, out_id, pthr = 5e-8, p = NULL) {
  if (is.null(p)) p <- mv_pair_input(exp_ids, out_id, pthr); d <- p$input$data; d <- d[d$mr_keep, , drop = FALSE]; K <- length(p$exposure_labels)
  fk <- vapply(seq_len(K), function(k) mean((d[[sprintf("beta_x_%d", k)]] / d[[sprintf("se_x_%d", k)]])^2), 0)
  list(ex_id = paste(unlist(exp_ids), collapse = "|"), out_id = out_id, exposure = p$exposure_label, outcome = p$outcome_label, instruments = nrow(d), mean_f = round(min(fk), 1), min_f = round(min(fk), 1),
       n_weak = 0L, palindromic = 0L, scale = p$input$outcome_scale, has_n = FALSE, simulated = isTRUE(p$simulated), anc_x = p$anc_x, anc_y = p$anc_y, ledger = p$ledger,
       ukb_both = FALSE, og_x = FALSE, og_y = FALSE, scale_x = "beta", steiger = list(ok = FALSE), i2gx = NA_real_, mv = TRUE, exposures = p$exposure_labels, mean_f_each = round(fk, 1))
}

# ---- report (Round 50) -------------------------------------------------------------------------------------------------------
# Round 53: how the covariance of the variant-exposure associations was set for the conditional F and Q_A.
.mv_cov_text <- function(st, labels, src = NULL) {
  K <- length(labels); if (K < 2) return("")
  if (identical(src, "individual-level data")) return("both allowed for covariance between the variant-exposure associations, computed from individual-level data in which the exposures were measured in the same participants (the covariance of each pair of single-variant regression slopes, from their residuals)")
  pr <- do.call(rbind, lapply(seq_len(K - 1), function(i) do.call(rbind, lapply((i + 1):K, function(j) data.frame(i = i, j = j, r = as.numeric(st[[sprintf("mv_r_%d_%d", i, j)]] %||% 0))))))
  if (all(pr$r == 0)) return("both assumed zero covariance between the variant-exposure associations, as for exposures from non-overlapping samples")
  sprintf("both allowed for covariance between the variant-exposure associations, taken as the phenotypic correlation of the exposures times the product of their standard errors (%s)",
          paste(sprintf("%s and %s r = %s", labels[pr$i], labels[pr$j], format(pr$r)), collapse = "; "))
}
MV_INTRO_TEXT <- "We performed a two-sample multivariable Mendelian randomization (MVMR) analysis using summary-level genetic associations, run in MR Studio, an application that writes every analysis step and setting into R scripts that can be re-run [Burgess 2023]."
#' The context of an MVMR pair for the report: its estimates (one direct effect per exposure and method), settings, ledger and data.
mv_context <- function(h, p) {
  pay <- pair_payload(h$dir, p$folder); est <- .est_df(pay)
  ce <- est[est$estimand == "direct_effect", , drop = FALSE]; ids <- unique(as.character(ce$exposure_id))   # in the order of the model
  labs <- as.character(unlist(p$exNames)); if (length(labs) != length(ids)) labs <- ids
  ce$exposure_label <- labs[match(as.character(ce$exposure_id), ids)]
  g <- function(m, n) Filter(function(d) identical(d$method, m) && identical(d$name, n), pay$diagnostics)
  list(mv = TRUE, p = p, pay = pay, est = est, ce = ce, methods = unlist(p$methods), simulated = isTRUE(p$simulated), settings = pay$settings %||% list(), ledger = p$ledger,
       cond_f = vapply(g("mvivw", "mv_cond_f"), function(d) as.numeric(d$value), 0), q = g("mvivw", "mv_q"), qa = g("mvivw", "mv_q_a"), labels = labs, cov_source = { n <- g("mvivw", "mv_pheno_cor"); if (!length(n)) n <- g("mvqhet", "mv_pheno_cor"); if (length(n)) n[[1]]$note else NULL }, vars = data.frame(), pthr = p$pthr %||% 5e-8)
}
mv_add_text <- function(txt, mvx) {
  st <- mvx[[1]]$settings; used <- unique(unlist(lapply(mvx, `[[`, "methods")))
  lib <- any(!vapply(mvx, `[[`, TRUE, "simulated"))
  txt$methods <- c(txt$methods, paste0("Multivariable MR (MVMR) estimated the direct effect of each exposure on the outcome, conditional on the other exposures in the model [Sanderson 2019]. ",
    if (lib) sprintf("Instruments were the union of the variants associated with any exposure at P < %s, clumped together (r2 < 0.001 within 10,000 kb, 1000 Genomes reference panel) using the smallest P over the exposures; the effects of every exposure at every instrument were taken from its full summary statistics and aligned to the same effect allele with TwoSampleMR::harmonise_data. ", format(mvx[[1]]$pthr, scientific = TRUE)) else "",
    "Methods (MendelianRandomization ", PIN_MENDELIANRANDOMIZATION, "): ", paste(c(if ("mvivw" %in% used) sprintf("MVMR-IVW (model: %s)", st$mv_model %||% "default"), if ("mvegger" %in% used) sprintf("MVMR-Egger, oriented to exposure %s [Rees 2017]", format(st$mv_orient %||% 1)),
      if ("mvmedian" %in% used) sprintf("the MVMR median with %s bootstrap iterations [Grant 2021]", format(st$mv_iter %||% 10000)), if ("mvlasso" %in% used) "MVMR-Lasso [Grant 2021]",
      if ("mvcml" %in% used) sprintf("MVMR-cML, constrained maximum likelihood with BIC selection of invalid variants%s [Lin 2023]", if (identical(st$mv_cml_dp, "ma")) "" else " and data perturbation"),
      if ("mvqhet" %in% used) sprintf("Q-statistic minimisation, which estimates the direct effects by minimising Q_A and so is robust to weak instruments, with an additive heterogeneity variance when the minimum exceeds its degrees of freedom and %s [Sanderson 2021]",
        if (identical(st$mv_qhet_ci, "bootstrap")) sprintf("bias-corrected and accelerated bootstrap intervals (%s iterations)", format(st$mv_qhet_iter %||% 1000)) else "jackknife standard errors (leaving out one variant at a time)")), collapse = "; "),
    ". Instrument strength was assessed with the conditional F statistic of each exposure (values below 10 indicate weak instruments), and heterogeneity with Cochran's Q of the MVMR-IVW fit and with Q_A, which also allows for the uncertainty of the variant-exposure associations [Sanderson 2021]; ", .mv_cov_text(st, mvx[[1]]$labels, mvx[[1]]$cov_source), "."))
  for (cx in mvx) {
    ce <- cx$ce; ex <- unique(ce$exposure_label); iv <- ce[ce$method == "mvivw", , drop = FALSE]
    line <- vapply(ex, function(x) { r <- ce[ce$exposure_label == x, , drop = FALSE]
      sprintf("%s: %s", x, paste(sprintf("%s %s", .short_mv(r$method), vapply(seq_len(nrow(r)), function(i) .ci_phrase(r[i, ], r$scale[i]), "")), collapse = "; ")) }, "")
    cf <- if (length(cx$cond_f)) paste(sprintf("%s %s", cx$labels[seq_along(cx$cond_f)], format(round(cx$cond_f, 1), nsmall = 1)), collapse = ", ") else ""
    q <- if (length(cx$q)) sprintf(" Cochran's Q of the MVMR-IVW fit was %.1f (P = %s)%s.", as.numeric(cx$q[[1]]$value), .fp(as.numeric(cx$q[[1]]$p)),
      if (length(cx$qa)) sprintf(" and Q_A was %.1f (P = %s)", as.numeric(cx$qa[[1]]$value), .fp(as.numeric(cx$qa[[1]]$p))) else "") else ""
    ic <- cx$est[cx$est$method == "mvegger" & cx$est$estimand == "pleiotropy_intercept", , drop = FALSE]
    led <- cx$ledger$mv
    txt$results <- c(txt$results, paste0(sprintf("MVMR of %s on %s (%s instruments). Direct effects: %s.", paste(ex, collapse = ", "), cx$p$outName, if (nrow(iv)) as.integer(iv$n_ivs[1]) else as.integer(led$final %||% NA), paste(line, collapse = ". ")),
      if (nzchar(cf)) sprintf(" Conditional F statistics: %s%s.", cf, if (any(cx$cond_f < 10, na.rm = TRUE)) "; at least one is below 10, so the estimates may be biased by weak instruments" else "") else "",
      q, if (nrow(ic)) sprintf(" The MVMR-Egger intercept was %s (P = %s).", .fnum(ic$estimate[1], 4), .fp(ic$p[1])) else "",
      if (!is.null(led$union) && !isTRUE(cx$simulated)) sprintf(" Of %d candidate instruments (%s), %d remained after joint clumping and %d had every exposure and the outcome.", as.integer(led$union), paste(sprintf("%s %d", led$exposures, as.integer(led$own)), collapse = ", "), as.integer(led$clumped), as.integer(led$final)) else ""))
    txt$legends <- c(txt$legends, sprintf("Figure. Direct effects of %s on %s from multivariable MR: for each exposure, the estimate and %s confidence interval of each method, conditional on the other exposures.", paste(ex, collapse = ", "), cx$p$outName, "95%"), "")
  }
  sim <- "SIMULATED DATA: at least one pair uses simulated demonstration traits. These results illustrate the software and must not be presented as scientific findings."
  if (any(vapply(mvx, `[[`, TRUE, "simulated")) && !(sim %in% txt$limitations)) txt$limitations <- c(txt$limitations, sim)
  txt$limitations <- c(txt$limitations, "Multivariable MR estimates are direct effects conditional on the exposures included in the model; they change if exposures are added or removed, and they assume that each variant affects the outcome only through the exposures in the model. The conditional F statistic assumes that the exposure samples do not overlap.")
  txt
}
.short_mv <- function(m) unname(c(mvivw = "MVMR-IVW", mvegger = "MVMR-Egger", mvmedian = "MVMR median", mvlasso = "MVMR-Lasso", mvqhet = "MVMR Q-min", mvcml = "MVMR-cML")[m])

# ---- Round 55: covariances of the variant-exposure associations from individual-level data (as MVMR::snpcov_mvmr, corrected) ----
#' For every variant j and pair of exposures (i, k), the covariance of the two single-variant regression slopes estimated in the same
#' people: cov(b_ij, b_kj) = sum(e_ij * e_kj) / ((n - 2) * sum((g_j - mean g_j)^2)), with e the residuals of exposure i (or k) on g_j and an
#' intercept. Its diagonal is exactly the squared standard error lm() reports. MVMR::snpcov_mvmr 0.4.9 uses the uncentred sum(g^2) * n
#' instead of (n - 2) * sum((g - mean)^2), which equals this only for centred genotypes; for 0/1/2 dosages it makes the variances too small.
#' Rows with a missing genotype or exposure are left out for that variant. Returns one row per variant with n and cov_x_i_j (i <= j).
mv_snpcov <- function(G, X) {
  G <- as.matrix(G); X <- as.matrix(X); K <- ncol(X); pr <- which(upper.tri(diag(K), diag = TRUE), arr.ind = TRUE); pr <- pr[order(pr[, 1], pr[, 2]), , drop = FALSE]
  rows <- lapply(seq_len(ncol(G)), function(j) {
    ok <- is.finite(G[, j]) & stats::complete.cases(X); g <- G[ok, j]; x <- X[ok, , drop = FALSE]; n <- length(g); gc <- g - mean(g); ssg <- sum(gc^2)
    if (n < 3 || ssg <= 0) return(c(n, rep(NA_real_, nrow(pr))))
    xc <- sweep(x, 2, colMeans(x)); e <- xc - outer(gc, as.vector(crossprod(gc, xc)) / ssg); V <- crossprod(e) / ((n - 2) * ssg)
    c(n, V[pr]) })
  out <- as.data.frame(do.call(rbind, rows)); names(out) <- c("n_cov", sprintf("cov_x_%d_%d", pr[, 1], pr[, 2]))
  cbind(data.frame(variant_id = colnames(G), stringsAsFactors = FALSE), out)
}

#' Read an individual-level file (one row per person): genotype dosage columns named by variant and one column per exposure. The
#' exposure columns are given in the order of the model; every other numeric column with values from 0 to 2 is taken as a genotype. Only the covariance
#' table is kept; the individual-level file is not copied into the study. A PLINK --recode A file (.raw) can be used as it is, with
#' the exposures added as columns.
mv_cov_from_file <- function(path, exposure_cols) {
  d <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  miss <- setdiff(exposure_cols, names(d)); if (length(miss)) stop(sprintf("Exposure column%s not in the file: %s.", if (length(miss) > 1) "s" else "", paste(miss, collapse = ", ")), call. = FALSE)
  num <- setdiff(names(d)[vapply(d, is.numeric, TRUE)], c(exposure_cols, "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))   # PLINK --recode A header columns are not genotypes
  dose <- vapply(num, function(k) { v <- d[[k]][is.finite(d[[k]])]; length(v) > 0 && min(v) >= 0 && max(v) <= 2 }, TRUE)
  gcols <- num[dose]   # genotype dosages lie in 0 to 2; other numeric columns (an ID, age) are ignored and reported
  if (!length(gcols)) stop("No genotype columns found: besides the exposures, the file needs one numeric column per variant (dosage 0 to 2), named by the variant ID.", call. = FALSE)
  out <- mv_snpcov(d[, gcols, drop = FALSE], d[, exposure_cols, drop = FALSE])
  # PLINK --recode A names a dosage column <variant>_<counted allele> (rs123_A); the variant ID is the part before it. The covariances
  # do not depend on which allele is counted (both slopes change sign together), so the allele is not needed.
  id <- sub("_[ACGTacgt]+$", "", out$variant_id); if (!anyDuplicated(id)) out$variant_id <- id
  attr(out, "meta") <- list(people = nrow(d), variants = length(gcols), exposures = exposure_cols, ignored = c(num[!dose], setdiff(names(d), c(num, exposure_cols)))); out
}

#' Add the covariance columns to an MVMR input by variant ID. Every variant of the model must have them; otherwise stop with the count.
mv_attach_cov <- function(input, cov) {
  d <- input$data; K <- sum(grepl("^beta_x_[0-9]+$", names(d))); need <- sprintf("cov_x_%d_%d", row(diag(K))[upper.tri(diag(K), diag = TRUE)], col(diag(K))[upper.tri(diag(K), diag = TRUE)])
  if (!all(need %in% names(cov))) stop(sprintf("The covariance table is for %s exposures but the model has %d.", format(round((sqrt(8 * sum(grepl("^cov_x_", names(cov))) + 1) - 1) / 2)), K), call. = FALSE)
  m <- match(d$variant_id, cov$variant_id); bad <- is.na(m) | !stats::complete.cases(cov[m, need, drop = FALSE])
  if (any(bad & d$mr_keep)) stop(sprintf("%d of the %d instruments have no covariance in the individual-level data (for example %s). Add them to the file or use the phenotypic correlations instead.",
                                         sum(bad & d$mr_keep), sum(d$mr_keep), paste(utils::head(d$variant_id[bad & d$mr_keep], 3), collapse = ", ")), call. = FALSE)
  for (k in need) d[[k]] <- cov[[k]][m]
  input$data <- d; input
}

#' Save a covariance table for a set of exposures in the study (02_datasets; summary-level only) and remember it in the study state.
mv_cov_save <- function(h, ex, cov, source_name) {
  key <- paste(unlist(ex), collapse = "|"); rel <- file.path("02_datasets", sprintf("mv_covariance_%s.csv", slugify(gsub("\\|", "_", key), 60L)))
  utils::write.csv(cov, file.path(h$dir, rel), row.names = FALSE); manifest_register(h$dir, rel, role = "mv_covariance")
  mt <- attr(cov, "meta") %||% list(); st <- Filter(function(x) !identical(x$key, key), ui_state_get(h$con, "mv_cov") %||% list())
  st[[length(st) + 1L]] <- list(key = key, ex = as.list(unlist(ex)), file = rel, source = source_name, people = mt$people, variants = mt$variants, exposure_cols = as.list(mt$exposures), ignored = as.list(mt$ignored))
  ui_state_set(h$con, "mv_cov", st); st
}
mv_cov_get <- function(h, ex) { key <- paste(unlist(ex), collapse = "|"); Filter(function(x) identical(x$key, key), ui_state_get(h$con, "mv_cov") %||% list()) }
mv_cov_clear <- function(h, ex) { key <- paste(unlist(ex), collapse = "|"); st <- Filter(function(x) !identical(x$key, key), ui_state_get(h$con, "mv_cov") %||% list()); ui_state_set(h$con, "mv_cov", st); st }
#' The MVMR input with the study's covariance table for these exposures merged in, if there is one.
mv_with_cov <- function(h, ex, input) {
  e <- mv_cov_get(h, ex); if (!length(e)) return(input)
  mv_attach_cov(input, utils::read.csv(file.path(h$dir, e[[1]]$file), stringsAsFactors = FALSE))
}
