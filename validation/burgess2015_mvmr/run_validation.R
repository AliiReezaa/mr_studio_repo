# Validation study 4: Burgess and Thompson 2015, multivariable MR of LDL-c, HDL-c and triglycerides on CHD (Am J Epidemiol 181:251-260,
# doi:10.1093/aje/kwu283, PMID 25632051). Data: the 28 variants of Waterworth et al. 2010 used in the paper, as bundled with
# MendelianRandomization 0.10.0 (ldlc, hdlc, trig, chdlodds and their standard errors). Run from the repository root:
#   Rscript --vanilla validation/burgess2015_mvmr/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "burgess2015_mvmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mvval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
suppressMessages(library(MendelianRandomization))

# The paper's likelihood-based results (text of the Example section): causal odds ratios of CHD per 30% reduction in each lipid
# (the lipids are log-transformed, so a 30% reduction is a change of log(0.7) and OR = exp(beta * log(0.7))).
paper <- data.frame(lipid = c("LDL-c", "HDL-c", "Triglycerides"), or = c(0.50, 1.22, 0.77), lo = c(0.40, 0.91, 0.68), hi = c(0.62, 1.63, 0.87))

d <- data.frame(variant_id = sprintf("waterworth_%02d", seq_along(ldlc)), ea = "A", oa = "G",
                beta_x_1 = ldlc, se_x_1 = ldlcse, beta_x_2 = hdlc, se_x_2 = hdlcse, beta_x_3 = trig, se_x_3 = trigse,
                beta_y = chdlodds, se_y = chdloddsse, mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- list(kind = "multi_iv", data = d, exposure_id = "LDL-c|HDL-c|Triglycerides", outcome_id = "CHD", outcome_scale = "log_or", exposure_scale = "beta")
h <- studio_new("Burgess 2015 MVMR"); out <- list(); k <- 0L
for (model in c("fixed", "random")) {
  k <- k + 1L
  r <- studio_run(h, sprintf("Lipids (%s)", model), "CHD", inp, c("mvivw", "mvegger", "mvmedian", "mvlasso"), idx = k, params = list(ci_level = 0.95, mv_model = model, mv_iter = 10000L))
  pl <- pair_payload(h$dir, r$pair_folder); e <- .est_df(pl); e <- e[e$estimand == "direct_effect", , drop = FALSE]
  direct <- MendelianRandomization::mr_mvivw(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse), by = chdlodds, byse = chdloddsse), model = model)
  iv <- e[e$method == "mvivw", , drop = FALSE]
  out[[k]] <- data.frame(model = model, lipid = paper$lipid, beta = iv$estimate, se = iv$se, or = exp(iv$estimate * log(0.7)), or_lo = exp(iv$ci_upper * log(0.7)), or_hi = exp(iv$ci_lower * log(0.7)),
                         paper_or = paper$or, paper_lo = paper$lo, paper_hi = paper$hi, max_diff_vs_package = max(abs(iv$estimate - direct@Estimate), abs(iv$se - direct@StdError)),
                         replay = pl$replay$status, all_ok = all(vapply(r$results, function(x) identical(x$status, "ok"), TRUE)), stringsAsFactors = FALSE)
  if (model == "random") others <- e[e$method != "mvivw", c("method", "exposure_id", "estimate", "ci_lower", "ci_upper", "p")]
}
# The paper's own estimator, reconstructed: the likelihood of its equation 1 with zero correlations (two-sample data). Profiling out
# the true variant-exposure associations gives sum_j (by_j - b.bx_j)^2 / v_j + log v_j with v_j = sy_j^2 + sum_k b_k^2 sx_jk^2, which
# allows for uncertainty in the lipid associations (MVMR-IVW does not). Maximised here (BFGS), with Wald intervals from the Hessian.
bx <- cbind(ldlc, hdlc, trig); sx <- cbind(ldlcse, hdlcse, trigse)
nll <- function(b) { v <- chdloddsse^2 + as.vector(sx^2 %*% b^2); r <- chdlodds - as.vector(bx %*% b); 0.5 * sum(r^2 / v + log(v)) }
ml <- stats::optim(res_start <- out[[1]]$beta, nll, method = "BFGS", hessian = TRUE, control = list(reltol = 1e-14)); mse <- sqrt(diag(solve(ml$hessian)))
res <- do.call(rbind, out)
res$lik_or <- rep(exp(ml$par * log(0.7)), 2); res$lik_lo <- rep(exp((ml$par + 1.96 * mse) * log(0.7)), 2); res$lik_hi <- rep(exp((ml$par - 1.96 * mse) * log(0.7)), 2)
# Criteria (set from the reconstruction above): the paper reports Bayesian posterior summaries from WinBUGS, which differ from the maximum
# likelihood estimate by up to 0.01 on the odds ratio scale on this data, and MVMR-IVW ignores the uncertainty in the lipid associations,
# which widens the gap to at most 0.015. Odds ratios must agree within 0.015 and fixed-effect interval limits within 0.035.
res$point_match <- abs(res$or - res$paper_or) <= 0.015
res$ci_within_0.035 <- abs(res$or_lo - res$paper_lo) <= 0.035 & abs(res$or_hi - res$paper_hi) <= 0.035
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
utils::write.csv(others, file.path(vdir, "results_robust_methods.csv"), row.names = FALSE)
print(res[, c("model", "lipid", "or", "or_lo", "or_hi", "lik_or", "lik_lo", "lik_hi", "paper_or", "paper_lo", "paper_hi", "point_match", "ci_within_0.035", "replay")], digits = 3)
fx <- res[res$model == "fixed", ]
ok <- all(res$all_ok) && all(res$replay == "identical") && all(res$max_diff_vs_package < 1e-12) && all(res$point_match) && all(fx$ci_within_0.035) && all(abs(res$lik_or - res$paper_or) <= 0.011)
cat(if (ok) "VALIDATION PASSED\n" else "VALIDATION FAILED\n"); unlink(tmp, recursive = TRUE); if (!ok) quit(status = 1)
