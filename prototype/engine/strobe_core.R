# STROBE-MR checklist (Round 20): Skrivankova et al. 2021, JAMA 326:1614 (PMID 34698778) and the explanation and elaboration,
# BMJ 375:n2233 (PMID 34702754, item list read from PMC8546498). 20 items with 30 sub-items (10d has parts i and ii).
# The "asks" column is MR Studio's own short paraphrase, not the checklist wording; the official checklist is cited in the file.
# Status: FILLED = MR Studio wrote it from the saved study; PARTLY = MR Studio wrote part of it and the authors complete it;
# AUTHOR = only the authors can report it. Nothing is marked FILLED unless the text or file it points to exists in the package.

STROBE_ITEMS <- data.frame(stringsAsFactors = FALSE, matrix(ncol = 3, byrow = TRUE, c(
  "Title and abstract", "1", "Name Mendelian randomization as the design when it is a main aim",
  "Introduction", "2", "Background: rationale, the exposure, why a causal effect is plausible, why MR suits the question",
  "Introduction", "3", "Objectives: aims and prespecified causal hypotheses",
  "Methods", "4a", "Design, population, setting and dates of each data source",
  "Methods", "4b", "Eligibility, selection, sample size and any power calculation",
  "Methods", "4c", "Measurement, quality control and selection of the genetic variants",
  "Methods", "4d", "How exposure, outcome and other variables were defined and measured",
  "Methods", "4e", "Ethics approval and consent",
  "Methods", "5", "The instrumental variable assumptions used",
  "Methods", "6a", "Scale, units and model for quantitative variables",
  "Methods", "6b", "How the variants were combined and weighted",
  "Methods", "6c", "The MR estimator and its statistics",
  "Methods", "6d", "How missing data were handled",
  "Methods", "6e", "How multiple testing was handled",
  "Methods", "7", "Methods or knowledge used to assess the assumptions",
  "Methods", "8", "Sensitivity and additional analyses",
  "Methods", "9a", "Software, packages, versions and settings",
  "Methods", "9b", "Pre-registration of the protocol",
  "Results", "10a", "Numbers at each stage and reasons for exclusion",
  "Results", "10b", "Summary of the exposure, outcome and other data",
  "Results", "10c", "Heterogeneity between studies in meta-analysed data sources",
  "Results", "10d-i", "Two-sample MR: why the variant-exposure associations apply to the outcome sample",
  "Results", "10d-ii", "Two-sample MR: overlap of people between the exposure and outcome studies",
  "Results", "11a", "Variant-exposure and variant-outcome associations",
  "Results", "11b", "MR estimates with their uncertainty on an interpretable scale",
  "Results", "11c", "Relative effects converted to absolute risk, if useful",
  "Results", "11d", "Plots of the results",
  "Results", "12a", "How valid the assumptions appear",
  "Results", "12b", "Further statistics such as heterogeneity between variants",
  "Results", "13a", "Sensitivity analyses for violations of the assumptions",
  "Results", "13b", "Other sensitivity or additional analyses",
  "Results", "13c", "Assessment of the direction of the effect",
  "Results", "13d", "Comparison with non-MR estimates, where relevant",
  "Results", "13e", "Further plots, such as leave-one-out",
  "Discussion", "14", "Key results in relation to the objectives",
  "Discussion", "15", "Limitations, including the assumptions and the direction and size of bias",
  "Discussion", "16a", "Meaning: cautious interpretation in the light of other evidence",
  "Discussion", "16b", "Mechanism and careful causal wording",
  "Discussion", "16c", "Clinical or policy relevance",
  "Discussion", "17", "Generalisability to other populations, timings and exposure levels",
  "Other information", "18", "Funding, including of the source studies",
  "Other information", "19", "Where the data and the analysis code can be obtained",
  "Other information", "20", "Conflicts of interest")))
names(STROBE_ITEMS) <- c("section", "item", "asks")

#' The checklist for a study, as a table: item, what it asks, status, what MR Studio provides and where.
#' ctxs = pair contexts (report_core), txt = report_text(), plan = the study plan, items = package items being built.
strobe_rows <- function(ctxs, txt, plan = list(), items = character()) {
  used <- unique(unlist(lapply(ctxs, `[[`, "methods"))); np <- length(ctxs); sim <- any(vapply(ctxs, `[[`, TRUE, "simulated"))
  pn <- function(cx) sprintf("%s -> %s", cx$p$exName, cx$p$outName)
  each <- function(f) paste(vapply(ctxs, function(cx) { v <- f(cx); if (length(v) && nzchar(v)) paste0(if (np > 1) paste0(pn(cx), ": ") else "", v) else "" }, ""), collapse = " ")
  has <- function(m) any(m %in% used)
  st <- steiger_any <- any(vapply(ctxs, function(cx) !is.null(steiger_facts(cx)), TRUE))
  ovl <- vapply(ctxs, function(cx) { o <- cx$overlap; if (is.null(o)) "not answered" else switch(o$answer, none = "no shared participants (analyst's statement)", some = sprintf("some shared participants%s (analyst's statement)", if (!is.null(o$pct)) sprintf(", about %s%%", o$pct) else ""), "not known") }, "")
  src <- function(s) if (!is.null(s$mv)) sprintf("OpenGWAS %s, each looked up at the joint instruments", paste(unlist(s$mv$exposures), collapse = ", ")) else if (!is.null(s$og)) sprintf("OpenGWAS %s (%s%s)", s$og$id, s$ancestry, if (nzchar(as.character(s$n))) sprintf(", N %s", format(as.numeric(s$n), big.mark = ",", scientific = FALSE)) else "") else if (length(s$file) && !is.na(s$file)) sprintf("file %s (%s)", s$file, s$ancestry) else if (is.null(s$file)) "Data library dataset" else "simulated demonstration trait"
  inst <- function(cx) { L <- cx$ledger
    sel <- if (!is.null(L$clump)) sprintf("P < %s, clumped on the analyst's computer with PLINK (r2 < %s within %s kb, 1000 Genomes %s)", format(L$threshold), L$clump$r2, format(L$clump$kb, big.mark = ","), L$clump$pop)
           else if (!is.null(cx$src_x$og)) sprintf("P < %s, %s (r2 < %s within %s kb, %s panel)", format(cx$src_x$og$params$pval), if (identical(cx$src_x$og$params$where, "server_live")) "clumped live on the OpenGWAS server" else if (ogw_stored_clumps(cx$src_x$og$params)) "from the stored, pre-computed clumps of OpenGWAS" else "clumped on the OpenGWAS server", cx$src_x$og$params$r2, format(cx$src_x$og$params$kb, big.mark = ","), cx$src_x$og$params$pop)
           else if (!is.null(L)) sprintf("P < %s from a file declared independent by the analyst", format(L$threshold)) else "simulated instruments"
    paste0(    sprintf("%s; %d instruments used, mean F %.1f (minimum %.1f, %d below 10)", sel, nrow(cx$kept), cx$mean_f, cx$min_f, cx$n_weak), { fx <- tryCatch(f_stats_text(L$fstats %||% cx$fstats), error = function(e) ""); if (nzchar(fx)) paste0(". ", sub("[.]$", "", fx)) else "" },
      if (!is.null(L$maf)) paste0(". ", sub("[.]$", "", maf_text(L$maf))) else "") }   # Round 109
  r <- function(item, status, provides, where = "") list(item = item, status = status, provides = provides, where = where)
  M <- "Draft text, Methods"; R <- "Draft text, Results"; Lm <- "Draft text, Limitations"
  ests <- c(ivw_mre = "IVW", divw = "debiased IVW", raps = "MR-RAPS", conmix = "contamination mixture", lasso = "MR-Lasso", egger_simex = "MR-Egger with SIMEX correction", egger = "MR-Egger", wm = "weighted median", smed = "simple median", mode = "weighted mode", smode = "simple mode", presso = "MR-PRESSO", radial = "radial IVW", radegger = "radial MR-Egger")
  robust <- intersect(c("egger", "wm", "smed", "mode", "smode", "raps", "conmix", "lasso", "egger_simex"), used); detect <- intersect(c("presso", "radial", "radegger"), used)
  rows <- list(
    r("1", "AUTHOR", "Name the design, for example 'a two-sample Mendelian randomization study', in the title and abstract."),
    r("2", "AUTHOR", "Rationale and why MR suits the question are the authors' to write."),
    r("3", if (nzchar(plan$aim %||% "")) "PARTLY" else "AUTHOR", if (nzchar(plan$aim %||% "")) sprintf("The aim recorded in MR Studio: \"%s\". State the prespecified causal hypothesis.", plan$aim) else "No aim was recorded in MR Studio. State the aims and the prespecified causal hypothesis."),
    r("4a", "PARTLY", sprintf("Two-sample MR on summary-level data. Sources: %s. Add the setting and dates of the source studies.", each(function(cx) sprintf("exposure %s; outcome %s", src(cx$src_x), src(cx$src_y)))), M),
    r("4b", "PARTLY", sprintf("Sample sizes as recorded with the data (%s). No power calculation is done by MR Studio: add one if needed.", each(function(cx) sprintf("median N exposure %s, outcome %s", if (length(cx$nx) == 1L && is.finite(cx$nx)) format(cx$nx, big.mark = ",") else "unknown", if (length(cx$ny) == 1L && is.finite(cx$ny)) format(cx$ny, big.mark = ",") else "unknown"))), M),
    r("4c", "FILLED", paste0("Instrument selection and quality control: ", each(inst), paste0(". Alleles harmonised with TwoSampleMR::harmonise_data; ", if (any(vapply(ctxs, function(cx) identical(cx$ledger$strand, "forward"), TRUE))) "palindromic variants kept where both datasets were declared forward strand (action 1), otherwise aligned by allele frequency or dropped." else "palindromic variants aligned by allele frequency or dropped.")), paste0(M, "; supplement workbook (variants of every pair)")),
    r("4d", "PARTLY", sprintf("Traits as named in the data sources (%s), effect sizes per unit of the exposure%s. Add the definitions of the traits and, for a disease, the case criteria.", each(function(cx) sprintf("%s; %s", cx$p$exName, cx$p$outName)), if (any(vapply(ctxs, function(cx) any(cx$ce$scale == "log_or"), TRUE))) ", log odds ratios for a binary outcome" else ""), M),
    r("4e", "AUTHOR", "Summary statistics only were used. State that the source studies had ethics approval and consent, citing them."),
    r("5", "AUTHOR", "State the three instrumental variable assumptions (relevance, independence, exclusion restriction) and any added for the sensitivity analyses. MR Studio does not write this sentence."),
    r("6a", "FILLED", "Effects are per unit of the exposure as provided by the source GWAS (log odds ratios for a binary outcome). Each method combines the variant-level ratio estimates (outcome effect divided by exposure effect) as described in its methods sentence.", M),
    r("6b", "FILLED", "Variants were combined by inverse-variance weighting of their ratio estimates; the weights of every method are stated in its methods sentence and in its script.", paste0(M, "; scripts (02_scripts)")),
    r("6c", "FILLED", sprintf("Estimators: %s, with their settings, confidence intervals and P values.", paste(ests[intersect(names(ests), used)], collapse = ", ")), M),
    r("6d", "PARTLY", { miss <- trimws(each(function(cx) { L <- cx$ledger; if (is.null(L)) "" else sprintf("%d instruments absent from the outcome%s", as.integer(L$missing_in_outcome), if (!is.null(L$by_proxy) && L$by_proxy > 0) sprintf(" (%d more found through proxies)", as.integer(L$by_proxy)) else "") }))
      paste0(if (nzchar(miss)) paste0("Missing variants: ", miss, ". ") else "No instrument was missing from the outcome data. ", "Effects or standard errors that are missing are never imputed; such variants are left out and counted. Add anything specific to your data.") }, M),
    if (any(vapply(ctxs, function(cx) !is.null(cx$p$med), TRUE))) r("6e", "FILLED", "One mediation analysis: its steps estimate parts of one effect and are not separate hypotheses, so no correction for multiple testing was applied.") else   # Round 98
    r("6e", if (np > 1) "PARTLY" else "FILLED", if (np > 1) sprintf("%d exposure-outcome pairs were tested. MR Studio applies no correction; a Bonferroni threshold would be P < %s. State the approach you used.", np, format(signif(0.05 / np, 3))) else "One exposure-outcome pair: no correction for multiple testing was needed."),
    r("7", "PARTLY", paste0("Instrument strength (F statistic), heterogeneity (Cochran's Q)", if (has("egger")) ", MR-Egger intercept for directional pleiotropy" else "", if (st) ", Steiger directionality test" else "", if (length(detect)) sprintf(", outlier detection (%s)", paste(ests[detect], collapse = ", ")) else "", ", the analyst's statement on sample overlap and the declared ancestry. Add prior knowledge on pleiotropy of the variants."), M),
    r("8", if (length(robust) || length(detect)) "FILLED" else "AUTHOR", if (length(robust) || length(detect)) sprintf("Sensitivity analyses: %s.", paste(ests[c(robust, detect)], collapse = ", ")) else "No sensitivity method was run. Add robust methods (for example MR-Egger, weighted median) on the Methods screen.", M),
    r("9a", "FILLED", sprintf("MR Studio %s; R %s; %s%s; every analysis is a generated R script with all arguments written out, re-run in a fresh R session (replay check).", if (exists("MRS_VERSION")) MRS_VERSION else "(version unknown)", paste(R.version$major, R.version$minor, sep = "."), { sp <- software_pkgs(unique(unlist(lapply(ctxs, function(x) x$methods)))); paste(c(sp$ver, if (length(sp$own)) paste("MR Studio's own implementation of", paste(vapply(sp$own, .own_name, ""), collapse = ", "))), collapse = "; ") }, if (any(vapply(ctxs, function(cx) !is.null(cx$ledger$clump), TRUE))) "; PLINK 1.9 for clumping" else ""), paste0(M, "; scripts (02_scripts); reproducibility capsule")),
    r("9b", "AUTHOR", if (isTRUE(plan$prereg)) "Pre-registration was switched on for this study, but MR Studio does not register protocols. Give the registry, date and identifier." else "State that the analysis was not pre-registered, or give the registry and date."),
    r("10a", "FILLED", paste0("Instrument flow: ", each(function(cx) { L <- cx$ledger; if (is.null(L)) sprintf("%d simulated instruments", nrow(cx$kept)) else sprintf("%s%d candidates, %d not in the outcome, %s, %d allele mismatches dropped, %d used", if (!is.null(L$clump)) sprintf("%d below P, %d removed by clumping, %d not in the LD panel; ", as.integer(L$below_threshold), as.integer(L$clump$clumped), as.integer(L$clump$not_in_panel)) else "", as.integer(L$instruments), as.integer(L$missing_in_outcome), if (identical(L$strand, "forward")) sprintf("%d palindromic kept (both datasets declared forward strand)", as.integer(L$kept_palindromic %||% 0L)) else sprintf("%d palindromic dropped", as.integer(L$dropped_palindromic)), as.integer(L$dropped_mismatch), as.integer(L$kept)) }), "."), paste0(M, "; Check screen; instrument_selection.csv where clumped here")),
    r("10b", "PARTLY", "Sample sizes and ancestry of each source are in the methods; the variant-level data are in the supplement workbook. Add descriptive data from the source publications.", M),
    r("10c", "AUTHOR", "Heterogeneity between the studies inside each GWAS meta-analysis is reported by the source publications; cite or summarise it."),
    r("10d-i", "PARTLY", sprintf("Declared ancestry: %s. Justify that the exposure associations apply to the outcome population.", each(function(cx) sprintf("exposure %s, outcome %s", cx$src_x$ancestry, cx$src_y$ancestry)))),
    r("10d-ii", if (all(ovl != "not answered")) "PARTLY" else "AUTHOR", sprintf("Sample overlap: %s. MR Studio cannot calculate overlap from summary data; give the evidence for your statement.", paste(if (np > 1) paste0(vapply(ctxs, pn, ""), ": ", ovl) else ovl, collapse = "; ")), Lm),
    r("11a", "FILLED", "Variant-exposure and variant-outcome associations of every instrument (beta, standard error, ratio, F) are in the supplement workbook and shown in the scatter plot.", "Supplement workbook; scatter plot"),
    r("11b", "FILLED", "Every estimate with its confidence interval and P value; odds ratios for a binary outcome.", paste0(R, "; report table")),
    r("11c", "AUTHOR", "MR Studio does not convert to absolute risk. Add it if it helps interpretation."),
    r("11d", if (any(c("figs", "report") %in% items)) "FILLED" else "PARTLY", "Forest plot across methods and scatter plot with fitted lines for every pair, with legends filled with the numbers.", "Figures folder; report"),
    r("12a", "PARTLY", paste0("Reported: instrument strength, Cochran's Q", if (has("egger")) ", the MR-Egger intercept" else "", if (st) ", the Steiger test" else "", ". Interpret how plausible the assumptions are."), R),
    r("12b", "FILLED", "Cochran's Q for heterogeneity between variants, with its P value. No E-value is calculated.", R),
    r("13a", if (length(robust)) "FILLED" else "AUTHOR", if (length(robust)) sprintf("Robust methods: %s.", paste(ests[robust], collapse = ", ")) else "No robust method was run.", R),
    r("13b", if (length(detect)) "FILLED" else "AUTHOR", if (length(detect)) sprintf("Outlier detection: %s; flagged variants are reported and were not removed from the main analysis.", paste(ests[detect], collapse = ", ")) else "No further analysis was run in MR Studio.", R),
    r("13c", if (st) "FILLED" else "AUTHOR", if (st) "Steiger directionality test for every pair where it could run." else "No direction test was run. Turn on the Steiger direction test in the IVW card, or state that direction was not tested.", if (st) R else ""),
    r("13d", "AUTHOR", "Compare with observational estimates from the literature, where relevant."),
    r("13e", if ("loo" %in% used) "FILLED" else "AUTHOR", if ("loo" %in% used) "Leave-one-out and funnel plots (supplementary figures) and the leave-one-out table in the supplementary tables." else "No leave-one-out analysis was run. Turn on Leave-one-out in the IVW card, or describe other plots.", if ("loo" %in% used) R else ""),
    r("14", "PARTLY", "The results text summarises every estimate. Relate it to the objectives.", R),
    r("15", "PARTLY", "A limitations paragraph is drafted from the diagnostics (instrument strength, heterogeneity, overlap, direction, ancestry). Discuss the likely direction and size of bias.", Lm),
    r("16a", "AUTHOR", "Interpret cautiously and in the light of other evidence."),
    r("16b", "AUTHOR", "Discuss mechanisms and keep causal wording careful."),
    r("16c", "AUTHOR", "Discuss clinical or policy relevance."),
    r("17", "AUTHOR", sprintf("Discuss generalisability. The data are of %s ancestry.", paste(unique(unlist(lapply(ctxs, function(cx) c(cx$src_x$ancestry, cx$src_y$ancestry)))), collapse = " and "))),
    r("18", "AUTHOR", "Report funding, including that of the source studies and consortia."),
    r("19", "FILLED", "Data availability statement names the sources; scripts, analysis data and the reproducibility capsule let others re-run every number.", "Draft text, Data availability; capsule"),
    r("20", "AUTHOR", "Declare conflicts of interest."))
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  out <- merge(STROBE_ITEMS, out, by = "item", sort = FALSE)
  out <- out[match(STROBE_ITEMS$item, out$item), c("section", "item", "asks", "status", "provides", "where")]
  out$page <- ""
  attr(out, "simulated") <- sim
  out$provides <- sci_thr_text(out$provides)   # Round 109
  rownames(out) <- NULL; attr(out, "simulated") <- sim; out
}

#' Write the checklist as a Word document (one table, a short guide and the citation of the official checklist).
strobe_write <- function(path, rows, study_name) {
  tab <- data.frame(Item = rows$item, Section = rows$section, `What it asks (paraphrased)` = rows$asks, Status = rows$status,
                    `What MR Studio provides` = rows$provides, Where = rows$where, `Page/line` = rows$page, check.names = FALSE, stringsAsFactors = FALSE)
  n <- table(factor(rows$status, c("FILLED", "PARTLY", "AUTHOR")))
  docx_write(path, Filter(Negate(is.null), list(
    list(type = "title", text = sprintf("STROBE-MR checklist: %s", study_name)),
    list(type = "note", text = sprintf("Prepared by MR Studio on %s. %d items filled from the saved study, %d partly, %d for the authors. Check every entry and add page and line numbers before submission.", format(Sys.Date()), n[["FILLED"]], n[["PARTLY"]], n[["AUTHOR"]])),
    if (isTRUE(attr(rows, "simulated"))) list(type = "note", text = "SIMULATED DATA: at least one pair uses simulated demonstration traits. This checklist shows how the app fills it in; it is not for submission."),
    list(type = "p", text = "FILLED: MR Studio wrote it from the saved results and settings. PARTLY: MR Studio wrote part of it; the authors complete it. AUTHOR: only the authors can report it. The second column is a short paraphrase; use the official checklist for the exact wording."),
    list(type = "table", data = tab, widths = c(6, 10, 21, 10, 35, 11, 9)),
    list(type = "p", text = paste("Official checklist:", ref_text(REFS$skrivankova2021), ref_text(REFS$skrivankova2021ee))))))
  invisible(n)
}
