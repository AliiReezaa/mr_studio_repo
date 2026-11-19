# ---- Round 63: bidirectional MR with constrained maximum likelihood (Xue and Pan 2022, PLoS Genet 18:e1010205) ----------------------
# BiDirMRcML and BiDirCDcML from BiDirectCausal (GitHub xue-hr/BiDirectCausal, MIT), with MRcML and MRCD (same authors), pinned.
# They need the variants of both directions at once: the instruments of X with their effects on Y, and the instruments of Y with their
# effects on X. In a bidirectional study these are the forward and the reverse pair, so the analysis is an across-pairs script
# (04_across_pairs/bidirectional_cml/<forward pair>.R) that reads both pairs' analysis data and runs both methods.
PIN_BIDIRECTCAUSAL <- "b463ff10621a700ad72d28337135e1a799229812"; PIN_MRCML <- "58a3481976093175cbb82a5e6bf06c10e3c576f2"; PIN_MRCD <- "11ca09599063ac236ba66c0307aa1186aebb0cad"

bidir_cml_script <- function(fwd, rev, xl, yl, seed, pthr) {
  q <- function(s) encodeString(s, quote = "\"")
  c("# MR Studio: bidirectional MR with constrained maximum likelihood (Xue and Pan 2022, PLoS Genet 18:e1010205).",
    sprintf("# Packages: BiDirectCausal commit %s; MRcML commit %s; MRCD commit %s (GitHub xue-hr).", PIN_BIDIRECTCAUSAL, PIN_MRCML, PIN_MRCD),
    "# Run from the study folder. Reads the analysis data of the forward pair (instruments of X, effects on Y) and of the reverse pair",
    "# (instruments of Y, effects on X), takes their union and estimates the effect in both directions at once.",
    sprintf("FORWARD <- %s   # %s on %s", q(fwd), xl, yl), sprintf("REVERSE <- %s   # %s on %s", q(rev), yl, xl),
    sprintf("SIG_CUTOFF <- %s   # P-value threshold that defines the instruments of each trait (as in the pairs)", format(pthr)), sprintf("seed <- %dL", as.integer(seed)),
    "for (p in c(\"BiDirectCausal\", \"MRcML\", \"MRCD\")) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf(\"Package %s is not installed (GitHub xue-hr/%s).\", p, p))",
    "suppressPackageStartupMessages({ library(MRcML); library(MRCD) })   # BiDirectCausal calls their functions without importing them",
    "rd <- function(f) { x <- utils::read.csv(Sys.glob(file.path(\"03_pairs\", f, \"01_data\", \"analysis_data_*.csv\"))[1], stringsAsFactors = FALSE); x[x$mr_keep %in% c(TRUE, \"TRUE\"), , drop = FALSE] }",
    "a <- rd(FORWARD); b <- rd(REVERSE)",
    "u <- rbind(data.frame(id = a$variant_id, bX = a$beta_x, sX = a$se_x, bY = a$beta_y, sY = a$se_y, nX = a$n_x, nY = a$n_y),",
    "           data.frame(id = b$variant_id, bX = b$beta_y, sX = b$se_y, bY = b$beta_x, sY = b$se_x, nX = b$n_y, nY = b$n_x))",
    "u <- u[!duplicated(u$id) & stats::complete.cases(u[, c(\"bX\", \"sX\", \"bY\", \"sY\")]), ]",
    "nX <- stats::median(u$nX, na.rm = TRUE); nY <- stats::median(u$nY, na.rm = TRUE)",
    "RNGkind(\"Mersenne-Twister\", \"Inversion\", \"Rejection\"); set.seed(seed)",
    "f1 <- BiDirectCausal::BiDirMRcML(u$bX, u$bY, u$sX, u$sY, nX, nY, sig.cutoff = SIG_CUTOFF, random.seed = seed)",
    "set.seed(seed); f2 <- BiDirectCausal::BiDirCDcML(u$bX, u$bY, u$sX, u$sY, nX, nY, sig.cutoff = SIG_CUTOFF, random.seed = seed)",
    "row <- function(m, f, v) do.call(rbind, lapply(c(\"XtoY\", \"YtoX\"), function(d) { b <- f[[sprintf(\"%s.est.%s\", d, v)]]; se <- f[[sprintf(\"%s.se.%s\", d, v)]]",
    "  data.frame(method = m, version = v, direction = d, estimate = b, se = se, ci_lower = b - 1.96 * se, ci_upper = b + 1.96 * se, p = 2 * stats::pnorm(-abs(b / se)), stringsAsFactors = FALSE) }))",
    "out <- do.call(rbind, c(lapply(c(\"S.DP\", \"S\", \"NoS.DP\", \"NoS\"), function(v) row(\"BiDirMRcML\", f1, v)), lapply(c(\"S.DP\", \"S\", \"NoS.DP\", \"NoS\"), function(v) row(\"BiDirCDcML\", f2, v))))",
    "out$variants <- nrow(u)",
    "dir.create(file.path(\"04_across_pairs\", \"bidirectional_cml\"), recursive = TRUE, showWarnings = FALSE)",
    "utils::write.csv(out, file.path(\"04_across_pairs\", \"bidirectional_cml\", paste0(FORWARD, \"_results.csv\")), row.names = FALSE)",
    "out")
}

#' Run the script for every forward/reverse couple of a study. Returns a list of results tables (one per couple), or an empty list.
bidir_cml_run <- function(h, pairs, pthr = 5e-8) {
  pr <- Filter(function(p) !is.null(p$folder) && !isTRUE(p$mv), pairs); out <- list()
  for (p in Filter(function(p) !isTRUE(p$rev), pr)) {
    r <- Filter(function(q) isTRUE(q$rev) && identical(q$ex, p$out) && identical(q$out, p$ex), pr); if (!length(r)) next
    rel <- file.path("04_across_pairs", "bidirectional_cml"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE)
    sf <- file.path(h$dir, rel, paste0(p$folder, ".R")); writeLines(bidir_cml_script(p$folder, r[[1]]$folder, p$exName, p$outName, .seed_of(paste("bidircml", p$folder)), p$pthr %||% pthr), sf)
    old <- setwd(h$dir); res <- tryCatch({ sys.source(sf, envir = new.env(parent = globalenv()), keep.source = FALSE)
      utils::read.csv(file.path(rel, paste0(p$folder, "_results.csv")), stringsAsFactors = FALSE) }, error = function(e) e, finally = setwd(old))
    if (inherits(res, "error")) { out[[p$folder]] <- structure(list(error = conditionMessage(res)), class = "bidir_cml_error"); next }
    manifest_register(h$dir, file.path(rel, paste0(p$folder, c(".R", "_results.csv"))), role = "bidirectional_cml")
    attr(res, "labels") <- c(x = p$exName, y = p$outName); out[[p$folder]] <- res
  }
  out
}

BIDIR_CML_METHODS_TEXT <- "Effects in both directions were also estimated jointly with bidirectional constrained maximum likelihood (BiDirMRcML) and with its causal-direction version on standardised effects (BiDirCDcML, per standard deviation), which use the instruments of both traits together and are robust to correlated and uncorrelated pleiotropy; the versions with instrument screening and data perturbation were the main results [Xue 2022]."
bidir_cml_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); g <- function(m, d) res[res$method == m & res$version == "S.DP" & res$direction == d, ][1, ]
  f <- function(r) sprintf("%.3f (95%% CI %.3f to %.3f; P = %s)", r$estimate, r$ci_lower, r$ci_upper, format(signif(r$p, 2)))
  sprintf("Bidirectional constrained maximum likelihood (%d variants): %s on %s %s, and %s on %s %s (BiDirMRcML); BiDirCDcML, which works on standardised effects (per SD of one trait on SD of the other), gave %s and %s.",
          res$variants[1], l[["x"]], l[["y"]], f(g("BiDirMRcML", "XtoY")), l[["y"]], l[["x"]], f(g("BiDirMRcML", "YtoX")), f(g("BiDirCDcML", "XtoY")), f(g("BiDirCDcML", "YtoX")))
}
