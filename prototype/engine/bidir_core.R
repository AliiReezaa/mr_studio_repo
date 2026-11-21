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

# ---- Round 72c: the reverse direction of a bidirectional study from OpenGWAS downloads ----------------------------------------------
# A study picks exposure A (instruments of A) and outcome B (B looked up at A's instruments). The reverse direction needs B's own
# instruments and A looked up at them, which these two downloads do not hold. When both come from OpenGWAS, the reverse pair is
# prepared from the same datasets with the same settings: B is downloaded as an exposure (P threshold, clumping, ancestry as for A)
# and A as an outcome for it (proxy settings as for B). Matching downloads already in the Data library are reused.
.lib_recipes <- function(root = studies_root()) {
  d <- list.dirs(library_dir(root), recursive = FALSE, full.names = TRUE); d <- d[file.exists(file.path(d, "recipe.json"))]
  r <- lapply(d, function(x) tryCatch(jsonlite::fromJSON(file.path(x, "recipe.json"), simplifyVector = FALSE), error = function(e) NULL)); names(r) <- basename(d); Filter(Negate(is.null), r)
}
.same_params <- function(a, b, keys) all(vapply(keys, function(k) identical(as.character(a[[k]] %||% ""), as.character(b[[k]] %||% "")), TRUE))
.free_name <- function(nm, root) { cand <- nm; i <- 2L; while (dir.exists(file.path(library_dir(root), cand))) { cand <- sprintf("%s %d", nm, i); i <- i + 1L }; cand }

#' Library ids for the reverse direction of the pair (exposure ex_id, outcome out_id). Returns list(ex, out, exName, outName, note).
bidir_reverse <- function(ex_id, out_id, jwt = NULL, offline = FALSE, root = studies_root(), get_exposure = ogw_get_exposure, get_outcome = ogw_get_outcome) {
  plain <- list(ex = out_id, out = ex_id, note = "")
  if (!startsWith(ex_id, "lib:") || !startsWith(out_id, "lib:")) return(plain)   # simulated traits: the same datasets serve both directions
  e <- library_dataset(sub("^lib:", "", ex_id), root)$recipe; o <- library_dataset(sub("^lib:", "", out_id), root)$recipe
  if (isTRUE(o$usable_as_exposure) && is.null(e$source$for_exposure) && isTRUE(e$usable_as_outcome)) return(plain)   # full summary statistics both ways
  es <- e$source; os <- o$source
  if (!identical(es$type, "opengwas") || !identical(os$type, "opengwas") || is.null(es$id) || is.null(os$id))
    stop(sprintf("The reverse direction needs \"%s\" as an exposure and \"%s\" as an outcome for it. Download both from OpenGWAS, or import the full summary statistics of both traits.", o$trait, e$trait), call. = FALSE)
  .say(sprintf("Reverse direction %s → %s: preparing %s as an exposure", o$trait, e$trait, o$trait), "Bidirectional MR tests each direction with its own instruments: the variants of the outcome trait become the instruments of the reverse direction.")
  rx <- ogw_ensure_exposure(os, o$trait, es$params, jwt, offline, root, get_exposure)
  ry <- ogw_ensure_outcome(es, e$trait, rx$name, os$params, jwt, offline, root, get_outcome)
  nx <- rx$name; ny <- ry$name; note <- c(rx$note, ry$note); note <- note[nzchar(note)]
  list(ex = paste0("lib:", nx), out = paste0("lib:", ny), exName = nx, outName = ny, note = paste(note, collapse = "; "))
}

#' The OpenGWAS dataset `src` (source block of a download) as an exposure with the instrument settings `px`; reused when present.
ogw_ensure_exposure <- function(src, label, px, jwt, offline, root, get_exposure = ogw_get_exposure) {
  kx <- c("pval", "clump", "r2", "kb", "pop", "where")
  hit <- Filter(function(r) identical(r$source$type, "opengwas") && identical(r$source$id, src$id) && isTRUE(r$usable_as_exposure) && .same_params(r$source$params, px, kx), .lib_recipes(root))
  if (length(hit)) { .say(sprintf("%s as an exposure is already in the Data library (same settings), so it is reused", src$trait %||% label)); return(list(name = names(hit)[1], note = "")) }
  .say(sprintf("Downloading the instruments of %s (%s) from OpenGWAS: P < %s, clumped", src$trait %||% label, src$id, format(as.numeric(px$pval %||% 5e-8))), "Instruments are the genome-wide significant variants of the trait, one per region.")
  nx <- .free_name(sprintf("%s (%s)", src$trait %||% label, src$id), root)
  nx <- get_exposure(src$id, jwt, as.numeric(px$pval %||% 5e-8), as.numeric(px$r2 %||% 0.001), as.numeric(px$kb %||% 10000), px$pop, nx, offline, root = root, where = px$where %||% "server")$trait
  list(name = nx, note = sprintf("downloaded \"%s\" as an exposure", nx))
}
#' The OpenGWAS dataset `src` looked up at the instruments of the library exposure `exposure_name`, with the proxy settings `py`.
ogw_ensure_outcome <- function(src, label, exposure_name, py, jwt, offline, root, get_outcome = ogw_get_outcome) {
  ky <- c("proxies", "rsq", "maf_threshold")
  hit <- Filter(function(r) identical(r$source$type, "opengwas") && identical(r$source$id, src$id) && identical(r$source$for_exposure, exposure_name) && .same_params(r$source$params, py, ky), .lib_recipes(root))
  if (length(hit)) { .say(sprintf("%s looked up for %s is already in the Data library, so it is reused", src$trait %||% label, exposure_name)); return(list(name = names(hit)[1], note = "")) }
  .say(sprintf("Looking up the instruments of %s in %s (%s) on OpenGWAS, with proxies where one is missing", exposure_name, src$trait %||% label, src$id), "Each instrument's effect on the outcome is needed; a proxy is a nearby variant in strong linkage (r² of at least 0.8) used when the instrument itself is missing.")
  ny <- .free_name(sprintf("%s for %s", src$trait %||% label, exposure_name), root)
  ny <- get_outcome(src$id, exposure_name, jwt, !identical(py$proxies, FALSE), as.numeric(py$rsq %||% 0.8), as.numeric(py$maf_threshold %||% 0.3), ny, offline, root = root)$trait
  list(name = ny, note = sprintf("downloaded \"%s\" as an outcome for %s", ny, exposure_name))
}
