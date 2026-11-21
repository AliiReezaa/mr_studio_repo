# ---- Round 70: CAUSE, causal analysis using summary effect estimates (Morrison et al. 2020, Nat Genet 52:740-747) ---------------
# PMID 32451458. Package cause (GitHub jean997/cause, GPL >= 2), pinned. From the full summary statistics of both traits it
# estimates the nuisance parameters (the joint distribution of the effects and rho, the correlation of the summary statistics that
# sample overlap brings) from a random sample of up to 1,000,000 variants, selects variants with exposure P < 1e-3 clumped at
# r2 < 0.01 (locally with PLINK against the app's 1000 Genomes panel) and compares a sharing model (correlated pleiotropy only) with
# a causal model. It reports the posterior of the causal effect gamma and the P value that the causal model fits better.
# Steps as in the package vignette: gwas_merge -> est_cause_params -> ld_clump -> cause.
PIN_CAUSE <- "6257e6c66e6c1518c8750d1e4d2b4d6e42644fa5"
CAUSE_MEMORY_GB <- function(rows) 1.5 + 0.65 * max(rows, 0, na.rm = TRUE) / 1e6   # Round 83: measured 6.65 GB peak (LDL 2.4 M and CHD 8.5 M variants, macOS); was 1.2 per million (11.7 GB), an estimate
CAUSE_SAFE_SHARE <- 0.85   # Round 85: CAUSE starts when its peak fits in 85% of the free memory (half, the rule for estimates, kept a 16 GB Mac from
                           # running it on complete GWAS; its estimate is measured, and it runs at low priority on two threads)
CAUSE_THREADS <- 2L   # RcppParallel threads; kept low so the computer stays usable

#' CAUSE input for one trait: snp beta se ea oa (the columns gwas_merge() is told to read).
cause_prepare <- function(d) {
  led <- list(rows = nrow(d))
  if (!("snp" %in% names(d))) stop("CAUSE needs rsIDs; this dataset identifies variants by position only.")
  out <- data.frame(snp = d$snp, beta = d$beta, se = d$se, ea = toupper(d$ea), oa = toupper(d$oa), stringsAsFactors = FALSE)
  out <- out[!is.na(out$snp) & is.finite(out$beta) & is.finite(out$se) & out$se > 0, , drop = FALSE]
  out <- out[!duplicated(out$snp), , drop = FALSE]; led$kept <- nrow(out)
  attr(out, "ledger") <- led; out
}

cause_script <- function(folder, files, labels, bfile, plink, seed, pthr = 1e-3, r2 = 0.01, nsample = 1e6, threads = CAUSE_THREADS) {
  q <- function(s) encodeString(s, quote = "\"")
  c("# MR Studio: CAUSE, Mendelian randomization accounting for correlated and uncorrelated pleiotropy (Morrison et al. 2020, Nat Genet 52:740-747).",
    sprintf("# Package: cause commit %s (GitHub jean997/cause, GPL >= 2).", PIN_CAUSE),
    "# Run from the study folder. Steps as in the package vignette: gwas_merge (alleles aligned, ambiguous variants removed),",
    "# est_cause_params (nuisance parameters and rho from a random sample of variants), ld_clump (exposure P < 1e-3, r2 < 0.01,",
    "# local PLINK), cause (sharing and causal models compared by their expected log predictive density, ELPD).",
    sprintf("PAIR <- %s", q(folder)), sprintf("TRAITS <- c(%s)   # exposure = %s, outcome = %s", paste(vapply(files, q, ""), collapse = ", "), labels[1], labels[2]),
    sprintf("BFILE <- %s; PLINK <- %s   # 1000 Genomes %s panel and PLINK of MR Studio", q(bfile), q(plink), basename(bfile)),   # Round 81: the panel used
    sprintf("seed <- %dL; P_THRESHOLD <- %s; R2_THRESHOLD <- %s; N_SAMPLE <- %s; THREADS <- %dL", as.integer(seed), format(pthr), format(r2), format(nsample, scientific = FALSE), as.integer(threads)),
    "if (!requireNamespace(\"cause\", quietly = TRUE)) stop(\"Package cause is not installed (GitHub jean997/cause).\")",
    "if (utils::packageVersion(\"loo\") >= \"2.9.0\") stop(\"CAUSE needs loo 2.8.0, the version MR Studio uses: from loo 2.9 loo_compare() returns a table that cause cannot read (it stops at the very end). Install it with remotes::install_version('loo', '2.8.0').\")",   # Round 83
    "suppressPackageStartupMessages(library(cause)); RcppParallel::setThreadOptions(numThreads = THREADS)",
    "invisible(loadNamespace(\"ieugwasr\"))   # ieugwasr draws a random number when it loads; loaded before the seed, every run is the same",
    "rd <- function(f) as.data.frame(data.table::fread(f))",
    "X <- gwas_merge(rd(TRAITS[1]), rd(TRAITS[2]), snp_name_cols = c(\"snp\", \"snp\"), beta_hat_cols = c(\"beta\", \"beta\"), se_cols = c(\"se\", \"se\"),",
    "                A1_cols = c(\"ea\", \"ea\"), A2_cols = c(\"oa\", \"oa\"))",
    "set.seed(seed); varlist <- with(X, sample(snp, size = min(N_SAMPLE, nrow(X)), replace = FALSE))",
    "params <- est_cause_params(X, varlist)",
    "cl <- ieugwasr::ld_clump(dat = data.frame(rsid = X$snp, pval = X$p1, stringsAsFactors = FALSE), clump_r2 = R2_THRESHOLD, clump_p = P_THRESHOLD,",
    "                         plink_bin = PLINK, bfile = BFILE)",
    "top_vars <- cl$rsid",
    "set.seed(seed); res <- cause(X = X, variants = top_vars, param_ests = params)",
    "s <- summary(res, ci_size = 0.95); g <- s$quants[[2]]   # causal model: rows median, lower, upper; columns gamma, eta, q",
    "e <- res$elpd",
    "out <- data.frame(gamma = g[1, \"gamma\"], gamma_lower = g[2, \"gamma\"], gamma_upper = g[3, \"gamma\"], eta = g[1, \"eta\"], q = g[1, \"q\"],",
    "  sharing_eta = s$quants[[1]][1, \"eta\"], sharing_q = s$quants[[1]][1, \"q\"], z = s$z, p = s$p, instruments = length(top_vars), rho = params$rho,",
    "  delta_elpd = e$delta_elpd[e$model1 == \"sharing\" & e$model2 == \"causal\"], se_delta_elpd = e$se_delta_elpd[e$model1 == \"sharing\" & e$model2 == \"causal\"],",
    "  stringsAsFactors = FALSE)",
    "od <- file.path(\"04_across_pairs\", \"cause\"); dir.create(od, recursive = TRUE, showWarnings = FALSE)",
    "utils::write.csv(out, file.path(od, paste0(PAIR, \"_results.csv\")), row.names = FALSE)",
    "utils::write.csv(e, file.path(od, paste0(PAIR, \"_elpd.csv\")), row.names = FALSE)",
    "writeLines(top_vars, file.path(od, paste0(PAIR, \"_variants.txt\")))",
    "out")
}

cause_pair_prepare <- function(h, pair, root = studies_root(), accept = NULL) {
  if (isTRUE(pair$mv) || is.null(pair$folder)) stop("CAUSE is for a pair with one exposure.")
  ids <- full_gwas_ids(c(pair$ex, pair$out), root); if (!all(startsWith(ids, "lib:"))) stop("CAUSE needs the full summary statistics of both traits, imported into the Data library.")
  if (is.na(plink_bin())) stop("CAUSE clumps its variants with PLINK, which is not available.")
  ds <- lapply(ids, function(i) library_dataset(sub("^lib:", "", i), root))
  lr <- ldref_resolve(lapply(ds, `[[`, "recipe"), "cause", accept)   # Round 81: the panel of the GWAS's ancestry, or the analyst's accepted choice
  rows <- vapply(ds, function(x) as.numeric(x$recipe$rows_kept %||% nrow(x$data)), 0)
  if (exists("mem_check")) { mc <- mem_check(CAUSE_MEMORY_GB(rows), "CAUSE", share = CAUSE_SAFE_SHARE); if (!mc$ok) stop(mc$msg, call. = FALSE) }
  rel <- file.path("02_datasets", "cause"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE); files <- character(2)
  for (k in 1:2) { tr <- ds[[k]]$recipe$trait; f <- file.path(rel, paste0(.ldsc_slug(tr), ".tsv.gz")); files[k] <- f
    if (!file.exists(file.path(h$dir, f))) { m <- cause_prepare(ds[[k]]$data)
      if (nrow(m) < 200000L) stop(sprintf("\"%s\" has %s usable variants; CAUSE needs genome-wide summary statistics.", tr, format(nrow(m), big.mark = ",")))
      data.table::fwrite(m, file.path(h$dir, f), sep = "\t")
      jf <- sub("\\.tsv\\.gz$", "_preparation.json", f); atomic_write_text(file.path(h$dir, jf), jsonlite::toJSON(c(list(trait = tr, dataset_sha256 = ds[[k]]$recipe$original$sha256), attr(m, "ledger")), auto_unbox = TRUE, pretty = TRUE))
      manifest_register(h$dir, c(f, jf), role = "cause_input") } }
  rd <- file.path("04_across_pairs", "cause"); dir.create(file.path(h$dir, rd), recursive = TRUE, showWarnings = FALSE); sf <- file.path(rd, paste0(pair$folder, ".R"))
  writeLines(cause_script(pair$folder, files, c(pair$exName, pair$outName), normalizePath(lr$bfile, winslash = "/", mustWork = FALSE),
                          plink_bin(), .seed_of(paste("cause", pair$folder))), file.path(h$dir, sf))
  ldref_save(h, pair, lr)   # Round 81: the LD reference (and any acceptance) beside the result
  sf
}

.cause_status_file <- function(dir, folder) file.path(dir, "04_across_pairs", "cause", paste0(folder, "_status.json"))
cause_job_read <- function(h, pair) { f <- .cause_status_file(h$dir, pair$folder); if (!file.exists(f)) return(NULL); j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  if (!is.null(j) && identical(j$state, "running") && !is.null(j$pid) && !.pid_alive(j$pid)) j <- list(state = "stopped", msg = "The run stopped (the app or the computer was closed). Run it again.", time = j$time)
  j }
cause_pair_start <- function(h, pair, root = studies_root(), accept = NULL) {
  j <- cause_job_read(h, pair); if (!is.null(j) && identical(j$state, "running")) stop("CAUSE is already running for this pair.")
  sf <- cause_pair_prepare(h, pair, root, accept = accept)
  script_job_start(h$dir, sf, .cause_status_file(h$dir, pair$folder), file.path(h$dir, "04_across_pairs", "cause", paste0(pair$folder, "_job.log")))
  invisible(cause_job_read(h, pair))
}
cause_pair_read <- function(h, pair) { f <- file.path(h$dir, "04_across_pairs", "cause", paste0(pair$folder, "_results.csv"))
  if (!file.exists(f)) return(NULL); r <- utils::read.csv(f, stringsAsFactors = FALSE); attr(r, "labels") <- c(x = pair$exName, y = pair$outName); r }

CAUSE_METHODS_TEXT <- "The causal effect was also assessed with CAUSE [Morrison 2020], which distinguishes a causal effect from correlated pleiotropy (a shared heritable factor) using variants with exposure P < 1e-3 (clumped at r2 < 0.01 against the 1000 Genomes European panel), with nuisance parameters, including the correlation that sample overlap brings, estimated from a random sample of up to one million variants. The causal model was compared with the sharing model by the expected log predictive density."
cause_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); r <- res[1, ]
  sprintf("CAUSE (%d variants): the causal model %s the sharing model (P = %s); %s on %s gamma %.2f (95%% credible interval %.2f to %.2f).",
          r$instruments, if (r$p < 0.05) "fitted better than" else "did not fit better than", format(signif(r$p, 2)), l[["x"]], l[["y"]], r$gamma, r$gamma_lower, r$gamma_upper)
}
