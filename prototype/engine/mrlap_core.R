# ---- Round 68: MRlap, IVW corrected for sample overlap, weak instruments and winner's curse (Mounier and Kutalik 2023) -----------
# Genet Epidemiol 47:314-331 (PMID 37036286). Package MRlap (GitHub n-mounier/MRlap, GPL-2), pinned. It needs the full summary
# statistics of both traits (it selects and prunes the instruments itself, by distance, and runs cross-trait LDSC to estimate the
# overlap), and the LDSC authors' LD scores and HapMap3 list (ldsc_core.R). The standard error of the corrected effect comes from a
# simulation (random draws), so the script sets the seed and replays exactly.
PIN_MRLAP <- "660f026864f8bfbbad5a8206bdff7d58f5d5d05b"
MRLAP_MEMORY_GB <- function(rows) 1.0 + 0.7 * max(rows, 0, na.rm = TRUE) / 1e6   # measured 1.8 GB on the package example (1.15 million variants)

#' MRlap's input for one trait: rsid chr pos a1 a2 beta se N. Positions come from the dataset, or from the 1000 Genomes panel's .bim
#' (`posmap`, rsID -> chr/pos) when the dataset has none; variants without a position are left out.
mrlap_prepare <- function(d, posmap = NULL) {
  led <- list(rows = nrow(d))
  if (!("snp" %in% names(d))) stop("MRlap needs rsIDs; this dataset identifies variants by position only.")
  if (!("n" %in% names(d)) || all(is.na(d$n))) stop("MRlap needs the sample size; give it when importing the dataset.")
  chr <- if ("chr" %in% names(d)) suppressWarnings(as.integer(sub("^chr", "", d$chr, ignore.case = TRUE))) else rep(NA_integer_, nrow(d))
  pos <- if ("pos" %in% names(d)) suppressWarnings(as.numeric(d$pos)) else rep(NA_real_, nrow(d))
  miss <- is.na(chr) | is.na(pos)
  if (any(miss) && !is.null(posmap)) { i <- match(d$snp[miss], posmap$snp); chr[miss] <- posmap$chr[i]; pos[miss] <- posmap$pos[i] }
  out <- data.frame(rsid = d$snp, chr = chr, pos = pos, a1 = d$ea, a2 = d$oa, beta = d$beta, se = d$se, N = d$n, stringsAsFactors = FALSE)
  k <- is.finite(out$chr) & is.finite(out$pos) & is.finite(out$beta) & is.finite(out$se) & out$se > 0 & is.finite(out$N)
  led$no_position <- sum(!(is.finite(out$chr) & is.finite(out$pos))); out <- out[k, , drop = FALSE]; led$kept <- nrow(out)
  attr(out, "ledger") <- led; out
}

#' Script for one pair. `pthr` is MRlap's instrument threshold (its default 5e-8).
mrlap_script <- function(folder, files, labels, ld, hm3, seed, pthr = 5e-8) {
  q <- function(s) encodeString(s, quote = "\"")
  c("# MR Studio: MRlap, IVW corrected for sample overlap, weak-instrument bias and winner's curse (Mounier and Kutalik 2023, Genet Epidemiol 47:314).",
    sprintf("# Package: MRlap commit %s (GitHub n-mounier/MRlap); cross-trait LDSC with GenomicSEM commit %s.", PIN_MRLAP, PIN_GENOMICSEM),
    "# Run from the study folder. The summary statistics in 02_datasets/mrlap were prepared by MR Studio from the Data library files.",
    sprintf("PAIR <- %s", q(folder)), sprintf("TRAITS <- c(%s)   # exposure = %s, outcome = %s", paste(vapply(files, q, ""), collapse = ", "), labels[1], labels[2]),
    sprintf("LD_SCORES <- %s; HM3 <- %s", q(ld), q(hm3)), sprintf("seed <- %dL; MR_THRESHOLD <- %s", as.integer(seed), format(pthr)),
    "if (!requireNamespace(\"MRlap\", quietly = TRUE)) stop(\"Package MRlap is not installed (GitHub n-mounier/MRlap).\")",
    "invisible(loadNamespace(\"ieugwasr\"))   # it draws a random number when it first loads; loaded before the seed, every run is the same",
    "set.seed(seed)   # the standard error of the corrected effect is simulated",
    "od <- file.path(\"04_across_pairs\", \"mrlap\"); dir.create(od, recursive = TRUE, showWarnings = FALSE)",
    "invisible(utils::capture.output(r <- MRlap::MRlap(exposure = TRAITS[1], exposure_name = \"exposure\", outcome = TRAITS[2], outcome_name = \"outcome\",",
    "  ld = LD_SCORES, hm3 = HM3, MR_threshold = MR_THRESHOLD, verbose = FALSE)))",
    "m <- r$MRcorrection; l <- r$LDSC",
    "out <- data.frame(quantity = c(\"observed_effect\", \"corrected_effect\", \"difference\"), estimate = c(m$observed_effect, m$corrected_effect, NA),",
    "  se = c(m$observed_effect_se, m$corrected_effect_se, NA), p = c(m$observed_effect_p, m$corrected_effect_p, m$p_difference),",
    "  z = c(m$observed_effect / m$observed_effect_se, m$corrected_effect / m$corrected_effect_se, m$test_difference), stringsAsFactors = FALSE)",
    "out$ci_lower <- out$estimate - 1.96 * out$se; out$ci_upper <- out$estimate + 1.96 * out$se; out$instruments <- m$m_IVs",
    "out$rg <- l$rg; out$int_crosstrait <- l$int_crosstrait; out$int_crosstrait_se <- l$int_crosstrait_se; out$h2_exposure <- l$h2_exp; out$h2_outcome <- l$h2_out",
    "utils::write.csv(out, file.path(od, paste0(PAIR, \"_results.csv\")), row.names = FALSE)",
    "writeLines(m$IVs, file.path(od, paste0(PAIR, \"_instruments.txt\")))",
    "out")
}

#' Prepare both datasets and write the script (fast); the run itself is started in the background by mrlap_pair_start().
mrlap_pair_prepare <- function(h, pair, root = studies_root(), ldsc = ldsc_ref_dir()) {
  if (isTRUE(pair$mv) || is.null(pair$folder)) stop("MRlap is for a pair with one exposure.")
  ids <- full_gwas_ids(c(pair$ex, pair$out), root); if (!all(startsWith(ids, "lib:"))) stop("MRlap needs the full summary statistics of both traits, imported into the Data library.")
  if (!ldsc_ref_status(ldsc)$official) stop("MRlap uses the European LD scores of LDSC; install them first (section above).")
  ds <- lapply(ids, function(i) library_dataset(sub("^lib:", "", i), root))
  if (!all(vapply(ds, function(x) toupper(x$recipe$facts$ancestry %||% "EUR"), "") == "EUR")) stop("MRlap uses European LD scores; both GWAS must be of European ancestry.")
  rows <- vapply(ds, function(x) as.numeric(x$recipe$rows_kept %||% nrow(x$data)), 0)
  if (exists("mem_check")) { mc <- mem_check(MRLAP_MEMORY_GB(rows), "MRlap"); if (!mc$ok) stop(mc$msg, call. = FALSE) }
  rel <- file.path("02_datasets", "mrlap"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE); files <- character(2); posmap <- NULL
  for (k in 1:2) { tr <- ds[[k]]$recipe$trait; f <- file.path(rel, paste0(.ldsc_slug(tr), ".tsv.gz")); files[k] <- f
    if (!file.exists(file.path(h$dir, f))) {
      d <- ds[[k]]$data
      if (is.null(posmap) && (!all(c("chr", "pos") %in% names(d)) || anyNA(d$pos)) && file.exists(file.path(ld_panel_dir(), "EUR.bim"))) {
        b <- data.table::fread(file.path(ld_panel_dir(), "EUR.bim"), select = c(1, 2, 4), col.names = c("chr", "snp", "pos")); posmap <- as.data.frame(b) }
      m <- mrlap_prepare(d, posmap)
      if (nrow(m) < 200000L) stop(sprintf("\"%s\" has %s usable variants; MRlap needs genome-wide summary statistics.", tr, format(nrow(m), big.mark = ",")))
      data.table::fwrite(m, file.path(h$dir, f), sep = "\t")
      jf <- sub("\\.tsv\\.gz$", "_preparation.json", f); atomic_write_text(file.path(h$dir, jf), jsonlite::toJSON(c(list(trait = tr, dataset_sha256 = ds[[k]]$recipe$original$sha256), attr(m, "ledger")), auto_unbox = TRUE, pretty = TRUE))
      manifest_register(h$dir, c(f, jf), role = "mrlap_input") } }
  rd <- file.path("04_across_pairs", "mrlap"); dir.create(file.path(h$dir, rd), recursive = TRUE, showWarnings = FALSE); sf <- file.path(rd, paste0(pair$folder, ".R"))
  writeLines(mrlap_script(pair$folder, files, c(pair$exName, pair$outName), normalizePath(ldsc_ld_folder("official", ldsc), winslash = "/"),
                          normalizePath(file.path(ldsc, "w_hm3.snplist"), winslash = "/"), .seed_of(paste("mrlap", pair$folder)), pair$pthr %||% 5e-8), file.path(h$dir, sf))
  sf
}

#' Run a study script in the background from the study folder; its state goes to `status` (running, done or error, with the process).
script_job_start <- function(study_dir, script, status, log, manifest = character(), role = "script") {
  code <- sprintf(paste0("sf <- %s; st <- function(s, m = '') writeLines(jsonlite::toJSON(list(state = s, msg = m, pid = Sys.getpid(), time = format(Sys.time(), '%%Y-%%m-%%d %%H:%%M:%%S')), auto_unbox = TRUE), sf)",
    "; invisible(tools::psnice(value = 15L)); st('running'); setwd(%s); tryCatch({ sys.source(%s, envir = new.env(parent = globalenv())); st('done') }, error = function(e) st('error', conditionMessage(e)))"),
    deparse(status), deparse(study_dir), deparse(script))
  writeLines(jsonlite::toJSON(list(state = "running", msg = "", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), status)
  rscript_background(c("--vanilla", "-e", shQuote(code)), log = log); invisible(TRUE)
}

.mrlap_status_file <- function(dir, folder) file.path(dir, "04_across_pairs", "mrlap", paste0(folder, "_status.json"))
mrlap_job_read <- function(h, pair) { f <- .mrlap_status_file(h$dir, pair$folder); if (!file.exists(f)) return(NULL); j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  if (!is.null(j) && identical(j$state, "running") && !is.null(j$pid) && !.pid_alive(j$pid)) j <- list(state = "stopped", msg = "The run stopped (the app or the computer was closed). Run it again.", time = j$time)
  j }
mrlap_pair_start <- function(h, pair, root = studies_root()) {
  j <- mrlap_job_read(h, pair); if (!is.null(j) && identical(j$state, "running")) stop("MRlap is already running for this pair.")
  sf <- mrlap_pair_prepare(h, pair, root)
  script_job_start(h$dir, sf, .mrlap_status_file(h$dir, pair$folder), file.path(h$dir, "04_across_pairs", "mrlap", paste0(pair$folder, "_job.log")))
  invisible(mrlap_job_read(h, pair))
}
mrlap_pair_read <- function(h, pair) { f <- file.path(h$dir, "04_across_pairs", "mrlap", paste0(pair$folder, "_results.csv"))
  if (!file.exists(f)) return(NULL); r <- utils::read.csv(f, stringsAsFactors = FALSE); attr(r, "labels") <- c(x = pair$exName, y = pair$outName); r }

MRLAP_METHODS_TEXT <- "IVW estimates were also corrected for weak-instrument bias and winner's curse, allowing for sample overlap estimated by cross-trait LD score regression, with MRlap [Mounier 2023], which selects its instruments from the full summary statistics (P < 5e-8, distance pruning 500 kb); the corrected and observed effects were compared with its difference test."
mrlap_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); g <- function(q) res[res$quantity == q, ][1, ]; o <- g("observed_effect"); k <- g("corrected_effect"); d <- g("difference")
  f <- function(r) sprintf("%.3f (95%% CI %.3f to %.3f)", r$estimate, r$ci_lower, r$ci_upper)
  sprintf("MRlap (%d instruments; cross-trait intercept %.3f): the IVW effect of %s on %s was %s and %s after correction for weak instruments, winner's curse and sample overlap (difference P = %s)%s.",
          o$instruments, o$int_crosstrait, l[["x"]], l[["y"]], f(o), f(k), format(signif(d$p, 2)), if (is.finite(d$p) && d$p < 0.05) "; the corrected effect should be preferred" else "; the observed IVW effect is not materially biased")
}
