# ---- Round 69: MR-APSS, MR accounting for pleiotropy and sample structure (Hu et al. 2022, PNAS 119:e2106858119) ---------------
# PMID 35787050. Package MRAPSS (GitHub YangLabHKUST/MR-APSS, MIT), pinned. From the full summary statistics of both traits it
# estimates the background (sample structure: the LDSC intercepts C; polygenic pleiotropy: Omega), selects instruments at a looser
# threshold (P < 5e-5, clumped locally with PLINK against the app's 1000 Genomes panel) and corrects for the selection bias that this
# threshold brings. Everything runs on this computer: format_data -> est_paras -> clump -> MRAPSS.
PIN_MRAPSS <- "85c18eaef8ccbc1b330d31db0b3984c2fd0f528f"
MRAPSS_MEMORY_GB <- function(rows) 1.0 + 0.3 * max(rows, 0, na.rm = TRUE) / 1e6   # Round 100: measured 3.0 GB peak in the shipped runtime (LDL 2.4 M and CHD 8.5 M variants, macOS; preparation 2.3 GB); was 0.6 per million (6.1 GB), an estimate
MRAPSS_SAFE_SHARE <- 0.85   # Round 100: measured on real data, like CAUSE and MRlap

#' MR-APSS input for one trait, in the column names format_data() reads by default: SNP A1 A2 b se n [freq].
mrapss_prepare <- function(d) {
  led <- list(rows = nrow(d))
  if (!("snp" %in% names(d))) stop("MR-APSS needs rsIDs; this dataset identifies variants by position only.")
  if (!("n" %in% names(d)) || all(is.na(d$n))) stop("MR-APSS needs the sample size; give it when importing the dataset.")
  out <- data.frame(SNP = d$snp, A1 = d$ea, A2 = d$oa, b = d$beta, se = d$se, n = d$n, stringsAsFactors = FALSE)
  if ("eaf" %in% names(d) && any(is.finite(d$eaf))) out$freq <- d$eaf
  out <- out[is.finite(out$b) & is.finite(out$se) & out$se > 0 & is.finite(out$n), , drop = FALSE]; led$kept <- nrow(out)
  attr(out, "ledger") <- led; out
}

mrapss_script <- function(folder, files, labels, ld, bfile, plink, seed, pthr = 5e-5) {
  q <- function(s) encodeString(s, quote = "\"")
  c("# MR Studio: MR-APSS, Mendelian randomization accounting for pleiotropy and sample structure (Hu et al. 2022, PNAS 119:e2106858119).",
    sprintf("# Package: MRAPSS commit %s (GitHub YangLabHKUST/MR-APSS, MIT).", PIN_MRAPSS),
    "# Run from the study folder. Steps as in the package: format_data (HapMap3 variants, MHC removed), est_paras (background: LDSC",
    "# intercepts C and polygenic pleiotropy Omega), clump (instruments at P < 5e-5, local PLINK), MRAPSS (with selection-bias correction).",
    sprintf("PAIR <- %s", q(folder)), sprintf("TRAITS <- c(%s)   # exposure = %s, outcome = %s", paste(vapply(files, q, ""), collapse = ", "), labels[1], labels[2]),
    sprintf("LD_SCORES <- %s   # %s", q(paste0(sub("/+$", "", ld), "/")), if (basename(sub("/+$", "", ld)) == "eur_w_ld_chr") "LDSC authors' European LD scores" else sprintf("LD scores computed by MR Studio from the 1000 Genomes %s panel", basename(sub("/+$", "", ld)))),
    sprintf("BFILE <- %s; PLINK <- %s   # 1000 Genomes %s panel and PLINK of MR Studio", q(bfile), q(plink), basename(bfile)),   # Round 81: the panel used
    sprintf("seed <- %dL; IV_THRESHOLD <- %s", as.integer(seed), format(pthr)),
    "if (!requireNamespace(\"MRAPSS\", quietly = TRUE)) stop(\"Package MRAPSS is not installed (GitHub YangLabHKUST/MR-APSS).\")",
    "suppressPackageStartupMessages(library(MRAPSS))   # its functions use the package's own data sets (w_hm3.snplist, MHC.SNPs)",
    "invisible(loadNamespace(\"ieugwasr\")); set.seed(seed)   # ieugwasr draws a random number when it loads; loaded first, every run is the same",
    "rd <- function(f) as.data.frame(data.table::fread(f))",
    "invisible(utils::capture.output({ x <- format_data(rd(TRAITS[1])); y <- format_data(rd(TRAITS[2]))",
    "  pa <- est_paras(dat1 = x, dat2 = y, trait1.name = \"exposure\", trait2.name = \"outcome\", ldscore.dir = LD_SCORES)",
    "  iv <- clump(pa$dat, IV.Threshold = IV_THRESHOLD, SNP_col = \"SNP\", pval_col = \"pval.exp\", clump_kb = 1000, clump_r2 = 0.001, bfile = BFILE, plink_bin = PLINK)",
    "  r <- MRAPSS(iv, exposure = \"exposure\", outcome = \"outcome\", C = pa$C, Omega = pa$Omega, Cor.SelectionBias = TRUE) }))",
    "z <- stats::qnorm(0.975); p <- as.numeric(r$pvalue)",
    "out <- data.frame(estimate = r$beta, se = r$beta.se, ci_lower = r$beta - z * r$beta.se, ci_upper = r$beta + z * r$beta.se, p = p, instruments = nrow(iv),",
    "  pi0 = r$pi0, sigma_sq = r$sigma.sq, tau_sq = r$tau.sq, iv_threshold_corrected = r$Threshold, c12 = pa$C[1, 2], stringsAsFactors = FALSE)",
    "od <- file.path(\"04_across_pairs\", \"mrapss\"); dir.create(od, recursive = TRUE, showWarnings = FALSE)",
    "utils::write.csv(out, file.path(od, paste0(PAIR, \"_results.csv\")), row.names = FALSE)",
    "utils::write.csv(iv, file.path(od, paste0(PAIR, \"_instruments.csv\")), row.names = FALSE)",
    "out")
}

mrapss_pair_prepare <- function(h, pair, root = studies_root(), ldsc = ldsc_ref_dir(), accept = NULL) {
  if (isTRUE(pair$mv) || is.null(pair$folder)) stop("MR-APSS is for a pair with one exposure.")
  ids <- full_gwas_ids(c(pair$ex, pair$out), root); if (!all(startsWith(ids, "lib:"))) stop("MR-APSS needs the full summary statistics of both traits, imported into the Data library.")
  if (is.na(plink_bin())) stop("MR-APSS clumps its instruments with PLINK, which is not available.")
  ds <- lapply(ids, function(i) library_dataset(sub("^lib:", "", i), root))
  lr <- ldref_resolve(lapply(ds, `[[`, "recipe"), "mrapss", accept, ldsc)   # Round 81: LD scores and panel of the GWAS's ancestry, or the analyst's accepted choice
  rows <- vapply(ds, function(x) as.numeric(x$recipe$rows_kept %||% nrow(x$data)), 0)
  if (exists("mem_check")) { mc <- mem_check(MRAPSS_MEMORY_GB(rows), "MR-APSS", share = MRAPSS_SAFE_SHARE); if (!mc$ok) stop(mc$msg, call. = FALSE) }
  rel <- file.path("02_datasets", "mrapss"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE); files <- character(2)
  for (k in 1:2) { tr <- ds[[k]]$recipe$trait; f <- file.path(rel, paste0(.ldsc_slug(tr), ".tsv.gz")); files[k] <- f
    if (!file.exists(file.path(h$dir, f))) { m <- mrapss_prepare(ds[[k]]$data)
      if (nrow(m) < 200000L) stop(sprintf("\"%s\" has %s usable variants; MR-APSS needs genome-wide summary statistics.", tr, format(nrow(m), big.mark = ",")))
      data.table::fwrite(m, file.path(h$dir, f), sep = "\t")
      jf <- sub("\\.tsv\\.gz$", "_preparation.json", f); atomic_write_text(file.path(h$dir, jf), jsonlite::toJSON(c(list(trait = tr, dataset_sha256 = ds[[k]]$recipe$original$sha256), attr(m, "ledger")), auto_unbox = TRUE, pretty = TRUE))
      manifest_register(h$dir, c(f, jf), role = "mrapss_input") } }
  rd <- file.path("04_across_pairs", "mrapss"); dir.create(file.path(h$dir, rd), recursive = TRUE, showWarnings = FALSE); sf <- file.path(rd, paste0(pair$folder, ".R"))
  writeLines(mrapss_script(pair$folder, files, c(pair$exName, pair$outName), normalizePath(lr$scores, winslash = "/"),
                           normalizePath(lr$bfile, winslash = "/", mustWork = FALSE), plink_bin(), .seed_of(paste("mrapss", pair$folder))), file.path(h$dir, sf))
  ldref_save(h, pair, lr)   # Round 81: the LD reference (and any acceptance) beside the result
  sf
}

.mrapss_status_file <- function(dir, folder) file.path(dir, "04_across_pairs", "mrapss", paste0(folder, "_status.json"))
mrapss_job_read <- function(h, pair) { f <- .mrapss_status_file(h$dir, pair$folder); if (!file.exists(f)) return(NULL); j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  if (!is.null(j) && identical(j$state, "running") && !is.null(j$pid) && !.pid_alive(j$pid)) j <- list(state = "stopped", msg = "The run stopped (the app or the computer was closed). Run it again.", time = j$time)
  j }
mrapss_pair_start <- function(h, pair, root = studies_root(), accept = NULL) {
  j <- mrapss_job_read(h, pair); if (!is.null(j) && identical(j$state, "running")) stop("MR-APSS is already running for this pair.")
  sf <- mrapss_pair_prepare(h, pair, root, accept = accept)
  script_job_start(h$dir, sf, .mrapss_status_file(h$dir, pair$folder), file.path(h$dir, "04_across_pairs", "mrapss", paste0(pair$folder, "_job.log")))
  invisible(mrapss_job_read(h, pair))
}
mrapss_pair_read <- function(h, pair) { f <- file.path(h$dir, "04_across_pairs", "mrapss", paste0(pair$folder, "_results.csv"))
  if (!file.exists(f)) return(NULL); r <- utils::read.csv(f, stringsAsFactors = FALSE); attr(r, "labels") <- c(x = pair$exName, y = pair$outName); r }

MRAPSS_METHODS_TEXT <- "The causal effect was also estimated with MR-APSS [Hu 2022], which accounts for pleiotropy and sample structure (population stratification, cryptic relatedness and sample overlap) estimated from the full summary statistics, selects instruments at P < 5e-5 (clumped at r2 < 0.001 within 1,000 kb against the 1000 Genomes European panel) and corrects for the resulting selection bias. Like MRlap, it works on the standardised scale (each effect is the Z statistic divided by the square root of the sample size), so its estimate is in standard deviations and cannot be compared in size with the IVW and other two-sample estimates."   # Round 101: scale stated
mrapss_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); r <- res[1, ]
  sprintf("MR-APSS (%d instruments at P < 5e-5, corrected for selection bias, pleiotropy and sample structure), on the standardised scale: %s on %s %.3f (95%% CI %.3f to %.3f; P = %s).",
          r$instruments, l[["x"]], l[["y"]], r$estimate, r$ci_lower, r$ci_upper, format(signif(r$p, 2)))
}
