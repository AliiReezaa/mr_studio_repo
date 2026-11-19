# ---- Round 65: LHC-MR, latent heritable confounder MR (Darrous et al. 2021, Nat Commun 12:7274) -----------------------------------
# lhcMR (GitHub LizaDarrous/lhcMR, MIT), pinned. Estimates the causal effect in both directions together with the effects of a latent
# heritable confounder, the direct heritabilities and the cross-trait intercept, from genome-wide summary statistics of both traits.
# Needs, besides the full GWAS, the authors' LD scores (UK10K, 4.65 million variants) and local-LD parameters, and the LD scores
# and HapMap3 list of LDSC (ldsc_core.R) for the cross-trait intercept used as a starting point. The authors' files are fetched
# from the links in the package README and checked by MD5.
# One correction to the package, applied by preparing the inputs: merge_sumstats() builds the LD table without the HLA region
# (chr6:28.5-33.5 Mb, as the paper and README state) but then joins the unfiltered table, so HLA variants stay in. MR Studio
# leaves them out of both traits' files before the call, which gives the documented behaviour.
PIN_LHCMR <- "aee7f348744dc7be8b4e7cd52d645ee225cd8b95"
LHCMR_REF <- list(
  ld = list(url = "https://drive.usercontent.google.com/download?id=1_B8DsGejRpeiW3XxNERWpmSQH2a8NWqS&export=download&confirm=t", file = "LDscores_filtered.csv", md5 = "2e8bffd59f56bfbfadf074222f54ecd1"),
  rho = list(url = "https://drive.usercontent.google.com/download?id=1PwsJ2KiCR2VDVOodEXnQu6Tec54a0HC1&export=download&confirm=t", file = "LD_GM2_2prm.csv", md5 = "75f5309d4dd58180cd610816073b5a59"))
LHCMR_M <- 1e7   # number of variants behind the authors' LD scores, the package default for their files
LHCMR_HLA <- c(28.5e6, 33.5e6)

lhcmr_ref_dir <- function() file.path(ldsc_ref_dir(), "lhcmr")
lhcmr_ref_status <- function(dir = lhcmr_ref_dir()) list(dir = dir, ld = file.exists(file.path(dir, LHCMR_REF$ld$file)), rho = file.exists(file.path(dir, LHCMR_REF$rho$file)),
                                                         package = requireNamespace("lhcMR", quietly = TRUE))
lhcmr_ref_fetch <- function(dir = lhcmr_ref_dir(), files = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (w in c("ld", "rho")) { r <- LHCMR_REF[[w]]; f <- file.path(dir, r$file)
    if (!is.null(files[[w]])) file.copy(files[[w]], f, overwrite = TRUE) else if (!(file.exists(f) && identical(unname(tools::md5sum(f)), r$md5))) utils::download.file(r$url, f, mode = "wb", quiet = TRUE)
    if (!identical(unname(tools::md5sum(f)), r$md5)) { unlink(f); stop(sprintf("The downloaded %s is not the authors' pinned file (MD5 differs).", r$file)) } }
  lhcmr_ref_status(dir)
}

#' Summary statistics of one trait in the form lhcMR reads (SNP CHR POS A1 A2 BETA SE TSTAT P N [MAF]): variants of the authors'
#' LD file (chromosome and position taken from it when the dataset has none), HLA left out. Returns the table with a ledger.
lhcmr_prepare <- function(d, ldpos) {
  led <- list(rows = nrow(d))
  if (!("snp" %in% names(d))) stop("LHC-MR needs rsIDs; this dataset identifies variants by position only.")
  if (!("n" %in% names(d)) || all(is.na(d$n))) stop("LHC-MR needs the sample size; give it when importing the dataset.")
  i <- match(d$snp, ldpos$RSID); d <- d[!is.na(i), , drop = FALSE]; i <- i[!is.na(i)]; led$in_ld_file <- nrow(d)
  chr <- ldpos$CHR[i]; pos <- ldpos$POS[i]; hla <- chr == 6 & pos >= LHCMR_HLA[1] & pos <= LHCMR_HLA[2]
  d <- d[!hla, , drop = FALSE]; chr <- chr[!hla]; pos <- pos[!hla]; led$hla_removed <- sum(hla)
  z <- d$beta / d$se
  out <- data.frame(SNP = d$snp, CHR = chr, POS = pos, A1 = d$ea, A2 = d$oa, BETA = d$beta, SE = d$se, TSTAT = z, P = 2 * stats::pnorm(-abs(z)), N = d$n, stringsAsFactors = FALSE)
  if ("eaf" %in% names(d) && any(is.finite(d$eaf))) out$MAF <- pmin(d$eaf, 1 - d$eaf)
  out <- out[is.finite(out$TSTAT) & is.finite(out$N), , drop = FALSE]; led$kept <- nrow(out); attr(out, "ledger") <- led; out
}

#' Script for one pair. `cores` > 1 uses parallel::mclapply (not on Windows).
lhcmr_script <- function(folder, files, labels, ldf, rhof, ldsc_dir, hm3, seed, cores, sp_pair = 50L, nblock = 200L, bfile = "", plink = "") {
  q <- function(s) encodeString(s, quote = "\"")
  c("# MR Studio: LHC-MR, latent heritable confounder Mendelian randomisation (Darrous et al. 2021, Nat Commun 12:7274).",
    sprintf("# Package: lhcMR commit %s (GitHub LizaDarrous/lhcMR); cross-trait intercept with GenomicSEM commit %s.", PIN_LHCMR, PIN_GENOMICSEM),
    "# Run from the study folder. The summary statistics in 02_datasets/lhcmr were prepared by MR Studio from the Data library files",
    "# (variants of the authors' LD file, HLA chr6:28.5-33.5 Mb left out as the paper states); each has a record of what was removed.",
    sprintf("PAIR <- %s", q(folder)), sprintf("TRAITS <- c(%s)   # X = %s, Y = %s", paste(vapply(files, q, ""), collapse = ", "), labels[1], labels[2]),
    sprintf("LD_FILE <- %s; RHO_FILE <- %s   # authors' LD scores and local-LD parameters", q(ldf), q(rhof)),
    sprintf("LDSC_LD <- %s; HM3 <- %s", q(ldsc_dir), q(hm3)),
    sprintf("seed <- %dL; CORES <- %dL; SP_PAIR <- %dL; N_BLOCK <- %dL", as.integer(seed), as.integer(cores), as.integer(sp_pair), as.integer(nblock)),
    "if (!requireNamespace(\"lhcMR\", quietly = TRUE)) stop(\"Package lhcMR is not installed (GitHub LizaDarrous/lhcMR).\")",
    "suppressPackageStartupMessages(library(lhcMR))   # the package calls fread, inner_join and %>% without prefixes",
    sprintf("BFILE <- %s; PLINK <- %s   # 1000 Genomes EUR panel and PLINK of MR Studio", q(bfile), q(plink)),
    "# The IVW starting points clump instruments with TwoSampleMR::clump_data, which by default asks the OpenGWAS server (a token is",
    "# needed). Here the same function clumps locally with PLINK against the same 1000 Genomes EUR panel (its bfile and plink_bin arguments).",
    "if (nzchar(BFILE) && nzchar(PLINK)) { cd <- TwoSampleMR::clump_data; utils::assignInNamespace(\"clump_data\", function(dat, ...) cd(dat, ..., bfile = BFILE, plink_bin = PLINK), \"TwoSampleMR\") }",
    "od <- normalizePath(file.path(\"04_across_pairs\", \"lhcmr\", PAIR), mustWork = FALSE); dir.create(od, recursive = TRUE, showWarnings = FALSE)",
    "inp <- lapply(TRAITS, function(f) as.data.frame(data.table::fread(f)))",
    "old <- setwd(od); on.exit(setwd(old), add = TRUE)   # lhcMR writes its tables into the working directory",
    "RNGkind(\"L'Ecuyer-CMRG\", \"Inversion\", \"Rejection\"); set.seed(seed)   # also reproducible with parallel::mclapply",
    "df <- merge_sumstats(inp, c(\"X\", \"Y\"), LD_FILE, RHO_FILE)",
    "sp <- calculate_SP(df, c(\"X\", \"Y\"), run_ldsc = TRUE, run_MR = TRUE, saveRFiles = TRUE, hm3 = HM3, ld = LDSC_LD, nStep = 2, SP_single = 3, SP_pair = SP_PAIR, SNP_filter = 10, nCores = CORES)",
    "lhc_mr(sp, c(\"X\", \"Y\"), paral_method = \"lapply\", nCores = CORES, nBlock = N_BLOCK)",
    "s <- utils::read.csv(\"SummarisedResults_X-Y.csv\", row.names = 1, check.names = FALSE)",
    "out <- data.frame(parameter = colnames(s), estimate = as.numeric(s[1, ]), se = as.numeric(s[2, ]), p = as.numeric(s[3, ]), stringsAsFactors = FALSE)",
    "out$ci_lower <- out$estimate - 1.96 * out$se; out$ci_upper <- out$estimate + 1.96 * out$se; out$variants <- nrow(df); out$variants_used <- nrow(sp$input.df_filtered)",
    "jk <- utils::read.csv(sprintf(\"%dblockJK_summary_X-Y.csv\", N_BLOCK)); out$bimodal <- jk$bimod[match(out$parameter, jk$Parameter)] %in% c(TRUE, \"TRUE\")",
    "out$axy_ivw <- sp$axy_MR; out$ayx_ivw <- sp$ayx_MR",
    "setwd(old); utils::write.csv(out, file.path(\"04_across_pairs\", \"lhcmr\", paste0(PAIR, \"_results.csv\")), row.names = FALSE)",
    "out")
}

# Each parallel worker holds a copy of the merged data (several GB with full GWAS), so the default stays gentle: at most 2 workers,
# fewer on small machines, 1 on Windows (no forking). MRSTUDIO_LHCMR_CORES overrides it. The run is also given low priority.
lhcmr_cores <- function() { e <- suppressWarnings(as.integer(Sys.getenv("MRSTUDIO_LHCMR_CORES", ""))); if (is.finite(e) && e >= 1L) return(e)
  if (.Platform$OS.type == "windows") 1L else max(1L, min(2L, parallel::detectCores() %/% 4L)) }

#' Prepare both datasets, write and run the LHC-MR script of a pair (long: one to several hours on the default 2 workers).
lhcmr_pair_run <- function(h, pair, root = studies_root(), dir = lhcmr_ref_dir(), ldsc = ldsc_ref_dir(), cores = lhcmr_cores(), sp_pair = 50L, nblock = 200L) {
  if (isTRUE(pair$mv) || is.null(pair$folder)) stop("LHC-MR is for a pair with one exposure.")
  ids <- c(pair$ex, pair$out); if (!all(startsWith(ids, "lib:"))) stop("LHC-MR needs the full summary statistics of both traits, imported into the Data library.")
  st <- lhcmr_ref_status(dir); if (!st$ld || !st$rho) stop("The LD files of the LHC-MR authors are not installed yet.")
  if (!ldsc_ref_status(ldsc)$official) stop("LHC-MR uses the European LD scores of LDSC for the cross-trait intercept; install them first.")
  ds <- lapply(ids, function(i) library_dataset(sub("^lib:", "", i), root))
  anc <- vapply(ds, function(x) toupper(x$recipe$facts$ancestry %||% "EUR"), ""); if (!all(anc == "EUR")) stop("LHC-MR's LD files are from Europeans (UK10K); both GWAS must be of European ancestry.")
  ldpos <- NULL; rel <- file.path("02_datasets", "lhcmr"); dir.create(file.path(h$dir, rel), recursive = TRUE, showWarnings = FALSE); files <- character(2)
  for (k in 1:2) { tr <- ds[[k]]$recipe$trait; f <- file.path(rel, paste0(.ldsc_slug(tr), ".tsv.gz")); files[k] <- f
    if (!file.exists(file.path(h$dir, f))) {
      if (is.null(ldpos)) { ldpos <- as.data.frame(data.table::fread(file.path(dir, LHCMR_REF$ld$file), select = c("chr", "pos", "rs"))); names(ldpos) <- c("CHR", "POS", "RSID") }   # the authors' file: chr, pos, rs
      m <- lhcmr_prepare(ds[[k]]$data, ldpos)
      if (nrow(m) < 500000L) stop(sprintf("\"%s\" has %s variants in the LHC-MR LD file; LHC-MR needs genome-wide summary statistics.", tr, format(nrow(m), big.mark = ",")))
      data.table::fwrite(m, file.path(h$dir, f), sep = "\t")
      jf <- sub("\\.tsv\\.gz$", "_preparation.json", f); atomic_write_text(file.path(h$dir, jf), jsonlite::toJSON(c(list(trait = tr, dataset_sha256 = ds[[k]]$recipe$original$sha256), attr(m, "ledger")), auto_unbox = TRUE, pretty = TRUE))
      manifest_register(h$dir, c(f, jf), role = "lhcmr_input") } }
  rd <- file.path("04_across_pairs", "lhcmr"); dir.create(file.path(h$dir, rd), recursive = TRUE, showWarnings = FALSE); sf <- file.path(h$dir, rd, paste0(pair$folder, ".R"))
  writeLines(lhcmr_script(pair$folder, files, c(pair$exName, pair$outName), normalizePath(file.path(dir, LHCMR_REF$ld$file), winslash = "/"), normalizePath(file.path(dir, LHCMR_REF$rho$file), winslash = "/"),
                          normalizePath(ldsc_ld_folder("official", ldsc), winslash = "/"), normalizePath(file.path(ldsc, "w_hm3.snplist"), winslash = "/"), .seed_of(paste("lhcmr", pair$folder)), cores, sp_pair, nblock,
                          bfile = if (file.exists(file.path(ld_panel_dir(), "EUR.bed"))) normalizePath(file.path(ld_panel_dir(), "EUR"), winslash = "/", mustWork = FALSE) else "", plink = if (is.na(plink_bin())) "" else plink_bin()), sf)
  old <- setwd(h$dir); res <- tryCatch({ utils::capture.output(sys.source(sf, envir = new.env(parent = globalenv()), keep.source = FALSE))
    utils::read.csv(file.path(rd, paste0(pair$folder, "_results.csv")), stringsAsFactors = FALSE) }, finally = setwd(old))
  manifest_register(h$dir, c(file.path(rd, paste0(pair$folder, c(".R", "_results.csv"))), file.path(rd, pair$folder, list.files(file.path(h$dir, rd, pair$folder)))), role = "lhcmr")
  attr(res, "labels") <- c(x = pair$exName, y = pair$outName); res
}

lhcmr_pair_read <- function(h, pair) { f <- file.path(h$dir, "04_across_pairs", "lhcmr", paste0(pair$folder, "_results.csv"))
  if (!file.exists(f)) return(NULL); r <- utils::read.csv(f, stringsAsFactors = FALSE); attr(r, "labels") <- c(x = pair$exName, y = pair$outName); r }

LHCMR_METHODS_TEXT <- "Causal effects in both directions were also estimated with LHC-MR, which models a latent heritable confounder of the two traits and uses genome-wide summary statistics rather than instruments only; standard errors are from a 200-block jackknife [Darrous 2021]."
lhcmr_text <- function(res) {
  if (!is.data.frame(res) || !nrow(res)) return(character())
  l <- attr(res, "labels"); g <- function(p) res[res$parameter == p, ][1, ]
  f <- function(r) sprintf("%.3f (95%% CI %.3f to %.3f; P = %s)", r$estimate, r$ci_lower, r$ci_upper, format(signif(r$p, 2)))
  tx <- g("tX"); ty <- g("tY")
  sprintf("LHC-MR (%s variants, every tenth used): %s on %s %s, and %s on %s %s; latent confounder effects %.3f on %s and %.3f on %s (P = %s and %s); cross-trait intercept %.3f.%s",
          format(res$variants[1], big.mark = ","), l[["x"]], l[["y"]], f(g("axy")), l[["y"]], l[["x"]], f(g("ayx")), tx$estimate, l[["x"]], ty$estimate, l[["y"]], format(signif(tx$p, 2)), format(signif(ty$p, 2)), g("iXY")$estimate,
          if (any(res$bimodal)) " The jackknife distribution of at least one parameter was bimodal, which the authors read as two competing solutions; interpret with care." else "")
}

# ---- Background run (one to two hours): separate R process, status in 04_across_pairs/lhcmr/<pair>_status.json ----------------------
.lhcmr_status_file <- function(dir, folder) file.path(dir, "04_across_pairs", "lhcmr", paste0(folder, "_status.json"))
lhcmr_job_read <- function(h, pair) { f <- .lhcmr_status_file(h$dir, pair$folder); if (!file.exists(f)) return(NULL); tryCatch(jsonlite::fromJSON(f), error = function(e) NULL) }
lhcmr_pair_start <- function(h, pair, root = studies_root(), engine = getOption("mrgui.engine_dir", ""), cores = lhcmr_cores()) {
  j <- lhcmr_job_read(h, pair); if (!is.null(j) && identical(j$state, "running")) stop("LHC-MR is already running for this pair.")
  sf <- .lhcmr_status_file(h$dir, pair$folder); dir.create(dirname(sf), recursive = TRUE, showWarnings = FALSE)
  pj <- file.path(dirname(sf), paste0(pair$folder, "_job.json")); writeLines(jsonlite::toJSON(pair[c("folder", "ex", "out", "exName", "outName")], auto_unbox = TRUE), pj)
  job <- file.path(dirname(sf), paste0(pair$folder, "_job.R"))
  writeLines(c("a <- commandArgs(TRUE); engine <- a[1]; sdir <- a[2]; root <- a[3]; pj <- a[4]; sf <- a[5]; cores <- as.integer(a[6])",
    "st <- function(state, msg = '') writeLines(jsonlite::toJSON(list(state = state, msg = msg, time = format(Sys.time(), '%Y-%m-%d %H:%M:%S')), auto_unbox = TRUE), sf)",
    "invisible(tools::psnice(value = 15L))   # low priority, so the computer stays responsive during the one to two hours",
    "tryCatch({ options(mrgui.engine_dir = engine); source(file.path(engine, 'studio_core.R')); load_contracts(find_contracts_dir(engine))",
    "  for (f in c('import_core.R', 'ld_core.R', 'ldsc_core.R', 'lhcmr_core.R')) source(file.path(engine, f))",
    "  pair <- jsonlite::fromJSON(pj); lhcmr_pair_run(list(dir = sdir), pair, root = root, cores = cores); st('done')",
    "}, error = function(e) st('error', conditionMessage(e)))"), job)
  writeLines(jsonlite::toJSON(list(state = "running", msg = "", time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), auto_unbox = TRUE), sf)
  env <- c(sprintf("MRGUI_LDSC_DIR=%s", ldsc_ref_dir()), sprintf("MRGUI_STUDIES_DIR=%s", root))
  system2(default_rscript(), c("--vanilla", shQuote(job), shQuote(engine), shQuote(h$dir), shQuote(root), shQuote(pj), shQuote(sf), cores), env = env,
          stdout = file.path(dirname(sf), paste0(pair$folder, "_job.log")), stderr = file.path(dirname(sf), paste0(pair$folder, "_job.log")), wait = FALSE)
  invisible(lhcmr_job_read(h, pair))
}

#' Download the authors' two LD files (about 500 MB) in a separate R process; progress is the size of the files.
lhcmr_fetch_start <- function(dir = lhcmr_ref_dir(), engine = getOption("mrgui.engine_dir", "")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE); sf <- file.path(dir, "fetch_status.json")
  writeLines(jsonlite::toJSON(list(state = "running", msg = ""), auto_unbox = TRUE), sf)
  code <- sprintf("sf <- %s; tryCatch({ source(file.path(%s, 'ldsc_core.R')); source(file.path(%s, 'lhcmr_core.R')); lhcmr_ref_fetch(%s); writeLines(jsonlite::toJSON(list(state = 'done', msg = ''), auto_unbox = TRUE), sf) }, error = function(e) writeLines(jsonlite::toJSON(list(state = 'error', msg = conditionMessage(e)), auto_unbox = TRUE), sf))",
                  deparse(sf), deparse(engine), deparse(engine), deparse(dir))
  system2(default_rscript(), c("--vanilla", "-e", shQuote(code)), stdout = FALSE, stderr = FALSE, wait = FALSE); invisible(TRUE)
}
lhcmr_status <- function(h = NULL, pair = NULL, dir = lhcmr_ref_dir()) {
  f <- file.path(dir, "fetch_status.json"); fj <- if (file.exists(f)) tryCatch(jsonlite::fromJSON(f), error = function(e) NULL) else NULL
  mb <- sum(file.size(file.path(dir, c(LHCMR_REF$ld$file, LHCMR_REF$rho$file))), na.rm = TRUE) / 1e6
  c(lhcmr_ref_status(dir), list(fetch = fj, fetched_mb = round(mb), ldsc_official = ldsc_ref_status()$official, job = if (!is.null(h) && !is.null(pair)) lhcmr_job_read(h, pair) else NULL))
}
