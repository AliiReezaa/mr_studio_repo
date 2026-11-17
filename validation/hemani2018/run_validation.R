# Validation study 2: Hemani et al. 2018, the MR-Base applied example (LDL cholesterol on coronary heart disease).
# Paper: eLife 7:e34408, doi:10.7554/eLife.34408 (PMID 29846171). Run on OpenGWAS (formerly MR-Base), so it tests MR Studio's
# download path: instruments with each of the three clumping choices (stored clumps, live server clumping, this computer), outcome
# lookup with proxies, harmonisation with allele frequencies, the generated scripts and the replay check.
#
#   Tier A (exact): the same OpenGWAS datasets through MR Studio's engine and through TwoSampleMR / ieugwasr directly.
#   Tier B (approximate): MR Studio against the numbers printed in the paper (2018 data; OpenGWAS may have changed since).
#
# Run from the repository root, in a Terminal where the token is set for that session only (never in .Renviron or a file):
#   export OPENGWAS_JWT='...'          # then:
#   Rscript --vanilla validation/hemani2018/run_validation.R
#   unset OPENGWAS_JWT
# The token is read from the environment, passed to the calls and never printed or written. Everything downloaded goes to a
# throw-away folder; only estimates and instrument rsID lists are written to validation/hemani2018/.
# Optional: MRV_EXPOSURE / MRV_OUTCOME choose other dataset IDs; MRV_LD_DIR the folder with EUR.bed/.bim/.fam (default: the app's panel).
# MRV_DRYRUN=1 runs everything on simulated data with real 1000 Genomes rsIDs and a fake OpenGWAS, no token and no internet.
repo <- normalizePath("."); vdir <- file.path(repo, "validation", "hemani2018"); eng <- file.path(repo, "prototype", "engine")
dry <- identical(Sys.getenv("MRV_DRYRUN"), "1")
jwt <- Sys.getenv("OPENGWAS_JWT", ""); if (dry) jwt <- "dry.run.token"
if (!nzchar(jwt)) stop("Set OPENGWAS_JWT in this Terminal session first (see the top of this script).", call. = FALSE)
X_ID <- Sys.getenv("MRV_EXPOSURE", "ieu-a-300"); Y_ID <- Sys.getenv("MRV_OUTCOME", "ieu-a-7")
LD_DIR <- Sys.getenv("MRV_LD_DIR", file.path(tools::R_user_dir("MRStudio", "data"), "ld_ref", "1kg.v3"))
stopifnot(all(file.exists(file.path(LD_DIR, paste0("EUR", c(".bed", ".bim", ".fam"))))))
PLINK <- file.path(eng, "tools", "plink", "plink"); stopifnot(file.exists(PLINK))
tmp <- tempfile("hemani-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"), MRGUI_LD_DIR = LD_DIR, MRGUI_PLINK = PLINK)
suppressMessages({ source(file.path(eng, "studio_core.R")); load_contracts(find_contracts_dir(eng))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R")) source(file.path(eng, f)) })
out_tag <- if (dry) "_dryrun" else ""; vout <- if (dry) tempdir() else vdir   # a dry run writes nothing into the repository

# ---- the reference half: TwoSampleMR and ieugwasr called directly (replaced by the same fake as the engine in a dry run) -------------
REF <- list(
  info = function(ids) as.data.frame(ieugwasr::gwasinfo(ids, opengwas_jwt = jwt)),
  instruments = function(id) TwoSampleMR::extract_instruments(id, opengwas_jwt = jwt),   # defaults: 5e-8, r2 0.001, 10,000 kb: stored clumps
  tophits = function(id, clump, force) as.data.frame(ieugwasr::tophits(id, clump = clump, force_server = force, opengwas_jwt = jwt)),
  outcome = function(snps, id) suppressMessages(TwoSampleMR::extract_outcome_data(snps, id, proxies = TRUE, rsq = 0.8, maf_threshold = 0.3, opengwas_jwt = jwt)))

if (dry) {   # simulated LDL-like exposure on real panel variants (pairs of neighbours spread over the genome), so local clumping is real
  bim <- utils::read.table(pipe(paste("awk '$2 ~ /^rs/ && (NR % 120000 < 3)' ", shQuote(file.path(LD_DIR, "EUR.bim")))), stringsAsFactors = FALSE)
  bim <- bim[grepl("^[ACGT]$", bim$V5) & grepl("^[ACGT]$", bim$V6), ]
  set.seed(2018); k <- nrow(bim); bx <- stats::rnorm(k, 0, 0.004) + ifelse(stats::runif(k) < 0.8, 0.06, 0); sx <- rep(0.004, k)
  by <- 0.37 * bx + stats::rnorm(k, 0, 0.01); sy <- rep(0.01, k)
  TH <- data.frame(rsid = bim$V2, chr = "19", position = bim$V4, ea = bim$V5, nea = bim$V6, eaf = stats::runif(k, 0.1, 0.9), beta = bx, se = sx, p = 2 * stats::pnorm(-abs(bx / sx)), n = 173082, id = X_ID, trait = "LDL cholesterol", stringsAsFactors = FALSE)
  TH <- TH[TH$p < 5e-8, ]; TH <- TH[order(TH$p), ]
  # Round 28: the strongest variants of six neighbouring triplets get P = 1e-200, as OpenGWAS stores them, so ties decide some leads
  tri <- utils::head(unique(match(TH$rsid, bim$V2) %/% 3L), 6L); TH$p[(match(TH$rsid, bim$V2) %/% 3L) %in% tri] <- 1e-200; TH <- TH[order(TH$p), ]
  cl <- ieugwasr::ld_clump_local(data.frame(rsid = TH$rsid, pval = TH$p), 10000, 0.001, 1, file.path(LD_DIR, "EUR"), PLINK); CL <- TH[TH$rsid %in% cl$rsid, ]
  CAT <- data.frame(id = c(X_ID, Y_ID), trait = c("LDL cholesterol", "Coronary heart disease"), population = "European", sex = "Males and Females", build = "HG19/GRCh37", author = "Sim", consortium = c("GLGC", "CARDIoGRAMplusC4D"),
    category = "Risk factor", subcategory = "Lipid", unit = "SD", note = "", pmid = "1", ncase = c(NA, 60801), ncontrol = c(NA, 123504), sample_size = c(173082, 184305), nsnp = 1e7, year = 2015, mr = 1, stringsAsFactors = FALSE)
  # fixed variants, whatever is asked: the weakest server instrument is missing, the two strongest come by proxy
  fake_out <- function(snps, id) { i <- match(snps, bim$V2); i <- i[!is.na(i)]; prox <- match(utils::head(CL$rsid, 2), bim$V2); i <- setdiff(i, match(utils::tail(CL$rsid, 1), bim$V2))
    data.frame(SNP = bim$V2[i], chr = "19", pos = bim$V4[i], beta.outcome = by[i], se.outcome = sy[i], samplesize.outcome = 184305, pval.outcome = 2 * stats::pnorm(-abs(by[i] / sy[i])),
      eaf.outcome = TH$eaf[match(bim$V2[i], TH$rsid)], effect_allele.outcome = bim$V5[i], other_allele.outcome = bim$V6[i], outcome = "Coronary heart disease", id.outcome = id, proxy.outcome = i %in% prox,
      target_snp.outcome = ifelse(i %in% prox, bim$V2[i], NA), proxy_snp.outcome = ifelse(i %in% prox, paste0(bim$V2[i], "p"), NA), mr_keep.outcome = TRUE, data_source.outcome = "igd", stringsAsFactors = FALSE) }
  REF <- list(info = function(ids) CAT[CAT$id %in% ids, ], tophits = function(id, clump, force) if (clump == 0) TH else CL,
    instruments = function(id) TwoSampleMR::format_data(data.frame(SNP = CL$rsid, beta = CL$beta, se = CL$se, effect_allele = CL$ea, other_allele = CL$nea, eaf = CL$eaf, pval = CL$p, samplesize = CL$n), type = "exposure"),
    outcome = function(snps, id) fake_out(snps, id))
  options(mrgui.ogw_backend = list(info = function(jwt, id = NULL) REF$info(id), tophits = function(id, p, r2, kb, pop, jwt, clump = 1, force = FALSE) REF$tophits(id, clump, force),
    outcome = function(snps, id, proxies, rsq, maf, jwt) REF$outcome(snps, id)))
}

# ---- 1. which datasets: the paper names the consortia, not the OpenGWAS IDs ---------------------------------------------------------
meta <- REF$info(c(X_ID, Y_ID))
cat("\nDatasets used (the paper reports GLGC with 214,370 subjects and CARDIoGRAMplusC4D with 195,813 subjects):\n")
print(meta[, intersect(c("id", "trait", "consortium", "author", "year", "population", "sample_size", "ncase", "ncontrol", "nsnp", "pmid"), names(meta))], row.names = FALSE)

# ---- 2. instruments: MR Studio, each clumping choice, against the reference -----------------------------------------------------------
cat("\nDownloading instruments for", X_ID, "with each clumping choice...\n")
# The catalogue gives "Mixed" ancestry for both datasets, so the reference population is chosen (EUR, the TwoSampleMR default).
ex <- list(stored = ogw_get_exposure(X_ID, jwt, pop = "EUR", where = "server"), live = ogw_get_exposure(X_ID, jwt, pop = "EUR", where = "server_live"), local = ogw_get_exposure(X_ID, jwt, pop = "EUR", where = "local"))
exd <- lapply(ex, function(r) library_dataset(r$trait)$data)
ref_stored <- REF$instruments(X_ID)
ref_live <- REF$tophits(X_ID, clump = 1, force = TRUE)
ref_cand <- REF$tophits(X_ID, clump = 0, force = FALSE); ref_cand <- ref_cand[order(ref_cand$p), ]
ref_local <- ieugwasr::ld_clump_local(data.frame(rsid = ref_cand$rsid, pval = ref_cand$p), clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = file.path(LD_DIR, "EUR"), plink_bin = PLINK)
ref_sets <- list(stored = ref_stored$SNP, live = ref_live$rsid, local = ref_local$rsid)
ties <- sum(duplicated(ref_cand$p[ref_cand$p < 5e-8]))
# Round 28, tie sensitivity reference: equal P broken by the larger |beta/SE|, given to ieugwasr as ranks (the same rule, written independently)
zc <- abs(ref_cand$beta / ref_cand$se); o <- order(ref_cand$p, -zc, method = "radix"); rk <- integer(nrow(ref_cand)); rk[o] <- seq_along(o)
ref_tie <- ieugwasr::ld_clump_local(data.frame(rsid = ref_cand$rsid, pval = rk / (nrow(ref_cand) + 1)), clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = file.path(LD_DIR, "EUR"), plink_bin = PLINK)
same_set <- function(a, b) length(a) == length(b) && setequal(a, b)
inst <- do.call(rbind, lapply(names(ex), function(k) { a <- exd[[k]]$snp; b <- ref_sets[[k]]
  data.frame(choice = k, app_n = length(a), ref_n = length(b), shared = length(intersect(a, b)), identical = same_set(a, b), stringsAsFactors = FALSE) }))
tsum <- ex$local$clump$ties; tie_run <- isTRUE(tsum$regions > 0)
if (tie_run) { ts0 <- lib_table_read(file.path(library_dir(), ex$local$trait), "instrument_selection_ties"); app_tie <- ts0$variant_id[ts0$status == "index"]
  inst <- rbind(inst, data.frame(choice = "local_ties", app_n = length(app_tie), ref_n = nrow(ref_tie), shared = length(intersect(app_tie, ref_tie$rsid)), identical = same_set(app_tie, ref_tie$rsid), stringsAsFactors = FALSE)) }
cat("\nInstruments, MR Studio against the reference (stored: TwoSampleMR::extract_instruments; live: ieugwasr::tophits force_server; local: ieugwasr::ld_clump_local on the same panel):\n"); print(inst, row.names = FALSE)
cat(sprintf("Candidates below 5e-8: %d; equal P values among them: %d. The primary analysis keeps PLINK's choice among them, as ieugwasr and the server do.\n", sum(ref_cand$p < 5e-8), ties))
cat(sprintf("Tie summary in MR Studio: %s candidates tied; ties decided %s lead variant(s). Kept: %s. Larger |beta/SE|: %s.\n", tsum$tied %||% 0, tsum$regions %||% 0, paste(unlist(tsum$out), collapse = " "), paste(unlist(tsum[["in"]]), collapse = " ")))
cat(sprintf("Between choices in MR Studio: stored and live share %d of %d/%d; live and local identical: %s.\n", length(intersect(exd$stored$snp, exd$live$snp)), nrow(exd$stored), nrow(exd$live), same_set(exd$live$snp, exd$local$snp)))
allsnp <- sort(unique(unlist(lapply(exd, `[[`, "snp"))))
allsnp <- sort(unique(c(allsnp, if (tie_run) app_tie)))
utils::write.csv(data.frame(rsid = allsnp, stored = allsnp %in% exd$stored$snp, live = allsnp %in% exd$live$snp, local = allsnp %in% exd$local$snp, local_ties = if (tie_run) allsnp %in% app_tie else NA),
  file.path(vout, paste0("instruments", out_tag, ".csv")), row.names = FALSE)

# ---- 3. outcome, harmonisation and estimates ------------------------------------------------------------------------------------------
METH <- c("ivw_mre", "egger", "wm", "smed", "mode")
TSMR <- c(ivw_mre = "mr_ivw", egger = "mr_egger_regression", wm = "mr_weighted_median", smed = "mr_simple_median", mode = "mr_weighted_mode")
key <- function(m, e) if (m == "egger" && e == "pleiotropy_intercept") "egger_int" else m
h <- studio_new("Hemani2018 LDL CHD"); rows <- list(); i <- 0L
CHOICES <- c(names(ex), if (tie_run) "local_ties")
for (k in CHOICES) {
  i <- i + 1L; cat("\nPair", i, "(", k, "): outcome lookup with proxies, harmonisation, methods, replay...\n")
  kx <- if (k == "local_ties") "local" else k
  oy <- ogw_get_outcome(Y_ID, ex[[kx]]$trait, jwt, proxies = TRUE)   # for local_ties: the same outcome download, reused
  if (k == "local_ties") { lp <- pair_input(paste0("lib:", ex$local$trait), paste0("lib:", oy$trait), ties = "strongest")
    rr <- studio_run_ties(h, paste0("lib:", ex$local$trait), paste0("lib:", oy$trait), METH, idx = i)
  } else { lp <- pair_input(paste0("lib:", ex[[k]]$trait), paste0("lib:", oy$trait))
    rr <- studio_run(h, lp$exposure_label, lp$outcome_label, lp$input, METH, idx = i) }
  pl <- pair_payload(h$dir, rr$pair_folder); led <- lp$ledger; d <- lp$input$data; dk <- d[d$mr_keep, ]
  # reference: the same instruments (the reference's own set for this choice) through TwoSampleMR
  rx <- if (k == "stored") ref_stored else { t <- if (k == "live") ref_live else if (k == "local_ties") ref_cand[ref_cand$rsid %in% ref_tie$rsid, ] else ref_cand[ref_cand$rsid %in% ref_local$rsid, ]
    TwoSampleMR::format_data(data.frame(SNP = t$rsid, beta = t$beta, se = t$se, effect_allele = t$ea, other_allele = t$nea, eaf = t$eaf, pval = t$p, samplesize = t$n), type = "exposure") }
  ry <- REF$outcome(rx$SNP, Y_ID)
  px <- ry$proxy.outcome %in% TRUE   # TwoSampleMR names a proxy row by the instrument it stands for (target), which the fake copies
  if (any(px)) cat(sprintf("  reference: %d proxy rows, named by the target rsID: %s\n", sum(px), all(ry$SNP[px] == ry$target_snp.outcome[px])))
  H <- suppressMessages(TwoSampleMR::harmonise_data(rx, ry, action = 2)); Hk <- H[H$mr_keep, ]
  m <- match(dk$variant_id, Hk$SNP)
  harm_same <- nrow(dk) == nrow(Hk) && !anyNA(m) && isTRUE(all.equal(dk$beta_x, Hk$beta.exposure[m])) && isTRUE(all.equal(dk$beta_y, Hk$beta.outcome[m]))
  if (!harm_same) { only_app <- setdiff(dk$variant_id, Hk$SNP); only_ref <- setdiff(Hk$SNP, dk$variant_id)
    cat("  harmonised data differ. Only in MR Studio:", paste(only_app, collapse = " "), "| only in the reference:", paste(only_ref, collapse = " "), "\n") }
  set.seed(1); R <- suppressMessages(TwoSampleMR::mr(H, method_list = unname(TSMR)))
  het <- TwoSampleMR::mr_heterogeneity(H, method_list = "mr_ivw"); egi <- TwoSampleMR::mr_pleiotropy_test(H)
  q_app <- TwoSampleMR::mr_ivw(dk$beta_x, dk$beta_y, dk$se_x, dk$se_y)
  for (e in pl$estimates) { mk <- key(e$method, e$estimand)
    ref <- if (mk == "egger_int") c(egi$egger_intercept, egi$se) else { r <- R[R$method == TwoSampleMR::mr_method_list()$name[match(TSMR[[mk]], TwoSampleMR::mr_method_list()$obj)], ]; c(r$b, r$se) }
    if (length(ref) < 2) ref <- c(NA_real_, NA_real_)
    rows[[length(rows) + 1L]] <- data.frame(choice = k, method = mk, app_est = e$estimate, app_se = e$se, app_n = e$n_ivs_used, ref_est = ref[1], ref_se = ref[2], ref_n = sum(H$mr_keep),
      harmonised_identical = harm_same, replay = pl$replay$status, q_app = q_app$Q, q_df_app = q_app$Q_df, q_ref = het$Q, q_df_ref = het$Q_df,
      instruments = led$instruments, in_outcome = led$in_outcome, by_proxy = led$by_proxy, missing_in_outcome = led$missing_in_outcome, flipped = led$flipped,
      palindromic_dropped = led$dropped_palindromic, mismatch_dropped = led$dropped_mismatch, kept = led$kept, stringsAsFactors = FALSE) }
}
res <- do.call(rbind, rows)

# ---- 4. the paper (Tier B): odds ratios per SD of LDL cholesterol, transcribed from the open-access text (PMC5976434) ---------------
P <- data.frame(method = c("ivw_mre", "egger", "wm", "mode", "egger_int"), paper_or = c(1.45, 1.85, 1.56, 1.68, NA), paper_lo = c(1.30, 1.48, 1.43, 1.56, NA), paper_hi = c(1.62, 2.32, 1.70, 1.80, NA),
  paper_b = c(NA, NA, NA, NA, -0.013), paper_se = c(NA, NA, NA, NA, 0.005), stringsAsFactors = FALSE)
res <- merge(res, P, by = "method", all.x = TRUE, sort = FALSE)
res$app_or <- ifelse(res$method == "egger_int", NA, exp(res$app_est)); res$app_lo <- ifelse(res$method == "egger_int", NA, exp(res$app_est - 1.96 * res$app_se)); res$app_hi <- ifelse(res$method == "egger_int", NA, exp(res$app_est + 1.96 * res$app_se))
res <- res[order(match(res$choice, CHOICES), match(res$method, c("ivw_mre", "egger", "egger_int", "wm", "smed", "mode"))), ]
utils::write.csv(res, file.path(vout, paste0("results", out_tag, ".csv")), row.names = FALSE)

cat("\n==== Tier A: MR Studio against TwoSampleMR on the same downloads (local_ties: the tie sensitivity analysis) ====\n")
for (k in CHOICES) { r <- res[res$choice == k, ]
  cat(sprintf("\n[%s] instruments %d, in outcome %d (%d by proxy, %d missing), flipped %d, palindromic dropped %d, mismatch dropped %d, kept %d; harmonised data identical: %s; replays: %s\n",
    k, r$instruments[1], r$in_outcome[1], r$by_proxy[1], r$missing_in_outcome[1], r$flipped[1], r$palindromic_dropped[1], r$mismatch_dropped[1], r$kept[1], r$harmonised_identical[1], paste(unique(r$replay), collapse = "/")))
  cat(sprintf("  %-10s app %8.4f (%.4f) n=%d   ref %8.4f (%.4f) n=%d   diff %.1e\n", r$method, r$app_est, r$app_se, r$app_n, r$ref_est, r$ref_se, r$ref_n, abs(r$app_est - r$ref_est)), sep = "")
  cat(sprintf("  Cochran's Q (IVW): app %.1f on %d df, ref %.1f on %d df\n", r$q_app[1], r$q_df_app[1], r$q_ref[1], r$q_df_ref[1])) }
det <- res$method %in% c("ivw_mre", "egger", "egger_int", "wm", "smed", "mode")
cat(sprintf("\nLargest |difference| in the estimate, all methods: %.2e; in the SE, deterministic methods (IVW, MR-Egger): %.2e\n", max(abs(res$app_est - res$ref_est)[det]),
  max(abs(res$app_se - res$ref_se)[res$method %in% c("ivw_mre", "egger", "egger_int")])))
cat("\n==== Tier B: MR Studio against the paper (2018; about 56 variants, inferred from Q = 122.5 with p = 4.72e-07) ====\n")
for (k in CHOICES) { r <- res[res$choice == k & res$method %in% P$method, ]
  cat(sprintf("\n[%s] Q %.1f on %d df (paper 122.5 on 55 df)\n", k, r$q_app[1], r$q_df_app[1]))
  cat(ifelse(r$method == "egger_int", sprintf("  %-10s paper %.3f (SE %.3f)                 app %.3f (SE %.3f)\n", r$method, r$paper_b, r$paper_se, r$app_est, r$app_se),
    sprintf("  %-10s paper OR %.2f (%.2f to %.2f)   app OR %.2f (%.2f to %.2f)\n", r$method, r$paper_or, r$paper_lo, r$paper_hi, r$app_or, r$app_lo, r$app_hi)), sep = "") }
cat("\nPaper outliers rs11065987, rs1250229, rs4530754 among the instruments:", paste(sapply(names(exd), function(k) sprintf("%s %s", k, paste(c("rs11065987", "rs1250229", "rs4530754") %in% exd[[k]]$snp, collapse = "/"))), collapse = "; "), "\n")
cat("\nWritten:", file.path(if (dry) vout else file.path("validation", "hemani2018"), paste0(c("results", "instruments"), out_tag, ".csv")), "\n")
unlink(tmp, recursive = TRUE)
