# Validation study 5: Rees, Wood and Burgess 2017, multivariable MR-Egger (Stat Med 36:4705-4718, doi:10.1002/sim.7492, PMID 28960498),
# Table 2: MVMR-IVW and MVMR-Egger (three orientations) for HDL-c, LDL-c and triglycerides on CHD with the 185 lipid variants of
# Do et al. 2013. The lipid associations are in validation/bowden2016/dodata.csv (figshare doi:10.6084/m9.figshare.1116328).
# The paper took CHD associations from the 2015 CARDIoGRAMplusC4D release; dodata.csv holds the earlier CARDIoGRAM associations.
# Part A runs the whole MR Studio path on dodata.csv (import, MVMR data assembly, generated scripts, replay check) and checks that it
# equals MendelianRandomization exactly. Part B replaces the CHD associations with the 2015 release (OpenGWAS ieu-a-7, Nikpay et al.)
# and compares with Table 2; it needs an OpenGWAS token in the environment variable OPENGWAS_JWT and is skipped without one.
#   Rscript --vanilla validation/rees2017_mvmr/run_validation.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "rees2017_mvmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("reesval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })

# Table 2 of the paper (log odds ratio of CHD per SD increase, 95% CI), transcribed from the open-access full text (PMC5725762).
T2 <- utils::read.csv(text = "analysis,lipid,est,lo,hi
mvivw,HDL-c,-0.039,-0.123,0.045
mvivw,LDL-c,0.375,0.292,0.457
mvivw,Triglycerides,0.173,0.063,0.283
egger_hdl,HDL-c,0.036,-0.063,0.134
egger_hdl,LDL-c,0.378,0.297,0.458
egger_hdl,Triglycerides,0.136,0.024,0.247
egger_hdl,Intercept,-0.005,-0.009,-0.001
egger_ldl,HDL-c,-0.034,-0.118,0.049
egger_ldl,LDL-c,0.420,0.318,0.522
egger_ldl,Triglycerides,0.194,0.081,0.308
egger_ldl,Intercept,-0.003,-0.007,0.001
egger_tg,HDL-c,-0.018,-0.102,0.066
egger_tg,LDL-c,0.350,0.267,0.433
egger_tg,Triglycerides,0.083,-0.045,0.211
egger_tg,Intercept,0.005,0.001,0.009", stringsAsFactors = FALSE)

r <- utils::read.csv(file.path(normalizePath(file.path("validation", "bowden2016")), "dodata.csv"), colClasses = "character")
logp <- function(x) { x <- trimws(x); m <- regmatches(x, regexec("^([0-9.]+)X10-([0-9]+)$", x, ignore.case = TRUE))
  vapply(seq_along(x), function(i) if (length(m[[i]]) == 3) log(as.numeric(m[[i]][2])) - as.numeric(m[[i]][3]) * log(10) else log(as.numeric(x[i])), 0) }
d <- data.frame(rsid = r$rsid, a1 = r$a1, a2 = r$a2, stringsAsFactors = FALSE)
for (k in c("ldlc", "hdlc", "tg", "chd")) { b <- as.numeric(r[[paste0(k, "beta")]]); lp <- logp(r[[paste0(k, "p")]])   # SE = |beta| / z(P/2), as in validation/bowden2016
  d[[paste0(k, "_b")]] <- b; d[[paste0(k, "_p")]] <- exp(lp); d[[paste0(k, "_se")]] <- abs(b) / -stats::qnorm(lp - log(2), log.p = TRUE) }
mp <- list(snp = "snp", ea = "ea", oa = "oa", beta = "beta", se = "se", p = "p")
fct <- function(nm, role, kind = "continuous", extra = list()) utils::modifyList(list(name = nm, role = role, kind = kind, scale = if (kind == "binary") "logor" else "beta", ea_ok = "yes", frq = "none", strand = "forward",
  build = "37", ancestry = "EUR", n = "", scope = "selected", sex = "both", ukb = FALSE, independent = "yes"), extra)
imp <- function(df, nm, role, kind = "continuous", extra = list()) { f <- tempfile(fileext = ".csv"); utils::write.csv(df, f, row.names = FALSE)
  v <- import_validate(f, ",", names(df), mp, fct(nm, role, kind, extra)); if (!isTRUE(v$ok)) stop("import refused: ", nm)
  invisible(import_save(f, basename(f), ",", names(df), mp, fct(nm, role, kind, extra))) }
# All four columns come from one table with the same effect allele, so the forward-strand statement is true and palindromic variants are
# kept (as in the paper, which used all variants). The lipids are imported with full summary statistics for all 185 variants ("either").
for (k in c("hdlc", "ldlc", "tg")) { s <- is.finite(d[[paste0(k, "_se")]]) & d[[paste0(k, "_se")]] > 0
  imp(data.frame(snp = d$rsid[s], ea = d$a1[s], oa = d$a2[s], beta = d[[paste0(k, "_b")]][s], se = d[[paste0(k, "_se")]][s], p = d[[paste0(k, "_p")]][s]), k, "either") }
oky <- is.finite(d$chd_se) & d$chd_se > 0
imp(data.frame(snp = d$rsid[oky], ea = d$a1[oky], oa = d$a2[oky], beta = d$chd_b[oky], se = d$chd_se[oky], p = d$chd_p[oky]), "CHD CARDIoGRAM", "outcome", "binary", list(ncase = "22233", ncontrol = "64762"))

# Instruments: every variant associated with any lipid. The paper used the 185 variants as given (they were selected by Do et al. as
# independent), so the joint clumping step is replaced by keeping every variant (an LD panel would drop none of these by design).
thr <- 5e-8; no_clump <- function(snp, p, pop) data.frame(variant_id = snp, status = "index", stringsAsFactors = FALSE)
lips <- c("lib:hdlc", "lib:ldlc", "lib:tg"); labs <- c("HDL-c", "LDL-c", "Triglycerides")
mv <- library_mv_input(lips, "lib:CHD CARDIoGRAM", pthr = thr, clump = no_clump)
cat(sprintf("Instruments: %d in the union at P < %g (paper: 185); %d with all three lipids and CHD.\n", mv$ledger$mv$union, thr, mv$ledger$mv$final))

run_set <- function(inp, tag) {
  h <- studio_new(paste("Rees 2017", tag)); out <- list(); k <- 0L
  for (an in c("mvivw", "egger_hdl", "egger_ldl", "egger_tg")) { k <- k + 1L
    m <- if (an == "mvivw") "mvivw" else "mvegger"; orient <- match(sub("egger_", "", an), c("hdl", "ldl", "tg"))
    rr <- studio_run(h, paste0("Lipids ", an), "CHD", inp, m, idx = k, params = list(ci_level = 0.95, mv_model = "random", mv_orient = if (is.na(orient)) 1L else orient))
    pl <- pair_payload(h$dir, rr$pair_folder); e <- .est_df(pl)
    dd <- inp$data[inp$data$mr_keep, ]; mi <- MendelianRandomization::mr_mvinput(bx = as.matrix(dd[, sprintf("beta_x_%d", 1:3)]), bxse = as.matrix(dd[, sprintf("se_x_%d", 1:3)]), by = dd$beta_y, byse = dd$se_y)
    ref <- if (m == "mvivw") { f <- MendelianRandomization::mr_mvivw(mi, model = "random"); c(f@Estimate) } else { f <- MendelianRandomization::mr_mvegger(mi, orientate = orient); c(f@Estimate, f@Intercept) }
    e$lipid <- ifelse(e$estimand == "pleiotropy_intercept", "Intercept", labs[match(e$exposure_id, lips)])
    out[[k]] <- data.frame(data = tag, analysis = an, lipid = e$lipid, est = e$estimate, lo = e$ci_lower, hi = e$ci_upper, n = e$n_ivs, vs_package = max(abs(e$estimate - ref)), replay = pl$replay$status, stringsAsFactors = FALSE)
  }
  x <- do.call(rbind, out); y <- merge(x, T2, by = c("analysis", "lipid"), suffixes = c("", "_paper")); y$diff <- y$est - y$est_paper; y[order(y$analysis, y$lipid), ]
}
A <- run_set(mv$input, "dodata CARDIoGRAM")
okA <- all(A$replay == "identical") && all(A$vs_package < 1e-12)
cat("\nPart A (MR Studio path on dodata.csv): replay identical and equal to MendelianRandomization:", okA, "\n")
print(A[, c("analysis", "lipid", "est", "est_paper", "diff", "n")], digits = 3)

B <- NULL; jwt <- Sys.getenv("OPENGWAS_JWT")
# Part B runs only when asked for explicitly (MRSTUDIO_VALIDATE_OPENGWAS=1), so that a token stored in the environment is never sent
# to OpenGWAS without the analyst choosing to.
if (identical(Sys.getenv("MRSTUDIO_VALIDATE_OPENGWAS"), "1") && nzchar(jwt)) {   # Part B: CHD from the 2015 CARDIoGRAMplusC4D release, as in the paper
  y <- as.data.frame(ieugwasr::associations(d$rsid, "ieu-a-7", proxies = 0, opengwas_jwt = jwt), stringsAsFactors = FALSE)
  y <- y[match(d$rsid, y$rsid), ]; ok <- !is.na(y$rsid)
  cat("Variants of the paper missing from ieu-a-7:", paste(d$rsid[!ok], collapse = ", "), "\n")
  imp(data.frame(snp = y$rsid[ok], ea = y$ea[ok], oa = y$nea[ok], beta = as.numeric(y$beta[ok]), se = as.numeric(y$se[ok]), p = as.numeric(y$p[ok])), "CHD CARDIoGRAMplusC4D 2015", "outcome", "binary", list(ncase = "60801", ncontrol = "123504"))
  mv2 <- library_mv_input(lips, "lib:CHD CARDIoGRAMplusC4D 2015", pthr = thr, clump = no_clump)
  B <- run_set(mv2$input, "ieu-a-7 CARDIoGRAMplusC4D 2015")
  cat("\nPart B (CHD from the 2015 release, as in the paper):\n"); print(B[, c("analysis", "lipid", "est", "est_paper", "diff", "n")], digits = 3)
} else cat("\nPart B skipped: set MRSTUDIO_VALIDATE_OPENGWAS=1 and OPENGWAS_JWT (a valid OpenGWAS token) to compare with the 2015 CHD release used in the paper.\n")
utils::write.csv(rbind(A, B), file.path(vdir, "results.csv"), row.names = FALSE)
# Criterion for Part B. Three of the paper's 185 variants (rs1998013, rs2247056, rs7422339) are absent from the OpenGWAS copy of the 2015
# release, with or without proxies. Leaving three variants out at random changes MVMR-IVW by up to 0.014, 0.015 and 0.022 (95th
# percentile of 2,000 draws; see README), so estimates must agree within 0.025 and the MVMR-Egger intercepts within 0.001.
okB <- is.null(B) || (all(B$replay == "identical") && all(B$vs_package < 1e-12) && all(abs(B$diff[B$lipid != "Intercept"]) <= 0.025) && all(abs(B$diff[B$lipid == "Intercept"]) <= 0.001))
cat(if (okA && okB) "VALIDATION PASSED" else "VALIDATION FAILED", if (is.null(B)) "(Part A only)" else "(Parts A and B)", "\n")
unlink(tmp, recursive = TRUE); if (!(okA && okB)) quit(status = 1)
