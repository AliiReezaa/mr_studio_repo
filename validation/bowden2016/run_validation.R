# Validation study 1: Bowden, Davey Smith, Haycock and Burgess 2016, Table 4 (lipid fractions on coronary artery disease).
# Paper: Genet Epidemiol 40:304-14, doi:10.1002/gepi.21965 (PMID 27061298). Data: Burgess, figshare doi:10.6084/m9.figshare.1116328
# (dodata.csv, CC BY 4.0, MD5 49485561fb97887a16805d1a1652b21c), from Supplementary Table 1 of Do et al. 2013, doi:10.1038/ng.2795.
# Runs the data through MR Studio's own engine (import, pairing and harmonisation, generated scripts, replay check) in a
# throw-away studies folder, and also through TwoSampleMR directly with every variant kept, then compares both with Table 4.
# Run from the repository root:  Rscript --vanilla validation/bowden2016/run_validation.R
repo <- normalizePath("."); vdir <- file.path(repo, "validation", "bowden2016"); args_dir <- file.path(repo, "prototype", "engine")
stopifnot(identical(unname(tools::md5sum(file.path(vdir, "dodata.csv"))), "49485561fb97887a16805d1a1652b21c"))
tmp <- tempfile("bowden-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R")) source(file.path(args_dir, f)) })

# The file gives betas and P values (written like 4X10-19), not standard errors. SE = |beta| / z(P/2), computed on the log scale so
# that P values below the smallest double (for example 1X10-292 is fine, but smaller would not be) do not become 0.
r <- utils::read.csv(file.path(vdir, "dodata.csv"), colClasses = "character")
logp <- function(x) { x <- trimws(x); m <- regmatches(x, regexec("^([0-9.]+)X10-([0-9]+)$", x, ignore.case = TRUE))
  vapply(seq_along(x), function(i) if (length(m[[i]]) == 3) log(as.numeric(m[[i]][2])) - as.numeric(m[[i]][3]) * log(10) else log(as.numeric(x[i])), 0) }
pnum <- function(x) as.numeric(sub("X10-", "e-", trimws(x), ignore.case = TRUE))
d <- data.frame(rsid = r$rsid, a1 = r$a1, a2 = r$a2, stringsAsFactors = FALSE)
for (k in c("ldlc", "hdlc", "tg", "chd")) { b <- as.numeric(r[[paste0(k, "beta")]]); lp <- logp(r[[paste0(k, "p")]])
  d[[paste0(k, "_b")]] <- b; d[[paste0(k, "_lp")]] <- lp; d[[paste0(k, "_p")]] <- pnum(r[[paste0(k, "p")]]); d[[paste0(k, "_se")]] <- abs(b) / -stats::qnorm(lp - log(2), log.p = TRUE) }
# A beta of exactly 0 with P = 1 has no recoverable SE: rs4942486 for CAD (an LDL instrument) and rs7422339 for triglycerides (not a TG instrument).
lips <- c("ldlc", "hdlc", "tg"); THR <- 1e-8   # Web Appendix 4 counts (73, 85, 47 and 61, 69, 31) are reproduced with P < 1e-8
sel <- function(ex, set) { oth <- setdiff(lips, ex); s <- is.finite(d[[paste0(ex, "_se")]]) & d[[paste0(ex, "_se")]] > 0
  if (set == "primary") s <- s & d[[paste0(ex, "_lp")]] < d[[paste0(oth[1], "_lp")]] & d[[paste0(ex, "_lp")]] < d[[paste0(oth[2], "_lp")]]; s }

mp <- list(snp = "snp", ea = "ea", oa = "oa", beta = "beta", se = "se", p = "p")
fct <- function(nm, role, kind = "continuous", extra = list()) utils::modifyList(list(name = nm, role = role, kind = kind, scale = if (kind == "binary") "logor" else "beta", ea_ok = "yes", frq = "none", strand = "unknown",
  build = "37", ancestry = "EUR", n = "", scope = "selected", sex = "both", ukb = FALSE, independent = "yes"), extra)
imp <- function(df, nm, role, kind = "continuous", extra = list()) { f <- tempfile(fileext = ".csv"); utils::write.csv(df, f, row.names = FALSE)
  v <- import_validate(f, ",", names(df), mp, fct(nm, role, kind, extra)); if (!isTRUE(v$ok)) stop("import refused: ", nm)
  invisible(import_save(f, basename(f), ",", names(df), mp, fct(nm, role, kind, extra))) }
oky <- is.finite(d$chd_se) & d$chd_se > 0
# Each trait is imported twice: strand not declared (the default), and declared forward strand (Round 26), which keeps palindromic
# variants. All four columns come from one table with the same effect allele, so the forward-strand statement is true here.
for (st in c("unknown", "forward")) { sfx <- if (st == "forward") " fwd" else ""
  imp(data.frame(snp = d$rsid[oky], ea = d$a1[oky], oa = d$a2[oky], beta = d$chd_b[oky], se = d$chd_se[oky], p = d$chd_p[oky]), paste0("CAD CARDIoGRAM", sfx), "outcome", "binary", list(ncase = "22233", ncontrol = "64762", strand = st))
  for (ex in lips) for (set in c("all", "primary")) { s <- sel(ex, set)
    imp(data.frame(snp = d$rsid[s], ea = d$a1[s], oa = d$a2[s], beta = d[[paste0(ex, "_b")]][s], se = d[[paste0(ex, "_se")]][s], p = d[[paste0(ex, "_p")]][s]), paste0(ex, " ", set, sfx), "exposure", extra = list(strand = st)) } }

# Table 4 of the paper (estimate, SE), transcribed from the open-access full text (PMC4849733).
T4 <- utils::read.csv(text = "exposure,set,method,est,se
ldlc,all,ivw,0.482,0.060
ldlc,all,smed,0.429,0.070
ldlc,all,wm,0.458,0.065
ldlc,all,egger,0.617,0.103
ldlc,all,egger_int,-0.009,0.005
ldlc,primary,ivw,0.470,0.055
ldlc,primary,smed,0.429,0.079
ldlc,primary,wm,0.457,0.065
ldlc,primary,egger,0.562,0.094
ldlc,primary,egger_int,-0.006,0.005
hdlc,all,ivw,-0.254,0.070
hdlc,all,smed,-0.267,0.090
hdlc,all,wm,-0.069,0.071
hdlc,all,egger,-0.013,0.115
hdlc,all,egger_int,-0.014,0.005
hdlc,primary,ivw,-0.137,0.066
hdlc,primary,smed,-0.224,0.085
hdlc,primary,wm,-0.066,0.065
hdlc,primary,egger,0.092,0.107
hdlc,primary,egger_int,-0.013,0.005
tg,all,ivw,0.416,0.081
tg,all,smed,0.512,0.101
tg,all,wm,0.516,0.084
tg,all,egger,0.422,0.140
tg,all,egger_int,-0.000,0.008
tg,primary,ivw,0.417,0.095
tg,primary,smed,0.565,0.105
tg,primary,wm,0.521,0.087
tg,primary,egger,0.464,0.155
tg,primary,egger_int,-0.004,0.009")

h <- studio_new("Bowden2016 Table4"); i <- 0L; rows <- list()
key <- function(m, e) if (m == "egger" && e == "pleiotropy_intercept") "egger_int" else if (m == "ivw_mre") "ivw" else m
for (ex in lips) for (set in c("all", "primary")) {
  i <- i + 1L; lp <- pair_input(paste0("lib:", ex, " ", set), "lib:CAD CARDIoGRAM", pthr = THR)
  rr <- studio_run(h, lp$exposure_label, lp$outcome_label, lp$input, c("ivw_mre", "egger", "wm", "smed"), idx = i)
  pl <- pair_payload(h$dir, rr$pair_folder); led <- lp$ledger
  i <- i + 1L; lf <- pair_input(paste0("lib:", ex, " ", set, " fwd"), "lib:CAD CARDIoGRAM fwd", pthr = THR)
  rf <- studio_run(h, lf$exposure_label, lf$outcome_label, lf$input, c("ivw_mre", "egger", "wm", "smed"), idx = i)
  pf <- pair_payload(h$dir, rf$pair_folder); fe <- function(k) { e <- Filter(function(e) key(e$method, e$estimand) == k, pf$estimates)[[1]]; c(e$estimate, e$se, e$n_ivs_used) }
  # Every variant kept (no harmonisation: both columns come from one table with the same effect allele), through TwoSampleMR directly.
  s <- sel(ex, set) & d[[paste0(ex, "_lp")]] < log(THR) & oky
  bx <- d[[paste0(ex, "_b")]][s]; sx <- d[[paste0(ex, "_se")]][s]; by <- d$chd_b[s]; sy <- d$chd_se[s]; pa <- TwoSampleMR::default_parameters()
  iv <- TwoSampleMR::mr_ivw(bx, by, sx, sy, pa); eg <- TwoSampleMR::mr_egger_regression(bx, by, sx, sy, pa)
  set.seed(1); wm <- TwoSampleMR::mr_weighted_median(bx, by, sx, sy, pa); set.seed(1); sm <- TwoSampleMR::mr_simple_median(bx, by, sx, sy, pa)
  allk <- list(ivw = c(iv$b, iv$se), egger = c(eg$b, eg$se), egger_int = c(eg$b_i, eg$se_i), wm = c(wm$b, wm$se), smed = c(sm$b, sm$se))
  for (e in pl$estimates) { k <- key(e$method, e$estimand); t4 <- T4[T4$exposure == ex & T4$set == set & T4$method == k, ]
    rows[[length(rows) + 1L]] <- data.frame(exposure = ex, set = set, method = k, paper_est = t4$est, paper_se = t4$se,
      app_est = round(e$estimate, 4), app_se = round(e$se, 4), app_n = e$n_ivs_used,
      fwd_est = round(fe(k)[1], 4), fwd_se = round(fe(k)[2], 4), fwd_n = fe(k)[3], fwd_palindromic_kept = lf$ledger$kept_palindromic, fwd_replay = pf$replay$status, all_kept_est = round(allk[[k]][1], 4), all_kept_se = round(allk[[k]][2], 4), all_kept_n = sum(s),
      candidates = led$instruments, missing_in_outcome = led$missing_in_outcome, palindromic_dropped = led$dropped_palindromic, replay = pl$replay$status) }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(vdir, "results.csv"), row.names = FALSE)
f3 <- function(e, s) sprintf("%6.3f (%.3f)", e, s)
cat(sprintf("%-5s %-8s %-10s  paper %s   app %s n=%d (pal.dropped %d)   app forward strand %s n=%d (pal.kept %d)   all kept %s n=%d   replays %s/%s\n", res$exposure, res$set, res$method,
  f3(res$paper_est, res$paper_se), f3(res$app_est, res$app_se), res$app_n, res$palindromic_dropped, f3(res$fwd_est, res$fwd_se), res$fwd_n, res$fwd_palindromic_kept,
  f3(res$all_kept_est, res$all_kept_se), res$all_kept_n, res$replay, res$fwd_replay), sep = "")
det <- res$method %in% c("ivw", "egger", "egger_int")
cat(sprintf("\nForward strand declared, deterministic methods: largest difference from TwoSampleMR with every variant kept %.2e; from Table 4 %.4f (estimate), %.4f (SE)\n",
  max(abs(res$fwd_est - res$all_kept_est)[det]), max(abs(res$fwd_est - res$paper_est)[det]), max(abs(res$fwd_se - res$paper_se)[det])))
cat(sprintf("\nLargest difference from Table 4, deterministic methods, every variant kept: estimate %.4f, SE %.4f\n",
  max(abs(res$all_kept_est - res$paper_est)[det]), max(abs(res$all_kept_se - res$paper_se)[det])))
