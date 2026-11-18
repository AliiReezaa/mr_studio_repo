# Validation study 7: Carter et al. 2019, BMJ 365:l1855 (doi:10.1136/bmj.l1855, PMID 31122926), two-sample MR mediation of the effect of
# education on coronary heart disease (CHD) by body mass index (BMI). Published (two-sample MR): total effect OR 0.63 (0.60 to 0.67) per SD of
# education (fixed-effect IVW); education on BMI -0.22 (-0.24 to -0.20) SD (IVW); proportion mediated by BMI 18% (14% to 23%) by the product
# of coefficients, with BMI on CHD from regression-based MVMR adjusted for education and delta-method intervals.
# Data, from OpenGWAS: education, Lee et al. 2018 without 23andMe (ieu-a-1239); BMI, Yengo et al. 2018 (ieu-b-40); CHD, CARDIoGRAMplusC4D
# 1000 Genomes (ieu-a-7). Runs only when asked for explicitly, with a valid token:
#   MRSTUDIO_VALIDATE_OPENGWAS=1 OPENGWAS_JWT=<token> Rscript --vanilla validation/carter2019_mediation/run_validation.R
# Add MRSTUDIO_VALIDATION_CACHE=<folder> to keep what OpenGWAS returns and rerun without querying again.
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "carter2019_mediation"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("medval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R", "mediation_core.R")) source(file.path(args_dir, f)) })
jwt <- Sys.getenv("OPENGWAS_JWT")
if (!identical(Sys.getenv("MRSTUDIO_VALIDATE_OPENGWAS"), "1") || !nzchar(jwt)) { cat("Skipped: set MRSTUDIO_VALIDATE_OPENGWAS=1 and OPENGWAS_JWT (a valid OpenGWAS token).\n"); quit(status = 0) }
options(ieugwasr_api = "https://api.opengwas.io/api/"); Sys.setenv(OPENGWAS_JWT = jwt)
EDU <- "ieu-a-1239"; BMI <- "ieu-b-40"; CHD <- "ieu-a-7"
# Optional cache of what OpenGWAS returns (MRSTUDIO_VALIDATION_CACHE=<folder>): each query is saved once and read back on later runs, so a
# rerun after a time-out or with the allowance used up does not query again. Summary data only; the token is never written.
cache_dir <- Sys.getenv("MRSTUDIO_VALIDATION_CACHE")
cached <- function(key, expr) {
  if (!nzchar(cache_dir)) return(expr)
  f <- file.path(cache_dir, paste0("carter2019_", gsub("[^A-Za-z0-9._-]", "_", key), ".rds"))
  if (file.exists(f)) return(readRDS(f))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE); v <- expr; saveRDS(v, f); v
}
ext_ins <- function(id, r2 = 0.001) cached(sprintf("ins_%s_r2_%g", id, r2), suppressMessages(TwoSampleMR::extract_instruments(id, r2 = r2)))
ext_out <- function(snps, out) cached(sprintf("out_%s_%s", out, digest::digest(sort(snps))), suppressMessages(TwoSampleMR::extract_outcome_data(snps, out, proxies = TRUE)))
paper <- list(total_or = c(0.63, 0.60, 0.67), edu_bmi = c(-0.22, -0.24, -0.20), prop = c(0.18, 0.14, 0.23))

# Instruments, as in the paper: genome-wide significant lead variants after clumping (r2 < 0.001, 10,000 kb, 1000 Genomes EUR) with the
# OpenGWAS clumping service, the server-side equivalent of TwoSampleMR::clump_data defaults. The paper used the 1271 education lead variants of
# the full Lee et al. sample (with 23andMe, r < 0.1); the public summary statistics exclude 23andMe, so fewer education variants reach P < 5e-8.
ie <- ext_ins(EDU)
pair <- function(ins, out) { o <- ext_out(ins$SNP, out)
  H <- suppressMessages(TwoSampleMR::harmonise_data(ins, o, action = 2)); H <- H[H$mr_keep, ]
  d <- data.frame(variant_id = H$SNP, ea = H$effect_allele.exposure, oa = H$other_allele.exposure, beta_x = H$beta.exposure, se_x = H$se.exposure, eaf_x = H$eaf.exposure,
                  beta_y = H$beta.outcome, se_y = H$se.outcome, eaf_y = H$eaf.outcome, p_x = H$pval.exposure, p_y = H$pval.outcome, n_x = H$samplesize.exposure, n_y = H$samplesize.outcome,
                  harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  list(kind = "pair_iv", data = d, exposure_id = H$id.exposure[1], outcome_id = out, outcome_scale = if (out == CHD) "log_or" else "beta", exposure_scale = "beta") }
in_tot <- pair(ie, CHD); in_xm <- pair(ie, BMI)
# MVMR of education and BMI on CHD: union of both instrument sets clumped together, both exposures looked up at every variant (TwoSampleMR's
# multivariable extraction, as the paper's "regression-based multivariable MR").
mx <- cached("mv_exposures", suppressMessages(TwoSampleMR::mv_extract_exposures(c(EDU, BMI)))); mo <- ext_out(mx$SNP, CHD); mh <- suppressMessages(TwoSampleMR::mv_harmonise_data(mx, mo))
ord <- match(c(EDU, BMI), colnames(mh$exposure_beta))
dm <- data.frame(variant_id = rownames(mh$exposure_beta), ea = "A", oa = "G", beta_x_1 = mh$exposure_beta[, ord[1]], se_x_1 = mh$exposure_se[, ord[1]], beta_x_2 = mh$exposure_beta[, ord[2]], se_x_2 = mh$exposure_se[, ord[2]],
                 beta_y = mh$outcome_beta, se_y = mh$outcome_se, mr_keep = TRUE, stringsAsFactors = FALSE)
in_mv <- list(kind = "multi_iv", data = dm, exposure_id = paste(EDU, BMI, sep = "|"), outcome_id = CHD, outcome_scale = "log_or", exposure_scale = "beta")

h <- studio_new("Carter 2019 mediation"); med <- list(x = EDU, ms = BMI, y = CHD)
r1 <- studio_run(h, "Education", "CHD", in_tot, "ivw_mre", idx = 1L, params = list(ci_level = 0.95, ivw_model = "fixed"))
r2 <- studio_run(h, "Education", "BMI", in_xm, "ivw_mre", idx = 2L, params = list(ci_level = 0.95, ivw_model = "fixed"))
r3 <- studio_run(h, "Education + BMI", "CHD", in_mv, "mvivw", idx = 3L, params = list(ci_level = 0.95, mv_model = "default"))
prs <- list(list(med = c(list(role = "total"), med), folder = r1$pair_folder, exName = "Education", outName = "CHD"),
            list(med = c(list(role = "xm", m = BMI), med), folder = r2$pair_folder, exName = "Education", outName = "BMI"),
            list(med = c(list(role = "direct", m = BMI), med), folder = r3$pair_folder, exName = "Education + BMI", outName = "CHD", exNames = list("Education", "BMI")))
res <- mediation_run(h, prs); g <- function(me, q) res[res$method == me & res$quantity == q, ][1, ]
tt <- res[res$quantity == "total", ]; xm <- g("two-step", "exposure on mediator"); pp <- g("product", "proportion mediated")
out <- data.frame(quantity = c("total effect OR", "education on BMI (SD)", "proportion mediated by BMI (product)"),
                  paper = c(paper$total_or[1], paper$edu_bmi[1], paper$prop[1]), paper_lo = c(paper$total_or[2], paper$edu_bmi[2], paper$prop[2]), paper_hi = c(paper$total_or[3], paper$edu_bmi[3], paper$prop[3]),
                  mr_studio = c(exp(tt$estimate), xm$estimate, pp$estimate), lo = c(exp(tt$ci_lower), xm$ci_lower, pp$ci_lower), hi = c(exp(tt$ci_upper), xm$ci_upper, pp$ci_upper))
out$within_paper_ci <- out$mr_studio >= out$paper_lo & out$mr_studio <= out$paper_hi
counts <- c(education_instruments = nrow(in_tot$data), education_in_bmi = nrow(in_xm$data), mvmr_variants = nrow(dm), paper_education_instruments = 1271L)
utils::write.csv(out, file.path(vdir, "results.csv"), row.names = FALSE); utils::write.csv(res, file.path(vdir, "mediation_results.csv"), row.names = FALSE)
print(out, digits = 3); print(counts)
reps <- vapply(list(r1, r2, r3), function(r) pair_payload(h$dir, r$pair_folder)$replay$status, "")

# 1. MR Studio against the packages on the same data: fixed-effect IVW (TwoSampleMR::mr_ivw_fe) for the total effect and education on BMI,
#    MVMR-IVW (MendelianRandomization::mr_mvivw) for the direct effects, and the mediation formulas recomputed here.
fe <- function(d) TwoSampleMR::mr_ivw_fe(d$beta_x, d$beta_y, d$se_x, d$se_y); ft <- fe(in_tot$data); fa <- fe(in_xm$data)
mv <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = cbind(dm$beta_x_1, dm$beta_x_2), bxse = cbind(dm$se_x_1, dm$se_x_2), by = dm$beta_y, byse = dm$se_y))
bb <- g("two-step", "mediator on outcome given exposure"); dd <- g("MVMR", "direct"); ip <- g("product", "indirect")
eq <- c(total = abs(tt$estimate - ft$b) + abs(tt$se - ft$se), edu_bmi = abs(xm$estimate - fa$b) + abs(xm$se - fa$se),
        mvmr = max(abs(c(dd$estimate, bb$estimate) - mv@Estimate), abs(c(dd$se, bb$se) - mv@StdError)),
        product = abs(ip$estimate - fa$b * mv@Estimate[2]) + abs(ip$se - sqrt(fa$b^2 * mv@StdError[2]^2 + mv@Estimate[2]^2 * fa$se^2)),
        proportion = abs(pp$estimate - fa$b * mv@Estimate[2] / ft$b))
# 2. How much education on BMI depends on the instrument set. The paper used the 1271 lead variants of the full Lee et al. sample (with
#    23andMe, independent at r < 0.1); the public statistics give fewer. Measured here: (a) the estimate with instruments clumped at r2 < 0.1,
#    closer to the paper's definition; (b) the spread of the estimate over 500 random subsets, each the size of the main set,
#    of the r2 < 0.1 instruments. Criterion for education on BMI: within 3 standard deviations of that spread of the paper's value.
big <- pair(ext_ins(EDU, 0.1), BMI)$data; set.seed(2019)   # one extra query only (OpenGWAS allowance)
sens <- rbind(c(r2 = 0.001, n = nrow(in_xm$data), b = fa$b, se = fa$se), c(r2 = 0.1, n = nrow(big), b = fe(big)$b, se = fe(big)$se))
sub <- vapply(1:500, function(i) fe(big[sample(nrow(big), nrow(in_xm$data)), ])$b, 0); sd_sub <- stats::sd(sub)
gap <- abs(xm$estimate - paper$edu_bmi[1]); within_sub <- gap <= 3 * sd_sub
print(round(eq, 14)); print(sens, digits = 3); cat(sprintf("education on BMI: gap to the paper %.3f; SD over instrument subsets %.4f (3 SD = %.3f)\n", gap, sd_sub, 3 * sd_sub))
utils::write.csv(data.frame(check = names(eq), abs_difference = eq), file.path(vdir, "package_equivalence.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(sens), file.path(vdir, "instrument_sensitivity.csv"), row.names = FALSE)
ok <- all(eq < 1e-10) && out$within_paper_ci[1] && out$within_paper_ci[3] && within_sub && all(reps == "identical")
cat(sprintf("\nreplay: %s\nVALIDATION %s\n", paste(reps, collapse = ", "), if (ok) "PASSED" else "FAILED (see README for the comparison)"))
