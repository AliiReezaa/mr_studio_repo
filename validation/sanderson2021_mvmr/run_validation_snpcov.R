# Validation study 6, part 4 (Round 55): covariances of the variant-exposure associations from individual-level data, as
# MVMR::snpcov_mvmr (Sanderson, Spiller and Bowden 2021, doi:10.1002/sim.9133, PMID 34338327). The MVMR package has no individual-level
# example data, so a cohort is simulated with a fixed seed: 20,000 people, 30 variants (dosages 0 to 2), three exposures measured in the
# same people and correlated through a shared confounder, and an outcome from a separate simulated sample. Run from the repository root:
#   Rscript --vanilla validation/sanderson2021_mvmr/run_validation_snpcov.R
args_dir <- normalizePath(file.path("prototype", "engine")); vdir <- normalizePath(file.path("validation", "sanderson2021_mvmr"))
options(mrgui.engine_dir = args_dir); tmp <- tempfile("mvval-"); root <- file.path(tmp, "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root, MRGUI_CONFIG_DIR = file.path(tmp, "cfg"))
suppressMessages({ source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
  for (f in c("import_core.R", "opengwas_core.R", "ld_core.R", "strobe_core.R", "feasibility_core.R", "report_core.R", "pdf_core.R", "outlier_core.R", "mv_core.R")) source(file.path(args_dir, f)) })
if (!requireNamespace("MVMR", quietly = TRUE)) stop("install the MVMR package to run this validation")

set.seed(55); n <- 20000L; J <- 30L; K <- 3L; maf <- stats::runif(J, 0.1, 0.5)
gen <- function(n) sapply(maf, function(p) stats::rbinom(n, 2, p))
A <- matrix(stats::rnorm(J * K, 0, 0.08), J, K); theta <- c(0.3, -0.2, 0.1)
G <- gen(n); colnames(G) <- sprintf("rs%07d", seq_len(J)); u <- stats::rnorm(n)
X <- G %*% A + cbind(u, 0.6 * u, -0.4 * u) + matrix(stats::rnorm(n * K), n, K); colnames(X) <- c("x_ldl", "x_hdl", "x_tg")
G2 <- gen(n); u2 <- stats::rnorm(n); Y <- (G2 %*% A + cbind(u2, 0.6 * u2, -0.4 * u2) + matrix(stats::rnorm(n * K), n, K)) %*% theta + u2 + stats::rnorm(n)
reg <- function(y, g) summary(stats::lm(y ~ g))$coefficients[2, 1:2]
bxs <- lapply(seq_len(K), function(k) t(vapply(seq_len(J), function(j) reg(X[, k], G[, j]), numeric(2)))); by <- t(vapply(seq_len(J), function(j) reg(Y, G2[, j]), numeric(2)))
f <- file.path(tmp, "cohort.csv"); data.table::fwrite(data.frame(id = seq_len(n), X, G), f)

# 1. MR Studio's covariances, read from the file as the app does (the id column is not a 0-2 dosage, so it is ignored).
cv <- mv_cov_from_file(f, colnames(X)); stopifnot(identical(cv$variant_id, colnames(G)), identical(attr(cv, "meta")$ignored, "id"))
lmse <- sapply(seq_len(K), function(k) bxs[[k]][, 2])
c1 <- max(abs(sqrt(as.matrix(cv[, sprintf("cov_x_%d_%d", 1:3, 1:3)])) / lmse - 1))   # diagonal against lm() standard errors
# 2. Against MVMR::snpcov_mvmr on centred genotypes (where the package formula is the OLS covariance up to n / (n - 2)).
Gc <- scale(G, scale = FALSE); pk <- MVMR::snpcov_mvmr(Gc, X); pkc <- t(vapply(pk, function(m) m[upper.tri(m, diag = TRUE)] * n / (n - 2), numeric(6)))
ours <- as.matrix(cv[, c("cov_x_1_1", "cov_x_1_2", "cov_x_2_2", "cov_x_1_3", "cov_x_2_3", "cov_x_3_3")]); c2 <- max(abs(ours / pkc - 1))
pk_raw <- MVMR::snpcov_mvmr(G, X); c2_raw <- stats::median(sqrt(vapply(pk_raw, function(m) m[1, 1], 0)) / lmse[, 1])   # the package on 0/1/2 dosages
# 3. The conditional F, Q_A and minimum Q with these covariances, against strength_mvmr and pleiotropy_mvmr given the same list.
gl <- lapply(seq_len(J), function(j) { m <- matrix(0, K, K); for (a in 1:K) for (b in 1:K) m[a, b] <- cv[[sprintf("cov_x_%d_%d", min(a, b), max(a, b))]][j]; m })
Fm <- MVMR::format_mvmr(BXGs = sapply(bxs, `[`, , 1), BYG = by[, 1], seBXGs = sapply(bxs, `[`, , 2), seBYG = by[, 2], RSID = colnames(G))
invisible(utils::capture.output({ rf <- unlist(MVMR::strength_mvmr(Fm, gencov = gl)); rq <- MVMR::pleiotropy_mvmr(Fm, gencov = gl); rqm <- MVMR::pleiotropy_mvmr(Fm, gencov = gl, estimator = "qmin") }))
d <- data.frame(variant_id = colnames(G), ea = "A", oa = "G", beta_x_1 = bxs[[1]][, 1], se_x_1 = bxs[[1]][, 2], beta_x_2 = bxs[[2]][, 1], se_x_2 = bxs[[2]][, 2],
                beta_x_3 = bxs[[3]][, 1], se_x_3 = bxs[[3]][, 2], beta_y = by[, 1], se_y = by[, 2], mr_keep = TRUE, stringsAsFactors = FALSE)
inp <- mv_attach_cov(list(kind = "multi_iv", data = d, exposure_id = "LDL-c|HDL-c|TG", outcome_id = "Y", outcome_scale = "beta", exposure_scale = "beta"), cv)
h <- studio_new("snpcov validation"); r <- studio_run(h, "Simulated lipids", "Y", inp, c("mvivw", "mvqhet"), idx = 1L, params = list(ci_level = 0.95, mv_r_1_2 = 0.9))
pl <- pair_payload(h$dir, r$pair_folder); dg <- pl$diagnostics; g <- function(m, nm) Filter(function(x) x$method == m && x$name == nm, dg)
cf <- vapply(g("mvivw", "mv_cond_f"), function(x) x$value, 0); qa <- g("mvivw", "mv_q_a")[[1]]$value; tau2 <- g("mvqhet", "mv_qhet_tau2")[[1]]$value
src <- g("mvivw", "mv_pheno_cor")[[1]]$note
# The package's qmin estimator (pleiotropy_mvmr(estimator = "qmin")) minimises Q_A without the additive tau2 that qhet_mvmr 0.4.9 and
# MR Studio add when the minimum exceeds its degrees of freedom. Here it does (tau2 > 0), so the minimum is compared on a second outcome
# with doubled standard errors, where the minimum is below its degrees of freedom and tau2 = 0: the case both estimators define the same way.
Fm2 <- Fm; Fm2$sebetaYG <- 2 * Fm$sebetaYG; invisible(utils::capture.output(rqm <- MVMR::pleiotropy_mvmr(Fm2, gencov = gl, estimator = "qmin")))
inp2 <- inp; inp2$data$se_y <- 2 * inp2$data$se_y; r2 <- studio_run(h, "Simulated lipids, wider outcome SE", "Y", inp2, "mvqhet", idx = 2L, params = list(ci_level = 0.95))
dg2 <- pair_payload(h$dir, r2$pair_folder)$diagnostics; tau2b <- Filter(function(x) x$name == "mv_qhet_tau2", dg2)[[1]]$value; qmin <- Filter(function(x) x$name == "mv_qhet_q", dg2)[[1]]$value
res <- data.frame(check = c("diagonal vs lm() SE^2 (max relative difference of the SE)", "vs snpcov_mvmr on centred genotypes x n/(n-2) (max relative)",
                            sprintf("conditional F %d vs strength_mvmr (relative)", 1:3), "Q_A vs pleiotropy_mvmr (relative)", "minimum Q vs pleiotropy_mvmr(qmin), outcome SE x 2 (relative; tau2 = 0)"),
                  value = c(c1, c2, abs(cf / rf - 1), abs(qa / rq$Qstat - 1), abs(qmin / rqm$Qstat - 1)), stringsAsFactors = FALSE)
res$pass <- res$value < c(1e-10, 1e-10, 1e-8, 1e-8, 1e-8, 1e-8, 1e-6)
ok <- all(res$pass) && tau2b == 0 && identical(src, "individual-level data") && identical(pl$replay$status, "identical") && all(vapply(r$results, function(x) identical(x$status, "ok"), TRUE))
utils::write.csv(res, file.path(vdir, "results_snpcov.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nsnpcov_mvmr on 0/1/2 dosages: its SE is %.3f x the lm() SE (median, exposure 1): the package uses the uncentred sum of g^2\n", c2_raw))
cat(sprintf("covariance source used: %s (the phenotypic correlation setting 0.9 is ignored); tau2 = %g (doubled outcome SE: %g); replay: %s\nVALIDATION %s\n", src, tau2, tau2b, pl$replay$status, if (ok) "PASSED" else "FAILED"))
if (!ok) quit(status = 1)
