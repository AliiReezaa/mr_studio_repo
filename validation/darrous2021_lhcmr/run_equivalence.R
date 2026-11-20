# Round 67: the adaptive, resumable runner equals lhcMR::lhc_mr() exactly on the same input (synthetic data: 3,000 variants,
# 3 starting points, 4 jackknife blocks), and a second run reads its saved fits instead of computing them. About 2.5 hours at low
# priority on a 2019 Intel Mac (each fit takes about 20 minutes whatever the data size). Run from the repository folder:
#   Rscript --vanilla validation/darrous2021_lhcmr/run_equivalence.R
suppressPackageStartupMessages(library(lhcMR)); source(file.path("prototype", "engine", "resource_core.R"))
eval(parse(text = LHCMR_RUNNER_SRC <- local({ e <- new.env(); sys.source(file.path("prototype", "engine", "lhcmr_runner.R"), e); e$LHCMR_RUNNER_SRC })))
set.seed(67); n <- 3000; z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .2, .2, 1), 2)) * 1.3
df <- data.frame(TSTAT.x = z[, 1], TSTAT.y = z[, 2], N.x = 3e5, N.y = 3e5, M_LOCAL = 5001, LDSC = runif(n, 5, 80), WEIGHT = runif(n, .1, 1), PIK = runif(n, .01, .3), SIGK = runif(n, .2, .9))
sp <- list(iX = 1.05, iY = 1.03, piX = 1e-4, piY = 1e-4, input.df_filtered = df,
           sp_mat = cbind(sp_h2X = c(.2, .3, .25), sp_h2Y = c(.3, .2, .4), sp_tX = c(.1, .2, .3), sp_tY = c(-.1, .1, .2), sp_axy = c(.05, 0, .1), sp_ayx = c(0, .05, -.05), sp_iXY = c(.1, .12, .08)))
w <- tempfile("lhceq-"); d1 <- file.path(w, "pkg"); d2 <- file.path(w, "run"); unlink(c(d1, d2), recursive = TRUE); dir.create(d1); dir.create(d2)
t0 <- Sys.time(); setwd(d1); invisible(capture.output(lhc_mr(sp, c("X", "Y"), paral_method = "lapply", nCores = 1, nBlock = 4))); t1 <- Sys.time()
setwd(d2); invisible(capture.output(lhcmr_adaptive(sp, c("X", "Y"), nBlock = 4L, ck = "units", max_workers = 2L, per_worker_gb = 0.5, mem = mem_info, poll = 1))); t2 <- Sys.time()
a <- read.csv(file.path(d1, "SummarisedResults_X-Y.csv"), row.names = 1); b <- read.csv(file.path(d2, "SummarisedResults_X-Y.csv"), row.names = 1)
cat("package", format(t1 - t0), "runner", format(t2 - t1), "\nmax difference in SummarisedResults:", max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE), "\n")
for (f in c("FullRes_X-Y.csv", "FullRes_X-Y_4blockJK.csv", "4blockJK_summary_X-Y.csv")) cat(f, identical(read.csv(file.path(d1, f)), read.csv(file.path(d2, f))), "\n")
setwd(d2); invisible(capture.output(lhcmr_adaptive(sp, c("X", "Y"), nBlock = 4L, ck = "units", max_workers = 2L, per_worker_gb = 0.5, mem = mem_info, poll = 1))); cat("resume from checkpoints:", format(Sys.time() - t2), "\n"); print(a, digits = 4)
