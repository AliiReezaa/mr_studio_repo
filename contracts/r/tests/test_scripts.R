# Tests for contracts v1.1 (draft): script emission, running the generated script, replay check.
# Run from the contracts/ folder:   Rscript r/tests/test_scripts.R <output_dir>
# Needs the TwoSampleMR package at the pinned commit (the generated scripts call it directly).
args <- commandArgs(trailingOnly = TRUE)
out_dir <- tempfile("run11-", tmpdir = if (length(args)) args[1] else tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c("r/canonical.R", "r/adapter.R", "r/security.R", "r/adapters/reference_adapters.R",
            "r/script.R", "r/adapters/reference_adapters_v11.R")) source(f)

n_pass <- 0L; n_fail <- 0L
expect <- function(desc, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("PASS ", desc, "\n") }
  else { n_fail <<- n_fail + 1L; cat("FAIL ", desc, "\n") }
}
near <- function(a, b, tol = 1e-12) isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))

# ---- same simulated harmonised table as test_contracts.R ---------------------------------------
set.seed(20261029)
k <- 30
bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE)
sx <- runif(k, 0.005, 0.012); sy <- runif(k, 0.008, 0.02)
by <- 0.3 * bx + rnorm(k, 0.002, 0.01) + rnorm(k, 0, sy)
dat <- data.frame(variant_id = sprintf("1:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G",
                  beta_x = bx, se_x = sx, eaf_x = 0.3, beta_y = by, se_y = sy, eaf_y = 0.3,
                  p_x = 2 * pnorm(-abs(bx / sx)), p_y = 2 * pnorm(-abs(by / sy)), n_x = 1e5, n_y = 1e5,
                  harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
dat$mr_keep[c(4, 17)] <- FALSE
input <- list(kind = "pair_iv", data = dat, exposure_id = "src:opengwas:ieu-a-300",
              outcome_id = "src:opengwas:ieu-a-7", outcome_scale = "log_or")
SEED <- 123456L

# ---- parameter schema (C13.1) ---------------------------------------------------------------
expect("param schema complete for ivw_mre", length(validate_param_schema(ADAPTER_IVW_MRE_V11)) == 0)
expect("param schema complete for egger", length(validate_param_schema(ADAPTER_EGGER_V11)) == 0)
bad <- ADAPTER_IVW_MRE_V11; bad$params$ci_level$evidence <- NULL
expect("schema check catches conventional parameter without evidence", length(validate_param_schema(bad)) == 1)
expect("v1.0.0 adapters still build without emit_script (additive)", is.null(ADAPTER_IVW_MRE$emit_script))

# ---- emission is deterministic and clean (C6.2 rules 1-3) ----------------------------------
spec <- list(exposure_id = input$exposure_id, outcome_id = input$outcome_id, outcome_scale = input$outcome_scale,
             data_file = "analysis_data_ivw_mre.csv", estimates_file = "ivw_mre_estimates.csv",
             diagnostics_file = "ivw_mre_diagnostics.csv")
p_ivw <- resolve_params(ADAPTER_IVW_MRE_V11, list())
s1 <- ADAPTER_IVW_MRE_V11$emit_script(spec, p_ivw, SEED)$method_script
s2 <- ADAPTER_IVW_MRE_V11$emit_script(spec, p_ivw, SEED)$method_script
expect("emit_script gives identical bytes twice", identical(charToRaw(s1), charToRaw(s2)))
expect("script has no date, time, user or absolute path",
       !grepl("Sys\\.time|Sys\\.Date|Sys\\.info|/home/|/Users/|[A-Z]:[/\\\\]|tempdir|getwd", s1))
expect("script calls the package function with every argument written out",
       grepl("TwoSampleMR::mr_ivw_mre\\(", s1) && grepl("parameters = TwoSampleMR::default_parameters\\(\\)", s1) &&
       grepl("ci_level <- 0.95", s1) && grepl("set.seed\\(123456\\)", s1) && grepl("RNGkind\\(\"Mersenne-Twister\"", s1))
expect("script contains no app code (no source() of engine files, no mrgui calls)",
       !grepl("mrgui|adapter\\.R|runner\\.R|canonical\\.R", s1))
expect("hint text appears as a comment next to the setting", grepl("# Confidence level: Width of the confidence interval", s1))
expect("different settings give a different script", !identical(s1,
       ADAPTER_IVW_MRE_V11$emit_script(spec, resolve_params(ADAPTER_IVW_MRE_V11, list(ci_level = 0.9)), SEED)$method_script))

# ---- run the generated scripts in a fresh R and compare with the in-process oracle -----------
pair <- file.path(out_dir, "P01_pair")
methods <- list(list(adapter = ADAPTER_IVW_MRE_V11, input = input, params = list(ci_level = 0.9)),
                list(adapter = ADAPTER_EGGER_V11, input = input, params = list()))
info <- materialise_pair(pair, methods, SEED)
expect("pair folder has data, scripts and master script",
       file.exists(file.path(pair, "01_data", "analysis_data_ivw_mre.csv")) &&
       file.exists(file.path(pair, "02_scripts", "00_run_all.R")) && file.exists(file.path(pair, "02_scripts", "01_ivw_mre.R")))
run <- run_in_dir(pair, run_pair_script(pair))
expect("master script runs in a fresh R (--vanilla) with exit status 0", run$status == 0L)
if (run$status != 0L) cat(utils::tail(run$console, 20), sep = "\n")
res_ivw <- read_method_results(pair, "ivw_mre"); res_eg <- read_method_results(pair, "egger")

ctx0 <- list(seed = 1L)
o_ivw <- ADAPTER_IVW_MRE$run(input, resolve_params(ADAPTER_IVW_MRE, list(ci_level = 0.9)), ctx0)
o_eg <- ADAPTER_EGGER$run(input, resolve_params(ADAPTER_EGGER, list()), ctx0)
cmp <- function(a, b) { num <- names(a)[vapply(a, is.numeric, logical(1))]; all(vapply(num, function(n) near(a[[n]], b[[n]]), logical(1))) }
expect("IVW script result equals oracle (1e-12), ci_level 0.9", cmp(res_ivw$estimates, o_ivw$estimates) && cmp(res_ivw$diagnostics, o_ivw$diagnostics))
expect("Egger script result equals oracle (1e-12), slope and intercept", nrow(res_eg$estimates) == 2 &&
       cmp(res_eg$estimates, o_eg$estimates) && cmp(res_eg$diagnostics, o_eg$diagnostics))
expect("result tables keep text columns (estimand, scale, ci_method)",
       identical(res_ivw$estimates$estimand, "causal_effect") && identical(res_ivw$estimates$ci_method, "wald") &&
       identical(res_eg$estimates$ci_method, c("wald_t", "wald_t")) && identical(res_ivw$estimates$scale, "log_or"))
expect("console output stored", file.exists(file.path(pair, "03_console", "run_all_console.txt")))

# ---- replay check (C6.5) ----------------------------------------------------------------------
stored <- list(ivw_mre = res_ivw, egger = res_eg)
rp <- replay_check(pair, c("ivw_mre", "egger"), stored)
expect("replay on the same computer: identical", identical(rp$status, "identical"))
expect("replay compares all 27 numeric values (IVW 10, Egger 17) with zero difference", rp$n_values_compared == 27L && identical(rp$max_rel_diff, 0))

# user edits the script (changes ci_level) -> replay must notice and report 'differs'
sp <- file.path(pair, "02_scripts", "01_ivw_mre.R")
orig <- readLines(sp); hash_orig <- sha256_file(sp)
writeLines(sub("ci_level <- 0.9$", "ci_level <- 0.8", orig), sp)
expect("edited script has a different hash", sha256_file(sp) != hash_orig)
rp2 <- replay_check(pair, c("ivw_mre", "egger"), stored)
expect("replay reports 'differs' for an edited script", identical(rp2$status, "differs") && rp2$max_rel_diff > 1e-3)
# a script that stops with an error -> 'failed'
writeLines(c(orig, "stop(\"broken on purpose\")"), sp)
rp3 <- replay_check(pair, c("ivw_mre", "egger"), stored)
expect("replay reports 'failed' when the script errors", identical(rp3$status, "failed"))
writeLines(orig, sp)
expect("restored script replays identical again", identical(replay_check(pair, c("ivw_mre", "egger"), stored)$status, "identical"))

# ---- guard against the wrong working directory ----------------------------------------------
wrong <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
          c("--vanilla", shQuote(file.path(pair, "02_scripts", "00_run_all.R"))), stdout = TRUE, stderr = TRUE))
expect("wrong working directory stops with a clear message", any(grepl("Wrong working directory", wrong)) && !is.null(attr(wrong, "status")))

# ---- numbers survive the CSV round trip (17 significant digits) ----------------------------------
rt <- tempfile(fileext = ".csv"); write_analysis_csv(input$data, rt)
back <- utils::read.csv(rt, stringsAsFactors = FALSE)
expect("analysis CSV round trip within 2e-16 relative",
       max(abs(back$beta_x - input$data$beta_x) / abs(input$data$beta_x)) < 2.3e-16 &&
       max(abs(back$se_y - input$data$se_y) / abs(input$data$se_y)) < 2.3e-16)

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
