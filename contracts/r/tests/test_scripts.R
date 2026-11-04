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
near <- function(a, b, tol = 1e-12) { na <- is.na(a); length(a) == length(b) && identical(na, is.na(b)) && isTRUE(all(abs(a - b)[!na] <= (tol * pmax(1, abs(b)))[!na])) }   # missing values must be missing in both

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
       grepl("TwoSampleMR::mr_ivw_mre,", s1) && grepl("fit <- fit_fun\\(", s1) && grepl("ivw_model <- \"mre\"", s1) && grepl("min_f <- 0", s1) && grepl("parameters = TwoSampleMR::default_parameters\\(\\)", s1) &&
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
# (files are written as raw bytes: writeLines() would turn line ends into CRLF on Windows)
put <- function(lines, path) writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), path)
sp <- file.path(pair, "02_scripts", "01_ivw_mre.R")
orig_raw <- readBin(sp, "raw", file.size(sp)); orig <- readLines(sp); hash_orig <- sha256_file(sp)
put(sub("ci_level <- 0.9$", "ci_level <- 0.8", orig), sp)
expect("edited script has a different hash", sha256_file(sp) != hash_orig)
rp2 <- replay_check(pair, c("ivw_mre", "egger"), stored)
expect("replay reports 'differs' for an edited script", identical(rp2$status, "differs") && rp2$max_rel_diff > 1e-3)
# a script that stops with an error -> 'failed'
put(c(orig, "stop(\"broken on purpose\")"), sp)
rp3 <- replay_check(pair, c("ivw_mre", "egger"), stored)
expect("replay reports 'failed' when the script errors", identical(rp3$status, "failed"))
writeBin(orig_raw, sp)
expect("restored script is byte-identical to the generated one", identical(sha256_file(sp), hash_orig))
expect("restored script replays identical again", identical(replay_check(pair, c("ivw_mre", "egger"), stored)$status, "identical"))

# ---- guard against the wrong working directory ----------------------------------------------
Sys.setenv(RHOME = R.home())   # the test calls Rscript directly
wrong <- suppressWarnings(system2(default_rscript(),
          c("--vanilla", shQuote(file.path(pair, "02_scripts", "00_run_all.R"))), stdout = TRUE, stderr = TRUE))
expect("wrong working directory stops with a clear message", any(grepl("Wrong working directory", wrong)) && !is.null(attr(wrong, "status")))

# ---- numbers survive the CSV round trip (17 significant digits) ----------------------------------
rt <- tempfile(fileext = ".csv"); write_analysis_csv(input$data, rt)
back <- utils::read.csv(rt, stringsAsFactors = FALSE)
expect("analysis CSV round trip within 2e-16 relative",
       max(abs(back$beta_x - input$data$beta_x) / abs(input$data$beta_x)) < 2.3e-16 &&
       max(abs(back$se_y - input$data$se_y) / abs(input$data$se_y)) < 2.3e-16)


# ---- Round 8: bootstrap-based robust methods: script equals oracle with the same seed, and replay is identical ----
pair8 <- file.path(out_dir, "P08_robust")
prm8 <- function(a, l) resolve_params(a, l)
m8 <- list(list(adapter = ADAPTER_WM_V11, input = input, params = list(nboot = 200L)),
           list(adapter = ADAPTER_SMED_V11, input = input, params = list(nboot = 200L)),
           list(adapter = ADAPTER_MODE_V11, input = input, params = list(nboot = 200L, phi = 0.5)),
           list(adapter = ADAPTER_SMODE_V11, input = input, params = list(nboot = 200L, phi = 0.5)))
materialise_pair(pair8, m8, SEED)
run8 <- run_in_dir(pair8, run_pair_script(pair8))
expect("robust methods: master script runs in a fresh R", run8$status == 0L)
stored8 <- list()
for (a in list(ADAPTER_WM_V11, ADAPTER_SMED_V11, ADAPTER_MODE_V11, ADAPTER_SMODE_V11)) {
  pr <- prm8(a, if (a$id %in% c("mode", "smode")) list(nboot = 200L, phi = 0.5) else list(nboot = 200L))
  o <- a$run(input, pr, list(seed = SEED)); r <- read_method_results(pair8, a$id); stored8[[a$id]] <- r
  expect(sprintf("%s: script result equals oracle (1e-12) with the same seed", a$id), cmp(r$estimates, o$estimates) && cmp(r$diagnostics, o$diagnostics))
}
rp8 <- replay_check(pair8, c("wm", "smed", "mode", "smode"), stored8)
expect("robust methods: replay is identical", identical(rp8$status, "identical"))

# ---- Round 10: outlier detectors. Planted outliers must be flagged; script equals oracle; replay identical ----
set.seed(77); k10 <- 14
bx10 <- runif(k10, 0.04, 0.12); sx10 <- runif(k10, 0.005, 0.01); sy10 <- runif(k10, 0.008, 0.015)
by10 <- 0.3 * bx10 + rnorm(k10, 0, sy10); by10[c(3, 9)] <- by10[c(3, 9)] + 0.06
dat10 <- data.frame(variant_id = sprintf("rs%d", seq_len(k10)), ea = "A", oa = "G", beta_x = bx10, se_x = sx10, eaf_x = 0.3,
                    beta_y = by10, se_y = sy10, eaf_y = 0.3, p_x = 1e-9, p_y = 0.5, n_x = 1e5, n_y = 1e5,
                    harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
input10 <- list(kind = "pair_iv", data = dat10, exposure_id = input$exposure_id, outcome_id = input$outcome_id, outcome_scale = "beta")
pair10 <- file.path(out_dir, "P10_detectors")
m10 <- list(list(adapter = ADAPTER_PRESSO_V11, input = input10, params = list(nbdist = 500L)),
            list(adapter = ADAPTER_RADIAL_V11, input = input10, params = list(alpha = 0.05, weights = "3")),
            list(adapter = ADAPTER_RADEGGER_V11, input = input10, params = list(alpha = 0.05)))
materialise_pair(pair10, m10, SEED)
run10 <- run_in_dir(pair10, run_pair_script(pair10))
expect("detectors: master script runs in a fresh R", run10$status == 0L)
stored10 <- list()
for (m in m10) {
  a <- m$adapter; pr <- resolve_params(a, m$params)
  o <- a$run(input10, pr, list(seed = SEED)); r <- read_method_results(pair10, a$id); stored10[[a$id]] <- r
  expect(sprintf("%s: script result equals oracle (1e-12)", a$id), cmp(r$estimates, o$estimates) && cmp(r$diagnostics, o$diagnostics))
  nm <- if (a$id == "presso") "presso_outliers" else "radial_outliers"
  flagged <- r$diagnostics$note[r$diagnostics$name == nm]
  expect(sprintf("%s: flags the planted outlier rs3 and rs9", a$id), !is.na(flagged) && all(c("rs3", "rs9") %in% strsplit(flagged, ";")[[1]]))
}
rp10 <- replay_check(pair10, c("presso", "radial", "radegger"), stored10)
expect("detectors: replay is identical", identical(rp10$status, "identical"))
expect("detector schemas complete", all(vapply(list(ADAPTER_PRESSO_V11, ADAPTER_RADIAL_V11, ADAPTER_RADEGGER_V11), function(a) length(validate_param_schema(a)) == 0, logical(1))))

# ---- Round 15: Steiger directionality (diagnostic-only) ----
pair15 <- file.path(out_dir, "P15_steiger")
m15 <- list(list(adapter = ADAPTER_STEIGER_V11, input = input10, params = list(r_xxo = 0.9, r_yyo = 0.95)))
materialise_pair(pair15, m15, SEED)
run15 <- run_in_dir(pair15, run_pair_script(pair15))
expect("steiger: master script runs in a fresh R", run15$status == 0L)
r15 <- read_method_results(pair15, "steiger"); o15 <- ADAPTER_STEIGER_V11$run(input10, resolve_params(ADAPTER_STEIGER_V11, m15[[1]]$params), list(seed = SEED))
expect("steiger: script result equals oracle (1e-12)", cmp(r15$diagnostics, o15$diagnostics) && nrow(r15$estimates) == 0L && nrow(o15$estimates) == 0L)
expect("steiger: replay is identical (header-only estimates file)", identical(replay_check(pair15, "steiger", list(steiger = r15))$status, "identical"))
expect("steiger: schema complete and flagged diagnostic-only", length(validate_param_schema(ADAPTER_STEIGER_V11)) == 0 && isTRUE(ADAPTER_STEIGER_V11$diagnostic_only) && !isTRUE(ADAPTER_IVW_MRE_V11$diagnostic_only))
expect("steiger: blocks a binary outcome without cases, controls and prevalence, and missing sample sizes", identical(ADAPTER_STEIGER_V11$check(modifyList(input10, list(outcome_scale = "log_or")), list())$status, "block") &&
       identical(ADAPTER_STEIGER_V11$check(modifyList(input10, list(data = transform(dat10, n_y = NA_real_))), list())$status, "block"))

# ---- Round 16: Steiger with a binary outcome and a binary exposure (liability scale, TwoSampleMR::get_r_from_lor; Lee 2012) ----
datb <- transform(dat10, eaf_y = 0.3, ncase_y = 20000, ncontrol_y = 80000, prevalence_y = 0.05)
inputb <- modifyList(input10, list(data = datb, outcome_scale = "log_or"))
expect("steiger binary: check passes when cases, controls, prevalence and frequency are given", identical(ADAPTER_STEIGER_V11$check(inputb, list())$status, "ok") || identical(ADAPTER_STEIGER_V11$check(inputb, list())$status, "warn"))
pair16 <- file.path(out_dir, "P16_steiger_binary")
m16 <- list(list(adapter = ADAPTER_STEIGER_V11, input = inputb, params = list(r_xxo = 1, r_yyo = 1)))
materialise_pair(pair16, m16, SEED)
run16 <- run_in_dir(pair16, run_pair_script(pair16))
expect("steiger binary: master script runs in a fresh R", run16$status == 0L)
r16 <- read_method_results(pair16, "steiger"); o16 <- ADAPTER_STEIGER_V11$run(inputb, resolve_params(ADAPTER_STEIGER_V11, m16[[1]]$params), list(seed = SEED))
expect("steiger binary: script result equals oracle (1e-12)", cmp(r16$diagnostics, o16$diagnostics))
expect("steiger binary: replay is identical", identical(replay_check(pair16, "steiger", list(steiger = r16))$status, "identical"))
rx16 <- TwoSampleMR::get_r_from_bsen(datb$beta_x, datb$se_x, datb$n_x)
ry16 <- TwoSampleMR::get_r_from_lor(datb$beta_y, datb$eaf_y, 20000, 80000, 0.05, model = "logit", correction = FALSE)
st16 <- TwoSampleMR::mr_steiger(rep(NA_real_, nrow(datb)), rep(NA_real_, nrow(datb)), datb$n_x, rep(100000, nrow(datb)), rx16, ry16, 1, 1)
g16 <- function(n) r16$diagnostics$value[r16$diagnostics$name == n]
expect("steiger binary: equals get_r_from_lor and mr_steiger called directly", abs(g16("steiger_r2_outcome") - st16$r2_out) < 1e-12 && abs(g16("steiger_r2_exposure") - st16$r2_exp) < 1e-12 &&
       abs(r16$diagnostics$p[r16$diagnostics$name == "steiger_direction"] - st16$steiger_test) < 1e-12 && g16("steiger_outcome_binary") == 1 && g16("steiger_exposure_binary") == 0)
expect("steiger binary: the script says liability scale and names the prevalence", grepl("get_r_from_lor", paste(readLines(file.path(pair16, "02_scripts", "01_steiger.R")), collapse = "\n"), fixed = TRUE) &&
       grepl("prevalence 0.05", r16$diagnostics$note[r16$diagnostics$name == "steiger_outcome_binary"], fixed = TRUE))
expect("steiger binary: blocked without prevalence, without frequency, without case counts",
       identical(ADAPTER_STEIGER_V11$check(modifyList(inputb, list(data = transform(datb, prevalence_y = NA_real_))), list())$codes, "FEAS_STEIGER_NEEDS_PREVALENCE") &&
       identical(ADAPTER_STEIGER_V11$check(modifyList(inputb, list(data = transform(datb, eaf_y = NA_real_))), list())$codes, "FEAS_STEIGER_NEEDS_EAF") &&
       identical(ADAPTER_STEIGER_V11$check(modifyList(inputb, list(data = transform(datb, ncase_y = NA_real_))), list())$codes, "FEAS_STEIGER_NEEDS_CASES"))
datbx <- transform(dat10, eaf_x = 0.3, ncase_x = 15000, ncontrol_x = 150000, prevalence_x = 0.1)
inputbx <- modifyList(input10, list(data = datbx, exposure_scale = "log_or"))
obx <- ADAPTER_STEIGER_V11$run(inputbx, resolve_params(ADAPTER_STEIGER_V11, list()), list(seed = SEED))
stbx <- TwoSampleMR::mr_steiger(rep(NA_real_, nrow(datbx)), rep(NA_real_, nrow(datbx)), rep(165000, nrow(datbx)), datbx$n_y,
         TwoSampleMR::get_r_from_lor(datbx$beta_x, datbx$eaf_x, 15000, 150000, 0.1, model = "logit", correction = FALSE), TwoSampleMR::get_r_from_bsen(datbx$beta_y, datbx$se_y, datbx$n_y), 1, 1)
expect("steiger binary exposure: equals mr_steiger called directly, and a binary exposure without counts is blocked",
       abs(obx$diagnostics$value[obx$diagnostics$name == "steiger_r2_exposure"] - stbx$r2_exp) < 1e-12 &&
       identical(ADAPTER_STEIGER_V11$check(modifyList(input10, list(exposure_scale = "log_or")), list())$codes, "FEAS_STEIGER_NEEDS_CASES"))
expect("steiger continuous: the two scale rows say continuous", all(r15$diagnostics$value[r15$diagnostics$name %in% c("steiger_exposure_binary", "steiger_outcome_binary")] == 0) &&
       sum(r15$diagnostics$name %in% c("steiger_exposure_binary", "steiger_outcome_binary")) == 2L)

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
