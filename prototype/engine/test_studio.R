# Headless tests for the screen logic and the Shiny server (no browser, no desktop shell).
# Run:  Rscript prototype/engine/test_studio.R        (from the repository root)
args_dir <- normalizePath(file.path("prototype", "engine"), winslash = "/")
source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir)); source(file.path(args_dir, "import_core.R")); source(file.path(args_dir, "liftover_core.R")); source(file.path(args_dir, "opengwas_core.R")); source(file.path(args_dir, "ld_core.R")); source(file.path(args_dir, "strobe_core.R")); source(file.path(args_dir, "feasibility_core.R")); source(file.path(args_dir, "report_core.R")); source(file.path(args_dir, "pdf_core.R")); source(file.path(args_dir, "outlier_core.R")); source(file.path(args_dir, "mv_core.R")); source(file.path(args_dir, "mediation_core.R")); source(file.path(args_dir, "bidir_core.R")); source(file.path(args_dir, "resource_core.R")); source(file.path(args_dir, "ldsc_core.R")); source(file.path(args_dir, "lhcmr_runner.R")); source(file.path(args_dir, "lhcmr_core.R")); source(file.path(args_dir, "mrlap_core.R")); source(file.path(args_dir, "mrapss_core.R")); source(file.path(args_dir, "cause_core.R")); source(file.path(args_dir, "fullgwas_core.R")); source(file.path(args_dir, "live_core.R")); source(file.path(args_dir, "bg_core.R")); source(file.path(args_dir, "app.R"))
root <- file.path(tempfile("studies-"), "MR Studio"); dir.create(root, recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root)

n_pass <- 0L; n_fail <- 0L
expect <- function(desc, ok) { if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("PASS ", desc, "\n") } else { n_fail <<- n_fail + 1L; cat("FAIL ", desc, "\n") } }
err_of <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

expect("studies root honours MRGUI_STUDIES_DIR", identical(studies_root(), root))
expect("study names are cleaned for the file system", identical(clean_study_name("  Omega-3: TG/risk? "), "Omega-3 TG risk") && identical(clean_study_name("..."), ""))
expect("empty study name is refused", grepl("Give the study a name", err_of(studio_new("   "))))
expect("a study in a cloud folder is refused", grepl("SYS_STUDY_FOLDER_REFUSED", err_of(studio_new("S", file.path(tempdir(), "OneDrive")))))
d <- demo_pair(args_dir)
expect("demo data: 40 simulated instruments with the contract columns",
       nrow(d$input$data) == 40L && all(c("beta_x", "se_x", "beta_y", "se_y", "mr_keep") %in% names(d$input$data)) && is.logical(d$input$data$mr_keep))

h <- studio_new("Test study")
expect("new study has the folder tree and appears in the list", all(dir.exists(file.path(h$dir, STUDY_FOLDERS))) && "Test study" %in% studies_list())
expect("creating the same name twice is refused", grepl("already exists", err_of(studio_new("Test study"))))
r <- studio_run(h, "Demo exposure", "Demo outcome", d$input)
expect("first run finished with both methods ok and an identical replay",
       !r$reused && all(vapply(r$results, function(x) x$status, "") == "ok") && identical(r$replay$status, "identical"))
res <- pair_results(h$dir, r$pair_folder)
expect("results table has IVW and Egger rows with estimates, intervals and P values",
       !is.null(res) && all(c("Method", "estimate", "ci_lower", "ci_upper", "p") %in% names(res)) && any(grepl("IVW", res$Method)) && any(grepl("Egger", res$Method)))
fr <- format_results(res)
expect("formatted table: readable headers, interval with level, P never printed as 0",
       identical(names(fr), c("Method", "Estimate", "Std. error", "Interval", "P value", "Variants", "Estimand")) &&
       all(grepl("^-?[0-9.]+ to -?[0-9.]+ \\(95%\\)$", fr$Interval[!is.na(res$ci_lower)])) && !any(fr$`P value` %in% c("0.000", "0")) &&
       "Pleiotropy intercept" %in% fr$Estimand)
expect("size formatting", identical(fmt_bytes(c(500, 2048, 3 * 1024^2)), c("500 B", "2.0 KB", "3.0 MB")))
expect("the IVW estimate recovers the simulated effect (0.25) within 0.1", abs(res$estimate[grepl("IVW", res$Method)][1] - 0.25) < 0.1)
expect("same labels give the same pair folder and the second run is a reuse",
       identical(studio_run(h, "Demo exposure", "Demo outcome", d$input)$pair_folder, r$pair_folder) && isTRUE(studio_run(h, "Demo exposure", "Demo outcome", d$input)$reused))
r2 <- studio_run(h, "Demo exposure", "Demo outcome", d$input, ci_level = 0.9)
expect("a changed confidence level re-runs (not a reuse) in the same pair folder", !isTRUE(r2$reused) && identical(r2$pair_folder, r$pair_folder))
r3 <- studio_run(h, "Other exposure", "Demo outcome", d$input, methods = "ivw_mre")
expect("a new pair gets the next index and one method can be run alone", grepl("^P02_", r3$pair_folder) && all(vapply(r3$results, function(x) x$status, "") == "ok") && length(r3$results) == 1L)
expect("only the chosen method was run for the new pair", identical(unique(pair_results(h$dir, r3$pair_folder)$Method), "IVW (multiplicative random effects)"))
expect("no method chosen is refused", grepl("at least one method", err_of(studio_run(h, "A", "B", d$input, character()))))
expect("file list comes from the manifest and is not empty", nrow(study_files(h$dir)) >= 11L)
expect("opening only works inside the studies root", grepl("Only study folders", err_of(open_in_file_manager(tempdir(), root))))
DBI::dbDisconnect(h$con)
h2 <- studio_open("Test study")
expect("reopening shows the same pairs and results", identical(studio_pair_folders(h2$dir), c(r$pair_folder, r3$pair_folder)) && !is.null(pair_results(h2$dir, r$pair_folder)))
DBI::dbDisconnect(h2$con)

# ---- new core functions (settings, catalogue, facts, overview, payload, diagnostics) --------------------------------------
Sys.setenv(MRGUI_CONFIG_DIR = tempfile("cfg-"))
st <- settings_read()
expect("settings default: paper theme (light), setup not done", identical(st$theme, "paper") && identical(st$setup_done, FALSE))
settings_set("theme", "rose"); expect("a theme choice is saved and read back", identical(settings_read()$theme, "rose"))
settings_set("theme", "peach"); expect("the Sage and Peach themes are accepted", identical(settings_read()$theme, "peach") && "sage" %in% SETTINGS_THEMES)
expect("an unknown theme is refused", grepl("theme", err_of(settings_set("theme", "neon")), ignore.case = TRUE))
expect("a wrong value type is refused", !is.na(err_of(settings_set("offline", "yes"))))
expect("a cloud-synchronised studies folder is refused", !is.na(err_of(settings_set("studies_root", file.path(tempdir(), "Dropbox", "x")))))
cat_ <- demo_catalogue()
expect("catalogue: 4 exposures and 4 outcomes, all simulated and EUR", length(cat_$exposures) == 4L && length(cat_$outcomes) == 4L && all(vapply(c(cat_$exposures, cat_$outcomes), function(x) x$anc == "EUR" && grepl("simulated", x$src), TRUE)))
f <- pair_facts("expA", "outA")
expect("pair facts: 40 instruments and a positive mean F", f$instruments == 40L && f$mean_f > 10 && f$scale == "beta")
expect("the same pair always gives the same data, a different outcome changes only the outcome", identical(demo_input("expA", "outA")$input$data, demo_input("expA", "outA")$input$data) &&
         identical(demo_input("expA", "outA")$input$data$beta_x, demo_input("expA", "outB")$input$data$beta_x) && !identical(demo_input("expA", "outA")$input$data$beta_y, demo_input("expA", "outB")$input$data$beta_y))
expect("an unknown trait is refused", grepl("Unknown", err_of(demo_input("zzz", "outA"))))
hh <- studio_new("Overview study"); study_plan_save(hh, list(kind = "basic", aim = "x", prereg = FALSE))
expect("the plan is saved as a file and in the database", file.exists(file.path(hh$dir, "01_study_plan", "study_plan.json")) && identical(ui_state_get(hh$con, "plan")$aim, "x"))
ov <- studies_overview(); expect("overview: a study without results is a DRAFT", ov[[which(vapply(ov, function(x) x$name, "") == "Overview study")]]$status == "DRAFT")
# two different traits whose labels slug to the same folder name must not share a pair folder
ra <- studio_run(hh, "Simulated exposure A", "Simulated outcome A", demo_input("expA", "outA")$input, idx = 1L)
rb <- studio_run(hh, "Simulated exposure B", "Simulated outcome A", demo_input("expB", "outA")$input, idx = 2L)
expect("two traits with the same 30-character slug get separate pair folders", !identical(ra$pair_folder, rb$pair_folder) && grepl("^P01_", ra$pair_folder) && grepl("^P02_", rb$pair_folder))
# method settings: each choice must give exactly what the TwoSampleMR function gives, and the label must say which was run
di <- demo_input("expA", "outA"); dd <- di$input$data; dd <- dd[dd$mr_keep, ]
ref <- function(f) suppressWarnings(f(dd$beta_x, dd$beta_y, dd$se_x, dd$se_y, TwoSampleMR::default_parameters()))
for (md in c("mre", "mre_floor", "fixed")) {
  rr <- studio_run(hh, paste("Model", md), "Simulated outcome A", di$input, "ivw_mre", params = list(ivw_model = md, ci_level = 0.9), idx = 10L + match(md, c("mre", "mre_floor", "fixed")))
  fn <- c(mre = TwoSampleMR::mr_ivw_mre, mre_floor = TwoSampleMR::mr_ivw, fixed = TwoSampleMR::mr_ivw_fe)[[md]]
  pe <- pair_payload(hh$dir, rr$pair_folder)$estimates[[1]]; fr <- ref(fn)
  expect(paste("IVW model", md, "matches the TwoSampleMR function and the 90% level"), abs(pe$se - fr$se) < 1e-9 && abs(pe$estimate - fr$b) < 1e-9 && abs(pe$ci_level - 0.9) < 1e-12 && identical(pe$method_name, IVW_MODEL_NAMES[[md]]))
}
# Round 30: the default IVW model is TwoSampleMR::mr_ivw (multiplicative random effects with the SE not below the fixed-effect SE), as in the Burgess guidelines
pdef <- pair_payload(hh$dir, ra$pair_folder); edef <- Filter(function(e) identical(e$method, "ivw_mre"), pdef$estimates)[[1]]; fdef <- ref(TwoSampleMR::mr_ivw)
expect("Round 30: with no setting the IVW model is TwoSampleMR::mr_ivw, the TwoSampleMR default, and is labelled as such",
       identical(pdef$settings$ivw_model, "mre_floor") && abs(edef$se - fdef$se) < 1e-9 && abs(edef$p - fdef$pval) < 1e-12 && identical(edef$method_name, "IVW (multiplicative random effects)") &&
       identical(resolve_params(ADAPTER_IVW_MRE_V11, list())$ivw_model, "mre_floor") && identical(ADAPTER_IVW_MRE_V11$params$ivw_model$values[1], "mre_floor"))
sigl <- summary(stats::lm(beta_y ~ -1 + beta_x, data = dd, weights = 1 / dd$se_y^2))$sigma
expect("Round 30: the default SE is never below the fixed-effect SE", edef$se >= ref(TwoSampleMR::mr_ivw_fe)$se - 1e-12 && (sigl >= 1 || abs(edef$se - ref(TwoSampleMR::mr_ivw_fe)$se) < 1e-12))
pst <- pair_payload(hh$dir, rr$pair_folder)$settings
expect("the payload reports the settings the pair was run with, read from its scripts", identical(pst$ivw_model, "fixed") && abs(pst$ci_level - 0.9) < 1e-12 && pst$min_f == 0)
big <- studio_run(hh, "F filter", "Simulated outcome A", di$input, c("ivw_mre", "egger"), params = list(min_f = 30), idx = 20L)
fk <- sum((dd$beta_x / dd$se_x)^2 >= 30); ps <- pair_payload(hh$dir, big$pair_folder)
expect("minimum F leaves out weak variants in both methods", fk < nrow(dd) && all(vapply(ps$estimates, function(e) e$n_ivs_used == fk, TRUE)))
weak <- di$input; weak$data$se_x <- weak$data$se_x * 100   # every variant now has F below 10
tooHigh <- tryCatch(studio_run(hh, "F too high", "Simulated outcome A", weak, "ivw_mre", params = list(min_f = 10), idx = 21L), error = function(e) e)
expect("a minimum F that leaves too few variants stops with a clear message", inherits(tooHigh, "error") || any(vapply(tooHigh$results, function(x) !identical(x$status, "ok"), TRUE)))
expect("an unknown model value is refused", inherits(tryCatch(studio_run(hh, "Bad", "Simulated outcome A", di$input, "ivw_mre", params = list(ivw_model = "nonsense"), idx = 22L), error = function(e) e), "error"))
pl <- pair_payload(hh$dir, ra$pair_folder)
expect("payload: estimates and diagnostics are arrays, variants and replay present", is.null(names(pl$estimates)) && is.null(names(pl$diagnostics)) && length(pl$variants) == 40L && identical(pl$replay$status, "identical"))
ov <- studies_overview(); expect("overview: a study with results is COMPLETE", ov[[which(vapply(ov, function(x) x$name, "") == "Overview study")]]$status == "COMPLETE")
dg <- run_diagnostics(root); expect("diagnostics: every item has a name, ok and detail", all(vapply(dg, function(x) all(c("name", "ok", "detail") %in% names(x)), TRUE)))
expect("diagnostics: the studies folder is writable", dg[[4]]$ok)
expect("folder size is counted", folder_size_bytes(hh$dir) > 1000)
DBI::dbDisconnect(hh$con)

# ---- instrument review, sensitivity branch, RStudio lookup, compare with default --------------------------------------------
it <- instrument_table("expA", "outA")
expect("instrument table: 40 rows, weight shares add up to 1, ratio is beta_y / beta_x", nrow(it) == 40L && abs(sum(it$share) - 1) < 1e-9 && all.equal(it$ratio, it$by / it$bx))
h3 <- studio_new("Sensitivity study"); base <- studio_run(h3, "Simulated exposure A", "Simulated outcome A", demo_input("expA", "outA")$input, idx = 1L)
top <- it$id[which.max(it$share)]
sens <- studio_run_without(h3, "expA", "outA", top, idx = 2L)
pb <- pair_payload(h3$dir, base$pair_folder); ps <- pair_payload(h3$dir, sens$pair_folder)
expect("sensitivity branch: own folder, one variant fewer, primary results untouched", grepl("^P02_", sens$pair_folder) && sens$kept == 39L &&
         pb$estimates[[1]]$n_ivs_used == 40L && ps$estimates[[1]]$n_ivs_used == 39L && identical(ps$replay$status, "identical"))
expect("sensitivity branch with an unknown variant is refused", grepl("not instruments", err_of(studio_run_without(h3, "expA", "outA", "rs0", idx = 3L))))
expect("sensitivity branch must leave at least 2 instruments", grepl("At least 2", err_of(studio_run_without(h3, "expA", "outA", it$id[1:39], idx = 3L))))
cd <- ci_default_rows(pb$estimates)
expect("compare with default: the 95% interval recomputed from the saved result equals the saved 95% interval", abs(cd[[1]]$lo - pb$estimates[[1]]$ci_lower) < 1e-9 && abs(cd[[1]]$hi - pb$estimates[[1]]$ci_upper) < 1e-9)
expect("compare with default: Egger uses the t quantile", abs(cd[[2]]$lo - pb$estimates[[2]]$ci_lower) < 1e-9)
h4 <- studio_new("Egger level"); r90 <- studio_run(h4, "A", "B", demo_input("expA", "outA")$input, ci_level = 0.9, idx = 1L); p90 <- pair_payload(h4$dir, r90$pair_folder)
expect("a changed confidence level reaches both methods", all(vapply(Filter(function(e) e$estimand == "causal_effect", p90$estimates), function(e) abs(e$ci_level - 0.9) < 1e-12, TRUE)))
DBI::dbDisconnect(h4$con)
expect("RStudio lookup returns a path or NA", is.character(rstudio_path()) && length(rstudio_path()) == 1L)
expect("open in RStudio refuses a pair without scripts", grepl("no script", err_of(open_in_rstudio(h3$dir, "P99_none"))))
pd_v <- file.path(h3$dir, "03_pairs", base$pair_folder)
expect("check files: a new pair folder holds VERIFY_IN_R.R, a README and an .Rproj", all(file.exists(file.path(pd_v, c("VERIFY_IN_R.R", "README_RERUN.md", "MR_Studio_pair.Rproj")))))
expect("check script: base R only, no setwd and no reference to the studies folder", {tx <- readLines(file.path(pd_v, "VERIFY_IN_R.R")); !any(grepl("^\\s*setwd\\(", tx)) && !any(grepl("mrgui|studio_core|contracts", tx))})
tab_md5 <- function(d) unname(tools::md5sum(list.files(file.path(d, "04_tables"), full.names = TRUE)))
md5_before <- tab_md5(pd_v)
run_verify <- function(dir, from = tempdir()) {
  old <- setwd(from); on.exit(setwd(old)); lib_old <- Sys.getenv("R_LIBS", unset = NA); Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)); rh_old <- Sys.getenv("RHOME", unset = NA); Sys.setenv(RHOME = R.home())
  on.exit(if (is.na(lib_old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = lib_old), add = TRUE)
  on.exit(if (is.na(rh_old)) Sys.unsetenv("RHOME") else Sys.setenv(RHOME = rh_old), add = TRUE)
  suppressWarnings(system2(default_rscript(), c("--vanilla", shQuote(file.path(dir, "VERIFY_IN_R.R"))), stdout = TRUE, stderr = TRUE))
}
vo <- run_verify(pd_v)
expect("check script: started from another folder, it finds the pair folder and reports that every number matches", any(grepl("ALL [0-9]+ NUMBERS MATCH", vo)))
expect("check script: the result is saved in a new run folder and the original results are untouched", {rf <- list.files(file.path(pd_v, "09_verification"), pattern = "^run_", full.names = TRUE); length(rf) == 1L && file.exists(file.path(rf, "verification_result.txt")) && file.exists(file.path(rf, "04_tables", "ivw_mre_estimates.csv")) && identical(tab_md5(pd_v), md5_before)})
src_old <- setwd(tempdir()); src_res <- tryCatch(suppressMessages(utils::capture.output(src_v <- source(file.path(pd_v, "VERIFY_IN_R.R"), echo = FALSE))), error = function(e) conditionMessage(e)); src_wd <- getwd(); setwd(src_old)
expect("check script: the Source button route (source()) works from another folder, matches, and leaves the working folder as it was", is.list(src_v$value) && isTRUE(src_v$value$match) && identical(normalizePath(src_wd), normalizePath(tempdir())))
expect("check script: the headline table shows the app estimate next to this run", any(grepl("MR_Studio", vo)) && any(grepl("This_run", vo)))
# a changed saved result must be reported
pd_alt <- file.path(tempfile("alt-"), "pair"); dir.create(dirname(pd_alt)); file.copy(pd_v, dirname(pd_alt), recursive = TRUE); pd_alt <- file.path(dirname(pd_alt), basename(pd_v)); unlink(file.path(pd_alt, "09_verification"), recursive = TRUE)
ef <- file.path(pd_alt, "04_tables", "ivw_mre_estimates.csv"); et <- utils::read.csv(ef, stringsAsFactors = FALSE); et$estimate[1] <- et$estimate[1] * 1.001; utils::write.csv(et, ef, row.names = FALSE)
vo2 <- run_verify(pd_alt)
expect("check script: a changed saved number is reported as a difference and names the file", any(grepl("DIFFERENCES FOUND", vo2)) && any(grepl("ivw_mre_estimates.csv", vo2)))
# an older pair folder (made before the check files existed) is upgraded when it is opened
unlink(file.path(pd_v, c("VERIFY_IN_R.R", "README_RERUN.md", "MR_Studio_pair.Rproj")))
vp <- verify_script_path(h3$dir, base$pair_folder)
expect("check files: an older pair folder gets them when opened, and the path stays inside the pair folder", file.exists(vp$script) && startsWith(vp$script, vp$pair) && identical(basename(vp$script), "VERIFY_IN_R.R"))
expect("check script: a pair folder outside the studies folder is refused", grepl("Only scripts", err_of(verify_script_path(h3$dir, base$pair_folder, root = file.path(tempdir(), "other_root")))))
# a missing package is explained with the exact install lines
pd_np <- file.path(tempfile("np-"), basename(pd_v)); dir.create(dirname(pd_np)); file.copy(pd_v, dirname(pd_np), recursive = TRUE)
sf <- list.files(file.path(pd_np, "02_scripts"), pattern = "^01_", full.names = TRUE); sx <- readLines(sf); sx <- sub("^# Package  : TwoSampleMR", "# Package  : NoSuchPackageXyz", sx); writeLines(sx, sf)
vo3 <- run_verify(pd_np)
expect("check script: a missing package gets the install lines, not an obscure error", any(grepl("remotes::install_github", vo3)) && any(grepl("NoSuchPackageXyz", vo3)))
expect("customise panel schema comes from the adapters", all(c("label", "hint", "min", "max", "default", "type") %in% names(param_schemas()$ivw_mre$ci_level)) && identical(names(param_schemas()$ivw_mre), c("ci_level", "ivw_model", "min_f")) && identical(names(param_schemas()$egger), c("ci_level", "min_f")) && identical(unlist(param_schemas()$ivw_mre$ivw_model$values), c("mre_floor", "mre", "fixed")))
DBI::dbDisconnect(h3$con)

# ---- data import ----------------------------------------------------------------------------------------------------------------
set.seed(1); ni <- 3000
g <- data.frame(MarkerName = paste0("rs", 100000 + seq_len(ni)), Chromosome = 1, BP_pos = 1e6 + seq_len(ni) * 100, Tested_Allele = sample(c("A", "C"), ni, TRUE), Other_Allele = sample(c("G", "T"), ni, TRUE),
                Frq = runif(ni, .05, .95), Effect = rnorm(ni, 0, .05), StdErr = runif(ni, .01, .03), `P-value` = runif(ni), Neff = 1e5, check.names = FALSE)
fi <- tempfile(fileext = ".tsv"); utils::write.table(g, fi, sep = "\t", quote = FALSE, row.names = FALSE)
hd <- import_read_head(fi, "bolt.tsv"); sg <- import_suggest(hd$columns, hd$rows); mp <- lapply(sg, `[[`, "column")
expect("import: headings in another spelling and order are all recognised", all(c("snp", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p", "n") %in% names(mp)) && identical(mp$p, "P-value") && identical(mp$beta, "Effect"))
g2 <- g; names(g2)[names(g2) == "StdErr"] <- "Stderrr"; fz <- tempfile(fileext = ".csv"); utils::write.csv(g2[, c(5:1, 6:10)], fz, row.names = FALSE)
h2i <- import_read_head(fz); s2 <- import_suggest(h2i$columns, h2i$rows)
expect("import: a misspelled heading gets a medium-confidence suggestion, and column order does not matter", identical(s2$se$column, "Stderrr") && identical(s2$se$confidence, "medium") && identical(s2$snp$column, "MarkerName"))
expect("import: the separator is detected (tab and comma)", identical(hd$sep, "\t") && identical(h2i$sep, ","))
fl <- list(frq = "eaf", scope = "selected", scale = "beta", n = "")
v <- import_validate(fi, hd$sep, hd$columns, mp, fl)
expect("import: a clean file passes with its counts", v$ok && v$stats$rows == ni && v$stats$kept == ni)
vb1 <- import_validate(fi, hd$sep, hd$columns, mp, c(fl, list(kind = "binary", ncase = "1200", ncontrol = "")))
vb2 <- import_validate(fi, hd$sep, hd$columns, mp, c(fl, list(kind = "binary")))
vb3 <- import_validate(fi, hd$sep, hd$columns, mp, c(fl, list(kind = "binary", ncase = "1200", ncontrol = "3400")))
expect("import: a binary trait needs both case counts as whole numbers, and says what is lost without them", !vb1$ok && vb2$ok && any(grepl("cases and controls", vb2$warnings, fixed = TRUE)) && vb3$ok && !any(grepl("cases and controls", vb3$warnings, fixed = TRUE)))
mpx <- mp; mpx$se <- NULL; mpx$p <- NULL; vx <- import_validate(fi, hd$sep, hd$columns, mpx, fl)   # Round 77: with a P value the SE is derived
expect("import: a missing required column (SE, with no P value to derive it from) blocks the import and names it", !vx$ok && any(grepl("standard error", vx$errors, ignore.case = TRUE)))
fm <- tempfile(fileext = ".tsv"); gm <- g; gm$Frq <- pmin(gm$Frq, 0.5); utils::write.table(gm, fm, sep = "\t", quote = FALSE, row.names = FALSE)
vm <- import_validate(fm, hd$sep, hd$columns, mp, list(frq = "eaf", scope = "selected", scale = "beta", n = ""))
expect("import: frequencies that look like MAF are flagged when the user says they are effect allele frequencies", any(grepl("minor allele frequency", vm$warnings)))
f_or <- tempfile(fileext = ".csv"); writeLines(c("SNP,A1,A2,OR,SE,P", "rs1,A,G,1.1,0.02,0.01", "rs2,C,T,-0.2,0.02,0.2", "rs3,A,C,0.9,0.03,0.4"), f_or)
h_or <- import_read_head(f_or, "or.csv"); mp_or <- lapply(import_suggest(h_or$columns, h_or$rows), `[[`, "column")
v_or <- import_validate(f_or, h_or$sep, h_or$columns, mp_or, list(frq = "none", scope = "selected", scale = "or", n = ""))
expect("import: odds ratios that are zero or negative are refused with a hint about log odds", !v_or$ok && any(grepl("log odds", v_or$errors)))
f_pal <- tempfile(fileext = ".csv"); writeLines(c("rsid,effect_allele,other_allele,beta,se", "rs1,A,T,0.1,0.02", "rs2,C,G,0.2,0.02", "rs3,A,G,0.1,0.03"), f_pal)
h_pal <- import_read_head(f_pal, "pal.csv"); mp_pal <- lapply(import_suggest(h_pal$columns, h_pal$rows), `[[`, "column")
v_pal <- import_validate(f_pal, h_pal$sep, h_pal$columns, mp_pal, list(frq = "none", scope = "selected", scale = "beta", n = ""))
expect("import: palindromic variants without a frequency column are reported", v_pal$ok && any(grepl("2 palindromic variants", v_pal$warnings)))
set.seed(5); f_un <- tempfile(fileext = ".csv"); writeLines(c("a,b,c,d,e", paste("rs", 1000 + 1:20, ",A,G,", round(runif(20, 0.005, 0.02), 4), ",", signif(runif(20, 0.2, 0.9), 3), sep = "")), f_un)
h_un <- import_read_head(f_un, "un.csv"); sg_un <- import_suggest(h_un$columns, h_un$rows)
expect("import: a column of small standard errors is not guessed to be the P value", !identical(sg_un$p$column, "d"))
expect("import: Z scores are refused without a Z score column, the sample size and a frequency (Round 77: accepted with them)", !import_validate(fi, hd$sep, hd$columns, mp, list(frq = "eaf", scope = "selected", scale = "z", n = ""))$ok)
fx_ <- tempfile(fileext = ".txt"); utils::write.table(data.frame(FID = 1, IID = 2, PHENO = 1, rs1 = 0, rs2 = 1), fx_, quote = FALSE, row.names = FALSE)
expect("import: individual-level data is refused before anything is read", grepl("individual-level", err_of(import_read_head(fx_))))
expect("import: VCF and zip files are refused with a clear message", { m <- err_of(import_read_head(fi, "x.vcf.gz")); grepl("not supported here", m) && grepl("Excel", m) && grepl("GWAS-VCF files have their own import", m) })
rec <- import_save(fi, "bolt.tsv", hd$sep, hd$columns, mp, list(name = "Import test", role = "exposure", frq = "eaf", scale = "beta", scope = "selected", n = ""))
sd <- file.path(library_dir(), "Import test")
expect("import: standard copy and recipe are saved, with the checksum of the original", file.exists(file.path(sd, "standardised.parquet")) && !file.exists(file.path(sd, "standardised.csv")) && identical(rec$original$sha256, sha256_file(fi)) && identical(rec$usable_in_runs, TRUE) && identical(rec$usable_as_exposure, FALSE))
expect("import: the original file is unchanged", identical(sha256_file(fi), rec$original$sha256))
expect("import: saving the same name twice is refused", grepl("already exists", err_of(import_save(fi, "bolt.tsv", hd$sep, hd$columns, mp, list(name = "Import test", frq = "eaf", scale = "beta", scope = "selected", n = "")))))
expect("import: the library lists the dataset", length(library_list()) == 1L && identical(library_list()[[1]]$trait, "Import test"))
std <- lib_table_read(sd, "standardised"); expect("import: the standard copy has the standard columns", all(c("snp", "ea", "oa", "beta", "se") %in% names(std)) && nrow(std) == ni)
# Round 35: library tables are Parquet; a dataset saved before as CSV is converted on first read, exactly, and the CSV removed.
lt <- file.path(tempdir(), "r35_libtable"); dir.create(lt, showWarnings = FALSE)
d35 <- data.frame(snp = c("rs1", "rs2", "rs3"), chr = c("1", "X", "22"), ea = c("A", "G", ""), oa = c("C", "T", "A"), beta = c(0.1234567890123456, -2e-17, NA), se = c(0.01, 0.02, 0.03), n = c(1000L, 2000L, NA), stringsAsFactors = FALSE)
utils::write.csv(d35, file.path(lt, "standardised.csv"), row.names = FALSE, na = "")
c35 <- utils::read.csv(file.path(lt, "standardised.csv"), stringsAsFactors = FALSE, colClasses = c(snp = "character", chr = "character", ea = "character", oa = "character"))
r35 <- lib_table_read(lt, "standardised")
expect("library: an old CSV table is read as before and converted to Parquet", isTRUE(all.equal(r35, c35, check.attributes = FALSE)) && file.exists(file.path(lt, "standardised.parquet")) && !file.exists(file.path(lt, "standardised.csv")))
expect("library: the converted table reads back identical, with character columns kept as text", isTRUE(all.equal(lib_table_read(lt, "standardised"), r35, check.attributes = FALSE)) && is.character(r35$chr) && identical(r35$ea[3], ""))
lib_table_write(d35, lt, "exact"); e35 <- lib_table_read(lt, "exact")
expect("library: Parquet keeps every number exactly (no rounding to 15 digits)", identical(e35$beta[1], 0.1234567890123456) && identical(e35$beta[2], -2e-17) && is.na(e35$beta[3]) && is.na(e35$n[3]))
expect("library: a missing table gives a clear message", grepl("is not in the Data library any more", tryCatch(lib_table_read(lt, "none"), error = conditionMessage)))
unlink(lt, recursive = TRUE)
template_save("Lab layout", hd$fingerprint, mp, list(role = "exposure"))
expect("template: the same headings propose the template, a different layout does not", identical(template_match(hd$fingerprint)$name, "Lab layout") && is.null(template_match(c("a", "b", "c"))))
expect("template: a near-identical set of headings still matches", !is.null(template_match(c(hd$fingerprint[-1], "extra"))) == (length(hd$fingerprint) / (length(hd$fingerprint) + 1) >= 0.8))
template_delete("Lab layout"); expect("template: deleting removes it", length(template_list()) == 0L)

# ---- the page and the Shiny server ---------------------------------------------------------------------------------------
html <- as.character(mr_ui_for(args_dir)(NULL))
expect("page: has the root element, the script and the footer the smoke test looks for", grepl('id="root"', html, fixed = TRUE) && grepl("MR GUI engine prototype", html, fixed = TRUE) && grepl('id="hello"', html, fixed = TRUE))
expect("page: no external script, only the font stylesheet is remote", !grepl("<script[^>]+src=\"https?://", html))
options(mrgui.engine_dir = args_dir)   # Round 83: background jobs load the engine from here
shiny::testServer(mr_server(args_dir), {
  session$setInputs(act = list(t = "init", nonce = 1))
  session$setInputs(act = list(t = "study_create", name = "", aim = "", kind = "basic", nonce = 2))
  expect("server: a blank study name does not crash and creates nothing", is.null(rv$h))
  session$setInputs(act = list(t = "study_create", name = "Server study", aim = "a", kind = "basic", prereg = FALSE, nonce = 3))
  expect("server: creating a study sets the current study", !is.null(rv$h) && grepl("Server study$", rv$h$dir))
  # Round 83: runs and reanalyses are background jobs; the test waits for them as the app does (passing on their messages)
  bg_wait <- function(n = 3000L) for (k in seq_len(n)) { later::run_now(0.1); if (!length(bg$jobs)) break }
  session$setInputs(act = list(t = "run_pair", i = 0, ex = "expA", out = "outB", methods = list("ivw_mre", "egger"), ci = 0.9, nonce = 4)); bg_wait()
  pp <- ui_state_get(rv$h$con, "pairs")
  expect("server: a pair run is recorded in the study database with its labels and folder", length(pp) == 1L && pp[[1]]$ex == "expA" && grepl("^P01_", pp[[1]]$folder))
  session$setInputs(act = list(t = "run_pair", i = 1, ex = "expB", out = "outB", methods = list("ivw_mre"), ci = 0.95, nonce = 5)); bg_wait()
  pf <- studio_pair_folders(rv$h$dir)
  expect("server: a second pair gets its own folder", length(pf) == 2L && grepl("^P02_", pf[2]))
  session$setInputs(act = list(t = "run_pair", i = 0, ex = "expA", out = "outB", methods = list("ivw_mre", "egger"), ci = 0.9, nonce = 6)); bg_wait()
  expect("server: rerunning a pair with the same settings reuses its folder", length(studio_pair_folders(rv$h$dir)) == 2L && length(ui_state_get(rv$h$con, "pairs")) == 2L)
  # Round 75f: the engine works through the queue itself (the window used to ask for each pair, so a locked screen stopped the run)
  pump <- function(n = 3000L) for (k in seq_len(n)) { later::run_now(0.1); if (is.null(rq$pairs) && !length(bg$jobs)) break }
  q2 <- list(list(i = 2, ex = "expA", out = "outA", methods = list("ivw_mre"), ci = 0.95), list(i = 3, ex = "expB", out = "outA", methods = list("ivw_mre"), ci = 0.95))
  session$setInputs(act = list(t = "run_all", pairs = q2, total = 4L, study = "Server study", nonce = 7)); pump()
  expect("Round 75f: a queued run finishes every pair by itself, with no request from the window between pairs", is.null(rq$pairs) && length(studio_pair_folders(rv$h$dir)) == 4L && rq$ok == 2L)
  q3 <- list(list(i = 4, ex = "expA", out = "outC", methods = list("ivw_mre"), ci = 0.95), list(i = 5, ex = "expB", out = "outC", methods = list("ivw_mre"), ci = 0.95))
  session$setInputs(act = list(t = "run_all", pairs = q3, total = 6L, nonce = 8)); session$setInputs(act = list(t = "run_pause", on = TRUE, nonce = 9)); pump(20L)
  n_paused <- length(studio_pair_folders(rv$h$dir))
  session$setInputs(act = list(t = "run_pause", on = FALSE, nonce = 10)); later::run_now(0.2)
  session$setInputs(act = list(t = "run_stop", nonce = 11)); pump()
  expect("Round 75f: Pause holds the queue between pairs, Resume continues it, and Stop ends it after the pair that is running",
         n_paused == 4L && is.null(rq$pairs) && length(studio_pair_folders(rv$h$dir)) == 5L)
  # Round 75g: a mediation step (exposure -> mediator) is rebuilt from its saved record; its plain labels name the mediator as the outcome,
  # which is not an outcome dataset (test 3: "TG was imported as an exposure only" in instrument review and the outlier analysis)
  med <- list(role = "xm", m = "expC", x = "expA", ms = c("expC", "expD"), y = "outA"); mi <- mediation_pair_input(med)
  r7 <- studio_run(rv$h, mi$exposure_label, mi$outcome_label, mi$input, c("ivw_mre", "egger"), idx = 6L, params = list(ci_level = 0.95))
  ui_state_set(rv$h$con, "pairs", c(ui_state_get(rv$h$con, "pairs") %||% list(), list(list(ex = "expA", out = "expC", med = med, folder = r7$pair_folder))))
  a7 <- list(folder = r7$pair_folder, ex = "expA", out = "expC", pthr = 5e-8)
  plain <- err_of(instrument_table("expA", "expC")); it7 <- instrument_table("expA", "expC", inp = pair_inp(a7))
  n_before <- length(studio_pair_folders(rv$h$dir))
  session$setInputs(act = list(t = "run_sens", folder = r7$pair_folder, ex = "expA", out = "expC", pthr = 5e-8, remove = list(it7$id[1]), why = "chosen by you", methods = list("ivw_mre", "egger"), nonce = 12)); bg_wait()
  expect("Round 75g: instrument review and sensitivity reruns of a mediation step use its own datasets (its plain labels cannot build it)",
         !is.na(plain) && nrow(it7) == sum(mi$input$data$mr_keep) && length(studio_pair_folders(rv$h$dir)) == n_before + 1L)
})
# ---- imported datasets in a run: harmonisation, ledger, refusals ----
set.seed(1); kk <- 60; snps <- sprintf("rs%d", 1000 + 1:kk); eaa <- sample(c("A", "C", "G", "T"), kk, TRUE); cmp <- c(A = "T", T = "A", C = "G", G = "C")
oaa <- vapply(eaa, function(x) sample(setdiff(c("A", "C", "G", "T"), x), 1), ""); bxx <- runif(kk, .05, .15) * sample(c(-1, 1), kk, TRUE); sxx <- runif(kk, .01, .015)
efq <- round(runif(kk, .1, .9), 3); byy <- 0.3 * bxx + rnorm(kk, 0, .01); syy <- runif(kk, .01, .02)
exd <- data.frame(snp = snps, ea = eaa, oa = oaa, eaf = efq, beta = bxx, se = sxx, p = 2 * pnorm(-abs(bxx / sxx)), n = 100000)
oud <- data.frame(snp = snps, ea = eaa, oa = oaa, eaf = efq, beta = byy, se = syy, p = 2 * pnorm(-abs(byy / syy)), n = 200000, stringsAsFactors = FALSE)
swr <- 1:10; tt <- oud$ea[swr]; oud$ea[swr] <- oud$oa[swr]; oud$oa[swr] <- tt; oud$beta[swr] <- -oud$beta[swr]; oud$eaf[swr] <- 1 - oud$eaf[swr]
palx <- (eaa == "A" & oaa == "T") | (eaa == "T" & oaa == "A") | (eaa == "C" & oaa == "G") | (eaa == "G" & oaa == "C")
str5 <- which(!palx & seq_len(kk) > 10)[1:5]; oud$ea[str5] <- cmp[oud$ea[str5]]; oud$oa[str5] <- cmp[oud$oa[str5]]
oud <- oud[-(55:58), ][sample(kk - 4), ]
fxx <- tempfile(fileext = ".csv"); foo <- tempfile(fileext = ".csv"); utils::write.csv(exd, fxx, row.names = FALSE); utils::write.csv(oud, foo, row.names = FALSE)
mpl <- list(snp = "snp", ea = "ea", oa = "oa", eaf = "eaf", beta = "beta", se = "se", p = "p", n = "n")
fct <- function(nm, role, ind, ea_ok = "yes") list(name = nm, role = role, kind = "continuous", scale = "beta", ea_ok = ea_ok, frq = "eaf", build = "37", ancestry = "EUR", n = "", scope = "selected", sex = "both", ukb = FALSE, independent = ind)
invisible(import_save(fxx, "ex.csv", ",", names(exd), mpl, fct("Trait X", "exposure", "yes")))
invisible(import_save(foo, "out.csv", ",", names(exd), mpl, fct("Trait Y", "outcome", "no")))
invisible(import_save(fxx, "ex2.csv", ",", names(exd), mpl, fct("Trait X2", "exposure", "no")))
expect("library catalogue lists only usable datasets by role", identical(vapply(library_catalogue()$exposures, `[[`, "", "id"), "lib:Trait X") && identical(vapply(library_catalogue()$outcomes, `[[`, "", "id"), "lib:Trait Y"))
lp <- pair_input("lib:Trait X", "lib:Trait Y", pthr = 1e-4); led <- lp$ledger; ld <- lp$input$data
expect("harmonisation ledger: 4 instruments missing in the outcome, every other one kept, flipped or dropped with a reason", led$missing_in_outcome == 4 && led$kept + led$dropped_palindromic + led$dropped_mismatch == led$in_outcome && led$flipped >= 10 && !lp$simulated)
expect("harmonised effects point the same way as the truth after allele alignment", mean(sign(ld$beta_x[ld$mr_keep]) == sign(ld$beta_y[ld$mr_keep])) > 0.9)
expect("an exposure that is not marked independent is refused with the reason", grepl("set to be clumped on this computer", tryCatch(pair_input("lib:Trait X2", "lib:Trait Y", 1e-4), error = conditionMessage)))
expect("a simulated and an imported trait cannot be paired", grepl("cannot be paired", tryCatch(pair_input("expA", "lib:Trait Y"), error = conditionMessage)))
expect("an exposure with no variant under the threshold is refused", grepl("Fewer than 2 variants", tryCatch(pair_input("lib:Trait X", "lib:Trait Y", 1e-300), error = conditionMessage)))
hl <- studio_new("Imported run"); rl <- studio_run(hl, lp$exposure_label, lp$outcome_label, lp$input, c("ivw_mre", "egger"), idx = 1L)
pll <- pair_payload(hl$dir, rl$pair_folder)
expect("an imported pair runs through both methods with an identical replay check", identical(pll$replay$status, "identical") && length(pll$estimates) == 3L && all(vapply(pll$estimates, function(e) e$n_ivs_used == led$kept, TRUE)))
ivw_l <- suppressWarnings(TwoSampleMR::mr_ivw_mre(ld$beta_x[ld$mr_keep], ld$beta_y[ld$mr_keep], ld$se_x[ld$mr_keep], ld$se_y[ld$mr_keep], TwoSampleMR::default_parameters()))
expect("the imported IVW estimate equals TwoSampleMR on the harmonised table", abs(pll$estimates[[1]]$estimate - ivw_l$b) < 1e-9)
# ---- Round 9: publication package from an imported pair ----
ui_state_set(hl$con, "pairs", list(list(ex = "lib:Trait X", out = "lib:Trait Y", exName = lp$exposure_label, outName = lp$outcome_label, folder = rl$pair_folder, methods = as.list(c("ivw_mre", "egger")), ci = 0.95, pthr = 1e-4, simulated = FALSE, ledger = lp$ledger, anc_x = "EUR", anc_y = "EUR")))
pk <- build_package(hl, c("scripts", "report", "figs", "tables", "text", "refs", "capsule"), 140)
pkd <- pk$dir; zl <- function(f) utils::unzip(file.path(pkd, f), list = TRUE)$Name
expect("package: every promised file exists", all(file.exists(file.path(pkd, c("report.docx", "report.html", "supplement_tables.xlsx", "draft_methods_results_limitations.docx", "references.bib", "references.ris", "reproducibility_capsule.zip", "README.txt", "figures/P1_forest.pdf", "figures/P1_scatter.pdf")))))
expect("package: the Word, Excel and zip files are valid archives", all(c("word/document.xml", "word/styles.xml") %in% zl("report.docx")) && "xl/worksheets/sheet1.xml" %in% zl("supplement_tables.xlsx") && "SESSION.txt" %in% zl("reproducibility_capsule.zip") && any(grepl("02_scripts/00_run_all.R$", zl("reproducibility_capsule.zip"))))
expect("package: the CRC-32 of 123456789 is the standard check value and a zip entry round-trips byte for byte", identical(crc32_raw(charToRaw("123456789")), -873187034L) && { zf <- tempfile(fileext = ".zip"); pay <- as.raw(c(0:255, 255:0)); zip_write(zf, list(`a/b.bin` = pay)); ex <- tempfile(); utils::unzip(zf, exdir = ex); identical(readBin(file.path(ex, "a", "b.bin"), "raw", 600L), pay) })
ivw1 <- Filter(function(e) e$method == "ivw_mre", pll$estimates)[[1]]; htm <- paste(readLines(file.path(pkd, "report.html"), warn = FALSE), collapse = "\n")
expect("package: the report quotes the saved IVW estimate and the harmonisation ledger", grepl(sprintf("beta = %.3f", ivw1$estimate), htm, fixed = TRUE) && grepl(sprintf("%d instruments were used", led$kept), htm, fixed = TRUE) && grepl(sprintf("%d were absent from the outcome file", led$missing_in_outcome), htm, fixed = TRUE))
expect("package: an imported pair has no simulated warning and says clumping was not done by the app", !grepl("SIMULATED DATA", htm, fixed = TRUE) && grepl("did not perform clumping", htm, fixed = TRUE))
bib <- paste(readLines(file.path(pkd, "references.bib"), warn = FALSE), collapse = "\n")
expect("package: references are complete and every citation in the text is listed", all(vapply(c("Hemani2018", "Burgess2013", "Bowden2015", "Burgess2023", "Skrivankova2021"), grepl, TRUE, x = bib, fixed = TRUE)) && !grepl("Bowden2016|Hartwig2017", bib))
expect("package: limitations name what was not tested", grepl("Sample overlap between the exposure and outcome studies was not tested", htm, fixed = TRUE) && grepl("Steiger directionality test was not run", htm, fixed = TRUE))
# ---- sample overlap: the analyst's declaration is saved with the study and used in the limitations --------------------------------------
ovh <- studio_new("Overlap study")
expect("overlap: nothing saved yet gives an empty answer list", length(overlap_get(ovh)) == 0L)
overlap_set(ovh, "lib:BMI", "lib:CHD", "some", 35, "  both use UK Biobank\n ")
ov1 <- overlap_get(ovh)[[overlap_key("lib:BMI", "lib:CHD")]]
expect("overlap: answer, percentage and note are saved, the note is cleaned", identical(ov1$answer, "some") && ov1$pct == 35 && identical(ov1$note, "both use UK Biobank"))
expect("overlap: a readable file is written and registered in the manifest", file.exists(file.path(ovh$dir, "01_study_plan", "sample_overlap.json")) && any(grepl("sample_overlap", utils::read.csv(file.path(ovh$dir, "MANIFEST.csv"), stringsAsFactors = FALSE)$role)))
expect("overlap: the answer survives closing and opening the study", {DBI::dbDisconnect(ovh$con); ovh <- studio_open(basename(ovh$dir)); identical(overlap_get(ovh)[[overlap_key("lib:BMI", "lib:CHD")]]$answer, "some")})
expect("overlap: an unknown answer value is refused", grepl("Choose one of", err_of(overlap_set(ovh, "lib:A", "lib:B", "maybe"))))
expect("overlap: a share outside 0 to 100 is refused", grepl("between 0 and 100", err_of(overlap_set(ovh, "lib:A", "lib:B", "some", 120))))
overlap_set(ovh, "lib:A", "lib:B", "none", 50, "")
expect("overlap: a share is dropped unless the answer is 'some'", is.null(overlap_get(ovh)[[overlap_key("lib:A", "lib:B")]]$pct))
mkc <- function(ex, out, ov, sim = FALSE) list(simulated = sim, p = list(ex = ex, exName = ex, out = out, outName = out), overlap = ov)
lim0 <- overlap_limitations(list(mkc("A", "B", NULL, sim = TRUE)))
expect("overlap text: simulated traits say their samples do not overlap (Round 102)", length(lim0) == 1L && grepl("do not overlap", lim0, fixed = TRUE) && !grepl("analyst stated", lim0))
lim1 <- overlap_limitations(list(mkc("X", "Y", list(answer = "none", note = "different biobanks"))))
expect("overlap text: 'none' is reported as the analyst's unverified statement", length(lim1) == 1L && grepl("analyst stated", lim1) && grepl("do not share participants", lim1) && grepl("did not verify", lim1) && grepl("different biobanks", lim1))
lim2 <- overlap_limitations(list(mkc("X", "Y", list(answer = "some", pct = 40, note = "UKB"))))
expect("overlap text: 'some' names the share, the bias and that nothing was corrected", grepl("share participants", lim2) && grepl("40%", lim2) && grepl("did not correct", lim2) && grepl("biases estimates", lim2))
lim3 <- overlap_limitations(list(mkc("X", "Y", NULL), mkc("P", "Q", list(answer = "unknown"))))
expect("overlap text: unanswered and 'I do not know' give the not-tested sentence naming the pairs", length(lim3) == 1L && grepl("was not tested", lim3) && grepl("X and Y", lim3) && grepl("P and Q", lim3))
lim4 <- overlap_limitations(list(mkc("X", "Y", list(answer = "none")), mkc("P", "Q", list(answer = "some"))))
expect("overlap text: pairs with different answers get one sentence per answer", length(lim4) == 2L)
fo <- pair_facts("expA", "outA"); expect("overlap: the pair facts carry the pair ids the answers are keyed by", is.character(fo$ex_id) && is.character(fo$out_id))
DBI::dbDisconnect(ovh$con)
pk2 <- build_package(hl, c("refs"), 140); expect("package: a second build never overwrites the first", !identical(pk2$dir, pkd) && file.exists(file.path(pkd, "report.html")))
expect("package: asking for a package before any result exists is refused", grepl("no saved results", tryCatch(build_package(studio_new("Empty report"), c("refs"), 140), error = conditionMessage)))
DBI::dbDisconnect(hl$con)
# "the effect size is for the other allele" must swap the labels in the saved copy
swd <- data.frame(snp = c("rs1", "rs2", "rs3"), ea = c("A", "C", "G"), oa = c("G", "T", "A"), beta = c(.1, .2, .3), se = c(.01, .01, .01)); fsw <- tempfile(fileext = ".csv"); utils::write.csv(swd, fsw, row.names = FALSE)
invisible(import_save(fsw, "sw.csv", ",", names(swd), list(snp = "snp", ea = "ea", oa = "oa", beta = "beta", se = "se"), modifyList(fct("Swap test", "outcome", "no", ea_ok = "swap"), list(frq = "none"))))
sd2 <- library_dataset("Swap test")$data
expect("import: answering 'no, the effect is for the other allele' swaps the allele labels", identical(sd2$ea, c("G", "T", "A")) && identical(sd2$oa, c("A", "C", "G")))
# ---- Round 8: weighted median, simple median, weighted mode ----
p8 <- pair_input("expA", "outC"); d8 <- p8$input$data; d8 <- d8[d8$mr_keep, ]
h8 <- studio_new("Robust methods")
r8 <- studio_run(h8, p8$exposure_label, p8$outcome_label, p8$input, c("ivw_mre", "wm", "smed", "mode"), idx = 1L, params = list(ci_level = 0.9, nboot = 300L, phi = 0.5))
pl8 <- pair_payload(h8$dir, r8$pair_folder); by8 <- function(m) Filter(function(e) e$method == m, pl8$estimates)[[1]]
expect("robust methods: all run and the replay check is identical", all(vapply(r8$results, function(x) identical(x$status, "ok"), TRUE)) && identical(pl8$replay$status, "identical"))
tp <- TwoSampleMR::default_parameters(); tp$nboot <- 300; tp$phi <- 0.5
dm <- function(f) f(d8$beta_x, d8$beta_y, d8$se_x, d8$se_y, tp)
expect("weighted and simple median point estimates equal TwoSampleMR", abs(by8("wm")$estimate - dm(TwoSampleMR::mr_weighted_median)$b) < 1e-12 && abs(by8("smed")$estimate - dm(TwoSampleMR::mr_simple_median)$b) < 1e-12)
expect("bootstrap standard errors are close to TwoSampleMR's with another seed", abs(by8("wm")$se / dm(TwoSampleMR::mr_weighted_median)$se - 1) < 0.35)
expect("robust methods use the confidence level and write the bootstrap size", by8("wm")$ci_level == 0.9 && abs((by8("wm")$ci_upper - by8("wm")$ci_lower) / 2 / by8("wm")$se - stats::qnorm(0.95)) < 1e-9 && any(vapply(pl8$diagnostics, function(g) identical(g$name, "bootstrap_replications") && g$value == 300, TRUE)))
expect("phi is read back from the scripts and changes the mode estimate", identical(pair_settings(h8$dir, r8$pair_folder)$phi, 0.5) && identical(pair_settings(h8$dir, r8$pair_folder)$nboot, 300) &&
       abs(by8("mode")$estimate - TwoSampleMR::mr_weighted_mode(d8$beta_x, d8$beta_y, d8$se_x, d8$se_y, modifyList(tp, list(phi = 1)))$b) > 0 &&
       abs(by8("mode")$estimate - TwoSampleMR::mr_weighted_mode(d8$beta_x, d8$beta_y, d8$se_x, d8$se_y, tp)$b) < 0.05)
r8b <- studio_run(h8, p8$exposure_label, p8$outcome_label, p8$input, c("wm"), idx = 1L, params = list(ci_level = 0.9, nboot = 300L, phi = 0.5))
expect("rerunning with the same settings reuses the saved results", isTRUE(r8b$reused))
expect("the parameter schemas of the new methods are complete", identical(names(param_schemas()$wm), c("ci_level", "nboot", "min_f")) && identical(names(param_schemas()$mode), c("ci_level", "nboot", "phi", "min_f")) && identical(names(param_schemas()$smed), c("ci_level", "nboot", "min_f")))
expect("a bad bootstrap size is refused by the adapter", grepl("nboot below minimum 100", tryCatch(studio_run(h8, "A", "B", p8$input, "wm", idx = 2L, params = list(nboot = 5L)), error = conditionMessage)))
ui_state_set(h8$con, "pairs", list(list(ex = "expA", out = "outC", exName = p8$exposure_label, outName = p8$outcome_label, folder = r8$pair_folder, methods = as.list(c("ivw_mre", "wm", "smed", "mode")), ci = 0.9, pthr = 5e-8, simulated = TRUE, ledger = NULL, anc_x = "EUR", anc_y = "EUR")))
pks <- build_package(h8, c("report", "text", "refs"), 89); hs <- paste(readLines(file.path(pks$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("package: simulated data carry a warning in the report and the README", grepl("SIMULATED DATA: not for publication", hs, fixed = TRUE) && pks$simulated && any(grepl("SIMULATED", readLines(file.path(pks$dir, "README.txt")))))
expect("package: the methods text uses the settings of the run (bootstrap size, bandwidth, confidence level)", grepl("300 replications", hs, fixed = TRUE) && grepl("phi = 0.5", hs, fixed = TRUE) && grepl("90% confidence intervals", hs, fixed = TRUE))
expect("package: weighted median, simple median and weighted mode are cited", grepl("Bowden2016", paste(readLines(file.path(pks$dir, "references.bib")), collapse = "\n"), fixed = TRUE) && grepl("Hartwig2017", paste(readLines(file.path(pks$dir, "references.bib")), collapse = "\n"), fixed = TRUE))
DBI::dbDisconnect(h8$con)
# ---- Round 10: simple mode and outlier detectors in the app ----
h10 <- studio_new("Detectors")
r10 <- studio_run(h10, p8$exposure_label, p8$outcome_label, p8$input, c("ivw_mre", "smode", "presso", "radial", "radegger"), idx = 1L,
                  params = list(ci_level = 0.95, nboot = 200L, phi = 0.5, nbdist = 500L, signif = 0.05, alpha = 0.01, weights = "3"))
pl10 <- pair_payload(h10$dir, r10$pair_folder); by10 <- function(m) Filter(function(e) e$method == m, pl10$estimates)[[1]]
expect("detectors: all methods run and replay identically", all(vapply(r10$results, function(x) identical(x$status, "ok"), TRUE)) && identical(pl10$replay$status, "identical"))
sm10 <- pair_settings(h10$dir, r10$pair_folder)
expect("detectors: their settings are read back from the scripts", identical(sm10$nbdist, 500) && identical(sm10$alpha, 0.01) && identical(sm10$weights, "3") && identical(sm10$signif, 0.05))
expect("simple mode equals TwoSampleMR with the same settings (within bootstrap noise of the point estimate)", abs(by10("smode")$estimate - TwoSampleMR::mr_simple_mode(d8$beta_x, d8$beta_y, d8$se_x, d8$se_y, modifyList(TwoSampleMR::default_parameters(), list(nboot = 200, phi = 0.5)))$b) < 1e-9)
rd <- RadialMR::format_radial(d8$beta_x, d8$beta_y, d8$se_x, d8$se_y, d8$variant_id); rv <- suppressWarnings(RadialMR::ivw_radial(rd, 0.01, 3, 1e-4, FALSE))
expect("radial IVW estimate equals RadialMR called directly", abs(by10("radial")$estimate - rv$coef[1, 1]) < 1e-12 && abs(by10("radial")$se - rv$coef[1, 2]) < 1e-12)
expect("detectors write the flagged variants as a list the screen can read", any(vapply(pl10$diagnostics, function(g) g$method == "presso" && g$name == "presso_outliers", TRUE)) && any(vapply(pl10$diagnostics, function(g) g$method == "radial" && g$name == "radial_outliers", TRUE)))

# ---- Round 15: Steiger directionality (diagnostic-only method) ----
h15 <- studio_new("Steiger")
r15 <- studio_run(h15, p8$exposure_label, p8$outcome_label, p8$input, c("ivw_mre", "steiger"), idx = 1L, params = list(ci_level = 0.95, r_xxo = 0.9))
pl15 <- pair_payload(h15$dir, r15$pair_folder); sg <- function(n) Filter(function(g) g$method == "steiger" && g$name == n, pl15$diagnostics)[[1]]
expect("steiger: runs next to IVW and replays identically", all(vapply(r15$results, function(x) identical(x$status, "ok"), TRUE)) && identical(pl15$replay$status, "identical"))
expect("steiger: gives diagnostics and no estimate", !any(vapply(pl15$estimates, function(e) e$method == "steiger", TRUE)) && length(Filter(function(g) g$method == "steiger", pl15$diagnostics)) == 8L)
st15 <- TwoSampleMR::mr_steiger(rep(NA_real_, nrow(d8)), rep(NA_real_, nrow(d8)), d8$n_x, d8$n_y, TwoSampleMR::get_r_from_bsen(d8$beta_x, d8$se_x, d8$n_x), TwoSampleMR::get_r_from_bsen(d8$beta_y, d8$se_y, d8$n_y), r_xxo = 0.9, r_yyo = 1)
expect("steiger: numbers equal TwoSampleMR::mr_steiger called directly", abs(sg("steiger_r2_exposure")$value - st15$r2_exp) < 1e-12 && abs(sg("steiger_r2_outcome")$value - st15$r2_out) < 1e-12 && abs(sg("steiger_direction_adjusted")$p - st15$steiger_test_adj) <= 1e-12 * max(1, st15$steiger_test_adj) + 1e-300)
expect("steiger: the setting is read back from the script", identical(pair_settings(h15$dir, r15$pair_folder)$r_xxo, 0.9))
res15 <- tryCatch(pair_results(h15$dir, r15$pair_folder), error = function(e) NULL)
expect("steiger: the results table skips it", !is.null(res15) && !any(grepl("Steiger", res15$Method)))
expect("steiger: refused for a binary outcome", inherits(tryCatch(studio_run(h15, p8$exposure_label, p8$outcome_label, modifyList(p8$input, list(outcome_scale = "log_or")), "steiger", idx = 2L), error = function(e) e), "error"))
expect("steiger: is listed as a method with a schema", "steiger" %in% METHOD_CHOICES && identical(names(param_schemas()$steiger), c("r_xxo", "r_yyo", "min_f")))
expect("steiger: pair facts say sample sizes are present", isTRUE(pair_facts("expA", "outA")$has_n))
ui_state_set(h15$con, "pairs", list(list(ex = "expA", out = "outA", exName = p8$exposure_label, outName = p8$outcome_label, folder = r15$pair_folder, methods = as.list(c("ivw_mre", "steiger")), ci = 0.95, pthr = 5e-8, simulated = TRUE)))
pk15 <- build_package(h15, c("scripts", "report", "text", "refs"), 140)
htm15 <- paste(readLines(file.path(pk15$dir, "report.html"), warn = FALSE), collapse = "\n"); bib15 <- paste(readLines(file.path(pk15$dir, "references.bib"), warn = FALSE), collapse = "\n")
expect("steiger: the report describes the test, quotes its numbers and cites Hemani 2017", grepl("Steiger directionality test, run together with the IVW analysis, compared", htm15, fixed = TRUE) && grepl(sprintf("explained %s of the variance in", .fnum(sg("steiger_r2_exposure")$value, 3)), htm15, fixed = TRUE) && grepl("Hemani2017", bib15, fixed = TRUE))
expect("steiger: the limitations no longer say directionality was not tested", !grepl("Steiger directionality test was not run", htm15, fixed = TRUE) && grepl("does not exclude confounding or reverse causation", htm15, fixed = TRUE))
# ---- Round 16: Steiger as part of IVW, with binary traits (liability scale) ----
expect("steiger16: pair facts report readiness for a continuous pair", identical(pair_facts("expA", "outA")$steiger$why, "") && !pair_facts("expA", "outA")$steiger$bin_y && identical(pair_facts("expA", "outA")$scale_x, "beta"))
pb <- p8; pb$input$outcome_scale <- "log_or"; pb$input$data$eaf_y <- 0.3; pb$binary <- list(x = NULL, y = list(ncase = 25000, ncontrol = 75000))
expect("steiger16: a continuous pair keeps its analysis file unchanged", identical(names(steiger_input(p8, list(y = 0.05))$data), names(p8$input$data)))
sib <- steiger_input(pb, list(y = 0.05))
expect("steiger16: a binary outcome gets case, control and prevalence columns", all(c("ncase_y", "ncontrol_y", "prevalence_y") %in% names(sib$data)) && all(sib$data$prevalence_y == 0.05) && !("ncase_x" %in% names(sib$data)))
rd16 <- steiger_readiness(pb)
expect("steiger16: readiness of a binary outcome does not depend on the prevalence yet", identical(rd16$why, "") && isTRUE(rd16$bin_y) && !isTRUE(rd16$bin_x))
pbn <- pb; pbn$binary$y$ncase <- NA_real_
expect("steiger16: a binary outcome without case counts is not possible, with a clear reason", grepl("cases and controls", steiger_readiness(pbn)$why, fixed = TRUE))
h16 <- studio_new("Steiger binary")
expect("steiger16: refused without a prevalence", grepl("prevalence", conditionMessage(tryCatch(studio_run(h16, p8$exposure_label, "Binary outcome", pb$input, c("ivw_mre", "steiger"), idx = 1L, steiger_in = steiger_input(pb, list())), error = function(e) e)), fixed = TRUE))
r16 <- studio_run(h16, p8$exposure_label, "Binary outcome", pb$input, c("ivw_mre", "steiger"), idx = 1L, steiger_in = sib)
pl16 <- pair_payload(h16$dir, r16$pair_folder); sg16 <- function(n) Filter(function(g) g$method == "steiger" && g$name == n, pl16$diagnostics)[[1]]
d16 <- sib$data[sib$data$mr_keep, ]
st16 <- TwoSampleMR::mr_steiger(rep(NA_real_, nrow(d16)), rep(NA_real_, nrow(d16)), d16$n_x, rep(1e5, nrow(d16)), TwoSampleMR::get_r_from_bsen(d16$beta_x, d16$se_x, d16$n_x),
          TwoSampleMR::get_r_from_lor(d16$beta_y, d16$eaf_y, 25000, 75000, 0.05, model = "logit", correction = FALSE), 1, 1)
expect("steiger16: binary run replays identically and equals TwoSampleMR called directly", identical(pl16$replay$status, "identical") && abs(sg16("steiger_r2_outcome")$value - st16$r2_out) < 1e-12 && sg16("steiger_outcome_binary")$value == 1)
expect("steiger16: IVW's own analysis file has no Steiger columns", !any(c("ncase_y", "prevalence_y") %in% names(utils::read.csv(file.path(h16$dir, "03_pairs", r16$pair_folder, "01_data", "analysis_data_ivw_mre.csv")))))
ui_state_set(h16$con, "pairs", list(list(ex = "expA", out = "outC", exName = p8$exposure_label, outName = "Binary outcome", folder = r16$pair_folder, methods = as.list(c("ivw_mre", "steiger")), ci = 0.95, pthr = 5e-8, simulated = TRUE, prev = list(y = 0.05))))
pk16 <- build_package(h16, c("report", "text", "refs"), 140)
htm16 <- paste(readLines(file.path(pk16$dir, "report.html"), warn = FALSE), collapse = "\n"); bib16 <- paste(readLines(file.path(pk16$dir, "references.bib"), warn = FALSE), collapse = "\n")
expect("steiger16: report explains the liability scale, names the prevalence and cites Lee 2012", grepl("liability scale", htm16, fixed = TRUE) && grepl("Binary outcome 0.05", htm16, fixed = TRUE) && grepl("Lee2012", bib16, fixed = TRUE) && !grepl("Lee2012", bib15, fixed = TRUE))
expect("steiger16: case counts come from the import facts or the OpenGWAS record", identical(.case_counts(list(facts = list(ncase = "1200", ncontrol = "3400")))$ncase, 1200) && identical(.case_counts(list(facts = list(), source = list(ncase = 60801, ncontrol = 123504)))$ncontrol, 123504) && is.na(.case_counts(list(facts = list()))$ncase))
expect("detector schemas are in the Customise list", identical(names(param_schemas()$presso), c("ci_level", "nbdist", "signif", "min_f")) && identical(names(param_schemas()$radial), c("ci_level", "alpha", "alpha_correction", "weights", "tol", "min_f")) && identical(names(param_schemas()$radegger), c("ci_level", "alpha", "alpha_correction", "weights", "min_f")))
ui_state_set(h10$con, "pairs", list(list(ex = "expA", out = "outC", exName = p8$exposure_label, outName = p8$outcome_label, folder = r10$pair_folder, methods = as.list(c("ivw_mre", "smode", "presso", "radial", "radegger")), ci = 0.95, pthr = 5e-8, simulated = TRUE, ledger = NULL, anc_x = "EUR", anc_y = "EUR")))
pk10 <- build_package(h10, c("report", "text", "refs", "figs"), 89); hs10 <- paste(readLines(file.path(pk10$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("package: detectors are described, their flagged variants reported, and nothing is said to be removed", grepl("MR-PRESSO (global test", hs10, fixed = TRUE) && grepl("flagged", hs10, fixed = TRUE) && grepl("not removed from the primary analysis", hs10, fixed = TRUE) && grepl("MRPRESSO (commit 3e3c92d7ed", hs10, fixed = TRUE))
expect("package: detector and mode references are cited", all(vapply(c("Verbanck2018", "Bowden2018", "Hartwig2017"), function(k) grepl(k, paste(readLines(file.path(pk10$dir, "references.bib")), collapse = "\n"), fixed = TRUE), TRUE)))
wi10 <- pair_input("expA", "outC"); sens10 <- studio_run_without(h10, "expA", "outC", wi10$input$data$variant_id[1], methods = c("ivw_mre", "presso", "radial", "radegger"), idx = 9L)
expect("sensitivity branches do not repeat the detectors", identical(sort(vapply(sens10$results, function(x) x$method, "")), "ivw_mre") || !any(c("presso", "radial", "radegger") %in% vapply(sens10$results, function(x) x$method, "")))
DBI::dbDisconnect(h10$con)
# ---- Round 11: Compare with default for every setting ----
h11 <- studio_new("Compare all")
r11 <- studio_run(h11, p8$exposure_label, p8$outcome_label, p8$input, c("ivw_mre", "egger", "wm", "mode", "radial"), idx = 1L,
                  params = list(ci_level = 0.9, nboot = 300L, phi = 0.5, alpha = 0.01, weights = "3"))
cp <- compare_defaults(h11$dir, r11$pair_folder, "expA", "outC", c("ivw_mre", "egger", "wm", "mode", "radial"))
cr <- function(m) Filter(function(x) x$method == m, cp$rows)[[1]]
h11b <- studio_new("Compare all defaults"); r11b <- studio_run(h11b, p8$exposure_label, p8$outcome_label, p8$input, c("ivw_mre", "egger", "wm", "mode", "radial"), idx = 1L, params = list(ci_level = 0.95))
pd11 <- pair_payload(h11b$dir, r11b$pair_folder); fd <- function(m) Filter(function(e) e$method == m, pd11$estimates)[[1]]
expect("compare all: methods with a changed setting beyond the confidence level are run again, the others are not", setequal(unlist(cp$rerun), c("wm", "mode", "radial")) && !cr("ivw_mre")$ran && !cr("egger")$ran && cr("wm")$ran)
expect("compare all: the default run equals a fresh run with default settings", all(vapply(c("wm", "mode", "radial"), function(m) abs(cr(m)$estimate - fd(m)$estimate) < 1e-12 && abs(cr(m)$lo - fd(m)$ci_lower) < 1e-12 && abs(cr(m)$hi - fd(m)$ci_upper) < 1e-12, TRUE)))
expect("compare all: a confidence-level-only difference gives the exact default interval (normal and t)", all(vapply(c("ivw_mre", "egger"), function(m) abs(cr(m)$lo - fd(m)$ci_lower) < 1e-9 && abs(cr(m)$hi - fd(m)$ci_upper) < 1e-9, TRUE)))
expect("compare all: each row says which setting differs", any(grepl("Bootstrap replications: 300 \\(default 1000\\)", unlist(cr("wm")$what))) && any(grepl("Weights: Modified second order", unlist(cr("radial")$what))) && any(grepl("Confidence level", unlist(cr("ivw_mre")$what))))
h11c <- studio_new("Compare radial"); r11c <- studio_run(h11c, p8$exposure_label, p8$outcome_label, p8$input, "radial", idx = 1L, params = list(ci_level = 0.9))
cc <- compare_defaults(h11c$dir, r11c$pair_folder, "expA", "outC", "radial")$rows[[1]]
expect("compare all: radial IVW with only the level changed gives its exact default interval (t with k - 1 degrees of freedom)", !cc$ran && abs(cc$lo - fd("radial")$ci_lower) < 1e-9 && abs(cc$hi - fd("radial")$ci_upper) < 1e-9)
DBI::dbDisconnect(h11c$con)
expect("compare all: the saved results and the scripts are untouched", identical(pair_settings(h11$dir, r11$pair_folder)$nboot, 300) && !any(grepl("compare", list.files(h11$dir, recursive = TRUE))))
expect("compare all: a pair with only default settings reports nothing to rerun", length(compare_defaults(h11b$dir, r11b$pair_folder, "expA", "outC", c("wm"))$rerun) == 0L)
DBI::dbDisconnect(h11$con); DBI::dbDisconnect(h11b$con)
# ---- Round 12: public GWAS from OpenGWAS (simulated OpenGWAS, no internet) ----
b64u <- function(x) gsub("[=\n]", "", chartr("+/", "-_", jsonlite::base64_enc(charToRaw(x))))
fake_jwt <- function(exp) paste(b64u('{"alg":"RS256","typ":"JWT"}'), b64u(sprintf('{"iss":"api.opengwas.io","sub":"someone@example.org","iat":1,"exp":%d}', as.integer(exp))), "c2lnbmF0dXJlLXNpZ25hdHVyZS1zaWduYXR1cmU", sep = ".")
tok_ok <- fake_jwt(as.numeric(Sys.time()) + 5 * 86400); tok_old <- fake_jwt(as.numeric(Sys.time()) - 86400)
expect("opengwas token: spaces, quotes and a Bearer prefix are tidied, text that is not a token is refused", identical(ogw_clean_token(paste0('  "Bearer ', tok_ok, '" ')), tok_ok) && is.null(ogw_clean_token("hello")) && is.null(ogw_clean_token("")) && is.null(ogw_clean_token(NULL)))
ti <- ogw_token_info(tok_ok)
expect("opengwas token: only the expiry is read (about 5 days left), the account is not returned", ti$set && ti$valid_shape && !ti$expired && abs(ti$days_left - 5) < 0.1 && !grepl("someone", paste(unlist(ti), collapse = " ")))
expect("opengwas token: an expired token is recognised and refused before any request", isTRUE(ogw_token_info(tok_old)$expired) && grepl("expired", err_of(ogw_catalogue_update(tok_old))))
expect("opengwas: no token gives a clear instruction, offline mode refuses before anything else", grepl("Paste your OpenGWAS access token", err_of(ogw_catalogue_update(NULL))) && grepl("Offline mode is on", err_of(ogw_catalogue_update(tok_ok, offline = TRUE))))

mk_snp <- function(i) sprintf("rs%d", 1000 + i)
set.seed(12); nx <- 15L; ex_t <- data.frame(rsid = mk_snp(1:nx), chr = rep(c("1", "2", "5"), 5), position = 1e6 * (1:nx), ea = rep(c("A", "C", "G"), 5), nea = rep(c("G", "T", "A"), 5), eaf = round(runif(nx, .15, .85), 3),
  beta = runif(nx, .06, .12) * sample(c(-1, 1), nx, TRUE), se = runif(nx, .005, .007), n = 300000, id = "sim-exp-1", trait = "Simulated body mass", stringsAsFactors = FALSE)
ex_t$p <- 2 * pnorm(-abs(ex_t$beta / ex_t$se))
ex_raw <- rbind(ex_t, transform(ex_t[1, ], rsid = "rs5000", ea = "AT", nea = "A"), transform(ex_t[2, ]), transform(ex_t[3, ], rsid = "rs5001", se = NA))
oy_beta <- 0.25 * ex_t$beta + rnorm(nx, 0, .004); oy_se <- runif(nx, .006, .01)
mk_out <- function(ix, proxy = FALSE) data.frame(SNP = if (proxy) mk_snp(ix) else mk_snp(ix), chr = ex_t$chr[ix], pos = ex_t$position[ix], beta.outcome = oy_beta[ix], se.outcome = oy_se[ix], samplesize.outcome = 230000, pval.outcome = 2 * pnorm(-abs(oy_beta[ix] / oy_se[ix])),
  eaf.outcome = ex_t$eaf[ix], effect_allele.outcome = ex_t$ea[ix], other_allele.outcome = ex_t$nea[ix], outcome = "Simulated heart disease", id.outcome = "sim-out-1", proxy.outcome = proxy,
  target_snp.outcome = if (proxy) mk_snp(ix) else NA_character_, proxy_snp.outcome = if (proxy) sprintf("rs%d", 7000 + ix) else NA_character_, stringsAsFactors = FALSE)
CATM <- data.frame(id = c("sim-exp-1", "sim-out-1", "sim-out-2", "sim-eas-1", "ukb-sim-1", "ukb-sim-2"), trait = c("Simulated body mass", "Simulated heart disease", "Simulated height", "Simulated East Asian trait", "UKB simulated A", "UKB simulated B"),
  population = c("European", "European", "European", "East Asian", "European", "European"), sex = "Males and Females", build = "HG19/GRCh37", author = c("Sim", "Sim", "Sim", "Sim", "Neale", "Neale"), consortium = c("GIANT", "", "", "BBJ", "", ""),
  category = "Risk factor", subcategory = "Anthropometric", unit = "SD", note = "", pmid = "1", ncase = c(NA, 30000, NA, NA, NA, 5000), ncontrol = c(NA, 200000, NA, NA, NA, 300000), sample_size = c(300000, 230000, 250000, 150000, 360000, 305000),
  nsnp = 1e7, year = 2020, mr = c(1, 1, 1, 0, 1, 1), stringsAsFactors = FALSE)
mock <- new.env(); mock$calls <- list(); mock$fail <- NULL
chk <- function(jwt) { if (!is.null(mock$fail)) stop(mock$fail); if (!identical(jwt, tok_ok)) stop("Status code from OpenGWAS API: 401 Message: must provide a token") }
BE <- list(
  info = function(jwt, id = NULL) { chk(jwt); mock$calls[[length(mock$calls) + 1L]] <- list("info", id); if (is.null(id)) CATM else { r <- CATM[CATM$id == id, ]; if (!nrow(r)) stop("Status code from OpenGWAS API: 404") else r } },
  tophits = function(id, p, r2, kb, pop, jwt) { chk(jwt); mock$calls[[length(mock$calls) + 1L]] <- list("tophits", id, p, r2, kb, pop); ex_raw },
  outcome = function(snps, id, proxies, rsq, maf, jwt) { chk(jwt); mock$calls[[length(mock$calls) + 1L]] <- list("outcome", id, length(snps), proxies, rsq, maf)
    direct <- setdiff(1:nx, c(4, 11, 15)); d <- mk_out(direct); if (isTRUE(proxies)) d <- rbind(d, mk_out(c(4, 11), TRUE)); d })
options(mrgui.ogw_backend = BE)
st <- ogw_catalogue_update(tok_ok)
expect("opengwas catalogue: saved on this computer with a date, and read back", st$have && st$n == 6L && nrow(ogw_catalogue_read()) == 6L && !is.null(ogw_catalogue_status()$fetched))
expect("opengwas catalogue: a wrong token gives the token message, not a technical one", grepl("did not accept the access token", err_of(ogw_catalogue_update(fake_jwt(as.numeric(Sys.time()) + 86400)))))
cat12 <- ogw_catalogue_read()
sr <- function(...) ogw_search(cat12, ...)$rows$id
expect("opengwas search: every word must match, the exact trait first, and an empty query returns nothing", identical(sr("heart"), "sim-out-1") && identical(sr("body mass"), "sim-exp-1") && length(sr("")) == 0L && length(sr("heart body")) == 0L && identical(sr("SIM-OUT-2"), "sim-out-2"))
expect("opengwas search: the MR-suitable filter, ancestry, size and binary/continuous filters work", !("sim-eas-1" %in% sr("simulated")) && identical(sr("simulated", population = "EAS", mr_only = FALSE), "sim-eas-1") &&
       !("sim-out-1" %in% sr("simulated", min_n = 240000)) && "sim-out-2" %in% sr("simulated", min_n = 240000) && setequal(sr("simulated", kind = "binary"), c("sim-out-1", "ukb-sim-2")) && "sim-exp-1" %in% sr("simulated", kind = "continuous"))
expect("opengwas search: larger studies come first among equal matches", identical(sr("simulated", kind = "continuous")[1], "sim-exp-1"))
pl <- ogw_search_payload(ogw_search(cat12, "heart")); expect("opengwas search: the payload gives the screen what it shows (binary kind, ancestry code, no empty fields)", pl$rows[[1]]$kind == "binary" && pl$rows[[1]]$pop == "EUR" && !any(vapply(pl$rows[[1]], function(v) length(v) == 1L && is.na(v), TRUE)))

lib0 <- length(library_list())
rx <- ogw_get_exposure("sim-exp-1", tok_ok)
th <- Filter(function(x) x[[1]] == "tophits", mock$calls)[[1]]
expect("opengwas exposure: clumping settings reach the server call (P 5e-8, r2 0.001, 10000 kb, European panel taken from the dataset)", th[[3]] == 5e-8 && th[[4]] == 0.001 && th[[5]] == 10000 && th[[6]] == "EUR")
expect("opengwas exposure: the 18 raw rows are cleaned to 15 (an indel, a repeat and a missing standard error dropped) and each drop is explained", rx$rows_in_file == 18L && rx$rows_kept == 15L && length(rx$warnings) >= 3L)
expect("opengwas exposure: saved as independent, with its source, settings, package versions and the catalogue facts", isTRUE(rx$usable_as_exposure) && !isTRUE(rx$usable_as_outcome) && identical(rx$facts$independent, "yes") && rx$facts$ancestry == "EUR" && identical(rx$source$id, "sim-exp-1") && rx$source$params$kb == 10000 && !is.null(rx$source$packages$TwoSampleMR) && rx$facts$kind == "continuous")
expect("opengwas exposure: it appears in the Data library and in the trait lists as OpenGWAS, as an exposure only", length(library_list()) == lib0 + 1L && any(vapply(catalogue_all()$exposures, function(r) identical(r$src, "OpenGWAS"), TRUE)) && !any(vapply(catalogue_all()$outcomes, function(r) identical(r$src, "OpenGWAS"), TRUE)))
n_calls <- length(mock$calls); rx2 <- ogw_get_exposure("sim-exp-1", tok_ok)
expect("opengwas exposure: asking again with the same settings reuses the saved copy and makes no request", isTRUE(rx2$reused) && length(mock$calls) == n_calls)
expect("opengwas exposure: the same name with other settings is refused instead of overwritten", grepl("already in the Data library", err_of(ogw_get_exposure("sim-exp-1", tok_ok, p = 5e-6))))
expect("opengwas exposure: out-of-range settings are refused with a sentence", grepl("between 1e-30 and 0.001", err_of(ogw_get_exposure("sim-exp-1", tok_ok, p = 0.5, name = "x"))) && grepl("r²", err_of(ogw_get_exposure("sim-exp-1", tok_ok, r2 = 0.9, name = "x"))) && grepl("window", err_of(ogw_get_exposure("sim-exp-1", tok_ok, kb = 1, name = "x"))))
expect("opengwas exposure: a dataset with no clumping panel (mixed ancestry) needs the user to choose one", { CATM$population[1] <- "Mixed"; mock$cat <- NULL; unlink(ogw_cache_file()); ogw_catalogue_update(tok_ok); r <- grepl("no clumping reference panel", err_of(ogw_get_exposure("sim-exp-1", tok_ok, name = "Mixed test"))); CATM$population[1] <- "European"; unlink(ogw_cache_file()); ogw_catalogue_update(tok_ok); r })

ro <- ogw_get_outcome("sim-out-1", rx$trait, tok_ok)
oc <- Filter(function(x) x[[1]] == "outcome", mock$calls)[[1]]
expect("opengwas outcome: only the 15 instruments are looked up, with proxies at r2 0.8 and the MAF limit 0.3", oc[[3]] == 15L && isTRUE(oc[[4]]) && oc[[5]] == 0.8 && oc[[6]] == 0.3)
expect("opengwas outcome: 14 of 15 found (12 direct, 2 through proxies), 1 lost, and the recipe and warnings say so", ro$lookup$found == 14L && ro$lookup$direct == 12L && ro$lookup$by_proxy == 2L && ro$lookup$missing == 1L && any(grepl("1 of 15 instruments are not in this outcome", unlist(ro$warnings))) && any(grepl("replaced by a proxy", unlist(ro$warnings))))
expect("opengwas outcome: binary outcome, usable as outcome only, tied to its exposure", ro$facts$kind == "binary" && isTRUE(ro$usable_as_outcome) && !isTRUE(ro$usable_as_exposure) && identical(ro$source$for_exposure, rx$trait))
od <- library_dataset(ro$trait)$data
expect("opengwas outcome: proxy variants keep the instrument's ID and record the proxy used", all(od$snp[od$proxy_snp != ""] %in% mk_snp(c(4, 11))) && setequal(od$proxy_snp[od$proxy_snp != ""], sprintf("rs%d", c(7004, 7011))))
ro_np <- ogw_get_outcome("sim-out-1", rx$trait, tok_ok, proxies = FALSE, name = "No proxies")
expect("opengwas outcome: with proxies off only the 12 direct hits are kept, and the warning says proxies were off", ro_np$lookup$found == 12L && any(grepl("proxies were off", unlist(ro_np$warnings))))

pp12 <- pair_input(paste0("lib:", rx$trait), paste0("lib:", ro$trait))
expect("opengwas pair: built by the same harmonisation as imported files; the ledger counts the proxies and what was lost", pp12$ledger$instruments == 15L && pp12$ledger$by_proxy == 2L && pp12$ledger$missing_in_outcome == 1L && pp12$ledger$kept == 14L && !isTRUE(pp12$ukb_both) && pp12$anc_x == "EUR" && pp12$anc_y == "EUR")
pf12 <- pair_facts(paste0("lib:", rx$trait), paste0("lib:", ro$trait))
expect("opengwas pair: the facts say which sides came from OpenGWAS, so the Check screen can word ancestry correctly", isTRUE(pf12$og_x) && isTRUE(pf12$og_y) && !isTRUE(pair_facts("expA", "outA")$og_x))
rhome_orig <- Sys.getenv("R_HOME"); Sys.setenv(R_HOME = "/keep/this/value"); .with_user_env(expect("programs started for the user do not inherit the bundled R settings", !nzchar(Sys.getenv("R_HOME"))))
expect("programs started for the user: the engine's own settings are put back afterwards", identical(Sys.getenv("R_HOME"), "/keep/this/value")); Sys.setenv(R_HOME = rhome_orig)
tmp_r <- tempfile(fileext = ".R"); writeLines("x", tmp_r)
expect("text editor: only scripts inside the studies folder can be opened", grepl("Only scripts", err_of(open_in_text_editor(tmp_r))))
h12 <- studio_new("OpenGWAS pair"); r12 <- studio_run(h12, pp12$exposure_label, pp12$outcome_label, pp12$input, c("ivw_mre", "egger"), idx = 1L)
expect("opengwas pair: the analysis runs from downloaded data, replay identical, IVW estimate near the planted 0.25", all(vapply(r12$results, function(x) x$status, "") == "ok") && identical(r12$replay$status, "identical") && abs(pair_results(h12$dir, r12$pair_folder)$estimate[1] - 0.25) < 0.08)
ui_state_set(h12$con, "pairs", list(list(ex = paste0("lib:", rx$trait), out = paste0("lib:", ro$trait), exName = pp12$exposure_label, outName = pp12$outcome_label, folder = r12$pair_folder, methods = as.list(c("ivw_mre", "egger")), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = pp12$ledger, anc_x = "EUR", anc_y = "EUR")))
pk12 <- build_package(h12, c("report", "refs"), 140); htm12 <- paste(readLines(file.path(pk12$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("opengwas report: the methods name the OpenGWAS datasets, the date, the clumping settings, the proxy rule and the package versions", grepl("obtained from the OpenGWAS database (dataset sim-exp-1", htm12, fixed = TRUE) && grepl("from the stored, pre-computed clumps of OpenGWAS (default settings, preclumped = 1)", htm12, fixed = TRUE) && grepl("within 10,000 kb, 1000 Genomes EUR", htm12, fixed = TRUE) && grepl("of at least 0.8", htm12, fixed = TRUE) && grepl("found only through proxy variants", htm12, fixed = TRUE) && !grepl("were were", htm12, fixed = TRUE) && grepl("sample size 300,000", htm12, fixed = TRUE))
expect("opengwas report: no wrong claim that the analyst declared independence, the limitations say the panel is a reference sample, data availability names the public source", !grepl("did not perform clumping", htm12, fixed = TRUE) && !grepl("rests on the analyst", htm12, fixed = TRUE) && grepl("reference sample", htm12, fixed = TRUE) && grepl("public OpenGWAS database (api.opengwas.io), dataset IDs: sim-exp-1", htm12, fixed = TRUE) && grepl("sim-out-1", htm12, fixed = TRUE) && !grepl("Add the public source", htm12, fixed = TRUE))
DBI::dbDisconnect(h12$con)
rx_b <- ogw_get_exposure("sim-exp-1", tok_ok, p = 5e-6, name = "Second exposure")
expect("opengwas pair: an outcome fetched for one exposure cannot be paired with another (its other variants were never downloaded)", grepl("downloaded only for the instruments of", err_of(pair_input(paste0("lib:", rx_b$trait), paste0("lib:", ro$trait)))))
expect("opengwas pair: a downloaded exposure with an imported outcome is allowed (same kind of ID)", is.null(tryCatch({ pair_input(paste0("lib:", rx$trait), paste0("lib:", ro_np$trait)); NULL }, error = function(e) e)))

ux <- ogw_get_exposure("ukb-sim-1", tok_ok, name = "UKB exposure"); uo <- ogw_get_outcome("ukb-sim-2", ux$trait, tok_ok, name = "UKB outcome")
expect("opengwas: UK Biobank datasets are flagged from the catalogue, and a pair of two of them carries the overlap flag", isTRUE(ux$facts$ukb) && isTRUE(uo$facts$ukb) && isTRUE(pair_input(paste0("lib:", ux$trait), paste0("lib:", uo$trait))$ukb_both) && isTRUE(pair_facts(paste0("lib:", ux$trait), paste0("lib:", uo$trait))$ukb_both))

BE$ping <- function() list(ok = TRUE, detail = "api.opengwas.io answered (HTTP 200)"); options(mrgui.ogw_backend = BE)
dg <- run_diagnostics(net = TRUE); dg0 <- run_diagnostics()
expect("diagnostics: the internet check is only run when asked for, and says what the server answered", length(dg) == length(dg0) + 1L && dg[[length(dg)]]$name == "Internet access to OpenGWAS" && isTRUE(dg[[length(dg)]]$ok) && grepl("HTTP 200", dg[[length(dg)]]$detail))
expect("opengwas errors: a network failure keeps the technical detail (without any token) so a failure can be diagnosed", grepl("Technical detail: .*certificate", ogw_friendly(simpleError(paste("SSL certificate problem: unable to get local issuer certificate", tok_ok)))) && !grepl("eyJ", ogw_friendly(simpleError(paste("SSL certificate problem", tok_ok)))))
fe <- function(m) ogw_friendly(simpleError(m))
expect("opengwas errors: each technical failure becomes a sentence that says what to do", grepl("access token", fe("Status code from OpenGWAS API: 401")) && grepl("limited how many requests", fe("Status code from OpenGWAS API: 429")) && grepl("Cannot reach OpenGWAS", fe("Could not resolve host: api.opengwas.io")) && grepl("server error", fe("Status code from OpenGWAS API: 503")) && grepl("does not know this dataset", fe("Status code from OpenGWAS API: 404")) && grepl("Cannot reach", fe("Timeout was reached")))
mock$fail <- paste("curl error", tok_ok, "Could not resolve host")
msg <- err_of(ogw_get_exposure("sim-eas-1", tok_ok, name = "Net down", pop = "EAS"))
expect("opengwas errors: a network failure during a download says so, saves nothing and never repeats the token", grepl("Cannot reach OpenGWAS", msg) && !grepl("eyJ", msg) && !dir.exists(file.path(library_dir(), "Net down")))
mock$fail <- paste("unexpected thing", tok_ok); msg <- err_of(ogw_catalogue_update(tok_ok))
expect("opengwas errors: an unknown error is shown shortened and with any token removed", grepl("unexpected answer", msg) && !grepl("eyJ", msg) && nchar(msg) < 400)
mock$fail <- NULL
expect("opengwas download: offline mode and a missing token are refused before anything is requested", { n0 <- length(mock$calls); r <- grepl("Offline mode", err_of(ogw_get_exposure("sim-exp-1", tok_ok, offline = TRUE, name = "z"))) && grepl("Paste your OpenGWAS", err_of(ogw_get_outcome("sim-out-1", rx$trait, NULL, name = "z"))); r && length(mock$calls) == n0 })
expect("opengwas download: asking for an unknown dataset ID gives the not-found sentence", grepl("does not know", err_of(ogw_get_exposure("nope-1", tok_ok, name = "q"))) && grepl("not a set of independent instruments|not in the Data library", err_of(ogw_get_outcome("sim-out-1", "no such dataset", tok_ok, name = "q"))))

shiny::testServer(mr_server(args_dir), {
  session$setInputs(act = list(t = "ogw_token_save", token = "not a token", nonce = 1)); expect("server: a text that is not a token is not stored", is.null(ogw_token_get()))
  session$setInputs(act = list(t = "ogw_token_save", token = paste0(" ", tok_ok, " "), nonce = 2)); expect("server: a pasted token is kept in memory only", identical(ogw_token_get(), tok_ok))
  session$setInputs(act = list(t = "ogw_search", q = "heart", pop = "", min_n = 0, mr_only = TRUE, kind = "any", nonce = 3)); expect("server: searching the saved catalogue does not crash", TRUE)
  session$setInputs(act = list(t = "ogw_dl_exposure", id = "sim-eas-1", p = 5e-8, r2 = 0.001, kb = 10000, pop = "EAS", name = "EAS exposure", nonce = 4)); expect("server: a download through the screen's action saves the dataset", dir.exists(file.path(library_dir(), "EAS exposure")))
  session$setInputs(act = list(t = "ogw_token_clear", nonce = 5)); expect("server: removing the token clears it", is.null(ogw_token_get()))
})
scan_for_token <- function(dirs) { fs <- unlist(lapply(dirs, function(d) list.files(d, recursive = TRUE, full.names = TRUE, all.files = TRUE))); fs <- fs[file.info(fs)$size > 0 & !file.info(fs)$isdir]; any(vapply(fs, function(f) { x <- readBin(f, "raw", file.size(f)); grepl(substr(tok_ok, 1, 40), rawToChar(replace(x, x == as.raw(0), as.raw(32))), fixed = TRUE) }, TRUE)) }
expect("opengwas: the token is in no file at all (studies folder, Data library, recipes, scripts, settings, catalogue cache)", !scan_for_token(c(root, settings_dir())))
ca_tmp <- tempfile(fileext = ".pem"); writeLines(c("-----BEGIN CERTIFICATE-----", rep(strrep("A", 70), 20), "-----END CERTIFICATE-----"), ca_tmp)
ca_old <- Sys.getenv(c("CURL_CA_BUNDLE", "SSL_CERT_FILE"), unset = NA)
ca_path <- ogw_ca_setup(candidates = c("/no/such/runner/ssl/cacert.pem", ca_tmp), force = TRUE)
expect("certificates: a missing build-machine path is skipped and the first file that exists is used", identical(ca_path, normalizePath(ca_tmp, winslash = "/")) && identical(Sys.getenv("CURL_CA_BUNDLE"), ca_path) && identical(Sys.getenv("SSL_CERT_FILE"), ca_path))
expect("certificates: with no file found nothing is changed and no error is raised", identical(ogw_ca_setup(candidates = "/no/such/file.pem", force = TRUE), "") && identical(Sys.getenv("CURL_CA_BUNDLE"), ca_path))
pem_good <- tempfile(fileext = ".pem"); writeLines(c("-----BEGIN CERTIFICATE-----", strrep("A", 70), rep(strrep("B", 70), 20), "-----END CERTIFICATE-----"), pem_good)
pem_bad <- tempfile(fileext = ".pem"); writeLines(rep("not a certificate file at all, only text", 40), pem_bad)
expect("certificates: a file without certificates in it, a folder and an empty name are not accepted", .ogw_ca_valid(pem_good) && !.ogw_ca_valid(pem_bad) && !.ogw_ca_valid(tempdir()) && !.ogw_ca_valid("") && !.ogw_ca_valid("/no/such.pem"))
vend <- file.path(args_dir, "certs", "cacert.pem")
expect("certificates: the copy shipped with the app exists and holds the usual 100+ authorities", file.exists(vend) && .ogw_ca_valid(vend) && sum(grepl("BEGIN CERTIFICATE", readLines(vend, warn = FALSE))) > 100)
options(mrgui.engine_dir = args_dir); expect("certificates: the shipped copy is among the places searched, after the system's", tail(ogw_ca_candidates(), 1) == file.path(args_dir, "certs", "cacert.pem"))
ogw_ca_setup(candidates = c("/no/such.pem", pem_bad, pem_good, vend), force = TRUE)
expect("certificates: only valid files are kept, in order", identical(.ogw_ca_state$list, normalizePath(c(pem_good, vend), winslash = "/")))
n_try <- 0L; r <- .ogw_do({ n_try <<- n_try + 1L; if (n_try == 1L) stop("Problem with the SSL CA cert (path? access rights?)") else "worked" })
expect("certificates: after a certificate error the next file is tried and the request is repeated", identical(r, "worked") && n_try == 2L && identical(.ogw_ca_state$path, normalizePath(vend, winslash = "/")))
n_try <- 0L; msg <- tryCatch(.ogw_do({ n_try <<- n_try + 1L; stop("Problem with the SSL CA cert") }), error = function(e) conditionMessage(e))
expect("certificates: when every file fails the user gets the normal sentence with technical detail, after exactly one try per file", n_try == 1L && grepl("Technical detail", msg))
n_try <- 0L; msg <- tryCatch(.ogw_do({ n_try <<- n_try + 1L; stop("HTTP 401 unauthorised") }), error = function(e) conditionMessage(e))
expect("certificates: an error that is not about certificates is never repeated", n_try == 1L)
ogw_ca_setup(candidates = c(pem_good, vend), force = TRUE)
for (nm in c("CURL_CA_BUNDLE", "SSL_CERT_FILE")) if (is.na(ca_old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(ca_old[[nm]]), nm))
ogw_ca_setup(force = TRUE)
options(mrgui.ogw_backend = NULL); ogw_token_set(NULL)
# ---- Round 19: clumping on this computer (PLINK 1.9 + LD reference panel) ----
options(mrgui.engine_dir = args_dir)
fix_dir <- file.path(args_dir, "test_fixtures", "ld_panel")   # synthetic panel: blocks A (rs91..) and B (rs92..), independent rs93.. and rs9400001
Sys.setenv(MRGUI_LD_DIR = fix_dir)
expect("clump: PLINK is bundled and runs", !is.na(plink_bin()) && grepl("PLINK v1.9", plink_version()))
expect("clump: panel status sees the installed ancestry only", isTRUE(ld_panel_status()$installed$EUR) && !isTRUE(ld_panel_status()$installed$EAS))
cand <- data.frame(snp = c(sprintf("rs91%05d", 1:4), sprintf("rs92%05d", 1:3), sprintf("rs93%05d", 1:5), "rs9400001", "rs9999999"),
                   p = c(1e-9, 1e-12, 1e-10, 1e-8, 1e-15, 1e-9, 1e-10, rep(1e-9, 5), 1e-9, 1e-20), stringsAsFactors = FALSE)
cl <- ld_clump_here(cand$snp, cand$p, 0.001, 10000, "EUR")
st_of <- function(id) cl$status[cl$variant_id == id]
expect("clump: the strongest variant of each LD block is kept, the rest are removed with their index variant",
       identical(sort(cl$variant_id[cl$status == "index"]), sort(c("rs9100002", "rs9200001", sprintf("rs93%05d", 1:5), "rs9400001"))) &&
       all(cl$index_variant[cl$variant_id %in% c("rs9100001", "rs9100003", "rs9100004")] == "rs9100002") && all(cl$index_variant[cl$variant_id %in% c("rs9200002", "rs9200003")] == "rs9200001"))
expect("clump: a variant missing from the panel is reported as such and not kept", identical(st_of("rs9999999"), "not_in_panel"))
ie <- ieugwasr::ld_clump_local(data.frame(rsid = cand$snp, pval = cand$p), clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = file.path(fix_dir, "EUR"), plink_bin = plink_bin())
expect("clump: keeps exactly the variants ieugwasr::ld_clump_local keeps", identical(sort(ie$rsid), sort(cl$variant_id[cl$status == "index"])))
expect("clump: a looser r2 (0.99) removes nothing in these blocks, 0.95 still removes the tightest pair", all(ld_clump_here(cand$snp, cand$p, 0.99, 10000, "EUR")$status[1:13] == "index") && sum(ld_clump_here(cand$snp, cand$p, 0.95, 10000, "EUR")$status == "clumped") == 1L)
expect("clump: a 0.05 kb window cannot reach variants 100 bp apart", all(ld_clump_here(cand$snp, cand$p, 0.001, 0.05, "EUR")$status[1:13] == "index"))
expect("clump: settings are recorded", identical(attr(cl, "settings")$pop, "EUR") && attr(cl, "settings")$r2 == 0.001 && grepl("PLINK", attr(cl, "settings")$plink))
expect("clump: clear messages for a missing ancestry, positions only and an unknown ancestry",
       grepl("not installed", err_of(ld_clump_here(cand$snp, cand$p, pop = "EAS"))) && grepl("rsIDs", err_of(ld_clump_here(c("1:1000", "1:2000"), c(1e-9, 1e-9)))) &&
       grepl("1000 Genomes ancestries", err_of(ld_clump_here(cand$snp, cand$p, pop = "Mixed"))))
# whole path: an imported genome-wide exposure file set to be clumped here, and an outcome file
td19 <- tempfile("clumpfiles-"); dir.create(td19)
set.seed(19); allsnp <- c(cand$snp, sprintf("rs8%06d", 1:30)); k19 <- length(allsnp)
pex <- c(cand$p, runif(30, 1e-4, 0.9)); bx19 <- runif(k19, 0.03, 0.08) * sample(c(-1, 1), k19, TRUE); sx19 <- abs(bx19) / stats::qnorm(pex / 2, lower.tail = FALSE)
utils::write.csv(data.frame(SNP = allsnp, effect_allele = "A", other_allele = "C", eaf = 0.35, beta = signif(bx19, 8), se = signif(sx19, 8), pval = pex, n = 300000), file.path(td19, "exp_full.csv"), row.names = FALSE)
utils::write.csv(data.frame(SNP = allsnp, effect_allele = "A", other_allele = "C", eaf = 0.35, beta = signif(0.25 * bx19 + rnorm(k19, 0, 0.005), 8), se = 0.006, pval = 0.01, n = 200000), file.path(td19, "out_full.csv"), row.names = FALSE)
imp19 <- function(f, facts) { hd <- import_read_head(file.path(td19, f)); sg <- import_suggest(hd$columns, hd$rows); mp <- lapply(sg, function(x) x$column); mp <- mp[!vapply(mp, is.null, TRUE)]; import_save(file.path(td19, f), f, hd$sep, hd$columns, mp, facts) }
b19 <- list(frq = "eaf", scale = "beta", ea_ok = "yes", build = "37", ancestry = "EUR", sex = "both", ukb = FALSE, n = "", kind = "continuous", scope = "genomewide")
hd19 <- import_read_head(file.path(td19, "exp_full.csv")); mp19 <- lapply(import_suggest(hd19$columns, hd19$rows), function(x) x$column); mp19 <- mp19[!vapply(mp19, is.null, TRUE)]
expect("clump import: an invalid r2 is refused before saving", !import_validate(file.path(td19, "exp_full.csv"), hd19$sep, hd19$columns, mp19, c(b19, list(role = "exposure", independent = "clump", clump_r2 = 2, clump_kb = 10000, clump_pop = "EUR")))$ok)
rx19 <- imp19("exp_full.csv", c(b19, list(name = "Full exposure", role = "exposure", independent = "clump", clump_r2 = 0.001, clump_kb = 10000, clump_pop = "EUR")))
imp19("out_full.csv", c(b19, list(name = "Full outcome", role = "outcome", independent = "no")))
expect("clump import: a full file set to be clumped is usable as an exposure", isTRUE(rx19$usable_as_exposure) && identical(rx19$facts$independent, "clump"))
p19 <- pair_input("lib:Full exposure", "lib:Full outcome", 5e-8)
expect("clump pair: 14 candidates below P, 5 removed for LD, 1 not in the panel, 8 instruments",
       p19$ledger$below_threshold == 14 && p19$ledger$clump$clumped == 5 && p19$ledger$clump$not_in_panel == 1 && p19$ledger$instruments == 8 && nrow(p19$input$data) == 8)
Sys.setenv(MRGUI_PLINK = file.path(tempdir(), "no-plink-here"))
p19b <- tryCatch(pair_input("lib:Full exposure", "lib:Full outcome", 5e-8), error = function(e) e)
Sys.unsetenv("MRGUI_PLINK")
expect("clump pair: the same inputs are read from the cache, without running PLINK again", !inherits(p19b, "error") && identical(p19b$selection$status, p19$selection$status))
h19 <- studio_new("Clumped here")
r19 <- studio_run(h19, p19$exposure_label, p19$outcome_label, p19$input, c("ivw_mre", "egger"), idx = 1L)
save_instrument_selection(h19$dir, r19$pair_folder, p19$selection)
sel_rel <- file.path("03_pairs", r19$pair_folder, "01_data", "instrument_selection.csv")
sel19 <- utils::read.csv(file.path(h19$dir, sel_rel), stringsAsFactors = FALSE)
expect("clump study: every candidate is saved with its status and index variant, and listed in the manifest",
       nrow(sel19) == 14 && all(c("variant_id", "status", "index_variant", "clump_r2", "clump_kb", "reference") %in% names(sel19)) && sel_rel %in% manifest_read(h19$dir)$relative_path && identical(r19$replay$status, "identical"))
ui_state_set(h19$con, "pairs", list(list(ex = "lib:Full exposure", out = "lib:Full outcome", exName = p19$exposure_label, outName = p19$outcome_label, folder = r19$pair_folder, methods = list("ivw_mre", "egger"), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = p19$ledger, anc_x = "EUR", anc_y = "EUR")))
pk19 <- build_package(h19, c("report", "text", "refs"), 140)
htm19 <- paste(readLines(file.path(pk19$dir, "report.html"), warn = FALSE), collapse = "\n"); bib19 <- paste(readLines(file.path(pk19$dir, "references.bib"), warn = FALSE), collapse = "\n")
expect("clump report: methods say clumped here with PLINK and the 1000 Genomes panel, with counts, and cite both",
       grepl("clumped with PLINK", htm19, fixed = TRUE) && grepl("5 variants were removed as correlated", htm19, fixed = TRUE) && grepl("Chang2015", bib19, fixed = TRUE) && grepl("1000Genomes2015", bib19, fixed = TRUE))
# ---- Round 33: allele-frequency check against the panel (synthetic panel: A allele at about 0.30, 0.40 and 0.35; files report EAF 0.35) ----
e19 <- p19$ledger$eaf_check
expect("Round 33: both datasets are checked against the EUR panel; with EAF 0.35 none differs by more than 0.2, and the counts add up",
       isTRUE(e19$x$run) && isTRUE(e19$y$run) && identical(e19$x$pop, "EUR") && e19$x$judged >= 1L && e19$x$flagged == 0L && e19$y$flagged == 0L &&
       e19$x$judged + e19$x$palindromic + e19$x$not_in_panel + e19$x$near_half == e19$x$with_eaf && e19$x$with_eaf == p19$ledger$kept && e19$x$max_difference < 0.2)
d19 <- p19$input$data; ds19 <- file.path(library_dir(), "Full exposure")
fl19 <- eaf_panel_check(d19, rep(0.85, nrow(d19)), ds19, "EUR", "Flipped")
expect("Round 33: an EAF given for the other allele (0.85 instead of about 0.3) is flagged for every judged variant and listed by ID; nothing is removed",
       isTRUE(fl19$run) && fl19$judged >= 1L && fl19$flagged == fl19$judged && setequal(fl19$flagged_ids, fl19$rows$variant_id[fl19$rows$judged]) && nrow(d19) == p19$ledger$kept)
expect("Round 33: variants with a panel minor allele frequency of 0.4 or more are not judged, palindromic variants are not judged",
       all(abs(pmin(fl19$rows$panel_eaf, 1 - fl19$rows$panel_eaf)[fl19$rows$judged]) < 0.4) && !any(fl19$rows$judged & ((fl19$rows$ea %in% c("A", "T") & fl19$rows$oa %in% c("A", "T")) | (fl19$rows$ea %in% c("C", "G") & fl19$rows$oa %in% c("C", "G")))))
expect("Round 33: the check is skipped with a reason for a mixed ancestry, a panel that is not installed, a dataset without EAF and a missing PLINK",
       grepl("not one of the five", eaf_panel_check(d19, d19$eaf_x, ds19, "Mixed", "T")$why) && grepl("East Asian (EAS) reference panel is not installed", eaf_panel_check(d19, d19$eaf_x, ds19, "EAS", "T")$why, fixed = TRUE) &&
       grepl("reports no effect-allele frequency", eaf_panel_check(d19, rep(NA_real_, nrow(d19)), ds19, "EUR", "T")$why) &&
       grepl("PLINK is not available", eaf_panel_check(d19, d19$eaf_x, ds19, "EUR", "T", bin = NA_character_)$why))
expect("Round 33: panel frequencies are cached in the dataset folder and read again without PLINK",
       length(list.files(file.path(ds19, "panel_freq"), pattern = "^EUR_.*\\.csv$")) >= 1L && nrow(ld_panel_freq_cached(ds19, d19$variant_id, "EUR", bin = "/no/plink")) >= 1L)
save_eaf_check(h19$dir, r19$pair_folder, p19$eaf_rows); ef_rel <- file.path("03_pairs", r19$pair_folder, "01_data", "allele_frequency_check.csv")
expect("Round 33: the check is saved per pair with one row per variant and dataset, and listed in the manifest",
       file.exists(file.path(h19$dir, ef_rel)) && nrow(utils::read.csv(file.path(h19$dir, ef_rel))) == 2L * p19$ledger$kept && ef_rel %in% manifest_read(h19$dir)$relative_path)
expect("Round 33: the Methods report the check of both datasets with the panel and the largest difference",
       grepl("The effect-allele frequencies of Full exposure were compared with the 1000 Genomes EUR reference panel", htm19, fixed = TRUE) && grepl("none differed by more than 0.2", htm19, fixed = TRUE) && grepl("The effect-allele frequencies of Full outcome", htm19, fixed = TRUE))
l33 <- p19$ledger; l33$eaf_check$x <- fl19[setdiff(names(fl19), "rows")]
ui_state_set(h19$con, "pairs", list(list(ex = "lib:Full exposure", out = "lib:Full outcome", exName = p19$exposure_label, outName = p19$outcome_label, folder = r19$pair_folder, methods = list("ivw_mre", "egger"), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = l33, anc_x = "EUR", anc_y = "EUR")))
t33 <- report_text(lapply(ui_state_get(h19$con, "pairs"), function(p) pair_context(h19, p)), "Clumped here")
expect("Round 33: flagged variants are named in the Methods (kept, listed in the file) and the Limitations advise a sensitivity analysis without them",
       grepl(sprintf("%d differed by more than 0.2, which suggests an allele or strand error", fl19$flagged), paste(t33$methods, collapse = " "), fixed = TRUE) &&
       grepl("Flagged variants were kept in the analysis and are listed in allele_frequency_check.csv", paste(t33$methods, collapse = " "), fixed = TRUE) &&
       grepl(sprintf("In one pair, %d variant", fl19$flagged), paste(t33$limitations, collapse = " "), fixed = TRUE) && grepl("a sensitivity analysis without", paste(t33$limitations, collapse = " "), fixed = TRUE))
# ---- Round 34: the analysis without the flagged variants, as a labelled sensitivity analysis ----
rm34 <- utils::head(sort(unique(fl19$flagged_ids)), 2L)   # 2 of the flagged variants, so that enough instruments remain
r34 <- studio_run_without(h19, "lib:Full exposure", "lib:Full outcome", rm34, c("ivw_mre", "egger"), 0.95, 2L, params = run_params(list()), pthr = 5e-8)
ui_state_set(h19$con, "sens", list(list(parent = r19$pair_folder, remove = as.list(rm34), folder = r34$pair_folder, label = r34$label, kept = r34$kept, kind = "eaf")))
t34 <- report_text(lapply(ui_state_get(h19$con, "pairs"), function(p) pair_context(h19, p)), "Clumped here")
e34 <- pair_payload(h19$dir, r34$pair_folder)$estimates; e34 <- Filter(function(e) e$method == "ivw_mre" && e$estimand == "causal_effect", e34)[[1]]
expect("Round 34: the run without the flagged variants leaves them out (and nothing else), and the primary analysis is unchanged",
       r34$kept == p19$ledger$kept - length(rm34) && e34$n_ivs_used == p19$ledger$kept - length(rm34) && identical(r19$replay$status, "identical") && !identical(r34$pair_folder, r19$pair_folder))
expect("Round 34: the Results report the analysis without the flagged variants, naming them, and the Limitations point to it instead of advising it",
       grepl(sprintf("Sensitivity analysis without the %d variant", length(rm34)), paste(t34$results, collapse = " "), fixed = TRUE) && grepl(rm34[1], paste(t34$results, collapse = " "), fixed = TRUE) &&
       grepl("compared with", paste(t34$results, collapse = " "), fixed = TRUE) && grepl("the Results report a sensitivity analysis without", paste(t34$limitations, collapse = " "), fixed = TRUE) &&
       !grepl("is advisable", paste(grep("allele frequency|effect-allele", t34$limitations, value = TRUE), collapse = " "), fixed = TRUE))
ui_state_set(h19$con, "sens", list())
# the download job, against a local copy of the fixture packed like the real archive
dl_src <- tempfile("panel-src-"); dir.create(dl_src); file.copy(file.path(fix_dir, c("EUR.bed", "EUR.bim", "EUR.fam")), dl_src)
tgz <- file.path(tempdir(), "fake-1kg.tgz"); old_wd <- setwd(dl_src); utils::tar(tgz, files = c("EUR.bed", "EUR.bim", "EUR.fam"), compression = "gzip", tar = "internal"); setwd(old_wd)
dl_dir <- tempfile("ldref-")
ld_panel_download_start("EUR", dir = dl_dir, url = paste0("file://", normalizePath(tgz, winslash = "/")), bytes = file.size(tgz), sha256 = digest::digest(file = tgz, algo = "sha256"))
for (i in 1:60) { js <- .ld_job_read(dl_dir); if (!is.null(js) && js$state %in% c("done", "error")) break; Sys.sleep(1) }
expect("panel download: the job downloads, checks the checksum, unpacks the chosen ancestry and removes the archive",
       identical(.ld_job_read(dl_dir)$state, "done") && isTRUE(ld_panel_status(dl_dir)$installed$EUR) && !file.exists(file.path(dl_dir, "1kg.v3.tgz.part")) && file.exists(file.path(dl_dir, "SOURCE.txt")))
bad_dir <- tempfile("ldref-bad-")
ld_panel_download_start("EUR", dir = bad_dir, url = paste0("file://", normalizePath(tgz, winslash = "/")), bytes = file.size(tgz), sha256 = strrep("0", 64))
for (i in 1:60) { js <- .ld_job_read(bad_dir); if (!is.null(js) && js$state %in% c("done", "error")) break; Sys.sleep(1) }
expect("panel download: a file with the wrong checksum is refused and deleted", identical(.ld_job_read(bad_dir)$state, "error") && grepl("SHA-256", .ld_job_read(bad_dir)$msg) && !isTRUE(ld_panel_status(bad_dir)$installed$EUR))
ld_panel_remove("EUR", dl_dir)
expect("panel: an ancestry can be removed", !isTRUE(ld_panel_status(dl_dir)$installed$EUR))
# ---- Round 23: an archive the user downloaded in a browser ----
wait_job <- function(d) { for (i in 1:60) { js <- .ld_job_read(d); if (!is.null(js) && js$state %in% c("done", "error")) break; Sys.sleep(1) }; .ld_job_read(d) }
tb <- file.size(tgz); pend <- "PANEL_SHA256_PENDING"
own <- tempfile("browser-"); dir.create(own); own_tgz <- file.path(own, "1kg.v3.tgz"); file.copy(tgz, own_tgz)
file.copy(tgz, file.path(own, "1kg.v3 (1).tgz.crdownload")); dir.create(file.path(own, "1kg.v3")); writeBin(as.raw(1:10), file.path(own, "1kg.v3.tar"))
expect("own archive: missing path, folder, unfinished download, unpacked .tar and wrong size each get their own sentence; the right file passes",
       grepl("No file was found", ld_archive_problem(file.path(own, "nope.tgz"), tb)) && grepl("folder, not the archive", ld_archive_problem(file.path(own, "1kg.v3"), tb)) &&
       grepl("unfinished browser download", ld_archive_problem(file.path(own, "1kg.v3 (1).tgz.crdownload"), tb)) && grepl("unpacked .tar", ld_archive_problem(file.path(own, "1kg.v3.tar"), tb)) &&
       grepl("incomplete or a different file", ld_archive_problem(own_tgz)) && is.null(ld_archive_problem(paste0('"', own_tgz, '"'), tb)))
fnd <- ld_archive_find(own, tb)
expect("own archive: files named 1kg.v3... are found with a verdict each, only the complete .tgz is usable",
       length(fnd) == 4L && sum(vapply(fnd, `[[`, TRUE, "ok")) == 1L && basename(Filter(function(x) x$ok, fnd)[[1]]$path) == "1kg.v3.tgz" && length(ld_archive_find(file.path(own, "none"), tb)) == 0L)
od <- tempfile("ldref-own-")
expect("own archive: a wrong file is refused before any job starts", grepl("incomplete", err_of(ld_panel_download_start("EUR", dir = od, bytes = tb + 1, sha256 = pend, local_file = own_tgz))) && !file.exists(.ld_job_file(od)))
ld_panel_download_start("EUR", dir = od, bytes = tb, sha256 = pend, local_file = own_tgz); jo <- wait_job(od); src_txt <- readLines(file.path(od, "SOURCE.txt"))
expect("own archive: checked and unpacked without downloading, the user's file is left in place, SHA-256 recorded (not called verified), sizes listed",
       identical(jo$state, "done") && identical(jo$mode, "local") && isTRUE(ld_panel_status(od)$installed$EUR) && file.exists(own_tgz) && !file.exists(file.path(od, "1kg.v3.tgz.part")) &&
       identical(jo$sha256, digest::digest(file = own_tgz, algo = "sha256")) && !isTRUE(jo$pinned) && any(grepl("browser", src_txt)) && any(grepl("recorded only", src_txt)) && any(grepl("^EUR: 13 variants", src_txt)))
expect("own archive: the installed panel clumps exactly as the fixture", { o <- Sys.getenv("MRGUI_LD_DIR"); Sys.setenv(MRGUI_LD_DIR = od); r <- identical(ld_clump_here(cand$snp, cand$p, 0.001, 10000, "EUR")$status, cl$status); Sys.setenv(MRGUI_LD_DIR = o); r })
od2 <- tempfile("ldref-own-bad-"); ld_panel_download_start("EUR", dir = od2, bytes = tb, sha256 = strrep("0", 64), local_file = own_tgz); jo2 <- wait_job(od2)
expect("own archive: when a SHA-256 is pinned and differs, the file is refused but not deleted", identical(jo2$state, "error") && grepl("SHA-256 differs", jo2$msg) && file.exists(own_tgz) && !isTRUE(ld_panel_status(od2)$installed$EUR))
dmg <- tempfile("dmg-"); dir.create(dmg); file.copy(file.path(fix_dir, c("EUR.bim", "EUR.fam")), dmg); bb <- readBin(file.path(fix_dir, "EUR.bed"), "raw", file.size(file.path(fix_dir, "EUR.bed"))); writeBin(bb[-length(bb)], file.path(dmg, "EUR.bed"))
dtgz <- file.path(tempdir(), "damaged-1kg.tgz"); old_wd <- setwd(dmg); utils::tar(dtgz, files = c("EUR.bed", "EUR.bim", "EUR.fam"), compression = "gzip", tar = "internal"); setwd(old_wd)
od3 <- tempfile("ldref-own-dmg-"); ld_panel_download_start("EUR", dir = od3, bytes = file.size(dtgz), sha256 = pend, local_file = dtgz); jo3 <- wait_job(od3)
expect("panel check: a .bed that does not match its .bim and .fam (one byte short) is refused and not installed", identical(jo3$state, "error") && grepl("does not match its .bim and .fam", jo3$msg) && !isTRUE(ld_panel_status(od3)$installed$EUR))
# ---- Round 24: panel files the user already has unpacked ----
fixp <- c(EUR = 300, EAS = 504, AFR = 661, AMR = 347, SAS = 489)     # the synthetic fixture has 300 people
unp <- tempfile("unpacked-"); dir.create(file.path(unp, "1000G_EUR"), recursive = TRUE); file.copy(file.path(fix_dir, c("EUR.bed", "EUR.bim", "EUR.fam")), file.path(unp, "1000G_EUR"))
sub_dir <- file.path(unp, "proj"); dir.create(sub_dir); file.copy(file.path(fix_dir, "EUR.bed"), file.path(sub_dir, "local489_EUR.bed")); file.copy(file.path(fix_dir, "EUR.bim"), file.path(sub_dir, "local489_EUR.bim"))
writeLines(readLines(file.path(fix_dir, "EUR.fam"))[1:250], file.path(sub_dir, "local489_EUR.fam"))
r_ok <- ld_unpacked_problem(file.path(unp, "1000G_EUR"), people = fixp)
expect("unpacked panel: a folder with EUR.bed/.bim/.fam is recognised as EUR with its full number of people", is.null(r_ok$problem) && identical(r_ok$pop, "EUR") && r_ok$people == 300L)
expect("unpacked panel: a subset made from the panel (fewer people), a missing .fam, a wrong ancestry choice and a non-.bed file are refused with a sentence",
       grepl("subset made from the panel", ld_unpacked_problem(file.path(sub_dir, "local489_EUR.bed"), people = fixp)$problem) &&
       { file.rename(file.path(sub_dir, "local489_EUR.fam"), file.path(sub_dir, "x.fam")); m <- ld_unpacked_problem(file.path(sub_dir, "local489_EUR.bed"), people = fixp)$problem; file.rename(file.path(sub_dir, "x.fam"), file.path(sub_dir, "local489_EUR.fam")); grepl(".fam file next to", m) } &&
       grepl("file name says EUR but EAS", ld_unpacked_problem(file.path(unp, "1000G_EUR", "EUR.bed"), "EAS", people = fixp)$problem) &&
       grepl("Give the .bed file", ld_unpacked_problem(file.path(unp, "1000G_EUR", "EUR.bim"), people = fixp)$problem))
cd24 <- tempfile("ldref-copy-")
ld_panel_copy_start(file.path(unp, "1000G_EUR", "EUR.bed"), dir = cd24, people = fixp); jc <- wait_job(cd24); sc <- readLines(file.path(cd24, "SOURCE_EUR.txt"))
expect("unpacked panel: copied in, checked, the originals untouched, SHA-256 of all three files recorded",
       identical(jc$state, "done") && identical(jc$mode, "folder") && isTRUE(ld_panel_status(cd24)$installed$EUR) && all(file.exists(file.path(unp, "1000G_EUR", c("EUR.bed", "EUR.bim", "EUR.fam")))) &&
       identical(unname(tools::md5sum(file.path(cd24, "EUR.bed"))), unname(tools::md5sum(file.path(fix_dir, "EUR.bed")))) && sum(grepl("^sha256 EUR\\.(bed|bim|fam): [0-9a-f]{64}$", sc)) == 3L &&
       identical(jc$sha256, digest::digest(file = file.path(fix_dir, "EUR.bed"), algo = "sha256")) && !dir.exists(file.path(cd24, "unpack")))
expect("unpacked panel: the copied panel clumps exactly as the fixture", { o <- Sys.getenv("MRGUI_LD_DIR"); Sys.setenv(MRGUI_LD_DIR = cd24); r <- identical(ld_clump_here(cand$snp, cand$p, 0.001, 10000, "EUR")$status, cl$status); Sys.setenv(MRGUI_LD_DIR = o); r })
dmg24 <- file.path(unp, "dmg"); dir.create(dmg24); file.copy(file.path(dmg, c("EUR.bed", "EUR.bim", "EUR.fam")), dmg24); cd24b <- tempfile("ldref-copy-dmg-")
ld_panel_copy_start(file.path(dmg24, "EUR.bed"), dir = cd24b, people = fixp); jcb <- wait_job(cd24b)
expect("unpacked panel: a damaged .bed is refused before copying and nothing is installed", identical(jcb$state, "error") && grepl("does not match its .bim and .fam", jcb$msg) && !isTRUE(ld_panel_status(cd24b)$installed$EUR) && !dir.exists(file.path(cd24b, "unpack")))
# ---- Round 29: the panel SHA-256 is pinned ----
expect("panel pin: the archive SHA-256 and the 15 file hashes are pinned (64 hex characters each)",
       grepl("^[0-9a-f]{64}$", LD_PANEL$sha256) && length(LD_FILE_SHA256) == 15L && all(grepl("^[0-9a-f]{64}$", LD_FILE_SHA256)) &&
       setequal(names(LD_FILE_SHA256), as.vector(outer(LD_PANEL$pops, c(".bed", ".bim", ".fam"), paste0))))
fixh <- vapply(c(".bed", ".bim", ".fam"), function(e) digest::digest(file = file.path(fix_dir, paste0("EUR", e)), algo = "sha256"), "")
fixpins <- setNames(fixh, paste0("EUR", c(".bed", ".bim", ".fam")))
cd29 <- tempfile("ldref-copy29-"); ld_panel_copy_start(file.path(unp, "1000G_EUR", "EUR.bed"), dir = cd29, people = fixp, pins = fixpins); j29 <- wait_job(cd29); s29 <- readLines(file.path(cd29, "SOURCE_EUR.txt"))
expect("panel pin: files copied in that equal the pinned files are marked as identical, in the job, SOURCE_EUR.txt and the status",
       identical(j29$state, "done") && isTRUE(j29$pinned) && any(grepl("identical to the files in the pinned MRC IEU 1kg.v3 archive", s29, fixed = TRUE)) &&
       isTRUE(ld_panel_matches_pin("EUR", cd29, files = fixpins)) && !isTRUE(ld_panel_matches_pin("EUR", cd29)))
expect("panel pin: the earlier copy (compared with the real pins) is recorded as differing, not refused, and the status says it is not the pinned panel",
       identical(jc$state, "done") && !isTRUE(jc$pinned) && any(grepl("SHA-256 differs from the files in the pinned MRC IEU 1kg.v3 archive", readLines(file.path(cd24, "SOURCE_EUR.txt")), fixed = TRUE)) &&
       isTRUE(ld_panel_status(cd24)$installed$EUR) && identical(ld_panel_status(cd24)$verified$EUR, FALSE))
ownsha <- digest::digest(file = own_tgz, algo = "sha256")
expect("panel pin: an archive install writes one record per ancestry, and it matches when the archive hash equals the pin",
       file.exists(file.path(od, "SOURCE_EUR.txt")) && isTRUE(ld_panel_matches_pin("EUR", od, pin = ownsha)) && !isTRUE(ld_panel_matches_pin("EUR", od)) && !isTRUE(ld_panel_matches_pin("EAS", od, pin = ownsha)))
old29 <- tempfile("ldref-old-"); dir.create(old29); writeLines(c("source: x", paste("sha256 of the archive:", LD_PANEL$sha256), "ancestries: EUR EAS"), file.path(old29, "SOURCE.txt"))
expect("panel pin: installs from before Round 29 (one SOURCE.txt) are judged by the archive hash and the ancestries it lists",
       isTRUE(ld_panel_matches_pin("EUR", old29)) && isTRUE(ld_panel_matches_pin("EAS", old29)) && !isTRUE(ld_panel_matches_pin("SAS", old29)))
ld_panel_remove("EUR", cd29)
expect("panel pin: removing an ancestry also removes its record", !file.exists(file.path(cd29, "SOURCE_EUR.txt")) && !isTRUE(ld_panel_status(cd29)$installed$EUR))
# ---- Round 22: OpenGWAS exposure clumped on this computer ----
set.seed(22); k22 <- nrow(cand)
th22 <- data.frame(rsid = cand$snp, chr = "1", position = 1e6 * seq_len(k22), ea = "A", nea = "C", eaf = 0.35, beta = runif(k22, 0.03, 0.08) * sample(c(-1, 1), k22, TRUE), se = NA_real_, p = cand$p, n = 300000, stringsAsFactors = FALSE)
th22$se <- abs(th22$beta) / stats::qnorm(th22$p / 2, lower.tail = FALSE)
th22 <- rbind(th22, transform(th22[1, ], rsid = "rs9100099", ea = "AT"))     # an indel: dropped before clumping
mock22 <- new.env(); mock22$calls <- list()
BE22 <- list(info = BE$info,
  tophits = function(id, p, r2, kb, pop, jwt, clump = 1, force = FALSE) { chk(jwt); mock22$calls[[length(mock22$calls) + 1L]] <- list("tophits", id, p, clump, force); if (clump == 0) th22 else th22[th22$rsid %in% c("rs9100002", "rs9200001", "rs9300001"), ] },
  outcome = function(snps, id, proxies, rsq, maf, jwt) { chk(jwt); mock22$calls[[length(mock22$calls) + 1L]] <- list("outcome", id, length(snps)); i <- match(snps, th22$rsid)
    data.frame(SNP = snps, chr = "1", pos = th22$position[i], beta.outcome = 0.3 * th22$beta[i] + rnorm(length(i), 0, 0.003), se.outcome = 0.008, samplesize.outcome = 230000, pval.outcome = 0.01, eaf.outcome = 0.35,
               effect_allele.outcome = "A", other_allele.outcome = "C", outcome = "Simulated heart disease", id.outcome = id, proxy.outcome = FALSE, target_snp.outcome = NA_character_, proxy_snp.outcome = NA_character_, stringsAsFactors = FALSE) })
options(mrgui.ogw_backend = BE22); ogw_token_set(tok_ok); mock$fail <- NULL; unlink(ogw_cache_file()); ogw_catalogue_update(tok_ok)
lx <- ogw_get_exposure("sim-exp-1", tok_ok, where = "local")
expect("opengwas local clump: the server is asked for unclumped hits (clump = 0) with the user's P threshold",
       { c1 <- Filter(function(x) x[[1]] == "tophits", mock22$calls)[[1]]; identical(c1[[4]], 0) && c1[[3]] == 5e-8 })
expect("opengwas local clump: 15 rows, 1 indel dropped, 14 candidates clumped here: 5 removed for LD, 1 not in the panel, 8 instruments (as ld_clump_here)",
       lx$rows_in_file == 15L && lx$clump$candidates == 14L && lx$clump$clumped == 5L && lx$clump$not_in_panel == 1L && lx$rows_kept == 8L && lx$clump$kept == 8L &&
       identical(sort(library_dataset(lx$trait)$data$snp), sort(cl$variant_id[cl$status == "index"])))
lxd <- file.path(library_dir(), lx$trait)
expect("opengwas local clump: saved with its own name, where = local, the candidates and every candidate's status kept in the dataset folder",
       grepl("clumped here$", lx$trait) && identical(lx$source$params$where, "local") && identical(lx$facts$independent, "yes") && isTRUE(lx$usable_as_exposure) &&
       nrow(lib_table_read(lxd, "candidates")) == 14L && setequal(lib_table_read(lxd, "instrument_selection")$status, c("index", "clumped", "not_in_panel")) && grepl("clumped on this computer", lx$source$note))
expect("opengwas local clump: a not-in-panel candidate is explained in the warnings", any(grepl("not in the reference panel", unlist(lx$warnings))))
n22 <- length(mock22$calls); lx2 <- ogw_get_exposure("sim-exp-1", tok_ok, where = "local")
expect("opengwas local clump: the same request reuses the saved copy without asking the server", isTRUE(lx2$reused) && length(mock22$calls) == n22)
sx22 <- ogw_get_exposure("sim-exp-1", tok_ok, name = "Server clumped 22")
expect("opengwas local clump: server clumping still sends clump = 1 and keeps its old settings (no where field)", is.null(sx22$source$params$where) && sx22$source$params$clump == 1 && is.null(sx22$clump))
expect("opengwas local clump: an ancestry whose panel is not installed is refused before anything is downloaded",
       { n0 <- length(mock22$calls); m <- err_of(ogw_get_exposure("sim-eas-1", tok_ok, where = "local", name = "EAS local")); grepl("East Asian \\(EAS\\) LD reference panel is not installed", m) && length(mock22$calls) == n0 })
expect("opengwas local clump: without PLINK it is refused before anything is downloaded",
       { Sys.setenv(MRGUI_PLINK = file.path(tempdir(), "no-plink-here")); n0 <- length(mock22$calls); m <- err_of(ogw_get_exposure("sim-exp-1", tok_ok, where = "local", name = "No plink")); Sys.unsetenv("MRGUI_PLINK"); grepl("PLINK was not found", m) && length(mock22$calls) == n0 })
lo <- ogw_get_outcome("sim-out-1", lx$trait, tok_ok, name = "Outcome for local")
expect("opengwas local clump: the outcome lookup sends only the 8 clumped instruments, not the 14 candidates", tail(Filter(function(x) x[[1]] == "outcome", mock22$calls), 1)[[1]][[3]] == 8L && lo$lookup$requested == 8L)
p22 <- pair_input(paste0("lib:", lx$trait), "lib:Outcome for local", 5e-8)
expect("opengwas local clump pair: the ledger shows the clumping (14 below P, 5 for LD, 1 not in panel, 8 kept) and the selection is passed on",
       p22$ledger$below_threshold == 14L && p22$ledger$clump$clumped == 5L && p22$ledger$clump$not_in_panel == 1L && p22$ledger$clump$kept == 8L && p22$ledger$instruments == 8L && nrow(p22$selection) == 14L)
p22s <- pair_input(paste0("lib:", lx$trait), "lib:Outcome for local", 5e-10)
fresh <- ld_clump_here(cand$snp[cand$p < 5e-10], cand$p[cand$p < 5e-10], 0.001, 10000, "EUR")
expect("opengwas local clump pair: a stricter study threshold cuts the saved selection, and the result equals clumping those candidates afresh",
       p22s$ledger$below_threshold == nrow(fresh) && p22s$ledger$clump$kept == sum(fresh$status == "index") && p22s$ledger$clump$clumped == sum(fresh$status == "clumped") &&
       setequal(p22s$input$data$variant_id, fresh$variant_id[fresh$status == "index"]) && p22s$ledger$instruments == sum(fresh$status == "index"))
h22 <- studio_new("OpenGWAS clumped here")
r22 <- studio_run(h22, p22$exposure_label, p22$outcome_label, p22$input, c("ivw_mre", "egger"), idx = 1L)
save_instrument_selection(h22$dir, r22$pair_folder, p22$selection)
ui_state_set(h22$con, "pairs", list(list(ex = paste0("lib:", lx$trait), out = "lib:Outcome for local", exName = p22$exposure_label, outName = p22$outcome_label, folder = r22$pair_folder, methods = list("ivw_mre", "egger"), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = p22$ledger, anc_x = "EUR", anc_y = "EUR")))
pk22 <- build_package(h22, c("report", "text", "refs", "strobe"), 140)
htm22 <- paste(readLines(file.path(pk22$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("opengwas local clump report: methods say downloaded unclumped (clump = 0) and clumped here with PLINK, with the outcome lookup, and never 'clumped on the OpenGWAS server'",
       grepl("downloaded from OpenGWAS without clumping", htm22, fixed = TRUE) && grepl("clumped with PLINK", htm22, fixed = TRUE) && grepl("extract_outcome_data", htm22, fixed = TRUE) &&
       grepl("5 variants were removed as correlated", htm22, fixed = TRUE) && !grepl("clumped on the OpenGWAS server with", htm22, fixed = TRUE) && identical(r22$replay$status, "identical") &&
       file.exists(file.path(h22$dir, "03_pairs", r22$pair_folder, "01_data", "instrument_selection.csv")))
# ---- Round 25: three places to clump, outcome names, ties at P = 0 ----
lv <- ogw_get_exposure("sim-exp-1", tok_ok, where = "server_live")
c25 <- tail(Filter(function(x) x[[1]] == "tophits", mock22$calls), 1)[[1]]
expect("clump choice: live server clumping asks OpenGWAS with force_server, is named '... live clump' and says so in its record",
       isTRUE(c25[[5]]) && identical(c25[[4]], 1) && grepl("live clump$", lv$trait) && identical(lv$source$params$where, "server_live") && grepl("force_server = TRUE", lv$source$note) && is.null(lv$clump))
expect("clump choice: stored clumps are recognised only for a plain server request with the default settings",
       ogw_stored_clumps(sx22$source$params) && !ogw_stored_clumps(lv$source$params) && !ogw_stored_clumps(lx$source$params) && !ogw_stored_clumps(modifyList(sx22$source$params, list(kb = 10001))) && grepl("stored, pre-computed clumps", sx22$source$note))
lo_lv <- ogw_get_outcome("sim-out-1", lv$trait, tok_ok); lo_sx <- ogw_get_outcome("sim-out-1", "Server clumped 22", tok_ok)
expect("outcome names: copies of one exposure (live, renamed, clumped here) get different default outcome names; the default-named exposure keeps the old ID tag",
       !identical(lo_lv$trait, lo_sx$trait) && grepl("for sim-exp-1 live clump$", lo_lv$trait) && grepl("for sim-exp-1 clumped here$", ogw_get_outcome("sim-out-1", lx$trait, tok_ok)$trait) && grepl("for Server clumped 22$", lo_sx$trait) &&
       { e0 <- .ogw_label("Simulated body mass", "sim-exp-1"); dn <- if (dir.exists(file.path(library_dir(), e0))) ogw_get_outcome("sim-out-1", e0, tok_ok)$trait else "for sim-exp-1"; grepl("for sim-exp-1$", dn) })
p25 <- pair_input(paste0("lib:", lv$trait), paste0("lib:", lo_lv$trait), 5e-8)
h25 <- studio_new("OpenGWAS live clump"); r25 <- studio_run(h25, p25$exposure_label, p25$outcome_label, p25$input, c("ivw_mre"), idx = 1L)
ui_state_set(h25$con, "pairs", list(list(ex = paste0("lib:", lv$trait), out = paste0("lib:", lo_lv$trait), exName = p25$exposure_label, outName = p25$outcome_label, folder = r25$pair_folder, methods = list("ivw_mre"), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = p25$ledger, anc_x = "EUR", anc_y = "EUR")))
pk25 <- build_package(h25, c("report", "strobe"), 140); htm25 <- paste(readLines(file.path(pk25$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("clump choice report: live server clumping is described as such (force_server), not as stored clumps", grepl("clumped live on the OpenGWAS server with ieugwasr::tophits (force_server = TRUE)", htm25, fixed = TRUE) && !grepl("pre-computed", htm25, fixed = TRUE))
# ties: two variants of LD block A with P = 0. Tie sensitivity order: the larger |beta/SE| leads, whichever comes first.
t1 <- ld_clump_here(c("rs9100001", "rs9100003", "rs9300001"), c(0, 0, 1e-9), z = c(5, 50, 6), tie_order = "strongest")
t2 <- ld_clump_here(c("rs9100001", "rs9100003", "rs9300001"), c(0, 0, 1e-9), z = c(50, 5, 6), tie_order = "strongest")
expect("ties at P = 0, strongest order: the variant with the larger |beta/SE| becomes the lead, and the tie-break is recorded",
       identical(t1$variant_id[t1$status == "index"][1], "rs9100003") && identical(t2$variant_id[t2$status == "index"][1], "rs9100001") && grepl("beta/SE", attr(t1, "settings")$tie_break))
# Round 28: the default order is that of ieugwasr::ld_clump_local and the OpenGWAS server (PLINK keeps the tied variant first in the panel)
s1 <- ld_clump_here(c("rs9100001", "rs9100003", "rs9300001"), c(0, 0, 1e-9), z = c(5, 50, 6)); s2 <- ld_clump_here(c("rs9100003", "rs9100001", "rs9300001"), c(0, 0, 1e-9))
ie28 <- ieugwasr::ld_clump_local(data.frame(rsid = c("rs9100001", "rs9100003", "rs9300001"), pval = c(0, 0, 1e-9)), clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, bfile = file.path(fix_dir, "EUR"), plink_bin = plink_bin())
expect("ties, default order (Round 28): PLINK keeps the tied variant first in the panel whatever the input order, beta/SE is ignored, and the result equals ieugwasr::ld_clump_local",
       identical(s1$variant_id[s1$status == "index"][1], "rs9100001") && "rs9100001" %in% s2$variant_id[s2$status == "index"] && !"rs9100003" %in% s2$variant_id[s2$status == "index"] && setequal(ie28$rsid, s1$variant_id[s1$status == "index"]) &&
       grepl("order of the reference panel", attr(s1, "settings")$tie_break) && identical(attr(s1, "settings")$tied, 2L))
ts28 <- ld_tie_summary(s1, t1)
expect("ties: the summary names the lead the tie order kept and the one with the larger |beta/SE| (1 region)", ts28$regions == 1L && identical(unlist(ts28$out), "rs9100001") && identical(unlist(ts28[["in"]]), "rs9100003") && ts28$tied == 2L)
expect("ties: the strongest order without beta/SE is refused", grepl("needs beta / SE", err_of(ld_clump_here(c("rs9100001", "rs9100003"), c(0, 0), tie_order = "strongest"))))
pu <- cand$p * (1 + seq_along(cand$p) / 1000); nz <- ld_clump_here(cand$snp, pu, 0.001, 10000, "EUR", z = rev(seq_along(pu)), tie_order = "strongest"); n0 <- ld_clump_here(cand$snp, pu, 0.001, 10000, "EUR")
expect("ties: when no P values are equal, the strongest order changes nothing", identical(nz$status, n0$status) && identical(nz$index_variant, n0$index_variant) && identical(attr(nz, "settings")$tie_break, "none needed"))
# ---- Round 28: tied P values from OpenGWAS (stored as 1e-200), the tie notice and the tie sensitivity analysis ----
th28 <- th22[th22$rsid != "rs9100099", ]
th28$p[th28$rsid %in% c("rs9100001", "rs9100002")] <- 1e-200; th28$se[th28$rsid == "rs9100001"] <- abs(th28$beta[th28$rsid == "rs9100001"]) / 20; th28$se[th28$rsid == "rs9100002"] <- abs(th28$beta[th28$rsid == "rs9100002"]) / 40
mock28 <- new.env(); mock28$calls <- list()
BE28 <- list(info = BE$info,
  tophits = function(id, p, r2, kb, pop, jwt, clump = 1, force = FALSE) { chk(jwt); if (clump == 0) th28 else th28[th28$rsid %in% c("rs9100001", "rs9200001", "rs9300001", "rs9300002"), ] },
  outcome = function(snps, id, proxies, rsq, maf, jwt) { chk(jwt); mock28$calls[[length(mock28$calls) + 1L]] <- snps; i <- match(snps, th28$rsid)
    data.frame(SNP = snps, chr = "1", pos = th28$position[i], beta.outcome = 0.3 * th28$beta[i] + rnorm(length(i), 0, 0.003), se.outcome = 0.008, samplesize.outcome = 230000, pval.outcome = 0.01, eaf.outcome = 0.35,
               effect_allele.outcome = "A", other_allele.outcome = "C", outcome = "Simulated heart disease", id.outcome = id, proxy.outcome = FALSE, target_snp.outcome = NA_character_, proxy_snp.outcome = NA_character_, stringsAsFactors = FALSE) })
options(mrgui.ogw_backend = BE28)
x28 <- ogw_get_exposure("sim-exp-1", tok_ok, where = "local", name = "Ties 28")
x28d <- library_dataset(x28$trait)$data$snp
expect("ties from OpenGWAS: clumped in the server's order (rs9100001, first in the panel, leads block A), the tie decided 1 region and the alternative selection is saved",
       "rs9100001" %in% x28d && !"rs9100002" %in% x28d && isTRUE(x28$clump$ties$regions == 1) && identical(unlist(x28$clump$ties[["in"]]), "rs9100002") && x28$clump$ties$tied >= 2 &&
       file.exists(file.path(library_dir(), x28$trait, "instrument_selection_ties.parquet")) && grepl("order of the reference panel", x28$clump$tie_break))
o28 <- ogw_get_outcome("sim-out-1", x28$trait, tok_ok, name = "Outcome for ties 28")
expect("ties from OpenGWAS: the outcome lookup also asks for the 1 alternative lead, which is not counted as an instrument",
       "rs9100002" %in% tail(mock28$calls, 1)[[1]] && length(tail(mock28$calls, 1)[[1]]) == length(x28d) + 1L && o28$lookup$requested == length(x28d) && o28$lookup$tie_alternatives == 1L && o28$lookup$tie_alternatives_found == 1L)
p28 <- pair_input(paste0("lib:", x28$trait), paste0("lib:", o28$trait), 5e-8)
q28 <- pair_input(paste0("lib:", x28$trait), paste0("lib:", o28$trait), 5e-8, ties = "strongest")
expect("ties pair: the primary analysis uses the server's lead; the ledger reports 1 tied region",
       "rs9100001" %in% p28$input$data$variant_id && !"rs9100002" %in% p28$input$data$variant_id && isTRUE(p28$ledger$ties$regions == 1) && is.null(p28$ledger$tie_order))
expect("ties pair, sensitivity: the lead with the larger |beta/SE| replaces it, every other instrument is the same, and the selection says so",
       "rs9100002" %in% q28$input$data$variant_id && !"rs9100001" %in% q28$input$data$variant_id && setequal(setdiff(p28$input$data$variant_id, "rs9100001"), setdiff(q28$input$data$variant_id, "rs9100002")) &&
       identical(q28$ledger$tie_order, "strongest") && identical(attr(q28$selection, "settings")$tie_order, "strongest"))
h28 <- studio_new("Ties 28")
r28 <- studio_run(h28, p28$exposure_label, p28$outcome_label, p28$input, c("ivw_mre", "egger"), idx = 1L)
rt28 <- studio_run_ties(h28, paste0("lib:", x28$trait), paste0("lib:", o28$trait), c("ivw_mre", "egger", "presso"), idx = 2L)
pt28 <- pair_payload(h28$dir, rt28$pair_folder)
expect("tie sensitivity run: its own labelled pair folder, detectors left out, replay identical, and the estimate differs from the primary",
       grepl("\\[ties: strongest\\]", rt28$label) && !identical(rt28$pair_folder, r28$pair_folder) && !any(vapply(pt28$estimates, function(e) identical(e$method, "presso"), TRUE)) &&
       identical(pt28$replay$status, "identical") && rt28$ties$regions == 1L &&
       { a <- Filter(function(e) e$method == "ivw_mre", pair_payload(h28$dir, r28$pair_folder)$estimates)[[1]]$estimate; b <- Filter(function(e) e$method == "ivw_mre", pt28$estimates)[[1]]$estimate; !isTRUE(all.equal(a, b)) })
ui_state_set(h28$con, "pairs", list(list(ex = paste0("lib:", x28$trait), out = paste0("lib:", o28$trait), exName = p28$exposure_label, outName = p28$outcome_label, folder = r28$pair_folder, methods = list("ivw_mre", "egger"), ci = 0.95, pthr = 5e-8, simulated = FALSE, ledger = p28$ledger, anc_x = "EUR", anc_y = "EUR")))
ui_state_set(h28$con, "sens", list(list(parent = r28$pair_folder, kind = "ties", remove = list(), folder = rt28$pair_folder, label = rt28$label, kept = rt28$kept, regions = rt28$ties$regions, swaps = list(out = rt28$ties$out, `in` = rt28$ties[["in"]]))))
pk28 <- build_package(h28, c("report", "text", "strobe"), 140); htm28 <- paste(readLines(file.path(pk28$dir, "report.html"), warn = FALSE), collapse = "\n")
expect("tie report: the Methods say PLINK kept the tied variant first in the reference panel, as ieugwasr and the server do, that the sensitivity analysis left out 1 and added 1 instrument (Round 29 wording), and the Results give the sensitivity IVW",
       grepl("PLINK kept the one that comes first in the reference panel, as in ieugwasr::ld_clump_local and OpenGWAS server clumping", htm28, fixed = TRUE) && grepl("A sensitivity analysis that kept the variant with the larger |beta/SE| instead left out 1 of the primary instruments and added 1 other", htm28, fixed = TRUE) && grepl("(1 primary instrument left out, 1 other added; ", htm28, fixed = TRUE) && !grepl("decided the lead variant", htm28, fixed = TRUE) &&
       grepl("Sensitivity analysis for tied P values", htm28, fixed = TRUE) && !grepl("were ordered by the absolute value of beta/SE", htm28, fixed = TRUE))
expect("tie sensitivity refusals: no tie that changes a lead, server clumping, simulated traits",
       grepl("No tied P value changed an instrument", err_of(pair_input(paste0("lib:", lx$trait), "lib:Outcome for local", 5e-8, ties = "strongest"))) &&
       grepl("clumped on the OpenGWAS server", err_of(pair_input("lib:Server clumped 22", paste0("lib:", lo_sx$trait), 5e-8, ties = "strongest"))) &&
       grepl("Simulated traits", err_of(pair_input("expA", "outA", ties = "strongest"))))
ps28 <- pair_input("lib:Server clumped 22", paste0("lib:", lo_sx$trait), 5e-8)
expect("ties with server clumping: the ledger says the ties could not be checked (server = TRUE) and counts the instruments that share a P value",
       isTRUE(ps28$ledger$ties$server) && is.numeric(ps28$ledger$ties$tied))
expect("Round 31: the ledger says whether the smallest P is itself shared, and the wording names a storage floor only then",
       is.logical(ps28$ledger$ties$floor_tied) && length(ps28$ledger$ties$floor_tied) == 1L &&
       identical(p_store_phrase(list(floor_tied = TRUE, floor = 0)), "OpenGWAS stores P values with a limited number of digits, and the smallest ones as the same value, 0") &&
       identical(p_store_phrase(list(floor_tied = FALSE, floor = 0)), "OpenGWAS stores P values with a limited number of digits") &&
       !grepl("below 0", p_store_phrase(list(floor_tied = TRUE, floor = 0)), fixed = TRUE))
expect("ties: an outcome downloaded before Round 28 is refused for the tie sensitivity analysis with a clear message",
       { rp <- file.path(library_dir(), o28$trait, "recipe.json"); rj <- jsonlite::fromJSON(rp, simplifyVector = FALSE); bak <- readLines(rp, warn = FALSE); rj$lookup[["tie_alternatives"]] <- NULL
         writeLines(jsonlite::toJSON(rj, auto_unbox = TRUE, pretty = TRUE, null = "null"), rp); m <- err_of(pair_input(paste0("lib:", x28$trait), paste0("lib:", o28$trait), 5e-8, ties = "strongest")); writeLines(bak, rp); grepl("Download the outcome again", m) })
options(mrgui.ogw_backend = NULL); ogw_token_set(NULL)
Sys.unsetenv("MRGUI_LD_DIR")

# ---- Round 20: STROBE-MR checklist ----
expect("strobe: the item list has 20 items with 30 sub-items (10d in two parts), 43 rows", nrow(STROBE_ITEMS) == 43L && !anyDuplicated(STROBE_ITEMS$item) &&
       length(unique(sub("[a-e](-i+)?$", "", STROBE_ITEMS$item))) == 20L)
pk20 <- build_package(h19, c("report", "text", "refs", "strobe"), 140)
c19 <- lapply(ui_state_get(h19$con, "pairs"), function(p) pair_context(h19, p))
rows19 <- strobe_rows(c19, report_text(c19, "Clumped here"), list(aim = "Does X cause Y?", prereg = FALSE), c("report", "strobe"))
g19 <- function(i, col = "status") rows19[[col]][rows19$item == i]
expect("strobe: every row has a known status and the rows follow the checklist order", identical(rows19$item, STROBE_ITEMS$item) && all(rows19$status %in% c("FILLED", "PARTLY", "AUTHOR")))
expect("strobe: instrument selection is filled with the local clumping and the counts", identical(g19("4c"), "FILLED") && grepl("clumped on the analyst's computer with PLINK", g19("4c", "provides"), fixed = TRUE) && grepl("5 removed by clumping", g19("10a", "provides"), fixed = TRUE))
expect("strobe: the recorded aim is quoted, one pair needs no multiple-testing correction, no direction test means AUTHOR",
       identical(g19("3"), "PARTLY") && grepl("Does X cause Y?", g19("3", "provides"), fixed = TRUE) && identical(g19("6e"), "FILLED") && identical(g19("13c"), "AUTHOR"))
expect("strobe: items only the authors can write stay AUTHOR, and 13e asks for leave-one-out when it was not run",
       all(vapply(c("1", "2", "4e", "5", "16a", "16b", "16c", "18", "20"), function(i) identical(g19(i), "AUTHOR"), TRUE)) && identical(g19("13e"), "AUTHOR") && grepl("Turn on Leave-one-out", g19("13e", "provides"), fixed = TRUE))
expect("strobe: plots count as filled only when figures are in the package", identical(g19("11d"), "FILLED") && identical(strobe_rows(c19, report_text(c19, "x"), list(), c("text"))$status[rows19$item == "11d"], "PARTLY"))
c15 <- lapply(list(list(ex = "expA", out = "outA", exName = p8$exposure_label, outName = p8$outcome_label, folder = r15$pair_folder, methods = list("ivw_mre", "steiger"), ci = 0.95, pthr = 5e-8, simulated = TRUE)), function(p) pair_context(h15, p))
r15s <- strobe_rows(c15, report_text(c15, "S"), list(), c("strobe"))
expect("strobe: a Steiger run fills 13c, and a simulated study is marked as such", identical(r15s$status[r15s$item == "13c"], "FILLED") && isTRUE(attr(r15s, "simulated")) && !isTRUE(attr(rows19, "simulated")))
sf <- file.path(pk20$dir, "strobe_mr_checklist.docx"); zl <- utils::unzip(sf, list = TRUE)$Name
td20 <- tempfile(); utils::unzip(sf, exdir = td20); dx <- paste(readLines(file.path(td20, "word", "document.xml"), warn = FALSE), collapse = "")
expect("strobe: the package holds the checklist as a Word file with every item and the official citations",
       file.exists(sf) && "word/document.xml" %in% zl && grepl("STROBE-MR checklist", dx, fixed = TRUE) && grepl("10d-ii", dx, fixed = TRUE) && grepl("34702754", dx, fixed = TRUE) && grepl("34698778", dx, fixed = TRUE))
expect("strobe: the package reports the counts, which add up to 43", !is.null(pk20$strobe) && pk20$strobe$filled + pk20$strobe$partly + pk20$strobe$author == 43L && !any(grepl("was not built", pk20$notes)))

# ---- Round 21: feasibility engine ----
fr <- feas_rules(); spec_ids <- utils::read.csv(file.path(dirname(dirname(args_dir)), "mr_feasibility_rules.csv"), stringsAsFactors = FALSE)$rule_id
known_tests <- c("divw_condition_below", "several_pairs", "ancestry_differs", "steiger_data", "iv_below_min", "iv_below_recommended", "mean_f_below", "i2gx_below", "i2gx_near", "overlap", "palindromic_dropped", "palindromic_kept", "binary_exposure", "presso_nbdist", "heavy_runtime", "q_significant", "egger_intercept", "steiger_reversed")
expect("feasibility: every rule id is unique and comes from the project's rule list, every active rule has a known test and a source",
       !anyDuplicated(fr$rule_id) && all(fr$rule_id %in% spec_ids) && all(fr$test[fr$active_v1 == "yes"] %in% known_tests) && all(nzchar(fr$source)) && all(spec_ids %in% fr$rule_id))
fx_base <- pair_facts("expA", "outA")
fz <- feasibility(list(fx_base))
expect("feasibility: the demo pair (40 strong instruments) leaves every method available, IVW without any rule (a method whose package is not in this runtime is locked for that reason only, Round 70b)", all(vapply(fz$methods, function(x) x$state %in% c("ok", "check") || all(vapply(x$reasons, `[[`, "", "rule") == "PKG"), TRUE)) && identical(fz$methods$ivw_mre$state, "ok") && !length(fz$methods$ivw_mre$reasons))
mk <- function(...) modifyList(fx_base, list(...))
f3 <- feasibility(list(mk(instruments = 3)))
expect("feasibility: 3 instruments lock MR-PRESSO and radial MR-Egger (R20) and warn for MR-Egger (R21); IVW stays usable with the R21b warning",
       f3$methods$presso$state == "locked" && f3$methods$radegger$state == "locked" && f3$methods$egger$state == "check" && "R21" %in% vapply(f3$methods$egger$reasons, `[[`, "", "rule") &&
       f3$methods$ivw_mre$state == "check" && identical(vapply(f3$methods$ivw_mre$reasons, `[[`, "", "rule"), "R21b"))
f4 <- feasibility(list(mk(instruments = 4)))
expect("Round 30: R21b warns for IVW below 4 instruments only, and says why (heterogeneity cannot be estimated reliably)",
       !"R21b" %in% vapply(f4$methods$ivw_mre$reasons, `[[`, "", "rule") && grepl("cannot be estimated reliably", f3$methods$ivw_mre$reasons[[1]]$msg) && grepl("3 instruments", f3$methods$ivw_mre$reasons[[1]]$msg))
fw <- feasibility(list(mk(mean_f = 8.2, i2gx = 0.8)))
rl <- function(fe, m) vapply(fe$methods[[m]]$reasons, `[[`, "", "rule")
expect("feasibility: mean F below 10 warns every estimator except debiased IVW and MR-RAPS (R23); I2GX below 0.9 warns the Egger methods only (R24)",
       all(vapply(setdiff(FEAS_CAUSAL, c("divw", "raps")), function(m) "R23" %in% rl(fw, m), TRUE)) && !("R23" %in% rl(fw, "raps")) && !("R23" %in% rl(fw, "steiger")) && "R24" %in% rl(fw, "egger") && "R24" %in% rl(fw, "radegger") && !("R24" %in% rl(fw, "ivw_mre")))
expect("feasibility: the F message does not claim that F above 10 means no bias", grepl("does not rule out bias", fw$methods$ivw_mre$reasons[[1]]$msg, fixed = TRUE))
expect("Round 32: the F message says the Egger, median and mode methods are biased more and names selection in the same exposure GWAS",
       grepl("MR-Egger and the weighted median and mode more than IVW", fw$methods$ivw_mre$reasons[[1]]$msg, fixed = TRUE) && grepl("winner's curse", fw$methods$ivw_mre$reasons[[1]]$msg, fixed = TRUE))
fn92 <- feasibility(list(mk(i2gx = 0.92))); fn89 <- feasibility(list(mk(i2gx = 0.89))); fn96 <- feasibility(list(mk(i2gx = 0.96)))
expect("Round 32: R24b notes radial MR-Egger for I2GX from 0.9 to below 0.95 only, as information; MR-Egger gets no note there and R24 takes over below 0.9",
       "R24b" %in% rl(fn92, "radegger") && !("R24" %in% rl(fn92, "radegger")) && !("R24b" %in% rl(fn92, "egger")) && fn92$methods$radegger$state == "ok" &&
       !("R24b" %in% rl(fn89, "radegger")) && "R24" %in% rl(fn89, "radegger") && !("R24b" %in% rl(fn96, "radegger")) &&
       grepl("about 11%", Filter(function(r) r$rule == "R24b", fn92$methods$radegger$reasons)[[1]]$msg, fixed = TRUE))
freal <- mk(simulated = FALSE, ex_id = "lib:A", out_id = "lib:B", ukb_both = FALSE)
fo1 <- feasibility(list(freal)); fo2 <- feasibility(list(freal), overlap = setNames(list(list(answer = "none")), "lib:A | lib:B"))
fo3 <- feasibility(list(freal), overlap = setNames(list(list(answer = "some", pct = 20)), "lib:A | lib:B"))
# Round 35: R24 says where the MR-Egger bias goes, by the overlap answer; R21 adds the median note for 3 or 4 instruments.
fe24 <- function(ans) { ov <- if (is.null(ans)) list() else setNames(list(list(answer = ans)), "lib:A | lib:B")
  r <- feasibility(list(modifyList(freal, list(i2gx = 0.8))), overlap = ov)$methods$egger$reasons; r[[which(vapply(r, `[[`, "", "rule") == "R24")]]$msg }
expect("Round 35: R24 says diluted toward the null without declared overlap, toward the observational association with some or unknown overlap",
       grepl("diluted toward the null", fe24(NULL)) && grepl("diluted toward the null", fe24("none")) &&
       grepl("toward the observational association", fe24("some")) && grepl("toward the observational association", fe24("unknown")) && !grepl("\\{", fe24("some")))
# Round 37: dIVW is not warned for weak-instrument bias (it removes it) but warns when its condition number is below 20.
fd1 <- feasibility(list(mk(mean_f = 8.2, instruments = 30))); fd2 <- feasibility(list(mk(mean_f = 4, instruments = 9)))
expect("Round 37: R23 skips debiased IVW; R23b warns when (mean F - 1) x sqrt(k) is below 20, not otherwise",
       !("R23" %in% rl(fd1, "divw")) && "R23" %in% rl(fd1, "ivw_mre") && !("R23b" %in% rl(fd1, "divw")) && "R23b" %in% rl(fd2, "divw") &&
       grepl("is 9.0 \\(below 20\\)", fd2$methods$divw$reasons[[which(rl(fd2, "divw") == "R23b")]]$msg) && !("R23b" %in% rl(fd2, "ivw_mre")))
expect("Round 37: debiased IVW is offered with its schema, citation and Methods sentence",
       "divw" %in% names(param_schemas()) && identical(METHOD_REFS$divw, "ye2021") && grepl("doi:10.1214/20-AOS2027.$", ref_text(REFS$ye2021)))
expect("Round 38: MR-RAPS is offered with its two settings, citation and a separate setting name from debiased IVW",
       all(c("raps_dispersion", "loss") %in% names(param_schemas()$raps)) && identical(METHOD_REFS$raps, "zhao2020") && grepl("doi:10.1214/19-AOS1866.$", ref_text(REFS$zhao2020)) &&
       "dispersion" %in% names(param_schemas()$divw) && !("dispersion" %in% names(param_schemas()$raps)))
expect("Round 38: run_params passes the enum settings of the new methods and refuses values outside their schema",
       identical(run_params(list(loss = "tukey", raps_dispersion = "none", dispersion = "none"))[c("loss", "raps_dispersion", "dispersion")], list(loss = "tukey", raps_dispersion = "none", dispersion = "none")) &&
       is.null(run_params(list(loss = "bogus"))$loss) && identical(run_params(list(ivw_model = "fixed"))$ivw_model, "fixed"))
expect("Round 40: contamination mixture and MR-Lasso are offered with their settings, citations and the R21 instrument-count warning",
       all(c("psi", "ci_step") %in% names(param_schemas()$conmix)) && "lasso_dist" %in% names(param_schemas()$lasso) && identical(METHOD_REFS$conmix, "burgess2020") && identical(METHOD_REFS$lasso, "rees2019") &&
       "R21" %in% rl(f3, "conmix") && "R21" %in% rl(f3, "lasso") && identical(run_params(list(psi = 2, ci_step = 0.001, lasso_dist = "t-dist"))[c("psi", "ci_step", "lasso_dist")], list(psi = 2, ci_step = 0.001, lasso_dist = "t-dist")))
expect("Round 41: MR-Egger with SIMEX correction is offered with its citation, R21 coverage and as the R24 recommendation",
       "simex_b" %in% names(param_schemas()$egger_simex) && identical(METHOD_REFS$egger_simex, "bowden2016i2") && "R21" %in% rl(f3, "egger_simex") &&
       grepl("SIMEX correction (sensitivity analysis)", feas_rules()$recommend[feas_rules()$rule_id == "R24"], fixed = TRUE) && !any(grepl("SIMEX correction (NOT BUILT)", feas_rules()$recommend, fixed = TRUE)) &&
       identical(run_params(list(simex_b = 300))$simex_b, 300L))
r21m <- function(n, m) { r <- feasibility(list(mk(instruments = n)))$methods[[m]]$reasons; r[[which(vapply(r, `[[`, "", "rule") == "R21")]]$msg }
expect("Round 35: R21 adds the median false-positive note for the weighted and simple median with 3 or 4 instruments only",
       grepl("about 10% instead of 5%", r21m(3, "wm")) && grepl("about 10% instead of 5%", r21m(4, "smed")) && !grepl("about 10%", r21m(5, "wm")) &&
       !grepl("about 10%", r21m(3, "egger")) && !grepl("\\{", r21m(5, "wm")))
expect("feasibility: sample overlap warns until answered 'none' (R28), and quotes a stated share", "R28" %in% rl(fo1, "ivw_mre") && !("R28" %in% rl(fo2, "ivw_mre")) && grepl("about 20%", fo3$methods$ivw_mre$reasons[[1]]$msg, fixed = TRUE))
expect("feasibility: different declared ancestries warn (R13)", "R13" %in% rl(feasibility(list(mk(anc_x = "EUR", anc_y = "EAS"))), "wm"))
fp <- feasibility(list(fx_base), params = list(nbdist = 30))
expect("feasibility: too few MR-PRESSO simulated datasets lock it with the fix (R31)", fp$methods$presso$state == "locked" && grepl("Raise the number of simulated datasets", fp$methods$presso$reasons[[which(rl(fp, "presso") == "R31")]]$msg, fixed = TRUE))
f2p <- feasibility(list(fx_base, pair_facts("expA", "outC")))
expect("feasibility: several pairs give the multiple-testing note once (R05b) and name the pair in pair-level messages", sum(rl(f2p, "ivw_mre") == "R05b") == 1L && grepl("Bonferroni", f2p$methods$ivw_mre$reasons[[1]]$msg, fixed = TRUE))
expect("feasibility: a binary exposure adds the interpretation note (R30), not a lock", "R30" %in% rl(feasibility(list(mk(scale_x = "log_or"))), "ivw_mre") && feasibility(list(mk(scale_x = "log_or")))$methods$ivw_mre$state == "ok")
# consistency with the adapters' own checks: a method is locked by the engine exactly when its adapter blocks
din <- demo_input("expA", "outA"); dd <- din$input$data[din$input$data$mr_keep, ]
agree <- TRUE
for (k in 2:6) {
  sub_in <- din$input; sub_in$data <- dd[seq_len(k), ]   # not modifyList: it would merge the data frame column by column
  fk <- feasibility(list(mk(instruments = k, mean_f = mean((dd$beta_x[1:k] / dd$se_x[1:k])^2), i2gx = feas_i2gx(dd$beta_x[1:k], dd$se_x[1:k]), steiger = steiger_readiness(list(input = sub_in, binary = NULL)))))
  for (m in FEAS_METHODS) { a <- param_adapters()[[m]]; blk <- identical(a$check(sub_in, resolve_params(a, list()))$status, "block")
    if (blk != identical(fk$methods[[m]]$state, "locked")) { agree <- FALSE; cat("  disagree:", m, "k =", k, "\n") } }
}
expect("feasibility: the engine and every adapter agree on what is blocked, for 2 to 6 instruments", agree)
expect("feasibility: results rules match the saved diagnostics (R25 heterogeneity)", {
  q <- Filter(function(d) d$method == "ivw_mre" && d$name == "cochran_q", pl10$diagnostics)[[1]]; adv <- feasibility_results(pl10)
  identical("R25" %in% vapply(adv, `[[`, "", "rule"), q$p < 0.05) && !is.null(pair_payload(h10$dir, r10$pair_folder)$advice) })
expect("feasibility: the Steiger result rule fires only when the direction is not clearly supported (R27)", {
  sg <- Filter(function(d) d$method == "steiger" && d$name == "steiger_direction", pl15$diagnostics)[[1]]
  identical("R27" %in% vapply(feasibility_results(pl15), `[[`, "", "rule"), !(sg$value == 1 && sg$p < 0.05)) })

# ---- Round 26: files without a sample-size or allele-frequency column (found while replicating Bowden 2016) ----
expect("an imported pair with no sample-size and no frequency column can be paired and run", {
  nx <- exd[, c("snp", "ea", "oa", "beta", "se", "p")]; ny <- exd[, c("snp", "ea", "oa", "beta", "se", "p")]; ny$beta <- 0.3 * nx$beta + rnorm(nrow(nx), 0, .01); ny$p <- 2 * pnorm(-abs(ny$beta / ny$se))
  f1 <- tempfile(fileext = ".csv"); f2 <- tempfile(fileext = ".csv"); utils::write.csv(nx, f1, row.names = FALSE); utils::write.csv(ny, f2, row.names = FALSE)
  m6 <- list(snp = "snp", ea = "ea", oa = "oa", beta = "beta", se = "se", p = "p")
  invisible(import_save(f1, "nx.csv", ",", names(nx), m6, modifyList(fct("NoN X", "exposure", "yes"), list(frq = "none"))))
  invisible(import_save(f2, "ny.csv", ",", names(ny), m6, modifyList(fct("NoN Y", "outcome", "no"), list(frq = "none"))))
  pn <- pair_input("lib:NoN X", "lib:NoN Y", pthr = 1e-4)
  pn$ledger$kept >= 2 && pn$ledger$kept + pn$ledger$dropped_palindromic + pn$ledger$dropped_mismatch == pn$ledger$in_outcome && all(is.na(pn$input$data$n_x)) })

# ---- Round 26: keep palindromic variants when both files are declared forward strand ----
sx26 <- exd[, c("snp", "ea", "oa", "beta", "se", "p")]; sy26 <- sx26; sy26$beta <- 0.3 * sx26$beta + rnorm(nrow(sx26), 0, .01); sy26$p <- 2 * pnorm(-abs(sy26$beta / sy26$se))
pal26 <- sum(paste0(sx26$ea, sx26$oa) %in% c("AT", "TA", "CG", "GC") & sx26$p < 1e-4)
m26 <- list(snp = "snp", ea = "ea", oa = "oa", beta = "beta", se = "se", p = "p")
imp26 <- function(df, nm, role, ind, strand) { f <- tempfile(fileext = ".csv"); utils::write.csv(df, f, row.names = FALSE)
  invisible(import_save(f, basename(f), ",", names(df), m26, modifyList(fct(nm, role, ind), list(frq = "none", strand = strand)))) }
imp26(sx26, "Fwd X", "exposure", "yes", "forward"); imp26(sy26, "Fwd Y", "outcome", "no", "forward"); imp26(sy26, "Unk Y", "outcome", "no", "unknown")
pf <- pair_input("lib:Fwd X", "lib:Fwd Y", pthr = 1e-4); pu <- pair_input("lib:Fwd X", "lib:Unk Y", pthr = 1e-4)
expect("strand: both files forward keeps every palindromic variant (action 1) and counts them", pal26 > 0 && identical(pf$ledger$strand, "forward") && pf$ledger$dropped_palindromic == 0 && pf$ledger$kept_palindromic == pal26 && pf$ledger$harmonise_action == 1L)
expect("strand: one file not declared forward drops them as before (action 2)", identical(pu$ledger$strand, "unknown") && pu$ledger$dropped_palindromic == pal26 && pu$ledger$kept_palindromic == 0 && pu$ledger$harmonise_action == 2L)
expect("strand: the import warning says palindromic variants are kept only when both files are forward strand", {
  f <- tempfile(fileext = ".csv"); utils::write.csv(sx26, f, row.names = FALSE); h6 <- import_read_head(f)
  any(grepl("unless both files are declared to be on the forward strand", import_validate(f, h6$sep, h6$columns, m26, modifyList(fct("w", "exposure", "yes"), list(frq = "none")))$warnings)) &&
  any(grepl("You declared the forward strand", import_validate(f, h6$sep, h6$columns, m26, modifyList(fct("w", "exposure", "yes"), list(frq = "none", strand = "forward")))$warnings)) })
h26 <- studio_new("Strand test"); r26 <- studio_run(h26, pf$exposure_label, pf$outcome_label, pf$input, c("ivw_mre"), idx = 1L)
ui_state_set(h26$con, "pairs", list(list(ex = "lib:Fwd X", out = "lib:Fwd Y", exName = pf$exposure_label, outName = pf$outcome_label, folder = r26$pair_folder, methods = list("ivw_mre"), ci = 0.95, pthr = 1e-4, simulated = FALSE, ledger = pf$ledger, anc_x = "EUR", anc_y = "EUR")))
c26 <- lapply(ui_state_get(h26$con, "pairs"), function(p) pair_context(h26, p)); t26 <- report_text(c26, "Strand test"); s26 <- strobe_rows(c26, t26, list(), c("strobe"))
expect("strand: Methods name action 1 and the kept palindromic variants, Limitations warn that the strand was not verified", {
  m <- paste(t26$methods, collapse = " "); l <- paste(t26$limitations, collapse = " ")
  grepl("action 1: the analyst declared that both datasets report alleles on the forward strand", m, fixed = TRUE) && grepl(sprintf("%d palindromic variants were kept", pal26), m, fixed = TRUE) &&
  !grepl("action 2", m, fixed = TRUE) && grepl("This was not verified", l, fixed = TRUE) })
expect("strand: STROBE 4c and 10a report the kept palindromic variants", grepl("declared forward strand (action 1)", s26$provides[s26$item == "4c"], fixed = TRUE) && grepl(sprintf("%d palindromic kept", pal26), s26$provides[s26$item == "10a"], fixed = TRUE))
expect("strand: the feasibility note R29b fires for kept palindromic variants and R29 does not", {
  fr <- feasibility(list(pair_facts("lib:Fwd X", "lib:Fwd Y", 1e-4)))$methods$ivw_mre$reasons; ids <- vapply(fr, `[[`, "", "rule"); "R29b" %in% ids && !("R29" %in% ids) })

# ---- Round 30: default IVW model (TwoSampleMR::mr_ivw) in the Methods, and the Limitations sentence below 4 variants ----
expect("Round 30: the Methods name the default IVW model (mr_ivw, SE not below fixed effect) and cite the guidelines", {
  m <- paste(t26$methods, collapse = " ")
  grepl("multiplicative random effects, with the residual standard error not allowed below 1", m, fixed = TRUE) && grepl("TwoSampleMR mr_ivw)", m, fixed = TRUE) &&
    grepl("Burgess 2013, Burgess 2023, Hemani 2018", m, fixed = TRUE) && !grepl("fewer than 4 variants", paste(t26$limitations, collapse = " ")) })
p3 <- pf$input; p3$data <- utils::head(p3$data[order(p3$data$se_x / abs(p3$data$beta_x)), ], 3)
h30 <- studio_new("Three variants"); r30 <- studio_run(h30, pf$exposure_label, pf$outcome_label, p3, c("ivw_mre"), idx = 1L)
ui_state_set(h30$con, "pairs", list(list(ex = "lib:Fwd X", out = "lib:Fwd Y", exName = pf$exposure_label, outName = pf$outcome_label, folder = r30$pair_folder, methods = list("ivw_mre"), ci = 0.95, pthr = 1e-4, simulated = FALSE, ledger = pf$ledger, anc_x = "EUR", anc_y = "EUR")))
c30 <- lapply(ui_state_get(h30$con, "pairs"), function(p) pair_context(h30, p)); t30 <- report_text(c30, "Three variants")
e30 <- pair_payload(h30$dir, r30$pair_folder)$estimates[[1]]
f30 <- suppressWarnings(TwoSampleMR::mr_ivw(p3$data$beta_x, p3$data$beta_y, p3$data$se_x, p3$data$se_y, TwoSampleMR::default_parameters()))
expect("Round 30: with 3 variants IVW still runs (equal to TwoSampleMR::mr_ivw) and the Limitations say heterogeneity could not be estimated reliably",
       e30$n_ivs_used == 3L && abs(e30$se - f30$se) < 1e-9 && abs(e30$estimate - f30$b) < 1e-9 &&
       grepl("In one pair, fewer than 4 variants were used, so the heterogeneity between the variants could not be estimated reliably", paste(t30$limitations, collapse = " "), fixed = TRUE))

# ---- Round 31: the Limitations name only the robust estimators that were run ----
lim_of <- function(meth, nm) { h <- studio_new(nm); r <- studio_run(h, pf$exposure_label, pf$outcome_label, pf$input, meth, idx = 1L)
  ui_state_set(h$con, "pairs", list(list(ex = "lib:Fwd X", out = "lib:Fwd Y", exName = pf$exposure_label, outName = pf$outcome_label, folder = r$pair_folder, methods = as.list(meth), ci = 0.95, pthr = 1e-4, simulated = FALSE, ledger = pf$ledger)))
  paste(report_text(lapply(ui_state_get(h$con, "pairs"), function(p) pair_context(h, p)), nm)$limitations, collapse = " ") }
l31a <- paste(t30$limitations, collapse = " "); l31b <- lim_of(c("ivw_mre", "wm"), "R31 median"); l31c <- lim_of(c("ivw_mre", "wm", "mode"), "R31 median mode"); l31d <- lim_of(c("ivw_mre", "smode"), "R31 mode")
expect("Round 31: with IVW (and MR-Egger) only, the Limitations do not mention median or mode estimators",
       !grepl("median", l31a, fixed = TRUE) && !grepl("mode estimators", l31a, fixed = TRUE) && grepl("cannot be proven from summary data.", l31a, fixed = TRUE))
expect("Round 31: a median estimator alone is described with its majority condition, a mode estimator alone with its plurality condition, both together with both",
       grepl("The median estimators relax the third assumption only if a majority of the instruments are valid.", l31b, fixed = TRUE) && !grepl("mode", l31b, fixed = TRUE) &&
       grepl("The mode estimators relax the third assumption only if a plurality of the instruments are valid.", l31d, fixed = TRUE) && !grepl("median", l31d, fixed = TRUE) &&
       grepl("The median and mode estimators relax the third assumption only under their own conditions (a majority or plurality of valid instruments).", l31c, fixed = TRUE))

# ---- Round 39: leave-one-out and single-variant analyses (part of IVW) ----
h39 <- studio_new("Leave one out"); wi39 <- pair_input("expA", "outA")
r39 <- studio_run(h39, wi39$exposure_label, wi39$outcome_label, wi39$input, c("ivw_mre", "egger", "loo"), idx = 1L)
p39 <- pair_payload(h39$dir, r39$pair_folder); k39 <- p39$estimates[[1]]$n_ivs_used
ui_state_set(h39$con, "pairs", list(list(ex = "expA", out = "outA", exName = wi39$exposure_label, outName = wi39$outcome_label, folder = r39$pair_folder, methods = as.list(c("ivw_mre", "egger", "loo")), ci = 0.95, pthr = 5e-8, simulated = TRUE, ledger = wi39$ledger, anc_x = "EUR", anc_y = "EUR")))
c39 <- lapply(ui_state_get(h39$con, "pairs"), function(p) pair_context(h39, p)); lf39 <- loo_facts(c39[[1]])
expect("Round 39: leave-one-out runs with IVW, replays identically and gives one row per variant", identical(p39$replay$status, "identical") && nrow(lf39) == k39 && all(is.finite(lf39$loo_b)) && all(is.finite(lf39$wr)))
d39 <- wi39$input$data[wi39$input$data$mr_keep, ]; j39 <- match(lf39$id[1], d39$variant_id)
m39 <- TwoSampleMR::mr_ivw(d39$beta_x[-j39], d39$beta_y[-j39], d39$se_x[-j39], d39$se_y[-j39])
expect("Round 39: the stored leave-one-out estimate equals mr_ivw without that variant", abs(lf39$loo_b[1] - m39$b) < 1e-12 && abs(lf39$wr[1] - d39$beta_y[j39] / d39$beta_x[j39]) < 1e-12)
t39 <- report_text(c39, "Leave one out")
expect("Round 39: Results text and two supplementary legends describe the analysis", any(grepl("^Leave-one-out analysis .* the IVW estimate ranged from", t39$results)) && length(t39$supp[[1]]) == 2L &&
       grepl("^Supplementary Figure S1\\. Leave-one-out", t39$supp[[1]][1]) && grepl("^Supplementary Figure S2\\. Funnel plot", t39$supp[[1]][2]) && any(grepl("leave-one-out analysis", t39$methods, fixed = TRUE)))
pk39 <- build_package(h39, c("figs", "text", "tables", "report"), 140); fd39 <- file.path(pk39$dir, "figures")
expect("Round 39: the package has the leave-one-out and funnel plots, the supplementary sheet and both in the HTML report",
       all(file.exists(file.path(fd39, c("P1_leave_one_out.pdf", "P1_funnel.pdf")))) && any(grepl("Pair 1 leave-one-out", readLines(utils::unzip(file.path(pk39$dir, "supplement_tables.xlsx"), "xl/workbook.xml", exdir = tempfile("x39")), warn = FALSE), fixed = TRUE)) &&
       any(grepl("Supplementary figures", readLines(file.path(pk39$dir, "report.html"), warn = FALSE), fixed = TRUE)))
s39 <- strobe_rows(c39, t39, list(aim = "x", prereg = FALSE), c("report", "strobe"))
expect("Round 39: STROBE-MR 13e is filled when leave-one-out ran", identical(s39$status[s39$item == "13e"], "FILLED"))
big <- data.frame(id = sprintf("rs%d", 1:80), loo_b = seq(0.1, 0.5, length.out = 80), loo_se = 0.01, loo_p = 0.01, wr = 0, wr_se = 1, wr_p = 1)
expect("Round 39: above 60 variants the plot shows the 30 whose removal changes the estimate most", nrow(.loo_rows(big, 0.1)) == 30L && min(.loo_rows(big, 0.1)$loo_b) > 0.3)
expect("Round 39: the leave-one-out switch is excluded from reduced-set sensitivity runs and is not a causal method", !("loo" %in% FEAS_CAUSAL) && "loo" %in% FEAS_METHODS)

# ---- Round 42: outlier module (FR-90 to FR-103) ----
d42 <- demo_pair(args_dir)$input; out42 <- d42$data$variant_id[c(3, 11, 27)]; hit42 <- d42$data$variant_id %in% out42
d42$data$beta_y[hit42] <- d42$data$beta_y[hit42] + 8 * d42$data$se_y[hit42]
h42 <- studio_new("Outliers"); inp42 <- list(input = d42, exposure_label = "Exposure 42", outcome_label = "Outcome 42")
r42 <- studio_run(h42, "Exposure 42", "Outcome 42", d42, c("ivw_mre", "egger", "wm", "radial", "radegger", "lasso"), idx = 1L, params = list(ci_level = 0.95, nboot = 200L))
p42 <- pair_payload(h42$dir, r42$pair_folder)
expect("Round 42: detectors and re-analysis methods run and replay identically", all(vapply(r42$results, function(x) identical(x$status, "ok"), TRUE)) && identical(p42$replay$status, "identical"))
ra42 <- Filter(function(g) g$method == "radial" && g$name == "radial_alpha", p42$diagnostics)[[1]]
expect("Round 42: radial detectors test each variant at 0.05 divided by the number of variants by default, written in the script",
       abs(ra42$value - 0.05 / 40) < 1e-12 && grepl("divided by 40", ra42$note) && identical(pair_settings(h42$dir, r42$pair_folder)$alpha_correction, "bonferroni"))
sg42 <- outlier_signals(p42)
expect("Round 42: planted outliers give a heterogeneity signal; a detector not run gives no signal", isTRUE(sg42$q) && "q" %in% outlier_triggered(sg42) && is.na(sg42$presso))
s42 <- outlier_sets(p42)
expect("Round 42: one set per detector that ran, plus the consensus set", setequal(attr(s42, "detectors_run"), c("radial", "radegger", "lasso")) && "consensus" %in% s42$set)
expect("Round 42: the radial IVW set holds the three planted outliers", all(out42 %in% s42$ids[[which(s42$set == "radial")]]))
expect("Round 42: the cap (20%) and floor (5) decide eligibility", all(s42$eligible == (s42$n_removed > 0 & s42$n_removed <= 8 & 40 - s42$n_removed >= 5)))
pays42 <- list(); k42 <- 1L
for (i in which(s42$eligible)) { k42 <- k42 + 1L; rr <- studio_run_reduced(h42, inp42, s42$ids[[i]], OUTLIER_REANALYSIS, 0.95, k42, params = list(ci_level = 0.95, nboot = 200L)); pays42[[s42$set[i]]] <- pair_payload(h42$dir, rr$pair_folder) }
sc42 <- outlier_scenarios(p42, s42, pays42); se42 <- outlier_select(sc42, outlier_triggered(sg42))
expect("Round 42: comparison has all instruments first, with Q, I2, Rucker Q', intercept and excess reduction", sc42$set[1] == "all" && all(c("q", "i2", "q_rucker", "int_p", "excess_reduction") %in% names(sc42)) && sc42$excess_reduction[1] == 0)
expect("Round 42: removing the planted outliers lowers Q by more than the number removed", sc42$excess_reduction[sc42$set == "radial"] > 3)
expect("Round 42: the rule picks a set that clears the Q signal (step 2) and never reads the causal estimate",
       se42$step == 2L && sc42$q_p[sc42$set == se42$set] >= 0.05 && identical(outlier_select(transform(sc42, ivw = -ivw, egger = 99, wm = NA), "q"), se42))
rel42 <- outlier_write(h42$dir, r42$pair_folder, p42, s42, sc42, se42)
fl42 <- utils::read.csv(file.path(h42$dir, rel42, "01_flagged", "outliers_flagged.csv")); cm42 <- utils::read.csv(file.path(h42$dir, rel42, "03_tables", "outlier_comparison.csv"))
expect("Round 42: flag ledger has one row per instrument with each detector's flag and the removal decision", nrow(fl42) == 40L && all(c("flagged_radial", "flagged_lasso", "n_detectors", "removed_in_selected_set") %in% names(fl42)) && all(fl42$removed_in_selected_set[fl42$variant_id %in% out42] == "yes"))
expect("Round 42: comparison table marks the selected set and says the P values are descriptive", sum(cm42$decision == "selected by the rule") == 1L && all(grepl("descriptive only", cm42$note)))
rec42 <- outlier_record(r42$pair_folder, sg42, outlier_triggered(sg42), sc42, se42, NULL); rec42$chosen <- "all"; outlier_write_choice(h42$dir, rec42)
cm42b <- utils::read.csv(file.path(h42$dir, rel42, "03_tables", "outlier_comparison.csv")); fl42b <- utils::read.csv(file.path(h42$dir, rel42, "01_flagged", "outliers_flagged.csv"))
expect("Round 42: keeping all instruments is recorded as the analyst's choice beside the rule's", cm42b$decision[cm42b$set == "all"] == "chosen by the analyst" && cm42b$decision[cm42b$set == se42$set] == "selected by the rule, not adopted" && all(fl42b$removed_in_selected_set == "no"))
expect("Round 42: the study manifest still checks clean after the module wrote its files", { mc <- manifest_check(h42$dir); all(mc$status[grepl("08_outlier_analysis", mc$relative_path)] == "ok") && !any(mc$status %in% c("missing", "modified")) })
mk42 <- function(set, n, elig, qp, ip, agr, er) data.frame(set = set, label = set, n_removed = n, eligible = elig, q_p = qp, int_p = ip, presso_p = NA_real_, agreement = agr, excess_reduction = er, stringsAsFactors = FALSE)
sy42 <- mk42(c("all", "a", "b", "c"), c(0, 2, 2, 5), c(TRUE, TRUE, TRUE, TRUE), c(0.001, 0.2, 0.3, 0.6), c(0.01, 0.01, 0.2, 0.3), c(NA, 1, 2, 1), c(0, 4, 3, 9))
expect("Round 42: with two signals both must clear; fewest removed wins", outlier_select(sy42, c("q", "intercept"))$set == "b")
expect("Round 42: equal counts go to higher detector agreement", outlier_select(sy42, "q")$set == "b" && grepl("detector agreement", outlier_select(sy42, "q")$reason))
expect("Round 42: no set clears the signal, so the largest positive excess reduction is chosen (step 3)", { s <- sy42; s$q_p <- 0.001; r <- outlier_select(s, "q"); r$set == "c" && r$step == 3L })
expect("Round 42: no improvement keeps all instruments", { s <- sy42; s$q_p <- 0.001; s$excess_reduction <- c(0, -1, -2, 0); outlier_select(s, "q")$set == "all" })
expect("Round 42: no signal keeps all instruments and removes nothing", outlier_select(sy42, character())$set == "all")
expect("Round 42: an over-cap set is shown but never chosen", { s <- sy42; s$eligible <- c(TRUE, FALSE, FALSE, TRUE); outlier_select(s, "q")$set == "c" })

rec42$chosen <- se42$set; rec42$review <- "reviewed"; rec42$chosen_row <- Filter(function(x) identical(x$set, se42$set), rec42$scenarios)[[1]]
ce42 <- data.frame(method = "ivw_mre", estimate = 0.3, ci_lower = 0.2, ci_upper = 0.4, p = 0.001, scale = "beta", stringsAsFactors = FALSE)
tx42 <- outlier_results_text(rec42, ce42, "beta")
expect("Round 42: Results text names the signal, the rule's set, the removed variants and the descriptive Q P", grepl("started by Cochran's Q", tx42) && grepl("pre-specified rule chose the ", tx42) && all(vapply(rec42$chosen_row$ids, function(i) grepl(i, tx42, fixed = TRUE), TRUE)) && grepl("descriptive", tx42))
expect("Round 42: Methods text states the cap, floor, Bonferroni level, rule and that the primary used all instruments", all(vapply(c("20%", "at least 5", "divided by the number of variants", "without using the causal estimates", "primary analysis used all instruments"), function(k) grepl(k, OUTLIER_METHODS_TEXT, fixed = TRUE), TRUE)))
expect("Round 42: the radial Methods phrase says Bonferroni unless the correction is off", grepl("Bonferroni", .radial_level_txt(list(alpha = 0.05))) && grepl("no correction", .radial_level_txt(list(alpha = 0.01, alpha_correction = "none"))))

ui_state_set(h42$con, "outliers", list(rec42)); back42 <- ui_state_get(h42$con, "outliers")[[1]]
expect("Round 42: the saved outlier record keeps small P values (study state is stored at full precision)", isTRUE(abs(back42$scenarios[[1]]$q_p - rec42$scenarios[[1]]$q_p) <= 1e-12 * max(1, rec42$scenarios[[1]]$q_p)) && rec42$scenarios[[1]]$q_p < 1e-4)

# ---- Round 43: dominant variant and leave-subset-out checks (FR-105), R25 and R26 advice ----
sf43 <- subset_facts(p39); w43 <- d39$beta_x^2 / d39$se_y^2
expect("Round 43: the strongest variant's share of the IVW weight is stored with its leave-one-out estimate", !is.null(sf43) && identical(sf43$dom_id, d39$variant_id[which.max(w43)]) && abs(sf43$dom_share - max(w43) / sum(w43)) < 1e-12 && is.finite(sf43$dom_without))
expect("Round 43: leave-subset-out runs with 10 or more variants and replays identically (fixed seed)", sf43$lso && sf43$lso_n == 200L && sf43$lso_lo <= sf43$lso_med && sf43$lso_med <= sf43$lso_hi && sf43$lso_min <= sf43$lso_lo && identical(p39$replay$status, "identical"))
sd43 <- as.integer(sub("^set\\.seed\\(([0-9]+)\\)$", "\\1", grep("^set\\.seed\\(", readLines(Sys.glob(file.path(h39$dir, "03_pairs", r39$pair_folder, "02_scripts", "*loo*.R"))[1]), value = TRUE)[1]))
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(sd43); dd43 <- .loo_core(d39, "mre_floor")
expect("Round 43: the stored quantiles equal a recomputation with the seed written in the script", is.finite(sd43) && abs(dd43$value[dd43$name == "lso_median"] - sf43$lso_med) < 1e-12)
lp43 <- .loo_core(d39[1:6, ], "mre_floor")
expect("Round 43: fewer than 10 variants give the dominant share but no leave-subset-out", "dominant_share" %in% lp43$name && !any(grepl("^lso_", lp43$name)))
expect("Round 43: Results text reports the strongest variant and the leave-subset-out spread", any(grepl("The strongest variant, .* carried .*% of the total IVW weight", t39$results)) && any(grepl("Leaving out random subsets", t39$results)))
fr43 <- feas_rules()
expect("Round 43: R25 names MR-RAPS, contamination mixture, MR-Lasso and the outlier analysis; R26 leaves out MR-RAPS (balanced pleiotropy only)",
       all(vapply(c("MR-RAPS", "contamination mixture", "MR-Lasso", "outlier analysis"), function(k) grepl(k, fr43$recommend[fr43$rule_id == "R25"], fixed = TRUE), TRUE)) && !grepl("MR-RAPS", fr43$recommend[fr43$rule_id == "R26"]) && grepl("MR-Lasso", fr43$recommend[fr43$rule_id == "R26"]))

# ---- Round 44: outlier comparison and agreement figures (FR-98) ----
ag44 <- outlier_agreement(rec42)
expect("Round 44: agreement matrix has one column per detector that ran and the planted outliers flagged by radial IVW", !is.null(ag44) && setequal(colnames(ag44$m), c("Radial IVW", "Radial MR-Egger", "MR-Lasso")) && all(ag44$m[intersect(out42, rownames(ag44$m)), "Radial IVW"]))
expect("Round 44: the removed column marks exactly the adopted set's variants", setequal(rownames(ag44$m)[ag44$removed], intersect(rownames(ag44$m), unlist(Filter(function(x) identical(x$set, rec42$chosen), rec42$scenarios)[[1]]$ids))))
f44 <- tempfile(fileext = ".pdf"); ok44 <- tryCatch({ grDevices::pdf(f44, 7, 3); draw_outlier_comparison(rec42, "beta"); draw_outlier_agreement(ag44); grDevices::dev.off(); TRUE }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
expect("Round 44: both figures draw without error", ok44 && file.size(f44) > 1000)
lg44 <- outlier_legends(rec42, list(p = list(exName = "X", outName = "Y")), 3L, ag44)
expect("Round 44: legends numbered from the given figure, saying the estimate was not used by the rule", length(lg44) == 3L && grepl("^Supplementary Figure S3\\.", lg44[1]) && grepl("^Supplementary Figure S4\\.", lg44[2]) && grepl("^Supplementary Figure S5\\. Radial", lg44[3]) && grepl("did not use", lg44[1]))

# ---- Round 45: bidirectional MR ----
rv45 <- pair_input("outB", "expB"); fw45 <- pair_input("expB", "outB")
expect("Round 45: the reverse demo pair uses the outcome trait's own instruments and swaps the labels", identical(rv45$exposure_label, "Simulated outcome B") && identical(rv45$outcome_label, "Simulated exposure B") && !length(intersect(rv45$input$data$variant_id, fw45$input$data$variant_id)))
h45 <- studio_new("Both ways"); m45 <- c("ivw_mre", "egger")
ra45 <- studio_run(h45, fw45$exposure_label, fw45$outcome_label, fw45$input, m45, idx = 1L); rb45 <- studio_run(h45, rv45$exposure_label, rv45$outcome_label, rv45$input, m45, idx = 2L)
ra45b <- pair_input("expA", "outA"); rc45 <- pair_input("outA", "expA")
rA45 <- studio_run(h45, ra45b$exposure_label, ra45b$outcome_label, ra45b$input, m45, idx = 3L); rB45 <- studio_run(h45, rc45$exposure_label, rc45$outcome_label, rc45$input, m45, idx = 4L)
pr45 <- function(ex, out, r, w) list(ex = ex, out = out, rev = w, exName = if (startsWith(ex, "exp")) DEMO_EXPOSURES$name[DEMO_EXPOSURES$id == ex] else DEMO_OUTCOMES$name[DEMO_OUTCOMES$id == ex], outName = if (startsWith(out, "exp")) DEMO_EXPOSURES$name[DEMO_EXPOSURES$id == out] else DEMO_OUTCOMES$name[DEMO_OUTCOMES$id == out], folder = r$pair_folder, methods = as.list(m45), ci = 0.95, pthr = 5e-8, simulated = TRUE)
ui_state_set(h45$con, "pairs", list(pr45("expB", "outB", ra45, FALSE), pr45("outB", "expB", rb45, TRUE), pr45("expA", "outA", rA45, FALSE), pr45("outA", "expA", rB45, TRUE)))
c45 <- lapply(ui_state_get(h45$con, "pairs"), function(p) pair_context(h45, p)); t45 <- report_text(c45, "Both ways")
expect("Round 45: forward and reverse pairs are matched", identical(bidir_pairs(c45), list(c(1L, 2L), c(3L, 4L))))
expect("Round 45: the reverse pair recovers the simulated reverse effect (0.10 for B) and none for A",
       { b <- c45[[2]]$ce; a <- c45[[4]]$ce; abs(b$estimate[b$method == "ivw_mre"] - 0.10) < 0.05 && b$p[b$method == "ivw_mre"] < 0.05 && a$p[a$method == "ivw_mre"] >= 0.05 })
expect("Round 45: Results say both directions for B and one direction for A", any(grepl("^Bidirectional analysis \\(Simulated exposure B .*Both directions gave P below 0.05", t45$results)) && any(grepl("^Bidirectional analysis \\(Simulated exposure A .*supports an effect of Simulated exposure A on Simulated outcome A but not", t45$results)))
expect("Round 45: Methods and Limitations describe bidirectional MR once", sum(t45$methods == BIDIR_METHODS_TEXT) == 1L && sum(t45$limitations == BIDIR_LIMITS_TEXT) == 1L)
expect("Round 45: verdict covers all four cases", identical(vapply(list(c(0.01, 0.5), c(0.5, 0.01), c(0.01, 0.01), c(0.5, 0.5)), function(x) bidir_verdict(x[1], x[2]), ""), c("forward", "reverse", "both", "neither")))
# ---- Round 46: PDF report ----
pk46 <- build_package(h45, c("report", "figs"))
pf46 <- file.path(h45$dir, grep("report\\.pdf$", pk46$files, value = TRUE))
expect("Round 46: the package has a PDF report next to the Word and HTML ones", length(pf46) == 1L && file.exists(pf46) && file.size(pf46) > 20000 && identical(readBin(pf46, "raw", 5), charToRaw("%PDF-")))
ff46 <- tempfile(fileext = ".pdf"); np46 <- pdf_write(ff46, list(list(type = "title", text = "T"), list(type = "p", text = paste(rep("word", 2000), collapse = " ")), list(type = "table", data = data.frame(a = 1:60, b = "x"))))
expect("Round 46: long text and tables flow onto further pages without error", file.exists(ff46) && np46 >= 3L && file.size(ff46) > 5000)

# ---- Round 47: radial plot and leave-one-out before and after outlier removal (FR-98) ----
lg47 <- outlier_legends(rec42, list(p = list(exName = "X", outName = "Y")), 5L, ag44, loo = TRUE)
expect("Round 47: legends are named and numbered in order: comparison, agreement, radial, leave-one-out", identical(names(lg47), c("comparison", "agreement", "radial", "loo")) && all(mapply(function(k, x) grepl(sprintf("^Supplementary Figure S%d\\.", k), x), 5:8, lg47)) && grepl("20% of them", lg47[["comparison"]], fixed = TRUE))
f47 <- tempfile(fileext = ".pdf"); ok47 <- tryCatch({ grDevices::pdf(f47, 7, 5); draw_outlier_radial(p42, rec42, "beta"); lf47 <- data.frame(loo_b = rnorm(30, 0.3, 0.01)); draw_outlier_loo(lf47, lf47[1:25, , drop = FALSE], 0.3, 0.29, "beta"); grDevices::dev.off(); TRUE }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); FALSE })
expect("Round 47: radial and before-and-after figures draw without error", ok47 && file.size(f47) > 1000)
rr47 <- studio_run_reduced(h42, inp42, out42, OUTLIER_REANALYSIS, 0.95, 20L, params = list(ci_level = 0.95, nboot = 200L), keep_loo = TRUE)
expect("Round 47: an outlier set can carry its own leave-one-out analysis (one row per remaining variant)", { lf <- loo_facts(list(pay = pair_payload(h42$dir, rr47$pair_folder))); !is.null(lf) && nrow(lf) == 37L })

# ---- Round 48: instrument strength per outlier set (FR-95c) ----
d48 <- d42$data[d42$data$mr_keep, ]
expect("Round 48: mean F and I2GX are recomputed for every set from the remaining variants", all(c("mean_f", "i2gx", "strength_flag") %in% names(sc42)) &&
       abs(sc42$mean_f[1] - mean((d48$beta_x / d48$se_x)^2)) < 1e-9 && { k <- !(d48$variant_id %in% sc42$ids[[which(sc42$set == "radial")]]); abs(sc42$i2gx[sc42$set == "radial"] - TwoSampleMR::Isq(abs(d48$beta_x[k]), d48$se_x[k])) < 1e-12 })
sw48 <- sc42; sw48$mean_f[2] <- 8; sw48$i2gx[3] <- 0.5; sw48 <- outlier_strength(sw48[, setdiff(names(sw48), c("mean_f", "i2gx", "strength_flag"))], p42)
expect("Round 48: a set below the R23 or R24 threshold is flagged but stays eligible", { f <- outlier_strength(transform(sc42[, setdiff(names(sc42), c("mean_f", "i2gx", "strength_flag"))]), list(variants = lapply(p42$variants, function(v) { v$sx <- as.numeric(v$sx) * 5; v }))); all(grepl("mean F below 10", f$strength_flag)) && identical(f$eligible, sc42$eligible) })

# ---- Round 49: genes and other traits of flagged variants (FR-106, FR-112); folder layout (FR-101) ----
expect("Round 49: the module folder has the FR-101 layout with a README, Excel comparison and the scenario index", all(file.exists(file.path(h42$dir, rel42, c("README.txt", "01_flagged/outliers_flagged.csv", "02_scenarios/scenarios.csv", "03_tables/outlier_comparison.csv", "03_tables/outlier_comparison.xlsx")))))
rs49 <- rec42; rs49$scenarios <- lapply(rs49$scenarios, function(x) { x$ids <- as.list(sub("^demo:([0-9]+):A:G$", "rs\\1", unlist(x$ids))); x })
fake49 <- list(genes = function(rs) data.frame(rsid = rs, chr = "1", pos = 100, gene = paste0("G", seq_along(rs)), stringsAsFactors = FALSE),
               phewas = function(snps, p, jwt) data.frame(rsid = rep(snps[1], 4), id = c("own-1", "t-1", "t-2", "t-2"), trait = c("Own GWAS", "Height", "Waist", "Waist"), p = c(1e-20, 1e-10, 1e-9, 1e-12), stringsAsFactors = FALSE))
lk49 <- outlier_lookup(rs49, "jwt", fake49, exclude = "own-1"); f49 <- lk49$rows[lk49$rows$variant_id == outlier_flagged_ids(rs49)[1], ]
expect("Round 49: lookup gives nearest gene and other traits, excluding the pair's own GWAS and counting each trait once", nrow(lk49$rows) == length(outlier_flagged_ids(rs49)) && f49$n_other_traits == 2L && grepl("^Waist \\(P = 1.00e-12\\); Height", f49$top_other_traits) && !grepl("Own GWAS", f49$top_other_traits))
expect("Round 49: without a token genes are still found and the note says traits need a token", { l <- outlier_lookup(rs49, NULL, fake49); all(!is.na(l$rows$nearest_gene)) && all(is.na(l$rows$n_other_traits)) && any(grepl("OpenGWAS token", l$notes)) })
expect("Round 49: simulated variants without rsIDs give a clear message", grepl("no rsIDs", err_of(outlier_lookup(rec42, "jwt", fake49))))
outlier_write_lookup(h42$dir, r42$pair_folder, lk49); fl49 <- utils::read.csv(file.path(h42$dir, rel42, "01_flagged", "outliers_flagged.csv"))
expect("Round 49: the lookup is saved on its own and merged into the flag ledger", file.exists(file.path(h42$dir, rel42, "01_flagged", "flagged_annotations.csv")) && all(c("nearest_gene", "n_other_traits", "top_other_traits") %in% names(fl49)))
rl49 <- rs49; rl49$lookup <- list(rows = lapply(seq_len(nrow(lk49$rows)), function(i) as.list(lk49$rows[i, ])), date = "2026-11-18")
expect("Round 49: the report lists each flagged variant with its gene and number of other traits", grepl("Flagged variants, with nearest gene and number of other traits", outlier_lookup_text(rl49)) && grepl("G1; 2 other traits", outlier_lookup_text(rl49), fixed = TRUE))

# ---- Round 50: multivariable MR ----
h50 <- studio_new("Multivariable"); p50 <- mv_pair_input(c("expA", "expB"), "outA"); m50 <- c("mvivw", "mvegger", "mvmedian", "mvlasso")
r50 <- studio_run(h50, p50$exposure_label, p50$outcome_label, p50$input, m50, idx = 1L, params = list(ci_level = 0.95, mv_iter = 1000L))
pl50 <- pair_payload(h50$dir, r50$pair_folder); e50 <- .est_df(pl50); d50 <- e50[e50$estimand == "direct_effect", ]
expect("Round 50: the four MVMR methods run as scripts and replay identically", all(vapply(r50$results, function(x) identical(x$status, "ok"), TRUE)) && identical(pl50$replay$status, "identical"))
expect("Round 50: one direct effect per exposure and method, plus the MVMR-Egger intercept", nrow(d50) == 8L && setequal(d50$exposure_id, c("demo:expA", "demo:expB")) && sum(e50$estimand == "pleiotropy_intercept") == 1L)
dd50 <- p50$input$data; mi50 <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = as.matrix(dd50[, c("beta_x_1", "beta_x_2")]), bxse = as.matrix(dd50[, c("se_x_1", "se_x_2")]), by = dd50$beta_y, byse = dd50$se_y), nx = c(DEMO_EXPOSURES$n[1], DEMO_EXPOSURES$n[2]))
expect("Round 50: MVMR-IVW equals MendelianRandomization::mr_mvivw called directly", max(abs(d50$estimate[d50$method == "mvivw"] - mi50@Estimate)) < 1e-12)
expect("Round 52: the conditional F statistics and Q_A follow Sanderson et al. 2021 (zero covariance), without needing sample sizes", {
  bx <- as.matrix(dd50[, c("beta_x_1", "beta_x_2")]); sx <- as.matrix(dd50[, c("se_x_1", "se_x_2")]); J <- nrow(bx)
  f1 <- { d <- coef(lm(bx[, 1] ~ 0 + bx[, 2])); sum((bx[, 1] - d * bx[, 2])^2 / (sx[, 1]^2 + d^2 * sx[, 2]^2)) / (J - 1) }
  qa <- sum((dd50$beta_y - bx %*% mi50@Estimate)^2 / (dd50$se_y^2 + sx^2 %*% mi50@Estimate^2))
  dg <- Filter(function(g) g$method == "mvivw", pl50$diagnostics); cf <- vapply(Filter(function(g) g$name == "mv_cond_f", dg), function(g) g$value, 0)
  length(cf) == 2L && abs(cf[1] - f1) < 1e-8 && abs(Filter(function(g) g$name == "mv_q_a", dg)[[1]]$value - qa) < 1e-8 })
r53 <- studio_run(h50, "MV correlated", p50$outcome_label, p50$input, "mvivw", idx = 2L, params = list(ci_level = 0.95, mv_r_1_2 = 0.3)); pl53 <- pair_payload(h50$dir, r53$pair_folder)
expect("Round 53: with a phenotypic correlation the conditional F and Q_A use covariance r x se1 x se2, as MVMR::phenocov_mvmr", {
  bx <- as.matrix(dd50[, c("beta_x_1", "beta_x_2")]); sx <- as.matrix(dd50[, c("se_x_1", "se_x_2")]); J <- nrow(bx); r <- 0.3
  d <- coef(lm(bx[, 1] ~ 0 + bx[, 2])); f1 <- sum((bx[, 1] - d * bx[, 2])^2 / (sx[, 1]^2 + d^2 * sx[, 2]^2 - 2 * d * r * sx[, 1] * sx[, 2])) / (J - 1)
  b <- mi50@Estimate; qa <- sum((dd50$beta_y - bx %*% b)^2 / (dd50$se_y^2 + b[1]^2 * sx[, 1]^2 + b[2]^2 * sx[, 2]^2 + 2 * b[1] * b[2] * r * sx[, 1] * sx[, 2]))
  dg <- Filter(function(g) g$method == "mvivw", pl53$diagnostics)
  identical(pl53$replay$status, "identical") && abs(Filter(function(g) g$name == "mv_cond_f", dg)[[1]]$value - f1) < 1e-8 && abs(Filter(function(g) g$name == "mv_q_a", dg)[[1]]$value - qa) < 1e-8 })
expect("Round 53: the correlation is read back from the script, and an invalid correlation matrix is blocked with the reason", {
  st <- pair_settings(h50$dir, r53$pair_folder); v <- ADAPTER_MVIVW_V11$check(mv_pair_input(c("expA", "expB", "expC"), "outA")$input, list(mv_r_1_2 = 0.9, mv_r_1_3 = 0.9, mv_r_2_3 = -0.9))
  isTRUE(all.equal(st$mv_r_1_2, 0.3)) && st$mv_r_1_3 == 0 && identical(v$status, "block") && identical(v$codes, "FEAS_INVALID_CORRELATION") })
expect("Round 53: the methods text states the covariance assumption, with the correlation when one is set", grepl("zero covariance", .mv_cov_text(list(), c("A", "B"))) &&
  grepl("A and B r = 0.3", .mv_cov_text(list(mv_r_1_2 = 0.3), c("A", "B")), fixed = TRUE))
r54 <- studio_run(h50, "MV Q-min", p50$outcome_label, p50$input, "mvqhet", idx = 3L, params = list(ci_level = 0.95, mv_r_1_2 = 0.3)); pl54 <- pair_payload(h50$dir, r54$pair_folder); e54 <- .est_df(pl54)
expect("Round 54: MVMR Q-minimisation runs as a script, replays identically and gives one direct effect per exposure with a jackknife SE", identical(pl54$replay$status, "identical") &&
  nrow(e54) == 2L && all(e54$estimand == "direct_effect") && all(is.finite(e54$se) & e54$se > 0) && all(e54$ci_method == "wald"))
expect("Round 54: the Q-minimisation estimate minimises Q_A (lower than at MVMR-IVW) and Q_A at the estimate is at most its degrees of freedom once tau2 is added", {
  bx <- as.matrix(dd50[, c("beta_x_1", "beta_x_2")]); sx <- as.matrix(dd50[, c("se_x_1", "se_x_2")]); r <- 0.3
  qa <- function(b, t2 = 0) sum((dd50$beta_y - bx %*% b)^2 / (dd50$se_y^2 + b[1]^2 * sx[, 1]^2 + b[2]^2 * sx[, 2]^2 + 2 * b[1] * b[2] * r * sx[, 1] * sx[, 2] + t2))
  dg <- Filter(function(g) g$method == "mvqhet", pl54$diagnostics); t2 <- Filter(function(g) g$name == "mv_qhet_tau2", dg)[[1]]$value; b <- e54$estimate
  qa(b, t2) <= qa(mi50@Estimate, t2) + 1e-9 && t2 >= 0 && qa(b, t2) <= nrow(bx) - 2 + 1e-6 })
expect("Round 55: the Q-minimisation script pins boot by its CRAN version (1.3-32, the conda-forge build for R 4.5), so the re-run helper can install it", {
  h <- readLines(Sys.glob(file.path(h50$dir, "03_pairs", r54$pair_folder, "02_scripts", "*mvqhet*.R"))[1], n = 12)
  any(grepl("^# Package  : MendelianRandomization version 0\\.10\\.0; boot version 1\\.3-32$", h)) })
r54b <- studio_run(h50, "MV Q-min boot", p50$outcome_label, p50$input, "mvqhet", idx = 4L, params = list(ci_level = 0.95, mv_qhet_ci = "bootstrap", mv_qhet_iter = 60L)); e54b <- .est_df(pair_payload(h50$dir, r54b$pair_folder))
expect("Round 54: the bootstrap option gives BCa intervals around the same estimate, and the methods text names the method and its uncertainty", all(e54b$ci_method == "bootstrap") &&
  all(e54b$ci_lower < e54b$estimate & e54b$estimate < e54b$ci_upper) && isTRUE(all.equal(e54b$estimate, .est_df(pair_payload(h50$dir, studio_run(h50, "MV Q-min 0", p50$outcome_label, p50$input, "mvqhet", idx = 5L, params = list(ci_level = 0.95))$pair_folder))$estimate)) && {
  tx <- report_text(list(mv_context(h50, list(folder = r54b$pair_folder, exNames = list("A", "B"), methods = list("mvqhet"), outName = "Y"))), "Q")
  any(grepl("Q-statistic minimisation", tx$methods)) && any(grepl("bootstrap intervals \\(60 iterations, boot 1\\.3-32\\)", tx$methods)) })
expect("Round 55: mv_snpcov gives the covariance of single-variant slopes whose diagonal is exactly lm()'s squared SE, and the file reader skips non-dosage columns and reads PLINK --recode A files", {
  set.seed(5); n <- 400; G <- sapply(1:3, function(j) rbinom(n, 2, 0.3)); colnames(G) <- c("rs1", "rs2", "rs3"); X <- cbind(G %*% c(.2, .1, 0) + rnorm(n), G %*% c(0, .1, .2) + rnorm(n))
  cv <- mv_snpcov(G, X); se <- summary(lm(X[, 2] ~ G[, 3]))$coefficients[2, 2]
  e2 <- resid(lm(X[, 1] ~ G[, 2])); f2 <- resid(lm(X[, 2] ~ G[, 2])); c12 <- sum(e2 * f2) / ((n - 2) * sum((G[, 2] - mean(G[, 2]))^2))
  f <- tempfile(fileext = ".csv"); utils::write.csv(data.frame(pid = 1001:1400, x1 = X[, 1], x2 = X[, 2], G), f, row.names = FALSE); fc <- mv_cov_from_file(f, c("x1", "x2"))
  fr <- tempfile(fileext = ".raw"); utils::write.table(data.frame(FID = 1:n, IID = 1:n, PAT = 0, MAT = 0, SEX = 1, PHENOTYPE = -9, x1 = X[, 1], x2 = X[, 2], rs1_A = G[, 1], rs2_G = G[, 2], rs3_T = G[, 3]), fr, sep = " ", row.names = FALSE, quote = FALSE)
  fp <- mv_cov_from_file(fr, c("x1", "x2"))
  identical(fp$variant_id, c("rs1", "rs2", "rs3")) && isTRUE(all.equal(fp$cov_x_1_2, cv$cov_x_1_2)) && setequal(attr(fp, "meta")$ignored, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")) &&
  abs(cv$cov_x_2_2[3] - se^2) < 1e-15 && abs(cv$cov_x_1_2[2] - c12) < 1e-15 && identical(fc$variant_id, c("rs1", "rs2", "rs3")) && identical(attr(fc, "meta")$ignored, "pid") && attr(fc, "meta")$people == 400L })
d55 <- p50$input$data; cv55 <- data.frame(variant_id = d55$variant_id, cov_x_1_1 = d55$se_x_1^2, cov_x_1_2 = 0.3 * d55$se_x_1 * d55$se_x_2, cov_x_2_2 = d55$se_x_2^2)
mv_cov_save(h50, c("expA", "expB"), cv55, "cohort.csv"); i55 <- mv_with_cov(h50, c("expA", "expB"), p50$input)
r55 <- studio_run(h50, "MV indiv", p50$outcome_label, i55, c("mvivw", "mvqhet"), idx = 6L, params = list(ci_level = 0.95)); pl55 <- pair_payload(h50$dir, r55$pair_folder)
expect("Round 55: a covariance table equal to 0.3 x se1 x se2 gives the same conditional F, Q_A and Q-minimisation as the phenotypic correlation 0.3, and is recorded as individual-level data", {
  v <- function(pl, m, n) vapply(Filter(function(g) g$method == m && g$name == n, pl$diagnostics), function(g) as.numeric(g$value), 0)
  identical(pl55$replay$status, "identical") && max(abs(v(pl55, "mvivw", "mv_cond_f") - v(pl53, "mvivw", "mv_cond_f"))) < 1e-10 && abs(v(pl55, "mvivw", "mv_q_a") - v(pl53, "mvivw", "mv_q_a")) < 1e-10 &&
    max(abs(.est_df(pl55)$estimate[.est_df(pl55)$method == "mvqhet"] - e54$estimate)) < 1e-8 && identical(Filter(function(g) g$name == "mv_pheno_cor" && g$method == "mvivw", pl55$diagnostics)[[1]]$note, "individual-level data") })
expect("Round 55: the covariance table is saved in 02_datasets and registered, a model with uncovered instruments is refused with the count, and clearing it returns to the correlations", {
  e <- mv_cov_get(h50, c("expA", "expB")); ok1 <- length(e) == 1L && file.exists(file.path(h50$dir, e[[1]]$file)) && startsWith(e[[1]]$file, "02_datasets/") && any(manifest_read(h50$dir)$relative_path == e[[1]]$file)
  bad <- tryCatch({ mv_attach_cov(p50$input, cv55[-(1:2), ]); "" }, error = function(x) conditionMessage(x))
  ok1 && grepl("^2 of the [0-9]+ instruments have no covariance", bad) && length(mv_cov_clear(h50, c("expA", "expB"))) == 0L && identical(mv_with_cov(h50, c("expA", "expB"), p50$input), p50$input) })
expect("Round 55: the methods text says the covariances came from individual-level data", grepl("computed from individual-level data", .mv_cov_text(list(), c("A", "B"), "individual-level data")))
# ---- Round 57: mediation ----
h57 <- studio_new("Mediation"); ms57 <- c("expC", "expD"); st57 <- c(list(list(role = "total")), lapply(ms57, function(m) list(role = "xm", m = m)), lapply(ms57, function(m) list(role = "direct", m = m)), list(list(role = "joint")))
pr57 <- lapply(seq_along(st57), function(k) { med <- c(st57[[k]], list(x = "expA", ms = ms57, y = "outA")); p <- mediation_pair_input(med); mv <- p$input$kind == "multi_iv"
  r <- studio_run(h57, p$exposure_label, p$outcome_label, p$input, if (mv) "mvivw" else c("ivw_mre", "egger"), idx = k, params = list(ci_level = 0.95))
  ex <- switch(med$role, total = "expA", xm = "expA", direct = c("expA", med$m), joint = c("expA", ms57)); out <- if (med$role == "xm") med$m else "outA"
  list(ex = ex, out = out, med = med, folder = r$pair_folder, exName = p$exposure_label, outName = p$outcome_label, exNames = as.list(p$exposure_labels), mv = mv, methods = as.list(if (mv) "mvivw" else c("ivw_mre", "egger")), replay = pair_payload(h57$dir, r$pair_folder)$replay$status) })
res57 <- mediation_run(h57, pr57); w57 <- demo_mediation_world("expA", ms57, "outA"); g57 <- function(m, me, q) res57[res57$mediator == m & res57$method == me & res57$quantity == q, ]
expect("Round 57: every mediation step runs as an ordinary pair and replays identically", all(vapply(pr57, function(p) identical(p$replay, "identical"), TRUE)) && length(pr57) == 6L)
expect("Round 57: the simulated mediation is recovered: indirect effect through C covers 0.4 x 0.2, through D covers 0, and the total is within 3 SE of 0.38 (one simulated draw)", {
  ic <- g57("Simulated exposure C", "product", "indirect"); id <- g57("Simulated exposure D", "product", "indirect"); tt <- res57[res57$quantity == "total", ]
  ic$ci_lower <= 0.08 && 0.08 <= ic$ci_upper && id$ci_lower <= 0 && 0 <= id$ci_upper && abs(tt$estimate - w57$truth$total) < 3 * tt$se })
expect("Round 57: the product, difference and proportion follow their formulas from the pair tables (delta method)", {
  e <- function(f, m) .est_df(pair_payload(h57$dir, f))[.est_df(pair_payload(h57$dir, f))$method == m, ]; tt <- e(pr57[[1]]$folder, "ivw_mre"); a <- e(pr57[[2]]$folder, "ivw_mre"); mv <- e(pr57[[4]]$folder, "mvivw")
  b <- mv[mv$exposure_id == "demo:expC", ]; d <- mv[mv$exposure_id == "demo:expA", ]; ip <- g57("Simulated exposure C", "product", "indirect"); pp <- g57("Simulated exposure C", "product", "proportion mediated"); pd <- g57("Simulated exposure C", "difference", "proportion mediated")
  abs(ip$estimate - a$estimate * b$estimate) < 1e-12 && abs(ip$se - sqrt(a$estimate^2 * b$se^2 + b$estimate^2 * a$se^2)) < 1e-12 && abs(pp$estimate - a$estimate * b$estimate / tt$estimate) < 1e-12 &&
    abs(pd$estimate - (1 - d$estimate / tt$estimate)) < 1e-12 && any(res57$mediator == "(all together)") })
expect("Round 57: the mediation script is written to 04_across_pairs, registered, and re-running it from the study folder gives the same table", {
  f <- file.path(h57$dir, "04_across_pairs", "mediation", "mediation.R"); old <- setwd(h57$dir); again <- sys.source(f, envir = new.env()); setwd(old)
  again <- utils::read.csv(file.path(h57$dir, "04_across_pairs", "mediation", "mediation_results.csv"), stringsAsFactors = FALSE)
  file.exists(f) && isTRUE(all.equal(again, res57, check.attributes = FALSE)) && all(c("04_across_pairs/mediation/mediation.R", "04_across_pairs/mediation/mediation_results.csv") %in% manifest_read(h57$dir)$relative_path) })
expect("Round 76: each mediation row carries its own scale (the exposure on a mediator in its units, proportions as proportions; test 3 labelled them log OR)", all(res57$scale[res57$quantity == "proportion mediated"] == "proportion") && all(res57$scale[res57$quantity != "proportion mediated"] == "beta"))
expect("Round 57: the check screen can describe every step, a mediator equal to the exposure is refused, and the results text reports the proportions", {
  f1 <- pair_facts("expA", "expC", 5e-8, mediation_pair_input(pr57[[2]]$med)); f2 <- mv_pair_facts(c("expA", "expC"), "outA", 5e-8, mediation_pair_input(pr57[[4]]$med))
  bad <- tryCatch({ mediation_pair_input(list(role = "total", x = "expA", ms = "expA", y = "outA")); "" }, error = function(e) conditionMessage(e)); tx <- mediation_results_text(res57, attr(res57, "spec"))
  f1$instruments == 50L && f2$instruments == 100L && grepl("must differ", bad) && any(grepl("^Through Simulated exposure C: .*proportion mediated [0-9]+%", tx)) && any(grepl("^All mediators together", tx)) })
expect("Round 57: the publication package text has the mediation methods, results, limitations and the Carter references", {
  ui_state_set(h57$con, "pairs", lapply(pr57, function(p) p[setdiff(names(p), "replay")])); pk57 <- build_package(h57, c("report", "text", "refs"), 140)
  hs <- paste(readLines(file.path(pk57$dir, "report.html"), warn = FALSE), collapse = "\n")
  grepl("Mediation was assessed with Mendelian randomization", hs) && grepl("Through Simulated exposure C", hs) && grepl("assume no interaction between the exposure and the mediator", hs) && grepl("Understanding the consequences of education inequality", hs) && grepl("Mendelian randomisation for mediation analysis", hs) })
# ---- Round 58: more MendelianRandomization estimators ----
h58 <- studio_new("MR methods"); p58 <- pair_input("expA", "outA"); m58 <- c("maxlik", "mbe", "cml", "robust_ivw", "pivw")
r58 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, m58, idx = 1L, params = list(ci_level = 0.95, mbe_iter = 1000L, cml_pert = 50L)); pl58 <- pair_payload(h58$dir, r58$pair_folder); e58 <- .est_df(pl58)
expect("Round 58: maximum likelihood, MBE, MR-cML, MR-Robust and pIVW run as scripts, replay identically and cover the simulated effect 0.30", identical(pl58$replay$status, "identical") &&
  all(vapply(r58$results, function(x) identical(x$status, "ok"), TRUE)) && setequal(e58$method, m58) && all(e58$ci_lower <= 0.3 & 0.3 <= e58$ci_upper))
expect("Round 58: maximum likelihood and MR-Robust equal MendelianRandomization called directly", { d <- p58$input$data; mi <- MendelianRandomization::mr_input(d$beta_x, d$se_x, d$beta_y, d$se_y)
  abs(e58$estimate[e58$method == "maxlik"] - MendelianRandomization::mr_maxlik(mi)@Estimate) < 1e-10 && abs(e58$estimate[e58$method == "robust_ivw"] - MendelianRandomization::mr_ivw(mi, robust = TRUE, penalized = TRUE)@Estimate) < 1e-10 })
expect("Round 58: MR-cML and MVMR-cML are blocked with the reason when sample sizes are missing", { q <- p58$input; q$data$n_x <- NA; m <- mv_pair_input(c("expA", "expB"), "outA")$input; m$data$n_y <- NULL
  identical(ADAPTER_CML_V11$check(q, resolve_params(ADAPTER_CML_V11, list()))$codes, "FEAS_NEEDS_N") && identical(ADAPTER_MVCML_V11$check(m, resolve_params(ADAPTER_MVCML_V11, list()))$codes, "FEAS_NEEDS_N") })
# ---- Round 59: TwoSampleMR estimators and diagnostics ----
m59 <- c("egger_boot", "pwm", "smode_nome", "wmode_nome", "rucker", "sign")
r59 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, m59, idx = 2L, params = list(ci_level = 0.95, eb_nboot = 200L, pwm_nboot = 200L, smn_nboot = 200L, wmn_nboot = 200L)); pl59 <- pair_payload(h58$dir, r59$pair_folder); e59 <- .est_df(pl59)
expect("Round 59: bootstrap Egger, penalised weighted median, NOME modes and Rucker run as scripts, replay identically and cover the simulated effect 0.30", identical(pl59$replay$status, "identical") &&
  all(vapply(r59$results, function(x) identical(x$status, "ok"), TRUE)) && { ce <- e59[e59$estimand == "causal_effect", ]; setequal(ce$method, setdiff(m59, "sign")) && all(ce$ci_lower <= 0.3 & 0.3 <= ce$ci_upper) })
expect("Round 59: the sign test gives a diagnostic only, Rucker reports Q, Q' and the selected model, and bootstrap Egger adds its intercept", {
  g <- function(n) Filter(function(x) x$name == n, pl59$diagnostics); !any(e59$method == "sign") && length(g("sign_concordance")) == 1L && all(c("q_ivw", "q_rucker", "q_diff", "rucker_model") %in% vapply(pl59$diagnostics, function(x) x$name, "")) &&
    any(e59$method == "egger_boot" & e59$estimand == "pleiotropy_intercept") })
expect("Round 59: the funnel figure has a precision-weighted density panel and its legend says so", { f <- tempfile(fileext = ".png"); lf <- data.frame(wr = rnorm(20, 0.3, 0.05), wr_se = runif(20, 0.02, 0.06))
  grDevices::png(f); .draw_funnel(lf, list()); grDevices::dev.off(); file.exists(f) && grepl("Right: density", paste(deparse(.report_text_uni), collapse = " ")) })
# ---- Round 60: MRMix ----
r60 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, "mrmix", idx = 3L, params = list(ci_level = 0.95)); pl60 <- pair_payload(h58$dir, r60$pair_folder); e60 <- .est_df(pl60)
expect("Round 60: MRMix runs as a script, replays identically, covers the simulated effect 0.30 and keeps its standardised estimate", identical(pl60$replay$status, "identical") && nrow(e60) == 1L &&
  e60$ci_lower <= 0.3 && 0.3 <= e60$ci_upper && any(vapply(pl60$diagnostics, function(g) identical(g$name, "mrmix_theta_std"), TRUE)))
expect("Round 60: MRMix is blocked with the reason without allele frequencies or sample sizes", { q <- p58$input; q$data$eaf_x <- NA; q$data$eaf_y <- NA; q2 <- p58$input; q2$data$n_y <- NA
  identical(ADAPTER_MRMIX_V11$check(q, list())$codes, "FEAS_NEEDS_EAF") && identical(ADAPTER_MRMIX_V11$check(q2, list())$codes, "FEAS_NEEDS_N") })
# ---- Round 61: MR-GRIP ----
r61 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, "grip", idx = 4L, params = list(ci_level = 0.95)); pl61 <- pair_payload(h58$dir, r61$pair_folder); e61 <- .est_df(pl61)
expect("Round 61: MR-GRIP runs as a script, replays identically, equals TwoSampleMR::mr_grip and is unchanged when an allele is recoded", {
  d <- p58$input$data; g <- TwoSampleMR::mr_grip(d$beta_x, d$beta_y, d$se_x, d$se_y, TwoSampleMR::default_parameters()); d2 <- d; d2$beta_x[1:5] <- -d2$beta_x[1:5]; d2$beta_y[1:5] <- -d2$beta_y[1:5]
  g2 <- TwoSampleMR::mr_grip(d2$beta_x, d2$beta_y, d2$se_x, d2$se_y, TwoSampleMR::default_parameters()); ce <- e61[e61$estimand == "causal_effect", ]
  identical(pl61$replay$status, "identical") && abs(ce$estimate - g$b) < 1e-10 && abs(g2$b - g$b) < 1e-12 && any(e61$estimand == "pleiotropy_intercept") })
# ---- Round 68b: IVW after Steiger filtering ----
in68b <- p58$input; uv <- which(in68b$data$variant_id %in% .usable(in68b$data)$variant_id)[1:3]
in68b$data$beta_y[uv] <- 6 * in68b$data$beta_x[uv]   # three variants that explain more variance in the outcome than in the exposure
r68b <- studio_run(h58, p58$exposure_label, p58$outcome_label, in68b, "ivw_steiger", idx = 6L, params = list(ci_level = 0.95)); pl68b <- pair_payload(h58$dir, r68b$pair_folder); e68b <- .est_df(pl68b)
expect("Round 68b: IVW after Steiger filtering equals TwoSampleMR::steiger_filtering then mr_ivw (3 variants removed), replays identically and lists the removed variants", {
  d <- .usable(in68b$data)
  td <- data.frame(SNP = d$variant_id, beta.exposure = d$beta_x, se.exposure = d$se_x, samplesize.exposure = d$n_x, pval.exposure = 2 * stats::pnorm(-abs(d$beta_x / d$se_x)),
                   beta.outcome = d$beta_y, se.outcome = d$se_y, samplesize.outcome = d$n_y, pval.outcome = 2 * stats::pnorm(-abs(d$beta_y / d$se_y)), units.exposure = NA, units.outcome = NA,   # no units: r2 from beta, se and n, as the app's Steiger test
                   id.exposure = "x", id.outcome = "y", exposure = "x", outcome = "y", mr_keep = TRUE, stringsAsFactors = FALSE)
  sf <- suppressMessages(TwoSampleMR::steiger_filtering(td)); k <- sf[sf$steiger_dir, ]; f <- TwoSampleMR::mr_ivw(k$beta.exposure, k$beta.outcome, k$se.exposure, k$se.outcome)
  dg <- pl68b$diagnostics; gv <- function(n) dg[vapply(dg, function(g) identical(g$name, n), TRUE)][[1]]
  identical(pl68b$replay$status, "identical") && abs(e68b$estimate - f$b) < 1e-10 && abs(e68b$se - f$se) < 1e-10 && gv("steiger_kept")$value == nrow(k) &&
    gv("steiger_removed")$value == 3 && setequal(strsplit(gv("steiger_removed")$note, " ")[[1]], sf$SNP[!sf$steiger_dir]) })
# ---- Round 62, Round 70e: BWMR, MR Studio's own implementation of the published algorithm ----
r62 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, "bwmr", idx = 5L, params = list(ci_level = 0.95)); pl62 <- pair_payload(h58$dir, r62$pair_folder); e62 <- .est_df(pl62)
sc62 <- readLines(Sys.glob(file.path(h58$dir, "03_pairs", r62$pair_folder, "02_scripts", "*_bwmr.R")))
expect("Round 70e: BWMR runs as a script that contains the whole algorithm and needs no package, replays identically and covers the simulated effect 0.30",
       identical(pl62$replay$status, "identical") && nrow(e62) == 1L && e62$ci_lower <= 0.3 && 0.3 <= e62$ci_upper &&
       any(grepl("^bwmr_mrs <- function", sc62)) && !any(grepl("BWMR::", sc62, fixed = TRUE)) && any(grepl("none beyond base R", sc62, fixed = TRUE)))
local({ eval(parse(text = BWMR_MRS_SRC)); set.seed(62); n <- 30; sx <- runif(n, .01, .03); g <- rnorm(n, .08, .04); sy <- runif(n, .01, .04); G <- .3 * g + rnorm(n, 0, sy); G[3] <- G[3] + .4
  f <- bwmr_mrs(g, G, sx, sy)
  expect("Round 70e: BWMR converges, down-weights the planted outlier and gives a variance larger than the mean-field one",
         isTRUE(f$converged) && f$weights[3] < 0.5 && all(f$weights[-3] > 0.5) && f$se > f$se_mfvb) })

# ---- Round 63: bidirectional cML ----
if (requireNamespace("BiDirectCausal", quietly = TRUE)) {
h63 <- studio_new("BiDir cML"); a63 <- pair_input("expB", "outB"); b63 <- pair_input("outB", "expB")
f63 <- studio_run(h63, a63$exposure_label, a63$outcome_label, a63$input, "ivw_mre", idx = 1L)$pair_folder; g63 <- studio_run(h63, b63$exposure_label, b63$outcome_label, b63$input, "ivw_mre", idx = 2L)$pair_folder
pr63 <- list(list(ex = "expB", out = "outB", folder = f63, exName = "B", outName = "outcome B", rev = FALSE), list(ex = "outB", out = "expB", folder = g63, exName = "outcome B", outName = "B", rev = TRUE))
res63 <- bidir_cml_run(h63, pr63)[[1]]
expect("Round 63: BiDirMRcML recovers both simulated directions (B on outcome B 0.25, back 0.10) from the forward and reverse pairs", { g <- function(d) res63[res63$method == "BiDirMRcML" & res63$version == "S.DP" & res63$direction == d, ]
  is.data.frame(res63) && g("XtoY")$ci_lower <= 0.25 && 0.25 <= g("XtoY")$ci_upper && g("YtoX")$ci_lower <= 0.10 && 0.10 <= g("YtoX")$ci_upper })
expect("Round 63: the bidirectional script is registered, re-running it gives the same table, and the text labels BiDirCDcML as standardised", {
  f <- file.path(h63$dir, "04_across_pairs", "bidirectional_cml", paste0(f63, ".R")); old <- setwd(h63$dir); sys.source(f, envir = new.env()); setwd(old)
  again <- utils::read.csv(file.path(h63$dir, "04_across_pairs", "bidirectional_cml", paste0(f63, "_results.csv")))
  isTRUE(all.equal(again, res63, check.attributes = FALSE)) && any(manifest_read(h63$dir)$relative_path == file.path("04_across_pairs", "bidirectional_cml", paste0(f63, ".R"))) && grepl("standardised", bidir_cml_text(res63)) })
} else cat("SKIP  Round 63 (bidirectional cML): package BiDirectCausal is not installed in this runtime\n")
# ---- Round 64: LD scores and LD score regression ----
local({
  set.seed(64); n <- 61L; m <- 40L; g <- matrix(rbinom(n * m, 2, 0.3), n, m); g[, 5] <- g[, 4]; g[3, 7] <- NA; g[, 9] <- 0   # a copy, a missing value, a monomorphic variant
  for (j in 2:m) if (j != 5 && j != 9) { k <- runif(n) < 0.6; g[k, j] <- g[k, j - 1] }   # LD between neighbours
  enc <- function(x) { code <- ifelse(is.na(x), 1L, c(3L, 2L, 0L)[x + 1L]); code <- c(code, rep(0L, (4 - n %% 4) %% 4)); q <- matrix(code, 4); as.raw(q[1, ] + 4L * q[2, ] + 16L * q[3, ] + 64L * q[4, ]) }
  pre <- file.path(tempdir(), "r64"); con <- file(paste0(pre, ".bed"), "wb"); writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con); for (j in 1:m) writeBin(enc(g[, j]), con); close(con)
  bp <- cumsum(sample(20000:60000, m)); writeLines(sprintf("22\trs%d\t0\t%d\tA\tG", 1:m, bp), paste0(pre, ".bim")); writeLines(sprintf("f%d i%d 0 0 0 -9", 1:n, 1:n), paste0(pre, ".fam"))
  d <- ldsc_l2_bfile(pre, wind = 2e5)
  keep <- setdiff(1:m, 9); gi <- g[, keep]; for (j in seq_len(ncol(gi))) gi[is.na(gi[, j]), j] <- mean(gi[, j], na.rm = TRUE)
  r2 <- suppressWarnings(stats::cor(gi))^2; u <- r2 - (1 - r2) / (n - 2); u[abs(outer(bp[keep], bp[keep], "-")) > 2e5] <- 0
  expect("Round 64: LD scores equal the LDSC estimator computed directly (corrected r2 summed within the window; missing set to the mean; monomorphic dropped)",
         nrow(d) == m - 1L && !("rs9" %in% d$SNP) && max(abs(d$L2 - rowSums(u))) < 1e-10)
  expect("Round 64: the decoded genotypes give the allele frequencies of the data", isTRUE(all.equal(d$MAF[1], min(mean(g[, 1]) / 2, 1 - mean(g[, 1]) / 2))))
  hm <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4", "rs5"), A1 = c("A", "C", "A", "G", "T"), A2 = c("G", "T", "C", "A", "G"), stringsAsFactors = FALSE)
  ds <- data.frame(snp = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6"), ea = c("A", "T", "A", "G", "T", "A"), oa = c("G", "C", "G", "A", "G", "C"), beta = c(0.2, 0.3, 0.1, 0.4, 0.5, 0.1), se = 0.1, eaf = c(0.3, 0.4, 0.2, 0.005, 0.5, 0.3), n = c(1000, 1000, 1000, 1000, 500, 1000), stringsAsFactors = FALSE)
  mu <- ldsc_munge(ds, hm); led <- attr(mu, "ledger")
  expect("Round 64: preparation keeps HapMap3 variants with matching alleles, signs Z for the list's A1, drops MAF < 0.01 and low N (munge_sumstats rules)",
         identical(mu$SNP, c("rs1", "rs2")) && isTRUE(all.equal(mu$Z, c(2, -3))) && led$hapmap3 == 5L && led$alleles_match == 4L && led$maf == 3L && led$kept == 2L)
  expect("Round 64: preparation refuses a dataset without rsIDs or without sample size",
         inherits(tryCatch(ldsc_munge(ds[, setdiff(names(ds), "snp")], hm), error = function(e) e), "error") && inherits(tryCatch(ldsc_munge(ds[, setdiff(names(ds), "n")], hm), error = function(e) e), "error"))
  sc <- ldsc_script("p1", c("02_datasets/ldsc/a.sumstats.gz", "02_datasets/ldsc/b.sumstats.gz"), c("A", "B"), "/ld", "official", c(NA, 0.4), c(NA, 0.1))
  expect("Round 64: the LDSC script parses, pins GenomicSEM and writes the liability-scale prevalences", !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_GENOMICSEM, sc)) && any(grepl("c(NA, 0.4)", sc, fixed = TRUE)))
  res <- data.frame(quantity = c("h2", "h2", "rg", "genetic_covariance", "intercept", "intercept", "cross_trait_intercept"), trait = c("A", "B", "both", "both", "A", "B", "both"),
                    estimate = c(0.44, 0.24, -0.17, -0.06, 1.3, 1.07, -0.128), se = c(0.02, 0.007, 0.018, 0.007, 0.057, 0.019, 0.013), p = c(1e-50, 1e-60, 3.7e-22, 1e-15, 1e-7, 1e-4, 2e-24),
                    ci_lower = -0.2, ci_upper = -0.14, scale = "observed", mean_chi2 = NA, snps = 1173569, stringsAsFactors = FALSE); attr(res, "labels") <- c(x = "A", y = "B")
  expect("Round 64: the results text gives heritabilities, genetic correlation and the cross-trait intercept with its reading", grepl("genetic correlation -0.170 (SE 0.018)", ldsc_text(res), fixed = TRUE) && grepl("overlapping samples", ldsc_text(res)))
  dd <- file.path(tempdir(), "ldsc_ref64"); dir.create(dd, showWarnings = FALSE)
  expect("Round 64: with nothing installed the status says so", { st <- ldsc_ref_status(dd); !st$hm3 && !st$official && !length(st$computed) })
})
# ---- Round 65: LHC-MR ----
local({
  ldpos <- data.frame(CHR = c(1, 1, 6, 6, 6), POS = c(1000, 2000, 28.4e6, 30e6, 33.6e6), RSID = paste0("rs", 1:5), stringsAsFactors = FALSE)
  d <- data.frame(snp = c(paste0("rs", 1:5), "rs9"), ea = "A", oa = "G", beta = c(0.02, -0.01, 0.03, 0.5, 0.01, 0.1), se = 0.01, eaf = c(0.3, 0.6, 0.2, 0.4, 0.1, 0.2), n = 1e5, stringsAsFactors = FALSE)
  m <- lhcmr_prepare(d, ldpos); led <- attr(m, "ledger")
  expect("Round 65: LHC-MR preparation keeps variants of the authors' LD file, takes chromosome and position from it and leaves out HLA (chr6 28.5-33.5 Mb)",
         identical(m$SNP, c("rs1", "rs2", "rs3", "rs5")) && identical(m$POS, c(1000, 2000, 28.4e6, 33.6e6)) && led$in_ld_file == 5L && led$hla_removed == 1L && isTRUE(all.equal(m$TSTAT, m$BETA / m$SE)) && isTRUE(all.equal(m$MAF[2], 0.4)))
  expect("Round 65: LHC-MR preparation refuses a dataset without sample size", inherits(tryCatch(lhcmr_prepare(d[, names(d) != "n"], ldpos), error = function(e) e), "error"))
  sc <- lhcmr_script("p1", c("02_datasets/lhcmr/a.tsv.gz", "02_datasets/lhcmr/b.tsv.gz"), c("A", "B"), "/ld.csv", "/rho.csv", "/ldsc", "/hm3", 7L, 1L, bfile = "/panel/EUR", plink = "/bin/plink")
  expect("Round 65: the LHC-MR script parses, pins lhcMR and GenomicSEM, sets the seed and clumps locally", !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_LHCMR, sc)) && any(grepl(PIN_GENOMICSEM, sc)) && any(grepl("set.seed(seed)", sc, fixed = TRUE)) && any(grepl("assignInNamespace(\"clump_data\"", sc, fixed = TRUE)) && any(grepl("mrstudio_local", sc, fixed = TRUE)) && which(grepl("loadNamespace(\"ieugwasr\")", sc, fixed = TRUE))[1] < which(grepl("set.seed(seed)", sc, fixed = TRUE))[1])
  res <- data.frame(parameter = c("axy", "ayx", "tX", "tY", "h2X", "h2Y", "iXY"), estimate = c(-0.128, -0.128, 0.032, -0.124, 0.256, 0.546, -0.103), se = 0.02, p = 1e-5, ci_lower = -0.2, ci_upper = 0, variants = 1e6, variants_used = 1e5, bimodal = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE); attr(res, "labels") <- c(x = "BMI", y = "Height")
  expect("Round 65: the LHC-MR text gives both directions, the confounder effects and the bimodality caution", { t <- lhcmr_text(res); grepl("BMI on Height -0.128", t) && grepl("confounder effects 0.032 on BMI", t) && grepl("bimodal", t) })
  expect("Round 65: with nothing installed the LHC-MR status says so", { st <- lhcmr_ref_status(file.path(tempdir(), "nolhc")); !st$ld && !st$rho })
})
# ---- Round 66: memory-aware heavy analyses ----
local({
  mi <- mem_info()
  expect("Round 66: the computer's total and free memory are read (macOS, Linux and Windows readers)", is.finite(mi$total_gb) && mi$total_gb > 0 && is.finite(mi$available_gb) && mi$available_gb > 0 && mi$available_gb <= mi$total_gb * 1.01)
  fake <- list(total_gb = 16, available_gb = 6)
  expect("Round 66: a run starts only when its peak fits in half of the free memory, and the refusal gives the numbers", mem_check(2.9, "X", fake)$ok && !mem_check(3.1, "X", fake)$ok && grepl("6.0 GB is free now", mem_check(3.1, "X", fake)$msg))
  cal <- list(base_gb = 1, main_gb_per_mvar = 1, per_worker_gb_per_mvar = 1, sp_fits = 50L, sp_minutes = 8, jk_fits = 200L, jk_minutes = 2, prep_hours = 0.2)
  o <- lhcmr_options(1e6, list(total_gb = 16, available_gb = 10), cal, max_workers = 8L); o2 <- lhcmr_options(1e6, list(total_gb = 16, available_gb = 5), cal, max_workers = 8L)
  expect("Round 66: LHC-MR offers only the worker counts whose estimated peak fits safely (none when even one does not)",
         (if (.Platform$OS.type == "windows") identical(o$workers, 1L) else identical(o$workers, c(1L, 2L))) && all(o$memory_gb <= 5) && nrow(o2) == 0L)
  expect("Round 66: the LHC-MR estimate grows with workers and variants, and time falls with workers", { a <- lhcmr_estimate(1e6, 1L, cal); b <- lhcmr_estimate(2e6, 2L, cal); b$memory_gb > a$memory_gb && lhcmr_estimate(1e6, 2L, cal)$hours < a$hours })
  expect("Round 66: LDSC and LD score computation have memory needs that grow with the data", ldsc_memory_gb(1e7) > ldsc_memory_gb(1e6) && LDSCORE_MEMORY_GB > 0)
  bl <- tempfile(fileext = ".log"); rscript_background(c("--vanilla", "-e", shQuote("cat(Sys.getenv('RHOME'), '|', Sys.getenv('MRGUI_T66'))")), log = bl, vars = list(MRGUI_T66 = "MR Studio test"))
  rd <- function(f) tryCatch(readLines(f, warn = FALSE), error = function(e) character())   # on Windows a file still open for writing may refuse to open
  for (w in 1:30) { if (file.exists(bl) && length(rd(bl))) break; Sys.sleep(1) }
  bx <- if (file.exists(bl)) paste(tryCatch(readLines(bl, warn = FALSE), error = function(e) ""), collapse = "") else ""
  ok66b <- grepl(R.home(), bx, fixed = TRUE) && grepl("| MR Studio test", bx, fixed = TRUE) && is.na(Sys.getenv("MRGUI_T66", unset = NA))
  if (!ok66b) cat("  background diagnostics: log exists =", file.exists(bl), "| log =", bx, "| R.home =", R.home(), "| Rscript =", default_rscript(), "\n")
  expect("Round 66: background R jobs get RHOME (needed by a moved runtime) and variables with spaces intact", ok66b)
  if (.Platform$OS.type != "windows") {
    sf <- tempfile(fileext = ".json"); writeLines('{"state":"running"}', sf)
    old <- Sys.getenv("MRSTUDIO_GUARD_FLOOR"); Sys.setenv(MRSTUDIO_GUARD_FLOOR = "2")
    job <- parallel::mcparallel(Sys.sleep(60)); pid <- job$pid
    glog <- tempfile(fileext = ".log"); mem_guard_start(pid, sf, every = 1, engine = if (exists("args_dir")) args_dir else getwd(), log = glog); for (w in 1:60) { if (grepl("stopped", paste(tryCatch(readLines(sf, warn = FALSE), error = function(e) ""), collapse = ""))) break; Sys.sleep(1) };   # up to 60 s: R starts slowly on CI runners
    Sys.setenv(MRSTUDIO_GUARD_FLOOR = old)
    st <- tryCatch(jsonlite::fromJSON(sf), error = function(e) list())
    gone <- !is.null(suppressWarnings(parallel::mccollect(job, wait = TRUE, timeout = 3)))   # a stopped child is collected at once
    if (!gone) { tools::pskill(pid); invisible(suppressWarnings(parallel::mccollect(job))) }
    ok66 <- identical(st$state, "stopped") && grepl("Stopped to protect the computer", st$msg) && gone
    if (!ok66) cat("  guard diagnostics: status =", paste(tryCatch(readLines(sf, warn = FALSE), error = function(e) ""), collapse = " "), "| gone =", gone, "| Rscript =", default_rscript(), "| log:", if (file.exists(glog)) paste(tryCatch(readLines(glog, warn = FALSE), error = function(e) ""), collapse = " / ") else "(none)", "\n")
    expect("Round 66: the memory guard stops the running analysis and writes why when free memory falls below its floor", ok66)
  }
})
# ---- Round 67: adaptive, resumable LHC-MR runner ----
local({
  re <- new.env(); eval(parse(text = LHCMR_RUNNER_SRC), re)
  X <- as.list(1:7); FUN <- function(x) { Sys.sleep(0.2); x^2 }; ck <- file.path(tempfile("ck67"), "units"); big <- function() list(total_gb = 16, available_gb = 12)
  r1 <- re$lhcmr_schedule(X, FUN, 1L, ck, max_workers = 3L, per_worker_gb = 1, poll = 0.2, mem = big)
  expect("Round 67: the scheduler returns what lapply returns, in order, and saves one checkpoint per unit", identical(r1, lapply(X, function(x) x^2)) && length(list.files(ck)) == 7L)
  t0 <- Sys.time(); r2 <- re$lhcmr_schedule(X, function(x) stop("must not be computed again"), 1L, ck, mem = big)
  expect("Round 67: a second run reads the checkpoints instead of computing (resume after a stop)", identical(r2, r1) && as.numeric(difftime(Sys.time(), t0, units = "secs")) < 5)
  r3 <- re$lhcmr_schedule(list(1, 2), FUN, 1L, ck, mem = big)
  expect("Round 67: a checkpoint is used only for the same input (inputs are keyed by hash)", identical(r3, list(1, 4)))
  sf <- tempfile(fileext = ".json"); calls <- 0L; seen <- character()
  lowfirst <- function() { calls <<- calls + 1L; if (calls <= 2L) list(total_gb = 16, available_gb = 0.4) else list(total_gb = 16, available_gb = 12) }
  ck2 <- file.path(tempfile("ck67"), "units"); r4 <- re$lhcmr_schedule(X[1:3], FUN, 1L, ck2, max_workers = 2L, per_worker_gb = 1, status = sf, poll = 0.2, pause = 0.2,
    mem = function() { m <- lowfirst(); seen <<- c(seen, tryCatch(jsonlite::fromJSON(sf)$state, error = function(e) "")); m })
  expect("Round 67: with too little memory for one worker the run pauses, then continues by itself and completes", identical(r4, list(1, 4, 9)) && "paused" %in% seen)
  if (.Platform$OS.type != "windows") {
    n <- 0L; dip <- function() { n <<- n + 1L; if (n %in% 3:4) list(total_gb = 16, available_gb = 1.5) else list(total_gb = 16, available_gb = 12) }   # memory falls during a batch
    ck3 <- file.path(tempfile("ck67"), "units"); st <- new.env(); slow <- function(x) { Sys.sleep(1); x + 100 }
    sf5 <- tempfile(fileext = ".json"); msgs <- character()
    r5 <- re$lhcmr_schedule(as.list(1:4), slow, 2L, ck3, max_workers = 3L, per_worker_gb = 1, floor_share = 0.15, status = sf5, poll = 0.2, pause = 0.2, state = st,
      mem = function() { msgs <<- c(msgs, tryCatch(jsonlite::fromJSON(sf5)$msg, error = function(e) "")); dip() })
    expect("Round 67: when free memory falls during a batch, the batch is stopped and redone with fewer workers, and the results are unchanged", isTRUE(all.equal(unlist(r5), 101:104)) && any(grepl("Free memory fell", msgs)))
  }
  hj <- list(dir = tempfile("h67")); pj <- list(folder = "p1"); dir.create(file.path(hj$dir, "04_across_pairs", "lhcmr"), recursive = TRUE)
  writeLines(jsonlite::toJSON(list(state = "running", msg = "", pid = 999999L, time = "x"), auto_unbox = TRUE), .lhcmr_status_file(hj$dir, "p1"))
  j1 <- lhcmr_job_read(hj, pj)
  writeLines(jsonlite::toJSON(list(state = "running", msg = "", pid = Sys.getpid(), time = "x"), auto_unbox = TRUE), .lhcmr_status_file(hj$dir, "p1"))
  expect("Round 67: a run whose process is gone is reported as stopped (resumable), a live one as running", identical(j1$state, "stopped") && grepl("continue from the fits already saved", j1$msg) && identical(lhcmr_job_read(hj, pj)$state, "running"))
  sc <- lhcmr_script("p1", c("a", "b"), c("A", "B"), "/ld", "/rho", "/ldsc", "/hm3", 7L, 2L)
  expect("Round 67: the study script embeds the runner and the memory reader, and no longer calls lhc_mr() directly", !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") &&
         any(grepl("^lhcmr_adaptive\\(sp", sc)) && any(grepl("^mem_info <-", sc)) && !any(grepl("^lhc_mr\\(", sc)) && !any(grepl("on.exit(", sc, fixed = TRUE)))
})
# ---- Round 68: MRlap ----
local({
  d <- data.frame(snp = c("rs1", "rs2", "rs3", "rs4"), chr = c("1", "chr2", NA, NA), pos = c(100, 200, NA, NA), ea = "A", oa = "G", beta = c(0.1, 0.2, 0.3, 0.4), se = 0.05, n = 1e5, stringsAsFactors = FALSE)
  pm <- data.frame(snp = "rs3", chr = 3L, pos = 300, stringsAsFactors = FALSE); m <- mrlap_prepare(d, pm); led <- attr(m, "ledger")
  expect("Round 68: MRlap input has rsid chr pos a1 a2 beta se N, positions from the dataset or the 1000 Genomes panel, variants without one left out",
         identical(names(m), c("rsid", "chr", "pos", "a1", "a2", "beta", "se", "N")) && identical(m$rsid, c("rs1", "rs2", "rs3")) && identical(m$chr, c(1L, 2L, 3L)) && led$no_position == 1L)
  sc <- mrlap_script("p1", c("02_datasets/mrlap/a.tsv.gz", "02_datasets/mrlap/b.tsv.gz"), c("A", "B"), "/ld", "/hm3", 5L)
  expect("Round 68: the MRlap script parses, pins MRlap and GenomicSEM, and loads ieugwasr before setting the seed (the corrected SE is simulated)",
         !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_MRLAP, sc)) && any(grepl(PIN_GENOMICSEM, sc)) &&
         which(grepl("loadNamespace(\"ieugwasr\")", sc, fixed = TRUE))[1] < which(grepl("set.seed(seed)", sc, fixed = TRUE))[1])
  res <- data.frame(quantity = c("observed_effect", "corrected_effect", "difference"), estimate = c(0.0856, 0.1146, NA), se = c(0.0398, 0.0528, NA), p = c(0.031, 0.03, 0.037), z = NA,
                    ci_lower = c(0.008, 0.011, NA), ci_upper = c(0.163, 0.218, NA), instruments = 39L, rg = 0.19, int_crosstrait = -0.0011, int_crosstrait_se = 0.0063, h2_exposure = 0.24, h2_outcome = 0.14, stringsAsFactors = FALSE); attr(res, "labels") <- c(x = "BMI", y = "SBP")
  expect("Round 68: the MRlap text gives the observed and corrected effects (Round 100: on the standardised scale, with the size of the change, no preference)", { t <- mrlap_text(res); grepl("39 instruments selected by MRlap", t) && grepl("0.115", t) && grepl("on the standardised scale", t) && grepl("a change of", t) && !grepl("should be preferred", t) })
  sd <- tempfile("st68"); dir.create(sd); writeLines(c("writeLines(getwd(), 'where.txt')"), file.path(sd, "ok.R")); writeLines("stop('planned failure')", file.path(sd, "bad.R"))
  s1 <- file.path(sd, "s1.json"); s2 <- file.path(sd, "s2.json"); script_job_start(sd, "ok.R", s1, file.path(sd, "l1.log")); script_job_start(sd, "bad.R", s2, file.path(sd, "l2.log"))
  rd <- function(f) tryCatch(jsonlite::fromJSON(f)$state, error = function(e) "")
  for (w in 1:60) { if (rd(s1) %in% c("done", "error") && rd(s2) %in% c("done", "error")) break; Sys.sleep(1) }
  expect("Round 68: a study script runs in the background from the study folder and reports done, or error with the reason",
         identical(rd(s1), "done") && file.exists(file.path(sd, "where.txt")) && identical(normalizePath(readLines(file.path(sd, "where.txt"))), normalizePath(sd)) &&
         identical(rd(s2), "error") && grepl("planned failure", tryCatch(jsonlite::fromJSON(s2)$msg, error = function(e) "")))
})
# ---- Round 69: MR-APSS ----
local({
  d <- data.frame(snp = c("rs1", "rs2", "rs3"), ea = "A", oa = "G", beta = c(0.1, NA, 0.3), se = 0.05, eaf = c(0.3, 0.4, 0.2), n = 1e5, stringsAsFactors = FALSE)
  m <- mrapss_prepare(d)
  expect("Round 69: MR-APSS input uses the column names format_data() reads (SNP A1 A2 b se n freq) and drops unusable rows",
         identical(names(m), c("SNP", "A1", "A2", "b", "se", "n", "freq")) && identical(m$SNP, c("rs1", "rs3")) && attr(m, "ledger")$kept == 2L)
  sc <- mrapss_script("p1", c("02_datasets/mrapss/a.tsv.gz", "02_datasets/mrapss/b.tsv.gz"), c("A", "B"), "/ld", "/panel/EUR", "/bin/plink", 9L)
  expect("Round 69: the MR-APSS script parses, pins MRAPSS, clumps locally with PLINK and loads ieugwasr before the seed",
         !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_MRAPSS, sc)) && any(grepl("bfile = BFILE, plink_bin = PLINK", sc, fixed = TRUE)) &&
         which(grepl("loadNamespace(\"ieugwasr\")", sc, fixed = TRUE))[1] <= which(grepl("set.seed(seed)", sc, fixed = TRUE))[1])
  res <- data.frame(estimate = 0.216, se = 0.026, ci_lower = 0.165, ci_upper = 0.268, p = 1.9e-16, instruments = 1227L, pi0 = 0.23, sigma_sq = 2e-6, tau_sq = 5e-7, iv_threshold_corrected = 2.2e-5, c12 = 0.05); attr(res, "labels") <- c(x = "BMI", y = "T2D")
  expect("Round 69: the MR-APSS text gives the estimate, the instruments and what was corrected", { t <- mrapss_text(res); grepl("1227 instruments at P < 5e-5", t) && grepl("BMI on T2D 0.216", t) })
})
# ---- Round 70b: a method whose package is missing is locked before the run ----
expect("Round 70b: a missing package is detected, an installed one is not", identical(.pkg_missing(c("stats", "NoSuchPackageMRS")), "NoSuchPackageMRS"))
expect("Round 70e: BWMR needs no package, so it is never locked for a missing one", { b <- feasibility(list())$methods$bwmr; is.null(b) || !identical(b$state, "locked") })
local({ uni <- list(data = data.frame(beta_x = 1, se_x = 1, beta_y = 1, se_y = 1)); mv <- list(data = data.frame(beta_x_1 = 1, se_x_1 = 1, beta_x_2 = 1, se_x_2 = 1, beta_y = 1, se_y = 1))
  a <- runnable_methods(c("mvivw", "ivw_mre", "egger"), uni); b <- runnable_methods(c("mvivw", "ivw_mre"), mv)
  expect("Round 70b: multivariable methods are left out of a one-exposure pair, and one-exposure methods out of a multivariable pair, with the reason",
         identical(a$methods, c("ivw_mre", "egger")) && grepl("two or more exposures", a$left_out[["mvivw"]]) && identical(b$methods, "mvivw") && grepl("one-exposure method", b$left_out[["ivw_mre"]])) })
expect("Round 72c: an ancestry with no reference panel keeps OpenGWAS's own word instead of a question mark", identical(.ogw_ancestry("Mixed"), "Mixed") && identical(.ogw_ancestry("European"), "EUR") && identical(.ogw_ancestry(NULL), "not reported"))
# ---- Round 72c: the reverse direction of a bidirectional study from OpenGWAS downloads ----
local({ root <- file.path(tempfile("bidirrev"), "MR Studio"); dir.create(library_dir(root), recursive = TRUE)
  mk <- function(name, src, usable_x, usable_y) { dd <- file.path(library_dir(root), name); dir.create(dd)
    lib_table_write(data.frame(snp = "rs1", ea = "A", oa = "G", beta = 0.1, se = 0.01), dd, "standardised")
    atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = name, usable_as_exposure = usable_x, usable_as_outcome = usable_y, source = src), auto_unbox = TRUE)) }
  px <- list(pval = 5e-8, clump = 1, r2 = 0.001, kb = 10000, pop = "EUR", where = "server_live"); py <- list(proxies = TRUE, rsq = 0.8, maf_threshold = 0.3)
  mk("LDL", list(type = "opengwas", role = "exposure", id = "ieu-a-300", trait = "LDL cholesterol", params = px), TRUE, FALSE)
  mk("CHD for LDL", list(type = "opengwas", role = "outcome", id = "ieu-a-7", trait = "Coronary heart disease", params = py, for_exposure = "LDL"), FALSE, TRUE)
  calls <- character()
  gx <- function(id, jwt, p, r2, kb, pop, name, offline, root, where) { calls <<- c(calls, sprintf("x %s %s %s %s", id, format(p), pop, where)); mk(name, list(type = "opengwas", role = "exposure", id = id, trait = "CHD", params = px), TRUE, FALSE); list(trait = name) }
  gy <- function(id, exposure_name, jwt, proxies, rsq, maf, name, offline, root) { calls <<- c(calls, sprintf("y %s for %s %s", id, exposure_name, rsq)); mk(name, list(type = "opengwas", role = "outcome", id = id, trait = "LDL", params = py, for_exposure = exposure_name), FALSE, TRUE); list(trait = name) }
  said <- list(); options(mrstudio.say = function(msg, why) said[[length(said) + 1L]] <<- list(msg = msg, why = why))
  r1 <- bidir_reverse("lib:LDL", "lib:CHD for LDL", "tok", root = root, get_exposure = gx, get_outcome = gy); n1 <- length(said)
  r2 <- bidir_reverse("lib:LDL", "lib:CHD for LDL", "tok", root = root, get_exposure = gx, get_outcome = gy); options(mrstudio.say = NULL)
  expect("Round 75d: the background work is narrated step by step, each step with a one-sentence explanation, and a reused download says so",
         n1 >= 3L && any(grepl("^Reverse direction", vapply(said, `[[`, "", "msg"))) && any(grepl("^Downloading the instruments", vapply(said, `[[`, "", "msg"))) &&
         all(nzchar(vapply(said[1:3], `[[`, "", "why"))) && any(grepl("reused", vapply(said[-seq_len(n1)], `[[`, "", "msg"))))
  expect("Round 72c: the reverse direction downloads the outcome trait as an exposure with the forward settings, and the exposure trait as an outcome for it",
         identical(calls, c("x ieu-a-7 5e-08 EUR server_live", sprintf("y ieu-a-300 for %s 0.8", r1$exName))) && identical(r1$ex, paste0("lib:", r1$exName)) && grepl("downloaded", r1$note))
  expect("Round 72c: a second preparation reuses the reverse downloads", length(calls) == 2L && identical(r2[c("ex", "out")], r1[c("ex", "out")]) && identical(r2$note, ""))
  mk("My file", list(type = "import"), TRUE, FALSE); mk("Outcome for My file", list(type = "import", for_exposure = "My file"), FALSE, TRUE)
  expect("Round 72c: data that is neither full summary statistics nor from OpenGWAS gets a clear message", grepl("needs .* as an exposure", err_of(bidir_reverse("lib:My file", "lib:Outcome for My file", root = root)))) })
# ---- Round 72c: multivariable data from OpenGWAS downloads ----
local({ root <- file.path(tempfile("ogwmv"), "MR Studio"); dir.create(library_dir(root), recursive = TRUE)
  set.seed(72); snps <- sprintf("rs%d", 101:140); b1 <- rnorm(40, 0, .05); b2 <- rnorm(40, 0, .05); by <- 0.3 * b1 + 0 * b2 + rnorm(40, 0, .002)
  truth <- list(`ieu-x1` = b1, `ieu-x2` = b2, `ieu-y` = by)
  mkx <- function(name, id, rows) { dd <- file.path(library_dir(root), name); dir.create(dd)
    lib_table_write(data.frame(snp = snps[rows], ea = "A", oa = "C", eaf = .3, beta = truth[[id]][rows], se = .002, p = 1e-12 / seq_along(rows), stringsAsFactors = FALSE), dd, "standardised")
    atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = name, usable_as_exposure = TRUE, usable_as_outcome = FALSE, facts = list(ancestry = "EUR", n = 1e5, kind = "continuous"),
      original = list(sha256 = name), source = list(type = "opengwas", role = "exposure", id = id, params = list(pval = 5e-8))), auto_unbox = TRUE)) }
  mkx("X1", "ieu-x1", 1:22); mkx("X2", "ieu-x2", 15:40)
  dd <- file.path(library_dir(root), "Y for X1"); dir.create(dd); lib_table_write(data.frame(snp = "rs101", ea = "A", oa = "C", beta = 0, se = 1), dd, "standardised")
  atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = "Y for X1", usable_as_exposure = FALSE, usable_as_outcome = TRUE, facts = list(ancestry = "EUR", kind = "continuous"),
    source = list(type = "opengwas", role = "outcome", id = "ieu-y", for_exposure = "X1", params = list(proxies = TRUE, rsq = 0.9, maf_threshold = 0.3))), auto_unbox = TRUE))
  asked <- list()
  bk <- list(clump = function(snps, p, r2, kb, pop, jwt) snps,
             outcome = function(q, id, proxies, rsq, maf, jwt) { asked[[length(asked) + 1L]] <<- list(id = id, proxies = proxies, rsq = rsq); i <- match(q, snps)
               data.frame(SNP = q, effect_allele.outcome = "A", other_allele.outcome = "C", eaf.outcome = .3, beta.outcome = truth[[id]][i], se.outcome = .002, pval.outcome = .01, samplesize.outcome = 1e5, stringsAsFactors = FALSE) })
  r1 <- ogw_mv_build(c("lib:X1", "lib:X2"), "lib:Y for X1", "tok", root = root, backend = bk, panel = FALSE)
  r2 <- ogw_mv_build(c("lib:X1", "lib:X2"), "lib:Y for X1", "tok", root = root, backend = bk, panel = FALSE)
  inp <- ogw_mv_input(c("lib:X1", "lib:X2"), "lib:Y for X1", root = root)$input$data; inp <- inp[inp$mr_keep, ]
  f <- stats::lm(beta_y ~ 0 + beta_x_1 + beta_x_2, data = inp, weights = 1 / inp$se_y^2)
  expect("Round 72c: OpenGWAS multivariable set: the union of both exposures' instruments, every exposure and the outcome looked up there, the outcome with its own proxy settings",
         nrow(inp) == 40L && r1$ledger$union == 40L && length(asked) == 3L && identical(asked[[3]]$id, "ieu-y") && isTRUE(asked[[3]]$proxies) && asked[[3]]$rsq == 0.9 && !isTRUE(asked[[1]]$proxies))
  expect("Round 72c: OpenGWAS multivariable set recovers the direct effects (0.3 and 0) and is reused while its downloads are unchanged",
         abs(coef(f)[[1]] - 0.3) < 0.02 && abs(coef(f)[[2]]) < 0.02 && isTRUE(r2$reused) && identical(r2$name, r1$name) && length(asked) == 3L)
  old <- Sys.getenv("MRGUI_STUDIES_DIR"); Sys.setenv(MRGUI_STUDIES_DIR = root); on.exit(Sys.setenv(MRGUI_STUDIES_DIR = old))
  expect("Round 72c: a multivariable pair of OpenGWAS downloads uses the prepared set", identical(mv_pair_input(c("lib:X1", "lib:X2"), "lib:Y for X1")$input$kind, "multi_iv")) })
# ---- Round 72c: mediation from OpenGWAS downloads ----
local({ root <- file.path(tempfile("ogwmed"), "MR Studio"); dir.create(library_dir(root), recursive = TRUE)
  set.seed(73); snps <- sprintf("rs%d", 201:240); bx <- rnorm(40, 0, .05); bm <- 0.5 * bx + rnorm(40, 0, .05); by <- 0.2 * bx + 0.4 * bm + rnorm(40, 0, .002)
  truth <- list(`ieu-x` = bx, `ieu-m` = bm, `ieu-y` = by)
  rec <- function(name, src, ux, uy, rows) { dd <- file.path(library_dir(root), name); dir.create(dd)
    lib_table_write(data.frame(snp = snps[rows], ea = "A", oa = "C", eaf = .3, beta = truth[[src$id]][rows], se = .002, p = 1e-12 / seq_along(rows), stringsAsFactors = FALSE), dd, "standardised")
    atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(list(trait = name, usable_as_exposure = ux, usable_as_outcome = uy, facts = list(ancestry = "EUR", n = 1e5, kind = "continuous"),
      original = list(sha256 = name), source = src), auto_unbox = TRUE)); list(trait = name) }
  px <- list(pval = 5e-8, clump = 1, r2 = 0.001, kb = 10000, pop = "EUR", where = "server_live"); py <- list(proxies = TRUE, rsq = 0.8, maf_threshold = 0.3)
  rec("X", list(type = "opengwas", role = "exposure", id = "ieu-x", params = px), TRUE, FALSE, 1:20)
  rec("M for X", list(type = "opengwas", role = "outcome", id = "ieu-m", for_exposure = "X", params = py), FALSE, TRUE, 1:20)
  rec("Y for X", list(type = "opengwas", role = "outcome", id = "ieu-y", for_exposure = "X", params = py), FALSE, TRUE, 1:20)
  gx <- function(id, jwt, p, r2, kb, pop, name, offline, root, where) rec(name, list(type = "opengwas", role = "exposure", id = id, params = px), TRUE, FALSE, 15:40)
  gy <- function(id, exposure_name, jwt, proxies, rsq, maf, name, offline, root) stop("not expected")
  bk <- list(clump = function(snps, p, r2, kb, pop, jwt) snps, outcome = function(q, id, proxies, rsq, maf, jwt) { i <- match(q, snps)
    data.frame(SNP = q, effect_allele.outcome = "A", other_allele.outcome = "C", eaf.outcome = .3, beta.outcome = truth[[id]][i], se.outcome = .002, pval.outcome = .01, samplesize.outcome = 1e5, stringsAsFactors = FALSE) })
  r <- mediation_prepare("lib:X", "lib:M for X", "lib:Y for X", "tok", root = root, get_exposure = gx, get_outcome = gy, backend = bk, panel = FALSE)
  old <- Sys.getenv("MRGUI_STUDIES_DIR"); Sys.setenv(MRGUI_STUDIES_DIR = root); on.exit(Sys.setenv(MRGUI_STUDIES_DIR = old))
  d <- mediation_pair_input(list(role = "direct", x = "lib:X", ms = list("lib:M for X"), y = "lib:Y for X", m = "lib:M for X", m_out = r$m_out[[1]], m_inst = r$m_inst[[1]]))$input$data; d <- d[d$mr_keep, ]
  f <- stats::lm(beta_y ~ 0 + beta_x_1 + beta_x_2, data = d, weights = 1 / d$se_y^2)
  expect("Round 72c: a mediator downloaded as an outcome is also downloaded as its own instruments, and the direct-effect step recovers 0.2 (exposure) and 0.4 (mediator)",
         identical(r$m_out[[1]], "lib:M for X") && grepl("^lib:", r$m_inst[[1]]) && r$m_inst[[1]] != "lib:M for X" && abs(coef(f)[[1]] - 0.2) < 0.03 && abs(coef(f)[[2]] - 0.4) < 0.03) })
# ---- Round 72e: MVMR heterogeneity review and reruns without chosen variants ----
local({ inp <- mv_pair_input(c("expA", "expB"), "outA"); h <- studio_new("MV review")
  r <- studio_run(h, inp$exposure_label, inp$outcome_label, inp$input, c("mvivw", "mvlasso"), idx = 1L)
  rv <- mv_review(h, r$pair_folder); q <- Filter(function(d) identical(d$method, "mvivw") && identical(d$name, "mv_q"), pair_payload(h$dir, r$pair_folder)$diagnostics)[[1]]
  expect("Round 72e: the per-variant contributions add up to Cochran's Q of MVMR-IVW, largest first", abs(rv$q_total - as.numeric(q$value)) < 1e-6 && !is.unsorted(rev(rv$rows$q)))
  s2 <- studio_run_reduced(h, inp, rv$rows$variant_id[1:2], "mvivw", 0.95, idx = 2L)
  expect("Round 72e: an MVMR pair reruns without chosen variants as a sensitivity analysis", identical(s2$results$mvivw$status, "ok") && s2$kept == sum(inp$input$data$mr_keep) - 2L) })
local({ root <- file.path(tempfile("fullgw"), "MR Studio"); dir.create(library_dir(root), recursive = TRUE)
  mk <- function(name, rec) { dd <- file.path(library_dir(root), name); dir.create(dd); lib_table_write(data.frame(snp = "rs1", ea = "A", oa = "G", beta = 0.1, se = 0.01), dd, "standardised")
    atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(c(list(trait = name), rec), auto_unbox = TRUE)) }
  mk("OG", list(usable_as_exposure = TRUE, usable_as_outcome = FALSE, rows_kept = 79, source = list(type = "opengwas", id = "ieu-a-300")))
  mk("OGY", list(usable_as_exposure = FALSE, usable_as_outcome = TRUE, rows_kept = 77, source = list(type = "opengwas", id = "ieu-a-7", for_exposure = "OG")))
  mk("TAB", list(usable_as_exposure = TRUE, usable_as_outcome = FALSE, rows_kept = 50, source = list(type = "import")))
  mk("FX", list(usable_as_exposure = TRUE, usable_as_outcome = TRUE, rows_kept = 2.4e6, source = list(type = "import")))
  mk("FY", list(usable_as_exposure = TRUE, usable_as_outcome = TRUE, rows_kept = 1.1e6, source = list(type = "import")))
  expect("Round 72e: complete-GWAS analyses are offered only when both traits are complete summary statistics, with the reason otherwise",
         identical(full_gwas_status("lib:OG", "lib:OGY", root)$why, "opengwas") && identical(full_gwas_status("lib:TAB", "lib:FY", root)$why, "instruments") &&
         isTRUE(full_gwas_status("lib:FX", "lib:FY", root)$ok) && identical(full_gwas_status("expA", "outA", root)$why, "simulated")) })
# ---- Round 72c: IVW after Steiger filtering without the prevalence of a binary outcome is left out, not failed ----
local({ d <- data.frame(variant_id = sprintf("rs%d", 1:12), ea = "A", oa = "G", beta_x = seq(.05, .16, length.out = 12), se_x = .01, eaf_x = .3, beta_y = seq(.05, .16, length.out = 12) * .3, se_y = .02, eaf_y = .3,
                        p_x = 1e-10, p_y = .01, n_x = 1e5, n_y = 2e5, ncase_y = 6e4, ncontrol_y = 1.4e5, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  inp <- list(kind = "pair_iv", data = d, exposure_id = "x", outcome_id = "y", outcome_scale = "log_or", exposure_scale = "beta")
  h <- studio_new("Steiger prevalence"); r <- studio_run(h, "X", "Y binary", inp, c("ivw_mre", "ivw_steiger"), idx = 1L, steiger_in = inp)
  expect("Round 72c: without the prevalence of a binary outcome, IVW after Steiger filtering is left out with the reason and IVW still runs",
         identical(names(r$results), "ivw_mre") && identical(r$results$ivw_mre$status, "ok") && grepl("prevalence", r$left_out[["ivw_steiger"]])) })
# ---- Round 70: CAUSE ----
local({
  d <- data.frame(snp = c("rs1", "rs2", "rs3", "rs1", NA), ea = c("a", "C", "G", "a", "T"), oa = c("g", "T", "A", "g", "C"), beta = c(0.1, NA, 0.3, 0.1, 0.2), se = 0.05, stringsAsFactors = FALSE)
  m <- cause_prepare(d)
  expect("Round 70: CAUSE input has the columns gwas_merge() is told to read, upper-case alleles, no duplicate, missing or unusable rows",
         identical(names(m), c("snp", "beta", "se", "ea", "oa")) && identical(m$snp, c("rs1", "rs3")) && identical(m$ea, c("A", "G")) && attr(m, "ledger")$kept == 2L)
  sc <- cause_script("p1", c("02_datasets/cause/a.tsv.gz", "02_datasets/cause/b.tsv.gz"), c("A", "B"), "/panel/EUR", "/bin/plink", 9L)
  expect("Round 70: the CAUSE script parses, pins cause, clumps locally with PLINK at P < 1e-3 and r2 < 0.01, and loads ieugwasr before the seed",
         !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_CAUSE, sc)) && any(grepl("plink_bin = PLINK, bfile = BFILE", sc, fixed = TRUE)) &&
         any(grepl("P_THRESHOLD <- 0.001; R2_THRESHOLD <- 0.01; N_SAMPLE <- 1000000", sc, fixed = TRUE)) &&
         which(grepl("loadNamespace(\"ieugwasr\")", sc, fixed = TRUE))[1] < which(grepl("set.seed(seed)", sc, fixed = TRUE))[1])
  res <- data.frame(gamma = 0.32, gamma_lower = 0.24, gamma_upper = 0.40, eta = -0.03, q = 0.03, sharing_eta = 0.39, sharing_q = 0.62, z = -4.337, p = 7.2e-6, instruments = 719L, rho = 0.065, delta_elpd = -7.075, se_delta_elpd = 1.631)
  attr(res, "labels") <- c(x = "LDL cholesterol", y = "CAD")
  expect("Round 70: the CAUSE text says whether the causal model fits better, with P and the posterior of gamma",
         { t <- cause_text(res); grepl("719 variants", t) && grepl("fitted better than the sharing model (P = 7.2e-06)", t, fixed = TRUE) && grepl("LDL cholesterol on CAD gamma 0.32 (95% credible interval 0.24 to 0.40)", t, fixed = TRUE) })
  res$p <- 0.4; expect("Round 70: a P value above 0.05 is reported as not fitting better", grepl("did not fit better than", cause_text(res)))
  expect("Round 70: the CAUSE reference is in the reference list", identical(REFS$morrison2020$pmid, "32451458"))
})
expect("Round 50: the demo recovers the simulated direct effects (A on outcome A 0.30, B none) with MVMR-IVW", { iv <- d50[d50$method == "mvivw", ]; a <- iv[iv$exposure_id == "demo:expA", ]; b <- iv[iv$exposure_id == "demo:expB", ]; a$ci_lower <= 0.3 && 0.3 <= a$ci_upper && b$ci_lower <= 0 && 0 <= b$ci_upper })
expect("Round 50: too few variants for the number of exposures are blocked with the reason", { q <- p50$input; q$data <- q$data[1:3, ]; v <- ADAPTER_MVEGGER_V11$check(q, list()); identical(v$status, "block") && grepl("needs at least 4 variants", paste(v$messages, collapse = " ")) })
expect("Round 50: one exposure, more than six, or mixing simulated and imported traits is refused", grepl("at least two exposures", err_of(mv_pair_input("expA", "outA"))) && grepl("at most six", err_of(mv_pair_input(rep("expA", 7), "outA"))) && grepl("cannot be mixed", err_of(mv_pair_input(c("expA", "lib:X"), "outA"))))
ui_state_set(h50$con, "pairs", list(list(ex = list("expA", "expB"), out = "outA", mv = TRUE, exName = p50$exposure_label, exNames = as.list(p50$exposure_labels), outName = p50$outcome_label, folder = r50$pair_folder, methods = as.list(m50), ci = 0.95, pthr = 5e-8, simulated = TRUE, ledger = p50$ledger)))
t50 <- report_text(list(mv_context(h50, ui_state_get(h50$con, "pairs")[[1]])), "MV")
expect("Round 50: report text gives the direct effects, conditional F, Q and intercept, with the MVMR references", any(grepl("^MVMR of Simulated exposure A, Simulated exposure B on Simulated outcome A", t50$results)) && any(grepl("Conditional F statistics: Simulated exposure A", t50$results)) && any(grepl("\\[Sanderson 2021\\]", t50$methods)) && any(grepl("SIMULATED DATA", t50$limitations)))
b50 <- build_package(h50, c("report", "figs", "tables"))
expect("Round 50: the package has the MVMR forest of direct effects and the PDF report", any(grepl("P1_mvmr_forest\\.png$", b50$files)) && any(grepl("report\\.pdf$", b50$files)) && !length(b50$notes))
f50 <- mv_pair_facts(c("expA", "expB"), "outA")
expect("Round 50: Check-screen facts for a multivariable pair", isTRUE(f50$mv) && f50$instruments == 60L && length(f50$mean_f_each) == 2L)
# library datasets: two exposures with full summary statistics, joint instruments, every exposure looked up at every instrument
set.seed(50); n5 <- 120; s5 <- sprintf("rs%d", 50000 + seq_len(n5)); g5 <- rep(1:3, each = 40)
b1 <- ifelse(g5 == 1, runif(n5, .05, .1), rnorm(n5, 0, .003)); b2 <- ifelse(g5 == 2, runif(n5, .05, .1), ifelse(g5 == 1, .3 * b1, rnorm(n5, 0, .003))); se5 <- rep(.008, n5)
y5 <- 0.3 * b1 + 0 * b2 + rnorm(n5, 0, .01)
w5 <- function(b, nm, role, ind) { f <- tempfile(fileext = ".csv"); utils::write.csv(data.frame(snp = s5, ea = "A", oa = "G", eaf = .3, beta = b, se = se5, p = 2 * pnorm(-abs(b / se5)), n = 100000), f, row.names = FALSE)
  invisible(import_save(f, paste0(nm, ".csv"), ",", c("snp", "ea", "oa", "eaf", "beta", "se", "p", "n"), mpl, fct(nm, role, ind))) }
w5(b1, "MV X1", "either", "yes"); w5(b2, "MV X2", "either", "yes"); w5(y5, "MV Y", "outcome", "no")
lm50 <- library_mv_input(c("lib:MV X1", "lib:MV X2"), "lib:MV Y", pthr = 1e-6, clump = function(snp, p, pop) data.frame(variant_id = snp, status = "index", stringsAsFactors = FALSE))
expect("Round 50: library MVMR data hold the union of both exposures' instruments with each exposure's effect at every one", lm50$ledger$mv$union == 80L && lm50$ledger$mv$final == 80L && all(c("beta_x_1", "beta_x_2", "n_x_1") %in% names(lm50$input$data)) && identical(lm50$input$kind, "multi_iv"))
expect("Round 50: library MVMR data recover the direct effects (0.3 and 0)", { d <- lm50$input$data; f <- MendelianRandomization::mr_mvivw(MendelianRandomization::mr_mvinput(bx = as.matrix(d[, c("beta_x_1", "beta_x_2")]), bxse = as.matrix(d[, c("se_x_1", "se_x_2")]), by = d$beta_y, byse = d$se_y)); abs(f@Estimate[1] - 0.3) < 0.05 && abs(f@Estimate[2]) < 0.05 })
m50x <- err_of(library_mv_input(c("lib:MV X1", "lib:Trait X"), "lib:MV Y", 1e-6, clump = function(snp, p, pop) data.frame(variant_id = snp, status = "index")))
expect("Round 50: an exposure that holds only its own instruments is refused for MVMR", grepl("holds only its own instruments|full summary statistics", m50x))

# ---- Round 51: MVMR keeps palindromic variants when every exposure is declared forward strand (found in the Rees 2017 validation) ----
w51 <- function(b, nm, role, ind, strand) { f <- tempfile(fileext = ".csv"); al <- ifelse(seq_len(n5) %% 4 == 0, "T", "G")
  utils::write.csv(data.frame(snp = s5, ea = "A", oa = al, beta = b, se = se5, p = 2 * pnorm(-abs(b / se5)), n = 100000), f, row.names = FALSE)
  invisible(import_save(f, paste0(nm, ".csv"), ",", c("snp", "ea", "oa", "beta", "se", "p", "n"), mpl[setdiff(names(mpl), "eaf")], utils::modifyList(fct(nm, role, ind), list(frq = "none", strand = strand)))) }
for (st in c("forward", "unknown")) { w51(b1, paste("P1", st), "either", "yes", st); w51(b2, paste("P2", st), "either", "yes", st); w51(y5, paste("PY", st), "outcome", "no", st) }
cl51 <- function(snp, p, pop) data.frame(variant_id = snp, status = "index", stringsAsFactors = FALSE)
f51 <- library_mv_input(c("lib:P1 forward", "lib:P2 forward"), "lib:PY forward", 1e-6, clump = cl51); u51 <- library_mv_input(c("lib:P1 unknown", "lib:P2 unknown"), "lib:PY unknown", 1e-6, clump = cl51)
expect("Round 51: declared forward strand keeps the palindromic (A/T) instruments in MVMR; undeclared, they are dropped without allele frequencies",
       f51$ledger$mv$final == 80L && isTRUE(f51$ledger$mv$palindromic_kept) && u51$ledger$mv$final == 60L && !isTRUE(u51$ledger$mv$palindromic_kept))

# ---- Round 75: complete GWAS inside the app (OpenGWAS GWAS-VCF downloads, VCF imports, free-space check) ----
local({ root <- file.path(tempfile("fg75"), "MR Studio"); dir.create(library_dir(root), recursive = TRUE)
  vf <- file.path(dirname(root), "t.vcf.gz"); con <- gzfile(vf, "w")
  set.seed(75); nn <- 300; bb <- round(rnorm(nn, 0, .02), 5); ss <- round(runif(nn, .005, .02), 5); lp <- round(-log10(2 * pnorm(-abs(bb / ss))), 4)
  body <- sprintf("%d\t%d\trs%d\t%s\t%s\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t%s:%s:%s:0.3:250000:rs%d", (seq_len(nn) %% 22) + 1L, seq_len(nn) * 1000L, seq_len(nn), "A", "G", bb, ss, lp, seq_len(nn))
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>", sprintf("##SAMPLE=<ID=ieu-a-300,TotalVariants=%d,StudyType=Continuous>", nn + 4),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tieu-a-300", body,
    "1\t9001\trs9001\tC\tT,G\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.02:5:0.2:250000:rs9001", "1\t9002\trs9002\tCT\tC\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.1:0.02:5:0.2:250000:rs9002",
    "1\t9003\t.\tG\tC\t.\tPASS\t.\tSE:ES:SS\t0.01:-0.05:250000", "1\t9004\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP:AF:SS:ID\t0.2:0.03:2:0.3:250000:rs1"), con); close(con)
  r <- vcf_read_gwas(vf, chunk = 100L)
  expect("Round 75: GWAS-VCF is read in chunks: effect of ALT (ea = ALT, oa = REF), P from LP, fields by FORMAT name in any order, with every row left out counted",
         r$stats$rows_read == nn + 4 && r$stats$rows_kept == nn + 1 && identical(unlist(r$stats$dropped), c(not_snv = 1, multiallelic = 1, missing = 0, duplicate = 1)) &&
         all(r$data$ea[1:nn] == "G") && all(r$data$oa[1:nn] == "A") && max(abs(r$data$p[1:nn] - 10^-lp)) < 1e-12 && r$data$beta[nn + 1] == -0.05 && r$data$snp[nn + 1] == "1:9003" && all(r$data$n == 250000))
  g <- file.path(dirname(root), "geno.vcf"); writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878", "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1"), g)
  nf <- file.path(dirname(root), "x.txt"); writeLines("a\tb", nf)
  expect("Round 75: a genotype VCF and a file that is not a VCF are refused with the reason", grepl("not a GWAS summary statistics", err_of(vcf_peek(g))) && grepl("not a VCF", err_of(vcf_peek(nf))))
  fb <- disk_free_bytes(root); expect("Round 75: the free space of the disk that holds the Data library is measured", is.finite(fb) && fb > 0)
  # the plan: size from the server, compared with the free space
  B <- getOption("mrgui.ogw_backend"); B$files <- function(id, jwt) { chk(jwt); list(list(paste0("https://example.org/files/", id, ".vcf.gz?sig=1"), paste0("https://example.org/files/", id, ".vcf.gz.tbi?sig=1"))) }
  pl <- function(free) fullgwas_plan("sim-exp-1", tok_ok, root = root, backend = B, size_fun = function(u) { if (!grepl("\\.vcf\\.gz\\?", u)) stop("wrong link"); 5e8 }, free_fun = function(p) free)
  p1 <- pl(1e9); p2 <- pl(2.5e9); p3 <- pl(1e10)
  expect("Round 75: before a download the space needed (download + standard copy) is compared with the free space: block, warn or ok, with the numbers",
         identical(c(p1$level, p2$level, p3$level), c("block", "warn", "ok")) && abs(p1$need_bytes - (5e8 + max(2.2 * 5e8, 60 * 1e7))) < 1 && grepl("Not enough free space", p1$msg) && grepl("only about", p2$msg))
  expect("Round 75: a download that does not fit cannot start", grepl("Not enough free space", err_of(fullgwas_download_start("sim-exp-1", "X", tok_ok, root = root, backend = B, plan = p1))))
  expect("Round 75: a dataset without a complete GWAS file says so", grepl("no complete summary statistics", err_of(ogw_vcf_link("sim-exp-1", tok_ok, list(files = function(id, jwt) list())))))
  # the job: download with resume (file:// link, half already downloaded), read, saved in the library; the original kept
  spec <- list(key = "ieu-a-300", mode = "opengwas", url = paste0("file://", normalizePath(vf)), bytes = file.size(vf), name = "LDL complete", keep_vcf = TRUE,
               facts = list(role = "either", kind = "continuous", ancestry = "EUR", independent = "clump", clump_pop = "EUR", scope = "genomewide"), source = list(type = "opengwas_vcf", id = "ieu-a-300", trait = "LDL cholesterol"),
               library = normalizePath(library_dir(root)))
  dir <- .fg_dl_dir("ieu-a-300", root); dir.create(dir, recursive = TRUE); writeLines(jsonlite::toJSON(spec, auto_unbox = TRUE, digits = NA), file.path(dir, "spec.json"))
  raw <- readBin(vf, "raw", file.size(vf)); writeBin(raw[seq_len(length(raw) %/% 2)], file.path(dir, "gwas.vcf.gz.part"))
  fullgwas_job_run(dir); js <- fullgwas_job_read("ieu-a-300", root); ds <- library_dataset("LDL complete", root)
  expect("Round 75: a download resumes from the part already downloaded, is checked by size, read and saved as a complete GWAS; the original VCF is kept with its SHA-256",
         identical(js$state, "done") && nrow(ds$data) == nn + 1 && isTRUE(ds$recipe$complete_gwas) && isTRUE(ds$recipe$usable_as_outcome) && isTRUE(ds$recipe$usable_as_exposure) &&
         identical(ds$recipe$source$sha256, digest::digest(file = vf, algo = "sha256")) && file.exists(file.path(ds$dir, "ieu-a-300.vcf.gz")) && !file.exists(file.path(dir, "gwas.vcf.gz.part")))
  expect("Round 75: the GWAS-VCF dataset is offered as an exposure (clumped on this computer) and as an outcome", identical(ds$recipe$facts$independent, "clump") && any(vapply(library_catalogue(root)$outcomes, function(x) x$name == "LDL complete", TRUE)))
  expect("Round 75: a finished dataset name cannot be downloaded again over it", grepl("already exists", err_of(fullgwas_job_start(spec, root))))
  # linking: the OpenGWAS instrument download of the same dataset uses the complete GWAS for the complete-GWAS analyses
  mk <- function(name, rec) { dd <- file.path(library_dir(root), name); dir.create(dd); lib_table_write(data.frame(snp = "rs1", ea = "A", oa = "G", beta = 0.1, se = 0.01), dd, "standardised")
    atomic_write_text(file.path(dd, "recipe.json"), jsonlite::toJSON(c(list(trait = name), rec), auto_unbox = TRUE)) }
  mk("LDL", list(usable_as_exposure = TRUE, usable_as_outcome = FALSE, rows_kept = 79, source = list(type = "opengwas", id = "ieu-a-300", trait = "LDL cholesterol")))
  mk("CHD", list(usable_as_exposure = FALSE, usable_as_outcome = TRUE, rows_kept = 77, source = list(type = "opengwas", id = "ieu-a-7", trait = "Coronary heart disease", for_exposure = "LDL")))
  mk("CHD complete", list(usable_as_exposure = TRUE, usable_as_outcome = TRUE, rows_kept = 9e6, complete_gwas = TRUE, source = list(type = "opengwas_vcf", id = "ieu-a-7")))
  rj <- jsonlite::fromJSON(file.path(ds$dir, "recipe.json"), simplifyVector = FALSE); rj$rows_kept <- 2.4e6; writeLines(jsonlite::toJSON(rj, auto_unbox = TRUE, null = "null", digits = NA), file.path(ds$dir, "recipe.json"))
  s1 <- full_gwas_status("lib:LDL", "lib:CHD", root)
  expect("Round 75: a pair of OpenGWAS downloads uses the complete GWAS of the same OpenGWAS datasets when they are downloaded, and says so",
         isTRUE(s1$ok) && identical(full_gwas_ids(c("lib:LDL", "lib:CHD"), root), c("lib:LDL complete", "lib:CHD complete")) && length(s1$linked) == 2L)
  unlink(file.path(library_dir(root), "CHD complete"), recursive = TRUE); s2 <- full_gwas_status("lib:LDL", "lib:CHD", root)
  expect("Round 75: without it, the pair names the OpenGWAS dataset to download", !isTRUE(s2$ok) && identical(s2$why, "opengwas") && length(s2$need) == 1L && identical(s2$need[[1]]$id, "ieu-a-7"))
  # a job whose process has gone is reported as stopped, ready to resume; discarding removes only the download folder
  d2 <- .fg_dl_dir("ieu-a-7", root); dir.create(d2, recursive = TRUE); .fg_write_status(d2, state = "running", step = "downloading", pid = 999999L, mode = "opengwas")
  expect("Round 75: a download whose process has gone (app or computer closed) is shown as stopped, to resume", identical(fullgwas_job_read("ieu-a-7", root)$state, "cancelled"))
  fullgwas_discard("ieu-a-300", root); expect("Round 75: removing a finished download clears it from the list and keeps the dataset", !dir.exists(dir) && dir.exists(ds$dir))
  # a VCF the user has, imported in the background (the separate R process)
  j <- vcf_import_start(vf, "My VCF", list(ancestry = "EUR"), root); t0 <- Sys.time()
  repeat { Sys.sleep(1); s <- fullgwas_job_read(j$key, root); if (!identical(s$state, "running") || difftime(Sys.time(), t0, units = "secs") > 90) break }
  expect("Round 75: a GWAS-VCF file is imported in a background R process, and the user's file is not changed or moved",
         identical(s$state, "done") && nrow(library_dataset("My VCF", root)$data) == nn + 1 && file.exists(vf) && identical(library_dataset("My VCF", root)$recipe$source$type, "vcf_file")) })

# ---- Round 75e: narration while the method scripts run; the watched run gives the same results as the plain one ----
local({ said <- list(); live_reset()
  pi <- pair_input("expA", "outA"); m <- c("ivw_mre", "egger", "wm")
  ha <- studio_new("Plain run"); ra <- studio_run(ha, pi$exposure_label, pi$outcome_label, pi$input, m, idx = 1L)
  options(mrstudio.say = function(msg, why) said[[length(said) + 1L]] <<- list(msg = msg, why = why), mrstudio.script_progress = function(dir, cn) live_scripts(dir, cn))
  hb <- studio_new("Watched run"); rb <- tryCatch(studio_run(hb, pi$exposure_label, pi$outcome_label, pi$input, m, idx = 1L), finally = options(mrstudio.say = NULL, mrstudio.script_progress = NULL))
  ea <- .est_df(pair_payload(ha$dir, ra$pair_folder)); eb <- .est_df(pair_payload(hb$dir, rb$pair_folder)); msgs <- vapply(said, `[[`, "", "msg")
  expect("Round 75e: a run watched for narration gives exactly the same estimates and replay result as the plain run",
         identical(ea$estimate, eb$estimate) && identical(ea$se, eb$se) && identical(pair_payload(hb$dir, rb$pair_folder)$replay$status, "identical"))
  expect("Round 75e: the Run screen is told when the method scripts start and when the replay check starts, each with why (fast methods finish before the first look)",
         any(grepl("^Running 3 method scripts", msgs)) && any(grepl("^Replay check: rerunning", msgs)) && all(nzchar(vapply(said[grepl("^Running 3 method scripts|^Replay check: rerunning", msgs)], `[[`, "", "why"))))   # Round 76: only the two start messages (a per-method replay line, said when the first look comes early, has no why by design)
  # a slow method: reported as still running after 45 s, once; a failed one is named
  d <- file.path(tempfile("live"), "pair"); dir.create(file.path(d, "02_scripts"), recursive = TRUE); dir.create(file.path(d, "04_tables")); dir.create(file.path(d, "03_console"))
  file.create(file.path(d, "02_scripts", c("00_run_all.R", "01_ivw_mre.R", "02_presso.R", "03_wm.R"))); said <<- list(); said <- list(); live_reset()
  options(mrstudio.say = function(msg, why) said[[length(said) + 1L]] <<- list(msg = msg, why = why))
  live_scripts(d); file.create(file.path(d, "04_tables", "ivw_mre_estimates.csv")); k <- paste(d, "run_all_console.txt"); .live[[k]]$t <- Sys.time() - 10; live_scripts(d)
  .live[[k]]$t <- Sys.time() - 50; live_scripts(d); live_scripts(d); writeLines("02_presso.R", file.path(d, "03_console", "failed_scripts.txt")); live_scripts(d)
  options(mrstudio.say = NULL); msgs <- vapply(said, `[[`, "", "msg")
  expect("Round 75e: a finished method is reported with its time, a slow one as still running (once, with why it is slow), a failed one by name",
         any(grepl("finished in 10 s \\(1 of 3\\)", msgs)) && sum(grepl("still running", msgs)) == 1L && grepl("refits the model many times", said[[which(grepl("still running", msgs))]]$why) &&
         any(grepl("MR-PRESSO.*stopped with an error", msgs)) && any(grepl("^Method 3 of 3", msgs))) })

# ---- Round 77: files from other sources (Excel, preambles, decimal commas, MAF and -log10 P headings, derived SE, Z scores, GRCh38) ----
local({
  td <- tempfile("imp77"); dir.create(td); set.seed(77)
  k <- 40; snp <- sprintf("rs%d", 5000 + seq_len(k)); ea <- rep(c("A", "C", "G", "T"), 10); oa <- rep(c("G", "T", "A", "C"), 10)
  b <- round(stats::runif(k, -0.1, 0.1), 6); se <- round(stats::runif(k, 0.01, 0.02), 6); p <- 2 * stats::pnorm(-abs(b / se)); f <- round(stats::runif(k, 0.05, 0.95), 4)
  mk <- function(name, lines) { fp <- file.path(td, name); writeLines(lines, fp); fp }
  head_of <- function(fp) { h <- import_read_head(fp); s <- import_suggest(h$columns, h$rows); list(h = h, s = s, hint = import_hints(h$columns, h$rows, s), map = lapply(s, `[[`, "column")) }
  facts0 <- list(name = "x", role = "either", kind = "continuous", scale = "beta", ea_ok = "yes", frq = "eaf", strand = "unknown", build = "37", ancestry = "EUR", n = "100000", scope = "selected", independent = "yes")
  # MAF heading: matched to the frequency column but marked as a minor allele frequency; an EAF heading in the same file wins
  a <- head_of(mk("maf.tsv", c("SNP\tA1\tA2\tMAF\tBETA\tSE\tP", sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s", snp, ea, oa, pmin(f, 1 - f), b, se, p))))
  a2 <- head_of(mk("both.tsv", c("SNP\tA1\tA2\tMAF\tEAF\tBETA\tSE\tP", sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", snp, ea, oa, pmin(f, 1 - f), f, b, se, p))))
  expect("Round 77: a MAF heading is marked as a minor allele frequency (not used as the effect allele frequency); an EAF heading wins over it",
         identical(a$hint$frq, "maf") && identical(a2$map$eaf, "EAF") && identical(a2$hint$frq, "eaf"))
  # -log10 P: hinted, converted, and a plain P column with values above 1 is refused with the reason
  ml <- mk("mlog.tsv", c("SNP\tA1\tA2\tBETA\tSE\tLOG10P", sprintf("%s\t%s\t%s\t%s\t%s\t%.10g", snp, ea, oa, b, se, -log10(p))))
  c1 <- head_of(ml); v1 <- import_validate(ml, c1$h$sep, c1$h$columns, c1$map, modifyList(facts0, list(p_scale = c1$hint$p_scale, frq = "none")))
  v2 <- import_validate(ml, c1$h$sep, c1$h$columns, c1$map, modifyList(facts0, list(p_scale = "p", frq = "none")))
  expect("Round 77: a -log10(P) column is recognised and converted; read as plain P it is refused with what to choose instead",
         identical(c1$hint$p_scale, "mlog10") && identical(c1$map$p, "LOG10P") && v1$ok && any(grepl("converted from -log10", unlist(v1$changes))) && !v2$ok && any(grepl("above 1", v2$errors)))
  # SE from the effect and P (Altman and Bland 2011), including a P value too small for double precision given as -log10(P)
  nose <- mk("nose.csv", c("rsid,effect_allele,other_allele,beta,pval", sprintf("%s,%s,%s,%s,%.17g", snp, ea, oa, b, p)))
  c2 <- head_of(nose); fa <- modifyList(facts0, list(name = "Derived SE", frq = "none", p_scale = "p"))
  r2 <- import_save(c2$h$path, "nose.csv", c2$h$sep, c2$h$columns, c2$map, fa); d2 <- lib_table_read(file.path(library_dir(), "Derived SE"), "standardised")
  lp <- -400 * log(10); zz <- stats::qnorm(lp - log(2), lower.tail = FALSE, log.p = TRUE)
  tiny <- .import_transform(data.frame(ea = "A", oa = "G", beta = "0.5", p = "400", stringsAsFactors = FALSE), list(p_scale = "mlog10", scale = "beta"))$d
  expect("Round 77: the SE is derived from the effect and its P value (exact for unrounded P; also for P below 1e-300 given as -log10 P); the change is recorded",
         isTRUE(all.equal(d2$se[match(snp, d2$snp)], se, tolerance = 1e-8)) && any(grepl("Altman and Bland", unlist(r2$changes))) && isTRUE(all.equal(tiny$se, 0.5 / zz)) && is.finite(tiny$se))
  # Z score with N and frequency (Zhu et al. 2016): b = z / sqrt(2f(1-f)(n + z^2)), se = 1 / sqrt(...)
  n <- 50000; zs <- stats::rnorm(k, 0, 4); den <- sqrt(2 * f * (1 - f) * (n + zs^2))
  zf <- mk("z.txt", c("ID REF ALT ALT_FREQ ZSCORE N", sprintf("%s %s %s %s %.12g %d", snp, oa, ea, f, zs, n)))
  c3 <- head_of(zf); m3 <- c3$map
  r3 <- import_save(c3$h$path, "z.txt", c3$h$sep, c3$h$columns, m3, modifyList(facts0, list(name = "From Z", scale = "z", frq = "eaf", n = "")))
  d3 <- lib_table_read(file.path(library_dir(), "From Z"), "standardised"); i3 <- match(snp, d3$snp)
  vb <- import_validate(c3$h$path, c3$h$sep, c3$h$columns, m3, modifyList(facts0, list(scale = "z", kind = "binary")))
  expect("Round 77: effect sizes and SEs are derived from Z scores with N and the frequency (Zhu 2016); ALT is the effect allele; refused for a binary trait",
         identical(m3$ea, "ALT") && identical(m3$oa, "REF") && identical(m3$z, "ZSCORE") && identical(c3$hint$scale, "z") && isTRUE(all.equal(d3$beta[i3], zs / den, tolerance = 1e-9)) && isTRUE(all.equal(d3$se[i3], 1 / den, tolerance = 1e-9)) &&
         !vb$ok && any(grepl("binary trait", vb$errors)))
  # Preamble, byte-order mark, semicolons with decimal commas, lower-case alleles, chr-prefixed chromosomes
  pre <- mk("pre.txt", c("﻿## produced by some pipeline", "## build GRCh37", "", "#CHROM;POS;ID;EA;NEA;BETA;SE", sprintf("chr1;%d;%s;%s;%s;%s;%s", 1000 + seq_len(k), snp, tolower(ea), tolower(oa), sub(".", ",", format(b, scientific = FALSE), fixed = TRUE), sub(".", ",", format(se, scientific = FALSE), fixed = TRUE))))
  c4 <- head_of(pre); v4 <- import_validate(c4$h$path, c4$h$sep, c4$h$columns, c4$map, modifyList(facts0, list(frq = "none")))
  d4 <- .import_transform(.import_read(c4$h$path, c4$h$sep, c4$h$columns, c4$map), facts0)$d
  expect("Round 77: comment lines, a byte-order mark, a #-heading, semicolons and decimal commas are handled; alleles upper-cased; chr1 written as 1",
         identical(c4$h$columns, c("CHROM", "POS", "ID", "EA", "NEA", "BETA", "SE")) && v4$ok && v4$stats$kept == k && isTRUE(all.equal(d4$beta, b)) && all(d4$chr == "1") && all(d4$ea == ea) && length(c4$h$notes) >= 3)
  # Excel: the results sheet of a workbook, with title, note and footnote rows
  x <- head_of(file.path(args_dir, "test_fixtures", "messy_gwas.xlsx"))
  vx <- import_validate(x$h$path, x$h$sep, x$h$columns, x$map, modifyList(facts0, list(p_scale = "p")))
  expect("Round 77: an Excel workbook is read from its results sheet (title, note and footnote rows left out), every column matched and 30 variants kept",
         identical(x$h$sheet, "Results") && identical(unname(unlist(x$map[c("snp", "ea", "oa", "eaf", "beta", "se", "p")])), c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value")) && vx$ok && vx$stats$kept == 30L && length(x$h$notes) == 3L)
  # GRCh38 positions: APOE rs7412 and rs429358 (GRCh37 19:45412079 and 19:45411941); a reverse-strand block complements the alleles
  bl <- liftover_blocks(); rb <- bl[bl$qstrand == "-" & bl$tchr == "1", ][1, ]
  g38 <- mk("g38.tsv", c("CHR\tBP\tA1\tA2\tBETA\tSE", "19\t44908822\tT\tC\t0.1\t0.01", "19\t44908684\tC\tT\t0.2\t0.01", sprintf("1\t%d\tA\tG\t0.3\t0.01", rb$tstart + 6), "1\t1\tA\tG\t0.4\t0.01"))
  old5 <- Sys.getenv("MRGUI_LD_DIR"); Sys.setenv(MRGUI_LD_DIR = tempfile("nopanel"))   # Round 79: without a panel, so the positions stay the IDs (no rsIDs added)
  c5 <- head_of(g38); Sys.setenv(MRGUI_LD_DIR = tempfile("nopanel")); r5 <- import_save(c5$h$path, "g38.tsv", c5$h$sep, c5$h$columns, c5$map, modifyList(facts0, list(name = "GRCh38 file", build = "38", frq = "none")))
  Sys.setenv(MRGUI_LD_DIR = old5)
  d5 <- lib_table_read(file.path(library_dir(), "GRCh38 file"), "standardised"); u5 <- import_validate(c5$h$path, c5$h$sep, c5$h$columns, c5$map, modifyList(facts0, list(build = "unknown")))
  expect("Round 77: GRCh38 positions are converted to GRCh37 (APOE variants at their known positions), reverse-strand alleles complemented, an unmappable one dropped, the chain recorded",
         identical(d5$pos[1:2], c(45412079, 45411941)) && identical(d5$ea[3], "T") && identical(d5$oa[3], "C") && nrow(d5) == 3L && identical(r5$liftover$sha256, LIFTOVER_CHAIN$sha256) && identical(r5$build_saved, "GRCh37") &&
         !u5$ok && any(grepl("genome build must be known", u5$errors)))
  # Frequencies from the 1000 Genomes panel (synthetic fixture panel), then Z scores to effects
  old_ld <- Sys.getenv("MRGUI_LD_DIR"); Sys.setenv(MRGUI_LD_DIR = file.path(args_dir, "test_fixtures", "ld_panel")); on.exit(Sys.setenv(MRGUI_LD_DIR = old_ld), add = TRUE)   # the synthetic fixture panel
  bim <- utils::read.table(file.path(Sys.getenv("MRGUI_LD_DIR"), "EUR.bim"), stringsAsFactors = FALSE)
  pz <- mk("pz.tsv", c("SNP\tEA\tOA\tZ\tN", sprintf("%s\t%s\t%s\t%s\t%d", bim$V2, bim$V5, bim$V6, round(stats::rnorm(nrow(bim), 0, 3), 4), 20000)))
  c6 <- head_of(pz); r6 <- tryCatch(import_save(c6$h$path, "pz.tsv", c6$h$sep, c6$h$columns, c6$map, modifyList(facts0, list(name = "Panel freq", scale = "z", frq = "none", eaf_fill = "panel", n = ""))), error = function(e) conditionMessage(e))
  d6 <- tryCatch(lib_table_read(file.path(library_dir(), "Panel freq"), "standardised"), error = function(e) NULL)
  expect("Round 77: missing frequencies are taken from the 1000 Genomes panel by rsID and allele, and then used to derive effects from Z scores",
         is.list(r6) && !is.null(d6) && nrow(d6) >= 1L && all(d6$eaf > 0 & d6$eaf < 1) && any(grepl("1000 Genomes EUR panel", unlist(r6$changes))) && isTRUE(all.equal(d6$se, 1 / sqrt(2 * d6$eaf * (1 - d6$eaf) * (20000 + (d6$beta / d6$se)^2)))))
})

# ---- Round 79: layouts of real sources (GWAS Catalog harmonised, FinnGen, Neale UKB, PLINK 2, SAIGE, papers' tables) ----
local({
  td <- tempfile("imp79"); dir.create(td); mk <- function(name, lines, eol = "\n") { fp <- file.path(td, name); con <- if (grepl("gz$", name)) gzfile(fp, "w") else file(fp, "wb"); writeLines(lines, con, sep = eol); close(con); fp }
  hs <- function(fp) { h <- import_read_head(fp); s <- import_suggest(h$columns, h$rows); list(h = h, s = s, map = lapply(s, `[[`, "column"), hint = import_hints(h$columns, h$rows, s, build_check = FALSE)) }
  tf <- function(x, facts = list()) { d <- .import_read(x$h$path, x$h$sep, x$h$columns, x$map); .import_transform(d, utils::modifyList(list(scale = "beta", p_scale = "p"), facts)) }
  # GWAS Catalog harmonised: hm_ prefixes; hm_beta wins over hm_odds_ratio; hm_rsid over variant_id
  a <- hs(mk("h.tsv", c("hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\thm_ci_lower\thm_ci_upper\thm_effect_allele_frequency\tvariant_id\tp_value\tstandard_error",
                         "1_100_A_G\trs1\t1\t100\tA\tG\t0.1\tNA\tNA\tNA\t0.3\trs1\t0.01\t0.02", "1_200_C_T\trs2\t1\t200\tC\tT\t-0.2\tNA\tNA\tNA\t0.4\trs2\t0.02\t0.03", "1_300_G_A\trs3\t1\t300\tG\tA\t0.05\tNA\tNA\tNA\t0.5\trs3\t0.3\t0.04")))
  expect("Round 79: GWAS Catalog harmonised headings (hm_) are recognised: effect hm_beta (not hm_odds_ratio), rsID hm_rsid, alleles and frequency",
         identical(unlist(a$map[c("snp", "chr", "pos", "ea", "oa", "beta", "eaf", "se", "p")]), c(snp = "hm_rsid", chr = "hm_chrom", pos = "hm_pos", ea = "hm_effect_allele", oa = "hm_other_allele", beta = "hm_beta", eaf = "hm_effect_allele_frequency", se = "standard_error", p = "p_value")) && is.null(a$hint$scale))
  # SAIGE: the effect allele is Allele2
  b <- hs(mk("saige.txt", c("CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 N BETA SE Tstat p.value", "1 100 rs1 A G 10 0.3 1000 0.1 0.02 5 1e-6", "1 200 rs2 C T 20 0.4 1000 -0.2 0.03 -6 1e-8", "1 300 rs3 G A 30 0.5 1000 0.05 0.04 1 0.2")))
  expect("Round 79: SAIGE output: Allele2 (whose frequency is AF_Allele2) is proposed as the effect allele, Allele1 as the other", identical(b$map$ea, "Allele2") && identical(b$map$oa, "Allele1") && identical(b$map$eaf, "AF_Allele2"))
  # PLINK 2: A1 can be REF; covariate rows are left out
  c2 <- hs(mk("p.glm.logistic", c("#CHROM\tPOS\tID\tREF\tALT\tA1\tA1_FREQ\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP", "1\t100\trs1\tA\tG\tG\t0.3\tADD\t1000\t1.2\t0.05\t3.6\t3e-4", "1\t200\trs2\tC\tT\tC\t0.4\tADD\t1000\t0.8\t0.06\t-3.7\t2e-4",
                                   "1\t200\trs2\tC\tT\tC\t0.4\tAGE\t1000\t1.01\t0.01\t1\t0.3", "1\t300\trs3\tG\tA\tA\t0.5\tADD\t1000\t1.05\t0.04\t1.2\t0.2")))
  t2 <- tf(c2, list(scale = c2$hint$scale))
  expect("Round 79: PLINK 2 --glm: A1 is the effect allele, the other allele is REF or (where A1 is REF) ALT; OR and LOG(OR)_SE read as log OR and SE; covariate rows left out",
         identical(c2$map$ea, "A1") && identical(c2$map$oa_alt, "ALT") && identical(c2$map$se, "LOG(OR)_SE") && identical(c2$hint$scale, "or") && nrow(t2$d) == 3L && identical(t2$d$oa, c("A", "T", "G")) && isTRUE(all.equal(t2$d$beta, log(c(1.2, 0.8, 1.05)))) && any(grepl("AGE", t2$changes)))
  # FinnGen: several rsIDs in one cell
  f <- hs(mk("fg.tsv.gz", c("#chrom\tpos\tref\talt\trsids\tpval\tmlogp\tbeta\tsebeta\taf_alt", "1\t100\tA\tG\trs1,rs99\t1e-6\t6\t0.1\t0.02\t0.3", "1\t200\tC\tT\trs2\t1e-8\t8\t-0.2\t0.03\t0.4", "1\t300\tG\tA\trs3\t0.2\t0.7\t0.05\t0.04\t0.5")))
  expect("Round 79: FinnGen: alt is the effect allele, af_alt its frequency, and a cell with several rsIDs keeps the first", identical(f$map$ea, "alt") && identical(f$map$eaf, "af_alt") && identical(tf(f)$d$snp, c("rs1", "rs2", "rs3")))
  # Neale lab UK Biobank: .bgz, variant chr:pos:ref:alt only
  n <- hs(mk("neale.tsv.bgz", c("variant\tminor_allele\tminor_AF\tn_complete_samples\tbeta\tse\ttstat\tpval", "1:100:A:G\tG\t0.3\t1000\t0.1\t0.02\t5\t1e-6", "1:200:C:T\tT\t0.4\t1000\t-0.2\t0.03\t-6\t1e-8", "1:300:G:A\tA\t0.45\t1000\t0.05\t0.04\t1\t0.2")))
  expect("Round 79: Neale lab files (.bgz, variant chr:pos:ref:alt): chromosome, position, REF (other) and ALT (effect) are read from the variant; it is not taken as an rsID",
         identical(unlist(n$map[c("chr", "pos", "ea", "oa")]), c(chr = "CHR", pos = "POS", ea = "ALT", oa = "REF")) && is.null(n$map$snp) && identical(tf(n)$d$pos, c(100, 200, 300)) && any(grepl("chr:pos:ref:alt", unlist(n$h$notes))))
  # A paper's table: combined columns and typeset P values
  sc <- .split_combined(c("Lead SNP", "Chr:position (hg19)", "Effect/other allele", "OR (95% CI)", "Beta (SE)", "P"),
                        data.frame(a = c("rs1", "rs2"), b = c("1:1,000", "2:2000"), c = c("t/c", "A/G"), d = c("1.12 (1.05–1.19)", "0.90 (0.80 to 1.01)"), e = c("−0.15 (0.01)", "0.2 (0.03)"), f = c("3.2 × 10−8", "0.04"), stringsAsFactors = FALSE))
  expect("Round 79: combined columns of published tables are split (chr:pos, effect/other allele, OR with its 95% CI, beta with its SE), and typeset P values are read",
         identical(sc$head, c("Lead SNP", "Chromosome", "Position (hg19)", "Effect allele", "Other allele", "OR", "OR CI lower", "OR CI upper", "Beta", "SE", "P")) && identical(sc$body[[3]], c("1000", "2000")) &&
         identical(sc$body[[4]], c("T", "A")) && identical(sc$body[[8]], c("1.19", "1.01")) && identical(sc$body[[9]], c("-0.15", "0.2")) && isTRUE(all.equal(.num_parse(c("3.2 × 10−8", "1.2x10^-5", "<1e-300", "5E-8")), c(3.2e-8, 1.2e-5, 1e-300, 5e-8))))
  # SE from a 95% CI (log scale for OR), and per-row SE from P where the SE cell is missing
  ci <- .import_transform(data.frame(ea = "A", oa = "G", beta = c("1.2", "0.8"), ci_lower = as.character(exp(log(c(1.2, 0.8)) - 1.959964 * 0.05)), ci_upper = as.character(exp(log(c(1.2, 0.8)) + 1.959964 * 0.05)), stringsAsFactors = FALSE), list(scale = "or"))$d
  pr <- .import_transform(data.frame(ea = "A", oa = "G", beta = c("0.1", "0.2"), se = c("0.02", "."), p = c("1e-6", as.character(2 * stats::pnorm(-0.2 / 0.04))), stringsAsFactors = FALSE), list(scale = "beta", p_scale = "p"))
  expect("Round 79: the SE is derived from a 95% confidence interval (on the log scale for odds ratios), and for single rows whose SE is missing, from the P value",
         isTRUE(all.equal(ci$se, c(0.05, 0.05), tolerance = 1e-6)) && isTRUE(all.equal(ci$beta, log(c(1.2, 0.8)))) && isTRUE(all.equal(pr$d$se, c(0.02, 0.04))) && any(grepl("1 variant", pr$changes)))
  w <- .import_transform(data.frame(ea = "A", oa = "G", beta = c("1.1", "0.95", "1.2", "1.05", "0.9"), se = "0.02", stringsAsFactors = FALSE), list(scale = "beta"))
  expect("Round 79: odds ratios in a column read as beta are flagged (all positive and centred on 1)", any(grepl("as odds ratios are", w$warnings)))
  # The 1000 Genomes panel (fixture): build detected from rsID positions (GRCh37, or GRCh38 converted), rsIDs added by position and alleles
  options(mrgui.engine_dir = args_dir); old79 <- Sys.getenv("MRGUI_LD_DIR"); Sys.setenv(MRGUI_LD_DIR = file.path(args_dir, "test_fixtures", "ld_panel")); on.exit(Sys.setenv(MRGUI_LD_DIR = old79), add = TRUE)   # the fixture panel
  b37 <- import_build_guess(list(list(ID = "rs9100001", CHR = "1", POS = "1000000"), list(ID = "rs9100002", CHR = "1", POS = "1000100"), list(ID = "rs9200001", CHR = "1", POS = "15000000"), list(ID = "rs9300001", CHR = "1", POS = "40000000")),
                            list(snp = list(column = "ID"), chr = list(column = "CHR"), pos = list(column = "POS")))
  bl <- liftover_blocks(); inv <- function(p37) { bb <- bl[bl$qchr == "1" & bl$qstrand == "+", ]; k <- which(p37 - 1 >= bb$qstart & p37 - 1 < bb$qstart + bb$size)[1]; bb$tstart[k] + (p37 - 1 - bb$qstart[k]) + 1 }
  b38 <- import_build_guess(lapply(c(1000000, 1000100, 15000000, 40000000), function(p) list(ID = c("rs9100001", "rs9100002", "rs9200001", "rs9300001")[match(p, c(1000000, 1000100, 15000000, 40000000))], CHR = "1", POS = as.character(inv(p)))),
                            list(snp = list(column = "ID"), chr = list(column = "CHR"), pos = list(column = "POS")))
  rr <- .rsid_from_panel(data.frame(chr = c("1", "1", "1"), pos = c(1000000, 1000100, 5), ea = c("C", "A", "A"), oa = c("A", "G", "C"), stringsAsFactors = FALSE))
  expect("Round 79: the genome build is detected from where rsIDs sit in the 1000 Genomes panel (GRCh37 directly, GRCh38 after conversion), and position-only variants get rsIDs by position and alleles",
         identical(b37$build, "37") && identical(b38$build, "38") && identical(rr$snp, c("rs9100001", NA, NA)) && rr$n == 1L)
})

# ---- Round 80: narration per kind of wait, in every kind of study; every finished method reported; OpenGWAS populations ----
local({
  # fast methods that finish between two looks are all reported
  d <- file.path(tempfile("live80"), "pair"); dir.create(file.path(d, "02_scripts"), recursive = TRUE); dir.create(file.path(d, "04_tables")); dir.create(file.path(d, "03_console"))
  file.create(file.path(d, "02_scripts", c("00_run_all.R", "01_ivw_mre.R", "02_egger.R", "03_wm.R", "04_mode.R", "05_presso.R"))); said <- list(); live_reset()
  options(mrstudio.say = function(msg, why) said[[length(said) + 1L]] <<- list(msg = msg, why = why))
  live_scripts(d); file.create(file.path(d, "04_tables", paste0(c("ivw_mre", "egger", "wm"), "_estimates.csv"))); live_scripts(d)
  file.create(file.path(d, "04_tables", "mode_estimates.csv")); writeLines("05_presso.R", file.path(d, "03_console", "failed_scripts.txt")); live_scripts(d)
  options(mrstudio.say = NULL); msgs <- vapply(said, `[[`, "", "msg")
  expect("Round 80: every method that finished between two looks is reported (three at once, then the next one and a failed one by name)",
         any(grepl("^Inverse-variance weighted.* finished in .*\\(1 of 5\\)", msgs)) && any(grepl("^Finished: MR-Egger, Weighted median \\(3 of 5\\)", msgs)) && any(grepl("^Weighted mode finished in", msgs)) && any(grepl("MR-PRESSO stopped with an error.*\\(5 of 5\\)", msgs)))
  # every engine action that narrates has its own box on the screen, the engine names the action in every message and marks its end
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n"); app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  grp <- regmatches(js, regexpr("const LIVE_GROUP = \\{[^}]*\\}", js)); keys <- regmatches(grp, gregexpr("[a-z_]+(?=: ')", grp, perl = TRUE))[[1]]
  acts <- c("facts", "mv_prepare", "med_prepare", "bidir_prepare", "run_all", "build_package", "ogw_dl_exposure", "ogw_dl_outcome", "run_sens", "run_outliers", "mediation", "bidir_cml", "instruments", "run_ties", "compare_default", "init", "diagnostics")
  expect("Round 80: each action that narrates (two-sample, multivariable, mediation and bidirectional studies, reruns, outliers, publishing, downloads) has its own box group; messages name their action and every action's end is sent",
         all(acts %in% keys) && grepl('act = getOption("mrstudio.act", "")', app, fixed = TRUE) && grepl('ev("live_end", act = t)', app, fixed = TRUE) && grepl('bg_start(rq$dir, "run_pair", a, act = "run_all")', app, fixed = TRUE) && !grepl("S\\.live\\b", js))
  # the screen's boxes, run in node: a run opens its own box (a new clock) after the Check box closed; news from another action does not land in it
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const cst = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, i + src.slice(i).search(/;\\r?\\n/) + 1); };",
      "global.document = { getElementById: () => null }; global.S = { lives: {}, _pin: {}, _top: {} }; let now = 1000; Date.now = () => now; global.setTimeout = () => 0; global.clearTimeout = () => 0;",
      "eval(cst('LIVE_GROUP') + cst('LIVE_PREP') + cst('mmss') + grab('liveOpen') + grab('liveEnd') + grab('liveAdd') + grab('liveFollow') + grab('liveTick'));",
      "liveOpen('check'); liveAdd('Pair 1 of 2', '', 'facts'); now = 200000; liveEnd('check');",
      "now = 380000; liveAdd('Running 31 method scripts', '', 'run_all'); now = 381000; liveAdd('Instrument review: loading', '', 'instruments');",
      "console.log([S.lives.run.t0, S.lives.run.lines.length, S.lives.check.open, S.lives.check.lines.length, S.lives.review.lines.length].join(' '));"), f)
    o <- strsplit(system2("node", f, stdout = TRUE), " ")[[1]]
    expect("Round 80: on the screen, a run gets its own box and clock (it showed the Check step's clock), the Check box stays closed, and news from another action goes to that action's box",
           identical(o, c("380000", "1", "false", "1", "1")))
  }
  # populations: the GWAS's own (OpenGWAS) in the Check, the ancestry rule and the draft text, not the panel used for clumping
  r <- list(source = list(type = "opengwas", population = "Mixed"), facts = list(ancestry = "EUR")); r2 <- list(source = list(type = "opengwas", population = "European"), facts = list(ancestry = "EUR"))
  expect("Round 80: an OpenGWAS download's population is OpenGWAS's own (Mixed stays Mixed although its instruments were clumped with the EUR panel; European is EUR)",
         identical(.gwas_pop(r), "Mixed") && identical(.gwas_pop(r2), "EUR") && is.null(.gwas_pop(list(source = list(type = "file"), facts = list(ancestry = "EUR")))))
})

# ---- Round 81: LD-based methods with the panel of the GWAS's ancestry; Mixed data only with the analyst's accepted choice ----
local({
  pd <- tempfile("panel81"); ld <- tempfile("ldsc81"); dir.create(pd); file.create(file.path(pd, c("EUR.bed", "EAS.bed")))
  full <- function(d) { dir.create(d, recursive = TRUE); file.create(file.path(d, c(paste0(1:22, ".l2.ldscore.gz"), paste0(1:22, ".l2.M_5_50")))) }
  full(file.path(ld, "eur_w_ld_chr")); writeLines("SNP A1 A2", file.path(ld, "w_hm3.snplist"))
  rc <- function(tr, anc, pop = NULL) list(trait = tr, facts = list(ancestry = anc), source = if (is.null(pop)) list(type = "file") else list(type = "opengwas_vcf", population = pop))
  eu <- rc("LDL", "EUR", "European"); eu2 <- rc("CHD", "EUR"); mx <- rc("CHD mixed", "", "Mixed"); mx2 <- rc("LDL mixed", "", "Mixed"); ea <- rc("BMI EAS", "EAS"); ea2 <- rc("T2D EAS", "EAS")
  r1 <- ldref_resolve(list(eu, eu2), "mrlap", NULL, ld, pd)
  expect("Round 81: two European GWAS use the LDSC authors' European LD scores, with nothing to accept and no Limitations sentence",
         identical(r1$source, "official") && is.null(r1$accepted) && !length(ldref_limitation(list(method = "mrlap", reference = r1$what, accepted = NULL), c("LDL", "CHD"))))
  i2 <- ldref_info(list(eu, mx), "cause", ld, pd); e2 <- tryCatch(ldref_resolve(list(eu, mx), "cause", NULL, ld, pd), error = conditionMessage)
  r2 <- ldref_resolve(list(eu, mx), "cause", list(pop = "EAS", ok = TRUE), ld, pd)
  expect("Round 81: a GWAS labelled Mixed needs the analyst's acceptance; with one European GWAS the European panel is used whatever was picked, and the acceptance is dated",
         isTRUE(i2$needs_accept) && i2$fixed && grepl("labelled Mixed", e2) && grepl("I accept this and take responsibility for it", e2, fixed = TRUE) &&
         identical(r2$pop, "EUR") && identical(basename(r2$bfile), "EUR") && isTRUE(r2$accepted$ok) && grepl("^\\d{4}-\\d\\d-\\d\\d", r2$accepted$at))
  r3 <- ldref_resolve(list(mx2, mx), "cause", list(pop = "EAS", ok = TRUE), ld, pd)
  e3 <- tryCatch(ldref_resolve(list(mx2, mx), "ldsc", list(pop = "EAS", ok = TRUE), ld, pd), error = conditionMessage)
  e3b <- tryCatch(ldref_resolve(list(mx2, mx), "cause", list(pop = "AFR", ok = TRUE), ld, pd), error = conditionMessage)
  expect("Round 81: for two Mixed GWAS the analyst's panel is used (East Asian for CAUSE); LD scores not computed yet and a panel not installed are refused with what to do",
         identical(r3$pop, "EAS") && identical(basename(r3$bfile), "EAS") && grepl("LD scores for EAS have not been computed yet.*about 30 minutes", e3) && grepl("AFR panel; install that panel first", e3b))
  full(file.path(ld, "computed", "EAS")); r4 <- ldref_resolve(list(ea, ea2), "mrapss", NULL, ld, pd)
  l4 <- ldref_limitation(list(method = "mrapss", reference = r4$what, accepted = NULL), c("BMI EAS", "T2D EAS"))
  expect("Round 81: two East Asian GWAS use LD scores computed from the East Asian panel and that panel, without acceptance; the Limitations say the scores were computed by MR Studio",
         identical(r4$source, "EAS") && grepl("computed/EAS$", r4$scores) && identical(basename(r4$bfile), "EAS") && is.null(r4$accepted) &&
         grepl("computed by MR Studio from the 1000 Genomes Phase 3 EAS panel \\(504 people\\)", l4) && grepl("less precise", l4))
  ref <- function(a, b, m) ldref_info(list(a, b), m, ld, pd)$refused
  expect("Round 81: different ancestries, an ancestry without a panel, and LHC-MR with Mixed or non-European data are refused with the reason; LHC-MR on two European GWAS is not",
         grepl("different ancestries \\(EUR and EAS\\)", ref(eu, ea, "ldsc")) && grepl("Middle Eastern ancestry, which has no 1000 Genomes panel", ref(eu, rc("X", "", "Middle Eastern"), "cause")) &&
         grepl("UK10K", ref(eu, mx, "lhcmr")) && grepl("UK10K", ref(ea, ea2, "lhcmr")) && !nzchar(ref(eu, eu2, "lhcmr")) && !nzchar(ref(eu, mx, "mrlap")))
  # saved beside the result; the draft's Limitations sentence names the datasets, the panel, the date and what may be biased
  h <- list(dir = tempfile("study81")); dir.create(h$dir); old <- get("manifest_register", globalenv()); assign("manifest_register", function(...) invisible(NULL), globalenv())
  on.exit(assign("manifest_register", old, globalenv()), add = TRUE)
  pr <- list(folder = "P01_LDL__to__CHD", exName = "LDL", outName = "CHD mixed"); ldref_save(h, pr, r2); back <- ldref_read(h, pr, "cause")
  lim <- ldref_limitation(back, c("LDL", "CHD mixed"))
  expect("Round 81: the LD reference and the acceptance are saved beside the result, and the Limitations sentence names the dataset labelled Mixed, the panel, the date and what may be biased",
         file.exists(file.path(h$dir, "04_across_pairs", "cause", "P01_LDL__to__CHD_ld_reference.json")) && identical(back$pop, "EUR") && isTRUE(back$accepted$ok) &&
         grepl("CAUSE for LDL and CHD mixed used the 1000 Genomes Phase 3 EUR panel", lim) && grepl("CHD mixed \\(Mixed\\)", lim) && grepl(substr(r2$accepted$at, 1, 10), lim, fixed = TRUE) && grepl("causal and sharing models", lim))
  # every method uses it: no European-only check is left in the four methods, and each saves its reference
  src <- paste(unlist(lapply(c("ldsc_core.R", "mrlap_core.R", "mrapss_core.R", "cause_core.R"), function(f) readLines(file.path(args_dir, f), warn = FALSE))), collapse = "\n")
  expect("Round 81: LDSC, MRlap, MR-APSS and CAUSE resolve their LD reference (no European-only check left) and save it with the result",
         !grepl("must be of European ancestry", src) && length(gregexpr("ldref_resolve\\(lapply\\(ds", src)[[1]]) == 4L && length(gregexpr("ldref_save\\(h, pair, lr\\)", src)[[1]]) == 4L)
  # Settings: the copy job's note follows the current pins (SAS showed "differs" while its files matched)
  pd2 <- tempfile("panel81b"); dir.create(pd2); file.create(file.path(pd2, paste0("SAS", c(".bed", ".bim", ".fam"))))
  writeLines(jsonlite::toJSON(list(state = "done", mode = "folder", pops = list("SAS"), pinned = FALSE), auto_unbox = TRUE), file.path(pd2, "download_status.json"))
  om <- get("ld_panel_matches_pin", globalenv()); assign("ld_panel_matches_pin", function(p, dir) TRUE, globalenv()); on.exit(assign("ld_panel_matches_pin", om, globalenv()), add = TRUE)
  s1 <- ld_panel_status(pd2); assign("ld_panel_matches_pin", function(p, dir) FALSE, globalenv()); s2 <- ld_panel_status(pd2)
  expect("Round 81: the panel copy note is judged with the current pins (it showed the flag written at install time)", isTRUE(s1$job$pinned) && identical(s2$job$pinned, FALSE))
  # the screen: a Mixed pair needs the tick before the Run button is offered; the panel used is shown beside the result
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const cst = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, i + src.slice(i).search(/;\\r?\\n/) + 1); };",
      "global.S = { ldscRef: { official: true, computed: [] } }; global.send = () => 0; global.esc = (x) => String(x);",
      "eval(cst('POPN') + grab('ldrefState') + grab('ldrefNote'));",
      "const av = { EUR: { panel: true, scores: true, official: true }, EAS: { panel: true, scores: false, official: false } };",
      "const R = { anc: { name: 'MRlap', pop: 'EUR', fixed: false, needs_accept: true, unknown: [{ trait: 'CHD', label: 'Mixed' }], refused: '', needs: { scores: true, panel: false }, avail: av, bias: 'the corrected effect' } }, pr = { folder: 'P01' };",
      "const a = ldrefState('mrlap', pr, R); S.ldacc['mrlap|P01'].ok = true; const b = ldrefState('mrlap', pr, R); S.ldacc['mrlap|P01'] = { pop: 'EAS', ok: true }; const c = ldrefState('mrlap', pr, R);",
      "const n = ldrefNote({ ldref: { reference: 'the LDSC authors\\' European LD scores', accepted: { at: '2026-11-21 19:30:00' } } });",
      "console.log(JSON.stringify([a.ok, /labelled Mixed/.test(a.html), /I accept this and take responsibility for it/.test(a.html), b.ok, b.accept, c.ok, /not been computed yet/.test(c.html), /2026-11-21/.test(n)]));"), f)
    o <- jsonlite::fromJSON(paste(system2("node", f, stdout = TRUE), collapse = ""), simplifyVector = FALSE)
    expect("Round 81: on the screen, a pair with Mixed data shows the warning and the tick, offers the run only after the tick, sends the panel and the acceptance, and shows the acceptance date beside the result",
           identical(o[[1]], FALSE) && isTRUE(o[[2]]) && isTRUE(o[[3]]) && isTRUE(o[[4]]) && identical(o[[5]], list(pop = "EUR", ok = TRUE)) && identical(o[[6]], FALSE) && isTRUE(o[[7]]) && isTRUE(o[[8]]))
  }
})

# Round 82: the European LD scores download in a background R process, with progress, a time limit per connection and resume
# (it ran inside the app: a stalled connection froze the app at 18 of 34 MB with no sign of progress)
local({
  pay <- tempfile(); set.seed(82); writeBin(as.raw(sample(0:255, 300000, TRUE)), pay); md5 <- unname(tools::md5sum(pay))
  srv <- tempfile(fileext = ".R"); pidf <- tempfile(); port <- httpuv::randomPort()
  writeLines(c("a <- commandArgs(TRUE); x <- readBin(a[1], 'raw', file.size(a[1])); writeLines(as.character(Sys.getpid()), a[3])",
    "hd <- list('Content-Type' = 'application/octet-stream')",
    "httpuv::runServer('127.0.0.1', as.integer(a[2]), list(call = function(req) { r <- if (is.null(req$HTTP_RANGE)) '' else req$HTTP_RANGE",
    "  if (nzchar(r) && !grepl('norange', req$QUERY_STRING)) { s <- as.numeric(sub('bytes=([0-9]+)-.*', '\\\\1', r)); list(status = 206L, headers = hd, body = x[(s + 1):length(x)]) }",
    "  else list(status = 200L, headers = hd, body = x) }))"), srv)
  rscript_background(c("--vanilla", shQuote(srv), shQuote(pay), port, shQuote(pidf)), log = tempfile())   # as the app starts R (a moved runtime needs RHOME; plain Rscript failed on the CI Macs)
  on.exit(if (file.exists(pidf)) tools::pskill(as.integer(readLines(pidf))), add = TRUE)
  url <- sprintf("http://127.0.0.1:%d/f.bin", port)
  up <- FALSE; for (i in 1:240) { if (isTRUE(tryCatch(curl::curl_fetch_memory(url)$status_code == 200, error = function(e) FALSE))) { up <- TRUE; break }; Sys.sleep(0.5) }
  dl <- function(...) tryCatch({ ldsc_download_resume(...); TRUE }, error = function(e) { cat("   download test:", conditionMessage(e), "\n"); FALSE })
  p1 <- tempfile(); ok1 <- dl(url, p1)
  p2 <- tempfile(); writeBin(readBin(pay, "raw", 100000), p2); first <- NA; ok2 <- dl(url, p2, progress = function(n) if (is.na(first)) first <<- n)
  p3 <- tempfile(); writeBin(as.raw(rep(7, 100000)), p3); ok3 <- dl(paste0(url, "?norange"), p3)
  e <- tryCatch(ldsc_download_resume(sprintf("http://127.0.0.1:%d/f.bin", httpuv::randomPort()), tempfile(), tries = 2L, wait = 0), error = function(e) conditionMessage(e))
  expect("Round 82: the LD score download continues where it stopped (Range), starts over when the server ignores Range, and stops with a clear message after its attempts",
         up && ok1 && ok2 && ok3 && identical(unname(tools::md5sum(p1)), md5) && identical(unname(tools::md5sum(p2)), md5) && isTRUE(first > 100000) && identical(unname(tools::md5sum(p3)), md5) &&
         grepl("after 2 attempts", e) && grepl("continue where it stopped", e))
  d <- tempfile(); dir.create(d); .ldsc_fetch_write(d, state = "running", pid = 999999L, done_bytes = 5, bytes = 10)
  s <- ldsc_fetch_read(d)
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n"); core <- paste(readLines(file.path(args_dir, "ldsc_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 82: the download runs in the background (the app no longer downloads it itself; LD score computation fetches its 5 MB list with the resumable download), and a job ended by closing the app reads as stopped",
         grepl("ldsc_fetch = { r <- tryCatch(ldsc_fetch_start()", app, fixed = TRUE) && !grepl("ldsc_ref_fetch()", app, fixed = TRUE) && !grepl('ldsc_ref_fetch("hm3", dir)\n', core, fixed = TRUE) &&
         identical(s$state, "stopped") && grepl("continue where it stopped", s$msg) && "fetch" %in% names(ldsc_status(d)))
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "global.esc = (x) => String(x); let sent = []; global.send = (a) => sent.push(a); let t = 0; global.setTimeout = (fn, ms) => { t = ms; return 1; }; global.clearTimeout = () => 0;",
      "eval(grab('ldscFetchUi') + grab('ldscPoll'));",
      "const a = ldscFetchUi({ fetch: { state: 'running', step: 'downloading', done_bytes: 9437184, bytes: 33357890 } }), b = ldscFetchUi({ fetch: { state: 'error', msg: 'The download stopped after 5 attempts' } }), c = ldscFetchUi({});",
      "global.S = { ldscRef: { fetch: { state: 'running' } } }; ldscPoll();",
      "console.log([/9\\.0 of 31\\.8 MB/.test(a) && !/button/.test(a), /Continue the download/.test(b) && /after 5 attempts/.test(b), /34 MB/.test(c), t].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE), collapse = "")
    expect("Round 82: on the screen, the download shows the megabytes received, offers to continue after a stop, and is followed every 2 seconds",
           identical(o, "true true true 2000"))
  }
})


# Round 83: long analyses run in a background process, so the Studies screen lists and opens studies while a study runs (the engine
# answered nothing for the minutes a pair took), a run stays with its own study, and the screen keeps a running study's state aside
local({
  h <- studio_new("R83 background"); DBI::dbDisconnect(h$con)
  old <- options(mrgui.engine_dir = args_dir); on.exit(options(old), add = TRUE)
  job <- bg_start(h$dir, "mediation", list(ci = 0.95), act = "mediation")
  run_now <- length(bg_running(h$dir)) == 1L && identical(studies_overview()[[which(vapply(studies_overview(), `[[`, "", "name") == "R83 background")]]$status, "RUNNING")
  for (i in 1:240) { r <- bg_poll(job, 0); if (isTRUE(r$done)) break; Sys.sleep(0.5) }
  ts <- vapply(r$lines, function(l) jsonlite::parse_json(l)$t, "")
  expect("Round 83: an analysis runs as a background job of its study (listed as RUNNING meanwhile); its messages arrive, as the screen would get them, and it ends done",
         run_now && isTRUE(r$done) && isTRUE(r$ok) && all(c("live", "mediation", "live_end") %in% ts) && !length(bg_running(h$dir)) &&
         identical(jsonlite::parse_json(r$lines[ts == "mediation"])$rows, list()))
  # a study opened while its job runs keeps the job's task running; with no job, a task left running by a crash goes back to pending
  h2 <- studio_open("R83 background"); db_upsert_task(h2$con, "t83", "method", "ivw_mre", "1"); DBI::dbExecute(h2$con, "UPDATE tasks SET state = 'running' WHERE task_id = 't83'"); DBI::dbDisconnect(h2$con)
  fake <- file.path(h$dir, ".engine", "jobs", "run_pair_fake"); dir.create(fake, recursive = TRUE); writeLines(jsonlite::toJSON(list(state = "running", pid = Sys.getpid()), auto_unbox = TRUE), file.path(fake, "status.json"))
  h3 <- studio_open("R83 background"); kept <- DBI::dbGetQuery(h3$con, "SELECT state FROM tasks WHERE task_id = 't83'")$state; DBI::dbDisconnect(h3$con)
  unlink(fake, recursive = TRUE); h4 <- studio_open("R83 background"); back <- DBI::dbGetQuery(h4$con, "SELECT state FROM tasks WHERE task_id = 't83'")$state; DBI::dbDisconnect(h4$con)
  gone <- file.path(h$dir, ".engine", "jobs", "x_gone"); dir.create(gone); file.create(file.path(gone, "events.jsonl")); writeLines('{"state":"running","pid":999999}', file.path(gone, "status.json")); g <- bg_poll(gone, 0)
  runner <- paste(readLines(file.path(dirname(args_dir), "..", "contracts", "r", "runner.R"), warn = FALSE), collapse = "\n")
  expect("Round 83: opening a study does not reset the running task of its background job (it still recovers tasks of a crash); a job whose process vanished ends with a message; the database waits for the other process's write",
         identical(kept, "running") && identical(back, "pending") && isTRUE(g$done) && any(grepl("stopped unexpectedly", g$lines)) &&
         regexpr("sqliteSetBusyHandler", runner) < regexpr("journal_mode=WAL", runner))
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n"); bodies <- paste(vapply(BG_ACTIONS, function(f) paste(deparse(f), collapse = "\n"), ""), collapse = "\n")
  expect("Round 83: the run of every kind of study (one pair at a time), the sensitivity, tie and outlier reanalyses, mediation, bidirectional MR-cML, LD score regression and the publication package run in the background; a run stays with its own study",
         setequal(names(BG_ACTIONS), c("run_pair", "run_sens", "run_outliers", "run_ties", "run_harm", "mediation", "bidir_cml", "build_package", "ldsc")) &&
         all(vapply(c("mediation", "bidir_cml", "run_sens", "run_outliers", "run_ties", "run_harm", "build_package"), function(k) grepl(sprintf('%s = bg_go("%s", a)', k, k), app, fixed = TRUE), TRUE)) &&
         grepl('if (identical(t, "ldsc")) return(bg_go("ldsc", a))', app, fixed = TRUE) && grepl('bg_start(rq$dir, "run_pair", a, act = "run_all")', app, fixed = TRUE) &&
         grepl("rq$dir <- rv$h$dir", app, fixed = TRUE) && grepl("A run is already going on in", app, fixed = TRUE) && !grepl("rv$h", bodies, fixed = TRUE) && !grepl("run_one(", app, fixed = TRUE))
  cz <- paste(readLines(file.path(args_dir, "cause_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 83: CAUSE's memory estimate covers the measured peak (7.35 GB for 8.5 million variants in the shipped runtime), and its script says it needs loo 2.8.0 before it starts (newer loo stops it at the very end)",
         CAUSE_MEMORY_GB(8479693) >= 7.35 && CAUSE_MEMORY_GB(8479693) < 8.5 && grepl("CAUSE needs loo 2.8.0", cz, fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const cst = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, i + src.slice(i).search(/;\\r?\\n/) + 1); };",
      "global.S = { lives: {} }; let sent = [], drawn = 0; global.prevOk = () => false; global.kindMethods = () => []; global.render = () => { if (!S._quiet) drawn++; }; global.liveTick = () => 0;",   # Round 90: a restored study restarts its live box
      "global.send = (t, p) => { if (S._quiet) { (S._deferred = S._deferred || []).push([t, p]); return; } sent.push(t); };",
      "global.handle = (m) => { if (m.t === 'pair_done') { S.run.done++; S.payloads[m.folder] = 1; } if (m.t === 'run_queue_done') { S.run.state = 'done'; send('mediation'); } render(); };",
      "const arr = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, src.indexOf('];', i) + 2); };",
      "eval(arr('STUDY_KEYS') + arr('PAIR_CACHES') + grab('studySnap') + grab('studyLoad') + grab('studyBusy') + grab('inStudy') + grab('route') + grab('applyStudy'));",
      "applyStudy({ name: 'A', dir: '/a', pairs: [], payloads: {} }); S.run = { state: 'running', done: 0, total: 2 }; S.ldsc = { P01: { rows: ['A'] } };",
      "applyStudy({ name: 'B', dir: '/b', pairs: [], payloads: {} }); const freshB = S.name === 'B' && !S.run && !Object.keys(S.ldsc).length;",
      "drawn = 0; route({ t: 'pair_done', _study: 'A', folder: 'P01' }); route({ t: 'run_queue_done', _study: 'A' }); const quietB = S.name === 'B' && !S.payloads.P01 && sent.length === 0;",
      "route({ t: 'pair_done', _study: 'Z', folder: 'P09' });",
      "const back = applyStudy({ name: 'A', dir: '/a', pairs: [], payloads: { P02: 2 } });",
      "console.log([freshB, quietB, !!S.notice && /\"A\": the run finished/.test(S.notice.msg), back === true && S.run.done === 1 && S.run.state === 'done' && S.payloads.P01 === 1 && S.payloads.P02 === 2 && S.ldsc.P01.rows[0] === 'A', sent.join(','), Object.keys(S.parked).length].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 83: on the screen, a running study is set aside when another is opened (the other starts clean, without the first's results), its news is filed into it without drawing, and opening it again shows it as it is, its requests sent then",
           identical(o, "true true true true mediation 0"))
  }
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "global.S = {}; let sent = 0; global.send = () => { sent++; }; global.STEPS = []; global.liveEndAll = () => 0; global.window = { scrollTo: () => 0 }; global.render = () => 0; global.outUnprepared = () => [];",   # Round 88: go() prepares two-sample outcomes
      "eval(grab('openStudy') + grab('go'));",
      "openStudy('A', 'auto'); openStudy('A', 'auto'); const once = sent === 1;",
      "const ok1 = !!S._opening && S._opening.nav === (S._nav || 0);",
      "go('outliers'); const late = !(S._opening && S._opening.nav === (S._nav || 0));",
      "console.log([once, ok1, late, S.screen].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    st <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
    expect("Round 83: a study opens once however often it is clicked, and an answer that arrives after the analyst moved on (clicks queued while the engine was busy) no longer jumps back to Results",
           identical(o, "true true true outliers") && grepl("if (!w || w.name !== m.name || w.nav !== (S._nav || 0)) break;", st, fixed = TRUE))
  }
})


# Round 84: opening a study is fast. Each script is read once for the settings (it was read again for every setting), a pair's results
# are kept while its files stay the same, their JSON is made once, and the status checks no longer load packages to see if they exist.
local({
  h <- studio_new("R84 fast open"); on.exit(DBI::dbDisconnect(h$con), add = TRUE)
  r <- studio_run(h, "expA", "outB", pair_input("expA", "outB")$input, c("ivw_mre", "egger", "wm"), idx = 1L, params = list(ci_level = 0.95))
  fo <- r$pair_folder; pd <- file.path(h$dir, "03_pairs", fo)
  b0 <- .pair_payload_build(h$dir, fo); p1 <- pair_payload(h$dir, fo)
  t0 <- Sys.time(); p2 <- pair_payload(h$dir, fo); fast <- as.numeric(Sys.time() - t0, units = "secs") < 0.05
  ef <- file.path(pd, "04_tables", "wm_estimates.csv"); keep <- readLines(ef); e <- utils::read.csv(ef); e$estimate[1] <- 123.25; Sys.sleep(1.1); utils::write.csv(e, ef, row.names = FALSE)
  p3 <- pair_payload(h$dir, fo); writeLines(keep, ef); Sys.sleep(1.1); Sys.setFileTime(ef, Sys.time()); p4 <- pair_payload(h$dir, fo)
  wm3 <- Filter(function(x) identical(x$method, "wm"), p3$estimates)[[1]]$estimate
  expect("Round 84: a pair's results are built once and kept while its files stay the same; a changed file builds them again; they equal a fresh build",
         identical(p1, b0) && identical(p2, b0) && fast && identical(wm3, 123.25) && identical(p4, b0))
  expect("Round 84: the results' JSON is made once and is what Shiny would send for them",
         identical(as.character(pair_payload_json(h$dir, fo)), as.character(shiny:::toJSON(b0))) && identical(pair_payload_json(h$dir, fo), pair_payload_json(h$dir, fo)))
  sc <- paste(readLines(file.path(args_dir, "studio_core.R"), warn = FALSE), collapse = "\n"); i <- regexpr("pair_settings <- function", sc, fixed = TRUE)
  ps <- substr(sc, i, i + regexpr("\n}\n", substr(sc, i, nchar(sc)), fixed = TRUE))
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  lh <- paste(readLines(file.path(args_dir, "lhcmr_core.R"), warn = FALSE), collapse = "\n"); ld <- paste(readLines(file.path(args_dir, "ldsc_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 84: settings read each script once; opening a study sends the kept JSON; the LHC-MR and LDSC status checks see whether a package is installed without loading it",
         lengths(regmatches(ps, gregexpr("readLines(", ps, fixed = TRUE))) == 1L && grepl("lns <- readLines(f, warn = FALSE)", ps, fixed = TRUE) &&
         grepl("pair_payload_json(h$dir, p$folder)", app, fixed = TRUE) && grepl("pair_payload_json(h$dir, sn$folder)", app, fixed = TRUE) &&
         grepl('package = nzchar(system.file(package = "lhcMR"))', lh, fixed = TRUE) && grepl('genomicsem = nzchar(system.file(package = "GenomicSEM"))', ld, fixed = TRUE) &&
         identical(pair_settings(h$dir, fo)$ci_level, 0.95))
  # the recipe of a dataset is read without its table (LHC-MR's options read a complete GWAS of millions of rows to know its variant count)
  lib <- library_list(); nm <- if (length(lib)) lib[[1]]$trait else NULL
  srcs <- paste(vapply(c("lhcmr_core.R", "app.R", "mv_core.R", "report_core.R"), function(f) paste(readLines(file.path(args_dir, f), warn = FALSE), collapse = "\n"), ""), collapse = "\n")
  expect("Round 84: a dataset's recipe is read without its table wherever only the recipe is needed (LHC-MR options took 10 to 14 s for a complete GWAS)",
         !is.null(nm) && identical(library_recipe(nm), library_dataset(nm)$recipe) && grepl("is not in the Data library any more", err_of(library_recipe("no such dataset"))) &&
         !grepl("library_dataset\\([^)]*\\)\\$recipe", srcs))
})


# Round 85: CAUSE starts when its measured peak fits in 85% of the free memory (half, the rule for estimated peaks, kept a 16 GB Mac from
# running it on complete GWAS: about 7 GB for 8.5 million variants needed 14 GB free); the other analyses keep the half rule
local({
  mi <- list(available_gb = 9.3, total_gb = 16); cz <- paste(readLines(file.path(args_dir, "cause_core.R"), warn = FALSE), collapse = "\n")
  a <- mem_check(CAUSE_MEMORY_GB(8479693), "CAUSE", mi, share = CAUSE_SAFE_SHARE); b <- mem_check(CAUSE_MEMORY_GB(8479693), "CAUSE", mi); c <- mem_check(CAUSE_MEMORY_GB(1.2e7), "CAUSE", mi, share = CAUSE_SAFE_SHARE)
  expect("Round 85: CAUSE on complete GWAS of 8.5 million variants starts with 9.3 GB free (85% rule), the half rule would refuse it, and a larger one is still refused with the memory it needs",
         isTRUE(a$ok) && !isTRUE(b$ok) && !isTRUE(c$ok) && grepl("starts only when about 12.4 GB is free", c$msg) && grepl('mem_check(CAUSE_MEMORY_GB(rows), "CAUSE", share = CAUSE_SAFE_SHARE)', cz, fixed = TRUE) &&
         identical(RES_SAFE_SHARE, 0.5))
  # Round 86: MRlap measured on the same complete GWAS (7.55 GB peak; the estimate from the package example said 6.9 GB)
  lp <- paste(readLines(file.path(args_dir, "mrlap_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 86: MRlap's memory estimate covers its measured peak on complete GWAS (7.55 GB for 8.5 million variants) and, measured, starts at 85% of the free memory",
         MRLAP_MEMORY_GB(8479693) >= 7.55 && MRLAP_MEMORY_GB(8479693) < 8.5 && isTRUE(mem_check(MRLAP_MEMORY_GB(8479693), "MRlap", list(available_gb = 9.6, total_gb = 16), share = MRLAP_SAFE_SHARE)$ok) &&
         grepl('mem_check(MRLAP_MEMORY_GB(rows), "MRlap", share = MRLAP_SAFE_SHARE)', lp, fixed = TRUE))
})

# Round 86: the draft text of a pair run with MR-Egger and another method that estimates an intercept stopped with an
# error ("the condition has length > 1"); bootstrapped MR-Egger and MR-GRIP have one too. It reports MR-Egger's own intercept
local({
  h <- studio_new("R86 two intercepts"); on.exit(DBI::dbDisconnect(h$con), add = TRUE)
  r <- studio_run(h, "expA", "outB", pair_input("expA", "outB")$input, c("ivw_mre", "egger", "egger_boot"), idx = 1L, params = list(ci_level = 0.95))
  ui_state_set(h$con, "pairs", list(list(ex = "expA", out = "outB", exName = r$exposure_label %||% "expA", outName = r$outcome_label %||% "outB", folder = r$pair_folder, methods = list("ivw_mre", "egger", "egger_boot"))))
  cx <- pair_context(h, ui_state_get(h$con, "pairs")[[1]]); n_ic <- sum(cx$est$estimand == "pleiotropy_intercept")
  eg <- cx$est$estimate[cx$est$estimand == "pleiotropy_intercept" & cx$est$method == "egger"]
  tx <- tryCatch(report_text(list(cx), "R86"), error = function(e) conditionMessage(e))
  expect("Round 86: the draft text of a pair with two intercept estimates (MR-Egger and bootstrapped MR-Egger) is written, with MR-Egger's own intercept",
         n_ic >= 2 && is.list(tx) && any(grepl(sprintf("The MR-Egger intercept was %s", .fnum(eg, 4)), unlist(tx$results), fixed = TRUE)))
})

# Round 87: lead variants that tie at the smallest P value when CAUSE clumps (ieugwasr then picks one arbitrarily) are stated beside the
# CAUSE result and in the Limitations of the draft text (the warning was only in the run's log)
local({
  d <- tempfile(); dir.create(file.path(d, "04_across_pairs", "cause"), recursive = TRUE); h <- list(dir = d); pr <- list(folder = "P01_X__to__Y", exName = "LDL", outName = "CHD")
  cd <- file.path(d, "04_across_pairs", "cause"); utils::write.csv(data.frame(gamma = 0.29, gamma_lower = 0.21, gamma_upper = 0.38, p = 4e-4, instruments = 862), file.path(cd, "P01_X__to__Y_results.csv"), row.names = FALSE)
  w1 <- "3 variants for Se6ATg share the smallest p-value (0), e.g. because of numerical underflow"; w2 <- "2 variants for Ab12Cd share the smallest p-value (0), e.g. because"
  writeLines(c("Refining grid", "Warning messages:", paste("1:", w1), paste("2:", w2)), file.path(cd, "P01_X__to__Y_job.log"))
  a <- cause_ties(h, pr); ta <- cause_ties_text(cause_pair_read(h, pr))
  writeLines(character(), file.path(cd, "P01_X__to__Y_ties.txt")); b <- cause_ties(h, pr)          # a new run with no tie: the empty file wins over the log
  writeLines("4 variants for Zz share the smallest p-value (1e-200), e.g. because p-values were capped", file.path(cd, "P01_X__to__Y_ties.txt")); cc <- cause_ties(h, pr); tc <- cause_ties_text(cause_pair_read(h, pr))
  expect("Round 87: the CAUSE tie warning is read from the run (or, for an older run, its log): regions and variants counted, worded for P = 0 and for capped P values; none when the run had none",
         identical(a$groups, 2L) && identical(a$variants, 5L) && grepl("5 variants in 2 regions shared the smallest P value (recorded as P = 0, a value too small to be stored)", ta, fixed = TRUE) &&
         grepl("CAUSE for LDL on CHD", ta, fixed = TRUE) && is.null(b) && identical(cc$variants, 4L) && grepl("P = 1e-200, perhaps capped", tc, fixed = TRUE) && grepl("of that region was chosen arbitrarily", tc, fixed = TRUE))
  cz <- paste(readLines(file.path(args_dir, "cause_core.R"), warn = FALSE), collapse = "\n"); rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 87: the CAUSE script keeps the clumping's tie warnings in its own file; the draft text's Limitations and the screen message carry the sentence",
         grepl("cl <- withCallingHandlers(ieugwasr::ld_clump(", cz, fixed = TRUE) && grepl('_ties.txt', cz, fixed = TRUE) && grepl("invokeRestart(\\\"muffleWarning\\\")", cz, fixed = TRUE) &&
         grepl("txt$limitations <- c(txt$limitations, unlist(lapply(cau, cause_ties_text)))", rc, fixed = TRUE) && grepl("ties = if (is.null(res)) \"\" else paste(cause_ties_text(res)", app, fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const i = src.indexOf(\"case 'cause':\"); const line = src.slice(i, src.indexOf('\\n', i)); const j = src.indexOf('function causeSection('); const body = src.slice(j, src.indexOf('\\nfunction ', j + 10));",
      "console.log([/ties: m\\.ties \\|\\| ''/.test(line), /R\\.ties \\? '<div role=\"note\"/.test(body)].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 87: on the screen, the CAUSE message keeps the tie sentence and the result shows it as a note", identical(o, "true true"))
  }
})


# Round 88: a two-sample study picks its outcome once. An OpenGWAS outcome download holds only the variants looked up for one exposure
# (so the other exposures' pairs were refused at Check); for each other exposure it is now looked up again at that exposure's instruments,
# with the same proxy settings, preferring a lookup the study already has, then one in the Data library
local({
  rt <- tempfile(); lib <- file.path(rt, "Data library"); dir.create(lib, recursive = TRUE)
  put <- function(nm, src, ex, out) { dir.create(file.path(lib, nm)); writeLines(jsonlite::toJSON(list(trait = nm, source = src, usable_as_exposure = ex, usable_as_outcome = out), auto_unbox = TRUE, null = "null"), file.path(lib, nm, "recipe.json")) }
  py <- list(proxies = TRUE, rsq = 0.8, maf_threshold = 0.3)
  for (x in c("LDL", "BMI", "HDL")) put(x, list(type = "opengwas", id = paste0("ieu-", x)), TRUE, FALSE)
  put("CHD for LDL", list(type = "opengwas", id = "ieu-a-7", trait = "Coronary heart disease", for_exposure = "LDL", params = py), FALSE, TRUE)
  put("CHD for BMI", list(type = "opengwas", id = "ieu-a-7", trait = "Coronary heart disease", for_exposure = "BMI", params = py), FALSE, TRUE)
  put("LDL full", list(type = "opengwas_vcf", id = "ieu-a-300"), TRUE, TRUE); put("My CHD", list(type = "import"), FALSE, TRUE)
  calls <- list(); fake <- function(id, exposure_name, jwt, proxies, rsq, maf, name, offline, root) { calls[[length(calls) + 1L]] <<- list(id = id, ex = exposure_name, proxies = proxies, rsq = rsq, maf = maf); put(name, list(type = "opengwas", id = id, for_exposure = exposure_name, params = py), FALSE, TRUE); list(trait = name) }
  old <- options(mrstudio.say = function(...) NULL); on.exit(options(old), add = TRUE)
  a <- pair_outcome_for("lib:LDL", "lib:CHD for LDL", "tok", root = rt, get_outcome = fake)
  b <- pair_outcome_for("lib:BMI", "lib:CHD for LDL", "tok", root = rt, get_outcome = fake, prefer = c("lib:CHD for LDL", "lib:CHD for BMI"))
  n0 <- length(calls); cc <- pair_outcome_for("lib:HDL", "lib:CHD for LDL", "tok", root = rt, get_outcome = fake)
  again <- pair_outcome_for("lib:HDL", "lib:CHD for BMI", "tok", root = rt, get_outcome = fake)          # the new lookup is reused
  d <- tryCatch(pair_outcome_for("lib:LDL full", "lib:CHD for BMI", "tok", root = rt, get_outcome = fake), error = function(e) conditionMessage(e))
  e <- pair_outcome_for("lib:HDL", "lib:My CHD", "tok", root = rt, get_outcome = fake)
  expect("Round 88: an OpenGWAS outcome serves other exposures: kept for its own, the study's own lookup preferred, otherwise looked up once with the same proxy settings and reused; refused with the reason for an exposure that is not an OpenGWAS instrument set; an imported outcome is left as it is",
         identical(a$out, "lib:CHD for LDL") && identical(b$out, "lib:CHD for BMI") && n0 == 0L && length(calls) == 1L &&
         identical(calls[[1]]$id, "ieu-a-7") && identical(calls[[1]]$ex, "HDL") && isTRUE(calls[[1]]$proxies) && identical(calls[[1]]$rsq, 0.8) && identical(calls[[1]]$maf, 0.3) &&
         grepl("^lib:Coronary heart disease for HDL", cc$out) && grepl("downloaded", cc$note) && identical(again$out, cc$out) &&
         grepl("was downloaded from OpenGWAS only for the instruments of \"BMI\"", d) && grepl("choose the outcome's complete GWAS or import it", d) && identical(e$out, "lib:My CHD"))
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 88: the engine prepares a two-sample study's outcomes on request, the study's own outcomes preferred",
         grepl("out_prepare = {", app, fixed = TRUE) && grepl("prefer = pref", app, fixed = TRUE) && grepl('ev("out_ready", rows = rows)', app, fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const isLib = (id) => String(id).indexOf('lib:') === 0; global.S = { pairs: [{ ex: 'lib:LDL', out: 'lib:CHD' }, { ex: 'lib:BMI', out: 'lib:CHD', outReady: true }, { ex: 'expA', out: 'outB' }, { ex: 'lib:X', out: 'lib:Y', rev: true }, { ex: ['lib:A', 'lib:B'], out: 'lib:Y', mv: true }] };",
      "eval(grab('outUnprepared')); const u = outUnprepared();",
      "const i = src.indexOf(\"case 'out_ready':\"); const blk = src.slice(i, src.indexOf(\"case 'bidir_ready':\", i));",
      "console.log([u.length === 1 && u[0].ex === 'lib:LDL', /kept once/.test(blk), /prepFailOut\\(\\)/.test(blk), /send\\('out_prepare'/.test(src), /out_prepare: 'check'/.test(src)].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 88: on the screen, a two-sample study's OpenGWAS pairs are prepared before the checks (narrated in the Check box); a pair that becomes another one's twin is kept once, and one that cannot be prepared leaves the study with the reason",
           identical(o, "true true true true true"))
  }
})


# Round 89 (the parked Round 77b): on the Data screen every dataset says what it is ready for in that role, and each list offers to import
# a file or find one on OpenGWAS, coming back to the study with the new dataset selected
local({
  pn <- tempfile(); dir.create(pn); file.create(file.path(pn, "EUR.bed"))
  rec <- function(type, ex = TRUE, out = FALSE, facts = list(), ...) list(trait = "T", source = c(list(type = type), list(...)), facts = facts, usable_as_exposure = ex, usable_as_outcome = out)
  L <- function(r, role, panel = pn) .lib_ready(r, role, panel)
  a <- L(rec("opengwas"), "exposure"); b <- L(rec("opengwas", params = list(where = "local", pop = "EUR")), "exposure")
  c1 <- L(rec("import", facts = list(independent = "clump", clump_pop = "EUR")), "exposure"); c2 <- L(rec("import", facts = list(independent = "clump", clump_pop = "EAS")), "exposure")
  d <- L(rec("opengwas_vcf", facts = list(ancestry = "EUR", clump_pop = "EUR")), "exposure"); e <- L(rec("import"), "exposure")
  f <- L(rec("opengwas", for_exposure = "LDL"), "outcome"); g <- L(rec("opengwas_vcf"), "outcome"); h <- L(rec("import", facts = list(scope = "genomewide")), "outcome"); i <- L(rec("import", facts = list(scope = "selected")), "outcome")
  expect("Round 89: each dataset's readiness for its role: OpenGWAS instruments (clumped on the server or here), clumped at Check (with a warning when its 1000 Genomes panel is missing), complete GWAS, an outcome looked up for one exposure, selected variants only",
         identical(a$txt, "Instruments, clumped on OpenGWAS") && grepl("clumped on this computer when downloaded (1000 Genomes EUR)", b$txt, fixed = TRUE) &&
         identical(c1$lvl, "ok") && identical(c2$lvl, "warn") && grepl("needs the 1000 Genomes EAS panel", c2$txt) && grepl("^Complete GWAS: clumped here at Check", d$txt) &&
         identical(e$txt, "Independent instruments, as imported") && grepl("^Looked up for LDL; looked up again for other OpenGWAS exposures", f$txt) && identical(f$lvl, "note") &&
         grepl("^Complete GWAS: any exposure", g$txt) && grepl("^Complete GWAS: any exposure", h$txt) && grepl("^Selected variants only", i$txt) &&
         grepl("ready = .lib_ready(r, role)$txt, ready_lvl = .lib_ready(r, role)$lvl", paste(readLines(file.path(args_dir, "studio_core.R"), warn = FALSE), collapse = "\n"), fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    fjs <- tempfile(fileext = ".js"); on.exit(unlink(fjs), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "global.esc = (x) => String(x); global.svg = () => ''; global.IC = {}; global.S = { out: [] };",
      "eval(grab('rowItem'));",
      "const w = rowItem({ id: 'lib:X', name: 'X', src: 'Imported', n: '1', anc: 'EUR', ready: 'Clumped here at Check: needs the 1000 Genomes EAS panel', ready_lvl: 'warn' }, 'out');",
      "const ok = rowItem({ id: 'lib:Y', name: 'Y', src: 'Imported', n: '1', anc: 'EUR', ready: 'Selected variants only', ready_lvl: 'note' }, 'out');",
      "const i = src.indexOf(\"case 'back_data':\"); const bd = src.slice(i, src.indexOf('break; }', i));",
      "console.log([/needs the 1000 Genomes EAS panel/.test(w) && /var\\(--warn\\)\">Clumped here/.test(w), /var\\(--mute\\)\">Selected variants only/.test(ok),",
      "  /lnk\\('Import a file', 'imp_from_data', key\\) \\+ lnk\\('Find on OpenGWAS', 'find_online', key\\)/.test(src), /S\\.imp\\.fromKey = v; if \\(v === 'out'\\) S\\.imp\\.facts\\.role = 'outcome'; else if \\(v === 'med'\\) S\\.imp\\.facts\\.role = 'either'/.test(src),",
      "  /'lib:' \\+ rec\\.trait/.test(bd) && /S\\.newFromStudy/.test(bd), /You came here from the study/.test(grab('vImport'))].join(' '));"), fjs)
    o <- paste(system2("node", fjs, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 89: on the screen, each dataset row shows its readiness (a missing panel in the warning colour); each list links to Import a file and Find on OpenGWAS; going back to the study selects the dataset just imported or downloaded",
           identical(o, "true true true true true true"))
  }
})


# Round 90: the live box of a study that was set aside (Round 83) stopped typing after two characters ("MR" of "MR-PRESSO is still
# running") once the analyst came back: news filed while it was away started typing there, and the typing timer was not restarted
local({
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const cst = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, i + src.slice(i).search(/;\\r?\\n/) + 1); };",
      "const arr = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, src.indexOf('];', i) + 2); };",
      "global.document = { getElementById: () => null }; let timer = 0; global.setTimeout = () => { timer++; return 1; }; global.clearTimeout = () => 0; global.send = () => 0; global.prevOk = () => false; global.kindMethods = () => []; global.render = () => 0;",
      "global.S = { lives: {}, _pin: {}, _top: {} };",
      "eval(cst('LIVE_GROUP') + cst('LIVE_PREP') + cst('mmss') + arr('STUDY_KEYS') + arr('PAIR_CACHES') + grab('liveOpen') + grab('liveAdd') + grab('liveFollow') + grab('liveTick') + grab('studySnap') + grab('studyLoad') + grab('studyBusy') + grab('inStudy') + grab('applyStudy'));",
      "applyStudy({ name: 'A', dir: '/a', pairs: [], payloads: {} }); S.run = { state: 'running' }; liveOpen('run'); liveAdd('Method 5 of 7: MR-PRESSO', '', 'run_all');",
      "applyStudy({ name: 'B', dir: '/b', pairs: [], payloads: {} });",
      "S._quiet = false; inStudy('A', () => liveAdd('MR-PRESSO is still running (45 s so far)', '', 'run_all'));",
      "const away = S.parked.A.lives.run, awayDone = away.typed === away.lines[away.lines.length - 1].msg.length;",
      "const t0 = timer; applyStudy({ name: 'A', dir: '/a', pairs: [], payloads: {} }); const L = S.lives.run, cur = L.lines[L.lines.length - 1];",
      "console.log([awayDone, L.typed === cur.msg.length, cur.msg === 'MR-PRESSO is still running (45 s so far)', timer > t0].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 90: news that arrives while a study is set aside is shown complete, and coming back to the study restarts the clocks of its live box (the line stopped at \"MR\")",
           identical(o, "true true true true"))
  }
  bg <- paste(readLines(file.path(args_dir, "bg_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 90: a background job sets a UTF-8 locale before it reads the engine (its messages showed \"LDL <e2><86><92> CHD\" when started without one)",
         grepl('!isTRUE(l10n_info()[["UTF-8"]])) for (lc in c("en_US.UTF-8", "C.UTF-8", "UTF-8"))', bg, fixed = TRUE) && regexpr("Sys.setlocale", bg) < regexpr("for (f in files) { source(", bg, fixed = TRUE) &&
         { j <- tempfile(); dir.create(j); writeBin(charToRaw('{"t":"live","msg":"LDL \u2192 CHD"}\n'), file.path(j, "events.jsonl")); writeLines('{"state":"done"}', file.path(j, "status.json")); ln <- bg_poll(j, 0)$lines; identical(Encoding(ln), "UTF-8") && grepl("\u2192", ln, fixed = TRUE) && grepl("\u2192", as.character(shiny:::toJSON(list(j = ln))), fixed = TRUE) })
})


# Round 91: faults the user found on screenshots of a finished run and its instrument review, and the same kinds of faults elsewhere
local({
  # P values: an underflowed P (0) read "0.00e+0" on the screen and would read "P = 0.00e+00" in the draft text and the advice texts
  rules <- utils::read.csv(file.path(args_dir, "feasibility_rules.csv"), stringsAsFactors = FALSE)
  expect("Round 91: an underflowed P value reads \"< 1e-300\" in the draft text and the advice texts (\"P < 1e-300\", never \"P = 0.00e+00\"); P values switch to scientific form below 0.001 as on the screen; no advice calls a built method NOT BUILT",
         identical(.fp(0), "< 1e-300") && identical(.pfix(list(a = "Q (P = < 1e-300)", b = list("P = < 1e-300")))$b[[1]], "P < 1e-300") &&
         identical(.feas_p(0), "< 1e-300") && identical(.feas_p(8.18e-4), "8.18e-04") && identical(.fill("Q P = {p}.", list(p = .feas_p(0))), "Q P < 1e-300.") &&
         !any(grepl("NOT BUILT", unlist(rules))))
  # the narration of a run: the running method's own duration, the replay's own "still running", and how the replay ended
  d <- tempfile(); dir.create(file.path(d, "02_scripts"), recursive = TRUE); dir.create(file.path(d, "04_tables")); dir.create(file.path(d, "03_console"))
  for (f in c("00_run_all.R", "01_ivw_mre.R", "02_presso.R")) file.create(file.path(d, "02_scripts", f))
  said <- character(); old <- options(mrstudio.say = function(msg, why = "") said <<- c(said, msg)); on.exit(options(old), add = TRUE); live_reset()
  live_scripts(d); file.create(file.path(d, "04_tables", "ivw_mre_estimates.csv")); live_scripts(d)
  k <- paste(d, "run_all_console.txt"); st <- .live[[k]]; st$t <- st$t - 400; .live[[k]] <- st; live_scripts(d)
  file.create(file.path(d, "04_tables", "presso_estimates.csv")); live_scripts(d)
  live_scripts(d, "replay_console.txt"); unlink(file.path(d, "04_tables", c("ivw_mre_estimates.csv", "presso_estimates.csv"))); live_scripts(d, "replay_console.txt")
  k2 <- paste(d, "replay_console.txt"); s2 <- .live[[k2]]; s2$t <- s2$t - 60; .live[[k2]] <- s2; live_scripts(d, "replay_console.txt")
  app <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 91: a run's narration gives the running method's own time (\"MR-PRESSO finished in 6 min 40 s\"), marks the replay's \"still running\" as the replay's, and ends with how the replay check ended",
         any(grepl("^MR-PRESSO is still running \\(started at [0-9]{2}:[0-9]{2}:[0-9]{2}\\)$", said)) && any(grepl("^MR-PRESSO finished in 6 min 4[0-9] s \\(2 of 2\\)$", said)) &&
         any(grepl("^Replay check: .* is still running", said)) && grepl("Replay check finished: all %s numbers identical to the first run. Saved.", app, fixed = TRUE))
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 91: figure PDFs are written with Cairo where it is available, which embeds their fonts as journals ask, and the screen reports it",
         grepl('if (.fig_cairo()) grDevices::cairo_pdf(f, width = w, height = h, family = "sans")', rc, fixed = TRUE) && grepl("fonts = .fig_cairo()", app, fixed = TRUE) &&
         grepl('opacity:0;pointer-events:none", "MR GUI engine prototype. Engine: "', app, fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    f <- tempfile(fileext = ".js"); on.exit(unlink(f), add = TRUE)
    writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath(file.path(args_dir, "ui", "studio.js")))),
      "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
      "const cst = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, i + src.slice(i).search(/;\\r?\\n/) + 1); };",
      "const line = (n) => { const i = src.indexOf('const ' + n + ' = '); return src.slice(i, src.indexOf('\\n', i)); };",
      "global.esc = (x) => String(x); global.mmss = (ms) => String(Math.round(ms / 1000)); global.S = { lives: {}, methods: ['ivw_mre', 'presso'], pairs: [{ folder: 'P01' }], payloads: { P01: { variants: Array.from({ length: 300 }, () => ({ keep: true })) } } };",
      "eval([line('fp'), line('pEq'), line('fcol'), cst('DETECTORS')].join('\\n').replace(/^const /gm, 'var ').replace(/;\\s*const /g, '; var '));",
      "var METHODS = [{ id: 'ivw_mre', built: true }, { id: 'presso', built: true }, { id: 'radial', built: true }, { id: 'radegger', built: true }, { id: 'lasso', built: true }]; var minInst = () => null; var pairFact = () => null;",
      "eval(grab('estPairSecs') + grab('etaParts') + grab('detMissing') + grab('vLive'));",
      "S.lives.run = { open: false, t0: 0, typed: 0, lines: [{ msg: 'Method 5 of 7: MR-PRESSO', t: 0 }, { msg: 'MR-PRESSO is still running (45 s so far)', t: 1 }, { msg: 'Replay check finished: all 940 numbers identical to the first run. Saved.', t: 2 }] };",
      "const box = vLive('Running the analysis.', '', false, 'run'); S.run = { state: 'done', times: [870] };",
      "const dm = detMissing(1), est = estPairSecs(['ivw_mre', 'presso'], S.pairs[0]);",
      "console.log([fp(0) === '< 1e-300', pEq(0) === '< 1e-300', pEq(0.012) === '= 0.012', src.indexOf(\"= ' + fp(\") < 0, fcol(1.09e-4) === '0.000109', fcol(-0.032) === '-0.032',",
      "  !/still running/.test(box) && /Finished\\./.test(box) && /all 940 numbers identical/.test(box), etaParts()[1] === 'The run took 14 min 30 s.',",
      "  /Only one detector was run/.test(dm) && /data-a=\"add_detectors\"/.test(dm) && /about 15 min/.test(dm), est > 800 && est < 1000,",
      "  /\\['ivw_mre', 'egger', 'wm', 'mode', 'smode'\\]/.test(src), /'units of ' \\+ esc\\(pair\\.outName\\)/.test(src), !/not built/i.test(src.replace(/=== 'NOT BUILT'/g, ''))].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 91: on the screen, an underflowed P reads \"< 1e-300\", table columns stay in decimals, a finished box drops its \"still running\" notes and says Finished, the time tile says how long the run took, the review offers the missing outlier detectors with a fair time, MR-PRESSO time follows its measurement, the recommended set holds all four detectors, the effect sentence names its units, and no screen says \"not built\"",
           identical(o, paste(rep("true", 13), collapse = " ")))
    if (!identical(o, paste(rep("true", 13), collapse = " "))) cat("   node said:", o, "\n")
  }
})


# Round 78: a trailing comment must never swallow code (Round 75d turned the rest of the Check screen's result handler into a comment, so the
# Check screen never finished in build-34). Every line of the screen code with code before "//" is scanned for statements after it.
local({
  L <- readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8")
  bad <- which(vapply(L, function(l) {
    m <- gregexpr("(?<![:'\"\\\\])//", l, perl = TRUE)[[1]]; if (m[1] < 0) return(FALSE)
    for (k in m) { before <- substr(l, 1, k - 1); q <- function(ch) lengths(regmatches(before, gregexpr(ch, before, fixed = TRUE)))
      if (q("'") %% 2 || q('"') %% 2 || q("`") %% 2) next
      after <- substr(l, k + 2, nchar(l)); return(grepl("\\S", before) && grepl("(break;|render\\(\\);|\\bS\\.[A-Za-z]+ = )", after)) }
    FALSE }, TRUE))
  expect("Round 78: no trailing comment in the screen code swallows statements (the bug that kept the Check screen from finishing)", !length(bad) || { cat("  lines:", bad, "\n"); FALSE })
})

# Round 80: the same check for the engine's R code (a trailing comment turned the UK Biobank overlap flag into a comment during Round 80)
local({
  fs <- c(list.files(args_dir, pattern = "\\.R$", full.names = TRUE), file.path(find_contracts_dir(args_dir), "r", c("script.R", "study.R")))
  fs <- fs[!grepl("test_studio\\.R$", fs)]
  hits <- unlist(lapply(fs, function(f) { L <- readLines(f, warn = FALSE, encoding = "UTF-8")
    k <- vapply(L, function(l) { p <- regexpr("#", l, fixed = TRUE); if (p < 2) return(FALSE); b <- substr(l, 1, p - 1)
      if (!grepl("\\S", b) || lengths(regmatches(b, gregexpr('"', b, fixed = TRUE))) %% 2 || lengths(regmatches(b, gregexpr("'", b, fixed = TRUE))) %% 2) return(FALSE)
      grepl("\\b[a-z_.]+ = (isTRUE|list|c\\(|as\\.|[a-z_.]+\\$|[a-z_.]+\\()|<-\\s*[a-zA-Z(]|\\)\\s*\\}\\s*,?\\s*$|^\\s*\\},", substr(l, p + 1, nchar(l)), perl = TRUE) }, TRUE)
    if (any(k)) paste0(basename(f), ":", which(k)) }))
  expect("Round 80: no trailing comment in the engine's R code swallows code", !length(hits) || { cat("  ", hits, "\n"); FALSE })
})

# Round 76: the run-time estimate counts on each pair only the methods it runs, at its own variants (test 3: 52 min shown, 26.5 min real)
if (nzchar(Sys.which("node"))) local({
  js <- tempfile(fileext = ".js"); on.exit(unlink(js))
  writeLines(c(sprintf("const src = require('fs').readFileSync(%s, 'utf8');", deparse(normalizePath("prototype/engine/ui/studio.js"))),
    "const grab = (n) => { const i = src.indexOf('function ' + n + '('); let d = 0, j = src.indexOf('{', i); for (;; j++) { if (src[j] == '{') d++; if (src[j] == '}') { d--; if (!d) break; } } return src.slice(i, j + 1); };",
    "const uni = 'ivw_mre egger wm divw bwmr grip ivw_steiger mrmix wmode_nome smode_nome pwm egger_boot rucker sign pivw maxlik raps mbe robust_ivw cml conmix lasso egger_simex smode mode smed radegger radial presso'.split(' '), mv = 'mvivw mvegger mvmedian mvlasso mvqhet mvcml'.split(' ');",
    "global.METHODS = uni.map((id) => ({ id })).concat(mv.map((id) => ({ id, mv: true })));",
    "global.S = { facts: [['a','y',78],['a','m1',75],['a','m2',75],['a|m1','y',142],['a|m2','y',117],['a|m1|m2','y',154]].map((f) => ({ ex_id: f[0], out_id: f[1], instruments: f[2], exposures: f[0].split('|') })) };",
    "S.pairs = S.facts.map((f) => ({ ex: f.ex_id.indexOf('|') >= 0 ? f.ex_id.split('|') : f.ex_id, out: f.out_id, mv: f.ex_id.indexOf('|') >= 0 }));",
    "eval(grab('minInst') + grab('pairFact') + grab('estPairSecs'));",
    "console.log(S.pairs.map((p) => estPairSecs(uni.concat(mv), p)).join(' '));"), js)
  e <- as.numeric(strsplit(system2("node", js, stdout = TRUE), " ")[[1]]); real <- c(348, 358, 333, 165, 125, 259)
  expect("Round 76: test 3 pairs estimated within 20% each and the study within 10% of the real 26.5 min (52 min before)",
         length(e) == 6L && all(abs(e / real - 1) < 0.2) && abs(sum(e) / sum(real) - 1) < 0.1)
})


# Round 92: the review build carries one version (0.9.0) in the installer, the screen and every package, so feedback names one build
local({
  shell <- file.path(args_dir, "..", "shell")
  if (file.exists(file.path(shell, "tauri.conf.json"))) {
    tv <- sub('.*"version": *"([^"]+)".*', "\\1", grep('"version"', readLines(file.path(shell, "tauri.conf.json"), warn = FALSE), value = TRUE)[1])
    cv <- sub('^version *= *"([^"]+)".*', "\\1", grep("^version", readLines(file.path(shell, "Cargo.toml"), warn = FALSE), value = TRUE)[1])
    expect("Round 92: the engine version agrees with the installer (tauri.conf.json) and the shell (Cargo.toml)", identical(tv, MRS_VERSION) && identical(cv, MRS_VERSION))
  }
  d <- run_diagnostics(tempfile("r92-"))
  expect("Round 92: Run diagnostics lists the MR Studio version first", identical(d[[1]]$name, "MR Studio") && identical(d[[1]]$detail, paste("version", mrs_version_text())) && startsWith(d[[1]]$detail, paste("version", MRS_VERSION)))   # Round 101: a build adds its number and commit
  htm <- paste(readLines(file.path(pkd, "report.html"), warn = FALSE), collapse = "\n")
  ses <- { z <- file.path(pkd, "reproducibility_capsule.zip"); n <- grep("SESSION.txt$", utils::unzip(z, list = TRUE)$Name, value = TRUE)[1]
           td <- tempfile("r92z-"); utils::unzip(z, n, exdir = td); readLines(file.path(td, n), warn = FALSE) }
  r9a <- paste(unlist(rows19[rows19$item == "9a", ]), collapse = " ")
  expect("Round 92: the report (methods and its note), the capsule's SESSION.txt and STROBE item 9a name the MR Studio version",
         grepl(paste0("run in MR Studio ", MRS_VERSION, ","), htm, fixed = TRUE) && grepl(paste0("Generated by MR Studio ", .mrs_ver_build(), " on"), htm, fixed = TRUE) &&
         any(ses == paste("MR Studio version:", .mrs_ver_build())) && grepl(paste0("MR Studio ", MRS_VERSION, "; R "), r9a, fixed = TRUE))
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 92: the screen keeps the version from the engine's init message and shows it in Settings",
         grepl("S.version = m.version || ''", js, fixed = TRUE) && grepl("' MR Studio ' + esc(S.version)", js, fixed = TRUE))
})


# Round 93: faults found in test H (simulated traits, offline, no token) on build-43
local({
  htm <- paste(readLines(file.path(pkd, "report.html"), warn = FALSE), collapse = "\n")
  ids <- regmatches(htm, gregexpr(' id="[^"]+"', htm))[[1]]
  expect("Round 93: every figure inlined in report.html keeps its own ids and references, so one figure's glyphs and clips never draw another's text",
         length(ids) > 50L && !anyDuplicated(ids) && grepl('href="#f2-', htm, fixed = TRUE) && !grepl('href="#glyph', htm, fixed = TRUE))
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 93: the scatter and funnel plots leave room for the y-axis title beside numbers such as -0.02 and for the x-axis title below",
         lengths(regmatches(rc, gregexpr("mar = c(4.4, 4.9, 0.8, 0.8), mgp = c(3.2, 0.6, 0)", rc, fixed = TRUE))) == 2L)
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 93: a pair of simulated traits gets the one box that says why the complete-GWAS analyses cannot run, not five sections marked Ready",
         grepl("if (st.ok) return", js, fixed = TRUE) && !grepl("st.why === 'simulated') return", js, fixed = TRUE) && grepl("The simulated traits are instrument tables of 40 variants", js, fixed = TRUE))
  expect("Round 93: Settings no longer shows the inactive automatic-updates switch, and the pre-registration switch no longer promises a plan writer",
         !grepl("'Automatic updates'", js, fixed = TRUE) && !grepl("plan writer is planned", js, fixed = TRUE) && grepl("you write and register the plan itself", js, fixed = TRUE))
  expect("Round 93: a finished live box drops its 'You can leave this screen; the work continues' footer",
         grepl("(closed ? '' : foot == null ?", js, fixed = TRUE))
})


# Round 94: faults found in tests A (LDL and BMI on CHD, OpenGWAS) and G (pause, stop) on build-43
local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 94: while paused the screen tells a pair still being calculated (FINISHING, 'Pausing after this pair') from one not yet started (QUEUED), and the bar does not count time for a pair that has not started",
         grepl("r.busy = true;", js, fixed = TRUE) && grepl("r.busy = false; r.i =", js, fixed = TRUE) && grepl("|| r.busy === false) return 0;", js, fixed = TRUE) &&
         grepl("(paused ? 'FINISHING' : 'RUNNING')", js, fixed = TRUE) && grepl("'Pausing after this pair'", js, fixed = TRUE) && !grepl("Paused after pair ' + r.done", js, fixed = TRUE))
  expect("Round 94: a stopped run says 'Analysis stopped' with how many pairs are saved, and marks the pairs it did not run NOT RUN, not QUEUED",
         grepl("stopped ? 'Analysis stopped'", js, fixed = TRUE) && grepl("(stopped ? 'NOT RUN' : 'QUEUED')", js, fixed = TRUE) &&
         grepl("are not kept in the study: to run them later, pick their traits again on the Data step", js, fixed = TRUE) && !grepl("you can run them later from the Data step", js, fixed = TRUE))
  expect("Round 94: the Methods panel shows the sample overlap the analyst stated at Check, not 'Not tested'",
         grepl("ovSummary(fc)", js, fixed = TRUE) && grepl("' (stated)'", js, fixed = TRUE) && !grepl("? 'Not tested' : 'None (simulated)'", js, fixed = TRUE))
  sn <- regmatches(js, regexpr("function steigerNotRun\\(\\) \\{[^}]*", js))
  expect("Round 94: a direction test skipped for a missing prevalence says so, instead of asking to turn on a switch that is on",
         length(sn) == 1L && regexpr("sp && !sp.ok", sn, fixed = TRUE) < regexpr("if (!S.steiger)", sn, fixed = TRUE))
  expect("Round 94: the package panel no longer says the fonts are not embedded when its own check says they are",
         grepl("(S.pkg.fonts ? 'Fonts are embedded in the PDF figures.'", js, fixed = TRUE))
  ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 94: the complete-GWAS box receives the downloads each pair needs, so its 'Download the complete GWAS' button appears",
         grepl("need = r$need %||% list(), linked = r$linked %||% list())", ap, fixed = TRUE))
  sc <- paste(readLines(file.path(args_dir, "studio_core.R"), warn = FALSE), collapse = "\n")
  rec <- list(source = list(type = "opengwas", population = "Mixed"), facts = list(ancestry = "EUR"))
  expect("Round 94: the Data screen badge shows the ancestry OpenGWAS lists (Mixed), not the panel the instruments were clumped with (EUR)",
         grepl("anc = .gwas_pop(r) %||% r$facts$ancestry", sc, fixed = TRUE) && identical(.gwas_pop(rec), "Mixed"))
  expect("Round 94: an automatic outcome name longer than 60 characters is cut after a whole word",
         identical(.word_cut("Coronary heart disease for Body mass index (ieu-a-2) live clump", 60L), "Coronary heart disease for Body mass index (ieu-a-2) live") &&
         identical(.word_cut("CHD for LDL", 60L), "CHD for LDL"))
})

local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 94: pairs reused unchanged do not count as measured times, so they no longer shrink the time left to seconds ('5 sec', then 'soon')",
         grepl("if (!m.reused) { r.times.push(", js, fixed = TRUE))
  expect("Round 94: the Home card of a study running in the background says so, not 'No results yet'",
         grepl("last.status === 'RUNNING' ? 'An analysis is running in the background. Open the study to follow it.'", js, fixed = TRUE))
})

local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 95: the Data panel of a bidirectional study says each pair is analysed in both directions and that the reverse is prepared at Check, and '1 pair' is singular",
         grepl("outcome pair' + (pairs === 1 ? '' : 's')", js, fixed = TRUE) && grepl("analysed in both directions (' + 2 * pairs + ' analyses; the reverse is prepared at Check)", js, fixed = TRUE))
  expect("Round 95: Check names a binary exposure (the reverse of a bidirectional study) and says its estimate is per unit of log odds, the effect of the liability, as a REVIEW item; instrument counts read 'each pair has at least N instruments'",
         grepl("x.scale_x === 'log_or'", js, fixed = TRUE) && grepl("per unit increase in the log odds of the exposure (multiply by 0.693 for a doubling of its odds)", js, fixed = TRUE) && grepl("(each pair has at least ' + minI + ' instruments)", js, fixed = TRUE))
})

# Round 95: test B (bidirectional LDL and CHD on OpenGWAS): MR-cML found no reverse pair and gave nothing, silently
local({
  root <- tempfile("r95-"); dir.create(root)
  mkr <- function(nm, id) { d <- file.path(library_dir(root), nm); dir.create(d, recursive = TRUE); jsonlite::write_json(list(trait = nm, source = list(type = "opengwas", id = id), facts = list(kind = "continuous")), file.path(d, "recipe.json"), auto_unbox = TRUE) }
  mkr("LDL live clump", "ieu-a-300"); mkr("CHD for LDL", "ieu-a-7"); mkr("CHD own instruments", "ieu-a-7"); mkr("LDL for CHD", "ieu-a-300")
  expect("Round 95: a forward pair and its reverse are recognised through their OpenGWAS dataset IDs although the reverse uses other downloads",
         identical(.bidir_trait("lib:CHD for LDL", root), .bidir_trait("lib:CHD own instruments", root)) && identical(.bidir_trait("lib:LDL live clump", root), "ogw:ieu-a-300") &&
         identical(.bidir_trait("expA", root), "expA"))
  h <- list(dir = tempfile("r95s-")); dir.create(h$dir)
  r <- bidir_cml_run(h, list(list(ex = "lib:LDL live clump", out = "lib:CHD for LDL", folder = "P01", exName = "LDL", outName = "CHD", rev = FALSE)), root = root)
  expect("Round 95: a forward pair without a reverse pair is reported as an error, not dropped in silence",
         inherits(r[["P01"]], "bidir_cml_error") && grepl("No reverse pair was found for LDL on CHD", r[["P01"]]$error, fixed = TRUE))
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 95: the report states the bidirectional design, reuses the saved MR-cML result, says the reverse direction was estimated when Steiger did not run, and gives the binary-exposure caveat with its reference",
         grepl("The study was bidirectional: each pair was also analysed in the reverse direction", rc, fixed = TRUE) && grepl("bidir_cml_run(h, pairs, reuse = TRUE)", rc, fixed = TRUE) &&
         grepl("The direction of the effect was examined instead by estimating each pair in both directions", rc, fixed = TRUE) && grepl("[Burgess and Labrecque 2018]", rc, fixed = TRUE) &&
         identical(REFS$burgess2018$pmid, "30039250") && grepl("if (length(.ctx_bin_x(ctxs))) \"burgess2018\"", rc, fixed = TRUE))
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 95: the screen pairs the two directions by OpenGWAS dataset ID too (both-directions panel and MR-cML), and shows an MR-cML error instead of nothing",
         grepl("tkey(q.ex) === tkey(pair.out) && tkey(q.out) === tkey(pair.ex)", js, fixed = TRUE) && grepl("S.bcmlErr = m.errors || {}", js, fixed = TRUE) &&
         grepl("Bidirectional MR-cML: ' + esc(er.error)", js, fixed = TRUE) && grepl("'bcml', 'bcmlErr'", js, fixed = TRUE))
})

# Round 96: test C (multivariable LDL, HDL, TG on CHD from OpenGWAS): the package failed and left empty folders
local({
  d <- tempfile("r96-"); dir.create(file.path(d, "P01", "05_figures"), recursive = TRUE); dir.create(file.path(d, "P01", "04_tables")); writeLines("x", file.path(d, "P01", "04_tables", "t.csv"))
  .rm_empty_dirs(d)
  expect("Round 96: folders an analysis did not use are left out of the package, folders with files stay",
         !dir.exists(file.path(d, "P01", "05_figures")) && file.exists(file.path(d, "P01", "04_tables", "t.csv")))
  sb <- paste(readLines(file.path(args_dir, "strobe_core.R"), warn = FALSE), collapse = "\n"); mv <- paste(readLines(file.path(args_dir, "mv_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 96: the STROBE checklist describes a multivariable OpenGWAS dataset and a dataset without a file, and a missing sample size, instead of stopping the package ('argument is of length zero')",
         grepl("if (!is.null(s$mv)) sprintf(\"OpenGWAS %s, each looked up at the joint instruments\"", sb, fixed = TRUE) && grepl("length(s$file) && !is.na(s$file)", sb, fixed = TRUE) &&
         grepl("length(cx$nx) == 1L && is.finite(cx$nx)", sb, fixed = TRUE))
  expect("Round 96: a multivariable report states the sample overlap answered at Check and the OpenGWAS dataset IDs under Data availability",
         grepl("overlap = tryCatch(overlap_get(h)[[ovk]]", mv, fixed = TRUE) && grepl("overlap_limitations(o)", mv, fixed = TRUE) && grepl("txt$availability <- unique(c(.avail_sentence(src), txt$availability))", mv, fixed = TRUE))
  expect("Round 96: a new study and a new pair get no folders that stay empty (05_report_and_supplement, 06_reproducibility, 05_figures, 06_method_reports, 08_outlier_analysis), and the MVMR forest legend is numbered and says why there is no scatter plot",
         !any(c("05_report_and_supplement", "06_reproducibility") %in% STUDY_FOLDERS) && !any(c("05_figures", "06_method_reports", "08_outlier_analysis") %in% PAIR_SUBFOLDERS) &&
         grepl("\"Figure %d. Direct effects of %s on %s from multivariable MR", mv, fixed = TRUE) && grepl("A scatter plot is not drawn for multivariable MR", mv, fixed = TRUE))
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 96: a pair that ran radial MR without an outlier analysis gets a numbered radial plot, supplementary figures are numbered in order, and the scripts folder carries every pair folder (outlier reruns too), as the capsule does",
         grepl("paste0(tag, \"_radial\")", rc, fixed = TRUE) && grepl("draw_outlier_radial(cx$pay, NULL, est$unit[1])", rc, fixed = TRUE) && grepl("so <- sn   # Round 44", rc, fixed = TRUE) &&
         !grepl("2L * sn", rc, fixed = TRUE) && grepl("list.dirs(file.path(h$dir, \"03_pairs\"), recursive = FALSE)", rc, fixed = TRUE))
})

# Round 97: test D (mediation BMI, TG and HDL on CHD from OpenGWAS)
local({
  mc <- paste(readLines(file.path(args_dir, "mediation_core.R"), warn = FALSE), collapse = "\n")
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  fc <- paste(readLines(file.path(args_dir, "feasibility_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 97: the total effect of a mediation study uses the outcome looked up for the exposure (an outcome download made for another exposure failed the step)",
         grepl("yt <- if (exists(\"pair_outcome_for\")) pair_outcome_for(x, y, jwt, offline, root, get_outcome)", mc, fixed = TRUE) && grepl("y_tot = yt$out", mc, fixed = TRUE) &&
         grepl("total = pair_input(x, med$y_tot %||% y, pthr)", mc, fixed = TRUE) && grepl("p.med.role === 'total' && r.y_tot && r.y_tot !== p.out", js, fixed = TRUE))
  expect("Round 97: Results open on the first pair with results when an earlier one failed, and the mediation section names a step without a result",
         grepl("const k = S.pairs.findIndex((q) => q.folder && S.payloads[q.folder])", js, fixed = TRUE) && grepl("Mediation cannot be computed: ", js, fixed = TRUE))
  expect("Round 97: going Home during Check and back keeps the choices and the Check in progress (preparation counts as busy; the same loaded study is not reloaded)",
         grepl("S.medBusy || S.revBusy || S.mvBusy || S.factsBusy || ld(S.ldsc)", js, fixed = TRUE) && grepl("S._kept = true; return false;", js, fixed = TRUE) &&
         grepl("if (S._kept) { S._kept = false; go(S._stepScreen || 'data'); break; }", js, fixed = TRUE))
  expect("Round 97: the mediation standard set adds weighted mode, MVMR-Egger and MVMR median; mediation steps get no multiple-testing note; a failed pair's reason is shown whole",
         grepl("['ivw_mre', 'egger', 'wm', 'mode', 'smode', 'mvivw', 'mvegger', 'mvmedian']", js, fixed = TRUE) && grepl("n_pairs = if (identical(kind, \"mediation\")) 1L else length(facts)", fc, fixed = TRUE) &&
         grepl("send('feas', { params: S.params, kind: S.kind })", js, fixed = TRUE) && !grepl(".slice(0, 220));", js, fixed = TRUE))
})

# Round 98: test D redo (the built package of a mediation study)
local({
  pk <- build_package(h57, c("scripts", "report", "figs", "tables", "text", "strobe", "refs", "capsule"), 140)
  hs <- paste(readLines(file.path(pk$dir, "report.html"), warn = FALSE), collapse = "\n")
  nums <- as.integer(regmatches(hs, gregexpr("(?<=Figure )[0-9]+(?=\\.)", hs, perl = TRUE))[[1]])
  wb <- paste(readLines(unz(file.path(pk$dir, "supplement_tables.xlsx"), "xl/workbook.xml"), warn = FALSE), collapse = "")
  cap <- utils::unzip(file.path(pk$dir, "reproducibility_capsule.zip"), list = TRUE)$Name
  expect("Round 98: a mediation package numbers its figures 1 to n without restarting, adds the mediation figure, sheet and scripts, and has no multiple-testing line",
         length(nums) > 3 && identical(nums, seq_along(nums)) && file.exists(file.path(pk$dir, "figures", "mediation_forest.pdf")) && grepl("Mediation of the effect of", hs) &&
         grepl("sheet name=\"Mediation\"", wb, fixed = TRUE) && "04_across_pairs/mediation/mediation.R" %in% cap &&
         dir.exists(file.path(pk$dir, "scripts_and_tables", "04_across_pairs")) && !grepl("No correction for multiple testing", hs, fixed = TRUE))
})
local({
  res <- data.frame(mediator = c("(none)", "M", "M", "M", "M"), method = c("univariable", "MVMR", "product", "difference", "two-step"), quantity = c("total", "direct", "indirect", "indirect", "exposure on mediator"),
                    estimate = c(0.4, 0.3, 0.05, 0.1, -0.2), ci_lower = c(0.3, 0.2, 0.02, -0.1, -0.3), ci_upper = c(0.5, 0.4, 0.08, 0.3, -0.1), scale = c("log_or", "log_or", "log_or", "log_or", "beta"), stringsAsFactors = FALSE)
  d <- mediation_forest_df(res)
  expect("Round 98: the mediation figure shows the total, direct and both indirect effects on the outcome scale, not the effect on the mediator",
         identical(d$label, c("Total effect", "M: direct", "M: indirect (product)", "M: indirect (difference)")) && identical(d$unit[1], "log odds ratio"))
  g <- function(id, ex, jwt, px, rsq, maf, name, offline, root) list(trait = name)
  r1 <- ogw_ensure_outcome(list(id = "ieu-a-7", trait = "Coronary heart disease"), "CHD", "Body mass index (ieu-a-2) live clump", list(), NULL, FALSE, root, g)
  r2 <- ogw_ensure_outcome(list(id = "ieu-a-7", trait = "Coronary heart disease"), "CHD", "BMI", list(), NULL, FALSE, root, g)
  expect("Round 98: an outcome looked up for a long exposure name is named with the dataset tags, not cut mid-word",
         identical(r1$name, "Coronary heart disease (ieu-a-7) for ieu-a-2 live clump") && identical(r2$name, "Coronary heart disease for BMI"))
  hf <- studio_new("Round 98 folders"); dir.create(file.path(hf$dir, "archive")); DBI::dbDisconnect(hf$con); hf2 <- studio_open("Round 98 folders")
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect("Round 98: a new study has no empty optional folders, opening a study removes one left empty, and a mediation step cites the dataset it used as outcome",
         !dir.exists(file.path(hf$dir, "02_datasets")) && !dir.exists(file.path(hf$dir, "archive")) && dir.exists(file.path(hf$dir, "03_pairs")) &&
         grepl("srcname(if (identical(p$med$role, \"xm\") && !is.null(p$med$m_out)) p$med$m_out else p$out", rc, fixed = TRUE))
  expect("Round 98: the activity box lists the newest entry first and mediation Results tabs say Pair, not Outcome",
         grepl("r.log.slice().reverse()", js, fixed = TRUE) && grepl("S.pairs.some((q) => !q.mv) ? 'Pair ' : 'Outcome '", js, fixed = TRUE))
  try(DBI::dbDisconnect(hf2$con), silent = TRUE)
})

# Round 99: test E (BMI on type 2 diabetes with an imported genome-wide outcome file)
local({
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  ic <- paste(readLines(file.path(args_dir, "import_core.R"), warn = FALSE), collapse = "\n")
  lc <- paste(readLines(file.path(args_dir, "live_core.R"), warn = FALSE), collapse = "\n")
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  a1 <- .avail_sentence(c("OpenGWAS:ieu-a-2", "Xue_T2D.gz")); a2 <- .avail_sentence(c("OpenGWAS:ieu-a-2", "OpenGWAS:ieu-a-7"))
  expect("Round 99: Data availability names OpenGWAS IDs and imported files separately (an OpenGWAS ID was listed as a file)",
         grepl("dataset IDs: ieu-a-2; and read from the files Xue_T2D.gz [Add the public source", a1, fixed = TRUE) && identical(a2, "Summary statistics were obtained from the public OpenGWAS database (api.opengwas.io), dataset IDs: ieu-a-2; ieu-a-7."))
  sg <- import_suggest(c("SNP", "A1", "A2", "frq_A1", "b", "se", "P"), list(list(SNP = "rs1", A1 = "T", A2 = "C", frq_A1 = "0.06", b = "0.001", se = "0.01", P = "0.9")))
  expect("Round 99: a frq_A1 column is recognised as the effect allele frequency with high confidence",
         identical(sg$eaf$column, "frq_A1") && identical(sg$eaf$confidence, "high"))
  expect("Round 99: an imported outcome states its build, cases and controls and that missing instruments had no proxies; the clumping limitation only with an imported exposure; no multiple-testing line for one pair",
         grepl("genome build GRCh%s", rc, fixed = TRUE) && grepl("were not replaced by proxy variants", rc, fixed = TRUE) && grepl("if (xfile) \"For imported exposure files", rc, fixed = TRUE) &&
         grepl("if (length(ctxs) > 1L) \"No correction for multiple testing", rc, fixed = TRUE))
  expect("Round 99: the import check says it looked at the first rows only, and the recipe records the rows of the whole file",
         grepl("This check looked at the first %s rows only; saving reads the whole file.", ic, fixed = TRUE) && grepl("rows_in_file = n_all", ic, fixed = TRUE))
  expect("Round 99: a long method says when it started instead of a count that stays at 45 s; one pair is not called all 1 pair; the complete-GWAS note says which trait is an OpenGWAS download",
         grepl("is still running (started at %s)", lc, fixed = TRUE) && !grepl("so far)\", if (replay)", lc, fixed = TRUE) && grepl("'All variants have F above 10 (the pair has '", js, fixed = TRUE) &&
         grepl("One trait of this pair is an OpenGWAS download", js, fixed = TRUE))
  expect("Round 99: a file imported through the Data library after leaving a study is selected in the study by its role when the analyst goes back",
         grepl("k = I && (I.fromKey || S.backKey || (rec && (rec.role === 'outcome' ? 'out' : rec.role === 'exposure' ? 'ex' : null)))", js, fixed = TRUE))
})

# Round 100: test F (LDL on CHD with complete OpenGWAS GWAS): MR-APSS in the shipped runtime, memory, and the report of genome-wide methods
local({
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n")
  ac <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  rs <- paste(readLines(file.path(args_dir, "resource_core.R"), warn = FALSE), collapse = "\n")
  ap <- paste(readLines(file.path(args_dir, "mrapss_core.R"), warn = FALSE), collapse = "\n")
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  ry <- file.path(args_dir, "..", "ci", "runtime.yml"); wf <- file.path(args_dir, "..", "..", ".github", "workflows", "prototype.yml")
  expect("Round 100: the runtime installs mvtnorm, which MRAPSS calls without declaring it; the build checks it and the app says so before starting",
         (!file.exists(ry) || any(trimws(readLines(ry, warn = FALSE)) == "- r-mvtnorm")) && (!file.exists(wf) || grepl("\"MRAPSS\", \"mvtnorm\")", paste(readLines(wf, warn = FALSE), collapse = "\n"), fixed = TRUE)) &&
         grepl("MR-APSS needs the package mvtnorm", ac, fixed = TRUE))
  expect("Round 100: MR-APSS memory is the measured peak (about 3.5 GB for 8.5 M variants) with the 85% share, so it starts with 9.9 GB free on a 16 GB computer",
         MRAPSS_MEMORY_GB(8.5e6) > 3 && MRAPSS_MEMORY_GB(8.5e6) < 4 && identical(MRAPSS_SAFE_SHARE, 0.85) &&
         isTRUE(mem_check(MRAPSS_MEMORY_GB(8.5e6), "MR-APSS", list(available_gb = 9.9, total_gb = 16), share = MRAPSS_SAFE_SHARE)$ok) &&
         grepl('mem_check(MRAPSS_MEMORY_GB(rows), "MR-APSS", share = MRAPSS_SAFE_SHARE)', ap, fixed = TRUE))
  expect("Round 100: on Windows free memory is the Available memory (with the standby cache), not FreePhysicalMemory alone",
         grepl("Win32_PerfFormattedData_PerfOS_Memory", rs, fixed = TRUE) && grepl("v[3] / 1073741824", rs, fixed = TRUE))
  expect("Round 100: a complete GWAS downloaded from OpenGWAS gives its population (the ancestry pill and the Methods were empty)",
         identical(.gwas_pop(list(source = list(type = "opengwas_vcf", population = "Mixed"))), "Mixed") && is.null(.gwas_pop(list(source = list(type = "file")))))
  res <- data.frame(quantity = c("observed_effect", "corrected_effect", "difference"), estimate = c(0.143, 0.149, 0.006), se = c(0.018, 0.018, 0.003), p = c(1e-16, 1e-16, 0.023),
                    ci_lower = c(0.109, 0.113, NA), ci_upper = c(0.178, 0.186, NA), instruments = 87L, int_crosstrait = 0.02); attr(res, "labels") <- c(x = "LDL", y = "CHD")
  t <- mrlap_text(res)
  expect("Round 100: MRlap is reported on the standardised scale with the size of the change, and its Methods give the scale, its own instruments and the simulation seed",
         grepl("on the standardised scale", t, fixed = TRUE) && grepl("a change of +4%", t, fixed = TRUE) && !grepl("preferred", t, fixed = TRUE) &&
         grepl("cannot be compared in size", MRLAP_METHODS_TEXT, fixed = TRUE) && grepl("1000 simulated draws", MRLAP_METHODS_TEXT, fixed = TRUE) && grepl("its own instruments", MRLAP_METHODS_TEXT, fixed = TRUE))
  lr <- function(icp) { d <- data.frame(quantity = c("h2", "h2", "rg", "intercept", "intercept", "cross_trait_intercept"), trait = c("LDL", "CHD", NA, "LDL", "CHD", NA),
                                         estimate = c(0.108, 0.066, 0.222, 1.004, icp, 0.021), se = c(0.014, 0.005, 0.054, 0.02, 0.008, 0.006), ci_lower = c(NA, NA, 0.116, NA, NA, NA), ci_upper = c(NA, NA, 0.327, NA, NA, NA),
                                         p = c(NA, NA, 3.8e-5, NA, NA, 4e-4), snps = 1173569); attr(d, "labels") <- c(x = "LDL", y = "CHD"); ldsc_text(d) }
  expect("Round 100: an LDSC intercept clearly below 1 is explained by genomic control; one near 1 is not",
         grepl("intercept below 1 for CHD usually means that genomic control was applied", lr(0.886), fixed = TRUE) && !grepl("genomic control", lr(0.99), fixed = TRUE))
  expect("Round 100: CAUSE states that its estimates depend on the seed, and a tie of P = 0 is worded plainly",
         grepl("depend on a seed recorded in the script", CAUSE_METHODS_TEXT, fixed = TRUE) && !grepl("too small for the computer", paste(deparse(cause_ties_text), collapse = " "), fixed = TRUE))
  e <- list(x = list(run = FALSE, trait = "LDL", why = "The ancestry of \"LDL\" () is not one of the five 1000 Genomes ancestries."))
  expect("Round 100: the allele-frequency sentence of an older ledger gives the OpenGWAS ancestry instead of empty brackets; OpenGWAS files are harmonised to the forward strand",
         grepl("its ancestry (Mixed) is not one", eaf_check_text(e, list(x = "Mixed")), fixed = TRUE) &&
         grepl("both datasets are OpenGWAS files", .strand_forward_text("(action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous)", list(), TRUE), fixed = TRUE) &&
         grepl("the analyst declared", .strand_forward_text("(action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous)", list(), FALSE), fixed = TRUE))
  expect("Round 100: Limitations no longer say clumping was not done when MR Studio clumped the file, name the methods that allow for overlap, and Data availability names OpenGWAS complete GWAS by ID",
         grepl("is.null(cx$src_x$og) && is.null(cx$ledger$clump)", rc, fixed = TRUE) && grepl("reported in the Results, allow%s for sample overlap", rc, fixed = TRUE) &&
         grepl("paste0(\"OpenGWAS:\", s$vcf$id, gwas_cite(s$vcf, \" \"), \" (complete summary statistics)\")", rc, fixed = TRUE) && grepl("of which 1 was absent", rc, fixed = TRUE))
  expect("Round 100: the supplement has a sheet for each genome-wide method that ran",
         grepl("Round 100: genome-wide methods in the supplement", rc, fixed = TRUE) && grepl("\"04_across_pairs\", m, \"*_results.csv\"", rc, fixed = TRUE))
  ok <- tryCatch({ .lhcmr_download(paste0("file://", file.path(tempdir(), "no_such_file.gz")), file.path(tempdir(), "lhc_x.gz"), tries = 1L); "no error" }, error = function(e) conditionMessage(e))
  expect("Round 100: the LHC-MR LD download allows time for a large file and fails in plain words",
         grepl("could not be downloaded from Google Drive after 1 tries", ok, fixed = TRUE) && !grepl("drive.usercontent", ok, fixed = TRUE) && grepl("max(7200, getOption(\"timeout\"))", paste(deparse(.lhcmr_download), collapse = " "), fixed = TRUE))
  expect("Round 100: the screen says it is waiting for memory instead of Ready, gives the start time instead of now, asks prevalence only for a binary trait, and words OpenGWAS strand",
         grepl("Waiting for memory: press Run again when more memory is free.", js, fixed = TRUE) && !grepl("time: 'now'", js, fixed = TRUE) && grepl("time: nowStamp()", js, fixed = TRUE) &&
         grepl("bx = !sb || !!sb.bin_x", js, fixed = TRUE) && grepl("which harmonises its files to the forward strand of the reference genome", js, fixed = TRUE))
})

local({   # Round 101: a full run keeps only the pairs of the current design (test E: a replaced outcome stayed in the package)
  h101 <- studio_new("Round 101 prune"); on.exit(DBI::dbDisconnect(h101$con), add = TRUE)
  ui_state_set(h101$con, "pairs", list(list(ex = "lib:BMI", out = "lib:T2D", folder = "P01_old"), list(ex = "lib:BMI", out = "lib:T2D v2", folder = "P02_new"),
                                       list(ex = list("lib:A", "lib:B"), out = "lib:C", mv = TRUE, folder = "P03_mv"), list(ex = "lib:T2D v2", out = "lib:BMI", rev = TRUE, folder = "P04_rev")))
  prune_pairs(h101$con, list(list(i = 0, ex = "lib:BMI", out = "lib:T2D v2", rev = FALSE), list(i = 1, ex = c("lib:A", "lib:B"), out = "lib:C", rev = FALSE), list(i = 2, ex = "lib:T2D v2", out = "lib:BMI", rev = TRUE)))
  kept <- vapply(ui_state_get(h101$con, "pairs"), `[[`, "", "folder")
  ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 101: a full run drops pairs left out of the design and keeps multivariable and reverse pairs",
         identical(kept, c("P02_new", "P03_mv", "P04_rev")) && grepl("if (identical(length(a$pairs), rq$total)) tryCatch(prune_pairs(rv$h$con, a$pairs)", ap, fixed = TRUE))
})

local({   # Round 101: an LHC-MR download error saved by an earlier version is shown in plain words
  d <- file.path(tempdir(), "lhc101"); dir.create(d, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(list(state = "error", msg = "download from 'https://drive.usercontent.google.com/x' failed"), auto_unbox = TRUE), file.path(d, "fetch_status.json"))
  st <- tryCatch(lhcmr_status(dir = d), error = function(e) NULL)
  expect("Round 101: an old raw LHC-MR download error is replaced by a plain message", !is.null(st) && grepl("^The LD files of the LHC-MR authors could not be downloaded", st$fetch$msg) && !grepl("drive.usercontent", st$fetch$msg, fixed = TRUE))
})

local({   # Round 101: wording found by tests F, E and 6
  w <- function(t, a) sprintf("The ancestry of \"%s\" (%s) is not one of the five 1000 Genomes ancestries.", t, a)
  e2 <- eaf_check_text(list(x = list(run = FALSE, trait = "LDL", why = w("LDL", "Mixed")), y = list(run = FALSE, trait = "CHD", why = "Its ancestry (Mixed) is not one of the five 1000 Genomes ancestries.")))
  expect("Round 101: allele frequencies not checked for an ancestry outside 1000 Genomes read as one sentence, without quotes or a capital after a colon",
         identical(e2, "The effect-allele frequencies of LDL and CHD were not compared with a reference panel, because their ancestry (Mixed) is not one of the five 1000 Genomes ancestries.") &&
         grepl("not compared with a reference panel (no allele frequency in the file).", eaf_check_text(list(x = list(run = FALSE, trait = "A", why = "No allele frequency in the file."))), fixed = TRUE))
  expect("Round 101: lists read A, B and C", identical(.and_join(c("MRlap", "MR-APSS", "CAUSE")), "MRlap, MR-APSS and CAUSE") && identical(.and_join(c("MRlap", "CAUSE")), "MRlap and CAUSE"))
  un <- list(list(trait = "LDL", label = "Mixed"), list(trait = "CHD", label = "Mixed")); acc <- list(at = "2026-11-23T10:00:00")
  rr <- lapply(c("ldsc", "mrlap", "mrapss"), function(m) list(method = m, reference = "the LDSC authors' European LD scores", unknown_ancestry = un, accepted = acc, labels = c("LDL", "CHD")))
  rr[[4]] <- list(method = "cause", reference = "the 1000 Genomes Phase 3 EUR panel (503 people)", unknown_ancestry = un, accepted = acc, labels = c("LDL", "CHD"))
  lim <- ldref_limitations(rr)
  expect("Round 101: the LD panels accepted for one pair on one day make one limitation sentence that names every method and its bias",
         length(lim) == 1L && grepl("^For LDL and CHD, LD score regression, MRlap and MR-APSS used the LDSC authors' European LD scores and CAUSE used the 1000 Genomes Phase 3 EUR panel", lim) &&
         grepl("the analyst chose these panels and accepted on 2026-11-23", lim, fixed = TRUE) && grepl("the corrected effect (MRlap)", lim, fixed = TRUE) && length(ldref_limitations(rr[1])) == 1L)
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n"); js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 101: MR-APSS states the standardised scale in the Methods, the Results and on screen, as MRlap does; Steiger limitation without a comma splice",
         grepl("works on the standardised scale", MRAPSS_METHODS_TEXT, fixed = TRUE) && grepl("on the standardised scale: %s on %s", paste(deparse(mrapss_text), collapse = " "), fixed = TRUE) &&
         lengths(regmatches(js, gregexpr("On the standardised scale (Z divided by the square root of N)", js, fixed = TRUE))) == 2L && grepl("only through the exposure and is sensitive", rc, fixed = TRUE) && !grepl("exposure, it is sensitive", rc, fixed = TRUE))
  expect("Round 101: the import asks for the type of trait instead of presetting continuous; Publish no longer promises to open the folder; Settings names the studies folder",
         grepl("kind: '',", js, fixed = TRUE) && grepl("Choose the type of trait: continuous, or binary (cases and controls).", js, fixed = TRUE) &&
         !grepl("The package folder opens when it is done", js, fixed = TRUE) && grepl("When it is done, press Show the package folder.", js, fixed = TRUE) && !grepl(">Local data library<", js, fixed = TRUE))
  d <- file.path(tempdir(), "bi101"); dir.create(d, showWarnings = FALSE); writeLines('{"build":"build-49","commit":"0123456789abcdef"}', file.path(d, "BUILD_INFO.json"))
  dg <- run_diagnostics(); wf <- paste(readLines(file.path(args_dir, "..", "..", ".github", "workflows", "prototype.yml"), warn = FALSE), collapse = "\n")
  expect("Round 101: the build number and commit come from BUILD_INFO.json, written by the workflow; diagnostics list the genome-wide method packages",
         identical(mrs_build(d), "build-49, commit 0123456") && identical(mrs_build(tempfile()), "") && grepl("BUILD_INFO.json", wf, fixed = TRUE) &&
         any(vapply(dg, function(x) identical(x$name, "Genome-wide method packages") && grepl("MRAPSS", x$detail, fixed = TRUE) && grepl("mvtnorm", x$detail, fixed = TRUE), TRUE)))
})

# ---- Round 102: Test 7 (scenario H, offline demo study) ----
local({
  ov <- overlap_limitations(list(list(simulated = TRUE)))
  expect("Round 102: a simulated pair says its samples do not overlap instead of 'not tested'", grepl("generated independently, so their samples do not overlap", ov, fixed = TRUE) && !grepl("not tested", ov, fixed = TRUE))
  mk <- function(p) list(pay = list(diagnostics = list(list(method = "ivw_mre", name = "cochran_q", value = 63.7, p = p))))
  h1 <- heterogeneity_limitation(list(mk(0.007))); h2 <- heterogeneity_limitation(list(mk(0.007), mk(0.4))); h0 <- heterogeneity_limitation(list(mk(0.4)))
  expect("Round 102: heterogeneity with Q P below 0.05 is written into the limitations, as the Results screen promises",
         startsWith(h1, "Cochran's Q showed more heterogeneity") && startsWith(h2, "In 1 of 2 pairs, Cochran's Q showed") && is.null(h0) && grepl("[Burgess 2023]", h1, fixed = TRUE))
  d <- tempfile("r102-"); dir.create(d); writeLines('{"build": "build-77", "commit": "abcdef0123456789"}', file.path(d, "BUILD_INFO.json"))
  op <- options(mrgui.engine_dir = d); v <- .mrs_ver_build(); options(op)
  expect("Round 102: the package note and the capsule's SESSION.txt name the build", identical(v, paste0(MRS_VERSION, " (build-77, commit abcdef0)")))
})

# ---- Round 103: every univariable method in one run (Part 1 on build-53) ----
local({
  est <- list(egger_boot = list(list(method = "egger_boot", estimand = "causal_effect", p = 0, ci_method = "wald"), list(method = "egger_boot", estimand = "pleiotropy_intercept", p = 0.476, ci_method = "wald")),
              pivw = list(list(method = "pivw", estimand = "causal_effect", p = 0, ci_method = "bootstrap")), divw = list(list(method = "divw", estimand = "causal_effect", p = 0, ci_method = "wald")))
  dg <- list(egger_boot = list(list(method = "egger_boot", name = "egger_boot_nboot", value = 500)))
  mk <- p_resolution_mark(est, dg); mk2 <- p_resolution_mark(list(pivw = list(list(method = "pivw", p = 0, ci_method = "wald"))), list())
  expect("Round 103: a bootstrap P (bootstrap MR-Egger, pIVW with bootstrap Fieller) carries its resolution, 1 / draws; other estimates do not",
         identical(mk$egger_boot[[1]]$p_resolution, 1 / 500) && identical(mk$egger_boot[[2]]$p_resolution, 1 / 500) && identical(mk$pivw[[1]]$p_resolution, 1 / 1000) &&
         is.null(mk$divw[[1]]$p_resolution) && is.null(mk2$pivw[[1]]$p_resolution) && identical(mk$egger_boot[[1]]$p, 0))
  r <- data.frame(estimate = 0.292, ci_lower = 0.194, ci_upper = 0.390, ci_level = 0.95, p = 0, p_resolution = 0.001)
  expect("Round 103: a bootstrap P of 0 reads \"< 0.001\" in the text and the tables, never \"< 1e-300\"; an underflowed P without a resolution keeps \"< 1e-300\"",
         identical(.fp_res(c(0, 0.5, 0), c(0.001, 0.001, NA)), c("< 0.001", "0.500", "< 1e-300")) && identical(.pfix(.ci_phrase(r, "beta")), "beta = 0.292 (95% CI 0.194 to 0.390; P < 0.001)") &&
         identical(.pfix("P = < 0.001"), "P < 0.001") && identical(.fp_res(numeric(0), NULL), character(0)))
  s1 <- software_sentence(c("ivw_mre", "pivw", "raps", "bwmr")); s2 <- software_sentence(c("ivw_mre", "egger"))
  expect("Round 103: the Methods name the package of every method that was run, and MR Studio's own BWMR, instead of 'All estimates come from TwoSampleMR'",
         grepl("TwoSampleMR (commit 5e93ee7d8e)", s1, fixed = TRUE) && grepl("MendelianRandomization (version ", s1, fixed = TRUE) && grepl("mr.raps (commit ", s1, fixed = TRUE) &&
         grepl("and with MR Studio's own implementation of BWMR", s1, fixed = TRUE) && !grepl("MendelianRandomization", s2, fixed = TRUE) && startsWith(s2, "The estimates were computed with TwoSampleMR") &&
         !grepl("All estimates come from", paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n"), fixed = TRUE))
  cx <- list(p = list(exName = "X", outName = "Y"), pay = list(diagnostics = list(list(method = "sign", name = "sign_concordance", value = 0.75, df = 40, p = 1.38e-6),
    list(method = "rucker", name = "rucker_model", value = NA, note = "B"), list(method = "rucker", name = "q_diff", value = 0.2466, df = 1, p = 0.6195))))
  mt <- method_result_text(cx)
  expect("Round 103: the sign concordance test and Rucker's choice are reported in the Results (they were described in the Methods only)",
         length(mt) == 2L && identical(mt[1], "Sign concordance test: 35 of 40 variants (88%) had associations with X and Y in the same direction (binomial test P = 1.38e-06).") &&
         identical(mt[2], "Rucker's model selection chose the random-effects IVW model (model B); the fall in Q from adding the MR-Egger intercept was 0.25 (P = 0.620).") &&
         length(method_result_text(list(p = cx$p, pay = list(diagnostics = list())))) == 0L)
  rc <- paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n"); sc <- paste(readLines(file.path(args_dir, "strobe_core.R"), warn = FALSE), collapse = "\n")
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 103: the report's scatter line, the Results card and the screen's scatter use MR-Egger's own intercept; the workbook has a Diagnostics sheet; STROBE 9a names the packages used",
         grepl('cx$est$estimand == "pleiotropy_intercept" & cx$est$method == "egger"]', rc, fixed = TRUE) && grepl("list(Diagnostics = dgs)", rc, fixed = TRUE) &&
         lengths(regmatches(js, gregexpr("e.method === 'egger' && e.estimand === 'pleiotropy_intercept'", js, fixed = TRUE))) == 3L &&
         !grepl("find((e) => e.estimand === 'pleiotropy_intercept')", js, fixed = TRUE) && grepl("software_pkgs(", sc, fixed = TRUE) && !grepl("TwoSampleMR commit", sc, fixed = TRUE))
  ids <- c("ivw_steiger", "bwmr", "grip", "mrmix", "egger_boot", "pwm", "smode_nome", "wmode_nome", "rucker", "maxlik", "mbe", "cml", "robust_ivw", "pivw", "divw", "raps", "conmix", "lasso", "egger_simex", "wm", "smed", "mode", "smode", "presso", "radial")
  mc <- regmatches(js, regexpr("const MCOL = \\{[^\n]*", js))
  expect("Round 103: the screen's scatter draws a coloured line for every method in its legend (it drew 8 of 27)",
         all(vapply(ids, function(i) grepl(paste0(i, ": '"), mc, fixed = TRUE), TRUE)) && grepl("filter((e) => ['ivw_mre', 'egger', 'radegger'].indexOf(e.method) < 0)", js, fixed = TRUE) && grepl("function methodResults(pay)", js, fixed = TRUE))
  if (nzchar(Sys.which("node"))) {
    gl <- function(n) { L <- strsplit(js, "\n", fixed = TRUE)[[1]]; sub("^const ", "var ", L[startsWith(L, paste0("const ", n, " = "))][1]) }
    f <- tempfile(fileext = ".js"); writeLines(c(gl("fp"), gl("fpr"), "console.log([fpr(0, 0.001) === '< 0.001', fpr(0) === '< 1e-300', fpr(0.0123, 0.001) === '0.012', fpr(0.0005, 0.001) === '< 0.001'].join(' '));"), f)
    o <- paste(system2("node", f, stdout = TRUE, stderr = TRUE), collapse = " ")
    expect("Round 103: on the screen a bootstrap P of 0 reads \"< 0.001\"", identical(o, "true true true true"))
  }
})

# ---- Round 104: every multivariable method in one run (Part 2 on build-54) ----
local({
  sp <- software_pkgs(c("mvivw", "mvegger", "mvqhet"))
  expect("Round 104: MVMR Q-minimisation is named as MR Studio's own implementation, not as MendelianRandomization or boot",
    identical(sp$own, "mvqhet") && identical(sp$ver, paste0("MendelianRandomization (version ", PIN_MENDELIANRANDOMIZATION, ")")))
  expect("Round 104: the software sentence with one package and one own implementation",
    identical(software_sentence(c("mvivw", "mvqhet")), paste0("The estimates were computed with MendelianRandomization (version ", PIN_MENDELIANRANDOMIZATION, ") and with MR Studio's own implementation of MVMR Q-minimisation")))
  m <- paste(t50$methods, collapse = " ")
  expect("Round 104: an MVMR-only report names MR Studio's version, the packages, the R version and the replay check",
    grepl(paste0("run in MR Studio ", .mrs_ver(), ","), m, fixed = TRUE) && grepl("The estimates were computed with MendelianRandomization", m, fixed = TRUE) &&
    grepl(R.version.string, m, fixed = TRUE) && grepl("the replay check was", m, fixed = TRUE) && !grepl("Methods (MendelianRandomization", m, fixed = TRUE))
  expect("Round 104: an MVMR-only report has a Data availability statement",
    length(t50$availability) >= 2L && AVAIL_SCRIPTS %in% t50$availability)
  src <- paste(readLines(file.path(args_dir, "mv_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 104: MVMR-Egger names the exposure it is oriented to; Q-minimisation is credited to MR Studio",
    grepl("oriented to %s [Rees 2017]", src, fixed = TRUE) && grepl("MR Studio's implementation of MVMR::qhet_mvmr", src, fixed = TRUE))
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  lc <- paste(readLines(file.path(args_dir, "live_core.R"), warn = FALSE), collapse = "\n")
  expect("Round 104: the direct effects table puts P on its own line; the replay check names every method when there are 10 or fewer",
    grepl("white-space:nowrap\">P ' + fpr(e.p, e.p_resolution)", js, fixed = TRUE) && grepl("length(ids) <= 10L ||", lc, fixed = TRUE))
})

# ---- Round 104b: no empty folder in a study, its package or its capsule, for every study kind ----
local({
  empt <- function(d) { x <- list.dirs(d, recursive = TRUE, full.names = TRUE); x <- x[!grepl("/[.]engine(/|$)", x)]
    x[vapply(x, function(p) !length(list.files(p, all.files = TRUE, no.. = TRUE)), TRUE)] }
  zempt <- function(z) { if (!file.exists(z)) return(character(0)); n <- utils::unzip(z, list = TRUE)$Name; d <- n[endsWith(n, "/")]
    d[vapply(d, function(x) !any(startsWith(n, x) & n != x), TRUE)] }
  full <- c("scripts", "report", "figs", "tables", "text", "strobe", "refs", "capsule")
  kinds <- list(two_sample = hl, multivariable = h50, mediation = h57)
  if (exists("h63")) { ui_state_set(h63$con, "pairs", lapply(pr63, function(p) c(p, list(methods = list("ivw_mre"))))); kinds$bidirectional <- h63 }
  for (k in names(kinds)) { h <- studio_open(basename(kinds[[k]]$dir)); study_plan_save(h, list(kind = k, aim = "", prereg = FALSE)); pk <- build_package(h, full, 140)   # the app saves the plan when a study is made
    e <- c(empt(h$dir), zempt(file.path(pk$dir, "reproducibility_capsule.zip")))
    expect(sprintf("Round 104b: a %s study, its full package and its capsule have no empty folder%s", gsub("_", "-", k), if (length(e)) paste0(" (", paste(sub(h$dir, "", e, fixed = TRUE), collapse = ", "), ")") else ""), !length(e)) }
})

# ---- Round 105: first launch (clean-install walkthrough as Dr. Chang) ----
local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  sh <- paste(readLines(file.path(args_dir, "..", "shell", "dist", "index.html"), warn = FALSE), collapse = "\n")
  expect("Round 105: the first Home screen says Welcome, not Welcome back, and says MR Studio made the studies folder and that it can be changed in Settings",
    grepl("h1(st.length ? 'Welcome back' : 'Welcome')", js, fixed = TRUE) && grepl("MR Studio saves your studies in a folder it made for you", js, fixed = TRUE) &&
    grepl("You can change it in Settings.", js, fixed = TRUE) && !grepl("Your studies are kept in", js, fixed = TRUE))
  expect("Round 105: a new install starts in the light Paper theme, in the settings, the screen and the loading screen",
    identical(SETTINGS_DEFAULT$theme, "paper") && grepl("theme: 'paper', settings: {}", js, fixed = TRUE) && grepl("THEMES[S.theme] || THEMES.paper", js, fixed = TRUE) &&
    grepl("--bg:#F5F3EE", sh, fixed = TRUE) && isFALSE(SETTINGS_DEFAULT$offline))
})

# ---- Round 106: Dr. Akbarzadeh walkthrough on build-56 ----
local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 106: Settings offers the system's folder picker for the studies folder; the chosen path fills the box and Save checks it as before",
    grepl("data-a=\"folder_pick\"", js, fixed = TRUE) && grepl("Choose folder...", js, fixed = TRUE) && grepl("case 'folder_picked'", js, fixed = TRUE) &&
    grepl("folder_pick = { r <- tryCatch(.with_user_env(pick_folder(", ap, fixed = TRUE) && is.function(pick_folder))
  if (!identical(.Platform$OS.type, "windows") && !identical(Sys.info()[["sysname"]], "Darwin"))
    expect("Round 106: where there is no folder picker the user is told to type the path", grepl("type the path", err_of(pick_folder()), fixed = TRUE))
  d <- tempfile("r106-"); dir.create(file.path(d, "02_scripts"), recursive = TRUE); dir.create(file.path(d, "04_tables")); dir.create(file.path(d, "03_console"))
  for (x in c("01_ivw_mre.R", "02_egger.R", "03_wmedian.R")) writeLines("x", file.path(d, "02_scripts", x))
  said <- character(); op <- options(mrstudio.say = function(msg, why = "") said <<- c(said, msg)); on.exit(options(op), add = TRUE)
  live_scripts(d, "replay_console.txt")
  for (m in c("ivw_mre", "egger")) writeLines("x", file.path(d, "04_tables", paste0(m, "_estimates.csv")))
  live_scripts(d, "replay_console.txt")
  expect("Round 106: the replay check names the quick methods that finished between two looks",
    any(grepl("^Replay check: finished .*\\(2 of 3\\)$", said)))
})

# ---- Round 107: quitting MR Studio during a run (scenario A walkthrough on build-56) ----
local({
  if (!exists("bg_orphan_runs")) source(file.path(args_dir, "bg_core.R"))
  st <- tempfile("r107-"); mk <- function(nm, state, t = "run_pair") { j <- file.path(.bg_jobs_dir(st), nm); dir.create(j, recursive = TRUE)
    saveRDS(list(t = t, a = list(ex = "lib:LDL", out = "lib:CHD"), act = "run_all"), file.path(j, "args.rds"))
    .bg_write_status(j, state = state, t = t, started = "2026-11-25 17:52:09", pid = Sys.getpid()); normalizePath(j, winslash = "/") }
  run <- mk("run_pair_a", "running"); mk("run_pair_b", "done"); mk("build_package_c", "running", "build_package")
  o <- bg_orphan_runs(st)
  expect("Round 107: a pair run still going from before MR Studio was closed is found, with its pair and start time; finished jobs and other kinds are not",
    length(o) == 1L && identical(o[[1]]$pair, "LDL -> CHD") && identical(o[[1]]$started, "17:52"))
  expect("Round 107: a run this window already watches is not reported again", !length(bg_orphan_runs(st, watched = vapply(o, `[[`, "", "job"))))
  ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n"); js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 107: opening the study says what is still finishing and that the rest must be run again; a new run is refused until it ends; its end is announced",
    grepl("bg_watch(o$job, name, quiet = TRUE, on_done = function(r) ev(\"bg_finished\", study = name))", ap, fixed = TRUE) &&
    grepl("is still finishing in the background. Wait until it ends", ap, fixed = TRUE) && grepl("An analysis is still finishing", js, fixed = TRUE) &&
    grepl("case 'bg_finished'", js, fixed = TRUE) && grepl("If you close MR Studio, the pair being calculated still finishes", js, fixed = TRUE))
})

# ---- Round 108: forest plot numbers ran under the scatter plot in a window that was not full size ----
local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 108: the forest plot number column fits a full line with a P value, the card wraps below the scatter plot before it is squeezed, and a very narrow window scrolls instead of hiding numbers",
    grepl("flex:none;width:330px;white-space:nowrap;text-align:right", js, fixed = TRUE) && grepl("margin:6px 342px 0 132px", js, fixed = TRUE) &&
    grepl("<div style=\"overflow-x:auto\"><div style=\"min-width:614px\">' + body + ax", js, fixed = TRUE) &&
    grepl("<section style=\"flex:1 1 620px;min-width:0", js, fixed = TRUE) && !grepl("width:236px", js, fixed = TRUE))
  ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 108: the pairs of a run are kept with the study when the run starts and come back when the study is reopened before they all finished",
    grepl('if (length(a$design)) tryCatch(ui_state_set(rv$h$con, "design", a$design)', ap, fixed = TRUE) && grepl('design = ui_state_get(h$con, "design") %||% list()', ap, fixed = TRUE) &&
    grepl("design: S.pairs.map((p) => ({ ex: p.ex, out: p.out", js, fixed = TRUE) && grepl("(m.design || []).filter((d) => !(m.pairs || []).some((p) => same(p, d)))", js, fixed = TRUE))
  h <- studio_new("Round 108 design")
  ui_state_set(h$con, "design", list(list(ex = "lib:LDL", out = "lib:CHD", rev = FALSE), list(ex = "lib:BMI", out = "lib:CHD", rev = FALSE)))
  expect("Round 108: the outcome download says it can be paired with other OpenGWAS exposures, as Check looks it up again for them",
    grepl("You can still pair this outcome with any other OpenGWAS exposure in a study: before Check, MR Studio looks it up again", js, fixed = TRUE) && !grepl("can only be paired with that exposure", js, fixed = TRUE))
  d <- ui_state_get(h$con, "design")
  expect("Round 108: a saved design reads back with every pair", length(d) == 2L && identical(d[[2]]$ex, "lib:BMI") && identical(d[[2]]$out, "lib:CHD"))
  DBI::dbDisconnect(h$con)
})

# ---- Round 109: further MendelianRandomization estimators (items 20 to 24) ----
h109 <- studio_new("Round 109 methods"); p109 <- pair_input("expA", "outA"); d109 <- p109$input$data[1:14, ]; q109 <- p109$input; q109$data <- d109
mi109 <- MendelianRandomization::mr_input(d109$beta_x, d109$se_x, d109$beta_y, d109$se_y)
r109 <- studio_run(h109, p109$exposure_label, "Outcome 14", q109, c("egger_rob", "hetpen", "clr"), idx = 1L, params = list(ci_level = 0.95)); pl109 <- pair_payload(h109$dir, r109$pair_folder); e109 <- .est_df(pl109)
expect("Round 109: robust / penalised MR-Egger, heterogeneity-penalised and CLR run as scripts and replay identically", identical(pl109$replay$status, "identical") && all(vapply(r109$results, function(x) identical(x$status, "ok"), TRUE)))
expect("Round 109: robust / penalised MR-Egger equals MendelianRandomization::mr_egger(robust, penalized)", abs(e109$estimate[e109$method == "egger_rob"] - MendelianRandomization::mr_egger(mi109, robust = TRUE, penalized = TRUE)@Estimate) < 1e-10)
expect("Round 109: heterogeneity-penalised equals mr_hetpen on the same grid", { r <- d109$beta_y / d109$beta_x; lo <- min(r) - 1; hi <- max(r) + 1
  f <- MendelianRandomization::mr_hetpen(mi109, CIMin = lo, CIMax = hi, CIStep = (hi - lo) / 2000); x <- e109[e109$method == "hetpen", ]
  abs(x$estimate - f@Estimate) < 1e-10 && abs(x$ci_lower - min(f@CILower)) < 1e-10 && abs(x$ci_upper - max(f@CIUpper)) < 1e-10 && is.na(x$se) })
expect("Round 109: CLR gives intervals as diagnostics and no point estimate, as mr_clr", { g <- Filter(function(x) x$name == "clr_interval", pl109$diagnostics)
  r <- d109$beta_y / d109$beta_x; lo <- min(r) - 1; hi <- max(r) + 1; f <- MendelianRandomization::mr_clr(mi109, nx = stats::median(d109$n_x), ny = stats::median(d109$n_y), CIMin = lo, CIMax = hi, CIStep = (hi - lo) / 2000)
  !any(e109$method == "clr") && length(g) == 1L && identical(g[[1]]$note, paste(sprintf("%.4g to %.4g", f@CLRlower, f@CLRupper), collapse = "; ")) })
expect("Round 109: heterogeneity-penalised is refused above 16 variants, CLR without sample sizes", { q <- p109$input; q2 <- q109; q2$data$n_x <- NA
  identical(ADAPTER_HETPEN_V11$check(q, resolve_params(ADAPTER_HETPEN_V11, list()))$codes, "FEAS_TOO_MANY_IVS") && identical(ADAPTER_CLR_V11$check(q2, resolve_params(ADAPTER_CLR_V11, list()))$codes, "FEAS_NEEDS_N") })
expect("Round 109: the package example of mr_clr (LDL on CHD) gives the help page intervals", { f <- MendelianRandomization::mr_clr(MendelianRandomization::mr_input(MendelianRandomization::ldlc, MendelianRandomization::ldlcse, MendelianRandomization::chdlodds, MendelianRandomization::chdloddsse), nx = 17723, ny = 60801)
  abs(f@CLRlower - 2.51) < 0.01 && abs(f@CLRupper - 3.74) < 0.01 })
mv109 <- mv_pair_input(c("expA", "expB"), "outA")
r109m <- studio_run(h109, "A and B", "outA", mv109$input, c("mvivwme", "mvgmm"), idx = 2L, params = list(ci_level = 0.95)); pl109m <- pair_payload(h109$dir, r109m$pair_folder); e109m <- .est_df(pl109m)
expect("Round 109: MVMR-IVW with measurement error and MVMR GMM run as scripts, replay identically and equal the package", { d <- mv109$input$data; d <- d[d$mr_keep, ]
  bx <- as.matrix(d[, c("beta_x_1", "beta_x_2")]); sx <- as.matrix(d[, c("se_x_1", "se_x_2")]); mi <- MendelianRandomization::mr_mvinput(bx = bx, bxse = sx, by = d$beta_y, byse = d$se_y)
  sd <- as.integer(sub("^seed <- ([0-9]+).*", "\\1", grep("^seed <- ", readLines(Sys.glob(file.path(h109$dir, "03_pairs", r109m$pair_folder, "02_scripts", "*mvivwme*.R"))), value = TRUE)[1])); a <- MendelianRandomization::mr_mvivwme(mi, correl.x = diag(2), seed = sd); b <- MendelianRandomization::mr_mvgmm(mi, nx = c(stats::median(d$n_x_1), stats::median(d$n_x_2)), ny = stats::median(d$n_y), cor.x = diag(2))
  identical(pl109m$replay$status, "identical") && all(vapply(r109m$results, function(x) identical(x$status, "ok"), TRUE)) &&
    max(abs(e109m$estimate[e109m$method == "mvivwme"] - a@Estimate)) < 1e-8 && max(abs(e109m$estimate[e109m$method == "mvgmm"] - b@Estimate)) < 1e-8 })
expect("Round 109: the new methods are named in the report, the UI list and the references", { js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  all(c("egger_rob", "hetpen", "clr", "mvivwme", "mvgmm") %in% names(METHOD_REFS)) && all(vapply(c("egger_rob", "hetpen", "clr", "mvivwme", "mvgmm"), function(m) grepl(paste0("id: '", m, "'"), js, fixed = TRUE), TRUE)) &&
    all(c("burgess2018hetpen", "wang2022clr", "patel2023") %in% names(REFS)) })
# ---- Round 109: items 2, 3, 5, 7, 10, 25, 26, 28, 29, 30, 31 ----
local({
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  expect("Round 109: the standard publication set is IVW, MR-Egger, weighted median, weighted mode and simple mode", grepl("['ivw_mre', 'egger', 'wm', 'mode', 'smode'];   // Round 109: the five methods", js, fixed = TRUE) && !grepl("'mode', 'presso', 'radial', 'radegger', 'lasso'];", js, fixed = TRUE))
  expect("Round 109: P cut points print in scientific form (5e-6, not 0.000005 or 5e-06); results such as P < 0.001 are untouched",
    identical(fmt_thr(c(0.000005, 5e-8, 1e-5, 0.05, 0.001)), c("5e-6", "5e-8", "1e-5", "0.05", "0.001")) && identical(sci_thr_text("variants below P 0.000005, P < 5e-06, P = 0.003, P < 0.001"), "variants below P 5e-6, P < 5e-6, P = 0.003, P < 0.001") &&
      grepl("fmtThr(l.threshold)", js, fixed = TRUE) && grepl("fmtThr(r.source.params.pval)", js, fixed = TRUE))
  expect("Round 109: Backspace, Delete, F5, Ctrl or Cmd+R and history keys do nothing outside text fields", grepl("if (blockNavKey(e)) return;", js, fixed = TRUE) && grepl("(k === 'Backspace' || k === 'Delete') && !editableTarget(e.target)", js, fixed = TRUE) && grepl("k === 'F5'", js, fixed = TRUE))
  expect("Round 109: the OpenGWAS window stays open when fewer than 3 instruments pass, and the next download replaces it", grepl("const OGW_FEW = 3;", js, fixed = TRUE) && grepl("replace: (P.replace && P.replace === P.name) ? P.replace : ''", js, fixed = TRUE) &&
    grepl("identical(replace, nm)) unlink(dest, recursive = TRUE)", paste(readLines(file.path(args_dir, "opengwas_core.R"), warn = FALSE), collapse = "\n"), fixed = TRUE))
  expect("Round 109: a binary outcome shows OR (95% CI) with the log OR on Results, in the report sentence and the supplement", grepl("(lor ? 'OR (' + ciL(lev) + ') ", js, fixed = TRUE) && grepl(" log OR' : ", js, fixed = TRUE) &&
    grepl("OR = 1.492 (95% CI 1.105 to 2.014; log OR 0.400; P = ", .ci_phrase(list(estimate = 0.4, ci_lower = 0.1, ci_upper = 0.7, ci_level = 0.95, p = 0.01), "log_or"), fixed = TRUE) &&
    grepl("beta = 0.400 (95% CI", .ci_phrase(list(estimate = 0.4, ci_lower = 0.1, ci_upper = 0.7, ci_level = 0.95, p = 0.01), "beta"), fixed = TRUE))
  expect("Round 109: SIMEX check row on Check and the instrument review table has a column for every method", grepl("Regression dilution of MR-Egger (I", js, fixed = TRUE) && grepl("oMs.map((m) => th(esc(MSHORT(m.id))", js, fixed = TRUE))
})
expect("Round 109: MAF below 0.01 is removed (exposure EAF, else outcome EAF), a variant with no EAF is kept and counted", {
  d <- data.frame(variant_id = paste0("rs", 1:5), eaf_x = c(0.005, 0.5, NA, NA, 0.995), eaf_y = c(0.5, 0.5, 0.004, NA, 0.5), mr_keep = TRUE, harmonise_action = "aligned", stringsAsFactors = FALSE)
  m <- maf_filter(d); identical(m$d$mr_keep, c(FALSE, TRUE, FALSE, TRUE, FALSE)) && m$led$removed == 3L && m$led$not_checked == 1L && identical(m$d$harmonise_action[1], "dropped_low_maf") &&
    grepl("could not be checked (kept)", maf_text(m$led), fixed = TRUE) })
expect("Round 109: the three F statistics follow their formulas, and say why when EAF or N is missing", {
  d <- data.frame(beta_x = c(0.1, 0.2, 0.05), se_x = c(0.01, 0.02, 0.01), eaf_x = c(0.3, 0.5, 0.1), n_x = 10000, stringsAsFactors = FALSE)
  f <- instrument_f_stats(d); r2 <- 2 * d$beta_x^2 * d$eaf_x * (1 - d$eaf_x); f2 <- r2 * (10000 - 2) / (1 - r2); R2 <- sum(r2); fa <- R2 * (10000 - 3 - 1) / (3 * (1 - R2))
  g <- instrument_f_stats(d[, c("beta_x", "se_x")])
  abs(f$f_z$mean - mean((d$beta_x / d$se_x)^2)) < 1e-9 && abs(f$f_r2$mean - mean(f2)) < 1e-9 && abs(f$f_all$value - fa) < 1e-9 && !is.null(g$f_r2$why) && !is.null(g$f_all$why) &&
    grepl("Overall F of the 3 instruments", f_stats_text(f), fixed = TRUE) && grepl("could not be computed", f_stats_text(g), fixed = TRUE) })
expect("Round 109: a library pair reports the MAF filter and the three F statistics in its ledger and at Check", {
  pf <- tryCatch(pair_facts("expA", "outA"), error = function(e) NULL); !is.null(pf) && is.character(pf$f_text) && nzchar(pf$f_text) })
expect("Round 109: SIMEX is added by itself when MR-Egger is chosen and I2GX is below 0.9, and not otherwise", {
  set.seed(1); d <- data.frame(beta_x = rnorm(20, 0.05, 0.01), se_x = 0.02, mr_keep = TRUE); d2 <- d; d2$se_x <- 0.001
  a <- simex_auto(c("ivw_mre", "egger"), d); b <- simex_auto(c("ivw_mre", "egger"), d2); cc <- simex_auto(c("ivw_mre", "wm"), d)
  identical(a$added, "egger_simex") && "egger_simex" %in% a$methods && is.null(b$added) && is.null(cc$added) })
expect("Round 109: each GWAS is cited with author, year and PMID in the report and on the Data card", identical(gwas_cite(list(author = "Xue A", year = 2018, pmid = "30054458")), "(Xue, 2018; PMID 30054458)") &&
  identical(gwas_cite(list(author = "", year = NA, pmid = "")), "") && grepl("pub = sub(", paste(deparse(library_catalogue), collapse = " "), fixed = TRUE))
expect("Round 109: the instrument review reruns every estimating method of the study, not the detectors", { p <- list(estimates = list(list(method = "ivw_mre", estimand = "causal_effect"), list(method = "presso", estimand = "causal_effect"), list(method = "cml", estimand = "causal_effect"), list(method = "mode", estimand = "causal_effect")))
  m <- outlier_methods(p); all(c("ivw_mre", "egger", "wm", "cml", "mode") %in% m) && !("presso" %in% m) &&
    all(c("m_cml", "m_cml_lo", "m_cml_hi") %in% names(outlier_stats(list(), c("ivw_mre", "cml")))) })
expect("Round 109: an outcome lookup says per variant why there was no proxy (or that the panel is needed to tell)", {
  b <- list(info = function(jwt, id = NULL) data.frame(id = "x-1", trait = "T", population = "European", sample_size = 1000, stringsAsFactors = FALSE))
  r <- proxy_reasons(c("rs1", "rs2"), dir = file.path(tempdir(), "no_panel")); is.null(r) && grepl("Why no proxy: ", paste(deparse(ogw_get_outcome), collapse = " "), fixed = TRUE) })
pnl <- file.path(Sys.getenv("HOME"), "Library", "Application Support", "org.R-project.R", "R", "MRStudio.before-clean-test", "ld_ref", "1kg.v3")
if (!file.exists(file.path(pnl, "EUR.bed"))) pnl <- ld_panel_dir()
real_pnl <- file.exists(file.path(pnl, "EUR.bed")) && file.info(file.path(pnl, "EUR.bed"))$size > 1e8   # the full 1000 Genomes panel, not the small test fixture
if (real_pnl && !is.na(plink_bin())) expect("Round 109: real case (VitD1 instruments lost in T2D): no r2 0.8 proxy, or proxies missing in the outcome", {
  r <- proxy_reasons(c("rs9328367", "rs80261862", "rs10469075"), 0.8, 0.3, "EUR", dir = pnl); if (is.null(r) || nrow(r) < 3L) cat("  proxy_reasons gave nothing; plink:", plink_bin(), "panel:", pnl, "\n") else cat("  ", paste(r$why, collapse = " | "), "\n")
  grepl("no other variant has r2 of at least 0.8", r$why[1], fixed = TRUE) && grepl("rs75145722", r$why[2], fixed = TRUE) && grepl("not in this outcome either", r$why[2], fixed = TRUE) && grepl("rs10469076", r$why[3], fixed = TRUE) && grepl("not in this outcome either", r$why[3], fixed = TRUE) })
# ---- Round 109: power (item 1) and the added figures (items 8, 9, 12, 13, 14) ----
expect("Round 109: power reproduces the mRnd worked example (N = 53,218 for 80% power) and rounds sample sizes up", {
  a <- mr_power_2sls(alpha = 0.05, byx = 1.3, bOLS = 1.41, R2xz = 0.01, varx = 1, vary = 116.6, power_target = 0.8); b <- mr_power_2sls(N = 53218, byx = 1.3, bOLS = 1.41, R2xz = 0.01, vary = 116.6)
  c1 <- mr_power_binary(N = 10000, K = 0.5, OR = 1.5, R2xz = 0.02, power_target = 0.8); y <- mr_power_yz(byx = 0.2, R2xz = 0.02, power_target = 0.8)
  identical(a$n_needed, 53218) && abs(b$power - 0.8) < 1e-4 && c1$n_needed == ceiling(c1$n_needed) && y$n_needed == ceiling(y$n_needed) && abs(b$f - (1 + 53218 * 0.01 / 0.99)) < 1e-9 })
expect("Round 109: power inputs are checked (alpha, R2, K, OR, negative residual variance)", grepl("alpha must be", err_of(mr_power_yz(N = 100, alpha = 2, byx = 0.1, R2xz = 0.1))) && grepl("share of cases", err_of(mr_power_binary(N = 100, K = 1.2, OR = 1.5, R2xz = 0.1))) &&
  grepl("not 1", err_of(mr_power_binary(N = 100, K = 0.2, OR = 1, R2xz = 0.1))) && grepl("negative residual variance", err_of(mr_power_2sls(N = 100, byx = 1, bOLS = 2, R2xz = 0.1, vary = 1))))
expect("Round 109: pair power fills N and R2 from the pair and writes the Methods sentence", { pw <- tryCatch(pair_power("expA", "outA", 0.3), error = function(e) e)
  !inherits(pw, "error") && is.finite(pw$power) && grepl("Brion et al. (mRnd", pw$text, fixed = TRUE) && is.finite(pw$r2xz) })
h109f <- studio_new("Round 109 figures"); p109f <- pair_input("expA", "outA")
r109f <- studio_run(h109f, p109f$exposure_label, p109f$outcome_label, p109f$input, c("ivw_mre", "egger", "wm", "mode", "egger_simex"), idx = 1L, params = list(ci_level = 0.95, nboot = 200L, simex_b = 50L))
ui_state_set(h109f$con, "pairs", list(list(ex = "expA", out = "outA", exName = p109f$exposure_label, outName = p109f$outcome_label, folder = r109f$pair_folder, methods = c("ivw_mre", "egger", "wm", "mode", "egger_simex"))))
pk109 <- tryCatch(build_package(h109f, c("report", "figures", "text"), 140), error = function(e) e)
expect("Round 109: the package has the MendelianRandomization forest, funnel, scatter and all-methods plots, the mode density, the MR-Egger plot and the SIMEX plot, named in the report", {
  !inherits(pk109, "error") && { f <- list.files(file.path(pk109$dir, "figures")); keys <- c("mr_forest_variants", "mr_forest_methods", "mr_funnel", "mr_plot", "mr_allmethods", "mode_density", "mregger", "mregger_simex", "knigge_mr_plot", "knigge_forest", "knigge_funnel_ivw", "knigge_funnel_egger")
    all(vapply(keys, function(k) any(grepl(paste0("_", k, "\\.pdf$"), f)), TRUE)) && grepl("as mrrobust mreggersimex", paste(readLines(file.path(pk109$dir, "report.html"), warn = FALSE), collapse = " "), fixed = TRUE) } })
# ---- Round 109: Cook's distance and studentized residuals (item 27) ----
h109o <- studio_new("Round 109 influence"); p109o <- pair_input("expA", "outA"); q109o <- p109o$input; q109o$data$beta_y[3] <- q109o$data$beta_y[3] + 0.5
r109o <- studio_run(h109o, p109o$exposure_label, "Outcome with an outlier", q109o, c("ivw_mre", "egger", "cook", "studres"), idx = 1L, params = list(ci_level = 0.95)); pl109o <- pair_payload(h109o$dir, r109o$pair_folder)
expect("Round 109: Cook's distance (4 / k) and studentized residuals (|r| > 2) flag the variants lm() flags, replay identically and feed the outlier analysis", {
  d <- q109o$data[q109o$data$mr_keep, ]; f <- stats::lm(beta_y ~ beta_x - 1, data = d, weights = 1 / d$se_y^2)
  c1 <- d$variant_id[stats::cooks.distance(f) > 4 / nrow(d)]; s1 <- d$variant_id[abs(stats::rstudent(f)) > 2]
  g <- function(n) { x <- Filter(function(z) identical(z$name, n), pl109o$diagnostics); if (length(x)) .o_ids(x[[1]]$note) else "missing" }
  st <- outlier_sets(pl109o)
  identical(pl109o$replay$status, "identical") && setequal(g("cook_outliers"), c1) && setequal(g("studres_outliers"), s1) && d$variant_id[3] %in% c1 && all(c("cook", "studres") %in% attr(st, "detectors_run")) })
# ---- Round 109: MR-Clust (item 15); runs where mrclust is installed (the build installs it), otherwise only the registration is checked ----
expect("Round 109: MR-Clust is registered (adapter, name, reference, screen) and left out with the reason where it is not installed", {
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n"); rm_ <- runnable_methods(c("ivw_mre", "mrclust"), pair_input("expA", "outA")$input)
  !is.null(param_adapters()$mrclust) && identical(METHOD_REFS$mrclust, "foley2021") && grepl("id: 'mrclust'", js, fixed = TRUE) &&
    (requireNamespace("mrclust", quietly = TRUE) || grepl("not installed", rm_$left_out[["mrclust"]], fixed = TRUE)) })
if (requireNamespace("mrclust", quietly = TRUE)) {
  expect("Round 109: MR-Clust reproduces the published SBP on CAD example (four clusters, one negative, seed 2000030885)", {
    s <- mrclust::SBP_CAD; d <- data.frame(variant_id = s$chr.pos, beta_x = s$bx, se_x = s$bxse, beta_y = s$by, se_y = s$byse, mr_keep = TRUE, stringsAsFactors = FALSE)
    r <- ADAPTER_MRCLUST_V11$run(list(data = d, exposure_id = "SBP", outcome_id = "CAD", outcome_scale = "log_or"), resolve_params(ADAPTER_MRCLUST_V11, list()), list(seed = 2000030885L))
    g <- r$diagnostics[grepl("^mrclust_cluster_", r$diagnostics$name), ]
    nrow(g) == 4L && all(abs(g$value - c(-0.071, 0.053, 0.113, 0.208)) < 5e-4) && identical(as.integer(g$df), c(4L, 45L, 9L, 4L)) && !nrow(r$estimates) &&
      r$diagnostics$value[r$diagnostics$name == "mrclust_null"] == 5 })
  h109c <- studio_new("Round 109 MR-Clust"); p109c <- pair_input("expA", "outA")
  r109c <- studio_run(h109c, p109c$exposure_label, p109c$outcome_label, p109c$input, c("ivw_mre", "mrclust"), idx = 1L, params = list(ci_level = 0.95)); pl109c <- pair_payload(h109c$dir, r109c$pair_folder)
  expect("Round 109: MR-Clust runs as a script in a study and replays identically", identical(pl109c$replay$status, "identical") && all(vapply(r109c$results, function(x) identical(x$status, "ok"), TRUE)))
}
# ---- Round 109: GSMR (item 6); the published example runs where gsmr2 is installed and MRGUI_GSMR_XMAT names gsmr_example.xmat.gz ----
expect("Round 109: GSMR is registered, needs the LD matrix, and is left out with the reason when the LD matrix cannot be built", {
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n"); sc <- paste(readLines(file.path(args_dir, "..", "..", "contracts", "r", "script.R"), warn = FALSE), collapse = "\n")
  h <- studio_new("Round 109 GSMR"); p <- pair_input("expA", "outA"); r <- tryCatch(studio_run(h, p$exposure_label, p$outcome_label, p$input, c("ivw_mre", "gsmr"), idx = 1L), error = function(e) e)
  !is.null(param_adapters()$gsmr) && identical(METHOD_REFS$gsmr, "zhu2018") && grepl("id: 'gsmr'", js, fixed = TRUE) && grepl("spec$ld_file <- sprintf(\"ld_matrix_%s.csv\", a$id)", sc, fixed = TRUE) &&
    !inherits(r, "error") && grepl("^left out,", r$left_out[["gsmr"]]) })
xm <- Sys.getenv("MRGUI_GSMR_XMAT")
if (requireNamespace("gsmr2", quietly = TRUE) && nzchar(xm) && file.exists(xm)) expect("Round 109: GSMR reproduces the gsmr2 vignette (LDL-C on CAD: bxy 0.4234802, SE 0.02221671; 39 weak, 5 in LD, 9 HEIDI outliers)", {
  e <- new.env(); utils::data("gsmr", package = "gsmr2", envir = e); g <- e$gsmr_data
  id <- scan(xm, what = "", nlines = 1, quiet = TRUE); X <- utils::read.table(xm, header = FALSE, skip = 2); sid <- intersect(g$SNP, id); gd <- g[match(sid, g$SNP), ]; o <- match(sid, id)
  ld <- stats::cor(X[, o]); dimnames(ld) <- list(id[o], id[o])
  d <- data.frame(variant_id = as.character(gd$SNP), beta_x = gd$bzx, se_x = gd$bzx_se, p_x = gd$bzx_pval, eaf_x = gd$a1_freq, n_x = gd$bzx_n, beta_y = gd$bzy, se_y = gd$bzy_se, p_y = gd$bzy_pval, n_y = gd$bzy_n, mr_keep = TRUE, stringsAsFactors = FALSE)
  r <- suppressWarnings(ADAPTER_GSMR_V11$run(list(data = d, ld = ld, exposure_id = "LDL-C", outcome_id = "CAD", outcome_scale = "log_or", exposure_scale = "beta"), resolve_params(ADAPTER_GSMR_V11, list(gsmr_std = "yes", gsmr_n_ref = 7703L)), list(seed = 1L)))
  abs(r$estimates$estimate[1] - 0.4234802) < 1e-6 && abs(r$estimates$se[1] - 0.02221671) < 1e-7 })
# ---- Round 109: harmonisation action 1, 2 or 3 as a sensitivity setting (item 4) ----
expect("Round 109: the harmonisation action can be chosen (1 keeps palindromic variants, 3 drops all of them) and each runs as its own labelled pair", {
  a1 <- tryCatch(pair_input("lib:Trait X", "lib:Trait Y", action = 1L), error = function(e) e); a3 <- tryCatch(pair_input("lib:Trait X", "lib:Trait Y", action = 3L), error = function(e) e)
  a0 <- tryCatch(pair_input("lib:Trait X", "lib:Trait Y"), error = function(e) e)
  ok <- !inherits(a1, "error") && !inherits(a3, "error") && !inherits(a0, "error") && identical(a1$ledger$harmonise_action, 1L) && identical(a3$ledger$harmonise_action, 3L) &&
    a1$ledger$dropped_palindromic == 0L && a3$ledger$dropped_palindromic >= a0$ledger$dropped_palindromic && grepl("must be 1, 2 or 3", err_of(pair_input("lib:Trait X", "lib:Trait Y", action = 4L)))
  h <- studio_new("Round 109 harmonisation"); r <- tryCatch(studio_run_harm(h, "lib:Trait X", "lib:Trait Y", 3L, c("ivw_mre"), 0.95, idx = 2L), error = function(e) e)
  ok && !inherits(r, "error") && grepl("[harmonisation: action 3]", r$label, fixed = TRUE) && identical(r$counts$dropped_palindromic, a3$ledger$dropped_palindromic) &&
    grepl("Simulated traits are already aligned", err_of(pair_input("expA", "outA", action = 1L))) })
local({ js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n"); ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  expect("Round 109: Results offers the other harmonisation actions, the comparison and the report name each one", grepl("function harmSection(pair)", js, fixed = TRUE) && grepl("case 'run_harm'", js, fixed = TRUE) &&
    grepl("x.kind === m.rec.kind", js, fixed = TRUE) && grepl("BG_ACTIONS$run_harm", ap, fixed = TRUE) && grepl("harmonised again with TwoSampleMR action", paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = " "), fixed = TRUE)) })
# ---- Round 109: GRAPPLE and PPMR (items 18, 19); full simulations are in private_notes/r109_methods/<pkg>/tests.R ----
expect("Round 109: GRAPPLE, MVMR GRAPPLE and PPMR are registered and left out with the reason where their package is not installed", {
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n"); rm_ <- runnable_methods(c("ivw_mre", "grapple"), pair_input("expA", "outA")$input)
  all(c("grapple", "mvgrapple", "ppmr") %in% names(param_adapters())) && identical(METHOD_REFS$ppmr, "yuan2020") && grepl("id: 'mvgrapple'", js, fixed = TRUE) &&
    (requireNamespace("GRAPPLE", quietly = TRUE) || grepl("not installed", rm_$left_out[["grapple"]], fixed = TRUE)) })
if (requireNamespace("GRAPPLE", quietly = TRUE)) {
  h109g <- studio_new("Round 109 GRAPPLE"); p109g <- pair_input("expA", "outA")
  r109g <- studio_run(h109g, p109g$exposure_label, p109g$outcome_label, p109g$input, c("ivw_mre", "grapple"), idx = 1L, params = list(ci_level = 0.95, gr_loss = "l2")); pl109g <- pair_payload(h109g$dir, r109g$pair_folder); e109g <- .est_df(pl109g)
  expect("Round 109: GRAPPLE runs as a script, replays identically and agrees with IVW within one SE on the demo pair", {
    g <- e109g[e109g$method == "grapple" & e109g$estimand == "causal_effect", ]; i <- e109g[e109g$method == "ivw_mre" & e109g$estimand == "causal_effect", ]
    identical(pl109g$replay$status, "identical") && nrow(g) == 1L && abs(g$estimate - i$estimate) < g$se })
}
if (requireNamespace("PPMR", quietly = TRUE)) expect("Round 109: PPMR gives an estimate with likelihood-ratio P values and no interval, as the package gives none", {
  set.seed(7); J <- 20; R <- 0.3^abs(outer(1:J, 1:J, "-")); bx <- as.vector(R %*% stats::rnorm(J, 0, 0.05)); by <- 0.3 * bx + as.vector(R %*% stats::rnorm(J, 0, 0.002))
  d <- data.frame(variant_id = paste0("rs", 1:J), beta_x = bx, se_x = 0.01, beta_y = by, se_y = 0.01, n_x = 10000, n_y = 10000, mr_keep = TRUE)
  ld <- R; dimnames(ld) <- list(d$variant_id, d$variant_id)
  r <- ADAPTER_PPMR_V11$run(list(data = d, ld = ld, exposure_id = "X", outcome_id = "Y", outcome_scale = "beta"), resolve_params(ADAPTER_PPMR_V11, list()), list(seed = 1L))
  nrow(r$estimates) == 1L && is.finite(r$estimates$estimate) && is.na(r$estimates$ci_lower) && is.finite(r$estimates$p) })
# ---- Round 109: mode estimates table (item 9) ----
expect("Round 109: the mode table gives simple, weighted and NOME mode estimates at phi 0.25, 0.5 and 1, equal to mr_mbe, and the workbook has the sheet", {
  p <- pair_input("expA", "outA"); d <- p$input$data[p$input$data$mr_keep, ]; k <- data.frame(id = d$variant_id, bx = d$beta_x, by = d$beta_y, sx = d$se_x, sy = d$se_y)
  m <- mode_table(k, 0.95, iterations = 200L); f <- MendelianRandomization::mr_mbe(MendelianRandomization::mr_input(k$bx, k$sx, k$by, k$sy), weighting = "weighted", stderror = "delta", phi = 0.5, seed = 12345L, iterations = 200L)
  w <- m[m$estimator == "Weighted mode" & m$phi == 0.5, ]
  nrow(m) == 12L && setequal(unique(m$estimator), c("Simple mode", "Weighted mode", "Simple mode (NOME)", "Weighted mode (NOME)")) && abs(w$estimate - f@Estimate) < 1e-12 && abs(w$se - f@StdError) < 1e-12 &&
    grepl("sheets$`Mode estimates` <- mtab", paste(readLines(file.path(args_dir, "report_core.R"), warn = FALSE), collapse = "\n"), fixed = TRUE) })
expect("Round 109: the MBE method with simple-mode weighting runs (mr_mbe calls it unweighted; Simple mode returned nothing before)", {
  p <- pair_input("expA", "outA"); r <- ADAPTER_MBE_V11$run(p$input, resolve_params(ADAPTER_MBE_V11, list(mbe_weighting = "simple", mbe_iter = 500L)), list(seed = 1L)); nrow(r$estimates) == 1L && is.finite(r$estimates$estimate) })
# ---- Round 109: PubMed ID at import (item 10) ----
expect("Round 109: an imported file can carry the PubMed ID of its GWAS; it is checked, saved and cited in the report and on the Data card", {
  r <- import_save(fxx, "ex.csv", ",", names(exd), mpl, modifyList(fct("Trait P", "exposure", "yes"), list(pmid = "PMID 30054458", author = "Xue A", year = 2018)))
  rec <- library_recipe("Trait P"); bad <- err_of(import_save(fxx, "ex.csv", ",", names(exd), mpl, modifyList(fct("Trait Q", "exposure", "yes"), list(pmid = "30054x58"))))
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n")
  identical(rec$facts$pmid, "30054458") && grepl("digits only", bad, fixed = TRUE) && grepl("PMID 30054458", gwas_cite_rec(rec), fixed = TRUE) && grepl("id=\"f_pmid\"", js, fixed = TRUE) })
# ---- Round 109: the further figures on the Results screen (items 8, 9, 11 to 14) ----
expect("Round 109: Results shows the package's further figures of a pair, drawn by the same code (SVG with caption for every figure)", {
  js <- paste(readLines(file.path(args_dir, "ui", "studio.js"), warn = FALSE), collapse = "\n"); ap <- paste(readLines(file.path(args_dir, "app.R"), warn = FALSE), collapse = "\n")
  pr <- ui_state_get(h109f$con, "pairs")[[1]]; cx <- pair_context(h109f, pr); td <- tempfile("figs-"); dir.create(file.path(td, "figures"), recursive = TRUE)
  xf <- extra_figures(cx, td, "P", 160); keys <- vapply(xf, `[[`, "", "key"); sv <- vapply(xf, function(x) { s <- grep("[.]svg$", x$files, value = TRUE); if (length(s)) paste(readLines(s[1], warn = FALSE), collapse = "") else "" }, "")
  all(c("mr_forest_variants", "mr_funnel", "mr_plot", "mode_density", "mregger", "mregger_simex", "knigge_mr_plot", "knigge_forest", "knigge_funnel_ivw", "knigge_funnel_egger") %in% keys) && all(grepl("<svg", sv, fixed = TRUE)) &&
    grepl("function figsSection(pair)", js, fixed = TRUE) && grepl("case 'pair_figs'", js, fixed = TRUE) && grepl("pair_figs = {", ap, fixed = TRUE) })
# ---- Round 109: the figures script of each pair reproduces the shipped figures ----
expect("Round 109: each pair gets 05_figures/make_figures.R (in the capsule), and running it redraws the pair's figures identical to the package", {
  pk <- build_package(h109f, c("figs", "capsule"), 140); pr <- ui_state_get(h109f$con, "pairs")[[1]]; fd <- file.path(h109f$dir, "03_pairs", pr$folder, "05_figures")
  cap <- utils::unzip(file.path(pk$dir, "reproducibility_capsule.zip"), list = TRUE)$Name
  td <- tempfile("mkfig-"); dir.create(td); file.copy(list.files(fd, full.names = TRUE), td)
  rb <- default_rscript(); old <- Sys.getenv(c("RHOME", "R_LIBS"), unset = NA); Sys.setenv(RHOME = R.home(), R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))   # as the app starts R (a moved runtime needs RHOME)
  wd <- getwd(); setwd(td); o <- suppressWarnings(system2(rb, c("--vanilla", "make_figures.R"), stdout = TRUE, stderr = TRUE)); setwd(wd)
  for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
  shipped <- list.files(file.path(pk$dir, "figures"), pattern = "^P1_.*[.]png$"); redrawn <- list.files(file.path(td, "figures"), pattern = "^P1_.*[.]png$")
  same <- vapply(shipped, function(f) file.exists(file.path(td, "figures", f)) && identical(unname(tools::md5sum(file.path(pk$dir, "figures", f))), unname(tools::md5sum(file.path(td, "figures", f)))), TRUE)
  if (!all(same) && nzchar(Sys.getenv("MRGUI_FIG_EVIDENCE"))) { ev <- Sys.getenv("MRGUI_FIG_EVIDENCE"); dir.create(ev, showWarnings = FALSE, recursive = TRUE); for (f in shipped[!same]) { file.copy(file.path(pk$dir, "figures", f), file.path(ev, paste0("package_", f)), overwrite = TRUE); file.copy(file.path(td, "figures", f), file.path(ev, paste0("script_", f)), overwrite = TRUE) }; file.copy(file.path(td, "figure_inputs.rds"), ev, overwrite = TRUE); writeLines(c(Sys.getlocale(), paste(names(options()[c("scipen", "digits", "OutDec")]), unlist(options()[c("scipen", "digits", "OutDec")])), paste("devices:", paste(names(grDevices::dev.list()), collapse = " ")), utils::capture.output(print(utils::sessionInfo()))), file.path(ev, "package_session.txt")) }   # evidence of a mismatch, kept only when MRGUI_FIG_EVIDENCE names a folder
  if (!all(same)) cat("  figures script output (last lines):", utils::tail(o, 8), sep = "\n    ")
  if (!all(same)) for (f in shipped[!same]) { if (!file.exists(file.path(td, "figures", f))) { cat("  not redrawn:", f, "\n"); next }; a <- png::readPNG(file.path(pk$dir, "figures", f)); b <- png::readPNG(file.path(td, "figures", f)); cat("  not identical:", f, if (identical(dim(a), dim(b))) sprintf("share of differing pixels %.5f, max diff %.3f", mean(abs(a - b) > 1e-6), max(abs(a - b))) else "size differs", "\n") }
  ck <- c(script = file.exists(file.path(fd, "make_figures.R")), capsule = any(grepl("05_figures/make_figures.R$", cap)), count = length(shipped) >= 10L, same_set = setequal(shipped, redrawn), identical = all(same))
  if (!all(ck)) cat("  failed checks:", names(ck)[!ck], "| shipped only:", setdiff(shipped, redrawn), "| redrawn only:", setdiff(redrawn, shipped), "| status:", attr(o, "status") %||% 0L, "\n")
  all(ck) })
# ---- Round 109a: R's cached text height is cleared before every figure device (labels were sometimes drawn about 11 px off) ----
expect("Round 109a: the cached height of M is cleared before each figure device and the report PDF, without leaving a device open", {
  fw <- paste(deparse(.fig_write), collapse = "\n"); pw <- paste(deparse(pdf_write), collapse = "\n")
  n0 <- length(grDevices::dev.list()); cur <- grDevices::dev.cur(); .fig_metric_reset(); same_dev <- identical(grDevices::dev.cur(), cur) && length(grDevices::dev.list()) == n0
  td <- tempfile("mreset-"); dir.create(td); est <- data.frame(label = c("IVW", "MR-Egger", "Weighted median"), estimate = c(0.3, 0.2, 0.25), ci_lower = c(0.2, 0.05, 0.15), ci_upper = c(0.4, 0.35, 0.35), col = c("#2b3fc7", "#c0392b", "#16a085"), ratio = FALSE, unit = "beta", stringsAsFactors = FALSE)
  f1 <- .fig_write(file.path(td, "a"), 140, 50, function() .draw_forest(est)); f2 <- .fig_write(file.path(td, "b"), 140, 50, function() .draw_forest(est))
  p1 <- grep("[.]png$", f1, value = TRUE); p2 <- grep("[.]png$", f2, value = TRUE)
  lengths(regmatches(fw, gregexpr(".fig_metric_reset()", fw, fixed = TRUE))) == 3L && grepl(".fig_metric_reset()", pw, fixed = TRUE) && same_dev &&
    length(p1) == 1L && identical(unname(tools::md5sum(p1)), unname(tools::md5sum(p2))) })
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail)); if (n_fail) quit(status = 1)
