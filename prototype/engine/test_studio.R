# Headless tests for the screen logic and the Shiny server (no browser, no desktop shell).
# Run:  Rscript prototype/engine/test_studio.R        (from the repository root)
args_dir <- normalizePath(file.path("prototype", "engine"), winslash = "/")
source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir)); source(file.path(args_dir, "import_core.R")); source(file.path(args_dir, "opengwas_core.R")); source(file.path(args_dir, "ld_core.R")); source(file.path(args_dir, "strobe_core.R")); source(file.path(args_dir, "feasibility_core.R")); source(file.path(args_dir, "report_core.R")); source(file.path(args_dir, "pdf_core.R")); source(file.path(args_dir, "outlier_core.R")); source(file.path(args_dir, "mv_core.R")); source(file.path(args_dir, "mediation_core.R")); source(file.path(args_dir, "bidir_core.R")); source(file.path(args_dir, "ldsc_core.R")); source(file.path(args_dir, "lhcmr_core.R")); source(file.path(args_dir, "app.R"))
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
expect("settings default: aurora theme, setup not done", identical(st$theme, "aurora") && identical(st$setup_done, FALSE))
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
mpx <- mp; mpx$se <- NULL; vx <- import_validate(fi, hd$sep, hd$columns, mpx, fl)
expect("import: a missing required column blocks the import and names it", !vx$ok && any(grepl("Standard error", vx$errors)))
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
expect("import: z-scores are refused for now", !import_validate(fi, hd$sep, hd$columns, mp, list(frq = "eaf", scope = "selected", scale = "z", n = ""))$ok)
fx_ <- tempfile(fileext = ".txt"); utils::write.table(data.frame(FID = 1, IID = 2, PHENO = 1, rs1 = 0, rs2 = 1), fx_, quote = FALSE, row.names = FALSE)
expect("import: individual-level data is refused before anything is read", grepl("individual-level", err_of(import_read_head(fx_))))
expect("import: VCF and zip files are refused with a clear message", grepl("not supported yet", err_of(import_read_head(fi, "x.vcf.gz"))))
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
shiny::testServer(mr_server(args_dir), {
  session$setInputs(act = list(t = "init", nonce = 1))
  session$setInputs(act = list(t = "study_create", name = "", aim = "", kind = "basic", nonce = 2))
  expect("server: a blank study name does not crash and creates nothing", is.null(rv$h))
  session$setInputs(act = list(t = "study_create", name = "Server study", aim = "a", kind = "basic", prereg = FALSE, nonce = 3))
  expect("server: creating a study sets the current study", !is.null(rv$h) && grepl("Server study$", rv$h$dir))
  session$setInputs(act = list(t = "run_pair", i = 0, ex = "expA", out = "outB", methods = list("ivw_mre", "egger"), ci = 0.9, nonce = 4))
  pp <- ui_state_get(rv$h$con, "pairs")
  expect("server: a pair run is recorded in the study database with its labels and folder", length(pp) == 1L && pp[[1]]$ex == "expA" && grepl("^P01_", pp[[1]]$folder))
  session$setInputs(act = list(t = "run_pair", i = 1, ex = "expB", out = "outB", methods = list("ivw_mre"), ci = 0.95, nonce = 5))
  pf <- studio_pair_folders(rv$h$dir)
  expect("server: a second pair gets its own folder", length(pf) == 2L && grepl("^P02_", pf[2]))
  session$setInputs(act = list(t = "run_pair", i = 0, ex = "expA", out = "outB", methods = list("ivw_mre", "egger"), ci = 0.9, nonce = 6))
  expect("server: rerunning a pair with the same settings reuses its folder", length(studio_pair_folders(rv$h$dir)) == 2L && length(ui_state_get(rv$h$con, "pairs")) == 2L)
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
expect("overlap text: simulated traits keep the old sentence", length(lim0) == 1L && grepl("was not tested", lim0) && !grepl("analyst stated", lim0))
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
expect("opengwas report: no wrong claim that the analyst declared independence, the limitations say the panel is a reference sample, data availability names the public source", !grepl("did not perform clumping", htm12, fixed = TRUE) && !grepl("rests on the analyst", htm12, fixed = TRUE) && grepl("reference sample", htm12, fixed = TRUE) && grepl("public OpenGWAS database (api.opengwas.io), dataset IDs: sim-exp-1; sim-out-1", htm12, fixed = TRUE) && !grepl("Add the public source", htm12, fixed = TRUE))
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
expect("feasibility: the demo pair (40 strong instruments) leaves every method available, IVW without any rule", all(vapply(fz$methods, `[[`, "", "state") %in% c("ok", "check")) && identical(fz$methods$ivw_mre$state, "ok") && !length(fz$methods$ivw_mre$reasons))
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
  any(grepl("Q-statistic minimisation", tx$methods)) && any(grepl("bootstrap intervals \\(60 iterations\\)", tx$methods)) })
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
# ---- Round 62: BWMR ----
r62 <- studio_run(h58, p58$exposure_label, p58$outcome_label, p58$input, "bwmr", idx = 5L, params = list(ci_level = 0.95)); pl62 <- pair_payload(h58$dir, r62$pair_folder); e62 <- .est_df(pl62)
expect("Round 62: BWMR runs as a script with the authors' R 4 fix, replays identically and covers the simulated effect 0.30", identical(pl62$replay$status, "identical") && nrow(e62) == 1L && e62$ci_lower <= 0.3 && 0.3 <= e62$ci_upper &&
  any(grepl("inherits(Sigma_hat", readLines(Sys.glob(file.path(h58$dir, "03_pairs", r62$pair_folder, "02_scripts", "*_bwmr.R"))), fixed = TRUE)))
# ---- Round 63: bidirectional cML ----
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
  expect("Round 65: the LHC-MR script parses, pins lhcMR and GenomicSEM, sets the seed and clumps locally", !inherits(tryCatch(parse(text = sc), error = function(e) e), "error") && any(grepl(PIN_LHCMR, sc)) && any(grepl(PIN_GENOMICSEM, sc)) && any(grepl("set.seed(seed)", sc, fixed = TRUE)) && any(grepl("assignInNamespace(\"clump_data\"", sc, fixed = TRUE)))
  res <- data.frame(parameter = c("axy", "ayx", "tX", "tY", "h2X", "h2Y", "iXY"), estimate = c(-0.128, -0.128, 0.032, -0.124, 0.256, 0.546, -0.103), se = 0.02, p = 1e-5, ci_lower = -0.2, ci_upper = 0, variants = 1e6, variants_used = 1e5, bimodal = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), stringsAsFactors = FALSE); attr(res, "labels") <- c(x = "BMI", y = "Height")
  expect("Round 65: the LHC-MR text gives both directions, the confounder effects and the bimodality caution", { t <- lhcmr_text(res); grepl("BMI on Height -0.128", t) && grepl("confounder effects 0.032 on BMI", t) && grepl("bimodal", t) })
  expect("Round 65: with nothing installed the LHC-MR status says so", { st <- lhcmr_ref_status(file.path(tempdir(), "nolhc")); !st$ld && !st$rho })
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

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail)); if (n_fail) quit(status = 1)
