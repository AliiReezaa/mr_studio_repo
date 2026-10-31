# Conformance and contract tests. Run from the contracts/ folder:
#   Rscript r/tests/test_contracts.R <output_dir>
args <- commandArgs(trailingOnly = TRUE)
out_dir <- tempfile("run-", tmpdir = if (length(args)) args[1] else tempdir())   # always fresh
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(mrgui.schema = "sql/schema_v1.sql")
for (f in c("r/canonical.R", "r/adapter.R", "r/runner.R", "r/security.R", "r/adapters/reference_adapters.R")) source(f)

n_pass <- 0L; n_fail <- 0L
expect <- function(desc, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("PASS ", desc, "\n") }
  else { n_fail <<- n_fail + 1L; cat("FAIL ", desc, "\n") }
}
near <- function(a, b, tol = 1e-12) isTRUE(all(abs(a - b) <= tol * pmax(1, abs(b))))

# ---- reference implementation from the pinned TwoSampleMR commit ------------
ref <- new.env()
exprs <- parse("ref/TwoSampleMR_mr_R_5e93ee7.R", keep.source = FALSE)
for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                     as.character(e[[2]]) %in% c("mr_ivw_mre", "mr_egger_regression")) eval(e, ref)
expect("reference functions loaded from pinned TwoSampleMR source",
       all(c("mr_ivw_mre", "mr_egger_regression") %in% ls(ref)))

# ---- simulated harmonised IV table (C5.1) ------------------------------------
set.seed(20261029)
k <- 30
bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE)
sx <- runif(k, 0.005, 0.012); sy <- runif(k, 0.008, 0.02)
by <- 0.3 * bx + rnorm(k, 0.002, 0.01) + rnorm(k, 0, sy)
dat <- data.frame(variant_id = sprintf("1:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G",
                  beta_x = bx, se_x = sx, eaf_x = 0.3, beta_y = by, se_y = sy, eaf_y = 0.3,
                  p_x = 2 * pnorm(-abs(bx / sx)), p_y = 2 * pnorm(-abs(by / sy)), n_x = 1e5, n_y = 1e5,
                  harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
dat$mr_keep[c(4, 17)] <- FALSE          # dropped palindromes must be ignored
input <- list(kind = "pair_iv", data = dat, exposure_id = "src:opengwas:ieu-a-300",
              outcome_id = "src:opengwas:ieu-a-7", outcome_scale = "log_or")
keep <- dat[dat$mr_keep, ]

ctx0 <- list(seed = 1L)
r_ivw <- ADAPTER_IVW_MRE$run(input, resolve_params(ADAPTER_IVW_MRE), ctx0)
t_ivw <- ref$mr_ivw_mre(keep$beta_x, keep$beta_y, keep$se_x, keep$se_y)
expect("IVW-MRE estimate matches TwoSampleMR", near(r_ivw$estimates$estimate, t_ivw$b))
expect("IVW-MRE SE matches TwoSampleMR", near(r_ivw$estimates$se, t_ivw$se))
expect("IVW-MRE P matches TwoSampleMR", near(r_ivw$estimates$p, t_ivw$pval))
expect("IVW-MRE Cochran Q matches TwoSampleMR", near(r_ivw$diagnostics$value, t_ivw$Q))
expect("IVW-MRE uses only mr_keep rows", r_ivw$estimates$n_ivs_used == 28L)

r_eg <- ADAPTER_EGGER$run(input, resolve_params(ADAPTER_EGGER), ctx0)
t_eg <- ref$mr_egger_regression(keep$beta_x, keep$beta_y, keep$se_x, keep$se_y, NULL)
e1 <- r_eg$estimates[r_eg$estimates$estimand == "causal_effect", ]
e2 <- r_eg$estimates[r_eg$estimates$estimand == "pleiotropy_intercept", ]
expect("Egger slope, SE, P match TwoSampleMR", near(e1$estimate, t_eg$b) && near(e1$se, t_eg$se) && near(e1$p, t_eg$pval))
expect("Egger intercept, SE, P match TwoSampleMR", near(e2$estimate, t_eg$b_i) && near(e2$se, t_eg$se_i) && near(e2$p, t_eg$pval_i))
expect("Egger Rucker Q matches TwoSampleMR", near(r_eg$diagnostics$value, t_eg$Q))

# C5.2: swapping the effect allele of some variants must not change any estimate
flip <- input; idx <- c(1, 5, 9, 22)
flip$data$beta_x[idx] <- -flip$data$beta_x[idx]; flip$data$beta_y[idx] <- -flip$data$beta_y[idx]
expect("IVW invariant to allele coding", near(ADAPTER_IVW_MRE$run(flip, resolve_params(ADAPTER_IVW_MRE), ctx0)$estimates$estimate, r_ivw$estimates$estimate, 1e-10))
expect("Egger invariant to allele coding (orientation inside adapter)",
       near(ADAPTER_EGGER$run(flip, resolve_params(ADAPTER_EGGER), ctx0)$estimates$estimate, r_eg$estimates$estimate, 1e-10))

# ---- feasibility verdicts (C6.2, C11) ----------------------------------------
two <- input; two$data <- keep[1:2, ]
expect("IVW runs with 2 IVs", ADAPTER_IVW_MRE$check(two, list())$status != "block")
expect("Egger blocked with 2 IVs (FEAS_TOO_FEW_IVS)", identical(ADAPTER_EGGER$check(two, list())$codes, "FEAS_TOO_FEW_IVS"))
weak <- input; weak$data$se_x <- abs(weak$data$beta_x) / 2      # F = 4
expect("weak instruments flagged DATA_WEAK_IVS", "DATA_WEAK_IVS" %in% ADAPTER_IVW_MRE$check(weak, list())$codes)
expect("undeclared parameter refused", inherits(try(resolve_params(ADAPTER_IVW_MRE, list(robust = TRUE)), silent = TRUE), "try-error"))
expect("out-of-range parameter refused", inherits(try(resolve_params(ADAPTER_IVW_MRE, list(ci_level = 1.2)), silent = TRUE), "try-error"))

# ---- canonical JSON, fingerprints, seeds (C8.3–C8.5) ------------------------
expect("canonical JSON: sorted keys, 17 digits, escapes",
       identical(.utf8_bytes(canonical_json(list(b = 1L, a = I(c(0.1, 2)), c = "\u00e9\n\""))),
                 .utf8_bytes("{\"a\":[0.10000000000000001,2],\"b\":1,\"c\":\"\u00e9\\u000a\\\"\"}")))
expect("canonical JSON refuses Inf", inherits(try(canonical_json(Inf), silent = TRUE), "try-error"))
expect("canonical JSON refuses named atomic vectors", inherits(try(canonical_json(c(a = 1)), silent = TRUE), "try-error"))
expect("key order does not change fingerprint",
       identical(task_fingerprint("method", list(x = "a", y = "b"), list(p = 1, q = 2), ADAPTER_IVW_MRE),
                 task_fingerprint("method", list(y = "b", x = "a"), list(q = 2, p = 1), ADAPTER_IVW_MRE)))
s1 <- task_seed("method:pair:A→B:ivw_mre")
expect("seed is deterministic and 31-bit", identical(s1, task_seed("method:pair:A→B:ivw_mre")) && s1 >= 0L)
expect("different tasks get different seeds", s1 != task_seed("method:pair:A→B:egger"))

# ---- write protocol, reuse and resume (C8.7, C8.8) --------------------------
root <- file.path(out_dir, "Test.mrgui")
con <- db_open(root)
blobs_in <- list(ivs = sha256_text(canonical_json(.df_rows(dat))))
tid <- "method:pair:src:opengwas:ieu-a-300→src:opengwas:ieu-a-7:ivw_mre"
r1 <- run_method_task(con, root, ADAPTER_IVW_MRE, input, list(), tid, blobs_in)
expect("first run executes and commits", !r1$reused && db_task(con, tid)$state == "done")
expect("stored blob passes integrity check", nchar(blob_get_text(root, r1$blob)) > 100)
r2 <- run_method_task(con, root, ADAPTER_IVW_MRE, input, list(), tid, blobs_in)
expect("same fingerprint is reused, not re-run", r2$reused)
r3 <- run_method_task(con, root, ADAPTER_IVW_MRE, input, list(ci_level = 0.9), tid, blobs_in)
expect("changed parameter changes fingerprint and re-runs", !r3$reused && r3$fingerprint != r1$fingerprint)
DBI::dbExecute(con, "UPDATE tasks SET state='running' WHERE task_id=?", params = list(tid))
db_recover(con)
expect("crash recovery resets running to pending", db_task(con, tid)$state == "pending")
orphan <- blob_put_text(root, "{\"orphan\":true}")      # crash after blob, before commit
expect("uncommitted blob is not referenced by any result",
       nrow(DBI::dbGetQuery(con, "SELECT 1 FROM results WHERE result_blob=?", params = list(orphan))) == 0)
eg_tid <- sub("ivw_mre$", "egger", tid)
r4 <- run_method_task(con, root, ADAPTER_EGGER, two, list(), eg_tid, list(ivs = "x"))
expect("blocked task stored with status blocked", r4$result$status == "blocked" && db_task(con, eg_tid)$state == "blocked")
expect("journal mode is WAL", DBI::dbGetQuery(con, "PRAGMA journal_mode")[[1]] == "wal")
DBI::dbDisconnect(con)

# Same task twice in fresh projects gives byte-identical estimates (determinism)
con2 <- db_open(file.path(out_dir, "Test2.mrgui"))
r5 <- run_method_task(con2, file.path(out_dir, "Test2.mrgui"), ADAPTER_IVW_MRE, input, list(), tid, blobs_in)
DBI::dbDisconnect(con2)
expect("re-run in a new project gives identical estimates",
       identical(canonical_json(.df_rows(r5$result$estimates)), canonical_json(.df_rows(r1$result$estimates))))

# ---- security helpers (C12) ---------------------------------------------------
expect("formula injection neutralised", identical(csv_safe_text(c("=HYPERLINK(1)", "LDL", "-1+1", NA)), c("'=HYPERLINK(1)", "LDL", "'-1+1", NA)))
expect("HTML escaped", identical(html_escape("<script>a&b</script>"), "&lt;script&gt;a&amp;b&lt;/script&gt;"))
expect("zip-slip refused", inherits(try(safe_path_join(out_dir, "../evil.txt"), silent = TRUE), "try-error"))
expect("absolute path refused", inherits(try(safe_path_join(out_dir, "/etc/passwd"), silent = TRUE), "try-error"))
expect("normal relative path allowed", grepl("panels/EUR.bed$", safe_path_join(out_dir, "panels/EUR.bed")))
tok <- "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ0ZXN0dXNlciJ9.c2lnbmF0dXJlX3BhcnQ"
expect("token redacted from logs", !grepl(tok, redact(paste("auth", tok), character()), fixed = TRUE))

# ---- files for JSON Schema validation ----------------------------------------
writeLines(result_to_json(r1$result), file.path(out_dir, "result_ivw.json"))
writeLines(result_to_json(r4$result), file.path(out_dir, "result_blocked.json"))
eg_ok <- run_method_task(db_open(file.path(out_dir, "Test3.mrgui")), file.path(out_dir, "Test3.mrgui"),
                         ADAPTER_EGGER, input, list(), eg_tid, blobs_in)
writeLines(result_to_json(eg_ok$result), file.path(out_dir, "result_egger.json"))
writeLines(manifest_json(ADAPTER_IVW_MRE), file.path(out_dir, "manifest_ivw_mre.json"))
writeLines(manifest_json(ADAPTER_EGGER), file.path(out_dir, "manifest_egger.json"))

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
