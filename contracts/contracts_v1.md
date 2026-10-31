# MR GUI — Engineering contracts (v1)

Contract version **1.0.0** · 2026-10-29 · status: draft for freeze at the end of the week-1 prototype. Companion to `requirements_spec.md` v1.5, `mr_inventory.json` v3.1 and `project_plan.md`.

These contracts are the fixed interfaces that let the data layer, feasibility engine, method adapters, task runner, UI and report builder be built in parallel. Anything not written here is an implementation detail and may change freely. Anything written here changes only through the versioning rule in C13.

Keywords MUST, MUST NOT, SHOULD and MAY are used in their usual sense. Machine-readable versions live in `contracts/` (JSON Schema, SQLite DDL, R reference code).

---

## C0. Principles

1. **Summary data only.** No object in these contracts can hold per-person rows.
2. **Store estimates on their native additive scale** (beta or log odds / log hazard). Odds ratios, per-SD conversions and wording are produced only at the presentation layer, so every number in the database is directly comparable and re-usable.
3. **Content addressing.** Every stored data file is named by the SHA-256 of its bytes. Nothing is overwritten; a pointer in SQLite says which file is current.
4. **Determinism.** A task with the same fingerprint MUST produce the same stored result on the same operating system, and results within the tolerance in C8.6 across operating systems.
5. **Single writer.** Only the main engine process writes to the project database. Workers compute and return.
6. **Adapters are pure.** A method adapter reads its input, writes only to its scratch folder, makes no network calls and leaves no global state.

## C1. Identifiers

| Object | ID form | Rule |
|---|---|---|
| Project | ULID (26 chars) | Created once; never reused. |
| Dataset | `src:<source>:<native_id>` e.g. `src:opengwas:ieu-a-7`, `src:finngen:R12_I9_CHD`, `src:file:<sha256-12>` | Stable across machines for public sources, so two users get the same ID for the same study. |
| Trait | EFO CURIE when mapped (`EFO:0001645`), else `trait:<slug>` | |
| Exposure–outcome pair | `pair:<exposure_dataset_id>→<outcome_dataset_id>` | Direction matters. |
| Method | `id` field from `mr_inventory.json` (e.g. `ivw_mre`, `egger`, `mr_presso`) | Adapter IDs MUST match an inventory ID. |
| Task | `<task_type>:<subject>:<method or step>` e.g. `method:pair:A→B:ivw_mre` | Human-readable, stable identity. Not a hash. |
| Fingerprint | SHA-256 hex (64) | Content identity of a task's inputs (C8.3). |
| Blob | SHA-256 hex of file bytes | File name in the blob store. |

## C2. Canonical variant table (summary statistics)

Every dataset, whatever its source, is converted on import to this table and stored as Parquet (zstd). Column names are fixed.

| Column | Type | Required | Meaning and rules |
|---|---|---|---|
| `variant_id` | string | yes | `chr:pos:ref:alt` on GRCh37, upper-case alleles. Primary key within a dataset. |
| `rsid` | string | no | dbSNP ID if known, `null` otherwise. |
| `chr` | string | yes | `1`–`22`, `X`. Y and MT are dropped on import in v1. |
| `pos_b37` | int32 | yes | 1-based position on GRCh37. |
| `pos_b38` | int32 | no | Kept for a later switch of canonical build. |
| `ea` | string | yes | Effect allele: the allele whose count the effect refers to. |
| `oa` | string | yes | Other allele. |
| `eaf` | float64 | no | Frequency of `ea` in the GWAS sample, in (0, 1). |
| `beta` | float64 | yes | Effect per copy of `ea`, on the scale in the dataset card (C3). Odds ratios are converted to log odds on import. |
| `se` | float64 | yes | Standard error of `beta`, > 0. |
| `p` | float64 | yes | Two-sided P value as reported. If only `beta`/`se` exist, computed from the normal distribution and flagged. |
| `log10p` | float64 | no | `-log10(p)` when reported, to keep precision for very small P values. |
| `n` | float64 | no | Per-variant sample size if reported. |
| `info` | float64 | no | Imputation quality if reported. |
| `qc_flags` | int32 | yes | Bit field (C2.1). 0 means clean. |

**C2.1 QC flag bits.** 1 = palindromic (A/T or C/G); 2 = EAF far from reference panel (|ΔEAF| above the threshold in the dataset card); 4 = P recomputed; 8 = lifted over from GRCh38; 16 = rsID mapped from position; 32 = multi-allelic site split; 64 = low INFO; 128 = duplicate removed. Flags never delete a row; deletion decisions are recorded in instrument sets (C4).

**C2.2 Import rules.** Rows with missing `beta`, `se`, alleles or position are dropped and counted in the import report. Alleles are normalised to upper case; indels are kept with their full allele strings. If the declared build is absent, the build is inferred and the inference recorded. Files that look like individual-level data (one row per person, genotype calls, more than a configured number of columns with sample-like names) are rejected (S-06).

## C3. Dataset card

One JSON document per dataset (`schemas/dataset_card.schema.json`). It travels with the variant table and is what the feasibility engine and the units passport read.

Required fields: `dataset_id`, `source` (`opengwas | gwas_catalog | finngen | panukb | ukb_ppp | eqtlgen | file`), `source_version`, `retrieved_at` (RFC 3339), `variant_table_blob` (SHA-256), `n_variants`, `trait_label`, `trait_id`, `trait_type` (`continuous | binary | count | molecular`), `effect_scale` (`beta | log_or | log_hr`), `units` (object, C3.1), `sample_size` (`n_total`, `n_cases`, `n_controls`, `n_effective`), `build_declared`, `build_used` (always `GRCh37` in v1), `ancestry` (C3.2), `cohorts` (list of cohort or biobank names, used for overlap), `full_sumstats` (bool: true if genome-wide, false if top hits only), `data_use_terms` (string or URL), `checksums` (original download SHA-256), `import_report` (counts of dropped and flagged rows).

**C3.1 Units passport.** `scale` (`sd | natural | log_transformed | inverse_normal | log_odds | unknown`), `unit_label` (e.g. "mg/dL"), `sd_estimated` (implied trait SD from C2 columns, with method), `sd_source` (`reported | estimated | unknown`), `user_confirmed` (bool).

**C3.2 Ancestry.** `declared` (list of `{group, proportion}` from metadata), `inferred` (list from allele frequencies, with `method` and `reference`), `classification` (`single | multi | unknown`), `main_group` (`EUR | EAS | AFR | AMR | SAS | MID | CSA | null`), `agreement` (`agree | disagree | inferred_only | declared_only`). The 90% single-ancestry cut-off is a parameter recorded here, not a constant in code.

## C4. Instrument set

The result of instrument selection for one exposure (or one set of exposures for multivariable methods), stored as a JSON card plus a Parquet variant list.

Card fields: `instrument_set_id` (fingerprint), `exposure_dataset_ids`, `selection` (`p_threshold`, `clump_r2`, `clump_kb`, `ld_panel` with ancestry and version, `region` for cis designs as `chr:start-end` on GRCh37, `method` = `clump | region | none`), `n_candidates`, `n_selected`, `dropped` (list of `{variant_id, reason}` with reason codes C11), `f_stat_mean`, `f_stat_min`, `r2_total`.

Per-variant columns: `variant_id`, `beta_x`, `se_x`, `eaf_x`, `f_stat` (approximate F = beta_x² / se_x²), `r2` (when EAF and units allow), `proxy_for` (variant it replaces, else `null`), `proxy_r2`.

## C5. Harmonised analysis inputs

Adapters never see raw datasets. They receive one of five input kinds, chosen by the method's `instrument_source` and `needs_ld` fields in the inventory.

| Input kind | Used by `instrument_source` | Contents |
|---|---|---|
| `pair_iv` | `clumped_independent`, `relaxed_threshold_independent` | Harmonised IV table for one exposure and one outcome. |
| `pair_iv_ld` | `correlated_variants…`, `allele_score_or_variants` | As `pair_iv` plus an LD correlation matrix for those variants, aligned to the same effect alleles. |
| `multi_iv` | multivariable methods | IV table with one `beta_x`/`se_x` column pair per exposure, plus the outcome, and optional LD matrix. |
| `region` | `cis_region_correlated`, `regional_full_sumstats` | All variants in a gene region for exposure and outcome, plus the region LD matrix. |
| `genome_wide_pair` | `genome_wide_full_sumstats` | Paths to both full harmonised variant tables plus LD scores; the adapter streams them. |

**C5.1 Harmonised IV table columns.** `variant_id`, `ea`, `oa`, `beta_x`, `se_x`, `eaf_x`, `beta_y`, `se_y`, `eaf_y`, `p_x`, `p_y`, `n_x`, `n_y`, `harmonise_action` (`aligned | flipped | palindrome_inferred | proxy`), `mr_keep` (bool).

**C5.2 Allele and sign convention.** Outcome effects are aligned to the exposure effect allele. The table is **not** re-oriented to make `beta_x` positive. Methods that need that orientation (MR-Egger and its variants, some mode and median implementations) MUST orient inside the adapter, by multiplying `beta_x`, `beta_y` by `sign(beta_x)`. This keeps one stored table valid for every method.

**C5.3 Palindromes.** Default handling is inference from allele frequency; palindromic variants with minor allele frequency above 0.42 are treated as ambiguous and set `mr_keep = FALSE` (the TwoSampleMR default). The threshold is a recorded parameter.

**C5.4 LD matrices.** Stored as Parquet (long format: `i`, `j`, `r`) or a dense float64 matrix blob with variant order recorded. Signs MUST be aligned to the same effect alleles as the IV table. The panel, ancestry and sample size of the panel are recorded. Matrices are checked for symmetry and unit diagonal; non-positive-definite matrices are repaired only by a recorded method (nearest positive-definite) and flagged.

## C6. Method adapter interface

Each method is a registered adapter object. The reference R implementation is in `contracts/r/adapter.R`.

**C6.1 Manifest (static).**

| Field | Meaning |
|---|---|
| `id` | Inventory method ID. |
| `adapter_version` | Semantic version of this adapter. Bump on any change that can alter results. |
| `label` | Name shown in the UI. |
| `input_kind` | One of the five in C5. |
| `estimands` | List of estimands it returns (C7.2), e.g. `causal_effect`, `pleiotropy_intercept`. |
| `requires` | `min_ivs` (hard minimum from the inventory), `needs_ld`, `needs_eaf`, `needs_n`, `needs_full_sumstats`, `min_exposures`, `max_exposures`. |
| `packages` | Named list of package names and exact versions or commit SHAs it calls. Must match the lockfile. |
| `stochastic` | TRUE if it uses random numbers (bootstrap, MCMC, simulation-based P values). |
| `supports_checkpoint` | TRUE if it can save and resume partial state. |
| `runtime_class` | `fast | moderate | slow | very_slow` from the inventory. |
| `params` | Parameter schema: name, type, default, allowed range. Only parameters listed here can be set. |

**C6.2 Functions.**

- `check(input, params)` → `Verdict`. Data-stage feasibility (the selection-stage rules run earlier from metadata). Returns `status` (`ok | warn | block`), `codes` (C11) and messages. MUST be cheap and side-effect free.
- `run(input, params, ctx)` → `MethodResult` (C7). Called only when `check` did not block.

**C6.3 Run context `ctx`.** `seed` (integer, C8.4), `scratch_dir` (empty folder the adapter may write to), `progress(fraction, message)`, `is_cancelled()`, `checkpoint_save(object)` and `checkpoint_load()` (no-ops unless `supports_checkpoint`), `log(level, message)`, `time_budget_s`.

**C6.4 Adapter rules.** An adapter MUST NOT: open network connections; read or write outside `scratch_dir` (inputs are passed as objects or read-only paths); change global options without restoring them; call `set.seed` itself (the runner sets it); print to the console as a way to return values; catch and hide errors from the underlying package. It MUST convert package errors and warnings into C11 codes and keep the original message text in `messages`.

## C7. Method result

One `MethodResult` per task (`schemas/method_result.schema.json`).

**C7.1 Envelope.** `contract_version`, `task_id`, `fingerprint`, `method_id`, `adapter_version`, `status` (`ok | warn | failed | blocked`), `codes`, `messages`, `estimates` (table C7.2), `diagnostics` (table C7.3), `variant_level_blob` (optional Parquet blob for per-variant output used by plots: leave-one-out, single-SNP, radial weights, outlier flags), `provenance` (C7.4).

**C7.2 Estimates table.** One row per estimand per exposure.

| Column | Meaning |
|---|---|
| `estimand` | `causal_effect`, `direct_effect`, `indirect_effect`, `total_effect`, `pleiotropy_intercept`, `proportion_mediated`, `colocalisation_pp`, `genetic_correlation`, `reverse_effect`. |
| `exposure_id`, `outcome_id` | Dataset IDs. |
| `scale` | `beta`, `log_or`, `log_hr`, `probability`, `correlation`, `proportion`. Causal effects are stored per unit of the exposure scale on the outcome scale. |
| `estimate`, `se` | `se` may be `null` when the method gives only an interval. |
| `ci_lower`, `ci_upper`, `ci_level` | As produced by the method (default 0.95). |
| `ci_method` | `wald` (normal quantile) `| wald_t` (t quantile, used when the method's P value is t-based, e.g. MR-Egger with k−2 df) `| profile | bootstrap | bayesian | robust | as_reported`. Asymmetric intervals are kept, never rebuilt from `se`. |
| `p` | As produced by the method, `null` if the method has none. |
| `n_ivs_used` | Instruments actually used after any outlier removal. |

**C7.3 Diagnostics table.** Long format: `name`, `value`, `df`, `p`, `note`. Standard names: `cochran_q`, `rucker_q`, `i2`, `i2_gx`, `egger_intercept`, `presso_global_p`, `presso_n_outliers`, `presso_distortion_p`, `steiger_direction`, `mean_f`, `min_f`, `n_removed_outliers`, `pp_h3`, `pp_h4`, `heidi_p`. Methods may add others prefixed with their ID.

**C7.4 Provenance.** `engine_version`, `r_version`, `os` (`macos-arm64 | macos-x64 | windows-x64 | linux-x64`), `packages` (name → version or SHA actually loaded), `seed`, `rng_kind`, `blas` (name), `threads`, `started_at`, `finished_at`, `elapsed_s`, `peak_memory_mb` where available.

**C7.5 Presentation rules (not stored).** Odds ratio = `exp(estimate)` with `exp` of the interval ends. Effects for binary exposures are also shown per doubling of the odds (multiply log-scale effect by ln 2). Per-SD conversion uses the units passport and is labelled as converted.

## C8. Task graph, fingerprints and resume

**C8.1 Task types.** `download`, `import`, `ancestry`, `instruments`, `harmonise`, `ld_matrix`, `ldsc`, `check`, `method`, `figure`, `table`, `report`. Tasks form a directed acyclic graph; a task depends only on blobs and cards produced by earlier tasks.

**C8.2 States.** `pending → running → done | failed | blocked | cancelled`. On start-up any task found `running` is set back to `pending` (the process that ran it is gone). `done` is final unless the fingerprint changes.

**C8.3 Fingerprint.** SHA-256 of the canonical JSON (C8.5) of: `contract_version`, `task_type`, the blob hashes of every input, the task parameters, the adapter ID and version (method tasks), the package versions listed in the adapter manifest, and the OS family. Secrets are never part of a fingerprint. A `done` task whose stored fingerprint equals the newly computed one is reused; otherwise it is re-run. This gives FR-17 (re-run only what changed) and FR-23 (identical after resume) from one rule.

**C8.4 Seeds and random numbers.** Before every task the runner calls `RNGkind("Mersenne-Twister", "Inversion", "Rejection")` and `set.seed(seed)`, where `seed` is the first 31 bits of SHA-256 of the task ID. The seed therefore depends only on what the task is, not on run order or worker. BLAS and OpenMP threads are fixed to 1 inside workers so floating-point summation order cannot vary. Parallelism happens between tasks, never inside one.

**C8.5 Canonical JSON.** UTF-8, object keys sorted by code point, no insignificant whitespace, integers written plainly, doubles written with C `%.17g` (17 significant digits; `-0` written as `0`), control characters escaped as `\u00XX`, `NA` and `NULL` written as `null`, booleans as `true`/`false`. Inf, NaN, factors and named atomic vectors are refused (use a named list for objects; wrap a length-1 vector in `I()` to force an array). Strings that are already valid UTF-8 are hashed byte-for-byte without re-encoding, so the result does not depend on the session locale. One function in the engine produces it; nothing else may hash JSON.

**C8.6 Reproducibility tolerances.** Same OS and same fingerprint: identical stored values. Across macOS and Windows: deterministic methods agree within relative difference 1e-8 on estimates and SE; stochastic methods within the Monte Carlo error the method reports, and the validation suite checks this on fixed datasets. A larger difference is a release-blocking bug.

**C8.7 Write protocol.** A finished task is committed in this order: (1) the worker writes each output to `blobs/tmp/<random>`, flushes it to disk, then renames it to `blobs/<aa>/<sha256>`; if that name already exists the temporary file is deleted, because identical content is already present; (2) the worker returns the blob hashes and the small JSON result to the main process; (3) the main process writes the task row, result JSON and state `done` in one SQLite transaction. A crash before step 3 leaves only unreferenced blobs, which a later clean-up removes. A crash can never leave a task marked `done` without its data.

**C8.8 Resume.** A reused task returns only its fingerprint and the pointer to its stored result; callers read the result from the store, never from memory. Resume = load the graph, reset `running` to `pending`, recompute fingerprints, skip `done` tasks with matching fingerprints, schedule the rest. Long tasks with `supports_checkpoint` reload their last checkpoint blob. Network tasks wait while offline; local tasks continue.

## C9. Project folder and database

```
<Project>.mrgui/
  project.sqlite            metadata, tasks, results (WAL mode)
  blobs/aa/<sha256>         content-addressed Parquet, JSON, PNG, PDF
  blobs/tmp/                temporary writes only
  exports/                  user-facing files (reports, figures, workbooks)
  logs/                     rotating text logs, secrets redacted
```

Reference data (LD panels, LD scores, catalogue) lives in a shared app-data folder with the same blob layout, so it is downloaded once per computer.

SQLite runs with `journal_mode=WAL`, `synchronous=FULL` (a power cut must not lose a committed task), `foreign_keys=ON`. Tables: `meta`, `datasets`, `instrument_sets`, `tasks`, `task_deps`, `results`, `blobs`, `ui_state`, `audit_log`. The DDL is in `contracts/sql/schema_v1.sql`. Schema changes use numbered migrations recorded in `meta`.

The project folder MUST NOT be inside a cloud-sync folder (iCloud Drive, OneDrive, Dropbox, Google Drive) or on a network drive, because sync tools can copy a SQLite database mid-write and corrupt it. The app detects these locations and refuses with an explanation; exports may be saved anywhere.

## C10. Desktop shell and engine connection

1. The shell starts the bundled R engine as a child process with a UTF-8 locale (`LANG=en_US.UTF-8` on macOS; R 4.2+ on Windows is UTF-8 natively), and every file reader in the engine declares `encoding = "UTF-8"`. It passes a random 256-bit secret and a one-time login token on **standard input**, never as command-line arguments or environment variables, which other local programs can read.
2. The engine binds an HTTP server to `127.0.0.1` on a port chosen by the operating system, then prints one JSON line to standard output: `{"ready": true, "port": <n>, "pid": <n>}`.
3. The shell opens its window at `http://127.0.0.1:<port>/auth?t=<one-time token>`. The engine checks the token once, sets an `HttpOnly; SameSite=Strict` session cookie holding the secret, and redirects to `/`. The token is then invalid.
4. Every HTTP request and WebSocket upgrade MUST carry the cookie, a `Host` header equal to `127.0.0.1:<port>`, and for WebSockets an `Origin` header equal to `http://127.0.0.1:<port>`. Anything else gets 403. The Host check stops DNS-rebinding attacks; the Origin check stops a web page open in the user's normal browser from connecting to the engine.
5. The shell watches the engine; if it exits unexpectedly the shell shows a restart screen, and the engine resumes from the database (C8.8). When the window closes, the shell asks the engine to shut down cleanly and kills it after a timeout.

Implementation note: Shiny serves through httpuv, and these checks need a small layer in front of Shiny's handlers for both HTTP and WebSocket traffic. Proving this layer works on both operating systems is a week-1 prototype exit criterion; if it cannot be done cleanly, the fallback is for the shell to proxy all traffic and attach the checks itself.

## C11. Codes for blocks, warnings, drops and errors

Codes are stable strings used in verdicts, results, instrument drop lists and logs. The UI maps each code to a plain-language message; the report maps it to a limitation sentence.

| Family | Examples |
|---|---|
| `FEAS_*` (block or warn before running) | `FEAS_TOO_FEW_IVS`, `FEAS_NEEDS_LD`, `FEAS_NEEDS_FULL_SUMSTATS`, `FEAS_ANCESTRY_MISMATCH`, `FEAS_MULTI_ANCESTRY_LD`, `FEAS_NEEDS_EAF`, `FEAS_NEEDS_N`, `FEAS_TOO_FEW_EXPOSURES`, `FEAS_TOO_MANY_IVS` |
| `DATA_*` (warnings about the data) | `DATA_WEAK_IVS`, `DATA_SAMPLE_OVERLAP`, `DATA_UNITS_UNKNOWN`, `DATA_BINARY_EXPOSURE`, `DATA_EAF_MISMATCH`, `DATA_HETEROGENEITY`, `DATA_LD_NOT_PD` |
| `DROP_*` (why a variant left an instrument set) | `DROP_NOT_IN_OUTCOME`, `DROP_PALINDROMIC_AMBIGUOUS`, `DROP_ALLELE_MISMATCH`, `DROP_CLUMPED`, `DROP_OUTLIER`, `DROP_STEIGER_REVERSE`, `DROP_NO_PROXY` |
| `RUN_*` (problems while running) | `RUN_PACKAGE_ERROR`, `RUN_PACKAGE_WARNING`, `RUN_NOT_CONVERGED`, `RUN_TIMEOUT`, `RUN_CANCELLED`, `RUN_OUT_OF_MEMORY` |
| `SYS_*` (environment) | `SYS_OFFLINE`, `SYS_DISK_FULL`, `SYS_CLOUD_FOLDER`, `SYS_TOKEN_MISSING`, `SYS_QUOTA_EXHAUSTED` |

## C12. Security model

**Assets:** the user's OpenGWAS token, their unpublished results, the integrity of results used in publications, and the user's computer.

**Threats and controls:**

| Threat | Control |
|---|---|
| Another program or a web page connects to the engine | C10: loopback only, secret cookie, Host and Origin checks, one-time login token. |
| Token leakage | Stored only in the OS credential store; read into memory just for network tasks; never in fingerprints, task parameters, logs, exports or diagnostics zips; logs pass through a redaction filter. |
| Malicious or malformed data file | Strict typed parsers with size limits; no evaluation of file content; rejection of individual-level shapes; all paths resolved and confined to the project or reference folder. |
| Archive extraction (`zip-slip`) | Every extracted path is checked to stay inside its target folder; links in archives are refused. |
| Injection into reports and the UI | Trait names, study titles and any text from metadata are HTML-escaped before display and in HTML/DOCX reports. |
| Formula injection in exported spreadsheets | Text cells that start with `=`, `+`, `-`, `@`, tab or carriage return are prefixed with an apostrophe in CSV/XLSX exports; numeric cells are written as numbers, never text. |
| Command injection via external tools | PLINK and other binaries are called with an argument vector (no shell), from their bundled absolute paths. |
| Tampered downloads | HTTPS only; allow-listed hosts (S-04); checksums verified where the source publishes them and recorded always. |
| Tampered app updates | Updates are verified with the Tauri updater signature (free key pair; private key held only in the build pipeline secrets) before installation. |
| Supply chain | Exact versions and commit SHAs in the lockfile; external binaries pinned by checksum; CI builds from a clean checkout. |
| Result integrity | Content-addressed blobs, fingerprints and the audit log let any number in a report be traced to its inputs and code version. |

## C13. Versioning and conformance

- `contract_version` follows semantic versioning. Adding an optional field is a minor change. Removing or renaming a field, changing a meaning, or changing canonical JSON or fingerprint rules is a major change and invalidates fingerprints, so affected tasks re-run.
- Every adapter ships with a conformance test: it runs on the fixed test datasets in the validation suite, checks the result against the schema, and compares estimates with the reference values from the underlying package within C8.6 tolerances. An adapter that fails its conformance test cannot be enabled in a release build.
- The reference adapter `ivw_mre` in `contracts/r/` shows the pattern every other adapter copies.

## Verification status (2026-10-29)

The reference code in `contracts/` was run on R 4.3.3 (Linux) under the C, C.UTF-8 and en_US.UTF-8 locales: 37 of 37 tests pass in each.

- `ivw_mre` and `egger` adapters reproduce TwoSampleMR (pinned commit `5e93ee7`, functions `mr_ivw_mre` and `mr_egger_regression` loaded from that source file) for estimate, SE, P, intercept and Q, to 1e-12 relative.
- Both adapters give the same estimates when the effect allele of some variants is swapped, confirming the C5.2 sign convention (Egger orients inside the adapter).
- Blocking (`FEAS_TOO_FEW_IVS`), weak-instrument warning, parameter validation, fingerprint reuse, re-run on changed parameters, crash recovery (`running` → `pending`), orphan blobs after a simulated crash, WAL mode and byte-identical re-runs in a fresh project all behave as specified.
- Security helpers: formula injection, HTML escaping, zip-slip and absolute paths, token redaction.
- Results, blocked results and adapter manifests validate against the JSON Schemas; a deliberately invalid result is rejected; the example dataset card validates.

Two defects were found and fixed while testing: `formatC(..., format = "g")` padded numbers with spaces (now C `%.17g`), and `enc2utf8()` mangled non-ASCII text in a C locale (now locale-independent, and the shell sets a UTF-8 locale). Not yet tested: macOS and Windows, conda-forge R, Parquet blobs, the C10 connection checks, and adapters other than `ivw_mre` and `egger`.

## Open items for the week-1 prototype

1. Prove the C10 request checks in front of Shiny on macOS and Windows, or switch to the shell-proxy fallback.
2. Confirm `synchronous=FULL` write throughput is acceptable with many small method tasks (expected yes; batch commits are the fallback).
3. Confirm single-threaded BLAS in the bundled R on both operating systems and measure cross-OS differences on the validation datasets.
4. Decide Parquet reader/writer (`arrow` vs `nanoparquet`) based on bundle size and Windows build reliability.
