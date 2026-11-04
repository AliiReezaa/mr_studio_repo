# MR GUI engineering contracts: v1.1.0 (DRAFT delta to v1.0.0)

Draft 2026-11-01. Not frozen. Read together with `contracts_v1.md` (v1.0.0, frozen) and `requirements_spec.md` v1.6 Section 4A (FR-52 to FR-117). This file lists only what changes. Anything not listed is unchanged.

## 0. Rules for this version
1. **Minor version, additive only.** Every new field is optional. Canonical JSON (C8.5) and the fingerprint composition (C8.3) do not change, so v1.0.0 fingerprints stay valid and the 37 existing contract tests must pass unchanged.
2. Everything stored for the new features lives in new optional fields or new task types, never inside the fingerprint inputs.
3. No method is enabled in a release build without all three: a parameter schema, a script emitter and a passing replay test (C13 additions).

## C2.1 QC flag bits (additive)
New: **256** other allele inferred from the reference panel; **512** EAF filled from the reference panel; **1024** SE or effect derived (from P, from OR limits, or from a z-score).
Not added: spec FR text proposed 2048 for "rsID mapped from position". Bit **16** already has exactly that meaning, so 2048 is dropped and the spec wording is corrected.

## C3 Dataset card (additive optional fields)
`sex` (`both | female | male | unknown`); `age_range` (`{min, max}` or null); `design` (`population | family | sibship | trio | unknown`); `socially_patterned` (bool or null); `selection_role` (`discovery | replication | selection_only | estimation_only | both | unknown`); `upload_recipe` (original file name, size, SHA-256, column mapping, facts form; present when `source = file`); `missing_ledger` (counts of rows kept, filled, derived and dropped, per policy rule).
`source` keeps the value `file` for user uploads (no new enum value, so v1.0.0 readers still work).

## C6 Method adapter (additive)
**C6.1 manifest, `params` entries gain** (all optional for old adapters, required for new ones):
`label`; `hint` (plain language, at most 280 characters, the same text is shown in the UI and written as a comment in the script); `package_arg` (the argument name in the R package); `kind` (`choice | int | float | bool | text`); `choices` or `range`; `default`; `conventional` (bool: only `true` parameters can be shown, so only conventionally adjustable arguments are offered); `advanced` (bool: hidden behind "More options"); `evidence` (PMID, DOI or package documentation reference for why it is conventional); `affects` (list of `estimate | se | ci | p | diagnostic | runtime`).
New manifest field **`column_requirements`**: list of `{column, level: required | recommended | optional, why, if_missing: {action: block | fill | degrade | ignore, code}}`. The importer and the feasibility engine read this list to tell the user what a missing column costs.

**C6.2 new function `emit_script(input_spec, params)`** returns a `ScriptBundle`: `method_script` (text), `helpers` (optional text), `packages` (name to version), `outputs` (list of `{file, role}`), `seed`. Rules:
1. Pure and deterministic: identical inputs give identical bytes. No date, time, user name, absolute path or random value in the text.
2. Every statistical step is a call to a public, versioned package function with **all arguments written out**, including defaults the user did not change, so a later change of a package default cannot alter results silently. No statistical code is hidden in app-only code. Sole exception: SMR and HEIDI (no public R implementation), inline in full, validated against the reference program, and marked `inline_code: true` in the manifest.
3. Reads only the CSV files of the pair's `01_data/` through relative paths; writes only inside the pair folder. A guard at the top stops with a clear message if the working directory is wrong.
4. First lines set `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`, `set.seed(seed)` and single-threaded BLAS and OpenMP, as in C8.4.
5. Results go to two CSV files per method in `04_tables/`, `<method_id>_estimates.csv` and `<method_id>_diagnostics.csv` (corrected from the first draft, which had one pair of fixed names; several methods share a pair folder), with the column names of C7.2 and C7.3 and doubles written with 17 significant digits. This is the only way back into the engine.
6. Format: fixed section headers, one comment line per argument (the hint text), a closing session-information block.

**C6.2 `run()` redefined for adapters with `emit_script`:** write the bundle to the scratch folder, execute it in a fresh `R --vanilla` process, read the two result CSVs. There is no second code path. For `ivw_mre` and `egger` the v1.0.0 in-process `run` is kept only as a test oracle (it is the code checked against TwoSampleMR).

**C6.5 replay check (engine function, not an adapter function).** When a pair finishes: copy the delivered files to a temporary folder, start a fresh `R --vanilla`, run `00_run_all.R`, compare every numeric output (estimate, SE, CI limits, P, Q, intercept, number of instruments) with the stored values. Same computer: identical to 15 significant digits. Across operating systems: relative difference 1e-8 (C8.6). Status values: `identical | within_tolerance | differs | failed | not_run`. Stochastic methods must use a fixed seed and one thread; if a method cannot be made identical it is labelled "reproducible to Monte Carlo error" with its tolerance stated.

## C7 Method result (additive optional fields)
`script_sha256`; `analysis_data_sha256`; `params_effective` (every parameter with its value and an `is_default` flag); `customised` (bool); `default_task_id` (the same method with all defaults, used by "Compare with default"); `user_modified` (bool, true when the script was edited and re-run, FR-87); `replay` = `{status, checked_at, r_version, os, max_rel_diff, tolerance, n_values_compared}`.
New standard diagnostic names: `outlier_scenario` (which detector's removal set was used) and `n_ivs_before_outlier_removal`.

## C8 Tasks (additive)
New task types: `script` (emit and run), `replay`, `outlier` (detect, compare scenarios, re-run), `user_script` (edited script; its inputs include the hash of the script bytes, so it has its own fingerprint and never overwrites the generated one). Replay results are stored with the result row and are not part of any fingerprint. Fingerprint composition and C8.5 are unchanged.

## C9 Project folder (replaced layout, same database rules)
Visible study folder chosen by the user (default `~/MR Studio`, folder picker, FR-52):
```
<Study>/
  01_study_plan/   02_datasets/   03_pairs/P01_<exp>__to__<out>/{01_data, 02_scripts, 03_console, 04_tables, 05_figures, 06_method_reports, 07_replay_check, 08_outlier_analysis}
  04_across_pairs/   05_report_and_supplement/   06_reproducibility/   archive/
  MANIFEST.csv
  .engine/         project.sqlite, blobs/, logs/   (hidden, unchanged from C9 v1.0.0)
```
Rules: every file in the numbered folders is a real file (no symlinks); written atomically (temporary file in the same folder, flush, rename); listed in `MANIFEST.csv` with `relative_path, sha256, bytes, role, task_id, created_by (app | user), version`. A file changed or deleted by the user is detected through the manifest and offered for regeneration, never overwritten silently (the new version is saved beside it). Re-running with changed settings moves earlier outputs to `archive/run_<datetime>/`. Names follow FR-54 (ASCII slug, at most 40 characters); the full path must stay below 240 characters (Windows limit).
The cloud-sync and network-drive refusal of v1.0.0 C9 still applies to the whole study folder. DDL migration 2: `results` gains the C7 columns; new tables `manifest` and `script_versions`.

## C11 New codes (proposed names)
`IMPORT_NO_LOCUS`, `IMPORT_NO_EA`, `IMPORT_NO_SE`, `IMPORT_NO_BETA`, `IMPORT_TOO_MANY_DROPS`, `IMPORT_UNKNOWN_BUILD`, `IMPORT_HEADER_AMBIGUOUS`; `SCRIPT_RUN_FAILED`, `SCRIPT_NOT_DETERMINISTIC`, `SCRIPT_USER_MODIFIED`; `REPLAY_FAILED`, `REPLAY_DIFFERS`, `REPLAY_TOLERANCE_ONLY`; `OUTLIER_CAP_REACHED`, `OUTLIER_TOO_FEW_IVS`, `OUTLIER_NO_AGREEMENT`; `DATA_SEX_MISMATCH`, `DATA_SELECTION_EQUALS_ESTIMATION`, `DATA_SOCIALLY_PATTERNED` (info), `FEAS_REPLICATION_OVERLAP`; `SYS_PATH_TOO_LONG`.

## C13 Conformance additions (every adapter)
1. Parameter schema valid: every default is valid, every choice runs, every `conventional` parameter has `evidence`.
2. `emit_script` determinism: two emissions give identical bytes.
3. Script result equals the in-process oracle within 1e-12 (reference adapters).
4. Replay check passes on the same OS (identical to 15 digits for deterministic methods).
5. CI: replay on macOS and Windows; cross-OS estimates within 1e-8.
6. Folder tests: slug collisions, path length on Windows, manifest equals files on disk, atomic write after a simulated crash.

## Open technical points (to settle while building the first two adapters)
1. `renv` with the conda-forge runtime: whether `renv::restore()` on a user's ordinary R reproduces our package versions (GitHub-only packages such as TwoSampleMR need a pinned commit).
2. Whether explicit default arguments differ between the pinned package versions on macOS and Windows.
3. Cost of the replay check for slow methods (run once at final run, not on every exploratory run).

## Decisions recorded (Mina, 2026-11-01)
- D1 Outside-app rerun: accepted. Users install the pinned package versions once (internet needed); cross-OS agreement is 1e-8, identical to 15 digits only on the same computer. The in-app replay check is the strict proof.
- D2 Output folder: default `~/MR Studio`; folders inside iCloud, OneDrive, Dropbox, Google Drive and network drives are refused.
- D3 Replay check: automatic for final runs; an optional switch skips it for quick exploratory runs. Never silently skipped for a finalised study.
- D4 No-shortcut rule adopted: a method is enabled in a release only with a parameter schema, a script emitter and a passing replay test.

## Verification status (2026-11-01, Linux, R 4.3.3, TwoSampleMR 0.7.13 at commit 5e93ee7)
Prototype code: `contracts/r/script.R`, `contracts/r/adapters/reference_adapters_v11.R`, tests `contracts/r/tests/test_scripts.R`; `adapter.R` gained two optional arguments only.
- The original 37 contract tests pass unchanged (the change is additive).
- 24 new tests pass: schema checks (hints, evidence, valid defaults); emit_script gives identical bytes twice and contains no date, user or absolute path; every argument written out; no app code in the script; script run in a fresh `R --vanilla` equals the in-process oracle within 1e-12 for `ivw_mre` (ci_level 0.9) and `egger` (slope and intercept); the replay check reports `identical` (27 numeric values, difference exactly 0), `differs` for an edited script, `failed` for a script that errors, `identical` again once restored; wrong working directory stops with a clear message; the analysis CSV round trip is within 2.3e-16 relative.
- Not yet tested: macOS and Windows (the CI runtime has no TwoSampleMR yet, see below), stochastic methods, `renv` restore on a user's computer, the engine running scripts through the task runner (the prototype runs them directly).

## Finding: the bundled R needs TwoSampleMR and its GitHub-only dependencies
TwoSampleMR 0.7.13 imports MRMix, MRPRESSO and RadialMR, none of which is installable from the standard registries for R 4.3 (RadialMR was no longer on CRAN; MRMix and MRPRESSO are GitHub-only). They had to be built from pinned GitHub sources (MRMix 56afdb2, MR-PRESSO 3e3c92d, RadialMR from GitHub HEAD at the time) and several others from Ubuntu binaries. The same will be true for the conda-forge runtime and for `renv` on a user's machine. Options to decide before M1 continues: (a) build these packages from pinned commits inside CI and ship them in the runtime (recommended; matches the lock-file idea), (b) have the scripts call minimal public functions only for the two reference methods, (c) vendor patched copies (not recommended: it would be hidden statistical code).

### Packaging decision applied (2026-11-01)
Option (a) was chosen: the pipeline installs TwoSampleMR, MRMix, MRPRESSO and RadialMR from pinned commits into the bundled runtime (`.github/workflows/prototype.yml`, step "Install the pinned MR packages into the runtime"); their CRAN dependencies come from conda-forge (`prototype/ci/runtime.yml`).
Pins: TwoSampleMR 5e93ee7d8e1f67581acb8134299759389d73efba, MRMix 56afdb2bc96760842405396f5d3f02e60e305039, MR-PRESSO 3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87, RadialMR 7ca35c052489001890347871e38998ca53a02b03.
Rehearsed on Linux with conda-forge R 4.5.3: all four install (pure R, no compiler needed), 37 of 37 original and 24 of 24 new tests pass. macOS and Windows are verified only when the pipeline runs.
Cost: the conda dependencies pull in pandoc (about 160 MB) through rmarkdown, which TwoSampleMR imports. Check the installer sizes after the first run before deciding whether to strip it.

### Verification update (2026-11-01, pipeline run 8)
Script generation, the run in a fresh R and the replay check passed on macOS Apple Silicon, macOS Intel and Windows (CI). Across the three systems the generated scripts and data files agree (scripts byte-identical except for one test artefact, fixed in the test), and result tables agree within 1.3e-15 relative. Two launcher issues found and fixed on the way: the bundled R launcher kept the build folder path, and R's launcher breaks when the install folder name contains a space (for example MR Studio.app); `relocate_macos.py` now patches both, and `script.R` passes RHOME to Rscript.

### Engine integration (2026-11-01, first version, Linux only so far)
`contracts/r/study.R` and `contracts/sql/schema_v2.sql` implement C9 and the C8 `script` task: study folder tree, refusal of cloud-sync and network folders (D2), path length limit (SYS_PATH_TOO_LONG), atomic writes, MANIFEST.csv with hashes and versions mirrored in the database, archive/run_<datetime>/ when settings change, SCRIPT_USER_MODIFIED instead of overwriting a user's edit (or `on_modified = "archive"`), automatic replay check with an optional switch (D3, recorded as `not_run`), C7 fields in the results row. 34 tests, run in a runtime moved into a folder with a space in its name. Not yet verified on macOS and Windows.
Reports (HTML, PDF, Word) are required in v1 (Mina, 2026-11-01), so pandoc stays in the runtime; a PDF engine is an open choice (a small one such as typst versus a LaTeX distribution).

### Diagnostic-only methods (2026-11-03, Round 15, additive, CONTRACT_VERSION stays 1.0.0)
Some methods test an assumption and give no causal estimate (the first is the Steiger directionality test, Hemani, Tilling, Davey Smith 2017). `mrgui_adapter()` gains the optional argument `diagnostic_only` (default FALSE; such an adapter has `estimands = character(0)`). A result of such a method carries `"diagnostic_only": true` and an empty `estimates` array; `method_result.schema.json` then requires at least one diagnostic instead of at least one estimate. The field is written only for diagnostic-only methods, so every existing result, fingerprint and saved study is unchanged. The generated script still writes a header-only `<id>_estimates.csv`, and the replay check compares the diagnostics.

### Steiger for binary traits and as part of IVW (2026-11-05, Round 16, additive, CONTRACT_VERSION stays 1.0.0)
The Steiger adapter goes to `adapter_version` 1.2.0. A pair input may carry the optional field `exposure_scale` (`beta` or `log_or`, default `beta`), next to the existing `outcome_scale`. For a binary side the Steiger adapter reads three extra columns of its own analysis data file, `ncase_<s>`, `ncontrol_<s>` and `prevalence_<s>` (s = x for the exposure, y for the outcome), and computes variance explained with `TwoSampleMR::get_r_from_lor(..., model = "logit", correction = FALSE)` (Lee et al. 2012, PMID 22714935); the presence of `ncase_<s>` is what the generated script reads. A continuous side is unchanged and gets no extra columns, so analysis files of continuous pairs and of all other methods are byte-identical to before. New verdict codes: `FEAS_STEIGER_NEEDS_CASES`, `FEAS_STEIGER_NEEDS_PREVALENCE`, `FEAS_STEIGER_NEEDS_EAF` (`FEAS_STEIGER_BINARY` is no longer used). Two diagnostics are added to every new Steiger result: `steiger_exposure_binary` and `steiger_outcome_binary` (value 1 or 0; the note states the method and, for a binary side, the cases, controls and prevalence). Results saved before this change have 6 Steiger diagnostics and are read as continuous.
On screen the method is shown as part of IVW; in the engine it is still the separate method `steiger`, so the contract and the replay check are unchanged.
