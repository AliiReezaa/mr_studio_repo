# MR GUI contracts v1 — reference files

Specification: `../contracts_v1.md` .

| Path | What it is |
|---|---|
| `r/canonical.R` | Canonical JSON, SHA-256 hashing, task fingerprints and seeds (C8.3–C8.5) |
| `r/adapter.R` | Adapter constructor, registry, parameter validation, verdicts, result rows, manifest and result serialisation (C6, C7) |
| `r/runner.R` | Blob store, project database, method-task runner, reuse and crash recovery (C8, C9) |
| `r/security.R` | Formula-injection guard, HTML escaping, safe path joining, secret redaction (C12) |
| `r/adapters/reference_adapters.R` | Reference adapters `ivw_mre` and `egger`: the pattern every adapter copies |
| `r/tests/test_contracts.R` | 37 conformance and contract tests |
| `sql/schema_v1.sql` | Project database schema |
| `schemas/*.schema.json` | JSON Schemas: method result, dataset card, adapter manifest |
| `examples/dataset_card_example.json` | Example dataset card (illustrative values) |
| `ref/TwoSampleMR_mr_R_5e93ee7.R` | TwoSampleMR `R/mr.R` at commit 5e93ee7, used as the conformance reference |

Run the tests from this folder (needs R with `digest`, `DBI`, `RSQLite`):

```
Rscript r/tests/test_contracts.R <output_dir>
```

Validate JSON outputs against the schemas with any JSON Schema 2020-12 validator.
