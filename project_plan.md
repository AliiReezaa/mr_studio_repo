# MR GUI — Phased project plan (v1)

Version 1 · 2026-10-29 · companion to requirements_spec.md v1.5 and working_policy_token_and_model.md. Hours are realistic focused effort (not optimistic), including about 15% contingency. Total about 700 h (range 550 to 900). Your own time is extra: about 120 to 180 h for decisions, review and macOS testing, plus about 30 to 50 h of Windows testing by your friend.

## Milestones

| Milestone | What you can do when it is done | Cumulative hours |
|---|---|---|
| M0 Prototype and contracts | Installer builds on macOS and Windows; one analysis (IVW) runs end to end; data, adapter and result contracts frozen | 60 to 90 |
| M1 MVP (publishable core) | Install, pick exposures/outcomes, feasibility engine, about 15 core methods, run with resume, basic export | 220 to 280 |
| M2 Full core | All 40 core analyses, STROBE-MR, draft text, journal-ready figures, supplement workbook | 450 to 560 |
| M3 Drug-target and bidirectional | Colocalisation, SMR/HEIDI in R, cis-MR, bidirectional and genetic correlation | 520 to 650 |
| M4 Hardening and release | Cross-OS validation, simulation calibration, installers, auto-update, documentation | 550 to 900 |

## Phases, model routing and parallelism

| Phase | Work | Hours | Model | Parallel with |
|---|---|---|---|---|
| 0a | Week-1 packaging prototype: Tauri vs Electron, bundled R on both OS, CI builds | 30 to 45 | Opus (design, debugging), Sonnet (scripts) | none (blocks everything) |
| 0b | Contracts: data schema, adapter interface, result object, task graph, security model | 15 to 25 | Opus | none |
| 1 | Data layer: catalogue, OpenGWAS token and cache, harmonisation, build and liftover, units passport, ancestry, PLINK and LD | 60 to 90 | Sonnet; Opus review of statistics | 2, 3 |
| 2 | Feasibility engine, rule set, threshold calibration by simulation | 40 to 60 | Opus (thresholds), Sonnet (code) | 1, 3 |
| 3 | UI: build the 10 designed screens in the real app, theme system | 80 to 120 | Sonnet | 1, 2 |
| 4 | Method adapters for 49 analyses from one reviewed template plus per-method config | 110 to 160 | Sonnet; Haiku for repetitive files; Opus for MR-Horse, CAUSE, SMR/HEIDI | after 0b |
| 5 | Resume and checkpointing | 30 to 45 | Opus (design), Sonnet | 4 |
| 6 | Publication support (FR-26 to FR-40) and figure system (FQ-1 to FQ-8) | 70 to 110 | Sonnet; Opus for statistical wording | 4 |
| 7 | Drug-target, colocalisation, bidirectional, genetic correlation | 40 to 60 | Opus for statistics, Sonnet for code | 6 |
| 8 | Validation against published results, cross-OS testing, bug fixing | 60 to 100 | Opus for statistics | after 4 to 7 |
| 9 | Installers, auto-update, Diagnose, user documentation | 30 to 50 | Sonnet | 8 |

## Rules for saving time and tokens

1. Freeze contracts before splitting work; rework in parallel branches is the biggest hidden cost.
2. Use GitHub Actions Windows and macOS runners so most testing needs no manual effort.
3. Generate adapters from one template plus a config table; review the template once on Opus.
4. Search project docs for the exact section instead of re-reading whole files.
5. Write tests with the code; build the MVP first and adjust later phases from real use.

## Dependencies and risks

- Packaging R on Windows is the largest technical risk. If the week-1 prototype fails, add about 40 h and fall back to a per-user R install.
- The 78 methods with no public code stay out of v1.
- Publisher figure widths must be verified against current author guidelines when presets are built.
- Heuristic thresholds (ancestry 90%, allele frequency difference 0.2, instrument counts) need simulation before release. Done in Rounds 29 and 32 (validation/calibration/README.md): all kept; the outlier module should use detectors at 0.05 / k.

## Additions 2026-10-31 and 2026-11-01 (DRAFT, requirements_spec.md Section 4A)

| Work | Extra hours (estimate) | Model | Placed in |
|---|---|---|---|
| Visible project tree, naming, manifest, backup (FR-52 to FR-55) | 10 to 15 | Sonnet | M1 |
| Script generation, rerun kit and automatic replay check (FR-56 to FR-60) | 25 to 40 | Opus for design and the stochastic-method cases, Sonnet for code | M1 (shapes every adapter) |
| Parameter schemas and Customise panel (FR-61 to FR-69) | 35 to 55 | Opus for choosing the exposed arguments, Haiku/Sonnet for authoring the 49 schemas | core methods M1, rest M2 |
| Missing-value policy and ledger (Section D) | within importer and Phase 1 | Opus for statistical rules | M1 |
| Import wizard and column requirements (FR-70 to FR-76) | 30 to 45 | Sonnet; Opus review of allele and EAF logic | M2 |
| Formatted script, analysis report, final analysis CSV (FR-80 to FR-83) | within script generation | Sonnet | M1 to M2 |
| Outlier detection, comparison and removal module (Part H, FR-90 to FR-104) | 30 to 50 | Opus for the selection rule and the simulation calibration, Sonnet for code | M2 (after the PRESSO, Radial and Lasso adapters) |
| Reports through a bundled document renderer, column-mapping templates, Open in RStudio, edit-and-rerun, Compare with default | within the items above, about 10 to 15 extra | Sonnet | M1 to M2 |
| Additional safeguards (Part I, FR-110 to FR-117): sex, winner's curse, instrument-trait screen, within-family, replication, control outcomes | 30 to 50 | Opus for rules and wording, Sonnet for code | M2 (card fields in M1) |
| Total | about +170 to 265 | | new overall range about 710 to 1150 h |

Contracts v1.1.0 (additive, needs confirmation) must be agreed before the adapters are generated.

