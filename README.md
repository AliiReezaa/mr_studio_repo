# MR Studio (MR GUI)

Free, local, point-and-click desktop app for Mendelian randomization with public GWAS summary data, for macOS and Windows.

Designed and developed by Ali Reza Soleymani Taloubaghi (all code and MR methodology). Methods and validation review: Dr. Danial Habibi (Babol University of Medical Sciences).

| Path | What it is |
|---|---|
| `requirements_spec.md` | Requirements specification |
| `project_plan.md` | Phased plan, milestones, effort estimate |
| `contracts_v1.md`, `contracts/` | Engineering contracts v1.0.0 (frozen) and reference code (data formats, adapters, results, resume, security) |
| `contracts_v1_1_draft.md` | Contracts v1.1.0 additive changes (draft) |
| `prototype/` | Week-1 packaging prototype: desktop shell, R engine, build pipeline |
| `mr_methods_inventory.md`, `mr_inventory.json`, `*.csv` | Inventory of MR methods, tools and feasibility rules |
| `data_limits_solutions.md` | Solutions for public-data limits |

Licence: MIT (see `LICENSE`). The installer also bundles R, R packages and PLINK, each under its own licence (`THIRD_PARTY_LICENSES.md`). How to cite: see `CITATION.cff`. Validation of every method: `validation/README.md`.
