# MR GUI — Verified Inventory of Mendelian Randomization Methods, Software & Resources

Version 3 · 2026-10-25 · adds **feasibility-requirement fields** and the **feasibility-rule set** used to enable/grey-out methods in the GUI. Machine-readable copies: `mr_methods.csv` (52 columns), `mr_feasibility_rules.csv`, `mr_tools_resources.csv`, `mr_inventory.json`. Data-limit solutions are in `data_limits_solutions.md`.

## 1. How this inventory was built
- **Seed list**: 11 PDFs in the MR_GUI folder, 5 GitHub repos (MRdmIVW, MRMU, XMR, LMM_MEC, MREILLS) and every chapter of the AndrewsLab UCSF MR tutorial.
- **Literature sweep**: PubMed title searches for MR methods papers (2013–2026) across methodology journals (AJHG, Nat Commun, Nat Genet, PLoS Genet, Genome Biol, HGG Adv, Sci Adv, PNAS, Bioinformatics, Brief Bioinform, Genet Epidemiol, Stat Med, Biometrics, SMMR, Biostatistics, IJE, AJE, Epidemiology, EJE, Cell Genom, Genome Res, eLife, Front Genet, Hum Genet, BMC Bioinform), plus reference lists of the source PDFs.
- **Reference check**: 194 of 210 method entries cite a PubMed record retrieved in this session (PMID + DOI). Entries marked ⚑ are statistics-journal or preprint papers not indexed in PubMed; their source is stated in the CSV column `ref_verification`.
- **Software check**: every listed GitHub repository returned HTTP 200 on 2026-10-25 (115 entries). CRAN packages were checked through the official `github.com/cran` mirror. Function names marked ✓ (67 entries) were confirmed in the package NAMESPACE, source file or tutorial code. Code quality/security was **not** audited — that is a later step before bundling anything into the GUI.
- **Not claimed**: completeness. MR methods appear monthly; the sweep favoured methodology journals and title keywords, so methods published only in applied journals or only on arXiv/bioRxiv may be missing.

## 2. Legend (columns used in the tables and CSV)
- **Data**: S = summary statistics; I = individual-level; S/I = either.
- **Design**: 1S one-sample; 2S two-sample; 3S three-sample (selection / exposure / outcome).
- **Addresses**: BP balanced pleiotropy · DP directional pleiotropy (InSIDE holds) · CP correlated pleiotropy (InSIDE violated) · UHP uncorrelated horizontal pleiotropy · WIV weak instruments · WC winner's curse · OVL sample overlap · LD correlated variants · POP population stratification · HET heterogeneous populations/datasets · MEAS measurement error · SEL selection / collider bias.
- **✓** after functions = function name confirmed. **⚑** = reference not in PubMed (see CSV for source).
- CSV-only columns useful for the GUI: `n_exposures`, `n_outcomes`, `instrument_scope`, `extra_inputs` (e.g. LD matrix, LDSC C/Ω matrices, individual data), `key_assumption`, `language`, `in_user_sources`.

## 3. Feasibility engine — schema and rules
Each method row now carries the fields the GUI needs to decide **enable / warn (amber) / block (grey-out with reason)**. Checks run in two stages: **selection stage** (metadata only, before any download) and **data stage** (after instrument extraction, clumping and harmonisation).

| Field | Meaning / allowed values |
|---|---|
| `multi_selection_behaviour` | pairwise_loop (univariable method repeated for every exposure–outcome pair) or joint_model (all selected traits in one model) |
| `runnable_tier` | r_package | python | external_binary | needs_wrapper | web_service | no_code | reference_only |
| `gui_role` | estimator, diagnostic, sensitivity, companion_analysis, network_inference, mediation, bias_correction, design_or_reference … |
| `data_required` | summary · individual · summary_or_individual · summary_plus_individual |
| `min/max_exposures, min/max_outcomes, min_traits_total` | counts for joint models ('any' = no upper limit) |
| `min_ivs_hard` | instruments below which the method cannot run (may be a formula in n_exposures) |
| `min_ivs_recommended` | below this the app warns rather than blocks |
| `max_ivs` | upper limit (e.g. heterogeneity-penalised method: 30) |
| `iv_rule_source` | package code (verified) · method definition · heuristic — validate in testing |
| `instrument_source` | clumped_independent · relaxed_threshold_independent · correlated_variants · cis_region_correlated · genome_wide_full_sumstats · allele_score |
| `needs_full_sumstats / needs_ld / needs_ldsc_matrices` | external inputs that must exist for the selected ancestry |
| `exposure_type / outcome_type` | e.g. molecular_QTL exposure, binary/count/mixed outcomes |
| `needs_eaf / needs_n` | allele frequency and sample size required |
| `population_rule / min_ancestries` | ancestry constraints (same ancestry required, ≥2 ancestries required …) |
| `overlap_handling` | assumes_none · models_overlap · corrects_overlap · n/a |
| `extra_data` | additional datasets needed (mediator, sex-stratified, family, covariate GWAS …) |
| `runtime_class` | fast · moderate · heavy (background job) |

**Instrument-count rules verified in package source code (2026-10-25):** TwoSampleMR — IVW, ML ≥2; MR-Egger, medians, modes ≥3; sign test ≥6; leave-one-out plot ≥3. MendelianRandomization — Egger and median require >2 variants; IVW defaults to random effects at ≥4; heterogeneity-penalised method refuses >30 variants (slow >25); MVMR-Lasso needs more valid IVs than exposures; MVMR-Egger t-df = nIV − n_exp − 1. MR-PRESSO — stops if nIV ≤ n_exp + 2 and requires NbDistribution > nIV. All other thresholds are labelled as definition-based or heuristic in `iv_rule_source`.

### Feasibility rules
| Stage | Rule | Condition | Action | Message shown | Basis |
|---|---|---|---|---|---|
| selection | R01 | data_required involves individual-level data | HIDE (out of scope) | Individual-level analyses are not supported in this app (summary data only). | design decision 2026-10-29 |
| selection | R02 | runnable_tier = no_code | BLOCK | Literature only — no public implementation available. | inventory impl_verification |
| selection | R03 | runnable_tier = reference_only / gui_role = design_or_reference | HIDE from run list | Shown as guidance only. | design decision |
| selection | R04 | gui_role = excluded_by_design | BLOCK | Requires an AI/LLM API key (excluded by project rule). | user requirement |
| selection | R05 | joint_model method and n_selected_exposures < min_exposures or > max_exposures | BLOCK | Method needs {min_exposures}–{max_exposures} exposures in one model; you selected {n}. | inventory |
| selection | R05b | pairwise_loop method with several exposures/outcomes selected | ALLOW (loop) + INFO | Will run separately for each of the {k} exposure–outcome pairs; multiple-testing correction applied. | design decision |
| selection | R06 | joint_model method and n_selected_outcomes < min_outcomes or > max_outcomes | BLOCK | Method needs {min_outcomes}–{max_outcomes} outcomes; you selected {n}. | inventory |
| selection | R07 | n_exposures + n_outcomes (+ mediators) < min_traits_total | BLOCK | Method needs at least {min_traits_total} traits. | inventory |
| selection | R08 | exposure_type = molecular_QTL and selected exposure is not a gene/protein/metabolite QTL dataset | BLOCK | cis-MR methods need a molecular (QTL) exposure with a defined gene region. | method definition |
| selection | R09 | needs_full_sumstats = Y and dataset provides only top hits | BLOCK | Needs full genome-wide summary statistics. | method definition |
| selection | R10 | needs_ld ≠ none and no LD reference for the selected ancestry | BLOCK | No matching LD reference for {ancestry}. | method definition |
| selection | R11 | needs_ldsc_matrices = Y and no LD scores for ancestry | BLOCK | LDSC parameters cannot be estimated for {ancestry}. | method definition |
| selection | R12 | min_ancestries ≥ 2 and datasets exist for only one ancestry | BLOCK | Cross-population methods need datasets from ≥2 ancestries. | method definition |
| selection | R13 | exposure and outcome GWAS ancestries differ and population_rule starts with same_ancestry | WARN (BLOCK if “required”) | Ancestry mismatch biases two-sample MR; consider XMR/TEMR/MR-EILLS or matched datasets. | MR methodology |
| selection | R14 | extra_data not available (e.g. sex-stratified, family, mediator, covariate GWAS) | BLOCK | Needs: {extra_data}. | inventory |
| selection | R15 | outcome_type/exposure_type mismatch (e.g. binary-only method with continuous outcome) | BLOCK | Method supports {outcome_type} outcomes only. | inventory |
| data | R20 | n_IV < min_ivs_hard (evaluated after clumping/harmonisation) | BLOCK | Needs ≥{min_ivs_hard} instruments; found {n_iv}. | iv_rule_source column |
| data | R21 | n_IV < min_ivs_recommended | WARN | Runs, but low power/unstable with {n_iv} instruments. | iv_rule_source column |
| data | R22 | n_IV > max_ivs | BLOCK (offer top-N subset) | Method limited to {max_ivs} instruments. | package code |
| data | R23 | mean F-statistic < 10 (univariable) or conditional F < 10 (MVMR) | WARN + RECOMMEND dIVW / pIVW / MR-RAPS / MRlap / MVMR-cML | Weak instruments detected. | Burgess 2011 (PMID 21414999); Sanderson 2021 (PMID 34338327) |
| data | R24 | I²GX < 0.9 for MR-Egger family | WARN + RECOMMEND SIMEX / Radial-Egger | Egger estimate may be diluted (NOME violation). | Bowden 2016 (PMID 27616674) |
| data | R25 | significant Cochran Q | WARN + RECOMMEND random-effects / robust methods | Heterogeneity suggests pleiotropy. | Bowden 2017 (PMID 28114746) |
| data | R26 | Egger intercept p < 0.05 | WARN | Directional pleiotropy; prefer methods robust to DP. | Bowden 2015 (PMID 26050253) |
| data | R27 | Steiger test indicates reverse direction for many IVs | WARN + offer Steiger filtering / bidirectional methods | Possible reverse causation. | Hemani 2017 (PMID 29149188) |
| data | R28 | estimated sample overlap > 0 (LDSC intercept or shared cohorts) and overlap_handling starts with assumes_none | WARN + RECOMMEND MRlap / MR-APSS / MR-SimSS / CAUSE | Sample overlap can bias toward the observational association. | Burgess 2016 (PMID 27625185) |
| data | R29 | palindromic SNPs without EAF | AUTO-DROP + INFO | Ambiguous palindromic variants removed. | TwoSampleMR harmonise_data |
| data | R30 | binary exposure | INFO | Estimate is per unit log-odds of exposure; ×0.693 gives per doubling of odds. Interpret via underlying continuous liability. | Burgess & Labrecque 2018 (PMID 30039250) |
| data | R31 | MR-PRESSO: NbDistribution ≤ n_IV | AUTO-FIX | NbDistribution raised above n_IV. | MR-PRESSO code |
| data | R32 | LD matrix not positive-definite / LD–GWAS mismatch | AUTO-FIX (prune / DENTIST QC) or BLOCK | LD reference inconsistent with GWAS. | Chen 2021 (PMID 34880243) |
| runtime | R40 | runtime_class = heavy | INFO + background job | Expected long runtime; runs in background with progress. | design decision |
| runtime | R41 | runnable_tier ∈ {external_binary, python, needs_wrapper} | CHECK dependency installed; else BLOCK with install hint | Dependency missing. | design decision |

Tier summary: r_package: 107, no_code: 78, reference_only: 9, external_binary (vendor website, not GitHub; not re-checked): 4, needs_wrapper: 4, web_service: 3, python: 3, r_package (source repository not located): 1, app_implemented (R re-implementation planned; official binaries differ by OS): 1.

## 4. Module overview (proposed GUI navigation)
| Code | Module | Entries |
|---|---|---|
| A | Core two-sample estimators | 15 |
| B | Pleiotropy-robust: median, mode, penalisation, outliers, mixtures | 27 |
| C | Genome-wide & correlated-pleiotropy models (InSIDE violation, sample structure) | 12 |
| D | Bayesian MR | 9 |
| E | Weak-instrument, winner's-curse, overlap & selection-bias corrections | 10 |
| F | Multivariable MR (MVMR) | 18 |
| G | Multi-outcome, bidirectional & causal-network MR | 18 |
| H | cis-MR, molecular-QTL (TWAS/PWAS) & drug-target MR | 23 |
| I | Colocalisation & fine-mapping (companion analyses) | 8 |
| J | Cross-population & multi-dataset MR | 4 |
| K | Non-linear, stratified & interaction MR | 8 |
| L | One-sample / individual-level MR | 19 |
| M | Mediation, time-varying & life-course MR | 7 |
| N | Collider, index-event & covariate-adjusted-GWAS bias | 9 |
| O | Diagnostics, sensitivity & visualisation | 10 |
| P | Study design, power & reporting | 7 |
| Q | Platforms, automation & AI | 6 |
| | **Total method/analysis entries** | **210** |

## 5. Methods by module
### A. Core two-sample estimators

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Wald ratio <br>*aka ratio estimate* | S/I | 2S | — | R · exp 1–1 · out 1–1 · pair · IV≥1 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_wald_ratio` ✓ | Lawlor 2008, Stat Med · PMID 17886233 · [doi](https://doi.org/10.1002/sim.3034) |
| Inverse-variance weighted (fixed effect) <br>*aka IVW-FE* | S | 2S | — | R · exp 1–1 · out 1–1 · pair · IV≥2 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_ivw_fe; mr_ivw(model="fixed")` ✓ | Burgess 2013, Genet Epidemiol · PMID 24114802 · [doi](https://doi.org/10.1002/gepi.21758) |
| IVW (multiplicative random effects) <br>*aka IVW-R; IVW-MRE* | S | 2S | BP | R · exp 1–1 · out 1–1 · pair · IV≥2 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_ivw_mre; mr_ivw(model="random")` ✓ | Bowden 2017, Stat Med · PMID 28114746 · [doi](https://doi.org/10.1002/sim.7221) |
| IVW / MR-Egger with correlated variants (GLS) <br>*aka correlated-IV IVW* | S | 2S | LD | R · exp 1–1 · out 1–1 · pair · IV≥2 · needs: LD | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_ivw(correl=); mr_egger(correl=)` | Burgess 2015, Stat Med · PMID 26661904 · [doi](https://doi.org/10.1002/sim.6835) |
| Maximum likelihood (two-sample) <br>*aka MaxLik* | S | 2S | MEAS | R · exp 1–1 · out 1–1 · pair · IV≥2 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_two_sample_ml; mr_maxlik` ✓ | Burgess 2013, Genet Epidemiol · PMID 24114802 · [doi](https://doi.org/10.1002/gepi.21758) |
| MR-Egger regression <br>*aka Egger* | S | 2S | DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_egger_regression; mr_egger` ✓ | Bowden 2015, Int J Epidemiol · PMID 26050253 · [doi](https://doi.org/10.1093/ije/dyv080) |
| MR-Egger (bootstrap) | S | 2S | DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_egger_regression_bootstrap` ✓ | Bowden 2015, Int J Epidemiol · PMID 26050253 · [doi](https://doi.org/10.1093/ije/dyv080) |
| Simple median | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_simple_median; mr_median(weighting="simple")` ✓ | Bowden 2016, Genet Epidemiol · PMID 27061298 · [doi](https://doi.org/10.1002/gepi.21965) |
| Weighted median <br>*aka WME* | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_weighted_median; mr_median(weighting="weighted")` ✓ | Bowden 2016, Genet Epidemiol · PMID 27061298 · [doi](https://doi.org/10.1002/gepi.21965) |
| Penalised weighted median <br>*aka PWM* | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_penalised_weighted_median; mr_median(weighting="penalized")` ✓ | Burgess S 2016, arXiv · arXiv:1606.03729 ⚑ |
| Mode-based estimate (simple / weighted; NOME variants) <br>*aka MBE; weighted mode; simple mode* | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_simple_mode; mr_weighted_mode; mr_simple_mode_nome; mr_weighted_mode_nome; mr_mbe` ✓ | Hartwig 2017, Int J Epidemiol · PMID 29040600 · [doi](https://doi.org/10.1093/ije/dyx102) |
| Debiased IVW <br>*aka dIVW* | S | 2S | WIV | R · exp 1–1 · out 1–1 · pair · IV≥2 | mr.divw; MendelianRandomization ([tye27/mr.divw](https://github.com/tye27/mr.divw), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr.divw; mr_divw` ✓ | Ye T 2021, Ann Stat · [doi](https://doi.org/10.1214/20-AOS2027) ⚑ |
| Modified debiased IVW <br>*aka mdIVW* | S | 2S | WIV; BP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 | — | Su 2024, Stat Med · PMID 39453381 · [doi](https://doi.org/10.1002/sim.10245) |
| Penalised IVW <br>*aka pIVW* | S | 2S | WIV; BP | R · exp 1–1 · out 1–1 · pair · IV≥2 | mr.pivw; MendelianRandomization ([siqixu/mr.pivw](https://github.com/siqixu/mr.pivw), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_pivw` ✓ | Xu 2022, Biometrics · PMID 35942938 · [doi](https://doi.org/10.1111/biom.13732) |
| Rerandomized IVW <br>*aka RIVW* | S | 2S | WC; WIV | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 · needs: extra data | — | Ma X 2023, Ann Stat · [doi](https://doi.org/10.1214/22-AOS2247) ⚑ |

### B. Pleiotropy-robust: median, mode, penalisation, outliers, mixtures

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Heterogeneity-penalised model averaging <br>*aka MR-HetPen* | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_hetpen` ✓ | Burgess 2018, Int J Epidemiol · PMID 29846613 · [doi](https://doi.org/10.1093/ije/dyy080) |
| MR-Robust (MM-estimation; penalised weights) <br>*aka robust IVW* | S | 2S | BP; outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_ivw(robust=TRUE, penalized=TRUE)` | Rees 2019, PLoS One · PMID 31545794 · [doi](https://doi.org/10.1371/journal.pone.0222362) |
| MR-Lasso | S | 2S | outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_lasso` ✓ | Rees 2019, PLoS One · PMID 31545794 · [doi](https://doi.org/10.1371/journal.pone.0222362) |
| MR-PRESSO (global, outlier, distortion tests) <br>*aka MRPRESSO* | S | 2S | outliers | R · exp 1–any · out 1–1 · pair · IV≥n_exposures + 3 | MRPRESSO; TwoSampleMR ([rondolab/MR-PRESSO](https://github.com/rondolab/MR-PRESSO), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_presso; run_mr_presso` ✓ | Verbanck 2018, Nat Genet · PMID 29686387 · [doi](https://doi.org/10.1038/s41588-018-0099-7) |
| Radial MR (radial IVW / radial Egger, outlier detection) <br>*aka Radial-MR* | S | 2S | outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 (2 for radial IVW only) | RadialMR; TwoSampleMR ([WSpiller/RadialMR](https://github.com/WSpiller/RadialMR), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `ivw_radial; egger_radial; plot_radial; mr_ivw_radial` ✓ | Bowden 2018, Int J Epidemiol · PMID 29961852 · [doi](https://doi.org/10.1093/ije/dyy101) |
| Contamination mixture <br>*aka MR-ConMix; ConMix* | S | 2S | BP; DP; multiple mechanisms | R · exp 1–1 · out 1–1 · pair · IV≥3 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_conmix` ✓ | Burgess 2020, Nat Commun · PMID 31953392 · [doi](https://doi.org/10.1038/s41467-019-14156-4) |
| MRMix (mixture model) <br>*aka MR-Mix; MIX* | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRMix; TwoSampleMR ([gqi/MRMix](https://github.com/gqi/MRMix), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `MRMix; run_mrmix` ✓ | Qi 2019, Nat Commun · PMID 31028273 · [doi](https://doi.org/10.1038/s41467-019-09432-2) |
| MR-Clust (clustering of variant-specific estimates) | S | 2S | multiple mechanisms | R · exp 1–1 · out 1–1 · pair · IV≥3 | mrclust ([cnfoley/mrclust](https://github.com/cnfoley/mrclust)) `mr_clust_em` ✓ | Foley 2021, Bioinformatics · PMID 32915962 · [doi](https://doi.org/10.1093/bioinformatics/btaa778) |
| MR-RAPS (robust adjusted profile score) <br>*aka RAPS* | S | 2S | WIV; BP; outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 | mr.raps; TwoSampleMR ([qingyuanzhao/mr.raps](https://github.com/qingyuanzhao/mr.raps), [MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr.raps; mr_raps` ✓ | Zhao Q 2020, Ann Stat · [doi](https://doi.org/10.1214/19-AOS1866) ⚑ |
| GRAPPLE | S | 3S | WIV; BP; multiple mechanisms | R · exp 1–any · out 1–1 · pair · IV≥3 · needs: extra data | GRAPPLE ([jingshuw/GRAPPLE](https://github.com/jingshuw/GRAPPLE)) `grappleRobustEst` ✓ | Wang 2021, PLoS Genet · PMID 34157017 · [doi](https://doi.org/10.1371/journal.pgen.1009575) |
| MR-cML (constrained maximum likelihood; cML-DP) <br>*aka MRcML; cML* | S | 2S | BP; DP; CP | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRcML; MendelianRandomization ([xue-hr/MRcML](https://github.com/xue-hr/MRcML), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_cML; mr_cML_DP` ✓ | Xue 2021, Am J Hum Genet · PMID 34214446 · [doi](https://doi.org/10.1016/j.ajhg.2021.05.014) |
| MIXIE (mixture of IVW & Egger regressions) | S | 2S | BP; DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | MIXIE ([ZhaotongL/MIXIE](https://github.com/ZhaotongL/MIXIE)) | Lin 2021, PLoS Genet · PMID 34793444 · [doi](https://doi.org/10.1371/journal.pgen.1009922) |
| IMRP (iterative MR & pleiotropy) | S | 2S | outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 | IMRP ([XiaofengZhuCase/IMRP](https://github.com/XiaofengZhuCase/IMRP)) | Zhu 2021, Bioinformatics · PMID 33226062 · [doi](https://doi.org/10.1093/bioinformatics/btaa985) |
| MRCIP | S | 2S | CP; idiosyncratic pleiotropy | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRCIP ([siqixu/MRCIP](https://github.com/siqixu/MRCIP)) | Xu 2021, Brief Bioinform · PMID 33704372 · [doi](https://doi.org/10.1093/bib/bbab019) |
| Constrained instrumental variables <br>*aka CIV* | S/I | 2S | DP (measured) | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 · needs: extra data | — | Jiang 2019, Genet Epidemiol · PMID 30635941 · [doi](https://doi.org/10.1002/gepi.22184) |
| MR-PATH (latent mixture of mechanisms) <br>*aka Path* | S | 2S | multiple mechanisms | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRPATH ([danieliong/MRPATH](https://github.com/danieliong/MRPATH)) | Iong D 2024, Ann Appl Stat · AOAS 18:966-990 ⚑ |
| Local-distribution MR <br>*aka Local* | S | 2S | invalid IVs | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Xu Z 2026, Stat Sin · Stat Sin 36:1-22 ⚑ |
| Welch-weighted Egger regression <br>*aka WWER* | S | 2S | CP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Brown 2021, Am J Hum Genet · PMID 34861175 · [doi](https://doi.org/10.1016/j.ajhg.2021.10.006) |
| MR-GRIP (allele-coding-invariant Egger) <br>*aka GRIPS* | S | 2S | DP | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_grip` ✓ | Dudbridge 2025, PLoS Genet · PMID 41468517 · [doi](https://doi.org/10.1371/journal.pgen.1011967) |
| MR-LOVA (latent outcome variable, stochastic EM) <br>*aka MRLOVA; MRLOVE* | S | 2S | BP; DP; CP (InSIDE test) | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRLOVA ([lamessad/MRLOVA](https://github.com/lamessad/MRLOVA)) | Amente 2025, Hum Genet · PMID 40214754 · [doi](https://doi.org/10.1007/s00439-025-02739-9) |
| dmIVW (debiased mixture IVW) / dmIVW-MA | S | 2S | WIV; BP; CP; MEAS | R · exp 1–1 · out 1–1 · pair · IV≥3 | MRdmIVW ([denglinsui/MRdmIVW](https://github.com/denglinsui/MRdmIVW)) `mr_dmIVW; plot_mr_dmIVW` ✓ | Deng 2026, Brief Bioinform · PMID 42704268 · [doi](https://doi.org/10.1093/bib/bbag468) |
| CMR (conditional-analysis two-stage MR) | S | 2S | LD-induced pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | — | Lin 2024, HGG Adv · PMID 39097774 · [doi](https://doi.org/10.1016/j.xhgg.2024.100339) |
| MR-MoE (mixture of experts method selection) | S | 2S | method selection | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_moe` ✓ | Hemani G 2017, bioRxiv · bioRxiv 173682 ⚑ |
| MR-TRYX (exploiting outliers to find pathways) <br>*aka TRYX* | S | 2S | outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: extra data | tryx ([explodecomputer/tryx](https://github.com/explodecomputer/tryx)) | Cho 2020, Nat Commun · PMID 32081875 · [doi](https://doi.org/10.1038/s41467-020-14452-4) |
| MR-SPI (selection of valid instruments + post-selection inference) | S | 2S | invalid IVs; WIV | R · exp 1–1 · out 1–1 · pair · IV≥3 | MR-SPI ([MinhaoYaooo/MR-SPI](https://github.com/MinhaoYaooo/MR-SPI)) | Yao 2024, Cell Genom · PMID 39637861 · [doi](https://doi.org/10.1016/j.xgen.2024.100700) |
| sisVIVE (some invalid, some valid IV estimator) | I | 1S | invalid IVs | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: indiv | sisVIVE ([cran/sisVIVE](https://github.com/cran/sisVIVE)) `sisVIVE` | Kang H 2016, J Am Stat Assoc · arXiv:1401.5755 ⚑ |
| TSHT (two-stage hard thresholding) / CIIV-type voting | S/I | 1S/2S | invalid IVs | R · exp 1–1 · out 1–1 · pair · IV≥3 | RobustIV ([zijguo/RobustIV](https://github.com/zijguo/RobustIV)) | Guo Z 2018, J R Stat Soc B · JRSSB 80(4):793-815; arXiv:1603.05224 ⚑ |

### C. Genome-wide & correlated-pleiotropy models (InSIDE violation, sample structure)

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| CAUSE | S | 2S | CP; UHP | R · exp 1–1 · out 1–1 · pair · IV≥n/a (uses genome-wide variants) · needs: full GWAS, LD | cause ([jean997/cause](https://github.com/jean997/cause)) `cause` ✓ | Morrison 2020, Nat Genet · PMID 32451458 · [doi](https://doi.org/10.1038/s41588-020-0631-4) |
| MR-APSS <br>*aka MRAPSS* | S | 2S | CP; OVL; POP; WC | R · exp 1–1 · out 1–1 · pair · IV≥n/a (uses genome-wide variants) · needs: full GWAS, LD, LDSC | MRAPSS ([YangLabHKUST/MR-APSS](https://github.com/YangLabHKUST/MR-APSS), [hxh0928/MR-APSS](https://github.com/hxh0928/MR-APSS)) `MRAPSS` ✓ | Hu 2022, Proc Natl Acad Sci U S A · PMID 35787050 · [doi](https://doi.org/10.1073/pnas.2106858119) |
| LHC-MR (latent heritable confounder MR) <br>*aka lhcMR* | S | 2S | CP; bidirectional | R · exp 1–1 · out 1–1 · pair · IV≥n/a (uses genome-wide variants) · needs: full GWAS, LD | lhcMR ([LizaDarrous/lhcMR](https://github.com/LizaDarrous/lhcMR)) `lhc_mr` ✓ | Darrous 2021, Nat Commun · PMID 34907193 · [doi](https://doi.org/10.1038/s41467-021-26970-w) |
| LCV (latent causal variable; genetic causality proportion) <br>*aka GCP* | S | 2S | CP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥n/a (uses genome-wide variants) · needs: full GWAS, LD | — | O'Connor 2018, Nat Genet · PMID 30374074 · [doi](https://doi.org/10.1038/s41588-018-0255-0) |
| MR-Corr2 | S | 2S | CP; LD | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | MR.Corr2 ([QingCheng0218/MR.Corr2](https://github.com/QingCheng0218/MR.Corr2)) | Cheng 2022, Bioinformatics · PMID 34499127 · [doi](https://doi.org/10.1093/bioinformatics/btab646) |
| MR-LDP | S | 2S | LD; pleiotropy | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | MR.LDP ([QingCheng0218/MR.LDP](https://github.com/QingCheng0218/MR.LDP)) | Cheng 2020, NAR Genom Bioinform · PMID 33575584 · [doi](https://doi.org/10.1093/nargab/lqaa028) |
| MR-CUE | S | 2S | CP | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | MR.CUE ([QingCheng0218/MR.CUE](https://github.com/QingCheng0218/MR.CUE)) | Cheng 2022, Nat Commun · PMID 36310177 · [doi](https://doi.org/10.1038/s41467-022-34164-1) |
| RMR-ICP | S | 2S | CP; idiosyncratic pleiotropy | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | RMR.ICP ([QingCheng0218/RMR.ICP](https://github.com/QingCheng0218/RMR.ICP)) | Cheng 2025, Brief Bioinform · PMID 41016013 · [doi](https://doi.org/10.1093/bib/bbaf508) |
| MRAID | S | 2S | CP; UHP; LD | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | MRAID ([yuanzhongshang/MRAID](https://github.com/yuanzhongshang/MRAID)) `MRAID` ✓ | Yuan 2022, Sci Adv · PMID 35235357 · [doi](https://doi.org/10.1126/sciadv.abl5744) |
| GSMR (with HEIDI-outlier) <br>*aka GCTA-GSMR* | S | 2S | outliers; LD | binary · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | GCTA (gsmr) `gcta --gsmr-file` | Zhu 2018, Nat Commun · PMID 29335400 · [doi](https://doi.org/10.1038/s41467-017-02317-2) |
| MR-MU (multiple unmeasured confounders) <br>*aka MRMU* | S | 2S | measured + unmeasured confounding; OVL | R · exp 1–any (exposure + measured confounders) · out 1–1 · pair · IV≥3 · needs: full GWAS, LD, LDSC | MRMU ([hxh0928/MRMU](https://github.com/hxh0928/MRMU)) `MRMU; MRMU_Input; MV_harmonise` ✓ | Hu X n/a, R package only · github.com/hxh0928/MRMU ⚑ |
| MRCI (reciprocal causation mixture model) | S | 2S | CP; bidirectional | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥n/a (uses genome-wide variants) · needs: full GWAS, LD, extra data | — | Liu 2023, Nat Commun · PMID 36854672 · [doi](https://doi.org/10.1038/s41467-023-36490-4) |

### D. Bayesian MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| BWMR (Bayesian weighted MR) | S | 2S | WIV; outliers | R · exp 1–1 · out 1–1 · pair · IV≥3 | BWMR ([jiazhao97/BWMR](https://github.com/jiazhao97/BWMR)) `BWMR` ✓ | Zhao 2020, Bioinformatics · PMID 31593215 · [doi](https://doi.org/10.1093/bioinformatics/btz749) |
| MR-Horse (horseshoe prior) <br>*aka MR-HORSE* | S | 2S | BP; DP; CP | R-script · exp 1–1 · out 1–1 · pair · IV≥3 | mr_horse.R (sourced script; R2jags/JAGS) ([AJ-Grant/mrhorse](https://github.com/AJ-Grant/mrhorse)) `mr_horse` ✓ | Grant 2024, Am J Hum Genet · PMID 38181732 · [doi](https://doi.org/10.1016/j.ajhg.2023.12.002) |
| MR-BMA (Bayesian model averaging, risk-factor selection) <br>*aka MVBMA* | S | 2S | correlated exposures | R-script · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 | demo_AMD scripts ([verena-zuber/demo_AMD](https://github.com/verena-zuber/demo_AMD)) | Zuber 2020, Nat Commun · PMID 31911605 · [doi](https://doi.org/10.1038/s41467-019-13870-3) |
| JAM-MR (pleiotropic loss Bayesian variable selection) | S | 2S | pleiotropy | R · exp 1–1 · out 1–1 · pair · IV≥3 | R2BGLiMS ([pjnewcombe/R2BGLiMS](https://github.com/pjnewcombe/R2BGLiMS)) `JAMMR` ✓ | Gkatzionis 2021, Stat Med · PMID 34155684 · [doi](https://doi.org/10.1002/sim.9109) |
| Bayesian MR-Egger with weakly informative priors <br>*aka BMRE* | S | 2S | DP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Schmidt 2018, Int J Epidemiol · PMID 29253155 · [doi](https://doi.org/10.1093/ije/dyx254) |
| Bayesian MR with multiple pleiotropic variants (horseshoe) | S/I | 2S | pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Berzuini 2020, Biostatistics · PMID 30084873 · [doi](https://doi.org/10.1093/biostatistics/kxy027) |
| Bayesian model averaging over pleiotropy models | S | 2S | BP; DP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Thompson 2017, Stat Med · PMID 28850703 · [doi](https://doi.org/10.1002/sim.7442) |
| Profile-likelihood BMA (MR-RAPS-based) | S | 2S | WIV; pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Shapland 2022, Stat Med · PMID 35060160 · [doi](https://doi.org/10.1002/sim.9320) |
| RBMR (robust Bayesian MR) | S | 2S | LD; idiosyncratic pleiotropy | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: LD | RBMR ([AnqiWang2021/RBMR](https://github.com/AnqiWang2021/RBMR)) | Wang 2022, Genet Epidemiol · PMID 35192729 · [doi](https://doi.org/10.1002/gepi.22445) |

### E. Weak-instrument, winner's-curse, overlap & selection-bias corrections

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| MRlap | S | 2S/overlap | OVL; WIV; WC | R · exp 1–1 · out 1–1 · pair · IV≥n/a (selects IVs internally) · needs: full GWAS, LD | MRlap ([n-mounier/MRlap](https://github.com/n-mounier/MRlap)) `MRlap` | Mounier 2023, Genet Epidemiol · PMID 37036286 · [doi](https://doi.org/10.1002/gepi.22522) |
| MR-SimSS (simulated sample splitting) | S | 2S/overlap | WC; OVL | R · exp 1–1 · out 1–1 · pair · IV≥inherits wrapped method · needs: full GWAS | mr.simss ([amandaforde/mr.simss](https://github.com/amandaforde/mr.simss)) | Forde 2026, PLoS Genet · PMID 42102131 · [doi](https://doi.org/10.1371/journal.pgen.1011949) |
| Participant-overlap bias assessment | S | overlap | OVL | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: extra data | — | Burgess 2016, Genet Epidemiol · PMID 27625185 · [doi](https://doi.org/10.1002/gepi.21998) |
| Weak-instrument robust tests (two-sample) | S | 2S | WIV | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 | — | Wang 2021, Biometrics · PMID 34213007 · [doi](https://doi.org/10.1111/biom.13524) |
| Conditional likelihood-ratio test & PC-GMM (weak cis instruments) <br>*aka mr_clr; mr_pcgmm* | S | 2S | WIV; LD | R · exp 1–any (mr_mvpcgmm) · out 1–1 · pair · IV≥2 · needs: LD | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_clr; mr_pcgmm; mr_mvpcgmm` ✓ | Patel 2023, Biometrics · PMID 37337418 · [doi](https://doi.org/10.1111/biom.13888) |
| Semi-empirical Bayes calibration for sex-specific MR | S | 2S | WIV | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 · needs: extra data | — | Huang 2025, Am J Hum Genet · PMID 40912238 · [doi](https://doi.org/10.1016/j.ajhg.2025.07.015) |
| Sex-stratified weak- & pleiotropy-robust MR | S/I | 1S/2S | WIV; pleiotropy cancellation | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 · needs: extra data | — | Karageorgiou 2023, Genet Epidemiol · PMID 36682072 · [doi](https://doi.org/10.1002/gepi.22512) |
| Instruments-for-selection (Heckman-type) adjustment | I | 1S | SEL | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | Gkatzionis 2024, Stat Med · PMID 39039030 · [doi](https://doi.org/10.1002/sim.10173) |
| Negative-control outcomes for selection bias | S/I | 1S/2S | SEL (detection) | ⛔ no code · exp 1–1 · out 1–1 · joint · IV≥— · needs: extra data | — | Gkatzionis 2026, Stat Med · PMID 42592747 · [doi](https://doi.org/10.1002/sim.70703) |
| Valid post-selection inference for two-sample MR | S | 2S | WC (selection) | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 | — | Wang 2026, Stat Med · PMID 42893902 · [doi](https://doi.org/10.1002/sim.70784) |

### F. Multivariable MR (MVMR)

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| MVMR-IVW <br>*aka multivariable IVW* | S | 2S | measured pleiotropy | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 | TwoSampleMR; MendelianRandomization; MVMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization), [WSpiller/MVMR](https://github.com/WSpiller/MVMR)) `mv_multiple; mv_ivw; mr_mvivw; ivw_mvmr` ✓ | Burgess 2015, Am J Epidemiol · PMID 25632051 · [doi](https://doi.org/10.1093/aje/kwu283) |
| MVMR-Egger <br>*aka MRMV-Egger* | S | 2S | DP | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_mvegger` ✓ | Rees 2017, Stat Med · PMID 28960498 · [doi](https://doi.org/10.1002/sim.7492) |
| MVMR-Median <br>*aka MRMV-WME* | S | 2S | BP; DP | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_mvmedian` ✓ | Grant 2021, Stat Med · PMID 34342032 · [doi](https://doi.org/10.1002/sim.9156) |
| MVMR-Lasso <br>*aka MVLasso* | S | 2S | outliers | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 valid IVs after penalisation | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_mvlasso` ✓ | Grant 2021, Stat Med · PMID 34342032 · [doi](https://doi.org/10.1002/sim.9156) |
| MVMR-Robust | S | 2S | outliers | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `(see package docs)` | Grant 2021, Stat Med · PMID 34342032 · [doi](https://doi.org/10.1002/sim.9156) |
| MVMR-cML <br>*aka MVcML* | S | 2S | BP; DP; CP | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | MVMR-cML; MendelianRandomization ([ZhaotongL/MVMR-cML](https://github.com/ZhaotongL/MVMR-cML), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `MVmr_cML; MVmr_cML_DP; mr_mvcML` ✓ | Lin 2023, Am J Hum Genet · PMID 36948188 · [doi](https://doi.org/10.1016/j.ajhg.2023.02.014) |
| MVMR-Horse <br>*aka MVHorse* | S | 2S | BP; DP; CP | R-script · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 | mr_horse.R (sourced script; R2jags/JAGS) ([AJ-Grant/mrhorse](https://github.com/AJ-Grant/mrhorse)) `mvmr_horse` ✓ | Grant 2024, Am J Hum Genet · PMID 38181732 · [doi](https://doi.org/10.1016/j.ajhg.2023.12.002) |
| MVMRmode (multivariable MBE & multivariable contamination mixture) <br>*aka MRMV-WBE; mv_mrmode* | S | 2S | BP; DP | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | MVMRmode `mv_mrmode` ✓ | Woolf 2024, PLoS One · PMID 38713711 · [doi](https://doi.org/10.1371/journal.pone.0291183) |
| Radial MVMR <br>*aka RMVMR* | S | 2S | outliers | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 | RMVMR ([WSpiller/RMVMR](https://github.com/WSpiller/RMVMR)) `format_rmvmr; ivw_rmvmr; strength_rmvmr; plot_rmvmr` ✓ | Spiller 2024, PLoS Genet · PMID 39680583 · [doi](https://doi.org/10.1371/journal.pgen.1011506) |
| Conditional F-statistic & Q_A heterogeneity for MVMR | S | 2S | WIV (conditional) | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: extra data | MVMR ([WSpiller/MVMR](https://github.com/WSpiller/MVMR)) `format_mvmr; strength_mvmr; pleiotropy_mvmr; qhet_mvmr` ✓ | Sanderson 2021, Stat Med · PMID 34338327 · [doi](https://doi.org/10.1002/sim.9133) |
| MRBEE (bias-correction estimating equations) | S | 2S | WIV; OVL; pleiotropy | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: full GWAS | MRBEE ([noahlorinczcomi/MRBEE](https://github.com/noahlorinczcomi/MRBEE)) | Lorincz-Comi 2024, HGG Adv · PMID 38582968 · [doi](https://doi.org/10.1016/j.xhgg.2024.100290) |
| LMM-MEC (linear mixed model with measurement-error correction) <br>*aka TS_LMM* | S | 2S | MEAS; LD | R-script · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: LD | LMM_MEC (R scripts) ([mingding-hsph/LMM_MEC](https://github.com/mingding-hsph/LMM_MEC)) `LMM_MEC.macro.r` | Ding 2026, Genet Epidemiol · PMID 41553776 · [doi](https://doi.org/10.1002/gepi.70032) |
| Robust MVMR (Wu, Kang & Ye) | S | 2S | WIV; pleiotropy | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 | — | Wu Y 2025, Biometrika · Biometrika 112:asaf053 ⚑ |
| Sparse PCA MVMR (highly correlated exposures) <br>*aka sPCA-MVMR* | S | 2S | correlated exposures | ⛔ no code · exp 3–any · out 1–1 · joint · IV≥n_exposures + 1 | — | Karageorgiou 2023, Elife · PMID 37074034 · [doi](https://doi.org/10.7554/eLife.80063) |
| PCA-based MVMR for shared clustered predictors | S | 2S | LD; correlated exposures | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: LD | — | Batool 2022, Genet Epidemiol · PMID 35638254 · [doi](https://doi.org/10.1002/gepi.22462) |
| Robust MVMR for highly correlated exposures (cML + Bayesian clustering) | S | 2S | CP; correlated exposures | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥n_exposures + 2 | — | Chan 2024, Am J Hum Genet · PMID 39106865 · [doi](https://doi.org/10.1016/j.ajhg.2024.07.007) |
| MR-PL (one-sample MVMR via PLS + Lasso) | I | 1S | high-dimensional exposures | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥1 · needs: indiv | — | Yang 2023, PLoS Genet · PMID 38150468 · [doi](https://doi.org/10.1371/journal.pgen.1011112) |
| Regularised selection of measured pleiotropic covariates | S | 2S | measured pleiotropy | ⛔ no code · exp 1–any · out 1–1 · joint · IV≥n_covariates + 2 · needs: extra data | — | Grant 2022, Biostatistics · PMID 33155035 · [doi](https://doi.org/10.1093/biostatistics/kxaa045) |

### G. Multi-outcome, bidirectional & causal-network MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| MR2 (multi-response MR) | S | 2S | shared/distinct exposures | R · exp 1–any · out 2–any · joint · IV≥n_exposures + 1 | MR2 ([lb664/MR2](https://github.com/lb664/MR2)) | Zuber 2023, Am J Hum Genet · PMID 37419091 · [doi](https://doi.org/10.1016/j.ajhg.2023.06.005) |
| MrDAG (Bayesian causal graphical model) | S | 2S | network | R · exp 1–any · out 2–any · joint · IV≥n_exposures + 1 | MrDAG ([lb664/MrDAG](https://github.com/lb664/MrDAG)) | Zuber 2025, Am J Hum Genet · PMID 40179887 · [doi](https://doi.org/10.1016/j.ajhg.2025.03.005) |
| METEOR (multi-outcome, data-adaptive) | S | 2S | OVL; CP | R · exp 1–any · out 2–any · joint · IV≥n_exposures + 1 · needs: full GWAS, LD, LDSC | METEOR ([yuanzhongshang/METEOR](https://github.com/yuanzhongshang/METEOR)) | Zhang 2026, Brief Bioinform · PMID 42407119 · [doi](https://doi.org/10.1093/bib/bbag364) |
| MRMO (two-stage multivariate MR on mixed outcomes) | I | 1S | — | ⛔ no code · exp 1–1 · out 2–any · joint · IV≥1 · needs: indiv | — | Deng 2023, Stat Methods Med Res · PMID 37338962 · [doi](https://doi.org/10.1177/09622802231181220) |
| BMRMO (Bayesian two-stage multivariate MR) | I | 1S | — | ⛔ no code · exp 1–1 · out 2–any · joint · IV≥1 · needs: indiv | — | Deng 2023, Stat Med · PMID 36998123 · [doi](https://doi.org/10.1002/sim.9721) |
| MR-AHC (agglomerative hierarchical clustering) <br>*aka MRAHC* | S | 2S | multiple mechanisms | R · exp 1–1 · out 2–any · joint · IV≥3 | MRAHC ([xiaoran-liang/MRAHC](https://github.com/xiaoran-liang/MRAHC)) | Liang 2024, Genet Epidemiol · PMID 39138631 · [doi](https://doi.org/10.1002/gepi.22582) |
| FusioMR (single- and multi-outcome; molecular exposures) | S | 2S | pleiotropy; few IVs | ⛔ no code · exp 1–1 · out 1–any · pair · IV≥1 | — | Kang 2026, Am J Hum Genet · PMID 41997142 · [doi](https://doi.org/10.1016/j.ajhg.2026.03.017) |
| Bidirectional MR-cML (BiDirectCausal) | S | 2S | CP; reverse causation | R · exp 1–1 · out 1–1 · pair · IV≥3 per direction · needs: extra data | BiDirectCausal ([xue-hr/BiDirectCausal](https://github.com/xue-hr/BiDirectCausal)) `BiDirMRcML; BiDirCDcML; BiDirCDMethod` ✓ | Xue 2022, PLoS Genet · PMID 35576237 · [doi](https://doi.org/10.1371/journal.pgen.1010205) |
| BayBiMR (Bayesian bidirectional MR) | S | 2S | CP; UHP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 per direction · needs: extra data | — | Chen 2026, Stat Med · PMID 42433205 · [doi](https://doi.org/10.1002/sim.70667) |
| Focusing framework for bi-directional causal tests | S | 2S | invalid IVs; reverse causation | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 per direction · needs: extra data | — | Li S 2024, J R Stat Soc B · JRSSB 87(2):529-548 ⚑ |
| MR2G (cyclic causal-network inference) | S | 2S | network; CP | R · exp 1–any · out 1–any · joint · IV≥3 per trait | MR2G ([ZhaotongL/MR2G](https://github.com/ZhaotongL/MR2G)) | Lin 2026, PLoS Genet · PMID 42189912 · [doi](https://doi.org/10.1371/journal.pgen.1012144) |
| Graph-MRcML (MR + network deconvolution) <br>*aka cMLgraph* | S | 2S | network; CP | web · exp 1–any · out 1–any · joint · IV≥3 per trait | Shiny app ([link](https://zhaotongl.shinyapps.io/cMLgraph/)) | Lin 2023, PLoS Genet · PMID 37200398 · [doi](https://doi.org/10.1371/journal.pgen.1010762) |
| DAG of phenotypes from GWAS summary statistics <br>*aka sumdag* | S | 2S | network | R · exp 1–any · out 1–any · joint · IV≥3 per trait | sumdag ([chunlinli/sumdag](https://github.com/chunlinli/sumdag)) | Zilinskas 2024, Biometrics · PMID 38470257 · [doi](https://doi.org/10.1093/biomtc/ujad039) |
| MRPC (causal networks with the MR principle) | I | 1S | network | R · exp 1–any · out 1–any · joint · IV≥1 · needs: indiv | MRPC ([cran/MRPC](https://github.com/cran/MRPC)) `MRPC` | Badsha 2019, Front Genet · PMID 31164902 · [doi](https://doi.org/10.3389/fgene.2019.00460) |
| cGAUGE (phenome-wide causal discovery) | S | 2S | network | ⛔ no code · exp 1–any · out 1–any · joint · IV≥3 per trait | — | Amar 2021, Nat Commun · PMID 33441555 · [doi](https://doi.org/10.1038/s41467-020-20516-2) |
| MR.RGM (Bayesian multivariate bidirectional MR networks) | S/I | 1S/2S | network; bidirectional | R · exp 1–any · out 1–any · joint · IV≥3 per trait | MR.RGM ([cran/MR.RGM](https://github.com/cran/MR.RGM), [bitansa/MR.RGM](https://github.com/bitansa/MR.RGM)) | Sarkar 2025, Bioinformatics · PMID 40128153 · [doi](https://doi.org/10.1093/bioinformatics/btaf130) |
| Bayesian network analysis with genetic anchors <br>*aka BN* | I | 1S | pleiotropy | ⛔ no code · exp 1–any · out 1–any · joint · IV≥1 · needs: indiv | — | Howey 2020, PLoS Genet · PMID 32119656 · [doi](https://doi.org/10.1371/journal.pgen.1008198) |
| Network MR (two-step) <br>*aka two-step MR* | S/I | 2S | mediation | ⛔ no code · exp 1–1 · out 1–1 · joint · IV≥2 per step · needs: extra data | — | Burgess 2014, Int J Epidemiol · PMID 25150977 · [doi](https://doi.org/10.1093/ije/dyu176) |

### H. cis-MR, molecular-QTL (TWAS/PWAS) & drug-target MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| SMR + HEIDI | S | 2S | linkage (HEIDI) | R · exp 1–1 · out 1–1 · pair · IV≥1 (top cis-QTL); HEIDI needs ≥3 variants · needs: LD, QTL exposure | SMR software `smr` | Zhu 2016, Nat Genet · PMID 27019110 · [doi](https://doi.org/10.1038/ng.3538) |
| OPERA (multi-omics pleiotropic association) | S | 2S | — | binary · exp 1–any (multiple omics layers) · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | OPERA software | Wu 2023, Cell Genom · PMID 37601976 · [doi](https://doi.org/10.1016/j.xgen.2023.100344) |
| TWMR (transcriptome-wide MR) | S | 2S | LD; pleiotropy (multi-gene) | R · exp 1–any (genes) · out 1–1 · pair · IV≥n_genes + 1 · needs: LD, QTL exposure | TWMR ([eleporcu/TWMR](https://github.com/eleporcu/TWMR)) | Porcu 2019, Nat Commun · PMID 31341166 · [doi](https://doi.org/10.1038/s41467-019-10936-0) |
| revTWMR (reverse TWMR) | S | 2S | reverse causation | R · exp 1–1 · out 1–any · pair · IV≥1 | TWMR ([eleporcu/TWMR](https://github.com/eleporcu/TWMR)) | Porcu 2021, Nat Commun · PMID 34561431 · [doi](https://doi.org/10.1038/s41467-021-25805-y) |
| MR-JTI (joint-tissue imputation + MR) | S/I | 2S | pleiotropy; LD | R · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure | MR-JTI ([gamazonlab/MR-JTI](https://github.com/gamazonlab/MR-JTI)) | Zhou 2020, Nat Genet · PMID 33020666 · [doi](https://doi.org/10.1038/s41588-020-0706-2) |
| PMR-Egger (probabilistic MR for TWAS) <br>*aka PMR* | S/I | 2S | DP; LD | R · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure | PMR ([yuanzhongshang/PMR](https://github.com/yuanzhongshang/PMR)) | Yuan 2020, Nat Commun · PMID 32737316 · [doi](https://doi.org/10.1038/s41467-020-17668-6) |
| MRLocus | S | 2S | allelic heterogeneity | R · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure | mrlocus ([mikelove/mrlocus](https://github.com/mikelove/mrlocus)) | Zhu 2021, PLoS Genet · PMID 33872308 · [doi](https://doi.org/10.1371/journal.pgen.1009455) |
| MR-link | I+S | 2S | pleiotropy; LD | Py · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | genome_integration ([adriaan-vd-graaf/genome_integration](https://github.com/adriaan-vd-graaf/genome_integration)) | van der Graaf 2020, Nat Commun · PMID 33004804 · [doi](https://doi.org/10.1038/s41467-020-18716-x) |
| MR-link-2 | S | 2S | pleiotropy; LD | Py · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure | mrlink2 ([adriaan-vd-graaf/mrlink2](https://github.com/adriaan-vd-graaf/mrlink2)) | van der Graaf 2025, Nat Commun · PMID 40610416 · [doi](https://doi.org/10.1038/s41467-025-60868-1) |
| cis-MRcML | S | 2S | CP; LD | R · exp 1–1 · out 1–1 · pair · IV≥2 · needs: LD, QTL exposure | cisMRcML ([ZhaotongL/cisMRcML](https://github.com/ZhaotongL/cisMRcML)) `cismr_cML_DP` ✓ | Lin 2024, Nat Commun · PMID 39025905 · [doi](https://doi.org/10.1038/s41467-024-50385-y) |
| PCA-IVW for fine-mapped correlated variants | S | 2S | LD | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 · needs: LD, QTL exposure | — | Burgess 2017, Genet Epidemiol · PMID 28944551 · [doi](https://doi.org/10.1002/gepi.22077) |
| Summary-data cis-MR methods (comparison & guidance) | S | 2S | LD; WIV | ref · exp 1–1 · out 1–1 · pair · IV≥2 · needs: LD, QTL exposure | — | Gkatzionis 2022, Genet Epidemiol · PMID 36273411 · [doi](https://doi.org/10.1002/gepi.22506) |
| cis-multivariable MR (dimension-reduced, overdispersion) | S | 2S | LD; heterogeneity | ⛔ no code · exp 1–any · out 1–1 · pair · IV≥n_exposures + 1 · needs: LD, QTL exposure | — | Patel 2024, Genet Epidemiol · PMID 38379245 · [doi](https://doi.org/10.1002/gepi.22551) |
| cis-MRBEE | S | 2S | WIV; pleiotropy | ⛔ no code · exp 1–any · out 1–1 · pair · IV≥n_exposures + 1 · needs: LD, QTL exposure | — | Yang 2025, Brief Bioinform · PMID 40471990 · [doi](https://doi.org/10.1093/bib/bbaf250) |
| EMIC (effective-median MR for causal genes) | S | 2S | pleiotropy | binary · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure | KGGSEE | Jiang 2022, Am J Hum Genet · PMID 35460606 · [doi](https://doi.org/10.1016/j.ajhg.2022.04.004) |
| MR-MtRobin (multi-tissue TWMR with invalid IVs) | S | 2S | invalid IVs | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | — | Gleason 2021, Genet Epidemiol · PMID 33834509 · [doi](https://doi.org/10.1002/gepi.22380) |
| mintMR (multi-context integrative MR) | S | 2S | tissue sharing | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | — | Lu 2024, Am J Hum Genet · PMID 39053459 · [doi](https://doi.org/10.1016/j.ajhg.2024.06.012) |
| MULTI (multi-tissue unified likelihood MR) | S | 2S | tissue sharing | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | — | Cheng 2026, Brief Bioinform · PMID 42520158 · [doi](https://doi.org/10.1093/bib/bbag414) |
| MT-TWMR (multi-tissue TWMR) | S | 2S | tissue sharing | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | — | Wang 2026, BMC Bioinformatics · PMID 42157043 · [doi](https://doi.org/10.1186/s12859-026-06451-x) |
| Drug-target MR framework | S | 2S | — | ⛔ no code · exp 1–1 · out 1–1 ·  · IV≥1 · needs: LD, QTL exposure | — | Schmidt 2020, Nat Commun · PMID 32591531 · [doi](https://doi.org/10.1038/s41467-020-16969-0) |
| MR-CARV (common + annotation-informed rare variants) | S | 2S | — | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: LD, QTL exposure, extra data | — | Zhang 2026, Brief Bioinform · PMID 41499219 · [doi](https://doi.org/10.1093/bib/bbaf649) |
| RARE (rare-variant multivariable MR) | S | 2S | — | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥1 · needs: LD, QTL exposure, extra data | — | Cheng 2025, Brief Bioinform · PMID 40370099 · [doi](https://doi.org/10.1093/bib/bbaf214) |
| hJAM (hierarchical joint analysis of marginal statistics) | S | 2S | LD | ⛔ no code · exp 1–any · out 1–1 · pair · IV≥n_exposures + 1 · needs: LD, QTL exposure | — | Jiang 2021, Am J Epidemiol · PMID 33404048 · [doi](https://doi.org/10.1093/aje/kwaa287) |

### I. Colocalisation & fine-mapping (companion analyses)

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| coloc (ABF colocalisation) | S | — | shared causal variant | R · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, extra data | coloc ([chr1swallace/coloc](https://github.com/chr1swallace/coloc)) `coloc.abf` ✓ | Giambartolomei 2014, PLoS Genet · PMID 24830394 · [doi](https://doi.org/10.1371/journal.pgen.1004383) |
| coloc-SuSiE (multiple causal variants) | S | — | multiple causal variants | R · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, LD | coloc; susieR ([chr1swallace/coloc](https://github.com/chr1swallace/coloc), [stephenslab/susieR](https://github.com/stephenslab/susieR)) `coloc.susie` ✓ | Wallace 2021, PLoS Genet · PMID 34587156 · [doi](https://doi.org/10.1371/journal.pgen.1009440) |
| HyPrColoc (multi-trait colocalisation) | S | — | many traits | R · exp 1–1 · out 1–any · pair · IV≥— · needs: full GWAS | hyprcoloc ([cnfoley/hyprcoloc](https://github.com/cnfoley/hyprcoloc)) `hyprcoloc` ✓ | Foley 2021, Nat Commun · PMID 33536417 · [doi](https://doi.org/10.1038/s41467-020-20885-8) |
| eCAVIAR | S | — | multiple causal variants | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, LD | — | Hormozdiari 2016, Am J Hum Genet · PMID 27866706 · [doi](https://doi.org/10.1016/j.ajhg.2016.10.003) |
| moloc (multiple-trait colocalisation) | S | — | ≥3 traits | R · exp 1–1 · out 1–1 · joint · IV≥— · needs: full GWAS | moloc ([clagiamba/moloc](https://github.com/clagiamba/moloc)) | Giambartolomei 2018, Bioinformatics · PMID 29579179 · [doi](https://doi.org/10.1093/bioinformatics/bty147) |
| SuSiE (fine-mapping) | S/I | — | — | R · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, LD | susieR ([stephenslab/susieR](https://github.com/stephenslab/susieR)) `susie; susie_rss` ✓ | Wang 2020, J R Stat Soc Series B Stat Methodol · PMID 37220626 · [doi](https://doi.org/10.1111/rssb.12388) |
| FINEMAP | S | — | — | binary · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, LD | FINEMAP | Benner 2016, Bioinformatics · PMID 26773131 · [doi](https://doi.org/10.1093/bioinformatics/btw018) |
| XMAP (cross-population fine-mapping) | S | — | confounding bias | R · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, LD, ≥2 ancestries | XMAP ([YangLabHKUST/XMAP](https://github.com/YangLabHKUST/XMAP)) | Cai 2023, Nat Commun · PMID 37898663 · [doi](https://doi.org/10.1038/s41467-023-42614-7) |

### J. Cross-population & multi-dataset MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| XMR (cross-population MR) | S | 2S (multi-ancestry) | HET; OVL; WIV | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: full GWAS, LD, LDSC, ≥2 ancestries | XMR ([YangLabHKUST/XMR](https://github.com/YangLabHKUST/XMR)) `fit_XMR` ✓ | Huang X 2026, medRxiv · [doi](https://doi.org/10.64898/2026.03.10.26348003) ⚑ |
| TEMR (trans-ethnic MR) | S | 2S (multi-ancestry) | HET | R · exp 1–1 · out 1–1 · pair · IV≥3 · needs: ≥2 ancestries | TEMR ([hhoulei/TEMR](https://github.com/hhoulei/TEMR)) | Hou 2024, Am J Hum Genet · PMID 39689714 · [doi](https://doi.org/10.1016/j.ajhg.2024.11.006) |
| MR-EILLS (invariance-based, multiple heterogeneous datasets) <br>*aka EILLS* | S | 2S (multi-dataset) | HET; CP; UHP | R · exp 1–any · out 1–1 · pair · IV≥3 · needs: ≥2 ancestries | MREILLS ([hhoulei/MREILLS](https://github.com/hhoulei/MREILLS)) `MREILLS; BOOT_MREILLS; CHO_lambda` ✓ | Hou 2025, Nat Commun · PMID 40826001 · [doi](https://doi.org/10.1038/s41467-025-62823-6) |
| TLMR (transfer-learning MR) | I | 1S | HET | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, ≥2 ancestries | — | Wei 2026, Brief Bioinform · PMID 41662355 · [doi](https://doi.org/10.1093/bib/bbag019) |

### K. Non-linear, stratified & interaction MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Non-linear MR: residual method (fractional polynomial, piecewise linear) <br>*aka nlmr; LACE* | I | 1S | non-linearity | R · exp 1–1 · out 1–1 · pair · IV≥1 (score) · needs: indiv, extra data | nlmr ([jrs95/nlmr](https://github.com/jrs95/nlmr)) `fracpoly_mr; piecewise_mr` ✓ | Staley 2017, Genet Epidemiol · PMID 28317167 · [doi](https://doi.org/10.1002/gepi.22041) |
| SUMnlmr (summarised-data non-linear MR) | I→S | 1S | non-linearity | R · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | SUMnlmr ([amymariemason/SUMnlmr](https://github.com/amymariemason/SUMnlmr)) `create_summary_data; frac_poly_summ_mr; piecewise_summ_mr` ✓ | Mason AM 2022, Int J Epidemiol · PMCID PMC9749704 ⚑ |
| Doubly-ranked stratification <br>*aka DRMR* | I | 1S | non-linearity | R · exp 1–1 · out 1–1 · pair · IV≥1 (score) · needs: indiv, extra data | DRMR ([HDTian/DRMR](https://github.com/HDTian/DRMR)) `Stratify; getSummaryInf; Smooth` ✓ | Tian 2023, PLoS Genet · PMID 37390109 · [doi](https://doi.org/10.1371/journal.pgen.1010823) |
| PolyMR (polynomial MR) | I | 1S | non-linearity | R · exp 1–1 · out 1–1 · pair · IV≥1 (score) · needs: indiv, extra data | PolyMR ([JonSulc/PolyMR](https://github.com/JonSulc/PolyMR)) `polymr; plot_polymr` ✓ | Sulc 2022, HGG Adv · PMID 35832928 · [doi](https://doi.org/10.1016/j.xhgg.2022.100124) |
| Quantile IV (machine-learning MR estimator) | I | 1S | non-linearity | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 (score) · needs: indiv, extra data | — | Legault 2025, Am J Hum Genet · PMID 40378846 · [doi](https://doi.org/10.1016/j.ajhg.2025.04.010) |
| Stratified MR framework for effect modifiers | I | 1S | effect modification | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 (score) · needs: indiv, extra data | — | Man 2025, Eur J Epidemiol · PMID 40072671 · [doi](https://doi.org/10.1007/s10654-025-01213-0) |
| Factorial MR (interactions) | I | 1S | interaction | ⛔ no code · exp 2–2 · out 1–1 · joint · IV≥1 · needs: indiv, extra data | — | Rees 2020, Int J Epidemiol · PMID 31369124 · [doi](https://doi.org/10.1093/ije/dyz161) |
| MR-GxE (gene-by-environment interaction) | I | 1S | pleiotropy correction | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | Spiller 2019, Int J Epidemiol · PMID 30462199 · [doi](https://doi.org/10.1093/ije/dyy204) |

### L. One-sample / individual-level MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Two-stage least squares <br>*aka 2SLS; TSLS* | I | 1S | — | R · exp 1–any · out 1–1 · pair · IV≥n_exposures (just-identified) · needs: indiv | ivreg; AER; OneSampleMR ([cran/ivreg](https://github.com/cran/ivreg), [cran/AER](https://github.com/cran/AER), [remlapmot/OneSampleMR](https://github.com/remlapmot/OneSampleMR)) `ivreg` ✓ | Burgess 2015, Stat Methods Med Res · PMID 26282889 · [doi](https://doi.org/10.1177/0962280215597579) |
| Limited-information maximum likelihood <br>*aka LIML* | I | 1S | WIV | ⛔ no code · exp 1–any · out 1–1 · pair · IV≥n_exposures (just-identified) · needs: indiv | — | Burgess 2015, Stat Methods Med Res · PMID 26282889 · [doi](https://doi.org/10.1177/0962280215597579) |
| Control function / two-stage residual inclusion <br>*aka 2SRI* | I | 1S | non-linear outcome models | R · exp 1–any · out 1–1 · pair · IV≥n_exposures (just-identified) · needs: indiv | ivtools; OneSampleMR ([cran/ivtools](https://github.com/cran/ivtools), [remlapmot/OneSampleMR](https://github.com/remlapmot/OneSampleMR)) `tsri (OneSampleMR); ivglm (ivtools)` ✓ | Burgess 2015, Stat Methods Med Res · PMID 26282889 · [doi](https://doi.org/10.1177/0962280215597579) |
| Two-stage predictor substitution <br>*aka 2SPS* | I | 1S | — | R · exp 1–any · out 1–1 · pair · IV≥n_exposures (just-identified) · needs: indiv | OneSampleMR ([remlapmot/OneSampleMR](https://github.com/remlapmot/OneSampleMR)) `tsps` ✓ | Burgess 2015, Stat Methods Med Res · PMID 26282889 · [doi](https://doi.org/10.1177/0962280215597579) |
| Generalised method of moments <br>*aka GMM* | I | 1S | — | R · exp 1–any · out 1–1 · pair · IV≥n_exposures (just-identified) · needs: indiv | ivtools; OneSampleMR; MendelianRandomization (mr_mvgmm, summary) ([cran/ivtools](https://github.com/cran/ivtools), [remlapmot/OneSampleMR](https://github.com/remlapmot/OneSampleMR)) `ivglm (ivtools; G-estimation/two-stage)` ✓ | Burgess 2015, Stat Methods Med Res · PMID 26282889 · [doi](https://doi.org/10.1177/0962280215597579) |
| Allele (genetic risk) score as single instrument <br>*aka GRS* | I | 1S | WIV | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Burgess 2015, Stat Med · PMID 26661904 · [doi](https://doi.org/10.1002/sim.6835) |
| MR for binary risk factor & outcome | I | 1S | — | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Allman 2021, Genet Epidemiol · PMID 33998053 · [doi](https://doi.org/10.1002/gepi.22387) |
| MR in the multivariate general linear model | I | 1S | — | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Allman 2021, Genet Epidemiol · PMID 34672390 · [doi](https://doi.org/10.1002/gepi.22435) |
| MR-BOIL (integrated likelihood, binary outcome) | I | 1S | — | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Shi 2023, Genet Epidemiol · PMID 36808763 · [doi](https://doi.org/10.1002/gepi.22520) |
| One-sample MR for count outcomes | I | 1S | CP; UHP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Liyanage 2024, Genet Epidemiol · PMID 39498871 · [doi](https://doi.org/10.1002/gepi.22602) |
| PLDMR (pleiotropy & LD adaptive MR) | I | 1S | LD; pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Wang 2021, Front Genet · PMID 34322150 · [doi](https://doi.org/10.3389/fgene.2021.634394) |
| MR-SPLIT (sample splitting + cross-fitting) | I | 1S | WIV; WC | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Shi 2024, PLoS Genet · PMID 39241053 · [doi](https://doi.org/10.1371/journal.pgen.1011391) |
| IV model averaging (LASSO/SCAD/MCP) | I | 1S | WIV | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Seng 2023, Stat Med · PMID 37476915 · [doi](https://doi.org/10.1002/sim.9819) |
| LASSO-type IV selection / jackknife IV | I | 1S | invalid & weak IVs | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Qasim 2024, Stat Methods Med Res · PMID 39544096 · [doi](https://doi.org/10.1177/09622802241281035) |
| Within-family (sibling / trio) MR | I/S (within-family GWAS) | 1S | POP; dynastic; assortative mating | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | Brumpton 2020, Nat Commun · PMID 32665587 · [doi](https://doi.org/10.1038/s41467-020-17117-4) |
| MR-Twin (family-based, robust to stratification) | I | 1S | POP | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | LaPierre 2023, Genome Res · PMID 37197991 · [doi](https://doi.org/10.1101/gr.277664.123) |
| Pleiotropy-robust MR (zero first-stage subgroup) | I | 1S | pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | van Kippersluis 2018, Int J Epidemiol · PMID 28338774 · [doi](https://doi.org/10.1093/ije/dyx002) |
| MR-GENIUS | I | 1S | invalid IVs | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Tchetgen Tchetgen EJ 2021, Stat Sci · [doi](https://doi.org/10.1214/20-STS802) ⚑ |
| MR-MiSTERI | I | 1S | invalid IVs | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Liu 2022, Biometrics · PMID 35950778 · [doi](https://doi.org/10.1111/biom.13735) |

### M. Mediation, time-varying & life-course MR

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Two-step MR (mediation) | S/I | 2S | mediation | ⛔ no code · exp 1–1 · out 1–1 · joint · IV≥2 per step · needs: extra data | — | Relton 2012, Int J Epidemiol · PMID 22422451 · [doi](https://doi.org/10.1093/ije/dyr233) |
| MR mediation: difference & product methods (MVMR-based) | S | 2S | mediation | R · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: extra data | TwoSampleMR; MVMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) | Carter 2021, Eur J Epidemiol · PMID 33961203 · [doi](https://doi.org/10.1007/s10654-021-00757-1) |
| Summary-data MR causal mediation (improved variance) | S | 2S | mediation | ⛔ no code · exp 1–1 · out 1–1 · joint · IV≥2 per step · needs: extra data | — | Lin 2025, Stat Med · PMID 39910926 · [doi](https://doi.org/10.1002/sim.10317) |
| Time-varying exposure MVMR | S | 2S | time-varying | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: extra data | — | Sanderson 2022, PLoS Genet · PMID 35849575 · [doi](https://doi.org/10.1371/journal.pgen.1010290) |
| Life-course MVMR (childhood vs adult) | S | 2S | time-varying | ⛔ no code · exp 2–any · out 1–1 · joint · IV≥n_exposures + 1 · needs: extra data | — | Richardson 2022, Nat Commun · PMID 35484151 · [doi](https://doi.org/10.1038/s41467-022-29932-y) |
| Structural mean models for life-course MR <br>*aka SMM* | I | 1S | time-varying | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv, extra data | — | Power 2026, Am J Epidemiol · PMID 39968723 · [doi](https://doi.org/10.1093/aje/kwaf029) |
| FA-MFMR (functional MVMR, longitudinal exposures) | I/S | — | correlated longitudinal exposures | R · exp 2–any · out 1–1 · joint · IV≥— · needs: extra data | FA-MFMR ([YQHuFD/FA-MFMR](https://github.com/YQHuFD/FA-MFMR)) | Chen 2026, Bioinformatics · PMID 42847660 · [doi](https://doi.org/10.1093/bioinformatics/btag725) |

### N. Collider, index-event & covariate-adjusted-GWAS bias

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Index-event bias adjustment (Dudbridge) | S | — | SEL (collider) | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, extra data | — | Dudbridge 2019, Nat Commun · PMID 30952951 · [doi](https://doi.org/10.1038/s41467-019-09381-w) |
| SlopeHunter | S | — | SEL (collider) | R · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, extra data | SlopeHunter ([Osmahmoud/SlopeHunter](https://github.com/Osmahmoud/SlopeHunter)) | Mahmoud 2022, Nat Commun · PMID 35110547 · [doi](https://doi.org/10.1038/s41467-022-28119-9) |
| IV-based collider-bias correction | S | — | SEL (collider) | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, extra data | — | Cai 2022, Genet Epidemiol · PMID 35583096 · [doi](https://doi.org/10.1002/gepi.22455) |
| Collider correction with multiple covariates via MVMR-cML | S | — | SEL (collider); CP | R · exp 1–any · out 1–1 · pair · IV≥n_exposures + 2 · needs: extra data | MVMR-cML ([ZhaotongL/MVMR-cML](https://github.com/ZhaotongL/MVMR-cML)) | Wang 2024, PLoS Genet · PMID 38648211 · [doi](https://doi.org/10.1371/journal.pgen.1011246) |
| MVMR with covariable-adjusted GWAS | S | 2S | collider (adjusted GWAS) | R · exp 1–any · out 1–1 · pair · IV≥n_exposures + 2 · needs: extra data | adjustedGWAS ([eleanorsanderson/adjustedGWAS](https://github.com/eleanorsanderson/adjustedGWAS)) | Gilbody 2025, Genet Epidemiol · PMID 39812504 · [doi](https://doi.org/10.1002/gepi.22606) |
| Unbiased MR with covariate-adjusted GWAS | S | 2S | collider (adjusted GWAS) | ⛔ no code · exp 1–any · out 1–1 · pair · IV≥n_exposures + 2 · needs: extra data | — | Wang 2025, HGG Adv · PMID 39891386 · [doi](https://doi.org/10.1016/j.xhgg.2025.100412) |
| Bayesian MR (BWMR / MR-Horse) for index-trait bias | S | — | SEL (index trait) | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥— · needs: full GWAS, extra data | — | Adhikari 2026, HGG Adv · PMID 41846316 · [doi](https://doi.org/10.1016/j.xhgg.2026.100592) |
| Index-event bias methods for disease progression (IPW, Heckman, SlopeHunter, MVMR) | S/I | — | SEL | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Zhang 2026, Am J Hum Genet · PMID 42372731 · [doi](https://doi.org/10.1016/j.ajhg.2026.06.004) |
| Collider-bias approach to apply two-sample methods to one-sample data | I | 1S | SEL | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥1 · needs: indiv | — | Barry 2021, PLoS Genet · PMID 34370750 · [doi](https://doi.org/10.1371/journal.pgen.1009703) |

### O. Diagnostics, sensitivity & visualisation

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Cochran's Q / Rücker's Q′ heterogeneity; Rücker model selection | S | 2S | heterogeneity | R · exp 1–1 · out 1–1 · pair · IV≥2 (IVW Q); 3 (Egger/Rücker Q′) | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_heterogeneity; mr_rucker; mr_rucker_bootstrap; mr_rucker_jackknife` ✓ | Bowden 2017, Stat Med · PMID 28114746 · [doi](https://doi.org/10.1002/sim.7221) |
| MR-Egger intercept (directional pleiotropy test) | S | 2S | DP (test) | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_pleiotropy_test` ✓ | Bowden 2015, Int J Epidemiol · PMID 26050253 · [doi](https://doi.org/10.1093/ije/dyv080) |
| I²GX & SIMEX correction for MR-Egger | S | 2S | NOME violation | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `Isq` ✓ | Bowden 2016, Int J Epidemiol · PMID 27616674 · [doi](https://doi.org/10.1093/ije/dyw220) |
| Modified (second-order) weights beyond NOME | S | 2S | NOME violation | R · exp 1–1 · out 1–1 · pair · IV≥3 | RadialMR ([WSpiller/RadialMR](https://github.com/WSpiller/RadialMR)) `ivw_radial(weights=)` | Bowden 2019, Int J Epidemiol · PMID 30561657 · [doi](https://doi.org/10.1093/ije/dyy258) |
| MR Steiger directionality test & Steiger filtering | S | 2S | reverse causation | R · exp 1–1 · out 1–1 · pair · IV≥1 · needs: extra data | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `directionality_test; steiger_filtering; mr_steiger` ✓ | Hemani 2017, PLoS Genet · PMID 29149188 · [doi](https://doi.org/10.1371/journal.pgen.1007081) |
| F-statistic / weak-instrument assessment | S/I | — | WIV | R · exp 1–1 · out 1–1 · pair · IV≥1 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `(compute from beta/se; add_rsq)` | Burgess 2011, Int J Epidemiol · PMID 21414999 · [doi](https://doi.org/10.1093/ije/dyr036) |
| Leave-one-out, single-SNP, scatter, forest, funnel & density plots | S | 2S | influential variants | R · exp 1–1 · out 1–1 · pair · IV≥3 | TwoSampleMR; MendelianRandomization ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR), [cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_leaveoneout; mr_singlesnp; mr_scatter_plot; mr_forest_plot; mr_funnel_plot; mr_density_plot; mr_loo; mr_plot` ✓ | Hemani 2018, Elife · PMID 29846171 · [doi](https://doi.org/10.7554/eLife.34408) |
| Powerful heterogeneity test (eigenvector combinations) | S | 2S | heterogeneity | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥3 | — | Wang 2024, Stat Med · PMID 39552275 · [doi](https://doi.org/10.1002/sim.10279) |
| Sensitivity statistics (robustness to residual stratification / pleiotropy) | S | 2S | POP; pleiotropy | ⛔ no code · exp 1–1 · out 1–1 · pair · IV≥2 | — | Cinelli 2022, Nat Commun · PMID 35232963 · [doi](https://doi.org/10.1038/s41467-022-28553-9) |
| Sign concordance test | S | 2S | — | R · exp 1–1 · out 1–1 · pair · IV≥6 | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `mr_sign` ✓ | - -, package documentation · - ⚑ |

### P. Study design, power & reporting

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| Three-sample genome-wide design (selection / exposure / outcome) | S | 3S | WC | ref · exp 1–1 · out 1–1 ·  · IV≥— · needs: extra data | mr.raps ([qingyuanzhao/mr.raps](https://github.com/qingyuanzhao/mr.raps)) | Zhao 2019, Int J Epidemiol · PMID 31298269 · [doi](https://doi.org/10.1093/ije/dyz142) |
| Power & sample-size calculation (mRnd) <br>*aka mRnd* | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— · needs: extra data | mRnd web app ([link](https://shiny.cnsgenomics.com/mRnd/)) | Brion 2013, Int J Epidemiol · PMID 24159078 · [doi](https://doi.org/10.1093/ije/dyt179) |
| STROBE-MR reporting guideline | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Skrivankova 2021, JAMA · PMID 34698778 · [doi](https://doi.org/10.1001/jama.2021.18236) |
| Guidelines for performing MR investigations (2023 update) | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Burgess 2023, Wellcome Open Res · PMID 32760811 · [doi](https://doi.org/10.12688/wellcomeopenres.15555.3) |
| MR primer (Nature Reviews Methods Primers) | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Sanderson 2022, Nat Rev Methods Primers · PMID 37325194 · [doi](https://doi.org/10.1038/s43586-021-00092-5) |
| Real-data benchmark of 16 two-sample MR methods | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Hu 2024, Am J Hum Genet · PMID 39059387 · [doi](https://doi.org/10.1016/j.ajhg.2024.06.016) |
| Comparison of robust MR methods (simulation) | — | — | — | ref · exp 1–1 · out 1–1 ·  · IV≥— | — | Slob 2020, Genet Epidemiol · PMID 32249995 · [doi](https://doi.org/10.1002/gepi.22295) |

### Q. Platforms, automation & AI

| Method | Data | Design | Addresses | Feasibility (tier · exposures · outcomes · multi-mode · min IVs · needs) | Implementation | Reference |
|---|---|---|---|---|---|---|
| MR-Base platform / TwoSampleMR | S | 2S | — | R · exp 1–1 · out 1–1 ·  · IV≥— | TwoSampleMR ([MRCIEU/TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR)) `extract_instruments; extract_outcome_data; harmonise_data; mr; mr_report` ✓ | Hemani 2018, Elife · PMID 29846171 · [doi](https://doi.org/10.7554/eLife.34408) |
| MendelianRandomization R package | S | 2S | — | R · exp 1–1 · out 1–1 ·  · IV≥— | MendelianRandomization ([cran/MendelianRandomization](https://github.com/cran/MendelianRandomization)) `mr_input; mr_mvinput; mr_allmethods` ✓ | Yavorska 2017, Int J Epidemiol · PMID 28398548 · [doi](https://doi.org/10.1093/ije/dyx034) |
| MR-PheWAS | S/I | — | hypothesis screening | ⛔ no code · exp 1–1 · out 2–any (phenome) · joint · IV≥— | — | Millard 2015, Sci Rep · PMID 26568383 · [doi](https://doi.org/10.1038/srep16645) |
| MRAgent (LLM-based automated MR agent) | S | 2S | — | Py · exp 1–1 · out 1–1 · pair · IV≥— | MRAgent ([xuwei1997/MRAgent](https://github.com/xuwei1997/MRAgent)) | Xu 2025, Brief Bioinform · PMID 40194554 · [doi](https://doi.org/10.1093/bib/bbaf140) |
| EpiGraphDB (MR evidence knowledge graph) | — | — | — | web · exp 1–1 · out 1–1 ·  · IV≥— | epigraphdb ([link](https://epigraphdb.org)) | Liu 2021, Bioinformatics · PMID 33165574 · [doi](https://doi.org/10.1093/bioinformatics/btaa961) |
| MRAD interactive platform (Alzheimer's disease) | — | — | — | web · exp 1–1 · out 1–1 ·  · IV≥— | web app ([link](https://gwasmrad.com/mrad/)) | Zhao 2024, Elife · PMID 39392298 · [doi](https://doi.org/10.7554/eLife.96224) |

## 6. Supporting tools & data resources (data acquisition, QC, LD, heritability)
| Resource | Role | Type | Link | Reference |
|---|---|---|---|---|
| MRC IEU OpenGWAS | GWAS summary-statistics database & API | data | [site](https://opengwas.io) | Elsworth B 2020, bioRxiv · [doi](https://doi.org/10.1101/2020.08.10.244293) ⚑ |
| ieugwasr | R client for OpenGWAS (instruments, LD clumping, LD matrix) | R package | [MRCIEU/ieugwasr](https://github.com/MRCIEU/ieugwasr) | - -, package documentation · - ⚑ |
| GWAS-VCF format & gwas2vcf | standard summary-statistics format / converter | format/tool | [MRCIEU/gwas2vcf](https://github.com/MRCIEU/gwas2vcf) | Lyon 2021, Genome Biol · PMID 33441155 · [doi](https://doi.org/10.1186/s13059-020-02248-0) |
| MungeSumstats | summary-statistics QC & standardisation | Bioconductor package | [neurogenomics/MungeSumstats](https://github.com/neurogenomics/MungeSumstats) | Murphy 2021, Bioinformatics · PMID 34601555 · [doi](https://doi.org/10.1093/bioinformatics/btab665) |
| MRInstruments | curated instrument sets (GWAS Catalog, QTLs) | R package | [MRCIEU/MRInstruments](https://github.com/MRCIEU/MRInstruments) | - -, package documentation · - ⚑ |
| NHGRI-EBI GWAS Catalog | GWAS associations & summary statistics | data | [site](https://www.ebi.ac.uk/gwas/) | Sollis 2023, Nucleic Acids Res · PMID 36350656 · [doi](https://doi.org/10.1093/nar/gkac1010) |
| FinnGen | biobank GWAS summary statistics | data | [site](https://www.finngen.fi/en) | Kurki 2023, Nature · PMID 36653562 · [doi](https://doi.org/10.1038/s41586-022-05473-8) |
| UK Biobank | individual-level biobank data | data | [site](https://www.ukbiobank.ac.uk) | Bycroft 2018, Nature · PMID 30305743 · [doi](https://doi.org/10.1038/s41586-018-0579-z) |
| GTEx v8 | tissue eQTL/sQTL | data | [site](https://gtexportal.org) | Consortium 2020, Science · PMID 32913098 · [doi](https://doi.org/10.1126/science.aaz1776) |
| eQTLGen | blood cis/trans-eQTL | data | [site](https://www.eqtlgen.org) | Võsa 2021, Nat Genet · PMID 34475573 · [doi](https://doi.org/10.1038/s41588-021-00913-z) |
| UKB-PPP | plasma pQTL (Olink) | data | — | Sun 2023, Nature · PMID 37794186 · [doi](https://doi.org/10.1038/s41586-023-06592-6) |
| deCODE plasma proteome | plasma pQTL (SomaScan) | data | — | Ferkingstad 2021, Nat Genet · PMID 34857953 · [doi](https://doi.org/10.1038/s41588-021-00978-w) |
| 1000 Genomes Phase 3 | LD reference panel | data | [site](https://www.internationalgenome.org) | Consortium 2015, Nature · PMID 26432245 · [doi](https://doi.org/10.1038/nature15393) |
| PhenoScanner / V2 | variant–trait lookup for pleiotropy screening | database | — | Kamat 2019, Bioinformatics · PMID 31233103 · [doi](https://doi.org/10.1093/bioinformatics/btz469) |
| LDlink / LDlinkR | LD proxies & LD matrices | web + R package | [CBIIT/LDlinkR](https://github.com/CBIIT/LDlinkR) | Machiela 2015, Bioinformatics · PMID 26139635; 32180801 · [doi](https://doi.org/10.1093/bioinformatics/btv402) |
| PLINK 1.9 / PLINK 2 | clumping, LD, genotype QC | CLI | [site](https://www.cog-genomics.org/plink/) | Purcell 2007, Am J Hum Genet · PMID 17701901; 25722852 · [doi](https://doi.org/10.1086/519795) |
| GCTA (incl. GSMR, COJO) | conditional analysis, GSMR | CLI | [site](https://yanglab.westlake.edu.cn/software/gcta/) | Yang 2010, Am J Hum Genet · PMID 21167468 · [doi](https://doi.org/10.1016/j.ajhg.2010.11.011) |
| LDSC (h², rg, intercepts) | heritability, genetic correlation, overlap parameters | CLI (Python) | [bulik/ldsc](https://github.com/bulik/ldsc) | Bulik-Sullivan 2015, Nat Genet · PMID 25642630; 26414676 · [doi](https://doi.org/10.1038/ng.3211) |
| LD Hub | web LDSC | web | — | Zheng 2016, Bioinformatics · PMID 27663502 · [doi](https://doi.org/10.1093/bioinformatics/btw613) |
| HDL (high-definition likelihood) | genetic correlation | R | [zhenin/HDL](https://github.com/zhenin/HDL) | Ning 2020, Nat Genet · PMID 32601477 · [doi](https://doi.org/10.1038/s41588-020-0653-y) |
| MTAG | multi-trait GWAS meta-analysis | CLI (Python) | [JonJala/mtag](https://github.com/JonJala/mtag) | Turley 2018, Nat Genet · PMID 29292387 · [doi](https://doi.org/10.1038/s41588-017-0009-4) |
| GenomicSEM | multivariate genetic SEM; munge/ldsc/hdl wrappers | R package | [GenomicSEM/GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) | Grotzinger 2019, Nat Hum Behav · PMID 30962613 · [doi](https://doi.org/10.1038/s41562-019-0566-x) |
| OneSampleMR | one-sample IV helpers | R package | [remlapmot/OneSampleMR](https://github.com/remlapmot/OneSampleMR) | - -, package documentation · - ⚑ |

## 7. Known gaps & cautions for the next design step
- **No public code located** for 78 entries (runnable_tier = no_code; e.g., mdIVW, RIVW, WWER, CMR, LCV, many one-sample and multi-tissue methods). GSMR, SMR/HEIDI, OPERA, FINEMAP and EMIC are distributed as binaries from the authors’ websites (tier external_binary, not re-checked). The GUI should only expose methods with a runnable, version-pinned implementation; the rest can stay in a "literature only" tier.
- **Script-only implementations** (not installable packages): MR-Horse/MVMR-Horse (`mr_horse.R`, needs JAGS), LMM-MEC (R scripts), MR-BMA (demo scripts). These need wrapping before GUI use.
- **Non-R tools** that need separate runtimes: GCTA-GSMR, SMR/HEIDI, OPERA, FINEMAP (C++ binaries); LDSC, MTAG (Python 2/3); MR-link/MR-link-2, MRAgent (Python).
- **Heavy external inputs**: MR-APSS, MRMU, XMR, METEOR and MRlap need LDSC-derived matrices; cis-MR, MR-Corr2, MR-LDP, MRAID, RBMR need LD matrices; one-sample methods need individual-level data (governance & security implications for a desktop app).
- **Unverified items**: penalised weighted median reference (arXiv:1606.03729) and MR-MoE (bioRxiv 173682) were not re-checked; MRMU has no paper yet (package only); XMR is a 2026 medRxiv preprint.
- **Benchmarks to drive a "method recommender"**: Hu et al. 2024 (16 methods on real data, PMID 39059387) and Slob & Burgess 2020 (PMID 32249995).

Sources: PubMed records cited by PMID/DOI above; publisher pages for Ann Stat / Stat Sci / JRSSB items; GitHub repository pages; the user's source PDFs and the UCSF MR tutorial.