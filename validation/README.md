# Validation of MR Studio

Every method in MR Studio is checked against a published result, against the reference package it wraps, or both. Each folder below
holds a script that reruns its check from the repository root, the data or the exact download it uses (with checksums), the numbers it
produced (`results.csv`) and a README that explains the criteria and any difference. Checks that need large downloads or an OpenGWAS
token run only on request (an environment variable named in each README).

Two kinds of check are used:

- **Reproduction**: MR Studio reproduces numbers printed in a paper or in the package's documentation (paper tables, vignettes,
  stored example results), within the precision printed or within the uncertainty stated.
- **Equivalence**: on the same input, MR Studio's result equals the reference package called directly, usually to the last digit,
  and the study script it writes replays identically.

| Folder | Methods | Reference | Data | Check | Result |
|---|---|---|---|---|---|
| `bowden2016` | IVW, MR-Egger, weighted and simple median | Bowden et al. 2016, Genet Epidemiol (PMID 27061298), Table 4 | Lipids and CAD, figshare (CC BY 4.0) | Reproduction | Passed |
| `hemani2018` | IVW and robust methods on OpenGWAS data, local and server clumping | Hemani et al. 2018, eLife (PMID 29846171) | LDL-c and CHD, OpenGWAS | Reproduction, equivalence | Passed |
| `burgess2015_mvmr` | MVMR-IVW (fixed and random effects) | Burgess and Thompson 2015, Am J Epidemiol (PMID 25632051) | 28 Waterworth et al. variants | Reproduction | Passed |
| `rees2017_mvmr` | MVMR-IVW, MVMR-Egger with orientation | Rees, Wood and Burgess 2017, Stat Med (PMID 28960498), Table 2 | 185 lipid variants (Do et al. 2013) | Reproduction | Passed |
| `sanderson2021_mvmr` | Conditional F, Q_A, Q-minimisation, phenotypic and individual-level covariances | Sanderson, Spiller and Bowden 2021, Stat Med (PMID 34338327); MVMR package | MVMR tutorial data (145 variants) | Equivalence | Passed (largest relative difference 2e-16) |
| `carter2019_mediation` | Two-step and multivariable MR mediation | Carter et al. 2019, BMJ (PMID 31122926) | Education, BMI and CHD, OpenGWAS | Reproduction, equivalence | Passed |
| `patel2023_mr_methods` | Maximum likelihood, mode-based estimate, MR-cML, MR-Robust, penalised IVW, MVMR-cML | Patel et al. 2023, Wellcome Open Res (PMID 37915953) | Waterworth et al. lipid example | Equivalence | Passed (every difference 0) |
| `twosamplemr_methods` | Bootstrap MR-Egger, penalised weighted median, NOME modes, Rucker selection, sign test | TwoSampleMR | Waterworth et al. lipid example | Equivalence | Passed (every difference 0) |
| `mrmix` | MRMix | Qi and Chatterjee 2019, Nat Commun (PMID 31028273) | Education on BMI and on CHD, OpenGWAS | Equivalence | Passed |
| `dudbridge2025_grip` | MR-GRIP | Dudbridge et al. 2025, PLoS Genet (PMID 41468517), Tables 5 and 6 | The paper's S1 Data (CC BY) | Reproduction, equivalence | Passed (every printed digit) |
| `zhao2020_bwmr` | BWMR | Zhao et al. 2020, Bioinformatics (PMID 31593215) | Package example data | Reproduction | Passed (vignette values exactly) |
| `xue2022_bidir_cml` | BiDirMRcML, BiDirCDcML | Xue and Pan 2022, PLoS Genet (PMID 35576237) | Simulated bidirectional data | Equivalence, recovery of known effects | Passed |
| `ldsc_neale_ukb` | LD scores (computed from 1000 Genomes), cross-trait LDSC | Bulik-Sullivan et al. 2015, Nat Genet (PMIDs 25642630, 26414676) | Neale lab UK Biobank height and BMI | Reproduction, equivalence | Passed (16 of 16 checks) |
| `mounier2023_mrlap` | MRlap | Mounier and Kutalik 2023, Genet Epidemiol (PMID 37036286) | The package's Example A and its stored result | Reproduction, equivalence | Passed |
| `darrous2021_lhcmr` | LHC-MR (adaptive, resumable runner) | Darrous et al. 2021, Nat Commun (PMID 34907193), Supplementary Data 2 | Neale lab UK Biobank BMI and height | Reproduction, equivalence | Passed (causal effects both ways and intercept within 1 SE of the paper) |
| `hu2022_mrapss` | MR-APSS | Hu et al. 2022, PNAS (PMID 35787050) | Neale lab UK Biobank BMI and height; package example | Reproduction, equivalence | Passed (both directions within 2 SE of the paper) |
| `morrison2020_cause` | CAUSE | Morrison et al. 2020, Nat Genet (PMID 32451458) | The package tutorial: LDL cholesterol and CAD, its stored inputs | Reproduction, equivalence | Passed (gamma identical; sharing vs causal ELPD within 1 Monte Carlo SD) |
| `calibration` | Feasibility thresholds (instrument counts, strength, allele frequency, ancestry, outliers, IVW models) | Simulation and real data | Simulated and OpenGWAS data | Calibration | See README |

Validation study numbers in the folder READMEs follow the order in which the studies were written; number 3 was not used.

## What "Passed" means

Each README states its criteria before giving the result. Where a criterion was changed after a first run, the README says so and why
(for example `mounier2023_mrlap`, part C, and `darrous2021_lhcmr`). Where a published value could not be reproduced exactly, the README
names the cause (a different instrument set, a newer package version, a smaller sample in public data) and measures its size.

## Running all the offline checks

    for s in bowden2016 burgess2015_mvmr rees2017_mvmr patel2023_mr_methods twosamplemr_methods dudbridge2025_grip zhao2020_bwmr xue2022_bidir_cml; do
      Rscript --vanilla validation/$s/run_validation.R; done

The others need downloads (UK Biobank summary statistics, LD scores) or an OpenGWAS token; their READMEs give the command.
