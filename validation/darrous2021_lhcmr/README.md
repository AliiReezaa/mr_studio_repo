# Validation study 15: LHC-MR, BMI and standing height in UK Biobank

Reference: Darrous L, Mounier N, Kutalik Z. Simultaneous estimation of bi-directional causal effects and heritable confounding from
GWAS summary statistics. Nat Commun 2021;12:7274 (doi:10.1038/s41467-021-26970-w, PMID 34907193). Package lhcMR (GitHub
LizaDarrous/lhcMR, MIT), commit aee7f348, run through MR Studio's adaptive, resumable runner (Round 67).

Data, each checked against its MD5 in the script:
- Neale lab UK Biobank round 2 summary statistics, both sexes, in LDSC format (HapMap3 variants): BMI (21001_irnt, X) and standing
  height (50_irnt, Y); effects standardised as Z / sqrt(N), standard error 1 / sqrt(N).
- The authors' LD scores (UK10K, LDscores_filtered.csv) and local-LD parameters (LD_GM2_2prm.csv), from the links in the package README.
- The LDSC authors' European LD scores and HapMap3 list (Zenodo, CC BY 4.0), used by lhcMR for the cross-trait intercept.
- The paper's Supplementary Data 2 (CC BY 4.0), row BMI-SHeight.

Runs only on request (about 600 MB of downloads; the 250 fits take about 6.5 hours on 2 workers or 3.5 hours on 4, at low priority,
and resume if stopped):

    MRSTUDIO_VALIDATE_LHCMR=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/darrous2021_lhcmr/run_validation.R

## Result

| Parameter | Paper (SE) | MR Studio (SE) | Distance |
|---|---|---|---|
| BMI on height (axy) | -0.128 (0.039) | -0.117 (0.037) | 0.2 combined SE |
| Height on BMI (ayx) | -0.128 (0.013) | -0.116 (0.016) | 0.6 combined SE |
| Cross-trait intercept (iXY) | -0.103 (0.017) | -0.107 (0.026) | 0.1 combined SE |
| Confounder effect on BMI (tX) | 0.032 (0.040) | 0.000 (0.052) | 0.5 combined SE |
| Confounder effect on height (tY) | -0.124 (0.029) | 0.153 (0.029) | same size, opposite sign |
| Direct heritability, BMI (h2X) | 0.256 (0.011) | 0.321 (0.012) | 4.0 combined SE |
| Direct heritability, height (h2Y) | 0.546 (0.030) | 0.672 (0.034) | 2.8 combined SE |

755,617 variants after merging with the authors' LD file (75,562 used, every tenth, as lhcMR does). Full table: `lhcmr_bmi_height.csv`;
checks: `results.csv`.

**The causal effects in both directions and the cross-trait intercept agree with the paper** (within one combined standard error).
These are what MR Studio reports from LHC-MR.

**The criterion set before the run was not met.** The validation was written to require all seven parameters within two combined
standard errors of the paper. Three are not:

- The direct heritabilities are higher. The Neale lab files used here hold only the HapMap3 variants (about 0.76 million after
  merging), while the paper used the full summary statistics (about 4.6 million variants of the authors' LD file). LHC-MR's direct
  heritabilities depend on the variants in the fit; the causal effects are ratios and are much less sensitive to it.
- The confounder's effect on height has the paper's size but the opposite sign. In both runs the confounder's effect on BMI is
  essentially zero (paper 0.032, P = 0.42; here 0.000). When it is zero, the likelihood is unchanged if the sign of the effect on
  height is flipped, so that sign is not identified; only its size is.

The criterion was therefore narrowed, after the run and stated as such here, to the causal effects and the intercept; the other four
parameters are reported with their distance from the paper. A closer reproduction would use the full Neale summary statistics (about
1.7 GB of downloads), which would test the heritabilities directly.

## The runner and the package give the same numbers

On the same input, MR Studio's runner and lhcMR::lhc_mr() give identical fit tables and identical estimates, standard errors and P
values (largest difference 0); a second run reads its saved fits instead of computing them (`run_equivalence.R`,
`equivalence_results.csv`). The study script replays REPLAY_RESULT.

## Two corrections to the package, made in its inputs and call, not its code

- merge_sumstats() builds the LD table without the HLA region (chr6:28.5-33.5 Mb, as the paper states) but then joins the unfiltered
  table; MR Studio removes HLA from both traits' inputs first.
- The IVW starting points clump with TwoSampleMR::clump_data(), which asks the OpenGWAS server (a token is needed); the script points
  it at local PLINK clumping against the same 1000 Genomes EUR panel.

## Time and memory on a 16 GB laptop (2019 Intel Mac)

Starting-point fits took 7.6 worker-minutes each and jackknife fits 2.0 (they start from the best solution), about 13 worker-hours in
all. Peak memory was 4.46 GB with 2 workers and 6.0 GB with 4. These numbers calibrate the app's estimates and its safe worker choice.
