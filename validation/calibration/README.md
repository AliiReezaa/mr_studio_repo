# Calibration of thresholds by simulation and real data (Rounds 29, 30 and 32)

This folder checks the three active heuristic rules of the feasibility engine (prototype/engine/feasibility_rules.csv) by simulation:

| Rule | Applies to | Threshold now | Source of the threshold |
|---|---|---|---|
| R21 | MR-Egger, weighted and simple median, weighted and simple mode, MR-PRESSO, radial IVW, radial MR-Egger | at least 10 instruments | rule of thumb |
| R23 | all causal estimates | mean F of at least 10 | Burgess 2011 (PMID 21414999) |
| R24 | MR-Egger, radial MR-Egger | I2GX of at least 0.9 | Bowden 2016 (PMID 27616674) |

No threshold was changed for R21, R23 or R24. One finding about IVW led to a change of the default IVW model (Round 30, see below).

## How to rerun

```
Rscript --vanilla validation/calibration/run_calibration.R          # full run, about 2 hours on 11 cores
Rscript --vanilla validation/calibration/run_calibration.R quick    # 8 repetitions, writes to a temporary folder
```

Full run on 2026-11-06: R 4.4.3, TwoSampleMR 0.7.13, RadialMR 1.2.5, 11 cores (parallel::mclapply, L'Ecuyer random numbers, fixed seeds per cell), 2 hours. Every method is called with the same default arguments the app uses (bootstrap 1000, phi 1). Radial IVW is computed directly by its regression; it equals RadialMR::ivw_radial to 1.39e-16 over 50 data sets (checked at the start of each run). MR-PRESSO is left out because it takes too long per data set.

## Data generation

Two-sample summary data with no sample overlap, 100,000 people in each GWAS, allele frequencies drawn from U(0.05, 0.5). For each variant the true exposure association is set by its strength mu (expected z statistic, F about mu squared), the observed exposure and outcome associations are drawn around the true values with their standard errors, and the outcome association is theta times the exposure association plus any pleiotropic effect.

## Criteria [design]

These are design choices, not values from the guidelines (Burgess et al. 2019, PMID 32760811, give no numbers):

1. False positive rate at theta = 0 of at most 0.075 for a nominal 0.05 test. With 500 repetitions the Monte Carlo SE of a rate near 0.05 is about 1.0 percentage point.
2. Dilution of MR-Egger toward the null of at most 10% (R24).
3. Relative bias of IVW of at most 10% (R23).

## Results

Full tables: results_R21_instruments.csv, results_R23_strength.csv, results_R24_i2gx.csv, results_R24_i2gx_bins.csv. Each row has the mean estimate, relative bias, empirical SD, mean SE, rejection rate at P < 0.05, 95% coverage and the Monte Carlo SE of the rejection rate.

### R21, number of instruments

k = 3, 4, 5, 6, 8, 10, 15, 20, 30 instruments, theta = 0, mu from U(6, 12), 500 repetitions per cell. Scenario "valid": no pleiotropy. Scenario "balanced": 20% of the instruments (at least 1) have a pleiotropic effect with mean 0 and SD 3 times the outcome SE.

False positive rate at P < 0.05:

| k | IVW (valid) | IVW (balanced) | Radial IVW (valid) | MR-Egger (valid) | Weighted median (valid) | Weighted median (balanced) |
|---|---|---|---|---|---|---|
| 3 | 0.194 | 0.114 | 0.054 | 0.000 | 0.024 | 0.090 |
| 4 | 0.152 | 0.124 | 0.054 | 0.000 | 0.038 | 0.058 |
| 5 | 0.116 | 0.098 | 0.042 | 0.000 | 0.032 | 0.052 |
| 6 | 0.124 | 0.084 | 0.068 | 0.008 | 0.026 | 0.030 |
| 8 | 0.092 | 0.090 | 0.050 | 0.008 | 0.032 | 0.056 |
| 10 | 0.086 | 0.064 | 0.054 | 0.012 | 0.020 | 0.026 |
| 15 | 0.056 | 0.054 | 0.046 | 0.022 | 0.028 | 0.042 |
| 20 | 0.064 | 0.058 | 0.050 | 0.026 | 0.028 | 0.034 |
| 30 | 0.060 | 0.044 | 0.052 | 0.034 | 0.028 | 0.032 |

Findings:

1. For the methods R21 covers, no false positive rate is above 0.075 at any k, except the weighted median at k = 3 in the balanced scenario (0.090). The mode methods and MR-Egger are conservative at small k (rates near 0, coverage 0.97 to 1.00), which means low power rather than false positives.
2. The precision part of the R21 message holds for MR-Egger and radial MR-Egger: their empirical SD is 4.3 and 5.4 times that of IVW at k = 10, and 12.1 times (MR-Egger) at k = 3. For the median and mode methods the SD is 0.7 to 2.0 times that of IVW at every k, so "imprecise" is mild for them.
3. The simulation shows no sharp change at 10 instruments: the false positive rates do not point to a different cut-off, and the SDs fall gradually with k. 10 stays a rule of thumb.
4. Not covered by R21, but found here: the app's IVW as it was then (TwoSampleMR::mr_ivw_mre, multiplicative random effects with no lower limit on the residual standard error) has a false positive rate of 0.194 at k = 3, 0.086 at k = 10 and at most 0.064 from k = 15 (valid scenario); coverage is 0.81 at k = 3. Radial IVW has the same estimate and, on average, the same SE, but takes its P value from a t distribution with k - 1 degrees of freedom; its rate stays between 0.042 and 0.068 at every k. The difference is the reference distribution.

### R23, instrument strength

Part A: 30 instruments with a fixed target F of 2, 3, 5, 10, 20 or 50, theta = 0.2 (bias) and theta = 0 (false positives), 500 repetitions; the mode methods are left out at theta = 0.

Relative bias at theta = 0.2:

| Target F | IVW | MR-Egger | Radial MR-Egger | Weighted median | Simple median | Weighted mode | Simple mode |
|---|---|---|---|---|---|---|---|
| 2 | -46.1% | -61.1% | -63.4% | -64.9% | -35.3% | -73.1% | -51.9% |
| 3 | -29.9% | -58.5% | -59.9% | -55.4% | -18.1% | -64.3% | -40.0% |
| 5 | -15.8% | -49.1% | -55.7% | -35.9% | -5.5% | -38.7% | -16.9% |
| 10 | -11.0% | -43.1% | -50.9% | -26.2% | -2.5% | -27.4% | -10.4% |
| 20 | -3.0% | -25.9% | -33.6% | -11.6% | 2.3% | -12.2% | -1.2% |
| 50 | -2.0% | -17.9% | -22.6% | -6.5% | 0.2% | -6.7% | -1.9% |

Part B: 2000 candidate variants with mu drawn from an exponential distribution (mean 1, 1.5, 2 or 3), selected at |z| > 5.45 (P < 5e-8) in the same exposure GWAS that gives the estimates, fast methods only, 500 repetitions.

| Mean mu | Instruments (mean) | Observed mean F | True mean F | IVW bias | IVW coverage |
|---|---|---|---|---|---|
| 1.0 | 13.8 | 42.5 | 32.6 | -13.3% | 0.90 |
| 1.5 | 66.6 | 50.6 | 43.9 | -9.2% | 0.82 |
| 2.0 | 148.0 | 59.6 | 54.4 | -6.3% | 0.79 |
| 3.0 | 343.6 | 80.7 | 77.1 | -3.4% | 0.81 |

Findings:

1. With no sample overlap, weak instruments bias every method toward the null and do not raise false positives: at theta = 0 the largest IVW rate is 0.092 at F = 3 and the rates are 0.046 to 0.070 from F = 5 upward.
2. At F = 10 the IVW bias is 11.0%, at the edge of criterion 3, and 3.0% at F = 20. The threshold of 10 fits IVW under criterion 3.
3. MR-Egger, radial MR-Egger, the weighted median and the weighted mode are biased far more: MR-Egger by 43% at F = 10 and still 18% at F = 50.
4. When instruments are selected at P < 5e-8 in the same exposure GWAS, every selected instrument has an observed F above 29.7 (5.45 squared), so R23 cannot fire, yet IVW is biased by up to 13% and coverage falls to 0.79 to 0.90. The bias comes from selection (winner's curse) and is not visible in F. The rule's note that "F above 10 does not rule out bias" is supported.

### R24, I2GX

30 instruments with mu = |10 + s * N(0, 1)| for spread s = 0.5, 1, 1.5, 2, 3, 4 and 6, theta = 0.2, 2000 repetitions per spread (14,000 data sets), grouped by the observed I2GX. Dilution is 1 minus the mean estimate divided by theta.

| Observed I2GX | Data sets | MR-Egger dilution | Radial MR-Egger dilution |
|---|---|---|---|
| below 0.5 | 17 | 20.7% | 62.0% |
| 0.5 to 0.6 | 60 | 44.1% | 99.2% |
| 0.6 to 0.7 | 401 | 25.4% | 73.4% |
| 0.7 to 0.8 | 2198 | 18.9% | 55.5% |
| 0.8 to 0.85 | 2446 | 16.2% | 47.2% |
| 0.85 to 0.9 | 2712 | 13.2% | 26.6% |
| 0.9 to 0.95 | 3175 | 7.5% | 10.6% |
| 0.95 and above | 2991 | 2.8% | 3.3% |

Findings:

1. For MR-Egger, dilution is 13.2% just below 0.9 and 7.5% just above it. The threshold of 0.9 meets criterion 2.
2. Radial MR-Egger is diluted more at the same I2GX: 10.6% between 0.9 and 0.95, just above criterion 2, and 3.3% above 0.95.
3. The two bins below 0.6 have few data sets (17 and 60), so their values are noisy.

## What this does not calibrate

1. MR-PRESSO under R21 (too slow to simulate here).
2. The 0.2 allele-frequency flip threshold and the 90% ancestry cut-off: calibrated in Round 32 (sections below). The allele-frequency check was built in Round 33; the ancestry classification is NOT BUILT.
3. The outlier-module rules in FR-104 (trigger, 20% cap, floor, selection rule): calibrated in Round 32 (section below) except the one-round rule; the module is NOT BUILT.
4. Sample overlap between the two GWAS, where weak instruments bias toward the observational estimate instead of the null.
5. Directional (non-balanced) pleiotropy under R21.

## IVW effects model with few variants (Round 30, changed)

Finding 4 of R21 was followed up with ivw_models.R (20,000 repetitions per cell, Monte Carlo SE about 0.15 percentage points; results_ivw_models.csv). The Burgess et al. guidelines (PMID 32760811) recommend multiplicative random effects as the default IVW model unless there are very few variants, because then the heterogeneity between variants cannot be estimated reliably. TwoSampleMR::mr_ivw (the IVW method TwoSampleMR::mr() runs by default) and MendelianRandomization::mr_ivw do not let the residual standard error fall below 1, so the standard error is never smaller than the fixed-effect one; TwoSampleMR::mr_ivw_mre, which MR Studio used before, has no such limit. The closed forms in the script equal MendelianRandomization::mr_ivw to 5.0e-16.

False positive rate at theta = 0:

| k | mr_ivw_mre (valid) | mr_ivw (valid) | mr_ivw_mre (balanced) | mr_ivw (balanced) | Fixed effect (balanced) |
|---|---|---|---|---|---|
| 3 | 0.187 | 0.038 | 0.147 | 0.075 | 0.323 |
| 4 | 0.148 | 0.037 | 0.117 | 0.064 | 0.278 |
| 5 | 0.118 | 0.036 | 0.099 | 0.058 | 0.247 |
| 10 | 0.081 | 0.040 | 0.066 | 0.056 | 0.255 |
| 15 | 0.071 | 0.040 | 0.061 | 0.057 | 0.260 |
| 30 | 0.058 | 0.039 | 0.055 | 0.055 | 0.263 |

Decision (Round 30): the default IVW model is now TwoSampleMR::mr_ivw for any number of variants; mr_ivw_mre and fixed effect stay available in Customise. MendelianRandomization switches to fixed effect below 4 variants; mr_ivw is kept there because fixed effect gives a false positive rate of 0.323 with 3 variants and balanced pleiotropy, against 0.075 for mr_ivw. A new feasibility rule, R21b, warns for IVW below 4 instruments that the heterogeneity cannot be estimated reliably, and the Limitations of the report say so.

## Options (no change made)

1. R21: keep 10. The simulation supports the message for MR-Egger and radial MR-Egger and shows no better cut-off.
2. IVW with few instruments: settled in Round 30 by following the guidelines (see the section above).
3. R23: kept at 10 for IVW. Done in Round 32: the message now says that MR-Egger and the weighted median and mode are biased more than IVW, and that instruments selected in the same exposure GWAS can be biased even with a high F (winner's curse).
4. R24: kept at 0.9 for MR-Egger and radial MR-Egger, the published threshold. Done in Round 32: a new information note, R24b, tells the user of radial MR-Egger that between 0.9 and 0.95 its estimate was diluted by about 11% here (MR-Egger by about 8%).

# Round 32: thresholds of checks that are not built yet

These three thresholds belong to checks the app does not have yet (marked NOT BUILT in the specification). They were calibrated now so that the numbers are ready when the checks are built. Criterion as above: bias of at most 10%, false positive rate of at most 0.075. Run on 2026-11-09 with R 4.5.3, TwoSampleMR 0.7.13, RadialMR 1.2.5, MRPRESSO 1.0, PLINK 1.9 (20260927) and the 1000 Genomes phase 3 panel MR Studio uses for clumping (503 Europeans).

## Allele-frequency flip threshold (|delta EAF| > 0.2): kept

Script: eaf_flip.R (real data). The effect-allele frequency reported in five OpenGWAS downloads (ieu-b-40 body mass index, ieu-a-8 coronary heart disease, ebi-a-GCST90029008 height, ieu-b-111 triglycerides, ebi-a-GCST006867 type 2 diabetes; all European) was compared with the frequency of the same allele in the European panel. 2,651 variants matched; 1 had other alleles and 330 were palindromic and were left out, leaving 2,320. Results: results_eaf_summary.csv, results_eaf_by_dataset.csv, results_eaf_thresholds.csv.

Correctly aligned variants differ from the panel by a median of 0.013; 99% differ by less than 0.062 and the largest difference is 0.181 (ieu-a-8). None differs by more than 0.2. If the dataset gave the frequency of the other allele (a flip), the difference becomes |1 - EAF - panel frequency|:

| Threshold | False flags (of 2,320) | Flips found, all | MAF below 0.3 | MAF 0.3 to 0.4 | MAF 0.4 to 0.5 |
|---|---|---|---|---|---|
| 0.1 | 4 (0.17%) | 87.7% | 100% | 100% | 53.8% |
| 0.15 | 1 (0.04%) | 80.9% | 100% | 100% | 28.4% |
| 0.2 | 0 | 72.9% | 100% | 95.0% | 2.9% |
| 0.25 | 0 | 66.9% | 100% | 74.1% | 0% |

Findings:

1. With 0.2, no correctly aligned European variant was flagged, and every flip with a minor allele frequency below 0.3 was found. 0.2 is kept.
2. Flips of variants with a minor allele frequency near 0.5 cannot be seen by any threshold, because both alleles then have similar frequencies. Such variants are left to the existing palindromic and allele checks.
3. The largest chance difference, 0.18, is close to 0.2, so the margin is small. These are European datasets against a European panel; with another or a mixed ancestry the differences are larger, and the check should only compare a dataset with the panel of its own ancestry.

## 90% single-ancestry cut-off: kept

Script: ancestry_cutoff.R (real frequencies, analytic bias). A GWAS adjusted for ancestry pools the within-group associations of a variant with weights proportional to each group's sample size times the variant's genotype variance. When instruments chosen in Europeans tag the causal variant less well in a share pi of the outcome sample (tagging ratio rho; rho = 0 is the worst case), the outcome association of variant j shrinks by w_Bj (1 - rho), where w_Bj = pi v_Bj / ((1 - pi) v_Aj + pi v_Bj). The formula was checked on 200,000 simulated people (estimate 0.929, formula 0.932). The genotype variances come from the 1000 Genomes frequencies of 448 body mass index instruments in all five panel ancestries, and the IVW weights from the body mass index on coronary heart disease pair. Results: results_ancestry_cutoff.csv, results_ancestry_check.csv.

IVW relative bias, worst case (rho = 0):

| Share of the other ancestry | AFR | AMR | EAS | SAS |
|---|---|---|---|---|
| 5% | -4.1% | -4.8% | -4.3% | -4.9% |
| 10% | -8.2% | -9.7% | -8.6% | -9.7% |
| 15% | -12.2% | -14.5% | -12.9% | -14.6% |
| 20% | -16.3% | -19.3% | -17.2% | -19.4% |

Findings:

1. Up to 10% of another ancestry, the IVW bias stays below 10% even when the instruments tag nothing in that ancestry; from 15% it exceeds 10%. The 90% cut-off meets the criterion and is kept.
2. With half the tagging (rho = 0.5) the bias at 10% is 4.1% to 4.9%; with rho = 0.8 it is 1.6% to 1.9%.
3. Single variants can be biased much more (up to 39% at 10%, AMR), so the cut-off is about the IVW estimate, not about any one variant.
4. Not covered: confounding by population structure, and ancestry differences in the exposure GWAS.

## Outlier-module rules (FR-90 to FR-99, validation FR-104)

Script: fr104_outliers.R (simulation). Six scenarios (valid; valid but heterogeneous, each variant with its own causal effect drawn from N(theta, 0.05 squared); balanced pleiotropy, 30% invalid; directional pleiotropy, 30% invalid; few strong outliers, 5%; many strong outliers, 30%), 10, 30, 100 and 200 instruments, theta = 0.2 and 0, 1000 repetitions each (48,000 data sets). The detectors are radial IVW and radial MR-Egger with first-order weights, computed directly and identical to RadialMR's flags on 30 data sets (0 differences); MR-PRESSO was run separately on 480 data sets. The trigger, candidate sets, cap, floor of 5 and selection rule follow FR-90 to FR-99; IVW is TwoSampleMR::mr_ivw, the app's default. Results: results_fr104_outliers.csv, results_fr104_presso.csv.

Detectors, share of valid variants flagged (valid scenario, theta = 0.2):

| Instruments | Radial IVW, alpha 0.05 | Radial IVW, alpha 0.05 / k | Radial MR-Egger, alpha 0.05 | Radial MR-Egger, alpha 0.05 / k | MR-PRESSO |
|---|---|---|---|---|---|
| 10 | 4.2% | 0.4% | 3.0% | 0.3% | 0.1% |
| 30 | 5.1% | 0.2% | 4.7% | 0.1% | 0.06% |
| 100 | 5.4% | 0.07% | 5.2% | 0.06% | not run |
| 200 | 5.4% | 0.02% | 5.2% | 0.02% | not run |

Every detector found strong outliers (alpha = +-8 SE) in 95% to 100% of cases at either level. Moderate pleiotropy is found less often: 45% to 77% at alpha 0.05, 24% to 46% at 0.05 / k, and 33% to 34% by MR-PRESSO.

IVW after the selection rule (cap 20%), coverage of the 95% interval at theta = 0.2 and false positive rate at theta = 0:

| Scenario, instruments | All instruments | Detectors at 0.05 | Detectors at 0.05 / k |
|---|---|---|---|
| Few strong outliers, 30 | 0.954 / 0.050 | 0.918 / 0.076 | 0.943 / 0.043 |
| Few strong outliers, 200 | 0.961 / 0.030 | 0.903 / 0.090 | 0.942 / 0.049 |
| Balanced pleiotropy, 100 | 0.937 / 0.060 | 0.882 / 0.108 | 0.903 / 0.087 |
| Directional pleiotropy, 30 | 0.943 / 0.053 | 0.910 / 0.088 | 0.847 / 0.141 |
| Directional pleiotropy, 200 | 0.939 / 0.057 | 0.939 / 0.058 | 0.844 / 0.157 |
| Valid, 200 | 0.939 / 0.039 | 0.938 / 0.041 | 0.939 / 0.040 |

Findings:

1. Detector level: at the package level alpha = 0.05, about 5% of valid variants are flagged at any number of instruments, so with strong true outliers the rule also removes valid variants and coverage falls to 0.90 to 0.92 (false positive rate up to 0.090). At the Bonferroni level 0.05 / k valid removals fall below 0.5% and coverage stays at 0.94 to 0.95. Recommended for the module: detectors at 0.05 / k (Bowden et al. 2018 discuss a Bonferroni level), the package default only for display.
2. Trigger (FR-90, Q P < 0.05 or Egger intercept P < 0.05): with valid instruments it fires in 8.7% (10 instruments) to 27% (200) of data sets, and in up to 54% when valid variants have heterogeneous effects, because heterogeneity grows with the number of instruments. With detectors at 0.05 / k the rule still removes something in only 2.0% to 3.0% (valid) and 3.9% to 9.4% (heterogeneous) of data sets, so the trigger can stay at 0.05.
3. Cap: 10%, 20% and 30% give nearly the same bias (within 1 percentage point). With 30% strong outliers the 20% cap keeps all instruments in every data set, and coverage is then 0.891 to 0.949, against 0.876 to 0.936 with a 30% cap. The 20% cap is kept.
4. Floor of 5: it never decided in these scenarios (with 10 instruments a 20% cap removes at most 2). It stays as a safety rule tied to the method minimums.
5. Directional pleiotropy: removing outliers made IVW worse, not better. With 30 to 200 instruments the false positive rate rose from 0.053 to 0.057 on all instruments to 0.134 to 0.157 after removal (0.05 / k), because the clearest invalid variants are removed and the remaining ones still shift the estimate, with a smaller standard error. This supports FR-100: the cleaned set is a labelled sensitivity analysis and never the primary analysis.
6. Not calibrated: the one-round rule (only one round was simulated), MR-Lasso (FR-92), and MR-PRESSO beyond 30 instruments.

How to rerun (from the repository root; the two real-data scripts need the Data library downloads and the 1000 Genomes panel):

```
Rscript --vanilla validation/calibration/fr104_outliers.R           # about 30 minutes on 11 cores
Rscript --vanilla validation/calibration/eaf_flip.R <Data library folder> <panel folder> <plink> validation/calibration
Rscript --vanilla validation/calibration/ancestry_cutoff.R <BMI dataset folder> <CHD dataset folder> <panel folder> <plink> validation/calibration
```

