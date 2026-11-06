# Calibration of three feasibility thresholds by simulation (Round 29)

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
2. The 0.2 allele-frequency threshold for palindromic variants and the 90% ancestry cut-off: no rule in the feasibility engine uses them yet.
3. The outlier-module rules in FR-104 (trigger, 20% cap, selection rule, one round): NOT BUILT.
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
3. R23: keep 10 for IVW. Optionally add to the message that MR-Egger and the median and mode methods are biased more, and that instruments selected in the same exposure GWAS can be biased even with a high F.
4. R24: keep 0.9 for MR-Egger. For radial MR-Egger either keep 0.9 (10.6% just above it) or use 0.95 for that method.
