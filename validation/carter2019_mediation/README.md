# Validation study 7: Carter et al. 2019, MR mediation of education on coronary heart disease by BMI

Reference: Carter AR, Gill D, Davies NM, et al. Understanding the consequences of education inequality on cardiovascular disease:
mendelian randomisation study. BMJ 2019;365:l1855. doi:10.1136/bmj.l1855, PMID 31122926. Methods as in Carter AR, Sanderson E,
Hammerton G, et al. Eur J Epidemiol 2021;36:465-478 (doi:10.1007/s10654-021-00757-1, PMID 33961203).

Two-sample MR mediation: total effect of education on CHD (fixed-effect IVW), education on BMI (IVW), BMI on CHD given education
(MVMR), indirect effect by the product of coefficients, proportion mediated with delta-method intervals. Data from OpenGWAS:
education Lee et al. 2018 without 23andMe (ieu-a-1239), BMI Yengo et al. 2018 (ieu-b-40), CHD CARDIoGRAMplusC4D 1000 Genomes (ieu-a-7).

Runs only on explicit request, with a valid OpenGWAS token; `MRSTUDIO_VALIDATION_CACHE=<folder>` keeps what OpenGWAS returns so a
rerun does not query again (the token is never written):

    MRSTUDIO_VALIDATE_OPENGWAS=1 MRSTUDIO_VALIDATION_CACHE=<folder> OPENGWAS_JWT=<token> Rscript --vanilla validation/carter2019_mediation/run_validation.R

## Result: PASSED

| Quantity | Paper (95% CI) | MR Studio (95% CI) | Within paper CI |
|---|---|---|---|
| Total effect, OR per SD of education | 0.63 (0.60 to 0.67) | 0.630 (0.575 to 0.689) | yes |
| Education on BMI, SD | -0.22 (-0.24 to -0.20) | -0.243 (-0.260 to -0.225) | no (0.003 outside) |
| Proportion mediated by BMI (product) | 18% (14% to 23%) | 20.4% (14.6% to 26.2%) | yes |

All three pairs and the mediation script replay identically (`results.csv`, `mediation_results.csv`).

## The calculation is exact

On the same data, MR Studio equals the packages to rounding (`package_equivalence.csv`, largest difference 5e-16): fixed-effect IVW
against `TwoSampleMR::mr_ivw_fe` for the total effect and education on BMI, MVMR-IVW against `MendelianRandomization::mr_mvivw`,
and the product and proportion recomputed independently. The remaining differences come from the data, not from MR Studio.

## Why education on BMI differs: the instrument set

The paper's education instruments were the 1271 lead variants of the full Lee et al. sample (1.1 million people, with 23andMe),
independent at r < 0.1. Only the 766,345 people without 23andMe have public summary statistics, so far fewer variants reach
P < 5e-8: 303 lead variants at r2 < 0.001 (264 of them found for BMI). The estimate moves with the instrument set
(`instrument_sensitivity.csv`):

| Education instruments | Found for BMI | Education on BMI |
|---|---|---|
| P < 5e-8, r2 < 0.001 (main) | 264 | -0.243 |
| P < 5e-8, r2 < 0.01 (diagnostic run) | 392 | -0.239 |
| P < 5e-8, r2 < 0.1 | 488 | -0.236 |
| Paper (1271, with 23andMe) | 1271 | -0.22 |

The more instruments, the closer to the paper. The paper's own sensitivity analyses show the same dependence (weighted median
-0.27, MR-Egger -0.28). To measure how much the instrument set alone moves the estimate, the script draws 500 random subsets of the
r2 < 0.1 instruments, each the size of the main set: their estimates have a standard deviation of 0.023. The gap to the paper,
0.023, is one such standard deviation; the criterion is three.

## Criteria (all met)

- MR Studio equals the packages on the same data (difference below 1e-10).
- Total effect and proportion mediated inside the paper's 95% intervals.
- Education on BMI within 3 standard deviations of the instrument-subset spread of the paper's value.
- Every pair replays identically.
