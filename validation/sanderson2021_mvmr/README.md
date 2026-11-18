# Validation study 6: Sanderson, Spiller and Bowden 2021, MVMR instrument strength and pleiotropy tests

Reference: Sanderson E, Spiller W, Bowden J. Testing and correcting for weak and pleiotropic instruments in two-sample multivariable
Mendelian randomization. Stat Med 2021;40:5434-5452. doi:10.1002/sim.9133, PMID 34338327. Data and reference values: the authors'
MVMR R package and its tutorial ("Multivariable MR Tutorial", rawdat_mvmr: 145 variants, LDL-c, HDL-c and triglycerides from GLGC on
systolic blood pressure from UK Biobank).

Run from the repository root (needs the MVMR package, which is not part of the app runtime):

    Rscript --vanilla validation/sanderson2021_mvmr/run_validation.R

## What is compared

MR Studio runs the tutorial data as an MVMR study pair (all four MVMR methods, generated script replayed) and its numbers are
compared with the MVMR package called directly with zero covariance (`strength_mvmr(gencov = 0)`, `pleiotropy_mvmr(gencov = 0)`,
`ivw_mvmr`). Criterion: the same formulas on the same data, so a relative difference below 1e-8 for every number.

| Quantity | MVMR package | MR Studio |
|---|---|---|
| Conditional F, LDL-c / HDL-c / triglycerides | 46.98 / 68.75 / 39.34 | 46.98 / 68.75 / 39.34 |
| Q_A (142 df) | 683.08 | 683.08 |
| MVMR-IVW estimates | -0.0218 / 0.0037 / 0.0256 | -0.0218 / 0.0037 / 0.0256 |
| MVMR-IVW standard errors | 0.0142 / 0.0103 / 0.0160 | 0.0142 / 0.0103 / 0.0160 |

Result: PASSED, largest relative difference 2e-16 (`results.csv`).

## Defect found and fixed (Round 52)

Before this round MR Studio reported the conditional F statistic from `MendelianRandomization::mr_mvivw(nx = ...)`. That is a
different statistic: it needs the exposure sample sizes (so it was missing when they were unknown) and its value depends on them.
On the tutorial data it gave 40.9 / 75.8 / 40.8 with n = 100,000 instead of the 47.0 / 68.8 / 39.3 of Sanderson et al., while the
methods text cited Sanderson et al. MR Studio now computes the Sanderson et al. conditional F itself (no sample sizes needed) and
also reports their Q_A, which allows for the uncertainty of the variant-exposure associations. Cochran's Q of the MVMR-IVW fit
(685.46 here, first-order weights) is still reported beside it; the two differ by design.

## Limitation

Both statistics assume zero covariance between the variant-exposure associations, which holds when the exposures come from
non-overlapping samples. The tutorial's lipids come from one GLGC sample; with the tutorial's illustrative phenotypic correlation
matrix the package gives other values (its `phenocov_mvmr` option). MR Studio states the zero-covariance assumption in the methods text.
