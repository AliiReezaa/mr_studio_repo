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

## Correlated exposures (Round 53)

When the exposures come from one GWAS (as the GLGC lipids here), the variant-exposure associations of a variant are correlated.
The MVMR package handles this with `phenocov_mvmr(pcor, se)`: covariance = phenotypic correlation x the product of the standard
errors. MR Studio takes the same correlations as MVMR-IVW settings (one per pair of exposures, default 0) and computes the
conditional F and Q_A with them. The estimates themselves do not change. With the tutorial's illustrative matrix
(LDL-HDL -0.1, LDL-TG -0.05, HDL-TG 0.2), compared with `strength_mvmr` and `pleiotropy_mvmr` given `gencov = phenocov_mvmr(...)`:

| Quantity | MVMR package | MR Studio |
|---|---|---|
| Conditional F, LDL-c / HDL-c / triglycerides | 48.88 / 70.52 / 40.33 | 48.88 / 70.52 / 40.33 |
| Q_A | 682.84 | 682.84 |

Identical (relative difference 0). A correlation matrix that is not positive definite (for example 0.95, -0.95, 0.95) is blocked
before the run with the reason.

## Q-statistic minimisation (Round 54)

Tutorial step 6: `qhet_mvmr(F.data, mvmrcormatrix)`. MR Studio's "MVMR Q-minimisation" method minimises the same Q_A (with the
same correlations), adds an additive heterogeneity variance tau2 when the minimum exceeds its degrees of freedom, as MVMR 0.4.9
does, and offers both of the package's uncertainty methods. Compared with the package (one core; the bootstrap with MR Studio's
script seed and 50 iterations to keep the run short):

| Quantity | MVMR package | MR Studio |
|---|---|---|
| Estimates, LDL-c / HDL-c / triglycerides | -0.02431 / 0.00897 / 0.00738 | -0.02431 / 0.00897 / 0.00738 |
| Jackknife SE | 0.01407 / 0.00941 / 0.01506 | 0.01407 / 0.00941 / 0.01506 |
| Bootstrap BCa interval | (-0.061, 0.021) / (-0.017, 0.029) / (-0.024, 0.028) | identical |

Estimates and SEs agree to a relative 1e-9 (the minimisation is an optimisation, so the criterion is 1e-6); the BCa limits, which
the package rounds to 3 decimals, are identical at that precision. The jackknife takes 13 s in MR Studio against 94 s in the
package, because the Q function is computed for all variants at once. MR Studio defaults to the jackknife: the bootstrap needs one
full minimisation per iteration (several minutes for 1000). The package marks qhet_mvmr as under development; MR Studio follows
version 0.4.9, the version validated here.

## Not covered

Correlations estimated
from individual-level data (`snpcov_mvmr`) are not supported; a correlation from the literature or from the GWAS (for example the
LD score regression intercept between the exposures) has to be entered.
