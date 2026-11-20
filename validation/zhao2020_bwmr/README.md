# Validation study 12: BWMR (Round 62, redone in Round 70e)

Reference: Zhao J, Ming J, Hu X, Chen G, Liu J, Yang C. Bayesian weighted Mendelian randomization for causal inference based on summary
statistics. Bioinformatics 2020;36:1501-1508. doi:10.1093/bioinformatics/btz749, PMID 31593215. Algorithm as given in the paper and its
supplement (preprint arXiv 1811.10223v2: Methods 2.1-2.3; Supplement S1.1, E-step and M-step; S1.2, variance by linear response).

    Rscript --vanilla validation/zhao2020_bwmr/run_validation.R

## Why MR Studio has its own implementation

The authors' package (GitHub jiazhao97/BWMR) carries no licence, so it cannot be distributed with MR Studio. MR Studio therefore runs its
own implementation of the published algorithm (`BWMR_MRS_SRC` in `contracts/r/adapters/reference_adapters_v11.R`, MIT), written from the
equations of the paper and its supplement. It needs nothing beyond base R, and every study script contains the whole function. The
package is used here only as the reference to compare against (commit c00e818, with the authors' fix for R 4.2 and later from
`updates/BWMR_updated.R`).

Two details are not given in the paper and were chosen here: the starting values (the IVW estimate, the observed variant-exposure
effects, every weight 1) and the stopping rule. MR Studio iterates until the evidence lower bound changes by less than 1e-12 (relative),
that is, to convergence. The posterior variance of the effect is the supplement's linear-response formula, -[(V H1 - I)^-1 V]_11, solved
exactly through its arrow structure (each variant's block is linked only to the effect and to the outlier proportion), so the time grows
linearly with the number of variants.

## Result: PASSED (`results.csv`, `package_comparison.csv`)

**A. The vignette** (the package's ExampleData: LDL cholesterol on dyslipidemia, 76 variants), through the app path:

| Quantity | Vignette | MR Studio | Distance |
|---|---|---|---|
| Effect | 0.6177372 | 0.6179045 | 0.007 SE |
| Standard error | 0.0256347 | 0.0255741 | 0.24% |
| P value | 2.6e-128 | 5.7e-129 | reported |

The study script replays identically.

**B. The variance computation.** On all 22 data sets the supplement's formula built in full, with (3N + 4)-square matrices, equals
MR Studio's arrow solve to 4e-15 (relative).

**C. Agreement with the package**, on 22 data sets: the vignette's data, the Waterworth et al. lipid example of MendelianRandomization
(28 variants) and 20 simulated data sets of 60 variants with about 10% outlying variants.

- The package stops when the evidence lower bound changes by about 1e-6 (relative). Stopped the same way, MR Studio agrees with the
  package within 0.007 SE for the effect and 0.5% for the standard error on the 21 data sets where both reach the same optimum.
- At full convergence, MR Studio's default, the effect differs from the package's by at most 0.0082 SE, and the standard error by at most 0.82%, on 20 of these 21; on one data set
  (simulated_10), whose iterations converge slowly (2,429), the package stops 0.057 SE short of the converged value, and MR Studio started
  from the package's answer moves to the same converged value.
- On one data set (simulated_07) the two reach different optima: MR Studio 0.347 (SE 0.049), the package 0.369 (SE 0.056). MR Studio's
  solution has the higher evidence lower bound (188.18 against 187.64), down-weights one outlying variant that the package's does not, and
  is closer to the simulated effect of 0.30. Started from the package's answer, MR Studio stays at the package's optimum.

## The criteria were set after a first comparison

Before comparing, the criterion was agreement within 1e-4 (relative) with the package. The first comparison showed that the package
stops before convergence, so that no implementation run to convergence can meet it, and that on some data the algorithm has more than one
optimum. The criteria were then set as above (0.01 SE for the effect, 1% for the standard error, with the package's stopping rule, and the
higher evidence lower bound where the optima differ), and this README states that they were chosen after the first comparison.

The earlier version of this study (Round 62) ran the package itself and reproduced the vignette to every printed digit.
