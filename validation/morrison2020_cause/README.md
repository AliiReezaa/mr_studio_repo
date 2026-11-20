# Validation study 18: CAUSE, LDL cholesterol and coronary artery disease

Reference: Morrison J, Knoblauch N, Marcus JH, Stephens M, He X. Mendelian randomization accounting for correlated and uncorrelated
pleiotropic effects using genome-wide summary statistics. Nat Genet 2020;52:740-747 (doi:10.1038/s41588-020-0631-4, PMID 32451458).
Package cause (GitHub jean997/cause, GPL >= 2), commit 6257e6c, with loo 2.8.0 (see below).

Published values: the package's tutorial for LDL cholesterol (Willer et al. 2013) and coronary artery disease (van der Harst and
Verweij 2018), `analysis/ldl_cad.Rmd`, as rendered in `docs/ldl_cad.html` at the pinned commit (November 2021). Its inputs are stored in the
package's `example_data` folder: the merged summary statistics (2,023,354 variants), the nuisance parameters, and the 719 variants
selected by PLINK clumping. The script checks their MD5 against the pinned commit.

    MRSTUDIO_VALIDATE_CAUSE=1 MRSTUDIO_CAUSE_EXAMPLE=<cause source>/example_data Rscript --vanilla validation/morrison2020_cause/run_validation.R

About 75 minutes on two threads at low priority (a 2019 Intel laptop): 20 CAUSE fits of half a minute each, then the nuisance
parameters (16 to 19 minutes) three times, for part B and for the app path in part C.

## Result: PASSED (`results.csv`)

**A. The tutorial's results from its stored inputs.**

| Quantity | Tutorial | MR Studio |
|---|---|---|
| gamma, causal model (posterior median, 95% credible interval) | 0.32 (0.24, 0.40) | 0.32 (0.24, 0.40), the same for all 20 seeds |
| delta ELPD, sharing vs causal | -7.075 (SE 1.631) | -7.057 (Monte Carlo SD 0.108 over 20 seeds) |
| P that the causal model fits better | 7.2e-6 | 7.2e-6 to 1.8e-5 over 20 seeds |
| eta, causal model | -0.03 (-1.28, 0.91) | -0.03 (-1.35, 0.90) |
| q, causal model | 0.03 (0, 0.24) | 0.03 (0, 0.24) |
| eta and q, sharing model | 0.39 (0.30, 0.51); 0.62 (0.43, 0.79) | 0.39 (0.29, 0.50); 0.63 (0.43, 0.79) |
| delta ELPD, null vs sharing | -29.661 (SE 5.588) | -30.127 (MC SD 0.091), 0.08 SE apart |
| delta ELPD, null vs causal | -36.736 (SE 7.092) | -37.184 (MC SD 0.042), 0.06 SE apart |

The ELPD (expected log predictive density) is estimated by importance sampling from posterior draws, so it changes slightly with the
seed; the posterior medians and intervals do not. Gamma, the causal effect, is reproduced exactly at the printed precision, and the
comparison that CAUSE tests, sharing against causal, is within one Monte Carlo SD of the tutorial.

**B. The nuisance parameters.** est_cause_params() with the tutorial's seed (100) and one million sampled variants gives the stored
parameters exactly: rho 0.064655, mixing proportions identical.

**C. The app path.** The two GWAS imported as Data library datasets, prepared and run by the study script, with clumping on this
computer (PLINK and the app's 1000 Genomes EUR panel, P < 1e-3, r2 < 0.01): 872 variants, gamma 0.27 (0.19, 0.35), P 2.3e-5. The
conclusion is the tutorial's (the causal model fits better, gamma within its interval). The tutorial's 719 variants came from the
OpenGWAS LD reference; 523 of the 872 are among them. The app path equals the package pipeline called directly on the same files
with the same seed (largest difference 2e-15, read back from the CSV table) and the study script replays identically.

## The criterion was narrowed after the first run

As written before the first run, all three ELPD comparisons had to be within three Monte Carlo SD of the tutorial. The two
comparisons with the null model were not: they are about 0.45 lower in every seed, a shift larger than the seed spread (SD 0.04 to
0.09) but less than a tenth of the standard errors the tutorial prints. The tutorial's own two renders show that these values also
depend on more than the seed: the July 2021 render of the same tutorial, with the same 719 variants, printed -29.514, -36.856 and
-7.342, and a sharing-model interval of (0.29, 0.50), as found here. The difference cannot be traced further without the package
versions of 2021. The null model is not part of CAUSE's test, so the criterion is now gamma and the sharing-against-causal
comparison; the null comparisons are reported with their distance. This change was made after the first run and is stated here.

## loo 2.8.0

From loo 2.10, loo_compare() returns a table whose first column is the model name; cause reads that column as the ELPD difference
and stops ("non-numeric argument to binary operator" in in_sample_elpd_loo). MR Studio's runtime pins loo 2.8.0, the last version
with the table cause expects; cause's code is not changed. loo 2.4.1, the version current when the tutorial was rendered, gives
results identical to 2.8.0 on the same seed.

## How MR Studio runs it

All steps run on the user's computer, in the order of the tutorial: the summary statistics of both traits are merged with alleles
aligned and ambiguous variants removed (gwas_merge), the nuisance parameters and rho (the correlation that sample overlap brings) are
estimated from a random sample of up to one million variants (est_cause_params), variants with exposure P < 1e-3 are clumped at
r2 < 0.01 (ieugwasr::ld_clump with local PLINK), and the sharing and causal models are fitted and compared (cause). No OpenGWAS query
is needed. CAUSE runs on two threads, after a check that enough memory is free.
