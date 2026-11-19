# Validation study 14: LD scores and cross-trait LD score regression on UK Biobank height and BMI

References: Bulik-Sullivan BK, Loh PR, Finucane HK, et al. LD Score regression distinguishes confounding from polygenicity in
genome-wide association studies. Nat Genet 2015;47:291-295 (doi:10.1038/ng.3211, PMID 25642630). Bulik-Sullivan B, Finucane HK,
Anttila V, et al. An atlas of genetic correlations across human diseases and traits. Nat Genet 2015;47:1236-1241
(doi:10.1038/ng.3406, PMID 26414676). The regression is GenomicSEM::ldsc (Grotzinger et al. 2019, Nat Hum Behav 3:513-525,
PMID 30962613), GitHub GenomicSEM/GenomicSEM commit da95d431.

Data, each checked against its MD5 in the script:
- Neale lab UK Biobank round 2 summary statistics in LDSC format, both sexes: standing height (50_irnt) and BMI (21001_irnt), and
  the lab's published LDSC table ukb31063_h2_all.02Oct2019.
- The LDSC authors' European LD scores (eur_w_ld_chr) and HapMap3 list (w_hm3.snplist), from their CC BY 4.0 copies on Zenodo
  (doi:10.5281/zenodo.8182036 and 10.5281/zenodo.7773502), since the original Broad Institute links are offline.
- The 1000 Genomes Phase 3 EUR panel that MR Studio uses for clumping.

Runs only on request (downloads about 130 MB; computing the LD scores takes 10 to 30 minutes):

    MRSTUDIO_VALIDATE_LDSC=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/ldsc_neale_ukb/run_validation.R

## Result: PASSED (16 of 16 checks, `results.csv`)

| Part | Check | Reference | MR Studio |
|---|---|---|---|
| A | LD scores vs corrected r2 sums from PLINK --r2 (200 variants, chr22) | 0 | largest difference 8e-6 (PLINK prints 6 digits) |
| B | Computed (1000G Phase 3) vs official (1000G Phase 1) LD scores | correlation 1 | 0.9951 |
| B | Mean computed / mean official LD score | 0.908 (share of HapMap3 variants in the panel) | 0.905 |
| C | App path vs GenomicSEM::ldsc on the same variants (h2, rg, covariance, intercepts) | 0 | 2e-15 |
| C | Replay of the study script | identical | identical |
| D | lambda GC, height / BMI | 2.5032 / 2.1946 | 2.5025 / 2.1958 |
| D | Intercept, height / BMI | 1.313 (0.052) / 1.105 (0.022) | 1.298 (0.057) / 1.074 (0.019) |
| D | SNP heritability, height / BMI | 0.485 (0.022) / 0.248 (0.008), partitioned | 0.441 (0.020) / 0.244 (0.007) |
| E | h2 height / h2 BMI / rg, computed vs official LD scores | 0.441 / 0.244 / -0.170 | 0.448 / 0.248 / -0.170 |

Genetic correlation of height and BMI: -0.170 (SE 0.018), P = 3.6e-22; cross-trait intercept -0.128 (SE 0.013), as expected
for two GWAS of the same UK Biobank participants (`ldsc_height_bmi.csv`).

## What the parts show

A. The LD score estimator is that of the LDSC software (bulik/ldsc 2fdeeb3, ldscore/ldscore.py): genotypes centred, missing
values set to the mean, scaled by their standard deviation, r2 corrected as r2 - (1 - r2)/(n - 2), summed over HapMap3 variants
within the window, the variant itself included. PLINK computes r2 independently; the sums agree to PLINK's printed precision.

B. The official scores were computed from 1000 Genomes Phase 1 (379 Europeans), MR Studio's from Phase 3 (503 Europeans), with the
same 1 cM windows. 9.2% of the HapMap3 variants are not in the app's panel, and since each LD score is a sum over variants, the
computed scores are lower by about that share (0.905 against 0.908). Their correlation is 0.995.

C. Both GWAS were put into the Data library with the alleles swapped (and the sign of Z flipped) for a random half of the variants,
then run through the app: preparation, study script, replay. The result equals GenomicSEM run directly on the Neale files once
those are restricted, as `munge_sumstats --merge-alleles` does, to HapMap3 variants with matching alleles (12,107 variants per file
are not in the list or have other alleles; the official LD score files hold more variants than the list, so GenomicSEM on the raw
files would use them).

D. The Neale lab's lambda GC is model free and agrees to the third decimal. Their intercepts and SNP heritabilities come from
partitioned LDSC with the baseline-LD model, which gives a higher heritability for height than the one-annotation model of
cross-trait LDSC; all values are within two combined standard errors (height h2 1.5, the others 1.1 or less).

E. With the LD scores computed from the panel, heritabilities move by less than one standard error and the genetic correlation is
unchanged to three decimals, so the computed scores can stand in for the official ones (and give LD scores for the other four
1000 Genomes ancestries, for which no official set exists).

## Not compared

The Python LDSC software itself was not run here. GenomicSEM's ldsc is the published R reimplementation; Part C checks the app
against it, and Part D against the published values of a Python LDSC pipeline.
