# Validation study 10: MRMix (Round 60)

Reference: Qi G, Chatterjee N. Mendelian randomization analysis using mixture models for robust and efficient estimation of causal
effects. Nat Commun 2019;10:1941. doi:10.1038/s41467-019-09432-2, PMID 31028273. Package MRMix 0.1.0 (GitHub gqi/MRMix, pinned to the
commit the runtime installs).

Data: real GWAS summary data, the education instruments of validation study 7 (Lee et al. 2018) with BMI (Yengo et al. 2018,
continuous) and with CHD (CARDIoGRAMplusC4D, binary). The script reads the OpenGWAS data saved by study 7 and sends no query:

    MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/mrmix/run_validation.R

## Scale

MRMix works on standardised associations: per SD of the genotype, and per SD of the trait for continuous traits (MRMix::standardize).
MR Studio reports MRMix on the same per-unit scale as the other methods: the standardised estimate times the median implied SD of the
outcome (continuous outcome only), divided by the median implied SD of the exposure, with SD = se x sqrt(n x 2p(1 - p)). The
standardised estimate, the proportion of variants without pleiotropy (pi0) and the conversion factor are kept as diagnostics. The
conversion is an approximation: the implied SD varies a little between variants.

## Result: PASSED

| Outcome | Variants | Standardised: MR Studio = TwoSampleMR::run_mrmix | Per unit (SE) |
|---|---|---|---|
| BMI (continuous) | 264 | -0.30 | -0.288 (0.119) |
| CHD (binary) | 102 | -0.24 | -0.228 (0.482) |

Identical to the package's reference wrapper on the same harmonised data; the per-unit value equals the formula above; replays
identical (`results.csv`). Fewer CHD variants: allele frequencies are missing for many of them in these data, and MRMix needs them.
