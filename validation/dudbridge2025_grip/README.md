# Validation study 11: MR-GRIP (Round 61)

Reference: Dudbridge F, Voller B, Woodward RM, et al. Getting to GRIPS with MR-Egger: modelling directional pleiotropy independently of
allele coding. PLoS Genet 2025;21:e1011967. doi:10.1371/journal.pgen.1011967, PMID 41468517. Implementation: TwoSampleMR::mr_grip (in the
pinned TwoSampleMR commit, as the paper says), regression of beta_y x beta_x on beta_x^2 with an intercept, weights 1 / (beta_x^2 se_y^2);
with the intercept fixed at zero it equals IVW. The weak-instrument adjusted version (b.adj) is the second option.

    Rscript --vanilla validation/dudbridge2025_grip/run_validation.R                      # Part A only
    MRSTUDIO_VALIDATE_OPENGWAS=1 MRSTUDIO_VALIDATION_CACHE=<folder> OPENGWAS_JWT=<token> Rscript --vanilla validation/dudbridge2025_grip/run_validation.R

## Part A: exact against the package (offline)

Waterworth lipid example (LDL-c on CHD, 28 variants): MR Studio equals TwoSampleMR::mr_grip for both versions (estimate and intercept;
2.7964 and 2.9043), replays identical (`results_partA.csv`).

## Part C: the paper's Tables 5 and 6, exactly (offline)

The paper's own input data, S1 Data Tables A (urate on CHD, 31 variants) and B (BMI on CAD, 96 variants), are saved here as
`s1_table_a.csv` and `s1_table_b.csv` (PLoS Genetics, CC BY 4.0; Dudbridge et al. 2025). Every printed number is reproduced: odds
ratios and limits to the 2 decimals printed, standard errors to the 3 significant digits printed (`results_partC.csv`).

| Table | Method | Paper = MR Studio |
|---|---|---|
| 5, urate on CHD | IVW | 1.11 (1.03 to 1.20), SE 0.0401 |
| 5 | MR-GRIP | 1.04 (0.96 to 1.13), SE 0.0417 |
| 5 | MR-GRIP, weak-instrument adjusted | 1.03 (0.95 to 1.12), SE 0.0421 |
| 6, BMI on CAD | IVW | 1.44 (1.26 to 1.63), SE 0.0648 |
| 6 | MR-GRIP | 1.54 (1.25 to 1.89), SE 0.105 |
| 6 | MR-GRIP, weak-instrument adjusted | 1.56 (1.07 to 2.26), SE 0.190 |

Result: PASSED.

## Part B (diagnostic): the same analysis with instruments from OpenGWAS

| Method | Paper OR (95% CI), SE | MR Studio OR, SE | Within paper CI |
|---|---|---|---|
| IVW | 1.44 (1.26 to 1.63), 0.065 | 1.62, 0.065 | yes |
| MR-GRIP | 1.54 (1.25 to 1.89), 0.105 | 1.69, 0.115 | yes |
| MR-GRIP, weak-instrument adjusted | 1.56 (1.07 to 2.26), 0.190 | 1.75, 0.163 | yes |

MR Studio equals TwoSampleMR on the same data (difference 0) and both pairs replay identically (`results.csv`). As in the paper, MR-GRIP
is above IVW. The estimates are higher than the paper's because the instruments differ: OpenGWAS gives 77 GIANT lead variants for
BMI (ieu-a-2, clumped at r2 < 0.001), while the paper used the 97 published Locke et al. variants (96 with CAD effects), listed in its
S1 Data. The standard errors agree closely, which points to the variant set rather than the calculation. With the paper's own variants (Part C) the
numbers are reproduced exactly, so the Part B differences come from the instrument set alone.
