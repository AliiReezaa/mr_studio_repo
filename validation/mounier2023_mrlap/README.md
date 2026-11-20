# Validation study 16: MRlap on the package's Example A

Reference: Mounier N, Kutalik Z. Bias correction for inverse variance weighting Mendelian randomization. Genet Epidemiol
2023;47:314-331 (doi:10.1002/gepi.22522, PMID 37036286). Package MRlap (GitHub n-mounier/MRlap, GPL-2), commit 660f0268.

Data: the package's own Example A, bundled with it: simulated BMI and systolic blood pressure, about 100,000 people each and no
sample overlap, 1.15 million variants; and the package's stored result of that example (data/A.RDS). LD scores and HapMap3 list:
the LDSC authors' European files (Zenodo, CC BY 4.0), as for LDSC. Nothing else is downloaded.

    MRSTUDIO_VALIDATE_MRLAP=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/mounier2023_mrlap/run_validation.R

About 20 minutes on one core (seven MRlap runs of about 2.3 minutes each, 1.8 GB peak).

## Result: PASSED (`results.csv`)

| Part | Check | Reference | MR Studio |
|---|---|---|---|
| A | MRlap called directly: observed and corrected effects, the 39 instruments, all LDSC and genetic-architecture values | stored Example A | identical (difference 0) |
| A | Standard error of the corrected effect (simulated by MRlap) | 0.05276 | 0.05271 (SD 0.00024 over 5 seeds) |
| B | App path (Data library, preparation, study script) vs stored Example A | stored Example A | largest difference 4e-7 |
| C | App path vs MRlap called directly on the same prepared files with the same seed | 0 | 3e-16 |
| D | Replay of the study script | identical | identical |

Observed IVW effect of BMI on SBP 0.0856 (SE 0.0398, 39 instruments); corrected for weak instruments and winner's curse 0.1146
(SE 0.053); difference P = 0.04 (`mrlap_example_a.csv`).

## Notes

- MRlap simulates the standard error of the corrected effect (random draws) and the package example sets no seed, so that value
  cannot be reproduced to the digit; its spread over seeds is measured instead, and the stored value lies within it. MR Studio's
  script sets the seed (after loading ieugwasr, which draws a random number when it loads), so its own results replay exactly.
- Part B differs from the stored result in the seventh decimal: the app gives MRlap the effect and standard error, from which MRlap
  computes Z, while the example file also carries a Z column rounded to five significant digits, which the package uses when present.
- Part C is at the precision of a number written to and read back from the app's CSV results table.
- The criterion for Part C was first set to exactly 0, which a CSV round trip cannot meet; it was changed to below 1e-10 after the
  first run, for that reason only.
