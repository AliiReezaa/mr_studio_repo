# Validation study 12: BWMR (Round 62)

Reference: Zhao J, Ming J, Hu X, Chen G, Liu J, Yang C. Bayesian weighted Mendelian randomization for causal inference based on summary
statistics. Bioinformatics 2020;36:1501-1508. doi:10.1093/bioinformatics/btz749, PMID 31593215. Package BWMR 0.1.1, GitHub jiazhao97/BWMR,
pinned to commit c00e818.

    Rscript --vanilla validation/zhao2020_bwmr/run_validation.R

The released function fails on R 4.2 and later (a `class(x) == "try-error"` comparison) and calls ggplot2 without importing it. The
authors give the fix in `updates/BWMR_updated.R`; it differs from the package only on that line. MR Studio's script applies that line to
the installed function at run time and attaches ggplot2. The package code is not copied into MR Studio, because it carries no licence.

Result: PASSED. On the package's ExampleData (LDL cholesterol on dyslipidemia, 76 variants), MR Studio reproduces the vignette's printed
result to every printed digit: effect 0.6177372, standard error 0.0256347, P = 2.642837e-128; replay identical (`results.csv`).
