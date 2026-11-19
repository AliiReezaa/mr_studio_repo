# Validation study 8: MendelianRandomization estimators added in Round 58

Reference: Patel A, Ye T, Xue H, et al. MendelianRandomization v0.9.0: updates to an R package for performing Mendelian
randomization analyses using summarized data. Wellcome Open Res 2023;8:449. doi:10.12688/wellcomeopenres.19995.2, PMID 37915953.

Methods: two-sample maximum likelihood (Burgess et al. 2013), mode-based estimate (Hartwig et al. 2017), MR-cML with data
perturbation (Xue et al. 2021, PMID 34214446), MR-Robust (robust regression with penalised weights; Rees et al. 2019, PMID 31545794),
penalised IVW (Xu et al. 2023, PMID 35942938) and MVMR-cML (Lin et al. 2023, PMID 36948188).

Data: the paper's worked example, the 28 Waterworth et al. variants for LDL-c, HDL-c and triglycerides with CHD (n = 17,723),
bundled with MendelianRandomization 0.10.0.

    Rscript --vanilla validation/patel2023_mr_methods/run_validation.R

## What can be compared, and the result

The paper prints the calls for these methods (for example `mr_cML(mr_input(ldlc, ldlcse, chdlodds, chdloddsse), n = 17723)`),
but their outputs were not printed in the text. So MR Studio, run as a study with replayed scripts, is compared with the package
called exactly as the paper shows (with MR Studio's seed for the stochastic methods). The only results the text does print, the
MVMR-IVW table, are compared with the paper.

| Method | MR Studio | Package | SE (both) |
|---|---|---|---|
| Maximum likelihood | 3.2249 | 3.2249 | 0.5690 |
| MBE | 2.9463 | 2.9463 | 0.9994 |
| MR-cML (DP) | 2.8506 | 2.8506 | 0.3930 |
| MR-Robust | 2.5759 | 2.5759 | 0.2515 |
| pIVW | 2.9319 | 2.9319 | 0.5276 |
| MVMR-cML, LDL-c / HDL-c / TG | 1.955 / -0.342 / 0.812 | identical | 0.555 / 0.690 / 0.271 |
| MVMR-IVW vs the paper's printed table | 1.925 / -0.590 / 0.723 | 1.925 / -0.590 / 0.723 | 0.439 / 0.555 / 0.230 |

Result: PASSED. Every difference is 0; both pairs replay identically (`results.csv`). The units are log odds of CHD per unit of
each lipid, as in the paper.
