# Validation study 4: Burgess and Thompson 2015, multivariable MR

Replication of the applied example in Burgess S, Thompson SG. *Multivariable Mendelian randomization: the use of pleiotropic genetic
variants to estimate causal effects.* Am J Epidemiol 2015;181:251-60, doi:10.1093/aje/kwu283 (PMID 25632051; full text checked in
PubMed Central, PMC4325677). The example estimates the direct effects of LDL cholesterol, HDL cholesterol and triglycerides on coronary
heart disease (CHD).

## Data

The 28 variants of Waterworth et al. 2010 used in the paper, as distributed with MendelianRandomization 0.10.0 (`ldlc`, `hdlc`, `trig`,
`chdlodds` and their standard errors). The lipids are log-transformed, so results are odds ratios of CHD per 30% reduction in each lipid,
OR = exp(beta x log(0.7)), as in the paper.

## How to run

From the repository root:

    Rscript --vanilla validation/burgess2015_mvmr/run_validation.R

The script builds MR Studio's multivariable input (contract input kind `multi_iv`) and runs MVMR-IVW, MVMR-Egger, the MVMR median and
MVMR-Lasso through the generated scripts with the replay check, once with the fixed-effect and once with the random-effects model. It also
reconstructs the paper's own estimator (the likelihood of its equation 1 with zero correlations, maximised numerically), to show where any
difference comes from. Output: `results.csv` (MVMR-IVW against the paper), `results_robust_methods.csv`.

## Result (MVMR-IVW, fixed effect)

| Lipid | Paper (likelihood, WinBUGS) | Paper's likelihood, reconstructed | MR Studio MVMR-IVW |
|---|---|---|---|
| LDL-c | 0.50 (0.40, 0.62) | 0.50 (0.40, 0.63) | 0.50 (0.41, 0.62) |
| HDL-c | 1.22 (0.91, 1.63) | 1.23 (0.92, 1.65) | 1.23 (0.94, 1.62) |
| Triglycerides | 0.77 (0.68, 0.87) | 0.77 (0.68, 0.87) | 0.77 (0.69, 0.86) |

Every script replays identically and MR Studio equals MendelianRandomization::mr_mvivw exactly (difference 0).

## Why the numbers differ slightly, and the criteria

The paper did not use MVMR-IVW: it fitted a likelihood that also allows for uncertainty in the lipid associations, in a Bayesian framework.
Reconstructed by maximum likelihood on the same data, that estimator gives the paper's numbers to within 0.01 (HDL-c 1.23 against 1.22, the
residual coming from posterior summaries of the Bayesian fit). MVMR-IVW, which ignores the uncertainty in the lipid associations, differs
from the paper by at most 0.014 on the odds ratio scale. Criteria (all met): odds ratios within 0.015 of the paper, fixed-effect interval limits
within 0.035, the reconstructed likelihood within 0.011, identical replay, and equality with the package.

The random-effects intervals are wider than the paper's (for example HDL-c 0.84 to 1.82) because Cochran's Q is 51.4 on 25 degrees of freedom
(residual standard error 1.43) and the paper's model has no term for excess heterogeneity. This is the expected behaviour of MR Studio's
default model, not an error.
