# Validation study 9: TwoSampleMR estimators and diagnostics added in Round 59

Methods: MR-Egger with bootstrap standard errors (Bowden et al. 2015), penalised weighted median (Burgess et al. 2017, Epidemiology
28:30-42, PMID 27749700), simple and weighted mode with NOME standard errors (Hartwig et al. 2017), Rucker model selection with its
bootstrap and jackknife versions (heterogeneity framework: Bowden et al. 2018, Am J Epidemiol 187:2681-2685, PMID 30188969), and the
sign concordance test (TwoSampleMR). The density plot of the single-variant estimates is drawn next to the funnel plot.

Data: the published Waterworth et al. lipid example (28 variants, LDL-c on CHD), as bundled with MendelianRandomization and used in
validation studies 4 and 8.

    Rscript --vanilla validation/twosamplemr_methods/run_validation.R

Criterion: MR Studio, run as study pairs with replayed scripts, equals TwoSampleMR (pinned commit 5e93ee7) called directly with the
same data, settings and seed. Result: PASSED, every difference 0, every replay identical (`results.csv`).

| Method | MR Studio = TwoSampleMR |
|---|---|
| MR-Egger (bootstrap) | 3.6851 |
| Penalised weighted median | 2.6769 |
| Simple mode (NOME) | 2.4056 |
| Weighted mode (NOME) | 2.9463 (equal to the MendelianRandomization mode-based estimate of study 8) |
| Rucker, selected model / bootstrap median / jackknife median | 2.8342 / 3.0451 / 3.1356 |
| Sign concordance | 0.185 |

Note: TwoSampleMR's jackknife reports the interval of the Rucker mean with its limits in reverse order (lower above upper); MR Studio
puts the limits in order.
