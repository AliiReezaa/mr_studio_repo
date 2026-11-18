# Validation study 5: Rees, Wood and Burgess 2017, multivariable MR-Egger

Replication of Table 2 of Rees JMB, Wood AM, Burgess S. *Extending the MR-Egger method for multivariable Mendelian randomization to correct
for both measured and unmeasured pleiotropy.* Stat Med 2017;36:4705-18, doi:10.1002/sim.7492 (PMID 28960498; full text checked in PubMed
Central, PMC5725762): MVMR-IVW and MVMR-Egger (with the variants oriented to the HDL-c, LDL-c or triglyceride increasing allele) for the
effects of HDL-c, LDL-c and triglycerides on CHD, log odds ratio per SD, with the 185 lipid variants of Do et al. 2013.

## Data

* Lipid associations: `validation/bowden2016/dodata.csv` (figshare doi:10.6084/m9.figshare.1116328), the 185 variants of Do et al. 2013.
  The file gives betas and P values; standard errors are derived as |beta| / z(P/2), as in validation 1.
* CHD: the paper used the 2015 CARDIoGRAMplusC4D release. `dodata.csv` holds the earlier CARDIoGRAM associations (Part A). Part B takes
  the 2015 release from OpenGWAS (`ieu-a-7`, Nikpay et al. 2015) for the same variants.

## How to run

From the repository root:

    Rscript --vanilla validation/rees2017_mvmr/run_validation.R

Part B needs an OpenGWAS token and runs only when asked for explicitly:

    MRSTUDIO_VALIDATE_OPENGWAS=1 OPENGWAS_JWT=<your token> Rscript --vanilla validation/rees2017_mvmr/run_validation.R

The script imports the three lipids (with full summary statistics, declared forward strand because all columns come from one table with
the same effect allele) and the outcome through MR Studio's importer, assembles the multivariable data with `library_mv_input` (the union of
the variants associated with any lipid at P < 5e-8; joint clumping is skipped because the paper used the 185 variants as given), and runs
MVMR-IVW and MVMR-Egger in three orientations through the generated scripts with the replay check. Output: `results.csv`.

## Result

Part A (earlier CARDIoGRAM release, 183 variants): every script replays identically and MR Studio equals MendelianRandomization exactly, but
the estimates differ from Table 2 by 0.05 to 0.15. Part B (2015 release, 182 variants):

| Analysis | Lipid | Paper | MR Studio | Difference |
|---|---|---|---|---|
| MVMR-IVW | HDL-c | -0.039 | -0.046 | -0.007 |
| MVMR-IVW | LDL-c | 0.375 | 0.388 | +0.013 |
| MVMR-IVW | Triglycerides | 0.173 | 0.161 | -0.012 |
| MVMR-Egger, HDL-c orientation | HDL-c | 0.036 | 0.034 | -0.002 |
| MVMR-Egger, HDL-c orientation | LDL-c | 0.378 | 0.391 | +0.013 |
| MVMR-Egger, HDL-c orientation | Triglycerides | 0.136 | 0.120 | -0.016 |
| MVMR-Egger, HDL-c orientation | Intercept | -0.005 | -0.006 | -0.001 |
| MVMR-Egger, LDL-c orientation | HDL-c | -0.034 | -0.041 | -0.007 |
| MVMR-Egger, LDL-c orientation | LDL-c | 0.420 | 0.436 | +0.016 |
| MVMR-Egger, LDL-c orientation | Triglycerides | 0.194 | 0.184 | -0.010 |
| MVMR-Egger, LDL-c orientation | Intercept | -0.003 | -0.003 | 0.000 |
| MVMR-Egger, TG orientation | HDL-c | -0.018 | -0.027 | -0.009 |
| MVMR-Egger, TG orientation | LDL-c | 0.350 | 0.365 | +0.015 |
| MVMR-Egger, TG orientation | Triglycerides | 0.083 | 0.082 | -0.001 |
| MVMR-Egger, TG orientation | Intercept | 0.005 | 0.004 | -0.001 |

## Problems found and how each was resolved

1. **A defect in MR Studio (fixed).** The first run kept only 177 of the variants. When the multivariable data were assembled, every
   exposure was looked up with TwoSampleMR::harmonise_data action 2, which drops palindromic variants when there is no allele frequency,
   even though the analyst had declared the files to be on the forward strand (the pairwise path respects that declaration, Round 26).
   `library_mv_input` now uses action 1 when every exposure is declared forward strand, records the choice in the ledger, and an engine
   test covers both cases. After the fix: 183 variants.
2. **Two variants without a recoverable standard error (data).** rs4942486 (CHD) and rs7422339 (triglycerides) have a beta of exactly 0 with
   P = 1 in `dodata.csv`, which gives no standard error. This is a limit of the deposited file, which has no standard errors.
3. **Different outcome release (data, confirmed).** Part A differs from the paper because the paper used the 2015 CHD release. Confirmed in
   two ways: on the same variants the univariable IVW estimate for HDL-c is -0.256 with `dodata.csv` (Bowden et al. 2016 report -0.254 on
   this file) against -0.130 in Rees Table 1, a gap of the same size as in Part A with no multivariable code involved; and with the 2015
   release (Part B) the gap shrinks to at most 0.016, with univariable IVW for HDL-c -0.135 (-0.232, -0.039) against the paper's -0.130
   (-0.227, -0.033).
4. **Three variants absent from the 2015 release on OpenGWAS (data, quantified).** rs1998013, rs2247056 and rs7422339 are not in `ieu-a-7`,
   with or without proxies, so Part B uses 182 of 185 variants. Leaving three variants out at random changes MVMR-IVW by up to 0.014, 0.015
   and 0.022 (95th percentile of 2,000 draws), which covers the remaining differences. Criteria for Part B (all met): estimates within 0.025,
   intercepts within 0.001, identical replay, equality with the package. Matching to three decimals would need the original release file
   with all 185 variants.
