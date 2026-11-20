# Validation study 17: MR-APSS, BMI and standing height in UK Biobank

Reference: Hu X, Zhao J, Lin Z, et al. Mendelian randomization for causal inference accounting for pleiotropy and sample structure
using genome-wide summary statistics. PNAS 2022;119:e2106858119 (doi:10.1073/pnas.2106858119, PMID 35787050). Package MRAPSS (GitHub
YangLabHKUST/MR-APSS, MIT), commit 85c18ea.

Published values: the authors' result file `Traits_MRAPSS.MRres` (repository YangLabHKUST/MRAPSS_RealDataAnalysis_reproduce), rows
Height_UKB -> BMI and BMI -> Height_UKB, from the UK Biobank summary statistics of Watanabe et al. 2019 (about 385,000 people).

Data here: the Neale lab UK Biobank summary statistics for BMI and standing height (about 360,000 people) in LDSC format, which holds
the HapMap3 variants that MR-APSS uses; effects Z / sqrt(N) with standard error 1 / sqrt(N). Background from the LDSC authors' European
LD scores; instruments clumped on this computer with PLINK and the 1000 Genomes EUR panel. Each download is checked by MD5.

    MRSTUDIO_VALIDATE_MRAPSS=1 MRSTUDIO_VALIDATION_CACHE=<folder> Rscript --vanilla validation/hu2022_mrapss/run_validation.R

About 30 minutes on one core.

## Result: PASSED (`results.csv`)

| Direction | Paper (Watanabe UKB) | MR Studio (Neale UKB) | Distance |
|---|---|---|---|
| Height on BMI | -0.0637 (SE 0.0135), 2,226 instruments | -0.0761 (SE 0.0104), 1,837 instruments | 0.7 combined SE |
| BMI on height | -0.0330 (SE 0.0602), 1,300 instruments | -0.1396 (SE 0.0478), 1,307 instruments | 1.4 combined SE |

Both directions agree with the paper within two combined standard errors. Height on BMI is clearly negative in both. For BMI on
height the sign agrees but the evidence differs: the paper's 95% interval includes zero (-0.151 to 0.085), MR Studio's does not
(-0.233 to -0.046); the two intervals overlap. The two analyses use different UK Biobank GWAS (Watanabe et al. against the Neale lab, with
different covariates and about 25,000 fewer people here) and therefore different instrument sets, so exact agreement is not expected.

On the same prepared files with the same seed, MR Studio equals the package pipeline called directly (`format_data`, `est_paras`,
`clump`, `MRAPSS`) to 3e-17 in both directions, and the study script replays identically.

## The package's bundled example

The package ships a BMI and type 2 diabetes example with precomputed background parameters (1,227 instruments). The pinned version
gives 0.216 (SE 0.026); version 1.0.0, which the paper used, gives 0.209 (SE 0.033) on the same example. The paper's own BMI and type
2 diabetes estimate, 0.326 (SE 0.056) with 1,296 instruments, came from the authors' formatted data with a different instrument set,
which is not distributed with the package, so it is not compared here.

## How MR Studio runs it

All steps run on the user's computer: the summary statistics of both traits are restricted to HapMap3 variants with the MHC region
removed (`format_data`), the background (LDSC intercepts for sample structure and the polygenic-pleiotropy covariance) is estimated
from them (`est_paras`), instruments are selected at P < 5e-5 and clumped at r2 < 0.001 within 1,000 kb (`clump`, local PLINK), and
MR-APSS corrects for the selection bias of that threshold. No OpenGWAS query is needed.
