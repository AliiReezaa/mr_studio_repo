# Validation study 2: Hemani et al. 2018, the MR-Base applied example

Status: **run on real OpenGWAS data on 2026-11-05 and again on 2026-11-06 (Round 28).** Tier A VERIFIED for all three clumping choices and for the tie sensitivity analysis. Tier B: IVW agrees with the paper.

Paper: Hemani G, Zheng J, Elsworth B, et al. *The MR-Base platform supports systematic causal inference across the human
phenome.* eLife 2018;7:e34408, doi:10.7554/eLife.34408 (PMID 29846171, open access PMC5976434). The applied example is the
effect of LDL cholesterol (GLGC) on coronary heart disease (CARDIoGRAMplusC4D), run on MR-Base, now called OpenGWAS.
Instruments: P < 5e-8, clumping r² 0.001, window 10,000 kb (the MR Studio defaults).

Reported in the text (odds ratio of CHD per SD of LDL cholesterol):

| Method | Paper |
|---|---|
| IVW, random effects | 1.45 (1.30 to 1.62) |
| MR-Egger | 1.85 (1.48 to 2.32); intercept -0.013 (SE 0.005), p = 0.020 |
| Weighted median | 1.56 (1.43 to 1.70) |
| Weighted mode | 1.68 (1.56 to 1.80) |
| Cochran's Q | 122.5, p = 4.72e-07 |

The number of variants is not printed in the text. Q = 122.5 with p = 4.72e-07 fits only 55 degrees of freedom, so about 56
variants were used (inferred). The paper names three outliers (rs11065987, rs1250229, rs4530754) and says none is palindromic.

## What is compared

- **Tier A, exact.** The same OpenGWAS datasets through MR Studio's engine and through TwoSampleMR and ieugwasr directly, for each
  clumping choice: stored clumps (against `TwoSampleMR::extract_instruments` with its defaults), live server clumping (against
  `ieugwasr::tophits(force_server = TRUE)`) and this computer (against `ieugwasr::ld_clump_local` on the same EUR panel and PLINK).
  Then the outcome lookup with proxies (r² 0.8), harmonisation (action 2, allele frequencies available), IVW (multiplicative random
  effects; TwoSampleMR mr_ivw since Round 30, equal to mr_ivw_mre here because Q/df is 4.0 to 5.2 for every clumping choice), MR-Egger, weighted median, simple median and weighted mode, Cochran's Q, and the replay check.
  Bootstrap SEs (median and mode methods) are expected to differ slightly, because the seeds differ; estimates should be identical.
- **Tier B, approximate.** MR Studio against the paper. Exact agreement is not expected: the paper used 2018 data and clumps, and
  OpenGWAS may have changed since.

## Dataset IDs

The paper names the consortia, not the OpenGWAS IDs. The defaults are `ieu-a-300` (LDL cholesterol, GLGC) and `ieu-a-7` (coronary
heart disease, CARDIoGRAMplusC4D); they are **not verified** as the paper's datasets. The script prints the catalogue record of both
(consortium, author, year, sample size, cases and controls) for comparison with the paper (GLGC 214,370 subjects,
CARDIoGRAMplusC4D 195,813 subjects). Other IDs can be given with `MRV_EXPOSURE` and `MRV_OUTCOME`.

## How to run

From the repository root, in a Terminal where the token is set for that session only (never in `.Renviron` or any file):

    export OPENGWAS_JWT='paste the token here'
    Rscript --vanilla validation/hemani2018/run_validation.R
    unset OPENGWAS_JWT

The token is read from the environment and never printed or written (checked: a rejected token gives the OpenGWAS 401 message,
and the token does not appear in the output). Everything downloaded goes to a throw-away folder. Only `results.csv` (estimates and
counts) and `instruments.csv` (rsIDs kept by each clumping choice, no effect sizes) are written here. Local clumping uses the EUR
panel installed in MR Studio (or `MRV_LD_DIR`) and the PLINK in `prototype/engine/tools/plink`.

`MRV_DRYRUN=1` runs the whole script on simulated effects at real 1000 Genomes variants with a fake OpenGWAS (no token, no
internet) and writes its tables to a temporary folder, not here. On 2026-11-05 the dry run gave identical instruments for all three
choices (78), identical harmonised data (73 kept, 2 by proxy, 4 palindromic dropped), estimates within 1.2e-15 of TwoSampleMR and
identical replays. Since Round 28 it gives six triplets of neighbouring variants P = 1e-200, so ties decide a lead; on 2026-11-06 the
tie sensitivity analysis also equalled its independent reference (within 1.5e-15).

The tie sensitivity analysis (`local_ties`) is checked against an independent reference: the candidates ranked by P and then by
the larger |beta/SE|, written by the script itself, and clumped with `ieugwasr::ld_clump_local`.

## Result (2026-11-06, Round 28; ieugwasr 1.2.0.9000, TwoSampleMR 0.7.13)

Datasets: `ieu-a-300` LDL cholesterol (GLGC, Willer 2013, N 173,082) and `ieu-a-7` CHD (CARDIoGRAMplusC4D, Nikpay 2015, 60,801
cases, 123,504 controls). Neither matches the sample sizes in the paper, so the paper may have used other datasets or versions.
The catalogue gives "Mixed" ancestry for both; MR Studio then refuses to choose a panel itself, so the script sets EUR (the
TwoSampleMR default).

**Tier A, MR Studio against TwoSampleMR / ieugwasr on the same downloads**

| Clumping | Instruments (app / reference) | Harmonised data | Estimates (all 5 methods) | Q | Replays |
|---|---|---|---|---|---|
| Stored clumps | 75 / 75, identical | identical (73 kept, 1 proxy, 2 palindromic dropped) | within 4.4e-16 | 285.5 on 72 df, identical | identical |
| Live server | 79 / 79, identical | identical (77 kept, 1 proxy, 2 palindromic dropped) | within 2.8e-16 | 363.8 on 76 df, identical | identical |
| This computer | 79 / 79, identical | identical (77 kept, 1 proxy, 2 palindromic dropped) | within 2.8e-16 | 363.8 on 76 df, identical | identical |
| Tie sensitivity (this computer, larger \|beta/SE\| among equal P) | 76 / 76, identical | identical (74 kept, 1 proxy, 2 palindromic dropped) | within 5.0e-16 | 379.9 on 73 df, identical | identical |

Bootstrap SEs of the median and mode methods differ in the third decimal (different seeds), as expected.
TwoSampleMR names a proxy row by the instrument it replaces (checked; the dry run had assumed it).

Tied P values: OpenGWAS stores the strongest GLGC P values as 1e-200, so 59 of the 3,023 candidates repeat another's P value
(102 candidates share a P value with at least one other). On 2026-11-05 MR Studio still broke these ties by |beta/SE| in its
primary analysis (Round 25), which gave 76 instruments instead of the server's 79. Since Round 28 the primary analysis keeps
PLINK's own choice among equal P values (the variant first in the reference panel), and local clumping equals
`ieugwasr::ld_clump_local` and live server clumping exactly (79 of 79). The tie-break is now the tie sensitivity analysis: it keeps
rs12740374 (SORT1, z = -36.6) instead of rs4970834 (z = -32.0) and rs7412 (APOE, z = -58.4) instead of rs7254892 (z = -40.8), and the
lead change frees or removes neighbours (rs8103315 in; rs413380, rs2965157, rs12721109, rs676388 out): 6 primary instruments out,
3 in, 76 instruments. The weighted mode is the most sensitive to this (0.255 against 0.540); IVW moves from 0.412 to 0.395.
Stored server clumps (75) cannot be checked for ties: the server does not say which tied variants it left out.

**Tier B, MR Studio against the paper** (odds ratio per SD of LDL cholesterol)

| Method | Paper (2018) | Stored | Live and this computer | Tie sensitivity |
|---|---|---|---|---|
| IVW | 1.45 (1.30 to 1.62) | 1.47 (1.36 to 1.60) | 1.51 (1.37 to 1.67) | 1.48 (1.35 to 1.63) |
| MR-Egger | 1.85 (1.48 to 2.32) | 1.49 (1.33 to 1.68) | 1.65 (1.42 to 1.93) | 1.55 (1.36 to 1.77) |
| Egger intercept | -0.013 (0.005) | -0.001 (0.004) | -0.007 (0.005) | -0.004 (0.004) |
| Weighted median | 1.56 (1.43 to 1.70) | 1.41 (1.29 to 1.54) | 1.49 (1.36 to 1.62) | 1.38 (1.27 to 1.50) |
| Weighted mode | 1.68 (1.56 to 1.80) | 1.29 (1.20 to 1.40) | 1.72 (1.47 to 2.01) | 1.29 (1.19 to 1.40) |
| Variants; Q | about 56; 122.5 | 73; 285.5 | 77; 363.8 | 74; 379.9 |

Live and this computer give the same estimates; their bootstrap intervals differ slightly (weighted mode on this computer: 1.72, 1.49 to 1.98).

The IVW estimate agrees with the paper with every clumping choice, and every method gives a harmful effect of LDL cholesterol.
The paper used fewer variants and found less heterogeneity and a clearer negative MR-Egger intercept; with today's datasets and
clumps these do not reproduce. The paper's outliers are among the instruments: all three with stored clumps, rs1250229 and
rs4530754 with live and local clumping and in the tie sensitivity analysis.

`results.csv` has every number above; `instruments.csv` lists which rsIDs each clumping choice kept.

## Attempt at a matched recent paper (2026-11-05)

To compare with a published study that used the same data and settings, Europe PMC full texts were searched for recent papers
naming `ieu-a-300` and `ieu-a-7`. The closest match was Bure et al. 2026, *Medicine* 105:e47880 (doi:10.1097/MD.0000000000047880,
PMID 41894294): LDL cholesterol (`ieu-a-300`) on CAD (`ieu-a-7`) with TwoSampleMR, P < 5e-8, r² < 0.001, a window given as
"10,000 base-pair", 113 variants, IVW 0.678496 (SE 0.042268), MR-Egger 0.410885 (SE 0.063061).
Run today through TwoSampleMR itself (0.7.13), with no MR Studio involved:

| Window | Clumping | Variants kept | IVW (SE) | MR-Egger (SE) |
|---|---|---|---|---|
| 10 kb | stored or live | 584 | 0.391 (0.019) | 0.483 (0.034) |
| 10,000 kb | stored | 73 | 0.388 (0.042) | 0.400 (0.060) |
| 10,000 kb | live | 77 | 0.412 (0.051) | 0.503 (0.078) |

Neither reading reproduces the paper's 113 variants or its IVW estimate, so the paper cannot be reproduced from its stated methods
even with the package it used; it was not used as a reference. A second paper, Lu et al. 2025 (MRanalysis, GigaScience,
PMC12616851), gives an IVW-type estimate of 0.4110 for the same pair (TwoSampleMR 0.6.8), close to today's live-clumping 0.4124 but
without the variant count or settings needed for an exact comparison.
