# Import of GWAS files from real sources

**Question.** When users bring their own files instead of OpenGWAS downloads, does MR Studio read every common layout correctly,
with as little help from the user as possible, and give the same MR result as the same data from OpenGWAS?

**Data.** Two datasets downloaded from OpenGWAS inside MR Studio (Data library): the LDL cholesterol instruments (ieu-a-300, 79
variants) and coronary heart disease (ieu-a-7) looked up at them. The reference pair (LDL on CHD) gives IVW 0.412377 (SE 0.023341,
77 variants after harmonisation).

**Corpus** (`make_corpus.R`). The same two datasets are written in 15 layouts that summary statistics come in, each as an exposure
file and an outcome file (30 files), following each source's documentation:

| # | Layout | What it tests |
|---|---|---|
| 01 | GWAS Catalog GWAS-SSF | GRCh38 positions; `rsid` and `variant_id` columns |
| 02 | GWAS Catalog harmonised | `hm_` prefixes; `hm_beta` and an empty `hm_odds_ratio`; GRCh38 |
| 03 | FinnGen | GRCh38; effect allele `alt`; cells with two rsIDs; `mlogp` |
| 04 | Neale lab UK Biobank | `.tsv.bgz`; only `variant` = chr:pos:ref:alt (no rsID); `minor_allele` |
| 05 | BOLT-LMM | `ALLELE1`/`ALLELE0`; two P columns |
| 06 | REGENIE | space separated; `LOG10P` only; `ALLELE0`/`ALLELE1`; `TEST` |
| 07 | PLINK 2 `--glm` | `A1` is REF in half the rows; covariate rows (`TEST` = AGE, SEX); `OR` and `LOG(OR)_SE` |
| 08 | SAIGE | the effect allele is `Allele2` (`AF_Allele2`) |
| 09 | METAL | lower-case alleles; Windows line endings |
| 10 | GIANT style | quoted CSV; `Tested_Allele`; a frequency heading naming the cohort |
| 11 | A paper's supplementary table (Excel) | title and note rows, footnote; `Effect/other allele` = "T/C"; `Beta (SE)`; `OR (95% CI)` = "0.992 (0.974–1.010)"; `Chr:position (hg19)`; P printed as "1.2 × 10−8" |
| 12 | European spreadsheet export | semicolons and decimal commas |
| 13 | Messy in-house text | `##` preamble; missing-value tokens; a missing SE; a duplicated row; an indel; `chr1`; Windows line endings |
| 14 | LDSC munged | Z scores and N only (exposure: frequencies from the 1000 Genomes panel; outcome: a binary trait, must be refused) |
| 15 | Estimate with 95% CI | no SE column |

**Check** (`run_import_check.R`). Each file is imported as a user does it: the app's column suggestions and pre-filled answers are
accepted, and only what a careful user knows from the source's documentation is answered (type of trait, cases and controls, role,
the genome build when the file has positions and the app could not detect it). Every answer that had to differ from the app's proposal
is counted (`user_changes`). Each imported exposure is paired with the reference outcome, each imported outcome with the reference
exposure, and IVW is compared with the reference pair. (The reference outcome is tied to the reference exposure's instruments; in the
temporary copy used here that tie is released, since the imported exposures hold the same variants.)

    Rscript validation/import_formats/make_corpus.R "<Data library>/LDL" "<Data library>/CHD for LDL v2" <corpus folder>
    Rscript validation/import_formats/run_import_check.R <corpus folder> "<Data library>/LDL" "<Data library>/CHD for LDL v2"

The 1000 Genomes panel (Settings, LD reference panel) must be installed for the build detection, the rsIDs of the Neale files and the
frequencies of the LDSC file.

**What the first run found** (each fixed in Round 79, each with an engine test):

- GWAS Catalog harmonised: `hm_` headings were not recognised and the empty odds-ratio column was taken as the effect; the outcome
  imported **without error but with odds ratios read as log odds ratios** (IVW 2.17 instead of 0.41).
- FinnGen: rows with two rsIDs in one cell were lost (70 of 77 variants).
- Neale lab: `.bgz` files were refused, and a file with only chr:pos:ref:alt IDs could not be paired with rsID datasets.
- PLINK 2: where A1 is REF, the other allele became the effect allele and the row was dropped (40 of 77 variants).
- SAIGE: Allele1 was proposed as the effect allele (it is Allele2).
- Papers' tables: combined columns were not split, so nothing could be imported.
- 95% CI without SE: the interval was not used, and an interval column was mistaken for the frequency.
- A single missing SE dropped the row instead of being derived from its P value.
- Every file with positions needed the genome build typed in.

**Result** (`results.csv`, run of 2026-11-21): see the table written by the script. All 29 importable files import; the LDSC outcome
(a binary trait from Z scores) is refused with the reason, as it must be. 23 files reproduce the reference IVW exactly (to 1e-8). The
others differ for stated reasons: no allele frequency in the file, so one palindromic variant cannot be aligned and is left out (Neale,
CI-only files); a confidence interval printed to three decimals (the paper's OR table); an effect derived from Z, N and panel
frequencies, per standard deviation (LDSC exposure); one Neale outcome variant not in the 1000 Genomes panel. With the fixes, the
user only answers what no file states: the type of trait, the numbers of cases and controls, and the build of position-only files
(Neale) or of files without rsIDs in the panel.
