# Demo files for testing the import wizard

All numbers in these files are simulated (random effects, random rsIDs). They are for testing the import screens only, not for any analysis. Each file has 300 rows unless stated, so in step 3 (the questions about the file) choose **Only selected variants** for the question about which variants the file holds, otherwise you will get the normal "you said genome-wide but the file has only 300 rows" warning.

How to use: Data library, Import a GWAS file, pick the file, then follow the four steps. Checked with the app's own import code; the outcomes below are what it produces.

## Files that should import cleanly

| File | What it tests | Expected |
|---|---|---|
| 01_clean_generic.csv | Plain headings (rsid, effect_allele, other_allele, eaf, beta, se, pval, n) | All 8 columns matched as HIGH. 300 rows read, 300 kept. 100 palindromic variants (they have a frequency, so no warning). |
| 02_bolt_lmm_style.tsv | BOLT-LMM style headings, tab separated, with chromosome and position | 9 columns matched. ALLELE0 (other allele) is MEDIUM, confirm it. Warning: no sample size, so type a total N in step 3 to remove it. |
| 02b_bolt_lmm_style_second_trait.tsv | Same headings as 02, different numbers | Import 02 first and tick "Remember this column mapping". On 02b the wizard should say TEMPLATE FOUND and fill the columns for you. |
| 03_gwas_catalog_style.tsv | GWAS Catalog style headings (variant_id, base_pair_location, standard_error, p_value) | All matched; position (base_pair_location) is MEDIUM. |
| 04_regenie_style.tsv | REGENIE style: ALLELE1 is the effect allele, LOG10P is not a P value | LOG10P is left unmatched on purpose (it is -log10 of P). No P column is fine. |
| 12_space_separated.txt | Space separated text | Separator detected, all matched. |
| 13_compressed_example.tsv.gz | Gzip-compressed file with METAL-style headings (MarkerName, Allele1, Freq1, Effect, StdErr, P-value) | Read without unzipping. All matched. |
| 18_misspelt_headings.csv | Headings with typos (effect_alelle, std_erorr, betaa) | Seven MEDIUM suggestions (the app says "looks similar to a known heading"). Check each, they are right. |
| 20_no_rsid_chr_pos_only.tsv | No rsID column, only chromosome and position | Allowed: chromosome plus position replace the variant ID. 300 rows kept. |

## Files that import with warnings

| File | What it tests | Expected |
|---|---|---|
| 05_odds_ratio_case_control.csv | Odds ratios | In step 3 choose **Odds ratio** (the SE column must be the SE of the log odds ratio, as the file says). Imports; values are saved as log odds. 100 palindromic variants and no frequency, so you get that warning. |
| 06_messy_values_warnings.csv | Bad values on purpose | Warnings: 1 effect size not a number, 2 variants with missing or non-positive SE, 1 P value outside 0 to 1, 2 insertions or deletions dropped, 3 repeated rows. 303 rows read, 293 kept. |
| 07_minor_allele_frequency.csv | Frequency column is a minor allele frequency | In step 3 choose **Minor allele frequency**. Warning that it will not be used as an effect allele frequency, plus the palindromic warning. |
| 19_palindromic_variants_no_frequency.csv | 120 A/T and C/G variants, no frequency | Warning that palindromic variants cannot be aligned without a frequency and will be left out at harmonisation. 120 rows kept. |

## Files the wizard should stop

| File | What it tests | Expected |
|---|---|---|
| 05b_log_odds_chosen_as_odds_ratio.csv | A column that holds log odds, but you choose "Odds ratio" | Error: odds ratios must be positive. Choose "Beta or log odds" and it imports. |
| 08_missing_standard_error.tsv | No standard error column | Error: no column chosen for "Standard error". It is required. |
| 09_individual_level_refused.tsv | One row per person (IID, FID, age, genotype columns) | Refused when the file is opened: individual-level data is not accepted by design. Nothing is read. |
| 10_zscore_only.tsv | Only Z scores | Errors: no effect size and no standard error. Z-scores are not supported yet. |
| 11_wrong_allele_column.csv | The "effect allele" column holds rsIDs | Error: more than half of the rows do not have single-letter alleles. |
| 14_unsupported_gwas_vcf.vcf | GWAS-VCF | Refused: file type not supported yet. |
| 15_too_few_columns.csv | Only 3 columns | Refused: fewer than 4 columns. |
| 16_duplicate_heading.csv | Two columns called BETA | Refused: two columns have the same heading. |

## File that needs your help

| File | What it tests | Expected |
|---|---|---|
| 17_unlabelled_columns.csv | Headings are col1 to col7 | Only variant, the two alleles and the frequency are guessed from the values, marked LOW ("please confirm"). Effect size (col5), standard error (col6) and P value (col7) must be chosen by hand from the drop-downs, and the wizard will not let you continue until you do. |

## Not covered yet

Decimal commas (for example 0,123), files above a few hundred MB, genome build detection and Z-score import are not supported. Imported datasets are saved to the Data library with their recipe, but cannot be used in a study yet because clumping and harmonisation are not built.

## A matching pair for a full study (files 23 and 24)
- `23_pair_exposure_clumped.csv`: 60 independent variants (simulated). Import it as an **exposure**, answer "Yes, clumped and independent", frequency column = effect allele frequency, "Only selected variants".
- `24_pair_outcome_full.csv`: the same variants in an outcome study (simulated, true effect 0.25), with 3 variants missing, 12 with swapped alleles and some on the opposite strand. Import it as an **outcome**.
- Then start a study, pick both on the Data screen and run. The Check screen lists what happened to every instrument. Expected: 60 variants, 3 not found, about 2 palindromic dropped, 13 flipped, about 50 kept, IVW close to 0.23.
