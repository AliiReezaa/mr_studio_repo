# MR GUI — Solutions for Public-GWAS Data Limits

Version 1 · 2026-10-25 · companion to `mr_methods_inventory.md` (v3). Each solution is something the app does automatically; the user only sees a badge, a warning or a recommendation. Everything here is free and runs locally — no AI keys.

Status labels: **[verified]** = function or behaviour checked in package source / documentation this session; **[standard]** = established statistical practice with the cited paper; **[design]** = our engineering choice.

---

## 1. Ancestry metadata that is coarse or missing

**Problem.** Population labels in OpenGWAS/GWAS Catalog are sometimes broad ("European"), missing, or describe multi-ancestry meta-analyses. A wrong label silently breaks LD-dependent methods and creates exposure–outcome population mismatch.

**Solution — "declared + inferred" ancestry with a confidence badge.**
1. **Declared ancestry**, in priority order: GWAS Catalog curated ancestry fields (broad ancestral category per study stage; Sollis 2023, PMID 36350656) → OpenGWAS `population` field → a small manual-override table that you and your friend curate. **[design]**
2. **Inferred ancestry from the data itself**: `bigsnpr::snp_ancestry_summary()` estimates ancestry proportions of any GWAS from its allele frequencies against 21 UK Biobank–derived reference groups (Privé 2022, PMID 35604078, [doi](https://doi.org/10.1093/bioinformatics/btac348)). **[verified: function exported in bigsnpr NAMESPACE]** Requires effect-allele frequency (EAF); if EAF is absent, the badge shows "declared only".
3. **Classification shown in the GUI:** *single-ancestry* (e.g., EUR 97%), *multi-ancestry meta-analysis*, or *unknown*, plus a flag when declared and inferred labels disagree. The ≥90%-one-group cut-off for "single-ancestry" is a **[design] heuristic** to tune.
4. **Consequences in the feasibility engine:** multi-ancestry or unknown datasets are blocked for methods with `needs_ld ≠ none` (no matching LD panel exists) and allowed with a warning for clumped-IV methods; exposure/outcome mismatch triggers rule R13 and offers XMR / TEMR / MR-EILLS when ≥2 ancestries are available.
5. **Matching LD reference:** 1000 Genomes Phase 3 super-populations (EUR, EAS, AFR, AMR, SAS; PMID 26432245) for clumping and LD matrices; Pan-UKB LD scores in LDSC-compatible format (`.l2.ldscore.gz`, `.M_5_50`) for six UK Biobank ancestry groups ([Pan-UKB LD docs](https://pan.ukbb.broadinstitute.org/docs/ld)) for LDSC-based methods. **[verified: documentation page]**
6. **Regional QC for cis-MR / colocalisation:** run DENTIST to detect and remove GWAS–LD-reference inconsistencies before LD-matrix methods (Chen 2021, PMID 34880243, [doi](https://doi.org/10.1038/s41467-021-27438-7)). **[standard]**

## 2. Units and scaling differ between datasets

**Problem.** One GWAS reports per-SD effects, another natural units, another log-odds or odds ratios; MVMR and result interpretation break when units are mixed or unknown.

**Solution — a "units passport" for every dataset, filled automatically and editable.**
1. **Record per dataset:** trait type (continuous / binary / count), effect type (beta, log-OR, OR, HR), scale (SD, natural unit, log-transformed, inverse-normal-transformed), N, ncase/ncontrol. **[design]**
2. **Auto-detect effect type:** if all effects are positive and centred near 1, treat as OR and convert to log-OR; confirm by checking that |β/SE| reproduces the reported P-value. **[design]**
3. **Detect whether a continuous trait is in SD units:** the implied trait SD ≈ median over SNPs of SE·√(2p(1−p)N) — values near 1 indicate SD units. TwoSampleMR exposes `estimate_trait_sd()` and `standardise_units()` for this. **[verified: both exported in TwoSampleMR NAMESPACE]**
4. **Standardise when needed:** β_std = z / √(2p(1−p)(N + z²)), SE_std = 1 / √(2p(1−p)(N + z²)) — the standard per-SD transformation used in SMR (Zhu 2016, PMID 27019110). **[standard]** Requires EAF and N; otherwise keep original units and display them.
5. **Binary exposures:** report the estimate per unit increase in log-odds and also per doubling of the odds (×ln 2 ≈ 0.693), with the caveat that the causal interpretation relates to the underlying continuous liability (Burgess & Labrecque 2018, PMID 30039250, [doi](https://doi.org/10.1007/s10654-018-0424-6)). **[standard]** → rule R30.
6. **GUI output always states units**, e.g. "OR 1.45 per 1-SD higher LDL-C". MVMR with exposures of unknown units is allowed but flagged amber. **[design]**

## 3. Mixed genome builds (GRCh37 vs GRCh38) and allele inconsistencies

**Problem.** Datasets use different builds, some lack rsIDs, and strand/allele coding differs.

**Solution — one internal standard, harmonised on entry.**
1. **Primary key = rsID + alleles**, as TwoSampleMR harmonisation does; positions are secondary. **[standard]**
2. **Infer the build** with `MungeSumstats::get_genome_builds()` and convert with `MungeSumstats::liftover()` (rtracklayer chain files) to one canonical build. **[verified: both exported in MungeSumstats NAMESPACE]** MungeSumstats also standardises columns and performs summary-statistics QC (Murphy 2021, PMID 34601555).
3. **Canonical build = GRCh37** for now, because 1000 Genomes Phase 3 and the common LDSC reference files are GRCh37; store GRCh38 coordinates alongside so a later switch is cheap. **[design]**
4. **Missing rsIDs:** map chr:pos:ref:alt to dbSNP rsIDs during munging (MungeSumstats). **[standard]**
5. **Allele/strand QC:** compare dataset EAF with the matched reference-panel frequency; flag variants with large differences as possible allele flips (|ΔEAF| > 0.2 is a **[design] heuristic**). Palindromic variants are resolved with EAF via `harmonise_data(action = 2)` and dropped when ambiguous (rule R29). **[standard]**
6. **Internal storage format:** GWAS-VCF (Lyon 2021, PMID 33441155) or an equivalent Parquet table with the same fields, so every method reads identical, harmonised data. **[design]**

## 4. OpenGWAS usage allowance can run out

**Problem.** Since 1 May 2024 every OpenGWAS API call needs a personal token, and each request is deducted from a free, periodically replenished allowance ([OpenGWAS announcement](https://blog.opengwas.io/posts/user-auth-spring-2024/)). Web downloads of full VCF files are limited to 20 datasets per 24 h, with links valid for 2 h ([example dataset page](https://opengwas.io/datasets/ieu-a-2)). **[verified: official pages]**

**Solution — "offline-first" data layer.**
1. **Local catalogue:** download the dataset metadata once and store it locally (SQLite/DuckDB); trait search and population filtering never call the API. Refresh weekly/monthly. **[design]**
2. **Cache every response and every file**, keyed by dataset ID + retrieval date; never request the same data twice. Full summary statistics are stored locally after first download; a download queue spreads large requests over days to respect the 20-per-day limit. **[design]**
3. **Local clumping and LD:** `ieugwasr::ld_clump()` runs locally when given a PLINK reference (`bfile`) and binary (`plink_bin`) instead of the API. **[verified: ieugwasr source]** LD matrices are computed locally with PLINK on 1000 Genomes. This removes the heaviest API usage.
4. **Polite API use:** batch SNP queries, exponential back-off on rate-limit responses, and an allowance indicator in the GUI. **[design]**
5. **Alternative free sources** for the same traits: GWAS Catalog summary-statistics downloads, FinnGen releases, Pan-UKB files. The catalogue records which sources hold each trait, so the app can fall back automatically. **[design]**
6. **Security:** each user keeps their own token, stored in the operating-system credential store (macOS Keychain / Windows Credential Manager, e.g. via the R `keyring` package) — never in code, config files or logs. **[design]**

## 5. Sample overlap can only be estimated

**Problem.** Exposure and outcome GWAS often share participants (UK Biobank especially); overlap with weak instruments biases two-sample MR toward the observational association (Burgess 2016, PMID 27625185).

**Solution — detect, quantify, avoid, then correct.**
1. **Cohort registry:** store the contributing cohorts/biobanks of each dataset; flag "likely overlap" when both contain the same cohort (e.g., both include UK Biobank). **[design]**
2. **Quantify from data:** the cross-trait LDSC intercept reflects sample overlap combined with phenotypic correlation (Bulik-Sullivan 2015, PMID 26414676). A non-zero intercept → overlap likely. These same LDSC quantities are the inputs that MR-APSS, MRlap, MRMU, XMR and METEOR use to model overlap. **[standard]**
3. **Show expected bias:** bias increases with the overlap fraction and with instrument weakness (lower F), so strong instruments keep it small (Burgess 2016). The GUI shows the mean F-statistic next to the overlap flag. **[standard]**
4. **Avoid first:** the app proposes non-overlapping alternatives automatically (e.g., UK Biobank exposure with FinnGen or consortium outcome; three-sample designs where an independent selection GWAS exists). **[design]**
5. **Correct if unavoidable:** rule R28 recommends overlap-aware methods — MRlap (PMID 37036286), MR-APSS (PMID 35787050), MR-SimSS (PMID 42102131), CAUSE (PMID 32451458) — and demotes methods with `overlap_handling = assumes_none` to amber.
6. **Report:** the overlap estimate, its source (registry vs LDSC) and the chosen correction are written into the results report. **[design]**

## 6. Messy trait names (needed for the exposure/outcome pickers)

**Problem.** The same trait appears under many labels across sources.

**Solution — ontology mapping without AI.**
1. Use EFO terms already attached to GWAS Catalog studies. **[standard]**
2. Map remaining free-text OpenGWAS trait names to EFO with EMBL-EBI's free ZOOMA / Ontology Lookup Service, then cache the mapping. **[design]**
3. Fuzzy string matching as fallback, plus a manual override table you curate.
4. The picker groups datasets by EFO term and shows, per dataset: source, ancestry badge, N / cases, build, units passport and full-summary-statistics availability — exactly the fields the feasibility engine reads.

---

### Sources
PubMed records: Sollis 2023 [doi](https://doi.org/10.1093/nar/gkac1010); Privé 2022 [doi](https://doi.org/10.1093/bioinformatics/btac348); 1000 Genomes 2015 [doi](https://doi.org/10.1038/nature15393); Chen 2021 [doi](https://doi.org/10.1038/s41467-021-27438-7); Zhu 2016 [doi](https://doi.org/10.1038/ng.3538); Burgess & Labrecque 2018 [doi](https://doi.org/10.1007/s10654-018-0424-6); Murphy 2021 [doi](https://doi.org/10.1093/bioinformatics/btab665); Lyon 2021 [doi](https://doi.org/10.1186/s13059-020-02248-0); Burgess 2016 [doi](https://doi.org/10.1002/gepi.21998); Bulik-Sullivan 2015 [doi](https://doi.org/10.1038/ng.3406); Mounier 2023 [doi](https://doi.org/10.1002/gepi.22522); Hu 2022 [doi](https://doi.org/10.1073/pnas.2106858119); Forde 2026 [doi](https://doi.org/10.1371/journal.pgen.1011949); Morrison 2020 [doi](https://doi.org/10.1038/s41588-020-0631-4).
Web: [OpenGWAS authentication announcement](https://blog.opengwas.io/posts/user-auth-spring-2024/) · [OpenGWAS dataset page](https://opengwas.io/datasets/ieu-a-2) · [Pan-UKB LD scores](https://pan.ukbb.broadinstitute.org/docs/ld). Package checks: bigsnpr, MungeSumstats, TwoSampleMR and ieugwasr source on GitHub (2026-10-25).
