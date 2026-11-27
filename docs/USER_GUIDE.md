# MR Studio user guide

MR Studio runs Mendelian randomization (MR) analyses with GWAS summary statistics on your own computer. You choose a question,
pick or import the data, and the app checks the data, offers only the methods the data can support, runs them, and builds a
publication package. Every analysis is written as an R script that you can rerun yourself.

This guide follows the seven steps of a study. For a worked example from start to finish, see the
[tutorial](TUTORIAL_LDL_CHD.md).

## Contents

1. [Install](#install)
2. [Before your first study](#before-your-first-study)
3. [The Data library](#the-data-library)
4. [A study in seven steps](#a-study-in-seven-steps)
5. [Analyses that need full summary statistics](#analyses-that-need-full-summary-statistics)
6. [Checking the results in R](#checking-the-results-in-r)
7. [Where things are kept](#where-things-are-kept)
8. [Troubleshooting](#troubleshooting)

## Install

Download the installer for your system from the project's Releases page on GitHub:

- **macOS**: the `.dmg` for Apple Silicon (M1 or later) or for Intel. Open it and drag MR Studio to Applications. The app is not yet
  signed with an Apple developer certificate, so the first time macOS refuses to open it. Open **System Settings > Privacy & Security**,
  find the message about MR Studio near the bottom and choose **Open Anyway**. This is needed only once.
- **Windows**: the `.exe` installer. If Windows SmartScreen stops it, choose **More info**, then **Run anyway**.

MR Studio includes everything it needs: R, the R packages of every method and PLINK. Nothing else has to be installed.

This guide describes MR Studio 0.9.0. The version is shown in Settings and in Run diagnostics, and is written into every report,
reproducibility capsule and STROBE-MR checklist.

## Before your first study

Open **Settings > Data and access**.

- **OpenGWAS access token.** Needed to search and download public GWAS from OpenGWAS. It is free: sign in at api.opengwas.io and
  copy your token. MR Studio keeps it in memory only, never in a file, so you paste it again after closing the app.
- **OpenGWAS catalogue.** The list of every public dataset. It is downloaded once and kept on your computer, so searching is instant.
- **1000 Genomes reference panel** (optional, 1.57 GB). Needed to clump on your own computer: your own GWAS files, or OpenGWAS
  downloads if you prefer local clumping. Without it, OpenGWAS downloads are clumped on the OpenGWAS server.

**Offline mode** (also in Settings) uses only data already on your computer. The simulated demo traits are always available, so
you can try every step without a token.

## The Data library

The Data library holds the traits on your computer. Each study picks its exposures and outcomes from it.

### Downloading from OpenGWAS

Choose **Find more on OpenGWAS**, search by trait name or dataset ID, and pick a dataset.

- **Use as exposure** downloads the independent instruments: the variants below the P-value threshold (5e-8 by default), with one
  variant kept per region (clumping). You choose where clumping happens:
  - **OpenGWAS, clumped live (recommended)**: the server clumps the variants for your threshold.
  - **OpenGWAS, stored clumps (fastest)**: the server sends clumps it stored in advance, which can keep a weaker variant as the lead
    of a region.
  - **This computer**: OpenGWAS sends every variant below the threshold and MR Studio clumps them with PLINK and the 1000 Genomes
    panel. You then see every candidate variant and why it was kept or removed.
- **Use as outcome** looks up the instruments of one exposure in the outcome GWAS. You choose that exposure. If an instrument is
  missing from the outcome, a proxy (a variant in strong linkage with it, r² at least 0.8 by default) is used instead.

- **Complete GWAS** downloads every variant of the dataset (its GWAS-VCF file, usually 0.2 to 1.5 GB). Before it starts, MR Studio
  shows the file size, the space it needs (the download and the standard copy made from it), the free space on your disk and the
  memory needed to read it. If the disk does not have enough space, the download cannot start; if little space would be left, you
  are warned. The download continues in the background and is listed under **Complete GWAS downloads** with its progress; **Stop**
  and **Resume** continue where it stopped, also after closing the app. A complete GWAS can be an outcome for any exposure, an
  exposure (clumped on this computer), and the data of the analyses that need full summary statistics. Complete downloads count
  against your OpenGWAS allowance.

An outcome downloaded with **Use as outcome** holds only the variants of the exposure you chose, so it can be paired only with that exposure. To
use the same outcome with another exposure, download it again for that exposure. The study's Check step tells you if a pair
cannot be used for this reason.

### Importing your own file

Choose **Import a GWAS file** and pick a text table (TSV, CSV or space separated, optionally `.gz`), or a GWAS-VCF file (`.vcf` or
`.vcf.gz`, the format of OpenGWAS and the GWAS Catalog). A VCF is read directly: you give only its name, the ancestry and the type of
trait (and the numbers of cases and controls if the file does not list them). MR Studio reads the columns,
shows what it recognised (variant ID, alleles, effect, standard error, P value, allele frequency, sample size), and asks about
anything unclear: for example whether effects are log odds ratios or odds ratios. It then saves the file in a standard layout.

When you import a file, you say whether it is:

- a set of **independent instruments** already (for example lead variants from a published table), or
- **full summary statistics**, which MR Studio clumps on your computer when it is used as an exposure (this needs rsIDs and the
  1000 Genomes panel), or
- used as an **outcome only**.

Files from other sources (for example a consortium website or a paper's supplementary table) often differ in layout. MR Studio
handles the common differences and lists every change on the **Check and save** step, before anything is saved:

- **Excel workbooks** (`.xlsx`, `.xls`): the sheet with the table is found (you can pick another), and title, note and footnote
  rows are left out.
- **Layout:** comment lines before the heading (`##`, or a heading that starts with `#`), a byte-order mark, decimal commas
  (`0,05`), and chromosome names such as `chr1` or `23`.
- **Headings:** case, spaces and punctuation do not matter, and common names are recognised (`MarkerName`, `A1`, `ALT`,
  `effect_allele`, `StdErr`, `P-value`, `Freq1`, `POS_hg38` and many others). Each match shows how sure it is; check those marked
  MEDIUM or LOW. Always check which allele the effect belongs to: the app shows your first row and asks.
- **A minor allele frequency** (`MAF`) is recognised as such and not used as the effect allele frequency.
- **-log10(P)** columns are recognised and converted.
- **A missing standard error** is derived from the effect and its P value (SE = |b| / z; Altman and Bland 2011). The P values must
  be two-sided; rounded P values give approximate standard errors.
- **Z scores only:** the effect and its standard error are derived per standard deviation of the trait from the Z score, the
  sample size and an allele frequency (Zhu et al. 2016). Not for binary traits.
- **No allele frequency:** frequencies can be taken from the 1000 Genomes panel of the file's ancestry (Settings, LD reference
  panel), matched by rsID and alleles.
- **GRCh38 positions** are converted to GRCh37 with the Ensembl chain file included in the app; variants that cannot be converted
  are dropped and counted. A file with positions but no rsIDs needs its genome build.

Your file is never changed. The recipe saved with the data records the file's checksum, the column mapping, your answers and every
change. For a mediator, choose **Any role (exposure, mediator or outcome)**.

## A study in seven steps

Choose **All studies**, then start a new study. The left side shows the seven steps. Your progress is saved in the study folder
after every step, so you can close the app and continue later.

### 1. Question

Give the study a name and choose what you want to find out:

| Question | Analysis |
|---|---|
| Does one thing cause another? | Two-sample MR of one or more exposures on one or more outcomes |
| Which of several factors really matters? | Multivariable MR (MVMR): the direct effects of two to six related exposures on each outcome |
| Does the effect act through something else? | Mediation: the total effect split into a direct part and parts through one or more mediators |
| Which way does the effect go? | Bidirectional MR: every pair in both directions |

The question decides which methods the study offers.

### 2. Data

Pick the exposures and outcomes from the Data library. Under each dataset a short line says what it is ready for in that role (for
example "Instruments, clumped on OpenGWAS", "Clumped here at Check", "Complete GWAS: any exposure; genome-wide methods" or "Looked up
for LDL"); a line in the warning colour means something is missing, such as the 1000 Genomes panel a dataset is clumped with. Each list
also offers "Import a file" and "Find on OpenGWAS": when you come back to the study, what you imported or downloaded is selected and
your other choices are kept. In a two-sample or bidirectional study, every exposure is paired with every
outcome. An OpenGWAS outcome download holds only the variants looked up for the exposure it was made for, so pick it once: before the
Check step MR Studio looks the same outcome up again for each other exposure that is an OpenGWAS instrument download (same proxy
settings), saves it in the Data library as "<outcome> for <exposure>" and reuses it afterwards. When you picked two downloads of the
same outcome, each exposure gets its own and the same pair is kept once. For an exposure clumped on this computer, use the outcome's
complete GWAS or an imported file. If you want different outcomes for different exposures, use one study per exposure.

Multivariable studies, and the multivariable step of a mediation study, look up each instrument's effect on every exposure. With
OpenGWAS downloads MR Studio does this itself on OpenGWAS before the Check step (as TwoSampleMR's `mv_extract_exposures` does); with
imported files, every exposure must be the full summary statistics. The outcome is looked up again at all joint instruments, so the
exposure an outcome download was made for does not matter.

### 3. Check

The app checks every pair before you choose methods, and shows each result with what it means:

- **Instrument strength** (F statistics), and how many instruments each pair has. F = (beta / SE)^2 is given for every exposure. The
  F from the variance explained (R2 from EAF) assumes the exposure is in SD units: for another unit R2 comes from beta, SE and N,
  and for a binary exposure it is not given.
- **What happened to each instrument**: found in the outcome or not, found through a proxy, dropped because the alleles could not be
  matched or because it is palindromic with an ambiguous strand.
- **Ancestry** of each GWAS, as OpenGWAS lists it or as you declared it.
- **Sample overlap**: whether the two GWAS share participants. This cannot be calculated from instruments alone, so you answer it
  (**Answer now**); the answer is reported in the draft text.
- **Units** and other data problems.

Items marked REVIEW do not stop the study. Each is reported as a limitation in the draft text. A pair with a data problem does not run.

### 4. Methods

The Methods screen shows the methods of your study's question. Each card says whether the method is **Available**, needs a
**Check** (it can run, with a note to read), or is **Locked**, with the reason (for example too few instruments). Nothing is
hidden.

**Use the standard publication set** selects the core set of a two-sample MR paper:

- inverse-variance weighted (IVW), with the Steiger direction test and leave-one-out analysis
- MR-Egger
- weighted median
- weighted mode
- simple mode

MR-PRESSO, radial MR and MR-Lasso, which the instrument review on Results uses, are not in the set; tick them when you want that review.

Each rests on a different assumption about pleiotropy. Add others as needed. **Customise** on a card changes its settings; the
settings are written into the script and the report.

A few methods use the full GWAS of both traits instead of instruments (LHC-MR, MRlap, MR-APSS, CAUSE). They run from the Results
page; see [below](#analyses-that-need-full-summary-statistics).

### 5. Run

Run analysis. Each pair is saved as soon as it finishes, so a closed laptop never loses finished work, and a resumed run gives the
same numbers as an uninterrupted one. Each method runs on its own: if one fails, the others still finish, and the activity log says
which failed and why.

The run works in the background, as do the reanalyses on the Results screen (leaving variants out, outliers, ties), mediation,
bidirectional MR-cML, LD score regression and the publication package. While it runs you can open Studies, look at other studies and
open the running one again: it shows its Run screen with the progress so far. The Studies list marks it RUNNING. One study runs at a
time; when a run ends while you are looking at another study, a note says so.

### 6. Results

For each pair:

- the estimate of every method, with its confidence interval, in a table and a forest plot;
- the scatter plot of the variant effects;
- diagnostics: heterogeneity (Cochran's Q), the MR-Egger intercept, the Steiger direction test, leave-one-out;
- **Which instruments drive the result?**: outlier analyses (MR-PRESSO, radial MR) and a rerun without chosen variants;
- the analyses that need full summary statistics, if you run them.
- **More figures** (Show the figures): the forest, funnel and scatter plots of MendelianRandomization, the all-methods plot, the
  mode-estimate density, the MR-Egger and SIMEX plots and Figures 5 to 8 of mendelianrandomization, as in the publication package,
  and an interactive scatter plot: hover a point to see the variant; pinch or drag a box to zoom, scroll with two fingers to move
  around while zoomed in, and double-click to reset. A multivariable pair shows
  the multivariable scatter plot (outcome associations against the fitted values), each variant's contribution to Cochran's Q
  with the Bonferroni cut-off, the radial MVMR plots without and with correction (when the optional RMVMR package is installed) and
  the interactive scatter plot. The plots of single-exposure pairs rest on each variant's ratio estimate, which a variant does not
  have with several exposures, so they are not drawn for a multivariable pair.

**Reading the result.** Check whether the methods agree in direction and roughly in size, whether the MR-Egger intercept is close
to zero, and whether heterogeneity is large. When the methods disagree, the result depends on which pleiotropy assumption holds,
and the report says so.

### 7. Publish

Build your publication package. Choose any of:

| Item | What it is |
|---|---|
| Full report | Summary table, methods, results, figures, limitations and references |
| Figures | Forest and scatter plots for every pair, at your journal's column width |
| Supplement workbook | Estimates, settings and the instruments of every pair |
| Draft Methods, Results and Limitations | Filled with your numbers and settings, ready to edit |
| STROBE-MR checklist | All 20 items, each marked as filled from your study, partly filled, or for you to write |
| Citations | Every method and tool named in the text, checked against PubMed |
| Scripts and data tables | The R scripts, the analysis data and every result table of every pair |
| Reproducibility capsule | Every pair folder with scripts, data, results and the replay check, plus the R version |

The draft text names each method and cites it. Read it, edit it, and add what only you know (for example why you chose the traits).

## Analyses that need full summary statistics

These run from a pair's Results page, in the background at low priority, after a check that enough memory is free:

| Analysis | What it adds | Needs |
|---|---|---|
| LD score regression | Heritability, genetic correlation, and the cross-trait intercept (a sign of shared participants) | Full GWAS of both traits; LD scores of their ancestry (European: 34 MB download; others computed from the panel) |
| LHC-MR | Effects in both directions with a heritable confounder | Full GWAS of both traits, both European; the LHC-MR LD files (524 MB); several hours |
| MRlap | IVW corrected for weak instruments, winner's curse and sample overlap | Full GWAS of both traits; LD scores of their ancestry |
| MR-APSS | Allows for pleiotropy and sample structure, with more instruments | Full GWAS of both traits; LD scores and 1000 Genomes panel of their ancestry |
| CAUSE | Tests a causal effect against correlated pleiotropy (a shared heritable factor) | Full GWAS of both traits; 1000 Genomes panel of their ancestry |

"Full GWAS" means genome-wide summary statistics with rsIDs and sample size in the Data library: downloaded with **Complete GWAS**
from OpenGWAS, or imported (a table or a GWAS-VCF file). A pair made from OpenGWAS instrument downloads uses the complete GWAS of the
same OpenGWAS datasets when you have downloaded them: the Results page offers the downloads it needs.

**Ancestry and LD reference.** These analyses need LD from the same population as the GWAS. Both GWAS must be of one ancestry, and
its 1000 Genomes panel must be installed (Settings, LD reference panel). For European GWAS the LDSC authors' LD scores are used; for
another ancestry, the LD score regression section computes LD scores from that panel (about 30 minutes, once for this computer).
When a GWAS is labelled Mixed by OpenGWAS, or has no ancestry reported, no panel matches it: the section shows a warning, you choose
one of the installed panels (European is proposed), and the analysis runs only after you tick **I accept this and take responsibility
for it**. Your choice and the date are saved beside the result and stated in the Limitations of the draft text. LHC-MR uses its
authors' European LD files, which cannot be replaced, so it needs two European GWAS.

## Checking the results in R

Each pair has its own folder in the study, with the analysis data (`01_data`), the R scripts (`02_scripts`), the console output
and the result tables. To rerun a pair in your own R:

1. Open the pair folder and open `VERIFY_IN_R.R` in RStudio (or `MR_Studio_pair.Rproj`).
2. Click **Source**.

The script copies the pair's data and scripts into a new folder (`09_verification/run_<date and time>`), runs them in a fresh R
session and compares every number with the results MR Studio saved. Your original results are never changed. If a package is
missing, it shows the lines to install it. `README_RERUN.md` in the same folder explains the package versions. After each run MR Studio does the same check itself (the
replay check), and the report states its result.

## Where things are kept

- **Studies and the Data library**: by default the `MR Studio` folder in your home folder (you can choose another in Settings). Each study is an ordinary folder that you can back
  up, share or archive. Do not keep it in a folder that a cloud service syncs while a run is going.
- **Inside a study**: `01_study_plan` (the question and settings), `02_datasets` (data prepared for the study, such as LD score inputs), `03_pairs` (one folder per
  analysis: `01_data`, `02_scripts`, `03_console`, `04_tables`, `07_replay_check`, and `08_outlier_analysis` once an outlier analysis is run), `04_across_pairs`
  (analyses that combine pairs, such as mediation, MR-cML, LD score regression, MRlap and CAUSE), `05_publication` (each package you build) and `archive` (earlier
  runs of a pair you ran again). Only `01_study_plan` and `03_pairs` are made with the study; the others appear when the analysis that uses them runs.
- **Reference data** (1000 Genomes panel, LD scores, the OpenGWAS catalogue): the app's data folder, shared by all studies.
- **Your OpenGWAS token**: in memory only.

## Troubleshooting

| What you see | What to do |
|---|---|
| "No token yet" or "token expired" | Paste a new token in Settings > Data and access (api.opengwas.io; tokens expire after a while). |
| An outcome "cannot be paired" with an exposure | It was downloaded for another exposure. Download it again for this exposure, or use one study per exposure. |
| A method is locked | The card gives the reason, for example too few instruments for that method. |
| A pair finished with some methods failed | The activity log names each failed method and its error; the others are saved. |
| Clumping on this computer is not offered | Install the 1000 Genomes panel in Settings > Data and access. |
| "Not enough free space" before a complete GWAS download | Free space on the disk shown, or choose another folder for studies and the Data library in Settings. |
| A complete GWAS download stopped | Choose **Resume** under Complete GWAS downloads; it continues where it stopped and asks OpenGWAS for a new link. |
| Slow or failing OpenGWAS requests | OpenGWAS may be busy. Try again later; finished downloads are kept. |

Each method's reference, and how MR Studio was checked against it, is listed in `validation/README.md`.
