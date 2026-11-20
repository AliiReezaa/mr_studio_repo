# Tutorial: LDL cholesterol and coronary heart disease

This tutorial runs a complete two-sample MR study with public data from OpenGWAS: the effect of LDL cholesterol on coronary heart
disease (CHD). It is the applied example of the MR-Base paper (Hemani et al. 2018, eLife 7:e34408), so you can compare your result
with a published one. It takes about 15 minutes, most of it downloads.

You need MR Studio installed and a free OpenGWAS token (see the [user guide](USER_GUIDE.md#before-your-first-study)).

| | Dataset | OpenGWAS ID |
|---|---|---|
| Exposure | LDL cholesterol, Global Lipids Genetics Consortium | `ieu-a-300` |
| Outcome | Coronary heart disease, CARDIoGRAMplusC4D | `ieu-a-7` |

## 1. Paste your token

Open **Settings > Data and access**, paste your OpenGWAS token and save. If the OpenGWAS catalogue is not on your computer yet,
download it there too (once, about a minute).

## 2. Download the exposure

1. Open the **Data library** and choose **Find more on OpenGWAS**.
2. Search for `ieu-a-300` and select the LDL cholesterol dataset.
3. Choose **Use as exposure**.
4. Keep the P-value threshold at 5e-8 and choose **OpenGWAS, clumped live (recommended)**.
5. Name it, for example `LDL cholesterol`, and download.

MR Studio downloads the genome-wide significant variants and keeps one per region (r² below 0.001 within 10,000 kb). The Data library
now lists the exposure with its number of instruments.

## 3. Download the outcome for this exposure

1. Choose **Find more on OpenGWAS** again and search for `ieu-a-7`.
2. Choose **Use as outcome**, and pick `LDL cholesterol` as the exposure whose instruments are looked up.
3. Keep **If an instrument is missing, use a proxy** on (minimum r² 0.8), name it, for example `CHD for LDL`, and download.

Only the LDL instruments are looked up in the CHD GWAS, so this outcome can be paired only with `LDL cholesterol`.

## 4. Create the study

1. Go to **All studies** and start a new study.
2. **Question**: choose **Does one thing cause another?** and give the study a name.
3. **Data**: tick `LDL cholesterol` as the exposure and `CHD for LDL` as the outcome. The study has one pair.

## 5. Check the data

The Check screen reports, for the pair:

- **Instrument strength**: the F statistic of every instrument. LDL instruments are strong, so weak-instrument bias is unlikely.
- **What happened to each instrument**: how many were found in the CHD GWAS, how many through a proxy, how many were dropped as
  palindromic or because the alleles could not be matched.
- **Ancestry**: the ancestry OpenGWAS lists for each GWAS (MR Studio shows it but cannot verify it).
- **Sample overlap**: choose **Answer now**. Compare the cohorts listed in the two consortia's papers and answer as well as you can;
  if some cohorts are in both, say so. Your answer goes into the draft text and its limitations.

Then continue to the methods.

## 6. Choose the methods

Choose **Use the standard publication set**. It selects IVW (with the Steiger direction test and leave-one-out), MR-Egger, weighted
median, weighted mode and MR-PRESSO. Leave the settings at their defaults. CHD is a yes/no trait, so the Steiger direction test asks
for its population prevalence (it puts the variance explained on the liability scale); enter a value from the literature and note
it in your methods. Run the analysis.

## 7. Read the results

CHD is a binary outcome, so the estimates are log odds ratios per standard deviation of LDL cholesterol; the report also gives
odds ratios. When this example was run through MR Studio (November 2026; 77 instruments with live clumping):

| Method | MR Studio, odds ratio of CHD per SD of LDL cholesterol | Hemani et al. 2018 |
|---|---|---|
| IVW (random effects) | 1.51 | 1.45 (1.30 to 1.62) |
| MR-Egger | 1.65 | 1.85 (1.48 to 2.32) |
| Weighted median | 1.49 | 1.56 (1.43 to 1.70) |
| Weighted mode | 1.72 | 1.68 (1.56 to 1.80) |

Your numbers can differ slightly: OpenGWAS updates its data and live clumping can choose different lead variants, and the paper
used about 56 variants (inferred from its Cochran's Q). The direction and size should agree: higher LDL cholesterol raises the risk
of CHD. The comparison in detail is in `validation/hemani2018`.

What to look at on the Results screen:

1. **Do the methods agree?** All estimates are above 1 and of similar size. Methods that rest on different pleiotropy assumptions
   giving the same answer makes the result more credible.
2. **The MR-Egger intercept.** Close to zero (here about -0.007, SE 0.005) means no clear sign of directional pleiotropy.
3. **Heterogeneity.** Cochran's Q is large (about 364 on 76 degrees of freedom): the variants' estimates vary more than chance
   allows. This is why IVW uses random effects, and why the outlier analyses matter.
4. **Which instruments drive the result?** MR-PRESSO and radial MR flag outlying variants; leave-one-out shows that no single
   variant decides the result. You can rerun without chosen variants as a sensitivity analysis.
5. **Steiger direction test.** The instruments explain more variance in LDL cholesterol than in CHD, as expected if the effect runs
   from LDL to CHD.

## 8. Build the publication package

Open **Publish**, choose the items you need (for example the full report, the figures, the draft text and the STROBE-MR checklist)
and build the package. Read the draft Methods, Results and Limitations: they contain your numbers, settings and citations, and the
limitations include every item the Check screen asked you to review, such as sample overlap.

## 9. Check it in R (optional)

Open the study folder, then `03_pairs`, then the pair folder. Open `VERIFY_IN_R.R` in RStudio and click **Source**. It reruns every
method in a fresh R session and confirms that the numbers are the same.

## Next steps

- Add more outcomes or exposures (download each outcome for its exposure, as in step 3).
- For multivariable MR, for example LDL cholesterol and triglycerides together, import the full summary statistics of both
  exposures and of the outcome, and choose **Which of several factors really matters?**
