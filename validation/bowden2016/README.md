# Validation study 1: Bowden et al. 2016, Table 4

Replication of the applied example in Bowden J, Davey Smith G, Haycock PC, Burgess S. *Consistent estimation in Mendelian
randomization with some invalid instruments using a weighted median estimator.* Genet Epidemiol 2016;40:304-14,
doi:10.1002/gepi.21965 (PMID 27061298). The example is the effect of LDL cholesterol, HDL cholesterol and triglycerides on
coronary artery disease. It has two variant sets per lipid: all variants associated with the lipid, and only the variants
more strongly associated with that lipid than with the other two ("primary").

## Data

`dodata.csv`, deposited by Stephen Burgess on figshare (doi:10.6084/m9.figshare.1116328, CC BY 4.0, MD5
`49485561fb97887a16805d1a1652b21c`). It holds 185 variants from Supplementary Table 1 of Do et al. 2013 (Nat Genet,
doi:10.1038/ng.2795): betas and P values for LDL-c, HDL-c and triglycerides (GLGC, SD units) and for CAD (CARDIoGRAM,
log odds). The file has no standard errors, so they are derived as |beta| / z(P/2), on the log scale.

Web Appendix 4 of the paper gives the number of variants per analysis (73, 85, 47; primary 61, 69, 31). These counts are
reproduced exactly with P < 1e-8. The threshold itself was lost from the extracted text of the paper; 1e-8 is inferred from
the counts.

## How to run

From the repository root:

    Rscript --vanilla validation/bowden2016/run_validation.R

The script works in a throw-away studies folder. It imports the four traits through MR Studio's importer, pairs and
harmonises them with `pair_input`, and runs IVW (multiplicative random effects; TwoSampleMR mr_ivw since Round 30, which on these data equals mr_ivw_mre because the variants are over-dispersed, so no number changed), MR-Egger, weighted median and simple
median through the generated scripts with the replay check, twice: with the strand left as "Not sure" (the default) and with
both files declared to be on the forward strand (Round 26). The second is true for this data, because all four traits come
from one table with the same effect allele. For comparison it also runs TwoSampleMR directly with every variant kept.
`results.csv` holds all four columns (paper, app, app forward strand, all kept).

## Result (2026-11-05)

- **Every variant kept:** the deterministic methods (IVW, MR-Egger slope and intercept) match Table 4 to within 0.0022 in
  the estimate and 0.0014 in the SE. HDL and triglycerides match to the third decimal. The LDL differences come from one
  variant, rs4942486: its CAD beta is 0 with P = 1, so no SE can be derived and it is left out (72 instead of 73).
  The median methods differ by up to 0.007, as expected from bootstrap SEs and from the one LDL variant.
- **Through MR Studio:** all 30 replay checks are identical. Triglycerides (primary: no palindromic variant) match Table 4
  exactly. MR Studio drops the palindromic variants (A/T or C/G) because the file has no allele frequencies:
  1 for LDL and triglycerides (all), 5 for HDL. For HDL this moves IVW from -0.254 to -0.312 and MR-Egger from -0.013 to
  -0.060. The conclusion of the paper still holds in MR Studio: the IVW effect of HDL is not supported by the weighted
  median (-0.081, SE 0.083) or by MR-Egger (-0.060, SE 0.146).
- **Through MR Studio, both files declared forward strand (Round 26):** the palindromic variants are kept, and the result
  equals TwoSampleMR with every variant kept exactly (difference 0 for IVW and MR-Egger). HDL IVW is -0.254 (0.070), as in
  Table 4. All 30 replay checks are identical.
- **Not compared:** the penalized weighted median (NOT BUILT in MR Studio).
- **Found and fixed on the way:** pairing an imported file without a sample-size column failed inside
  `TwoSampleMR::format_data` (Round 26, test added).

The two datasets share participants (ADVANCE, deCODE, LURIC: 9.3% of GLGC, 41.9% of CARDIoGRAM, Web Appendix 4).
