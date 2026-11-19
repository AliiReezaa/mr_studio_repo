# Validation study 13: bidirectional MR with constrained maximum likelihood (Round 63)

Reference: Xue H, Pan W. Robust inference of bi-directional causal relationships in presence of correlated pleiotropy with GWAS summary
data. PLoS Genet 2022;18:e1010205. doi:10.1371/journal.pgen.1010205, PMID 35576237. Packages from GitHub xue-hr, pinned: BiDirectCausal
(BiDirMRcML, BiDirCDcML; MIT), MRcML and MRCD. BiDirectCausal calls MRcML and MRCD without importing them, so the script attaches both.

In a bidirectional study, MR Studio writes 04_across_pairs/bidirectional_cml/<forward pair>.R, which reads the forward pair (instruments of
X with their effects on Y) and the reverse pair (instruments of Y with their effects on X), takes their union and runs both methods in
all four versions (with or without instrument screening, with or without data perturbation; screening with perturbation is the main
result). BiDirCDcML works on standardised effects (per SD), BiDirMRcML on the per-unit scale of the other methods.

    Rscript --vanilla validation/xue2022_bidir_cml/run_validation.R

Result: PASSED. (1) The script equals the package called directly on the same union, for all 16 numbers (largest difference 1e-16).
(2) On the simulated traits B (0.25 forward, 0.10 back; 80 variants) BiDirMRcML gives 0.274 and 0.118, whose intervals cover both.
The paper's applications use genome-wide data for 16 traits that are not bundled, so they are not reproduced here.
