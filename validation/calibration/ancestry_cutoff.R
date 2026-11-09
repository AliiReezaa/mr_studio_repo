# Calibration of the 90% single-ancestry cut-off (data_limits_solutions.md item 3, contracts C3.2), Round 32.
# Question: if an outcome GWAS is called "single-ancestry" when at least 90% of its sample is one group, how large can the bias be
# when the remaining share comes from another ancestry? The cut-off is a labelling rule; the classification itself is NOT BUILT.
#
# Mechanism studied: the instruments are chosen in a European exposure GWAS and tag a causal variant with correlation r_A in
# Europeans. In the other ancestry they tag it with r_B = rho * r_A (rho = 1: same tagging; rho = 0: no tagging, the worst case).
# A GWAS adjusted for ancestry estimates the variant's association within each group and pools them with weights proportional to
# n_k * 2 p_k (1 - p_k), the group's sample size times the variant's genotype variance. With a share pi from ancestry B, the outcome
# association of variant j is scaled by 1 - w_Bj (1 - rho), where w_Bj = pi v_Bj / ((1 - pi) v_Aj + pi v_Bj), while the exposure
# association (European) is not. The ratio estimate of variant j is therefore biased by -w_Bj (1 - rho), and IVW by the IVW-weighted
# mean of these. check_formula() confirms the pooling formula on simulated genotypes.
# Real data: v_A and v_B from the 1000 Genomes frequencies of the BMI instruments (EUR and each other panel ancestry), IVW weights
# from the BMI (ieu-b-40) on coronary heart disease (ieu-a-8) pair. Not covered: confounding by population structure, and allele
# frequency or LD differences in the exposure GWAS.
#
# Usage: Rscript ancestry_cutoff.R <exposure dataset folder> <outcome dataset folder> <panel folder> <plink program> <output folder>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
ex <- utils::read.csv(file.path(args[1], "standardised.csv"), stringsAsFactors = FALSE, colClasses = c(snp = "character", ea = "character", oa = "character"))
ou <- utils::read.csv(file.path(args[2], "standardised.csv"), stringsAsFactors = FALSE, colClasses = c(snp = "character", ea = "character", oa = "character"))
panel <- args[3]; plink <- args[4]; out <- args[5]; work <- tempfile("anc"); dir.create(work)
pops <- c("EUR", "AFR", "AMR", "EAS", "SAS")

# ---- 1. check of the pooling formula on simulated genotypes (two groups, one variant, ancestry-adjusted regression) ----
check_formula <- function(n = 200000, pi = 0.1, pA = 0.3, pB = 0.6, rho = 0.4, seed = 32) {
  set.seed(seed); grp <- runif(n) < pi; p <- ifelse(grp, pB, pA)
  g <- rbinom(n, 2, p); gs <- (g - 2 * p) / sqrt(2 * p * (1 - p))                 # standardised tag genotype within its group
  rr <- ifelse(grp, rho * 0.9, 0.9)                                                 # tag-causal correlation: 0.9 in A, rho * 0.9 in B
  cz <- rr * gs + sqrt(1 - rr^2) * rnorm(n)                                         # causal variant (standardised), correlated with the tag
  c_raw <- cz * sqrt(2 * p * (1 - p))                                               # back on the allele-count scale within group
  y <- 1 * c_raw + rnorm(n)                                                         # large effect, so sampling error is small
  fitted <- coef(lm(y ~ g + grp))[["g"]]
  vA <- 2 * pA * (1 - pA); vB <- 2 * pB * (1 - pB); wB <- pi * vB / ((1 - pi) * vA + pi * vB)
  c(simulated_ratio = fitted / 0.9, formula_ratio = 1 - wB * (1 - rho))
}
chk <- check_formula()
cat(sprintf("pooling check: simulated %.4f, formula %.4f\n", chk[1], chk[2]))

# ---- 2. the analysed pair: harmonised variants with their IVW weights ----
m <- merge(ex[, c("snp", "ea", "oa", "beta", "se")], ou[, c("snp", "ea", "oa", "se")], by = "snp", suffixes = c("_x", "_y"))
m <- m[toupper(m$ea_x) == toupper(m$ea_y) & toupper(m$oa_x) == toupper(m$oa_y) | toupper(m$ea_x) == toupper(m$oa_y) & toupper(m$oa_x) == toupper(m$ea_y), ]
m$w <- m$beta^2 / m$se_y^2                                                          # first-order IVW weight of the ratio estimate
writeLines(m$snp, file.path(work, "snps.txt"))

# ---- 3. genotype variance of each instrument in each ancestry (1000 Genomes, PLINK --freq), aligned to the EUR allele ----
fr <- list()
for (pp in pops) {
  st <- system2(plink, c("--bfile", shQuote(file.path(panel, pp)), "--extract", shQuote(file.path(work, "snps.txt")), "--freq", "--out", shQuote(file.path(work, pp))), stdout = FALSE, stderr = FALSE)
  stopifnot(st == 0L)
  f <- utils::read.table(file.path(work, paste0(pp, ".frq")), header = TRUE, stringsAsFactors = FALSE, colClasses = c(SNP = "character", A1 = "character", A2 = "character"))
  fr[[pp]] <- f[is.finite(f$MAF), c("SNP", "A1", "MAF")]
}
v <- Reduce(function(a, b) merge(a, b, by = "SNP"), lapply(pops, function(pp) { f <- fr[[pp]]; d <- data.frame(SNP = f$SNP, x = 2 * f$MAF * (1 - f$MAF)); names(d)[2] <- pp; d }))
m <- merge(m, v, by.x = "snp", by.y = "SNP")
cat(sprintf("%d instruments with frequencies in all five ancestries\n", nrow(m)))

# ---- 4. IVW relative bias for each other ancestry, share pi and tagging ratio rho ----
grid <- expand.grid(other = setdiff(pops, "EUR"), pi = c(0.02, 0.05, 0.1, 0.15, 0.2, 0.3), rho = c(0, 0.5, 0.8), stringsAsFactors = FALSE)
grid$ivw_bias <- mapply(function(o, p, r) { wB <- p * m[[o]] / ((1 - p) * m$EUR + p * m[[o]]); -sum(m$w * wB * (1 - r)) / sum(m$w) }, grid$other, grid$pi, grid$rho)
grid$worst_variant_bias <- mapply(function(o, p, r) { wB <- p * m[[o]] / ((1 - p) * m$EUR + p * m[[o]]); -max(wB * (1 - r)) }, grid$other, grid$pi, grid$rho)
grid$pooling_weight_B <- mapply(function(o, p) { wB <- p * m[[o]] / ((1 - p) * m$EUR + p * m[[o]]); sum(m$w * wB) / sum(m$w) }, grid$other, grid$pi)
grid$instruments <- nrow(m)
utils::write.csv(grid, file.path(out, "results_ancestry_cutoff.csv"), row.names = FALSE)
utils::write.csv(data.frame(check = "pooling formula on 200,000 simulated people (pi 0.1, pA 0.3, pB 0.6, rho 0.4)", simulated_ratio = chk[1], formula_ratio = chk[2]),
                 file.path(out, "results_ancestry_check.csv"), row.names = FALSE)
print(grid[grid$rho == 0, ])
unlink(work, recursive = TRUE)
