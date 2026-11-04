# Builds the small synthetic LD reference panel used by the clumping tests (test_fixtures/ld_panel/EUR.bed/.bim/.fam).
# Synthetic on purpose: the tests must know exactly which variants are in LD. Not real people, no licence restrictions.
# Run from the repository root with PLINK 1.9 available:  Rscript prototype/engine/test_fixtures/make_ld_fixture.R <path to plink>
args <- commandArgs(TRUE); plink <- if (length(args)) args[1] else "plink"
set.seed(20261105); n <- 300
base_a <- rbinom(n, 2, 0.30); base_b <- rbinom(n, 2, 0.40)
noisy <- function(g, flip = 0.03) { k <- runif(length(g)) < flip; g[k] <- rbinom(sum(k), 2, 0.35); g }
snps <- rbind(
  data.frame(id = sprintf("rs91%05d", 1:4), chr = 1, pos = 1000000 + (0:3) * 100, block = "A"),   # LD block A
  data.frame(id = sprintf("rs92%05d", 1:3), chr = 1, pos = 15000000 + (0:2) * 100, block = "B"),   # LD block B, 14 Mb away (outside a 10,000 kb window)
  data.frame(id = sprintf("rs93%05d", 1:5), chr = 1, pos = c(40, 60, 80, 100, 120) * 1e6, block = "indep"),   # 20 Mb apart, beyond any 10,000 kb window
  data.frame(id = "rs9400001", chr = 2, pos = 3e6, block = "indep"))
geno <- sapply(seq_len(nrow(snps)), function(i) switch(snps$block[i], A = noisy(base_a), B = noisy(base_b), rbinom(n, 2, 0.35)))
al <- function(g) ifelse(g == 0, "C C", ifelse(g == 1, "A C", "A A"))
ped <- cbind(sprintf("F%03d", 1:n), sprintf("I%03d", 1:n), 0, 0, 1, -9, apply(geno, 2, al))
tmp <- tempfile(); utils::write.table(ped, paste0(tmp, ".ped"), quote = FALSE, row.names = FALSE, col.names = FALSE)
utils::write.table(data.frame(snps$chr, snps$id, 0, format(snps$pos, scientific = FALSE, trim = TRUE)), paste0(tmp, ".map"), quote = FALSE, row.names = FALSE, col.names = FALSE)
out <- file.path("prototype", "engine", "test_fixtures", "ld_panel", "EUR")
system2(plink, c("--file", tmp, "--make-bed", "--out", out, "--silent"))
unlink(paste0(out, c(".log", ".nosex")))
cat("Wrote", out, ".bed/.bim/.fam\n")
