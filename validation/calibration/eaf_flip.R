# Calibration of the allele-frequency flip threshold (|delta EAF| > 0.2, data_limits_solutions.md item 5), Round 32.
# Real data: the effect-allele frequencies (EAF) reported in OpenGWAS downloads are compared with the frequency of the
# same allele in the 1000 Genomes European panel that MR Studio uses for clumping. Correctly aligned variants show how
# far the two differ by chance and by population difference (false flags); the same variants with the EAF reported for
# the other allele (1 - EAF) show which flips a threshold can find (detection).
#
# Usage: Rscript eaf_flip.R <data library folder> <panel folder with EUR.bed/.bim/.fam> <plink program> <output folder>
# The datasets are only read. One dataset per OpenGWAS ID is used (the first one found), only EUR datasets with an EAF column.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
lib <- args[1]; panel <- file.path(args[2], "EUR"); plink <- args[3]; out <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- tempfile("eafflip"); dir.create(work)

# ---- 1. one dataset per OpenGWAS ID, EUR, with an EAF column ----
sets <- list()
for (d in sort(list.dirs(lib, recursive = FALSE))) {
  rj <- file.path(d, "recipe.json"); sf <- file.path(d, "standardised.csv")
  if (!file.exists(rj) || !file.exists(sf)) next
  r <- jsonlite::fromJSON(rj, simplifyVector = FALSE)
  id <- r$source$params$id %||% r$source$id %||% NULL
  if (!identical(r$source$type, "opengwas") || is.null(id) || !identical(r$facts$ancestry, "EUR") || id %in% names(sets)) next
  x <- utils::read.csv(sf, stringsAsFactors = FALSE, colClasses = c(snp = "character", ea = "character", oa = "character"))
  if (!"eaf" %in% names(x)) next
  x <- x[is.finite(x$eaf) & x$eaf > 0 & x$eaf < 1 & grepl("^rs[0-9]+$", x$snp), c("snp", "ea", "oa", "eaf")]
  x$dataset <- id; sets[[id]] <- x
}
d <- do.call(rbind, sets); rownames(d) <- NULL
cat(sprintf("%d datasets (%s), %d variants\n", length(sets), paste(names(sets), collapse = ", "), nrow(d)))

# ---- 2. panel frequencies of these variants (PLINK 1.9 --freq) ----
writeLines(unique(d$snp), file.path(work, "snps.txt"))
st <- system2(plink, c("--bfile", shQuote(panel), "--extract", shQuote(file.path(work, "snps.txt")), "--freq", "--out", shQuote(file.path(work, "frq"))), stdout = FALSE, stderr = FALSE)
stopifnot(st == 0L)
fq <- utils::read.table(file.path(work, "frq.frq"), header = TRUE, stringsAsFactors = FALSE, colClasses = c(SNP = "character", A1 = "character", A2 = "character"))
fq <- fq[is.finite(fq$MAF), ]
m <- merge(d, fq[, c("SNP", "A1", "A2", "MAF", "NCHROBS")], by.x = "snp", by.y = "SNP")

# ---- 3. align: panel frequency of the dataset's effect allele; drop allele mismatches and palindromic variants ----
up <- function(a) toupper(a)
same <- up(m$ea) == m$A1 & up(m$oa) == m$A2; swap <- up(m$ea) == m$A2 & up(m$oa) == m$A1
pal <- (up(m$ea) %in% c("A", "T") & up(m$oa) %in% c("A", "T")) | (up(m$ea) %in% c("C", "G") & up(m$oa) %in% c("C", "G"))
m$p_panel <- ifelse(same, m$MAF, ifelse(swap, 1 - m$MAF, NA_real_))
n_all <- nrow(m); n_mismatch <- sum(!(same | swap)); n_pal <- sum(pal & (same | swap))
a <- m[(same | swap) & !pal, ]
a$maf <- pmin(a$p_panel, 1 - a$p_panel)
a$d_ok <- abs(a$eaf - a$p_panel)            # correctly aligned: chance and population difference
a$d_flip <- abs((1 - a$eaf) - a$p_panel)    # same variant with the EAF of the other allele (a flip)
cat(sprintf("matched %d, allele mismatch %d, palindromic %d, analysed %d\n", n_all, n_mismatch, n_pal, nrow(a)))

# ---- 4. summaries ----
q <- stats::quantile(a$d_ok, c(0.5, 0.9, 0.95, 0.99, 0.999, 1))
summ <- data.frame(datasets = paste(names(sets), collapse = "; "), variants_matched = n_all, allele_mismatch = n_mismatch, palindromic = n_pal, analysed = nrow(a),
                   median = q[1], p90 = q[2], p95 = q[3], p99 = q[4], p999 = q[5], max = q[6], panel_chromosomes = stats::median(a$NCHROBS))
by_set <- do.call(rbind, lapply(split(a, a$dataset), function(s) data.frame(dataset = s$dataset[1], analysed = nrow(s), median = stats::median(s$d_ok),
  p99 = stats::quantile(s$d_ok, 0.99), max = max(s$d_ok), above_0.1 = mean(s$d_ok > 0.1), above_0.2 = mean(s$d_ok > 0.2))))
bins <- cut(a$maf, c(0, 0.1, 0.2, 0.3, 0.4, 0.5), include.lowest = TRUE)
thr <- do.call(rbind, lapply(c(0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3), function(t) {
  det <- tapply(a$d_flip > t, bins, mean); nb <- tapply(a$d_flip, bins, length)
  data.frame(threshold = t, false_flag_rate = mean(a$d_ok > t), false_flags = sum(a$d_ok > t), detection_all = mean(a$d_flip > t),
             detection_maf_0_0.1 = det[1], detection_maf_0.1_0.2 = det[2], detection_maf_0.2_0.3 = det[3], detection_maf_0.3_0.4 = det[4], detection_maf_0.4_0.5 = det[5],
             n_maf_bins = paste(nb, collapse = "/"))
}))
rownames(summ) <- rownames(by_set) <- rownames(thr) <- NULL
utils::write.csv(summ, file.path(out, "results_eaf_summary.csv"), row.names = FALSE)
utils::write.csv(by_set, file.path(out, "results_eaf_by_dataset.csv"), row.names = FALSE)
utils::write.csv(thr, file.path(out, "results_eaf_thresholds.csv"), row.names = FALSE)
print(summ); print(by_set); print(thr)
unlink(work, recursive = TRUE)
