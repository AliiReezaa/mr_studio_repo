# Import formats: writes one exposure and one outcome dataset in the layouts that GWAS summary statistics come in from real sources,
# all made from the same two datasets of the Data library (an exposure's instruments and the outcome looked up at them), so every
# import must reproduce the same pair. Run from the repository root:
#   Rscript validation/import_formats/make_corpus.R <exposure dataset folder> <outcome dataset folder> <output folder>
# Each layout follows the documentation of its source: column names, allele conventions, genome build, number formats, quirks.
args <- commandArgs(TRUE)
ex_dir <- if (length(args) >= 1) args[1] else file.path(path.expand("~"), "MR Studio", "Data library", "LDL")
out_dir <- if (length(args) >= 2) args[2] else file.path(path.expand("~"), "MR Studio", "Data library", "CHD for LDL v2")
dest <- if (length(args) >= 3) args[3] else file.path(tempdir(), "import_corpus")
source("prototype/engine/liftover_core.R"); options(mrgui.engine_dir = "prototype/engine")
dir.create(dest, recursive = TRUE, showWarnings = FALSE); dest <- normalizePath(dest)
rd <- function(d) { x <- as.data.frame(nanoparquet::read_parquet(file.path(d, "standardised.parquet"))); x[order(as.numeric(x$chr), x$pos), c("snp", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p", "n")] }
X <- rd(ex_dir); Y <- rd(out_dir)

# GRCh37 -> GRCh38 by inverting the bundled GRCh38 -> GRCh37 chain blocks (only to write the GRCh38 files; the app converts them back)
to38 <- function(chr, pos) {
  b <- liftover_blocks(); out <- data.frame(pos = rep(NA_real_, length(pos)), rev = FALSE)
  for (i in seq_along(pos)) {
    bb <- b[b$qchr == chr[i], ]; q <- pos[i] - 1
    qq <- ifelse(bb$qstrand == "+", q, bb$qsize - 1 - q)   # the position on the block's strand
    k <- which(qq >= bb$qstart & qq < bb$qstart + bb$size); if (length(k) != 1L) next
    out$pos[i] <- bb$tstart[k] + (qq[k] - bb$qstart[k]) + 1; out$rev[i] <- bb$qstrand[k] == "-"
  }
  out
}
comp <- function(a) chartr("ACGT", "TGCA", a)
# A minimal .xlsx (every cell as text), written with the system zip: no Excel package is needed to make the test file
write_xlsx <- function(m, path) {
  colref <- function(j) { s <- ""; while (j > 0) { r <- (j - 1) %% 26; s <- paste0(LETTERS[r + 1], s); j <- (j - 1) %/% 26 }; s }
  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE)
  rows <- vapply(seq_len(nrow(m)), function(i) { v <- m[i, ]; k <- which(nzchar(v))
    paste0('<row r="', i, '">', paste0('<c r="', vapply(k, colref, ""), i, '" t="inlineStr"><is><t>', esc(v[k]), '</t></is></c>', collapse = ""), '</row>') }, "")
  td <- tempfile("xlsx"); dir.create(file.path(td, "_rels"), recursive = TRUE); dir.create(file.path(td, "xl", "_rels"), recursive = TRUE); dir.create(file.path(td, "xl", "worksheets"))
  ns <- "http://schemas.openxmlformats.org/"
  writeLines(paste0('<?xml version="1.0" encoding="UTF-8"?><Types xmlns="', ns, 'package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>'), file.path(td, "[Content_Types].xml"), useBytes = TRUE)
  writeLines(paste0('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="', ns, 'package/2006/relationships"><Relationship Id="rId1" Type="', ns, 'officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'), file.path(td, "_rels", ".rels"))
  writeLines(paste0('<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="', ns, 'spreadsheetml/2006/main" xmlns:r="', ns, 'officeDocument/2006/relationships"><sheets><sheet name="Table S2" sheetId="1" r:id="rId1"/></sheets></workbook>'), file.path(td, "xl", "workbook.xml"))
  writeLines(paste0('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="', ns, 'package/2006/relationships"><Relationship Id="rId1" Type="', ns, 'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>'), file.path(td, "xl", "_rels", "workbook.xml.rels"))
  con <- file(file.path(td, "xl", "worksheets", "sheet1.xml"), "wb"); writeLines(enc2utf8(paste0('<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="', ns, 'spreadsheetml/2006/main"><sheetData>', paste(rows, collapse = ""), '</sheetData></worksheet>')), con, useBytes = TRUE); close(con)
  old <- setwd(td); on.exit(setwd(old)); unlink(path); utils::zip(path, c("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/_rels/workbook.xml.rels", "xl/worksheets/sheet1.xml"), flags = "-X -q")
}
p_txt <- function(p) formatC(p, format = "g", digits = 6)
wtab <- function(d, f, sep = "\t", gz = FALSE, eol = "\n", quote = FALSE, na = "NA", dec = ".") {
  con <- if (gz) gzfile(file.path(dest, f), "w") else file(file.path(dest, f), "wb"); on.exit(close(con))
  utils::write.table(d, con, sep = sep, quote = quote, row.names = FALSE, na = na, eol = eol, dec = dec)
}
made <- list(); reg <- function(file, role, layout, note) made[[length(made) + 1L]] <<- data.frame(file = file, role = role, layout = layout, note = note)
both <- function(fun) { fun(X, "exposure"); fun(Y, "outcome") }
tag <- function(role) if (role == "exposure") "x" else "y"

# 1 GWAS Catalog summary statistics format (GWAS-SSF, GRCh38): effect and other allele, rsid column, variant_id chr_pos_other_effect
both(function(d, role) { l <- to38(d$chr, d$pos); ea <- ifelse(l$rev, comp(d$ea), d$ea); oa <- ifelse(l$rev, comp(d$oa), d$oa)
  o <- data.frame(chromosome = d$chr, base_pair_location = l$pos, effect_allele = ea, other_allele = oa, beta = d$beta, standard_error = d$se,
                  effect_allele_frequency = d$eaf, p_value = p_txt(d$p), variant_id = paste(d$chr, l$pos, oa, ea, sep = "_"), rsid = d$snp, n = d$n)[!is.na(l$pos), ]
  f <- sprintf("01_gwas_catalog_ssf_%s.tsv.gz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "GWAS Catalog GWAS-SSF (GRCh38)", "GRCh38 positions, rsid column") })
# 2 GWAS Catalog harmonised (hm_ columns, GRCh38), odds ratio and its interval for the outcome
both(function(d, role) { l <- to38(d$chr, d$pos); ea <- ifelse(l$rev, comp(d$ea), d$ea); oa <- ifelse(l$rev, comp(d$oa), d$oa)
  o <- data.frame(hm_variant_id = paste(d$chr, l$pos, oa, ea, sep = "_"), hm_rsid = d$snp, hm_chrom = d$chr, hm_pos = l$pos, hm_other_allele = oa, hm_effect_allele = ea,
                  hm_beta = d$beta, hm_odds_ratio = if (role == "outcome") exp(d$beta) else NA, hm_ci_lower = NA, hm_ci_upper = NA, hm_effect_allele_frequency = d$eaf, hm_code = 10,
                  variant_id = d$snp, p_value = p_txt(d$p), standard_error = d$se)[!is.na(l$pos), ]
  f <- sprintf("02_gwas_catalog_harmonised_%s.h.tsv.gz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "GWAS Catalog harmonised (hm_ columns, GRCh38)", "hm_ prefixes; empty CI columns") })
# 3 FinnGen (GRCh38): #chrom pos ref alt rsids ... pval mlogp beta sebeta af_alt; the effect allele is alt; X written as 23
both(function(d, role) { l <- to38(d$chr, d$pos); alt <- ifelse(l$rev, comp(d$ea), d$ea); ref <- ifelse(l$rev, comp(d$oa), d$oa)
  o <- data.frame(`#chrom` = sub("^X$", "23", d$chr), pos = l$pos, ref = ref, alt = alt, rsids = ifelse(seq_len(nrow(d)) %% 9 == 0, paste0(d$snp, ",rs999", seq_len(nrow(d))), d$snp),
                  nearest_genes = "GENE1", pval = p_txt(d$p), mlogp = -log10(pmax(d$p, 1e-300)), beta = d$beta, sebeta = d$se, af_alt = d$eaf, check.names = FALSE)[!is.na(l$pos), ]
  f <- sprintf("03_finngen_%s.tsv.gz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "FinnGen (GRCh38)", "effect allele = alt; some rsids hold two IDs") })
# 4 Neale lab UK Biobank (GRCh37, .tsv.bgz): variant chr:pos:ref:alt, beta for alt, no rsID column
both(function(d, role) {
  o <- data.frame(variant = paste(d$chr, d$pos, d$oa, d$ea, sep = ":"), minor_allele = ifelse(d$eaf < 0.5, d$ea, d$oa), minor_AF = pmin(d$eaf, 1 - d$eaf), low_confidence_variant = "false",
                  n_complete_samples = d$n, AC = round(2 * d$n * d$eaf), ytx = 0, beta = d$beta, se = d$se, tstat = d$beta / d$se, pval = p_txt(d$p))
  f <- sprintf("04_neale_ukb_%s.tsv.bgz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "Neale lab UK Biobank (GRCh37)", "variant chr:pos:ref:alt only; .bgz; minor allele columns") })
# 5 BOLT-LMM: ALLELE1 is the effect allele, A1FREQ its frequency
both(function(d, role) {
  o <- data.frame(SNP = d$snp, CHR = d$chr, BP = d$pos, GENPOS = 0, ALLELE1 = d$ea, ALLELE0 = d$oa, A1FREQ = d$eaf, INFO = 0.99, BETA = d$beta, SE = d$se, P_BOLT_LMM_INF = p_txt(d$p), P_BOLT_LMM = p_txt(d$p))
  f <- sprintf("05_bolt_lmm_%s.txt.gz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "BOLT-LMM", "ALLELE1/ALLELE0, two P columns") })
# 6 REGENIE: space separated, ALLELE1 is the effect allele, LOG10P (-log10 P), no P column
both(function(d, role) {
  o <- data.frame(CHROM = d$chr, GENPOS = d$pos, ID = d$snp, ALLELE0 = d$oa, ALLELE1 = d$ea, A1FREQ = d$eaf, INFO = 1, N = d$n, TEST = "ADD", BETA = d$beta, SE = d$se, CHISQ = (d$beta / d$se)^2,
                  LOG10P = -log10(pmax(d$p, 1e-300)), EXTRA = "NA")
  f <- sprintf("06_regenie_%s.regenie", tag(role)); wtab(o, f, sep = " "); reg(f, role, "REGENIE", "space separated; LOG10P only; ALLELE0/ALLELE1") })
# 7 PLINK 2 --glm: A1 is the effect allele and can be REF or ALT; covariate rows (TEST = AGE, SEX); logistic output gives OR and LOG(OR)_SE
both(function(d, role) { sw <- seq_len(nrow(d)) %% 2 == 0; ref <- ifelse(sw, d$ea, d$oa); alt <- ifelse(sw, d$oa, d$ea)
  base <- data.frame(`#CHROM` = d$chr, POS = d$pos, ID = d$snp, REF = ref, ALT = alt, A1 = d$ea, A1_FREQ = d$eaf, TEST = "ADD", OBS_CT = d$n, check.names = FALSE)
  o <- if (role == "outcome") data.frame(base, OR = exp(d$beta), `LOG(OR)_SE` = d$se, Z_STAT = d$beta / d$se, P = p_txt(d$p), ERRCODE = ".", check.names = FALSE)
       else data.frame(base, BETA = d$beta, SE = d$se, T_STAT = d$beta / d$se, P = p_txt(d$p), ERRCODE = ".", check.names = FALSE)
  cv <- o[1:4, ]; cv$TEST <- rep(c("AGE", "SEX"), 2); cv[[if (role == "outcome") "OR" else "BETA"]] <- if (role == "outcome") 1.01 else 0.01
  o <- rbind(o, cv); o <- o[order(o$ID, o$TEST != "ADD"), ]
  f <- sprintf("07_plink2_glm_%s.%s", tag(role), if (role == "outcome") "glm.logistic.hybrid" else "glm.linear"); wtab(o, f); reg(f, role, "PLINK 2 --glm", "A1 equals REF in half the rows; covariate rows; OR + LOG(OR)_SE for the outcome") })
# 8 SAIGE: Allele2 is the effect allele (AF_Allele2), Allele1 the reference
both(function(d, role) {
  o <- data.frame(CHR = d$chr, POS = d$pos, MarkerID = d$snp, Allele1 = d$oa, Allele2 = d$ea, AC_Allele2 = round(2 * d$n * d$eaf), AF_Allele2 = d$eaf, imputationInfo = 1, N = d$n,
                  BETA = d$beta, SE = d$se, Tstat = d$beta / d$se, `p.value` = p_txt(d$p), check.names = FALSE)
  f <- sprintf("08_saige_%s.txt", tag(role)); wtab(o, f, sep = " "); reg(f, role, "SAIGE", "effect allele = Allele2") })
# 9 METAL: lower-case alleles, Windows line endings, Direction column
both(function(d, role) {
  o <- data.frame(MarkerName = d$snp, Allele1 = tolower(d$ea), Allele2 = tolower(d$oa), Freq1 = d$eaf, FreqSE = 0.01, MinFreq = d$eaf, MaxFreq = d$eaf, Effect = d$beta, StdErr = d$se,
                  `P-value` = p_txt(d$p), Direction = "++-+", check.names = FALSE)
  f <- sprintf("09_metal_%s.tbl", tag(role)); wtab(o, f, eol = "\r\n"); reg(f, role, "METAL", "lower-case alleles; CRLF line endings") })
# 10 GIANT-style CSV with quotes
both(function(d, role) {
  o <- data.frame(SNP = d$snp, CHR = d$chr, POS = d$pos, Tested_Allele = d$ea, Other_Allele = d$oa, Freq_Tested_Allele_in_HRS = d$eaf, BETA = d$beta, SE = d$se, P = p_txt(d$p), N = d$n)
  f <- sprintf("10_giant_style_%s.csv", tag(role)); wtab(o, f, sep = ",", quote = TRUE); reg(f, role, "GIANT style (quoted CSV)", "Tested_Allele; frequency heading names the cohort") })
# 11 A paper's supplementary table (Excel): title and note rows, combined columns, typeset P values, footnote
both(function(d, role) {
  ptxt <- function(p) { e <- floor(log10(p)); m <- p / 10^e; ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.1f × 10−%d", m, -e)) }
  est <- if (role == "outcome") sprintf("%.3f (%.3f–%.3f)", exp(d$beta), exp(d$beta - 1.959964 * d$se), exp(d$beta + 1.959964 * d$se)) else sprintf("%.4f (%.4f)", d$beta, d$se)
  hd <- c("Locus", "Lead SNP", "Chr:position (hg19)", "Effect/other allele", "EAF", if (role == "outcome") "OR (95% CI)" else "Beta (SE)", "P")
  rows <- rbind(c(sprintf("Supplementary Table 2. Lead variants for %s", if (role == "outcome") "coronary heart disease" else "LDL cholesterol"), rep("", 6)),
                c("Effects per allele; P values two-sided.", rep("", 6)), hd,
                cbind(paste0("Locus ", seq_len(nrow(d))), d$snp, paste0(d$chr, ":", d$pos), paste0(d$ea, "/", d$oa), sprintf("%.3f", d$eaf), est, ptxt(d$p)),
                c("SNP, single-nucleotide polymorphism; EAF, effect allele frequency.", rep("", 6)))
  f <- sprintf("11_supplementary_table_%s.xlsx", tag(role)); write_xlsx(rows, file.path(dest, f)); reg(f, role, "Paper supplementary table (Excel)", "combined allele and estimate columns; typeset P values") })
# 12 European CSV: semicolons and decimal commas
both(function(d, role) {
  o <- data.frame(SNP = d$snp, CHR = d$chr, POS = d$pos, EA = d$ea, NEA = d$oa, EAF = d$eaf, BETA = d$beta, SE = d$se, P = d$p)
  f <- sprintf("12_european_%s.csv", tag(role)); wtab(o, f, sep = ";", dec = ","); reg(f, role, "European spreadsheet export", "semicolons; decimal commas") })
# 13 Messy text: a preamble, missing-value tokens, a duplicated row, an indel, a row without rsID, trailing spaces, Windows line endings
both(function(d, role) {
  o <- data.frame(rsid = d$snp, chr = paste0("chr", d$chr), pos = d$pos, effect_allele = d$ea, other_allele = d$oa, eaf = as.character(d$eaf), beta = as.character(d$beta), se = as.character(d$se), pval = p_txt(d$p))
  o$eaf[3] <- "NA"; o$se[5] <- "."; o$beta[7] <- paste0(o$beta[7], "  ")
  o <- rbind(o, o[2, ], data.frame(rsid = "rs999999901", chr = "chr1", pos = 1234567, effect_allele = "TA", other_allele = "T", eaf = "0.2", beta = "0.01", se = "0.01", pval = "0.3"))
  L <- c("## produced by an in-house pipeline", "## columns: see the README", paste(names(o), collapse = "\t"), apply(o, 1, paste, collapse = "\t"))
  f <- sprintf("13_messy_%s.txt", tag(role)); con <- file(file.path(dest, f), "wb"); writeLines(L, con, sep = "\r\n"); close(con); reg(f, role, "Messy in-house text", "preamble; NA tokens; duplicate; indel; CRLF") })
# 14 LDSC munged sumstats: Z and N only (A1 the effect allele); continuous exposure only (an outcome on a binary trait cannot be derived)
both(function(d, role) {
  o <- data.frame(SNP = d$snp, A1 = d$ea, A2 = d$oa, Z = d$beta / d$se, N = d$n)
  f <- sprintf("14_ldsc_sumstats_%s.sumstats.gz", tag(role)); wtab(o, f, gz = TRUE); reg(f, role, "LDSC munged (Z only)", if (role == "outcome") "binary trait: must be refused" else "Z scores; frequencies from the 1000 Genomes panel") })
# 15 Effect with a confidence interval and no SE
both(function(d, role) {
  o <- if (role == "outcome") data.frame(SNP = d$snp, EA = d$ea, OA = d$oa, OR = exp(d$beta), LCI = exp(d$beta - 1.959964 * d$se), UCI = exp(d$beta + 1.959964 * d$se), P = p_txt(d$p))
       else data.frame(SNP = d$snp, EA = d$ea, OA = d$oa, BETA = d$beta, L95 = d$beta - 1.959964 * d$se, U95 = d$beta + 1.959964 * d$se, P = p_txt(d$p))
  f <- sprintf("15_ci_no_se_%s.tsv", tag(role)); wtab(o, f); reg(f, role, "Estimate with 95% CI, no SE", "SE from the confidence interval") })

utils::write.csv(do.call(rbind, made), file.path(dest, "corpus.csv"), row.names = FALSE)
cat(sprintf("%d files written to %s\n", length(made), dest))
