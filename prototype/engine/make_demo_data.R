# Builds demo/demo_pair_iv.csv: SIMULATED instruments for the built-in demo pair. Not real GWAS data.
# Run once from the engine folder: Rscript make_demo_data.R   (the CSV is committed; this file documents how it was made)
set.seed(20261101); k <- 40
bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE)
sx <- runif(k, 0.005, 0.012); sy <- runif(k, 0.008, 0.02)
by <- 0.25 * bx + rnorm(k, 0.001, 0.006) + rnorm(k, 0, sy)
d <- data.frame(variant_id = sprintf("demo:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G",
                beta_x = round(bx, 6), se_x = round(sx, 6), eaf_x = round(runif(k, 0.1, 0.9), 3),
                beta_y = round(by, 6), se_y = round(sy, 6), eaf_y = NA_real_, p_x = NA_real_, p_y = NA_real_,
                n_x = 100000, n_y = 100000, harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
d$eaf_y <- d$eaf_x
d$p_x <- signif(2 * pnorm(-abs(d$beta_x / d$se_x)), 6); d$p_y <- signif(2 * pnorm(-abs(d$beta_y / d$se_y)), 6)
write.csv(d, "demo/demo_pair_iv.csv", row.names = FALSE, quote = FALSE)
