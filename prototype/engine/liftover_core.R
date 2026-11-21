# Round 77: GRCh38 to GRCh37 positions for imported files (MR Studio works on GRCh37). Base R only: reads a UCSC-format chain file
# and maps single positions, as liftOver does for a one-base interval. Positions that fall in a gap, or in more than one block, are not
# mapped (reported, never guessed). Variants on a reverse-strand block get their alleles complemented.
# Chain file: Ensembl assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz (Ensembl data are free to use without restriction),
# kept to the primary chromosomes 1-22, X, Y and MT (3,258 of 22,364 chains).
LIFTOVER_CHAIN <- list(file = "GRCh38_to_GRCh37.chain.gz", sha256 = "466a62474ab53ab08ee26279958c05e0133bce89aef3c0525c63f9b36a6ee0ee",
                       source = "Ensembl, https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz (sha256 of the full file 6bee5bfd...), primary chromosomes only")

liftover_chain_path <- function() {
  d <- getOption("mrgui.engine_dir", getOption("mrstudio.engine_dir", if (exists("args_dir")) args_dir else "prototype/engine"))
  file.path(d, "data", "liftover", LIFTOVER_CHAIN$file)
}

.lift_env <- new.env()
#' The chain as one table of ungapped blocks: target (GRCh38) chromosome and 0-based start, size, query (GRCh37) chromosome, 0-based start on
#' the query strand, the strand and the query chromosome length.
liftover_blocks <- function(path = liftover_chain_path()) {
  if (!is.null(.lift_env$b) && identical(.lift_env$path, path)) return(.lift_env$b)
  if (!file.exists(path)) stop("The GRCh38 to GRCh37 chain file is missing from this installation (", path, ").", call. = FALSE)
  if (!identical(digest::digest(file = path, algo = "sha256"), LIFTOVER_CHAIN$sha256)) stop("The GRCh38 to GRCh37 chain file of this installation is damaged (its checksum does not match). Reinstall MR Studio.", call. = FALSE)
  L <- readLines(gzfile(path)); L <- L[nzchar(L)]
  h <- startsWith(L, "chain"); cid <- cumsum(h); hd <- do.call(rbind, strsplit(L[h], " ", fixed = TRUE))
  bl <- strsplit(L[!h], "[ \t]+"); bc <- cid[!h]
  size <- as.numeric(vapply(bl, `[`, "", 1)); dt <- as.numeric(vapply(bl, function(x) if (length(x) > 1) x[2] else "0", "")); dq <- as.numeric(vapply(bl, function(x) if (length(x) > 2) x[3] else "0", ""))
  t0 <- as.numeric(hd[bc, 6]); q0 <- as.numeric(hd[bc, 11])
  # start of each block = chain start + the sizes and gaps of the blocks before it in the same chain
  prev_t <- stats::ave(size + dt, bc, FUN = function(v) c(0, cumsum(v)[-length(v)])); prev_q <- stats::ave(size + dq, bc, FUN = function(v) c(0, cumsum(v)[-length(v)]))
  b <- data.frame(tchr = hd[bc, 3], tstart = t0 + prev_t, size = size, qchr = hd[bc, 8], qstart = q0 + prev_q, qstrand = hd[bc, 10], qsize = as.numeric(hd[bc, 9]), stringsAsFactors = FALSE)
  b <- b[order(b$tchr, b$tstart), ]; rownames(b) <- NULL
  .lift_env$b <- b; .lift_env$path <- path; b
}

.chr_norm <- function(x) { x <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE)); x[x == "23"] <- "X"; x[x == "24"] <- "Y"; x[x %in% c("M", "25", "26")] <- "MT"; x }

#' Map 1-based GRCh38 positions to GRCh37. Returns chr, pos (NA when not mapped), strand ("+" or "-") and why a position was not mapped.
liftover_38_to_37 <- function(chr, pos, path = liftover_chain_path()) {
  b <- liftover_blocks(path); chr <- .chr_norm(chr); pos <- as.numeric(pos); n <- length(pos)
  out <- data.frame(chr = rep(NA_character_, n), pos = rep(NA_real_, n), strand = rep(NA_character_, n), why = rep("", n), stringsAsFactors = FALSE)
  out$why[is.na(pos) | pos < 1] <- "no position"
  for (cc in unique(chr[!is.na(pos)])) {
    i <- which(chr == cc & !is.na(pos) & pos >= 1); bb <- b[b$tchr == cc, , drop = FALSE]
    if (!nrow(bb)) { out$why[i] <- "chromosome not in GRCh38"; next }
    x <- pos[i] - 1; k <- findInterval(x, bb$tstart); en <- bb$tstart + bb$size
    hit <- k > 0; hit[hit] <- x[hit] < en[k[hit]]
    cmx <- cummax(en); amb <- k > 1; amb[amb] <- cmx[k[amb] - 1] > x[amb]   # an earlier block covers it too
    out$why[i[!hit]] <- "in a gap between GRCh38 and GRCh37"; out$why[i[hit & amb]] <- "maps to more than one place"
    ok <- hit & !amb; j <- k[ok]; off <- x[ok] - bb$tstart[j]; q <- bb$qstart[j] + off
    plus <- bb$qstrand[j] == "+"
    out$chr[i[ok]] <- bb$qchr[j]; out$strand[i[ok]] <- bb$qstrand[j]
    out$pos[i[ok]] <- ifelse(plus, q + 1, bb$qsize[j] - q)   # a reverse-strand block counts from the end of the chromosome
  }
  out
}

.complement <- function(a) chartr("ACGTacgt", "TGCAtgca", a)
