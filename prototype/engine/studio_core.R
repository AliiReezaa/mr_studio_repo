# MR Studio engine: study logic used by the screen. No Shiny code in this file, so it can be tested without a browser.
# Depends on the contracts v1.1 R code (loaded by load_contracts()).

CONTRACT_FILES <- c("canonical.R", "adapter.R", "security.R", "runner.R", file.path("adapters", "reference_adapters.R"),
                    "script.R", file.path("adapters", "reference_adapters_v11.R"), "study.R")

#' Find the contracts folder: MRGUI_CONTRACTS_DIR, the installed app (next to the engine), or the repository.
find_contracts_dir <- function(engine_dir) {
  cands <- c(Sys.getenv("MRGUI_CONTRACTS_DIR", ""), file.path(engine_dir, "..", "contracts"),
             file.path(engine_dir, "..", "..", "contracts"))
  for (d in cands[nzchar(cands)]) if (file.exists(file.path(d, "r", "study.R"))) return(normalizePath(d, winslash = "/"))
  stop("The contracts folder was not found. Looked in: ", paste(cands[nzchar(cands)], collapse = "; "))
}

load_contracts <- function(contracts_dir) {
  for (f in CONTRACT_FILES) source(file.path(contracts_dir, "r", f), local = .GlobalEnv)
  options(mrgui.schema = file.path(contracts_dir, "sql", "schema_v2.sql"))
  invisible(contracts_dir)
}

#' Why instruments clumped on the OpenGWAS server can share a P value (Round 31). OpenGWAS stores P values with a
#' limited number of digits; when the smallest P is itself shared, it is named as the stored floor (for example 0).
p_store_phrase <- function(ties) {
  if (isTRUE(ties$floor_tied) && is.finite(ties$floor %||% NA_real_)) sprintf("OpenGWAS stores P values with a limited number of digits, and the smallest ones as the same value, %s", format(ties$floor))
  else "OpenGWAS stores P values with a limited number of digits"
}

#' Default place for studies (decision D2). MRGUI_STUDIES_DIR overrides it (tests, support).
studies_root <- function() {
  r <- Sys.getenv("MRGUI_STUDIES_DIR", "")
  if (nzchar(r)) return(r)
  s <- settings_read()$studies_root
  if (is.character(s) && length(s) == 1L && nzchar(s)) s else file.path(path.expand("~"), "MR Studio")
}

# ---- user settings (theme, toggles, studies folder): one small JSON file, remembered per user (UI-1) ----------------
settings_dir <- function() {
  d <- Sys.getenv("MRGUI_CONFIG_DIR", "")
  if (nzchar(d)) d else tools::R_user_dir("MRStudio", "config")
}
SETTINGS_DEFAULT <- list(theme = "aurora", setup_done = FALSE, offline = FALSE, update = TRUE, awake = TRUE, notify = TRUE,
                         studies_root = "")
SETTINGS_THEMES <- c("aurora", "paper", "graphite", "ocean", "rose", "lilac", "sage", "peach")   # Sage and Peach added 2026-11-05
settings_read <- function() {
  f <- file.path(settings_dir(), "settings.json")
  s <- SETTINGS_DEFAULT
  if (file.exists(f)) {
    u <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) list())
    for (k in names(s)) if (!is.null(u[[k]]) && length(u[[k]]) == 1L && identical(class(u[[k]]), class(s[[k]]))) s[[k]] <- u[[k]]
  }
  if (!s$theme %in% SETTINGS_THEMES) s$theme <- "aurora"
  s
}
settings_set <- function(key, value) {
  if (!key %in% names(SETTINGS_DEFAULT)) stop("Unknown setting: ", key)
  if (!identical(class(value), class(SETTINGS_DEFAULT[[key]]))) stop("Wrong type for setting ", key)
  if (identical(key, "theme") && !value %in% SETTINGS_THEMES) stop("Unknown theme.")
  if (identical(key, "studies_root") && nzchar(value)) {
    prob <- study_path_problem(value); if (!is.null(prob)) stop("SYS_STUDY_FOLDER_REFUSED: ", prob, ". Choose a folder on this computer that is not synchronised.")
    check_path_length(value)
  }
  s <- settings_read(); s[[key]] <- value
  dir.create(settings_dir(), recursive = TRUE, showWarnings = FALSE)
  atomic_write_text(file.path(settings_dir(), "settings.json"), jsonlite::toJSON(s, auto_unbox = TRUE, pretty = TRUE))
  invisible(s)
}

#' Study folder names are what the user sees in the file manager: keep them readable, remove only what a file system refuses.
clean_study_name <- function(x) {
  x <- trimws(gsub("[\\\\/:*?\"<>|[:cntrl:]]", " ", x %||% ""))
  x <- gsub("\\s+", " ", x)
  x <- sub("^\\.+", "", x); x <- sub("[ .]+$", "", x)
  substr(x, 1L, 60L)
}
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

studies_list <- function(root = studies_root()) {
  if (!dir.exists(root)) return(character())
  d <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  d[file.exists(file.path(root, d, ".engine", "project.sqlite"))]
}

#' Create a study folder under `root` and open its database.
studio_new <- function(name, root = studies_root()) {
  nm <- clean_study_name(name)
  if (!nzchar(nm)) stop("Give the study a name.")
  dir <- file.path(root, nm)
  if (dir.exists(dir)) stop("A folder named \"", nm, "\" already exists in ", root, ". Choose another name, or open that study.")
  study <- study_create(dir)
  list(dir = study, con = db_open(file.path(study, ".engine")))
}

studio_open <- function(name, root = studies_root()) {
  dir <- file.path(root, name)
  if (!file.exists(file.path(dir, ".engine", "project.sqlite"))) stop("This folder is not a study: ", dir)
  study <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  con <- db_open(file.path(study, ".engine"))
  db_recover(con)    # a task left "running" by a crash goes back to pending (contracts C8)
  list(dir = study, con = con)
}

studio_pair_folders <- function(study) {
  d <- list.dirs(file.path(study, "03_pairs"), recursive = FALSE, full.names = FALSE)
  sort(d[grepl("^P[0-9]+_", d)])
}

#' The same exposure and outcome labels always map to the same pair folder (so a second run is a reuse, not a duplicate).
studio_pair_index <- function(study, exposure_label, outcome_label) {
  existing <- studio_pair_folders(study)
  key <- sub("^P[0-9]+_", "", pair_folder_name(1L, exposure_label, outcome_label))
  hit <- which(sub("^P[0-9]+_", "", existing) == key)
  if (length(hit)) return(as.integer(sub("^P([0-9]+)_.*", "\\1", existing[hit[1]])))
  length(existing) + 1L
}

# ---- demo data (SIMULATED; the real public GWAS download is a later milestone) ------------------------------------------
demo_pair <- function(engine_dir) {
  f <- file.path(engine_dir, "demo", "demo_pair_iv.csv")
  if (!file.exists(f)) stop("Demo data file is missing: ", f)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$mr_keep <- as.logical(d$mr_keep)
  list(exposure_label = "Demo exposure", outcome_label = "Demo outcome",
       input = list(kind = "pair_iv", data = d, exposure_id = "demo:simulated-exposure",
                    outcome_id = "demo:simulated-outcome", outcome_scale = "beta", exposure_scale = "beta"))
}

# ---- demo catalogue: SIMULATED traits so every screen can be tested before the real catalogue exists --------------------
DEMO_EXPOSURES <- data.frame(id = c("expA", "expB", "expC", "expD"), name = paste("Simulated exposure", LETTERS[1:4]),
                             n = c(173082, 681275, 757601, 204402), stringsAsFactors = FALSE)
DEMO_OUTCOMES <- data.frame(id = c("outA", "outB", "outC", "outD"), name = paste("Simulated outcome", LETTERS[1:4]),
                            n = c(184305, 392423, 440328, 487511), stringsAsFactors = FALSE)
DEMO_EFFECT <- matrix(c(0.30, 0.00, -0.15, 0.10,
                        0.00, 0.25,  0.00, 0.05,
                        0.20, 0.00,  0.00, 0.35,
                        0.00, 0.10, -0.25, 0.00), nrow = 4, byrow = TRUE, dimnames = list(DEMO_EXPOSURES$id, DEMO_OUTCOMES$id))

demo_catalogue <- function() {
  mk <- function(d) lapply(seq_len(nrow(d)), function(i) list(id = d$id[i], name = d$name[i], src = "Demo, simulated", n = format(d$n[i], big.mark = ","), anc = "EUR"))
  list(exposures = mk(DEMO_EXPOSURES), outcomes = mk(DEMO_OUTCOMES))
}
#' Datasets imported by the user, offered next to the simulated traits. An exposure must be marked as already independent, or set to be clumped on this computer (Round 19).
library_catalogue <- function(root = studies_root()) {
  recs <- tryCatch(library_list(root), error = function(e) list())
  mk <- function(r, role) list(id = paste0("lib:", r$trait), name = r$trait, src = if (identical(r$source$type, "opengwas")) "OpenGWAS" else "Imported", n = if (!is.null(r$facts$n) && nzchar(as.character(r$facts$n))) format(as.numeric(r$facts$n), big.mark = ",") else "unknown", anc = r$facts$ancestry %||% "?")
  list(exposures = lapply(Filter(function(r) isTRUE(r$usable_as_exposure), recs), mk, role = "exposure"),
       outcomes = lapply(Filter(function(r) isTRUE(r$usable_as_outcome), recs), mk, role = "outcome"))
}
catalogue_all <- function(root = studies_root()) {
  d <- demo_catalogue(); l <- library_catalogue(root)
  list(exposures = c(l$exposures, d$exposures), outcomes = c(l$outcomes, d$outcomes))
}
.seed_of <- function(x) sum(utf8ToInt(x) * seq_along(utf8ToInt(x))) %% 2147483647

#' One simulated pair of 40 instruments. The instruments belong to the exposure; the outcome noise belongs to the pair.
demo_input <- function(exp_id, out_id) {
  ei <- match(exp_id, DEMO_EXPOSURES$id); oi <- match(out_id, DEMO_OUTCOMES$id)
  if (is.na(ei) || is.na(oi)) stop("Unknown demo trait.")
  k <- 40L
  set.seed(.seed_of(exp_id))
  bx <- runif(k, 0.02, 0.12) * sample(c(-1, 1), k, replace = TRUE); sx <- runif(k, 0.005, 0.012); eaf <- round(runif(k, 0.1, 0.9), 3)
  set.seed(.seed_of(paste(exp_id, out_id)))
  sy <- runif(k, 0.008, 0.02); by <- DEMO_EFFECT[ei, oi] * bx + rnorm(k, 0.001, 0.006) + rnorm(k, 0, sy)
  d <- data.frame(variant_id = sprintf("demo:%d:A:G", seq_len(k) * 1000), ea = "A", oa = "G", beta_x = round(bx, 6), se_x = round(sx, 6),
                  eaf_x = eaf, beta_y = round(by, 6), se_y = round(sy, 6), eaf_y = eaf, p_x = signif(2 * pnorm(-abs(bx / sx)), 6),
                  p_y = signif(2 * pnorm(-abs(by / sy)), 6), n_x = DEMO_EXPOSURES$n[ei], n_y = DEMO_OUTCOMES$n[oi],
                  harmonise_action = "aligned", mr_keep = TRUE, stringsAsFactors = FALSE)
  list(exposure_label = DEMO_EXPOSURES$name[ei], outcome_label = DEMO_OUTCOMES$name[oi],
       input = list(kind = "pair_iv", data = d, exposure_id = paste0("demo:", exp_id), outcome_id = paste0("demo:", out_id), outcome_scale = "beta", exposure_scale = "beta"))
}

.vkey <- function(d) {
  if ("snp" %in% names(d) && any(nzchar(d$snp))) return(list(key = d$snp, kind = "rsID"))
  if (all(c("chr", "pos") %in% names(d))) return(list(key = paste0(d$chr, ":", d$pos), kind = "position"))
  stop("The dataset has neither variant IDs nor chromosome and position.")
}

#' A pair of imported datasets, made ready for the methods: instruments are the exposure variants with P below the threshold (the file must
#' already hold independent variants), looked up in the outcome and harmonised with TwoSampleMR::harmonise_data (action 2: palindromic
#' variants are aligned by allele frequency or dropped; action 1, keeping them, when both files are declared forward strand). Every instrument that is lost on the way is counted in the `ledger`.
#' Round 28: `ties` = "strongest" builds the tie sensitivity analysis instead: in every region where equal P values decided the lead
#' variant, the variant with the larger |beta/SE| is used. The ledger's `ties` entry says how many regions that changes.
library_pair_input <- function(exp_id, out_id, pthr = 5e-8, root = studies_root(), ties = c("server", "strongest")) {
  ties <- match.arg(ties)
  ex <- library_dataset(sub("^lib:", "", exp_id), root); ou <- library_dataset(sub("^lib:", "", out_id), root)
  if (!isTRUE(ex$recipe$usable_as_exposure)) stop("\"", ex$recipe$trait, "\" is neither marked as a set of independent instruments nor set to be clumped on this computer, so it cannot be an exposure. Import it again and answer the question about independent variants. Import a file that holds only clumped, genome-wide significant variants and answer \"Yes\" to the question about independent variants.")
  if (!isTRUE(ou$recipe$usable_as_outcome)) stop("\"", ou$recipe$trait, "\" was imported as an exposure only.")
  fe <- ou$recipe$source$for_exposure
  if (!is.null(fe) && !identical(fe, ex$recipe$trait)) stop("\"", ou$recipe$trait, "\" was downloaded only for the instruments of \"", fe, "\", so it cannot be paired with \"", ex$recipe$trait, "\". Download this outcome again for the exposure you want (Data library).")
  kx <- .vkey(ex$data); ky <- .vkey(ou$data)
  if (!identical(kx$kind, ky$kind)) stop("One dataset identifies variants by rsID and the other only by chromosome and position, so they cannot be matched. Import both with the same kind of ID.")
  X <- ex$data; X$key <- kx$key; Y <- ou$data; Y$key <- ky$key
  if (!("p" %in% names(X))) X$p <- 2 * stats::pnorm(-abs(X$beta / X$se))
  led <- list(candidates = nrow(X))
  X <- X[is.finite(X$p) & X$p < pthr, , drop = FALSE]; led$instruments <- nrow(X)
  sel <- NULL
  if (identical(ex$recipe$facts$independent, "clump")) {   # Round 19: clump the candidates on this computer (PLINK + 1000 Genomes)
    f <- ex$recipe$facts; num <- function(v, d) { v <- suppressWarnings(as.numeric(v %||% d)); if (length(v) == 1L && is.finite(v)) v else d }
    if (!identical(kx$kind, "rsID")) stop("\"", ex$recipe$trait, "\" identifies variants by position only. Clumping needs rsIDs, so it cannot be clumped on this computer.")
    if (!nrow(X)) stop("No variant in \"", ex$recipe$trait, "\" reaches P < ", format(pthr, scientific = TRUE), ", so there is nothing to clump.")
    cr2 <- num(f$clump_r2, 0.001); ckb <- num(f$clump_kb, 10000); cpop <- f$clump_pop %||% ld_pop_for(f$ancestry)
    sel <- ld_clump_cached(ex$dir, X$key, X$p, cr2, ckb, cpop)   # Round 28: PLINK's own order for equal P, as ieugwasr::ld_clump_local
    alt <- if (anyDuplicated(X$p[is.finite(X$p)])) ld_clump_cached(ex$dir, X$key, X$p, cr2, ckb, cpop, z = X$beta / X$se, tie_order = "strongest") else NULL
    led$ties <- ld_tie_summary(sel, alt)
    if (identical(ties, "strongest")) sel <- alt %||% sel
    X <- X[X$key %in% sel$variant_id[sel$status == "index"], , drop = FALSE]
    led$below_threshold <- led$instruments; led$instruments <- nrow(X)
    led$clump <- c(attr(sel, "settings"), list(clumped = sum(sel$status == "clumped"), not_in_panel = sum(sel$status == "not_in_panel"), kept = nrow(X)))
  }
  if (!is.null(ex$recipe$clump) && identical(ex$recipe$source$params$where, "local")) {   # Round 22: an OpenGWAS exposure clumped on this computer at download
    if (lib_table_exists(ex$dir, "instrument_selection")) {
      s0 <- lib_table_read(ex$dir, "instrument_selection")
      # Clumping keeps variants in order of P, so the clumped set of the candidates below a stricter threshold is exactly the
      # clumped set restricted to that threshold: the selection is cut, not clumped again.
      sel <- s0[is.finite(s0$p) & s0$p < pthr, , drop = FALSE]; rownames(sel) <- NULL
      st0 <- ex$recipe$clump[intersect(c("r2", "kb", "pop", "panel", "plink", "command", "tie_break", "tie_order"), names(ex$recipe$clump))]
      if (lib_table_exists(ex$dir, "instrument_selection_ties")) {   # Round 28: the same candidates with equal P broken by |beta/SE|
        a0 <- lib_table_read(ex$dir, "instrument_selection_ties")
        alt <- a0[is.finite(a0$p) & a0$p < pthr, , drop = FALSE]; rownames(alt) <- NULL
        led$ties <- ld_tie_summary(sel, alt)
        if (identical(ties, "strongest")) {
          if (!lib_table_exists(ex$dir, "candidates")) stop("The downloaded candidates of \"", ex$recipe$trait, "\" are missing, so the tie sensitivity analysis cannot be built. Download the exposure again.")
          cd <- lib_table_read(ex$dir, "candidates"); cd$key <- cd$snp
          X <- cd[cd$key %in% alt$variant_id[alt$status == "index"], intersect(names(X), names(cd)), drop = FALSE]
          sel <- alt; st0$tie_break <- "equal P ordered by |beta/SE|"; st0$tie_order <- "strongest"
        }
      } else if (identical(ties, "strongest")) stop("\"", ex$recipe$trait, "\" has no tied P values that change a lead variant, so there is no tie sensitivity analysis to run.")
      attr(sel, "settings") <- st0
      led$below_threshold <- nrow(sel); led$instruments <- nrow(X)
      led$clump <- c(attr(sel, "settings"), list(clumped = sum(sel$status == "clumped"), not_in_panel = sum(sel$status == "not_in_panel"), kept = sum(sel$status == "index"), at_download = TRUE))
    }
  }
  if (is.null(sel) && identical(ex$recipe$source$type, "opengwas")) {   # Round 28: clumped on the server; the ties cannot be checked here
    pt <- X$p[is.finite(X$p)]   # Round 31: floor_tied says whether the smallest P itself is shared (a storage floor such as 0 or 1e-200)
    led$ties <- list(server = TRUE, tied = .n_tied(pt), floor = if (length(pt)) min(pt) else NA_real_, floor_tied = length(pt) > 0 && sum(pt == min(pt)) > 1L)
  }
  if (identical(ties, "strongest")) {
    if (is.null(led$ties)) stop(if (is.null(sel)) paste0("MR Studio did not clump the instruments of \"", ex$recipe$trait, "\" (the file was declared to hold independent variants), so there are no ties to break.") else paste0("No two candidates of \"", ex$recipe$trait, "\" share a P value, so there is no tie sensitivity analysis to run."))
    if (isTRUE(led$ties$server)) stop("The instruments of \"", ex$recipe$trait, "\" were clumped on the OpenGWAS server, so MR Studio cannot see which tied variants were left out. Download the exposure with clumping on this computer to run the tie sensitivity analysis.")
    if (!isTRUE(led$ties$regions > 0)) stop("No tied P value changed an instrument of \"", ex$recipe$trait, "\", so the tie sensitivity analysis would repeat the primary analysis.")
    lk <- ou$recipe$lookup
    if (identical(ou$recipe$source$type, "opengwas") && !is.null(lk) && is.null(lk[["tie_alternatives"]]))
      stop("\"", ou$recipe$trait, "\" was downloaded before MR Studio looked up the tie alternatives. Download the outcome again (under another name) to run the tie sensitivity analysis.")
    led$tie_order <- "strongest"
  }
  if (nrow(X) < 2L) stop("Fewer than 2 variants in \"", ex$recipe$trait, "\" reach P < ", format(pthr, scientific = TRUE), ". At least 2 instruments are needed.")
  Y <- Y[Y$key %in% X$key, , drop = FALSE]; led$in_outcome <- nrow(Y); led$missing_in_outcome <- led$instruments - nrow(Y); led$by_proxy <- if ("proxy_snp" %in% names(Y)) sum(!is.na(Y$proxy_snp) & nzchar(Y$proxy_snp)) else 0L
  if (nrow(Y) < 2L) stop("Fewer than 2 instruments were found in \"", ou$recipe$trait, "\".")
  fmt <- function(D, type) {
    D2 <- data.frame(SNP = D$key, beta = D$beta, se = D$se, effect_allele = D$ea, other_allele = D$oa, stringsAsFactors = FALSE)
    if ("eaf" %in% names(D) && any(is.finite(D$eaf))) D2$eaf <- D$eaf
    if ("p" %in% names(D)) D2$pval <- D$p
    if ("n" %in% names(D)) D2$samplesize <- D$n
    suppressMessages(TwoSampleMR::format_data(D2, type = type, snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
      eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize"))   # always a name: format_data skips a column that is absent, but fails on NULL (Round 26)
  }
  # Round 26: palindromic variants are kept (action 1) only when the analyst declared, for both files, that alleles are on the forward strand.
  fwd <- identical(ex$recipe$facts$strand, "forward") && identical(ou$recipe$facts$strand, "forward")
  H <- suppressWarnings(suppressMessages(TwoSampleMR::harmonise_data(fmt(X, "exposure"), fmt(Y, "outcome"), action = if (fwd) 1 else 2)))
  if (!nrow(H)) stop("No instrument could be matched to the outcome.")
  oy <- stats::setNames(Y$ea, Y$key)
  flipped <- toupper(H$effect_allele.exposure) != toupper(oy[H$SNP])
  act <- ifelse(H$mr_keep, ifelse(flipped, "flipped", "aligned"), ifelse(H$palindromic %in% TRUE, "dropped_palindromic", "dropped_allele_mismatch"))
  col <- function(nm) if (nm %in% names(H)) H[[nm]] else rep(NA_real_, nrow(H))
  d <- data.frame(variant_id = H$SNP, ea = toupper(H$effect_allele.exposure), oa = toupper(H$other_allele.exposure), beta_x = H$beta.exposure, se_x = H$se.exposure,
                  eaf_x = col("eaf.exposure"), beta_y = H$beta.outcome, se_y = H$se.outcome, eaf_y = col("eaf.outcome"), p_x = col("pval.exposure"), p_y = col("pval.outcome"),
                  n_x = col("samplesize.exposure"), n_y = col("samplesize.outcome"), harmonise_action = act, mr_keep = H$mr_keep %in% TRUE, stringsAsFactors = FALSE)
  d <- d[order(d$p_x), , drop = FALSE]; rownames(d) <- NULL
  led$flipped <- sum(d$harmonise_action == "flipped"); led$dropped_palindromic <- sum(d$harmonise_action == "dropped_palindromic")
  led$dropped_mismatch <- sum(d$harmonise_action == "dropped_allele_mismatch"); led$kept <- sum(d$mr_keep); led$threshold <- pthr; led$id_kind <- kx$kind
  led$strand <- if (fwd) "forward" else "unknown"; led$harmonise_action <- if (fwd) 1L else 2L
  led$kept_palindromic <- sum(H$palindromic[match(d$variant_id, H$SNP)] %in% TRUE & d$mr_keep)
  if (led$kept < 2L) stop("Fewer than 2 instruments remain after harmonisation.")
  # Round 33: allele-frequency check of both datasets against the panel of their ancestry (reported, nothing removed)
  ec <- list(x = tryCatch(eaf_panel_check(d, d$eaf_x, ex$dir, ex$recipe$facts$ancestry, ex$recipe$trait), error = function(e) list(run = FALSE, why = conditionMessage(e), trait = ex$recipe$trait)),
             y = tryCatch(eaf_panel_check(d, d$eaf_y, ou$dir, ou$recipe$facts$ancestry, ou$recipe$trait), error = function(e) list(run = FALSE, why = conditionMessage(e), trait = ou$recipe$trait)))
  eaf_rows <- do.call(rbind, lapply(ec, `[[`, "rows")); led$eaf_check <- lapply(ec, function(z) z[setdiff(names(z), "rows")])
  sc <- if (identical(ou$recipe$facts$kind, "binary")) "log_or" else "beta"
  scx <- if (identical(ex$recipe$facts$kind, "binary")) "log_or" else "beta"
  list(exposure_label = ex$recipe$trait, outcome_label = ou$recipe$trait, simulated = FALSE, ledger = led, selection = sel, eaf_rows = eaf_rows,
       anc_x = ex$recipe$facts$ancestry %||% "?", anc_y = ou$recipe$facts$ancestry %||% "?", ukb_both = isTRUE(ex$recipe$facts$ukb) && isTRUE(ou$recipe$facts$ukb),
       source_x = ex$recipe$source$id %||% NULL, source_y = ou$recipe$source$id %||% NULL,
       binary = list(x = if (identical(scx, "log_or")) .case_counts(ex$recipe) else NULL, y = if (identical(sc, "log_or")) .case_counts(ou$recipe) else NULL),
       input = list(kind = "pair_iv", data = d, exposure_id = exp_id, outcome_id = out_id, outcome_scale = sc, exposure_scale = scx))
}

#' Save the clumping record of a pair (every candidate, kept or not, with its index variant) in its 01_data folder and the manifest.
save_instrument_selection <- function(study_dir, pair_folder, sel) {
  if (is.null(sel)) return(invisible(NULL))
  rel <- file.path("03_pairs", pair_folder, "01_data", "instrument_selection.csv")
  st <- attr(sel, "settings")
  out <- cbind(sel, clump_r2 = st$r2, clump_kb = st$kb, reference = paste(st$panel, st$pop), plink = st$plink %||% NA_character_)
  utils::write.csv(out, file.path(study_dir, rel), row.names = FALSE, na = "")
  manifest_register(study_dir, rel, role = "instrument_selection")
  invisible(rel)
}

#' Round 33: save the allele-frequency check of a pair (one row per variant and dataset) in its 01_data folder and the manifest.
save_eaf_check <- function(study_dir, pair_folder, rows) {
  if (is.null(rows) || !nrow(rows)) return(invisible(NULL))
  rel <- file.path("03_pairs", pair_folder, "01_data", "allele_frequency_check.csv")
  utils::write.csv(cbind(rows, threshold = EAF_FLIP_THRESHOLD, reference = LD_PANEL$name), file.path(study_dir, rel), row.names = FALSE, na = "")
  manifest_register(study_dir, rel, role = "allele_frequency_check")
  invisible(rel)
}

#' Numbers of cases and controls of a binary dataset: from the import facts, else from the OpenGWAS record (NA when unknown).
.case_counts <- function(rec) {
  num <- function(v) { v <- suppressWarnings(as.numeric(v %||% NA)); if (length(v) == 1L && is.finite(v) && v >= 1) v else NA_real_ }
  list(ncase = num(rec$facts$ncase %||% rec$source$ncase), ncontrol = num(rec$facts$ncontrol %||% rec$source$ncontrol))
}

#' The input the Steiger test reads: the pair input plus, for each binary side, its numbers of cases and controls and the population
#' prevalence entered by the analyst (columns ncase_*, ncontrol_*, prevalence_*). A continuous side gets no extra columns, so the
#' analysis file of a continuous pair is unchanged. `prevalence` is list(x = , y = ) (NULL or NA when not entered).
steiger_input <- function(p, prevalence = list()) {
  inp <- p$input; d <- inp$data
  for (s in c("x", "y")) {
    b <- p$binary[[s]]; if (is.null(b)) next
    pv <- suppressWarnings(as.numeric(prevalence[[s]] %||% NA)); if (length(pv) != 1L) pv <- NA_real_
    d[[paste0("ncase_", s)]] <- b$ncase; d[[paste0("ncontrol_", s)]] <- b$ncontrol; d[[paste0("prevalence_", s)]] <- pv
  }
  inp$data <- d; inp
}

#' Why the Steiger test cannot run for a pair whatever prevalence is entered (NULL when the data allow it), and which sides are binary.
steiger_readiness <- function(p) {
  probe <- steiger_input(p, list(x = 0.5, y = 0.5))   # a placeholder prevalence, only to test everything else
  v <- ADAPTER_STEIGER_V11$check(probe, list())
  list(bin_x = !is.null(p$binary$x), bin_y = !is.null(p$binary$y),
       why = if (identical(v$status, "block")) paste(v$messages, collapse = " ") else "")
}

#' One pair as the methods need it, from the simulated catalogue or from imported datasets (never a mixture).
pair_input <- function(exp_id, out_id, pthr = 5e-8, ties = "server") {
  le <- startsWith(exp_id, "lib:"); lo <- startsWith(out_id, "lib:")
  if (le != lo) stop("A simulated trait cannot be paired with an imported one. Choose two simulated traits or two imported datasets.")
  if (le) return(library_pair_input(exp_id, out_id, pthr, ties = ties))
  if (!identical(ties, "server")) stop("Simulated traits have no tied P values, so there is no tie sensitivity analysis to run.")
  p <- demo_input(exp_id, out_id); p$simulated <- TRUE; p$anc_x <- "EUR"; p$anc_y <- "EUR"; p$ledger <- NULL; p
}

#' What the engine can say about a pair before running it (Methods and Check screens). Real numbers from the data table.
pair_facts <- function(exp_id, out_id, pthr = 5e-8) {
  p <- pair_input(exp_id, out_id, pthr); d <- p$input$data; kept <- d[d$mr_keep, ]
  Fs <- (kept$beta_x / kept$se_x)^2
  pal <- sum((kept$ea == "A" & kept$oa == "T") | (kept$ea == "T" & kept$oa == "A") | (kept$ea == "C" & kept$oa == "G") | (kept$ea == "G" & kept$oa == "C"))
  list(ex_id = exp_id, out_id = out_id, exposure = p$exposure_label, outcome = p$outcome_label, instruments = nrow(kept), mean_f = round(mean(Fs), 1), min_f = round(min(Fs), 1),
       n_weak = sum(Fs < 10), palindromic = pal, scale = p$input$outcome_scale,
       has_n = nrow(kept) > 0 && all(c("n_x", "n_y") %in% names(kept)) && all(is.finite(kept$n_x) & is.finite(kept$n_y) & kept$n_x > 2 & kept$n_y > 2),
       simulated = isTRUE(p$simulated), anc_x = p$anc_x, anc_y = p$anc_y, ledger = p$ledger, ukb_both = isTRUE(p$ukb_both), og_x = !is.null(p$source_x), og_y = !is.null(p$source_y),
       scale_x = p$input$exposure_scale %||% "beta", steiger = steiger_readiness(p),
       i2gx = if (exists("feas_i2gx")) feas_i2gx(kept$beta_x, kept$se_x) else NA_real_)
}

METHOD_CHOICES <- c("Inverse-variance weighted, multiplicative random effects" = "ivw_mre", "Debiased IVW (dIVW)" = "divw", "MR-RAPS (robust adjusted profile score)" = "raps", "MR-Egger regression" = "egger", "Weighted median" = "wm", "Simple median" = "smed", "Weighted mode" = "mode", "Simple mode" = "smode", "MR-PRESSO" = "presso", "Radial IVW" = "radial", "Radial MR-Egger" = "radegger", "MR Steiger directionality" = "steiger")

param_adapters <- function() list(ivw_mre = ADAPTER_IVW_MRE_V11, divw = ADAPTER_DIVW_V11, raps = ADAPTER_RAPS_V11, egger = ADAPTER_EGGER_V11, wm = ADAPTER_WM_V11, smed = ADAPTER_SMED_V11, mode = ADAPTER_MODE_V11, smode = ADAPTER_SMODE_V11, presso = ADAPTER_PRESSO_V11, radial = ADAPTER_RADIAL_V11, radegger = ADAPTER_RADEGGER_V11, steiger = ADAPTER_STEIGER_V11, loo = ADAPTER_LOO_V11)

studio_run <- function(h, exposure_label, outcome_label, input, methods = c("ivw_mre", "egger"), ci_level = 0.95, idx = NULL, params = NULL, steiger_in = NULL) {
  params <- c(params[setdiff(names(params), "ci_level")], list(ci_level = if (!is.null(params$ci_level)) params$ci_level else ci_level))
  if (!length(methods)) stop("Choose at least one method.")
  exposure_label <- trimws(exposure_label); outcome_label <- trimws(outcome_label)
  if (!nzchar(exposure_label) || !nzchar(outcome_label)) stop("Give the exposure and the outcome a name.")
  ad <- param_adapters()
  if (is.null(steiger_in)) steiger_in <- input   # a continuous pair needs nothing extra; a binary one gets steiger_input()
  if ("steiger" %in% methods) { v <- ad$steiger$check(steiger_in, list()); if (identical(v$status, "block")) stop(paste(v$messages, collapse = " ")) }   # the script would stop on missing sample sizes but could not tell a binary outcome
  spec <- lapply(methods, function(m) list(adapter = ad[[m]], input = if (identical(m, "steiger")) steiger_in else input,
                                           params = params[intersect(names(params), names(ad[[m]]$params))]))
  if (is.null(idx)) idx <- studio_pair_index(h$dir, exposure_label, outcome_label)
  r <- run_script_task(h$con, h$dir, idx, exposure_label, outcome_label, spec)
  r$pair_folder <- basename(r$pair_dir)
  r
}

METHOD_NAMES <- c(ivw_mre = "IVW (multiplicative random effects)", divw = "Debiased IVW", raps = "MR-RAPS", egger = "MR-Egger", wm = "Weighted median", smed = "Simple median", mode = "Weighted mode", smode = "Simple mode", presso = "MR-PRESSO", radial = "Radial IVW", radegger = "Radial MR-Egger", steiger = "MR Steiger directionality", loo = "Leave-one-out and single-variant analyses")
IVW_MODEL_NAMES <- c(mre_floor = "IVW (multiplicative random effects)", mre = "IVW (multiplicative random effects, no fixed-effect floor)", fixed = "IVW (fixed effect)")

#' The settings a pair was run with, read back from the settings lines of its scripts (the scripts are the record).
pair_settings <- function(study, pair_folder) {
  out <- list(ci_level = 0.95, min_f = 0, ivw_model = "mre", nboot = 1000, phi = 1, nbdist = 1000, signif = 0.05, alpha = 0.05, weights = "1", tol = 1e-4, r_xxo = 1, r_yyo = 1, dispersion = "balanced", raps_dispersion = "balanced", loss = "huber")
  for (f in Sys.glob(file.path(study, "03_pairs", pair_folder, "02_scripts", "[0-9]*_*.R"))) {
    for (k in names(out)) {
      ln <- grep(paste0("^", k, " <- "), readLines(f, warn = FALSE), value = TRUE)
      if (length(ln)) { v <- sub(paste0("^", k, " <- "), "", ln[1]); out[[k]] <- if (k %in% c("ivw_model", "weights", "dispersion", "raps_dispersion", "loss")) gsub("\"", "", v) else as.numeric(v) }
    }
  }
  out
}

#' The name of a method as it was actually run, read from the settings line of its script (so a reopened study is labelled correctly).
method_label <- function(study, pair_folder, m) {
  if (!identical(m, "ivw_mre")) return(METHOD_NAMES[[m]])
  f <- Sys.glob(file.path(study, "03_pairs", pair_folder, "02_scripts", "*_ivw_mre.R"))
  if (!length(f)) return(METHOD_NAMES[[m]])
  ln <- grep("^ivw_model <- ", readLines(f[1], warn = FALSE), value = TRUE)
  k <- if (length(ln)) gsub("^ivw_model <- \"|\".*$", "", ln[1]) else "mre"
  if (k %in% names(IVW_MODEL_NAMES)) IVW_MODEL_NAMES[[k]] else METHOD_NAMES[[m]]
}

#' Estimates of a pair as one table, read from the delivered CSV files (works for a reopened study too).
pair_results <- function(study, pair_folder) {
  tdir <- file.path(study, "03_pairs", pair_folder, "04_tables")
  out <- list()
  for (m in names(METHOD_NAMES)) {
    f <- file.path(tdir, paste0(m, "_estimates.csv"))
    if (!file.exists(f)) next
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (!nrow(d)) next   # a diagnostic-only method (Steiger) has a header-only estimates file
    keep <- intersect(c("estimand", "estimate", "se", "ci_lower", "ci_upper", "ci_level", "p", "n_ivs_used"), names(d))
    out[[m]] <- cbind(Method = method_label(study, pair_folder, m), d[, keep, drop = FALSE])
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}

fmt_num <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "fg", flag = "#"))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, formatC(p, digits = 2, format = "e"), formatC(p, digits = 3, format = "f")))
fmt_bytes <- function(b) ifelse(b < 1024, paste0(b, " B"), ifelse(b < 1024^2, sprintf("%.1f KB", b / 1024), sprintf("%.1f MB", b / 1024^2)))

#' The results table as the screen shows it: readable headers, interval in one column, P values that never print as 0.
format_results <- function(res) {
  est <- c(causal_effect = "Causal effect", pleiotropy_intercept = "Pleiotropy intercept")
  lvl <- if ("ci_level" %in% names(res)) res$ci_level else rep(0.95, nrow(res))
  data.frame(Method = res$Method, Estimate = fmt_num(res$estimate), `Std. error` = fmt_num(res$se),
             Interval = ifelse(is.na(res$ci_lower), "", sprintf("%s to %s (%g%%)", fmt_num(res$ci_lower), fmt_num(res$ci_upper), 100 * lvl)),
             `P value` = fmt_p(res$p), Variants = as.character(res$n_ivs_used),
             Estimand = unname(ifelse(res$estimand %in% names(est), est[res$estimand], res$estimand)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

pair_replay <- function(study, pair_folder) {
  f <- file.path(study, "03_pairs", pair_folder, "07_replay_check", "replay_check.json")
  if (!file.exists(f)) return(NULL)
  jsonlite::fromJSON(f)
}

study_files <- function(study) {
  m <- manifest_read(study)
  if (!nrow(m)) return(m)
  m <- m[!duplicated(m$relative_path, fromLast = TRUE), c("relative_path", "role", "bytes"), drop = FALSE]
  rownames(m) <- NULL
  m
}

#' Run a program that belongs to the user (an editor, RStudio) without the engine's own R settings. The bundled R sets R_HOME and
#' library paths; a program started from here would otherwise inherit them and could pick up the bundled R instead of its own.
.with_user_env <- function(expr) {
  nm <- c("R_HOME", "R_SHARE_DIR", "R_INCLUDE_DIR", "R_DOC_DIR", "R_LIBS", "R_LIBS_SITE", "R_LIBS_USER", "R_ARCH", "R_PLATFORM", "R_PAPERSIZE", "R_ENVIRON", "R_PROFILE")
  old <- Sys.getenv(nm, unset = NA); on.exit({ for (n in nm) if (!is.na(old[[n]])) do.call(Sys.setenv, stats::setNames(list(old[[n]]), n)) }, add = TRUE)
  Sys.unsetenv(nm); force(expr)
}

#' Open a folder in Finder or Explorer. Only paths inside the studies root are opened.
open_in_file_manager <- function(path, root = studies_root()) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  r <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!startsWith(p, paste0(r, "/")) && !identical(p, r)) stop("Only study folders can be opened from here.")
  if (.Platform$OS.type == "windows") utils::shell.exec(normalizePath(p))
  else .with_user_env(system2(if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open", shQuote(p), wait = FALSE))
  invisible(TRUE)
}

#' Open a script in the computer's plain text editor (TextEdit, Notepad). Used when RStudio is not installed. The script is not
#' handed to whatever program is registered for .R files, because on many computers that is R itself, which cannot edit it well.
open_in_text_editor <- function(path, root = studies_root()) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  r <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!startsWith(p, paste0(r, "/"))) stop("Only scripts inside the studies folder can be opened from here.")
  sys <- Sys.info()[["sysname"]]
  .with_user_env(
    if (sys == "Darwin") system2("open", c("-t", shQuote(p)), wait = FALSE)
    else if (.Platform$OS.type == "windows") system2("notepad.exe", shQuote(normalizePath(p, winslash = "\\")), wait = FALSE)
    else system2("xdg-open", shQuote(p), wait = FALSE))
  invisible(p)
}


# ---- UI state kept in the study database (plan, pairs) --------------------------------------------------------------------
ui_state_get <- function(con, key) {
  r <- DBI::dbGetQuery(con, "SELECT value_json FROM ui_state WHERE key = ?", params = list(key))
  if (!nrow(r)) NULL else jsonlite::fromJSON(r$value_json, simplifyVector = FALSE)
}
ui_state_set <- function(con, key, value) {
  DBI::dbExecute(con, "INSERT OR REPLACE INTO ui_state(key, value_json, updated_at) VALUES (?, ?, ?)",
                 params = list(key, as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null")), format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  invisible(TRUE)
}

# ---- Sample overlap: the analyst's own statement about the two cohorts (not a calculation) -------------------------------------
# Overlap between the exposure and outcome studies cannot be calculated from the summary statistics the app reads, so the Check screen
# asks the analyst and records the answer in the study, where the reports read it. The answer is a declaration and is reported as one.
OVERLAP_ANSWERS <- c("none", "some", "unknown")
overlap_key <- function(exp_id, out_id) paste(exp_id, out_id, sep = " | ")
overlap_get <- function(h) { ov <- ui_state_get(h$con, "overlap"); if (is.null(ov) || !is.list(ov)) list() else ov }
overlap_set <- function(h, exp_id, out_id, answer, pct = NULL, note = "") {
  if (!nzchar(exp_id) || !nzchar(out_id)) stop("The pair is missing.")
  if (!(answer %in% OVERLAP_ANSWERS)) stop("Choose one of: the studies do not share people, they share some people, or I do not know.")
  pct <- suppressWarnings(as.numeric(pct)); if (length(pct) != 1L || !is.finite(pct) || answer != "some") pct <- NA_real_
  if (!is.na(pct) && (pct < 0 || pct > 100)) stop("The share of overlapping people must be between 0 and 100 percent.")
  note <- trimws(substr(gsub("[[:cntrl:]]+", " ", as.character(note %||% "")), 1, 500))
  ov <- overlap_get(h)
  ov[[overlap_key(exp_id, out_id)]] <- list(exposure = exp_id, outcome = out_id, answer = answer, pct = if (is.na(pct)) NULL else pct, note = note,
                                           saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  rel <- "01_study_plan/sample_overlap.json"
  atomic_write_text(file.path(h$dir, rel), jsonlite::toJSON(unname(ov), auto_unbox = TRUE, pretty = TRUE, null = "null"))
  manifest_register(h$dir, rel, role = "sample_overlap")
  ui_state_set(h$con, "overlap", ov)
  invisible(ov)
}

#' Save the question page of a new study (what the user typed) as a readable file and in the database.
study_plan_save <- function(h, plan) {
  plan$saved_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  rel <- "01_study_plan/study_plan.json"
  atomic_write_text(file.path(h$dir, rel), jsonlite::toJSON(plan, auto_unbox = TRUE, pretty = TRUE))
  manifest_register(h$dir, rel, role = "study_plan")
  ui_state_set(h$con, "plan", plan)
  invisible(plan)
}

age_label <- function(t) {
  d <- as.numeric(difftime(Sys.time(), t, units = "days"))
  if (d < 1) "today" else if (d < 2) "yesterday" else if (d < 14) paste(floor(d), "days ago") else if (d < 60) paste(floor(d / 7), "weeks ago") else format(t, "%b %Y")
}

#' Recent studies for Home and Studies: status comes from the files (results present means complete).
studies_overview <- function(root = studies_root()) {
  nm <- studies_list(root)
  out <- lapply(nm, function(n) {
    dir <- file.path(root, n)
    pf <- list.dirs(file.path(dir, "03_pairs"), recursive = FALSE, full.names = FALSE); pf <- pf[grepl("^P[0-9]+_", pf)]
    has_res <- any(vapply(pf, function(p) length(Sys.glob(file.path(dir, "03_pairs", p, "04_tables", "*_estimates.csv"))) > 0, logical(1)))
    nmeth <- length(unique(unlist(lapply(pf, function(p) sub("_estimates\\.csv$", "", basename(Sys.glob(file.path(dir, "03_pairs", p, "04_tables", "*_estimates.csv"))))))))
    mt <- file.mtime(file.path(dir, ".engine", "project.sqlite"))
    list(name = n, status = if (has_res) "COMPLETE" else "DRAFT", pairs = length(pf), methods = nmeth, when = age_label(mt), mtime = as.numeric(mt))
  })
  out[order(-vapply(out, function(x) x$mtime, 0))]
}

#' Everything the Results screen needs for one pair, read back from the delivered files.
pair_payload <- function(study, pair_folder) {
  pd <- file.path(study, "03_pairs", pair_folder)
  rd <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  est <- list(); dg <- list()
  for (m in names(METHOD_NAMES)) {
    e <- rd(file.path(pd, "04_tables", paste0(m, "_estimates.csv")))
    if (!is.null(e)) est[[m]] <- lapply(seq_len(nrow(e)), function(i) c(list(method = m, method_name = method_label(study, pair_folder, m)), as.list(e[i, , drop = FALSE])))
    g <- rd(file.path(pd, "04_tables", paste0(m, "_diagnostics.csv")))
    if (!is.null(g)) dg[[m]] <- lapply(seq_len(nrow(g)), function(i) c(list(method = m), as.list(g[i, , drop = FALSE])))
  }
  v <- rd(Sys.glob(file.path(pd, "01_data", "analysis_data_*.csv"))[1])
  vars <- if (!is.null(v)) lapply(seq_len(nrow(v)), function(i) list(id = v$variant_id[i], bx = v$beta_x[i], by = v$beta_y[i], sx = v$se_x[i], sy = v$se_y[i], keep = isTRUE(as.logical(v$mr_keep[i])))) else list()
  out <- list(folder = pair_folder, estimates = unname(unlist(est, recursive = FALSE)), diagnostics = unname(unlist(dg, recursive = FALSE)), variants = vars, settings = pair_settings(study, pair_folder),
       replay = pair_replay(study, pair_folder))
  if (exists("feasibility_results")) out$advice <- feasibility_results(out)   # rules that need results (R25 to R27)
  out
}

folder_size_bytes <- function(path) {
  f <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  if (!length(f)) 0 else sum(file.size(f), na.rm = TRUE)
}

#' Self-check for Settings > Run diagnostics. Each item: name, ok, detail.
run_diagnostics <- function(root = studies_root(), net = FALSE) {
  chk <- function(name, expr) { r <- tryCatch(expr, error = function(e) list(ok = FALSE, detail = conditionMessage(e))); c(list(name = name), r) }
  Filter(Negate(is.null), list(
    chk("R engine", list(ok = TRUE, detail = R.version.string)),
    chk("TwoSampleMR package", list(ok = requireNamespace("TwoSampleMR", quietly = TRUE), detail = if (requireNamespace("TwoSampleMR", quietly = TRUE)) paste("version", as.character(utils::packageVersion("TwoSampleMR"))) else "not installed")),
    chk("Script runner (Rscript)", { r <- default_rscript(); list(ok = file.exists(r), detail = if (file.exists(r)) "found" else paste("missing:", r)) }),
    chk("Studies folder is writable", { dir.create(root, recursive = TRUE, showWarnings = FALSE); f <- tempfile(".mrstudio-check-", root); ok <- tryCatch({ writeLines("x", f); TRUE }, error = function(e) FALSE, warning = function(w) FALSE); unlink(f); list(ok = ok, detail = root) }),
    chk("Studies folder is not synchronised", { p <- study_path_problem(root); list(ok = is.null(p), detail = if (is.null(p)) "local folder" else p) }),
    chk("PLINK (clumping on this computer)", { v <- if (exists("plink_version")) plink_version() else NA_character_; list(ok = !is.na(v), detail = if (is.na(v)) "not found in this installation" else v) }),
    if (isTRUE(net)) chk("Internet access to OpenGWAS", ogw_ping())
  ))
}

# ---- Open in RStudio (FR-86) and Compare with default (FR-66) ------------------------------------------------------------
#' Where RStudio is installed, or NA. Checked in the usual places only.
rstudio_path <- function() {
  sys <- Sys.info()[["sysname"]]
  cand <- if (sys == "Darwin") c("/Applications/RStudio.app", "/Applications/Posit/RStudio.app", file.path(Sys.getenv("HOME"), "Applications", "RStudio.app"))
          else if (sys == "Windows") c(file.path(Sys.getenv("ProgramFiles"), "RStudio", "rstudio.exe"), file.path(Sys.getenv("ProgramFiles"), "Posit", "RStudio", "rstudio.exe"), file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "RStudio", "rstudio.exe"))
          else unname(Sys.which("rstudio"))
  hit <- cand[nzchar(cand) & file.exists(cand)]
  if (length(hit)) hit[[1]] else NA_character_
}

#' Path of the verification script of a pair (VERIFY_IN_R.R). The files are written when missing, so pair folders made by older builds work too.
verify_script_path <- function(study, pair_folder, root = studies_root()) {
  pd <- file.path(study, "03_pairs", pair_folder)
  if (!file.exists(file.path(pd, "02_scripts", "00_run_all.R"))) stop("This pair has no script yet. Run the analysis first.")
  p <- normalizePath(pd, winslash = "/", mustWork = TRUE); r <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!startsWith(p, paste0(r, "/"))) stop("Only scripts inside the studies folder can be opened from here.")
  ensure_verify_files(p)
  list(script = file.path(p, "VERIFY_IN_R.R"), pair = p)
}

#' Open the check script of a pair (VERIFY_IN_R.R) in RStudio. It finds its own folder and runs on a copy, so nothing is set, copied or changed by the user.
open_in_rstudio <- function(study, pair_folder, root = studies_root()) {
  v <- verify_script_path(study, pair_folder, root)
  rs <- rstudio_path()
  if (is.na(rs)) return(list(opened = FALSE, script = v$script, pair = v$pair))
  .with_user_env(if (Sys.info()[["sysname"]] == "Darwin") system2("open", c("-a", shQuote(rs), shQuote(v$script)), wait = FALSE)
                 else system2(rs, shQuote(v$script), wait = FALSE))
  list(opened = TRUE, script = v$script, pair = v$pair)
}

#' The interval a result would have had at the default level (95%), recomputed from the saved estimate and standard error with the
#' same interval rule the method used. Only used for the Compare with default panel; the saved results are never changed.
ci_default_rows <- function(estimates, level = 0.95) {
  out <- list()
  for (e in estimates) {
    if (!identical(e$estimand, "causal_effect") || is.na(e$se)) next
    k <- if (e$method %in% c("radial", "presso")) 1 else 2   # degrees of freedom: variants used minus 2 (Egger types) or minus 1 (radial IVW, MR-PRESSO)
    q <- if (identical(e$ci_method, "wald_t")) stats::qt(1 - (1 - level) / 2, df = max(1, e$n_ivs_used - k)) else stats::qnorm(1 - (1 - level) / 2)
    out[[length(out) + 1L]] <- list(method = e$method, lo = e$estimate - q * e$se, hi = e$estimate + q * e$se, level = level)
  }
  out
}

#' Descriptive table of the instruments of a pair for the outlier review. Nothing here is an estimator: the single-variant ratio is
#' beta_y / beta_x, the weight share is the first-order inverse-variance weight (beta_x^2 / se_y^2) as a share of the total, and F is
#' (beta_x / se_x)^2. Detection of outliers is left to the packages named in the spec (not built yet).
instrument_table <- function(exp_id, out_id, pthr = 5e-8) {
  d <- pair_input(exp_id, out_id, pthr)$input$data; d <- d[d$mr_keep, ]
  w <- d$beta_x^2 / d$se_y^2
  data.frame(id = d$variant_id, bx = d$beta_x, by = d$beta_y, ratio = d$beta_y / d$beta_x, share = w / sum(w), f = (d$beta_x / d$se_x)^2, stringsAsFactors = FALSE)
}

#' Sensitivity branch: the same pair analysed again without the chosen variants, as its own pair folder (the primary results stay as they are).
studio_run_without <- function(h, exp_id, out_id, remove, methods = c("ivw_mre", "egger"), ci_level = 0.95, idx, params = NULL, pthr = 5e-8) {
  inp <- pair_input(exp_id, out_id, pthr)
  d <- inp$input$data; hit <- d$variant_id %in% remove
  if (!all(remove %in% d$variant_id)) stop("Some of the chosen variants are not instruments of this pair.")
  d$mr_keep[hit] <- FALSE; inp$input$data <- d
  kept <- sum(d$mr_keep)
  if (kept < 2L) stop("At least 2 instruments must remain.")
  methods <- setdiff(methods, c("presso", "radial", "radegger", "steiger", "loo"))   # detectors are not repeated on a reduced set: they are what suggests the sensitivity analysis
  if (kept < 3L) methods <- setdiff(methods, c("egger", "wm", "smed", "mode", "smode"))
  lab <- sprintf("%s [without %d]", inp$exposure_label, length(remove))
  r <- studio_run(h, lab, inp$outcome_label, inp$input, methods, ci_level, idx = idx, params = params)
  r$kept <- kept; r$label <- lab; r
}


#' Tie sensitivity analysis (Round 28): the same pair analysed again with the variant of larger |beta/SE| as the lead of every region
#' where equal P values decided it, as its own pair folder (the primary results stay as they are). Detectors are not repeated.
studio_run_ties <- function(h, exp_id, out_id, methods = c("ivw_mre", "egger"), ci_level = 0.95, idx, params = NULL, pthr = 5e-8) {
  inp <- pair_input(exp_id, out_id, pthr, ties = "strongest")
  methods <- setdiff(methods, c("presso", "radial", "radegger", "steiger", "loo"))
  kept <- sum(inp$input$data$mr_keep)
  if (kept < 3L) methods <- setdiff(methods, c("egger", "wm", "smed", "mode", "smode"))
  if (!length(methods)) methods <- "ivw_mre"
  lab <- sprintf("%s [ties: strongest]", inp$exposure_label)
  r <- studio_run(h, lab, inp$outcome_label, inp$input, methods, ci_level, idx = idx, params = params)
  r$kept <- kept; r$label <- lab; r$ties <- inp$ledger$ties; r$selection <- inp$selection; r$ledger <- inp$ledger; r
}

#' Compare with default settings, for every setting (FR-66). The same pair is analysed again with the default of each setting a method
#' offers, in a throw-away study folder, and set beside the saved result. A method whose only changed setting is the confidence level is
#' not run again: its estimate and standard error cannot change, so its default interval is recomputed from them.
#' The saved results are never touched.
compare_defaults <- function(study, pair_folder, exp_id, out_id, methods, pthr = 5e-8) {
  methods <- setdiff(methods, "steiger")   # no estimate to compare
  st <- pair_settings(study, pair_folder); ads <- param_adapters()
  dflt <- function(m, k) ads[[m]]$params[[k]]$default
  same <- function(a, b) if (is.character(b)) identical(as.character(a), b) else isTRUE(all.equal(as.numeric(a), as.numeric(b)))
  shown <- function(m, k, v) { sp <- ads[[m]]$params[[k]]; if (identical(sp$type, "enum") && !is.null(sp$value_labels)) unname(sp$value_labels[[as.character(v)]]) else format(v) }
  pay <- pair_payload(study, pair_folder); rows <- list(); rerun <- character(); diffs <- list()
  for (m in methods) {
    if (is.null(ads[[m]])) next
    d <- Filter(function(k) !identical(k, "ci_level") && !is.null(st[[k]]) && !same(st[[k]], dflt(m, k)), names(ads[[m]]$params))
    diffs[[m]] <- vapply(d, function(k) sprintf("%s: %s (default %s)", ads[[m]]$params[[k]]$label, shown(m, k, st[[k]]), shown(m, k, dflt(m, k))), "")
    if (length(d)) rerun <- c(rerun, m)
  }
  ci_changed <- !same(st$ci_level, 0.95)
  if (length(rerun)) {
    root <- tempfile("compare-"); dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    hs <- studio_new("compare", root); on.exit(try(DBI::dbDisconnect(hs$con), silent = TRUE), add = TRUE, after = FALSE)
    inp <- pair_input(exp_id, out_id, pthr)
    r <- studio_run(hs, inp$exposure_label, inp$outcome_label, inp$input, rerun, ci_level = 0.95, idx = 1L, params = list(ci_level = 0.95))
    dp <- pair_payload(hs$dir, r$pair_folder)
  }
  for (e in pay$estimates) {
    if (!identical(e$estimand, "causal_effect")) next
    m <- e$method
    if (m %in% rerun) { f <- Filter(function(x) identical(x$method, m) && identical(x$estimand, "causal_effect"), dp$estimates)
      if (!length(f)) next; f <- f[[1]]
      rows[[length(rows) + 1L]] <- list(method = m, ran = TRUE, estimate = f$estimate, lo = f$ci_lower, hi = f$ci_upper, p = f$p, n = f$n_ivs_used, what = as.list(unname(diffs[[m]])))
    } else {
      cd <- ci_default_rows(list(e))
      rows[[length(rows) + 1L]] <- list(method = m, ran = FALSE, estimate = e$estimate, lo = if (length(cd)) cd[[1]]$lo else e$ci_lower, hi = if (length(cd)) cd[[1]]$hi else e$ci_upper, p = e$p, n = e$n_ivs_used,
                                        what = if (ci_changed) list(sprintf("Confidence level: %s%% (default 95%%)", format(round(st$ci_level * 100, 1)))) else list())
    }
  }
  list(rows = rows, rerun = as.list(rerun))
}
