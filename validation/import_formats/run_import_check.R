# Import formats: imports every file of the corpus (make_corpus.R) the way a user does in "Import a GWAS file", pairs it with the
# reference dataset it was made from, and compares the pair with the reference pair. Run from the repository root:
#   Rscript validation/import_formats/run_import_check.R <corpus folder> <exposure dataset folder> <outcome dataset folder>
# The user is played honestly: the app's column suggestions and pre-filled answers are accepted; then only what a careful user knows
# from the source's documentation is answered (type of trait, cases and controls, role, which allele is the effect allele). Every
# answer that differs from what the app proposed is counted in user_changes, so the table shows where the app needs the user's help.
args <- commandArgs(TRUE)
corpus <- normalizePath(args[1]); ex_ref <- args[2]; out_ref <- args[3]
root <- file.path(tempdir(), "MR Studio"); dir.create(file.path(root, "Data library"), recursive = TRUE)
Sys.setenv(MRGUI_STUDIES_DIR = root)
args_dir <- normalizePath("prototype/engine"); options(mrgui.engine_dir = args_dir)
source(file.path(args_dir, "studio_core.R")); load_contracts(find_contracts_dir(args_dir))
for (f in c("import_core.R", "liftover_core.R", "opengwas_core.R", "ld_core.R")) source(file.path(args_dir, f))
file.copy(c(ex_ref, out_ref), file.path(root, "Data library"), recursive = TRUE)
EX <- basename(ex_ref); OUT <- basename(out_ref)
# The reference outcome was downloaded for the reference exposure's instruments; the imported copies of that exposure hold the same
# variants, so in this temporary copy the outcome is released from that one exposure (the app ties it, rightly, to avoid other exposures).
rp <- file.path(root, "Data library", OUT, "recipe.json"); rr <- jsonlite::fromJSON(rp, simplifyVector = FALSE); rr$source$for_exposure <- NULL; rr$lookup <- NULL
writeLines(jsonlite::toJSON(rr, auto_unbox = TRUE, pretty = TRUE, null = "null"), rp)
ivw <- function(d) { d <- d[d$mr_keep, ]; w <- 1 / d$se_y^2; b <- sum(w * d$beta_x * d$beta_y) / sum(w * d$beta_x^2); c(b = b, se = 1 / sqrt(sum(w * d$beta_x^2)), n = nrow(d)) }
ref <- ivw(pair_input(paste0("lib:", EX), paste0("lib:", OUT))$input$data)

# What the source documentation tells a careful user. The allele question is answered from the source's convention (for SAIGE the effect
# allele is Allele2): ea_ok "yes" means the app already proposed the right effect allele. (only these may differ from the app's proposal)
know <- list(
  x = list(ea_ok = "yes", role = "exposure", kind = "continuous", independent = "yes", scope = "selected", ancestry = "EUR"),
  y = list(ea_ok = "yes", role = "outcome", kind = "binary", ncase = "60801", ncontrol = "123504", independent = "no", scope = "selected", ancestry = "EUR"))
per_file <- list(   # facts the documentation of each source states
  "01" = list(build = "38"), "02" = list(build = "38"), "03" = list(build = "38"), "04" = list(build = "37"),
  "07y" = list(scale = "or"), "11y" = list(scale = "or"), "14" = list(scale = "z", eaf_fill = "panel"), "15y" = list(scale = "or"))
cc <- utils::read.csv(file.path(corpus, "corpus.csv"), stringsAsFactors = FALSE)
res <- lapply(seq_len(nrow(cc)), function(i) {
  f <- cc$file[i]; r <- substr(f, nchar(sub("\\..*$", "", f)), nchar(sub("\\..*$", "", f))); id <- substr(f, 1, 2); nm <- sprintf("imp %s %s", id, r)
  out <- data.frame(file = f, layout = cc$layout[i], role = cc$role[i], imported = FALSE, user_changes = "", mapping = "", problem = "", kept = NA, n_pair = NA, ivw = NA, ivw_se = NA, ref_ivw = ref[["b"]], diff = NA, stringsAsFactors = FALSE)
  tryCatch({
    h <- import_read_head(file.path(corpus, f), f); sg <- import_suggest(h$columns, h$rows); hn <- import_hints(h$columns, h$rows, sg)
    mp <- lapply(sg, `[[`, "column")
    out$mapping <- paste(sprintf("%s=%s", names(mp), unlist(mp)), collapse = " ")
    facts <- list(name = nm, role = "exposure", kind = "continuous", scale = "beta", ea_ok = "yes", frq = if (!is.null(mp$eaf)) "eaf" else "none", p_scale = "p", eaf_fill = "none", strand = "unknown", build = "unknown",
                  ancestry = "EUR", n = "", scope = "genomewide", sex = "both", independent = "no")
    facts <- utils::modifyList(facts, hn)                                    # what the app proposes
    want <- utils::modifyList(know[[r]], c(per_file[[id]] %||% list(), per_file[[paste0(id, r)]] %||% list()))
    if (!(facts$build %in% c("37", "38")) && is.null(want$build) && !is.null(mp$chr) && !is.null(mp$pos)) want$build <- "37"   # asked only when the file has positions
    ch <- names(want)[vapply(names(want), function(k) !identical(as.character(facts[[k]]), as.character(want[[k]])), TRUE)]
    out$user_changes <- paste(setdiff(ch, c("role", "kind", "ncase", "ncontrol", "independent", "scope", "ancestry")), collapse = " ")   # the questions every user answers
    facts <- utils::modifyList(facts, want)
    v <- import_validate(h$path, h$sep, h$columns, mp, facts)
    if (!v$ok) stop(paste(v$errors, collapse = " "))
    rec <- import_save(h$path, f, h$sep, h$columns, mp, facts, orig_path = file.path(corpus, f), prep_notes = unlist(h$notes))
    out$imported <- TRUE; out$kept <- rec$rows_kept
    p <- if (r == "x") pair_input(paste0("lib:", nm), paste0("lib:", OUT)) else pair_input(paste0("lib:", EX), paste0("lib:", nm))
    e <- ivw(p$input$data); out$n_pair <- e[["n"]]; out$ivw <- e[["b"]]; out$ivw_se <- e[["se"]]; out$diff <- e[["b"]] - ref[["b"]]
  }, error = function(e) out$problem <<- conditionMessage(e))
  out
})
res <- do.call(rbind, res)
utils::write.csv(res, file.path(dirname(sys.frame(1)$ofile %||% "validation/import_formats/x"), "results.csv"), row.names = FALSE)
cat(sprintf("reference IVW %.6f (SE %.6f, %d variants)\n", ref[["b"]], ref[["se"]], ref[["n"]]))
for (i in seq_len(nrow(res))) with(res[i, ], cat(sprintf("%-38s %-5s kept %-4s pair %-4s ivw %-9s diff %-10s changes [%s]\n   map: %s\n%s", file, imported, kept, n_pair, signif(ivw, 6), signif(diff, 3), user_changes, mapping, if (nzchar(problem)) paste0("   PROBLEM: ", problem, "\n") else "")))
