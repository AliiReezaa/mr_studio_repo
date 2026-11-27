# Publication package for MR Studio (Round 9). Base R only: no packages beyond the ones the engine already needs.

# Round 96: folders an analysis did not use (figures, method reports, outlier analysis) are left out of the package
.rm_empty_dirs <- function(dir) { d <- list.dirs(dir, recursive = TRUE, full.names = TRUE); d <- d[order(-nchar(d))]
  for (x in d) if (!identical(normalizePath(x), normalizePath(dir)) && !length(list.files(x, all.files = TRUE, no.. = TRUE))) unlink(x, recursive = TRUE); invisible(dir) }
.mrs_ver <- function() if (exists("MRS_VERSION")) MRS_VERSION else "(version unknown)"
.mrs_ver_build <- function() if (exists("mrs_version_text")) mrs_version_text() else .mrs_ver()
# Everything written here is a DRAFT for the author to check: numbers come from the saved result tables of the study,
# texts are assembled from the settings recorded in the generated scripts, and every limitation that applies is stated.
# Nothing here is an estimator.

# ---------- small helpers ----------
.xe <- function(s) { s <- as.character(s); s[is.na(s)] <- ""; s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE); s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("\"", "&quot;", s, fixed = TRUE) }
.fnum <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
.nth <- function(x, i) if (i <= length(x)) x[[i]] else NULL   # Round 50: lists that cover only the univariable pairs
.fp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-300, "< 1e-300", ifelse(p < 0.001, formatC(p, digits = 2, format = "e"), formatC(p, digits = 3, format = "f"))))   # Round 91: an underflowed P (0) is "< 1e-300", never "0.00e+00"
.fp_res <- function(p, res = NA_real_) { res <- rep_len(as.numeric(res %||% NA_real_), length(p)); out <- .fp(p); lo <- !is.na(p) & is.finite(res) & p < res; out[lo] <- paste("<", format(res[lo])); out }   # Round 103: a bootstrap P of 0 is "< 0.001" (1 / draws)
.pfix <- function(x) if (is.character(x)) gsub("= < ", "< ", x, fixed = TRUE) else if (is.list(x)) lapply(x, .pfix) else x
# Round 42: the radial detectors' level, as written in their scripts.
.radial_level_txt <- function(st) if (identical(st$alpha_correction %||% "bonferroni", "none")) sprintf("alpha = %s per variant, no correction", format(st$alpha %||% 0.05)) else sprintf("alpha = %s divided by the number of variants (Bonferroni)", format(st$alpha %||% 0.05))
.pct <- function(l) paste0(format(round(l * 100, 1), nsmall = 0), "%")

# ---------- minimal ZIP writer (stored, no compression) ----------
.crc_tab <- local({
  tab <- integer(256); poly <- -306674912L   # 0xEDB88320 as a signed 32-bit integer
  for (n in 0:255) { cc <- n; for (k in 1:8) cc <- if (bitwAnd(cc, 1L) != 0L) bitwXor(poly, bitwShiftR(cc, 1L)) else bitwShiftR(cc, 1L); tab[n + 1L] <- cc }
  tab })
crc32_raw <- function(x) {
  cc <- -1L; b <- as.integer(x)
  for (i in seq_along(b)) cc <- bitwXor(.crc_tab[bitwAnd(bitwXor(cc, b[i]), 255L) + 1L], bitwShiftR(cc, 8L))
  bitwXor(cc, -1L)
}
.u16 <- function(x) writeBin(as.integer(x), raw(), size = 2L, endian = "little")
.u32 <- function(x) writeBin(as.integer(x), raw(), size = 4L, endian = "little")
.dos <- function(t = Sys.time()) { lt <- as.POSIXlt(t); list(time = lt$hour * 2048L + lt$min * 32L + lt$sec %/% 2L, date = max(0L, (lt$year + 1900L - 1980L)) * 512L + (lt$mon + 1L) * 32L + lt$mday) }
#' entries: named list, name = path inside the archive, value = raw vector.
zip_write <- function(path, entries) {
  d <- .dos(); con <- file(path, "wb"); on.exit(close(con)); off <- 0L; cen <- list()
  for (nm in names(entries)) {
    data <- entries[[nm]]; nmr <- charToRaw(enc2utf8(nm)); crc <- crc32_raw(data); n <- length(data)
    hdr <- c(as.raw(c(0x50, 0x4b, 0x03, 0x04)), .u16(20L), .u16(2048L), .u16(0L), .u16(d$time), .u16(d$date), .u32(crc), .u32(n), .u32(n), .u16(length(nmr)), .u16(0L), nmr)
    writeBin(hdr, con); writeBin(data, con)
    cen[[length(cen) + 1L]] <- c(as.raw(c(0x50, 0x4b, 0x01, 0x02)), .u16(20L), .u16(20L), .u16(2048L), .u16(0L), .u16(d$time), .u16(d$date), .u32(crc), .u32(n), .u32(n),
                                 .u16(length(nmr)), .u16(0L), .u16(0L), .u16(0L), .u16(0L), .u32(0L), .u32(off), nmr)
    off <- off + length(hdr) + n
  }
  cd <- do.call(c, cen); writeBin(cd, con)
  writeBin(c(as.raw(c(0x50, 0x4b, 0x05, 0x06)), .u16(0L), .u16(0L), .u16(length(cen)), .u16(length(cen)), .u32(length(cd)), .u32(off), .u16(0L)), con)
  invisible(path)
}
.txt <- function(x) charToRaw(enc2utf8(paste(x, collapse = "\n")))
.file_raw <- function(f) readBin(f, "raw", file.info(f)$size)

# ---------- minimal DOCX writer ----------
# blocks: list of list(type = "h1"|"h2"|"p"|"title"|"note"|"table"|"image", ...)
docx_write <- function(path, blocks) {
  media <- list(); rels <- character(); body <- character()
  para <- function(txt, style = NULL, bold = FALSE, italic = FALSE) sprintf("<w:p>%s<w:r>%s<w:t xml:space=\"preserve\">%s</w:t></w:r></w:p>",
    if (!is.null(style)) sprintf("<w:pPr><w:pStyle w:val=\"%s\"/></w:pPr>", style) else "",
    if (bold || italic) sprintf("<w:rPr>%s%s</w:rPr>", if (bold) "<w:b/>" else "", if (italic) "<w:i/>" else "") else "", .xe(txt))
  for (b in blocks) {
    body <- c(body, switch(b$type,
      title = para(b$text, "Title"), h1 = para(b$text, "Heading1"), h2 = para(b$text, "Heading2"), p = para(b$text), note = para(b$text, italic = TRUE),
      table = {
        d <- b$data; cell <- function(v, hdr = FALSE) sprintf("<w:tc><w:tcPr><w:tcW w:w=\"0\" w:type=\"auto\"/></w:tcPr><w:p><w:r>%s<w:t xml:space=\"preserve\">%s</w:t></w:r></w:p></w:tc>", if (hdr) "<w:rPr><w:b/></w:rPr>" else "", .xe(v))
        wts <- if (!is.null(b$widths)) b$widths else vapply(names(d), function(n) max(nchar(n) + 3, min(24, max(nchar(as.character(d[[n]])), na.rm = TRUE))), 1); wd <- round(9298 * wts / sum(wts))   # widths: optional relative column widths
        cell <- function(v, w, hdr = FALSE) sprintf("<w:tc><w:tcPr><w:tcW w:w=\"%d\" w:type=\"dxa\"/></w:tcPr><w:p><w:pPr><w:spacing w:after=\"0\"/></w:pPr><w:r><w:rPr>%s<w:sz w:val=\"18\"/></w:rPr><w:t xml:space=\"preserve\">%s</w:t></w:r></w:p></w:tc>", w, if (hdr) "<w:b/>" else "", .xe(v))
        rows <- c(sprintf("<w:tr>%s</w:tr>", paste(vapply(seq_along(d), function(j) cell(names(d)[j], wd[j], TRUE), ""), collapse = "")),
                  vapply(seq_len(nrow(d)), function(i) sprintf("<w:tr>%s</w:tr>", paste(vapply(seq_along(d), function(j) cell(d[[j]][i], wd[j]), ""), collapse = "")), ""))
        paste0("<w:tbl><w:tblPr><w:tblStyle w:val=\"TableGrid\"/><w:tblW w:w=\"9298\" w:type=\"dxa\"/><w:tblLayout w:type=\"fixed\"/></w:tblPr><w:tblGrid>", paste(sprintf("<w:gridCol w:w=\"%d\"/>", wd), collapse = ""), "</w:tblGrid>", paste(rows, collapse = ""), "</w:tbl>", "<w:p/>")
      },
      image = {
        k <- length(media) + 1L; media[[sprintf("word/media/image%d.png", k)]] <- .file_raw(b$file); rid <- sprintf("rId%d", 100L + k)
        rels <- c(rels, sprintf("<Relationship Id=\"%s\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/image\" Target=\"media/image%d.png\"/>", rid, k))
        cx <- round(b$width_mm * 36000); cy <- round(b$height_mm * 36000)
        sprintf("<w:p><w:r><w:drawing><wp:inline distT=\"0\" distB=\"0\" distL=\"0\" distR=\"0\"><wp:extent cx=\"%d\" cy=\"%d\"/><wp:docPr id=\"%d\" name=\"Figure %d\" descr=\"%s\"/><a:graphic xmlns:a=\"http://schemas.openxmlformats.org/drawingml/2006/main\"><a:graphicData uri=\"http://schemas.openxmlformats.org/drawingml/2006/picture\"><pic:pic xmlns:pic=\"http://schemas.openxmlformats.org/drawingml/2006/picture\"><pic:nvPicPr><pic:cNvPr id=\"%d\" name=\"image%d.png\"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed=\"%s\"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x=\"0\" y=\"0\"/><a:ext cx=\"%d\" cy=\"%d\"/></a:xfrm><a:prstGeom prst=\"rect\"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>", cx, cy, k, k, .xe(b$alt %||% "Figure"), k, k, rid, cx, cy)
      }, ""))
  }
  doc <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\" xmlns:wp=\"http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing\"><w:body>", paste(body, collapse = ""), "<w:sectPr><w:pgSz w:w=\"11906\" w:h=\"16838\"/><w:pgMar w:top=\"1440\" w:right=\"1304\" w:bottom=\"1440\" w:left=\"1304\" w:header=\"708\" w:footer=\"708\" w:gutter=\"0\"/></w:sectPr></w:body></w:document>")
  sty <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><w:styles xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\"><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii=\"Calibri\" w:hAnsi=\"Calibri\" w:cs=\"Calibri\"/><w:sz w:val=\"22\"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after=\"120\" w:line=\"276\" w:lineRule=\"auto\"/></w:pPr></w:pPrDefault></w:docDefaults><w:style w:type=\"paragraph\" w:default=\"1\" w:styleId=\"Normal\"><w:name w:val=\"Normal\"/></w:style><w:style w:type=\"paragraph\" w:styleId=\"Title\"><w:name w:val=\"Title\"/><w:basedOn w:val=\"Normal\"/><w:pPr><w:spacing w:after=\"240\"/></w:pPr><w:rPr><w:b/><w:sz w:val=\"40\"/></w:rPr></w:style><w:style w:type=\"paragraph\" w:styleId=\"Heading1\"><w:name w:val=\"heading 1\"/><w:basedOn w:val=\"Normal\"/><w:pPr><w:keepNext/><w:spacing w:before=\"360\" w:after=\"120\"/><w:outlineLvl w:val=\"0\"/></w:pPr><w:rPr><w:b/><w:sz w:val=\"30\"/></w:rPr></w:style><w:style w:type=\"paragraph\" w:styleId=\"Heading2\"><w:name w:val=\"heading 2\"/><w:basedOn w:val=\"Normal\"/><w:pPr><w:keepNext/><w:spacing w:before=\"240\" w:after=\"80\"/><w:outlineLvl w:val=\"1\"/></w:pPr><w:rPr><w:b/><w:sz w:val=\"25\"/></w:rPr></w:style><w:style w:type=\"table\" w:styleId=\"TableGrid\"><w:name w:val=\"Table Grid\"/><w:tblPr><w:tblBorders><w:top w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/><w:left w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/><w:bottom w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/><w:right w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/><w:insideH w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/><w:insideV w:val=\"single\" w:sz=\"4\" w:color=\"999999\"/></w:tblBorders><w:tblCellMar><w:left w:w=\"80\" w:type=\"dxa\"/><w:right w:w=\"80\" w:type=\"dxa\"/></w:tblCellMar></w:tblPr></w:style></w:styles>"
  ct <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/><Default Extension=\"png\" ContentType=\"image/png\"/><Override PartName=\"/word/document.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml\"/><Override PartName=\"/word/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml\"/></Types>"
  rr <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"word/document.xml\"/></Relationships>"
  dr <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/>", paste(rels, collapse = ""), "</Relationships>")
  zip_write(path, c(list(`[Content_Types].xml` = .txt(ct), `_rels/.rels` = .txt(rr), `word/document.xml` = .txt(doc), `word/styles.xml` = .txt(sty), `word/_rels/document.xml.rels` = .txt(dr)), media))
}

# ---------- minimal XLSX writer ----------
.col_letter <- function(i) { s <- ""; while (i > 0) { r <- (i - 1L) %% 26L; s <- paste0(intToUtf8(65L + r), s); i <- (i - 1L) %/% 26L }; s }
xlsx_write <- function(path, sheets) {
  nm <- vapply(names(sheets), function(n) substr(gsub("[\\[\\]:*?/\\\\]", " ", n), 1, 31), ""); nm <- make.unique(nm, sep = " ")
  sheet_xml <- function(d) {
    cells <- function(vals, r) paste0("<row r=\"", r, "\">", paste(vapply(seq_along(vals), function(j) {
      v <- vals[[j]]; ref <- paste0(.col_letter(j), r)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) sprintf("<c r=\"%s\"><v>%s</v></c>", ref, format(v, digits = 15, scientific = FALSE, trim = TRUE))
      else if (is.logical(v) && length(v) == 1L && !is.na(v)) sprintf("<c r=\"%s\" t=\"b\"><v>%d</v></c>", ref, as.integer(v))
      else if (length(v) == 1L && !is.na(v)) sprintf("<c r=\"%s\" t=\"inlineStr\"><is><t xml:space=\"preserve\">%s</t></is></c>", ref, .xe(v)) else sprintf("<c r=\"%s\"/>", ref)
    }, ""), collapse = ""), "</row>")
    rows <- c(cells(as.list(names(d)), 1L), vapply(seq_len(nrow(d)), function(i) cells(lapply(d, function(col) col[[i]]), i + 1L), ""))
    cw <- vapply(seq_along(d), function(j) min(60, max(10, nchar(names(d)[j]) + 3, suppressWarnings(max(nchar(as.character(d[[j]])), na.rm = TRUE)) + 2)), 1)
    colsxml <- paste0("<cols>", paste(sprintf("<col min=\"%d\" max=\"%d\" width=\"%.1f\" customWidth=\"1\"/>", seq_along(cw), seq_along(cw), cw), collapse = ""), "</cols>")
    paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><sheetViews><sheetView workbookViewId=\"0\"><pane ySplit=\"1\" topLeftCell=\"A2\" activePane=\"bottomLeft\" state=\"frozen\"/></sheetView></sheetViews>", colsxml, "<sheetData>", paste(rows, collapse = ""), "</sheetData></worksheet>")
  }
  n <- length(sheets)
  ct <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"xml\" ContentType=\"application/xml\"/><Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>", paste(sprintf("<Override PartName=\"/xl/worksheets/sheet%d.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>", seq_len(n)), collapse = ""), "</Types>")
  rr <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/></Relationships>"
  wb <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><sheets>", paste(sprintf("<sheet name=\"%s\" sheetId=\"%d\" r:id=\"rId%d\"/>", .xe(nm), seq_len(n), seq_len(n)), collapse = ""), "</sheets></workbook>")
  wr <- paste0("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?><Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">", paste(sprintf("<Relationship Id=\"rId%d\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet%d.xml\"/>", seq_len(n), seq_len(n)), collapse = ""), "</Relationships>")
  ent <- c(list(`[Content_Types].xml` = .txt(ct), `_rels/.rels` = .txt(rr), `xl/workbook.xml` = .txt(wb), `xl/_rels/workbook.xml.rels` = .txt(wr)),
           setNames(lapply(sheets, function(d) .txt(sheet_xml(d))), sprintf("xl/worksheets/sheet%d.xml", seq_len(n))))
  zip_write(path, ent)
}

# ---------- references (verified against PubMed records) ----------
REFS <- list(
  hemani2018 = list(pmid = "29846171", key = "Hemani2018", authors = "Hemani G, Zheng J, Elsworth B, et al.", title = "The MR-Base platform supports systematic causal inference across the human phenome", journal = "eLife", year = "2018", vol = "7", pages = "e34408", doi = "10.7554/eLife.34408"),
  burgess2013 = list(pmid = "24114802", key = "Burgess2013", authors = "Burgess S, Butterworth A, Thompson SG.", title = "Mendelian randomization analysis with multiple genetic variants using summarized data", journal = "Genet Epidemiol", year = "2013", vol = "37", pages = "658-665", doi = "10.1002/gepi.21758"),
  bowden2015 = list(pmid = "26050253", key = "Bowden2015", authors = "Bowden J, Davey Smith G, Burgess S.", title = "Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression", journal = "Int J Epidemiol", year = "2015", vol = "44", pages = "512-525", doi = "10.1093/ije/dyv080"),
  bowden2016 = list(pmid = "27061298", key = "Bowden2016", authors = "Bowden J, Davey Smith G, Haycock PC, Burgess S.", title = "Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator", journal = "Genet Epidemiol", year = "2016", vol = "40", pages = "304-314", doi = "10.1002/gepi.21965"),
  hartwig2017 = list(pmid = "29040600", key = "Hartwig2017", authors = "Hartwig FP, Davey Smith G, Bowden J.", title = "Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption", journal = "Int J Epidemiol", year = "2017", vol = "46", pages = "1985-1998", doi = "10.1093/ije/dyx102"),
  burgess2023 = list(pmid = "32760811", key = "Burgess2023", authors = "Burgess S, Davey Smith G, Davies NM, et al.", title = "Guidelines for performing Mendelian randomization investigations: update for summer 2023", journal = "Wellcome Open Res", year = "2023", vol = "4", pages = "186", doi = "10.12688/wellcomeopenres.15555.3"),
  skrivankova2021 = list(pmid = "34698778", key = "Skrivankova2021", authors = "Skrivankova VW, Richmond RC, Woolf BAR, et al.", title = "Strengthening the reporting of observational studies in epidemiology using Mendelian randomisation (STROBE-MR): explanation and elaboration", journal = "JAMA", year = "2021", vol = "326", pages = "1614-1621", doi = "10.1001/jama.2021.18236"))
METHOD_REFS <- list(ivw_steiger = "hemani2017", bwmr = "zhao2020", grip = "dudbridge2025", mrmix = "qi2019", egger_boot = "bowden2015", pwm = "burgess2017sens", smode_nome = "hartwig2017", wmode_nome = "hartwig2017", rucker = "bowden2018het", sign = "hemani2018", maxlik = "burgess2013", mbe = "hartwig2017", cml = "xue2021", robust_ivw = "rees2019robust", pivw = "xu2023", cook = "cook1977", studres = "cook1977", mrclust = "foley2021", gsmr = "zhu2018", grapple = "wang2021grapple", mvgrapple = "wang2021grapple", ppmr = "yuan2020", egger_rob = "rees2019robust", hetpen = "burgess2018hetpen", clr = "wang2022clr", mvivwme = "patel2023", mvgmm = "patel2023", mvcml = "lin2023", mvivw = c("sanderson2019", "sanderson2021"), mvegger = "rees2017", mvmedian = "grant2021", mvlasso = "grant2021", mvqhet = "sanderson2021", ivw_mre = "burgess2013", divw = "ye2021", raps = "zhao2020", conmix = "burgess2020", lasso = "rees2019", egger_simex = "bowden2016i2", egger = "bowden2015", wm = "bowden2016", smed = "bowden2016", mode = "hartwig2017", smode = "hartwig2017", presso = "verbanck2018", radial = "bowden2018", radegger = "bowden2018", steiger = "hemani2017")
REFS$bowden2016i2 <- list(pmid = "27616674", key = "Bowden2016I2", authors = "Bowden J, Del Greco M F, Minelli C, Davey Smith G, Sheehan NA, Thompson JR.", title = "Assessing the suitability of summary data for two-sample Mendelian randomization analyses using MR-Egger regression: the role of the I2 statistic", journal = "Int J Epidemiol", year = "2016", vol = "45", pages = "1961-1974", doi = "10.1093/ije/dyw220")
REFS$burgess2020 <- list(pmid = "31953392", key = "Burgess2020", authors = "Burgess S, Foley CN, Allara E, Staley JR, Howson JMM.", title = "A robust and efficient method for Mendelian randomization with hundreds of genetic variants", journal = "Nat Commun", year = "2020", vol = "11", pages = "376", doi = "10.1038/s41467-019-14156-4")
# Round 50: multivariable MR (PMIDs checked in PubMed on 2026-11-18).
REFS$sanderson2019 <- list(pmid = "30535378", key = "Sanderson2019", authors = "Sanderson E, Davey Smith G, Windmeijer F, Bowden J.", title = "An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings", journal = "Int J Epidemiol", year = "2019", vol = "48", pages = "713-727", doi = "10.1093/ije/dyy262")
REFS$sanderson2021 <- list(pmid = "34338327", key = "Sanderson2021", authors = "Sanderson E, Spiller W, Bowden J.", title = "Testing and correcting for weak and pleiotropic instruments in two-sample multivariable Mendelian randomization", journal = "Stat Med", year = "2021", vol = "40", pages = "5434-5452", doi = "10.1002/sim.9133")
REFS$carter2021 <- list(pmid = "33961203", key = "Carter2021", authors = "Carter AR, Sanderson E, Hammerton G, Richmond RC, Davey Smith G, Heron J, Taylor AE, Davies NM, Howe LD.", title = "Mendelian randomisation for mediation analysis: current methods and challenges for implementation", journal = "Eur J Epidemiol", year = "2021", vol = "36", pages = "465-478", doi = "10.1007/s10654-021-00757-1")
REFS$carter2019 <- list(pmid = "31122926", key = "Carter2019", authors = "Carter AR, Gill D, Davies NM, et al.", title = "Understanding the consequences of education inequality on cardiovascular disease: mendelian randomisation study", journal = "BMJ", year = "2019", vol = "365", pages = "l1855", doi = "10.1136/bmj.l1855")
REFS$bulik2015a <- list(pmid = "25642630", key = "BulikSullivan2015a", authors = "Bulik-Sullivan BK, Loh PR, Finucane HK, et al.", title = "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies", journal = "Nat Genet", year = "2015", vol = "47", pages = "291-295", doi = "10.1038/ng.3211")
REFS$bulik2015b <- list(pmid = "26414676", key = "BulikSullivan2015b", authors = "Bulik-Sullivan B, Finucane HK, Anttila V, et al.", title = "An atlas of genetic correlations across human diseases and traits", journal = "Nat Genet", year = "2015", vol = "47", pages = "1236-1241", doi = "10.1038/ng.3406")
REFS$grotzinger2019 <- list(pmid = "30962613", key = "Grotzinger2019", authors = "Grotzinger AD, Rhemtulla M, de Vlaming R, et al.", title = "Genomic structural equation modelling provides insights into the multivariate genetic architecture of complex traits", journal = "Nat Hum Behav", year = "2019", vol = "3", pages = "513-525", doi = "10.1038/s41562-019-0566-x")
REFS$darrous2021 <- list(pmid = "34907193", key = "Darrous2021", authors = "Darrous L, Mounier N, Kutalik Z.", title = "Simultaneous estimation of bi-directional causal effects and heritable confounding from GWAS summary statistics", journal = "Nat Commun", year = "2021", vol = "12", pages = "7274", doi = "10.1038/s41467-021-26970-w")
REFS$mounier2023 <- list(pmid = "37036286", key = "Mounier2023", authors = "Mounier N, Kutalik Z.", title = "Bias correction for inverse variance weighting Mendelian randomization", journal = "Genet Epidemiol", year = "2023", vol = "47", pages = "314-331", doi = "10.1002/gepi.22522")
REFS$morrison2020 <- list(pmid = "32451458", key = "Morrison2020", authors = "Morrison J, Knoblauch N, Marcus JH, Stephens M, He X", title = "Mendelian randomization accounting for correlated and uncorrelated pleiotropic effects using genome-wide summary statistics", journal = "Nat Genet", year = "2020", vol = "52", pages = "740-747", doi = "10.1038/s41588-020-0631-4")
REFS$hu2022 <- list(pmid = "35787050", key = "Hu2022", authors = "Hu X, Zhao J, Lin Z, et al.", title = "Mendelian randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics", journal = "Proc Natl Acad Sci U S A", year = "2022", vol = "119", pages = "e2106858119", doi = "10.1073/pnas.2106858119")
REFS$xue2022 <- list(pmid = "35576237", key = "Xue2022", authors = "Xue H, Pan W.", title = "Robust inference of bi-directional causal relationships in presence of correlated pleiotropy with GWAS summary data", journal = "PLoS Genet", year = "2022", vol = "18", pages = "e1010205", doi = "10.1371/journal.pgen.1010205")
REFS$zhao2020 <- list(pmid = "31593215", key = "Zhao2020", authors = "Zhao J, Ming J, Hu X, Chen G, Liu J, Yang C.", title = "Bayesian weighted Mendelian randomization for causal inference based on summary statistics", journal = "Bioinformatics", year = "2020", vol = "36", pages = "1501-1508", doi = "10.1093/bioinformatics/btz749")
REFS$dudbridge2025 <- list(pmid = "41468517", key = "Dudbridge2025", authors = "Dudbridge F, Voller B, Woodward RM, Saxby KL, Frayling TM, Pilling LC, Bowden J.", title = "Getting to GRIPS with MR-Egger: modelling directional pleiotropy independently of allele coding", journal = "PLoS Genet", year = "2025", vol = "21", pages = "e1011967", doi = "10.1371/journal.pgen.1011967")
REFS$qi2019 <- list(pmid = "31028273", key = "Qi2019", authors = "Qi G, Chatterjee N.", title = "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects", journal = "Nat Commun", year = "2019", vol = "10", pages = "1941", doi = "10.1038/s41467-019-09432-2")
REFS$burgess2017sens <- list(pmid = "27749700", key = "Burgess2017", authors = "Burgess S, Bowden J, Fall T, Ingelsson E, Thompson SG.", title = "Sensitivity analyses for robust causal inference from Mendelian randomization analyses with multiple genetic variants", journal = "Epidemiology", year = "2017", vol = "28", pages = "30-42", doi = "10.1097/EDE.0000000000000559")
REFS$bowden2018het <- list(pmid = "30188969", key = "Bowden2018", authors = "Bowden J, Hemani G, Davey Smith G.", title = "Invited commentary: detecting individual and global horizontal pleiotropy in Mendelian randomization, a job for the humble heterogeneity statistic?", journal = "Am J Epidemiol", year = "2018", vol = "187", pages = "2681-2685", doi = "10.1093/aje/kwy185")
REFS$xue2021 <- list(pmid = "34214446", key = "Xue2021", authors = "Xue H, Shen X, Pan W.", title = "Constrained maximum likelihood-based Mendelian randomization robust to both correlated and uncorrelated pleiotropic effects", journal = "Am J Hum Genet", year = "2021", vol = "108", pages = "1251-1269", doi = "10.1016/j.ajhg.2021.05.014")
REFS$rees2019robust <- list(pmid = "31545794", key = "Rees2019b", authors = "Rees JMB, Wood AM, Dudbridge F, Burgess S.", title = "Robust methods in Mendelian randomization via penalization of heterogeneous causal estimates", journal = "PLoS One", year = "2019", vol = "14", pages = "e0222362", doi = "10.1371/journal.pone.0222362")
REFS$xu2023 <- list(pmid = "35942938", key = "Xu2023", authors = "Xu S, Wang P, Fung WK, Liu Z.", title = "A novel penalized inverse-variance weighted estimator for Mendelian randomization with applications to COVID-19 outcomes", journal = "Biometrics", year = "2023", vol = "79", pages = "2184-2195", doi = "10.1111/biom.13732")
REFS$burgess2018hetpen <- list(pmid = "29846613", key = "Burgess2018", authors = "Burgess S, Zuber V, Gkatzionis A, Foley CN.", title = "Modal-based estimation via heterogeneity-penalized weighting: model averaging for consistent and efficient estimation in Mendelian randomization when a plurality of candidate instruments are valid", journal = "Int J Epidemiol", year = "2018", vol = "47", pages = "1242-1254", doi = "10.1093/ije/dyy080")   # Round 109
REFS$wang2022clr <- list(pmid = "34213007", key = "Wang2022", authors = "Wang S, Kang H.", title = "Weak-instrument robust tests in two-sample summary-data Mendelian randomization", journal = "Biometrics", year = "2022", vol = "78", pages = "1699-1713", doi = "10.1111/biom.13524")   # Round 109
REFS$cook1977 <- list(pmid = "", key = "Cook1977", authors = "Cook RD.", title = "Detection of influential observation in linear regression", journal = "Technometrics", year = "1977", vol = "19", pages = "15-18", doi = "10.1080/00401706.1977.10489493")   # Round 109 (not in PubMed)
REFS$foley2021 <- list(pmid = "32915962", key = "Foley2021", authors = "Foley CN, Mason AM, Kirk PDW, Burgess S.", title = "MR-Clust: clustering of genetic variants in Mendelian randomization with similar causal estimates", journal = "Bioinformatics", year = "2021", vol = "37", pages = "531-541", doi = "10.1093/bioinformatics/btaa778")   # Round 109
REFS$zhu2018 <- list(pmid = "29335400", key = "Zhu2018", authors = "Zhu Z, Zheng Z, Zhang F, et al.", title = "Causal associations between risk factors and common diseases inferred from GWAS summary data", journal = "Nat Commun", year = "2018", vol = "9", pages = "224", doi = "10.1038/s41467-017-02317-2")   # Round 109
REFS$wang2021grapple <- list(pmid = "34157017", key = "Wang2021", authors = "Wang J, Zhao Q, Bowden J, et al.", title = "Causal inference for heritable phenotypic risk factors using heterogeneous genetic instruments", journal = "PLoS Genet", year = "2021", vol = "17", pages = "e1009575", doi = "10.1371/journal.pgen.1009575")   # Round 109
REFS$yuan2020 <- list(pmid = "32737316", key = "Yuan2020", authors = "Yuan Z, Zhu H, Zeng P, et al.", title = "Testing and controlling for horizontal pleiotropy with probabilistic Mendelian randomization in transcriptome-wide association studies", journal = "Nat Commun", year = "2020", vol = "11", pages = "3861", doi = "10.1038/s41467-020-17668-6")   # Round 109
REFS$patel2023 <- list(pmid = "37915953", key = "Patel2023", authors = "Patel A, Ye T, Xue H, et al.", title = "MendelianRandomization v0.9.0: updates to an R package for performing Mendelian randomization analyses using summarized data", journal = "Wellcome Open Res", year = "2023", vol = "8", pages = "449", doi = "10.12688/wellcomeopenres.19995.2")   # Round 109
REFS$lin2023 <- list(pmid = "36948188", key = "Lin2023", authors = "Lin Z, Xue H, Pan W.", title = "Robust multivariable Mendelian randomization based on constrained maximum likelihood", journal = "Am J Hum Genet", year = "2023", vol = "110", pages = "592-605", doi = "10.1016/j.ajhg.2023.02.014")
REFS$rees2017 <- list(pmid = "28960498", key = "Rees2017", authors = "Rees JMB, Wood AM, Burgess S.", title = "Extending the MR-Egger method for multivariable Mendelian randomization to correct for both measured and unmeasured pleiotropy", journal = "Stat Med", year = "2017", vol = "36", pages = "4705-4718", doi = "10.1002/sim.7492")
REFS$grant2021 <- list(pmid = "34342032", key = "Grant2021", authors = "Grant AJ, Burgess S.", title = "Pleiotropy robust methods for multivariable Mendelian randomization", journal = "Stat Med", year = "2021", vol = "40", pages = "5813-5830", doi = "10.1002/sim.9156")
REFS$rees2019 <- list(pmid = "31545794", key = "Rees2019", authors = "Rees JMB, Wood AM, Dudbridge F, Burgess S.", title = "Robust methods in Mendelian randomization via penalization of heterogeneous causal estimates", journal = "PLoS One", year = "2019", vol = "14", pages = "e0222362", doi = "10.1371/journal.pone.0222362")
REFS$zhao2020 <- list(pmid = "", key = "Zhao2020", authors = "Zhao Q, Wang J, Hemani G, Bowden J, Small DS.", title = "Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score", journal = "Ann Stat", year = "2020", vol = "48", pages = "1742-1769", doi = "10.1214/19-AOS1866")
REFS$ye2021 <- list(pmid = "", key = "Ye2021", authors = "Ye T, Shao J, Kang H.", title = "Debiased inverse-variance weighted estimator in two-sample summary-data Mendelian randomization", journal = "Ann Stat", year = "2021", vol = "49", pages = "2079-2100", doi = "10.1214/20-AOS2027")
REFS$hemani2017 <- list(pmid = "29149188", key = "Hemani2017", authors = "Hemani G, Tilling K, Davey Smith G.", title = "Orienting the causal relationship between imprecisely measured traits using GWAS summary data", journal = "PLoS Genet", year = "2017", vol = "13", pages = "e1007081", doi = "10.1371/journal.pgen.1007081")
REFS$skrivankova2021ee <- list(pmid = "34702754", key = "Skrivankova2021EE", authors = "Skrivankova VW, Richmond RC, Woolf BAR, et al.", title = "Strengthening the reporting of observational studies in epidemiology using mendelian randomisation (STROBE-MR): explanation and elaboration", journal = "BMJ", year = "2021", vol = "375", pages = "n2233", doi = "10.1136/bmj.n2233")
REFS$chang2015 <- list(pmid = "25722852", key = "Chang2015", authors = "Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ.", title = "Second-generation PLINK: rising to the challenge of larger and richer datasets", journal = "Gigascience", year = "2015", vol = "4", pages = "7", doi = "10.1186/s13742-015-0047-8")
REFS$kgp2015 <- list(pmid = "26432245", key = "1000Genomes2015", authors = "1000 Genomes Project Consortium, Auton A, Brooks LD, et al.", title = "A global reference for human genetic variation", journal = "Nature", year = "2015", vol = "526", pages = "68-74", doi = "10.1038/nature15393")
REFS$burgess2018 <- list(pmid = "30039250", key = "BurgessLabrecque2018", authors = "Burgess S, Labrecque JA.", title = "Mendelian randomization with a binary exposure variable: interpretation and presentation of causal estimates", journal = "Eur J Epidemiol", year = "2018", vol = "33", pages = "947-952", doi = "10.1007/s10654-018-0424-6")
REFS$lee2012 <- list(pmid = "22714935", key = "Lee2012", authors = "Lee SH, Goddard ME, Wray NR, Visscher PM.", title = "A better coefficient of determination for genetic profile analysis", journal = "Genet Epidemiol", year = "2012", vol = "36", pages = "214-224", doi = "10.1002/gepi.21614")
REFS$verbanck2018 <- list(pmid = "29686387", key = "Verbanck2018", authors = "Verbanck M, Chen CY, Neale B, Do R.", title = "Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases", journal = "Nat Genet", year = "2018", vol = "50", pages = "693-698", doi = "10.1038/s41588-018-0099-7")
REFS$bowden2018 <- list(pmid = "29961852", key = "Bowden2018", authors = "Bowden J, Spiller W, Del Greco M F, et al.", title = "Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression", journal = "Int J Epidemiol", year = "2018", vol = "47", pages = "1264-1278", doi = "10.1093/ije/dyy101")
.has_pmid <- function(r) !is.null(r$pmid) && nzchar(r$pmid)   # Round 37: journals outside PubMed (for example Ann Stat) have no PMID
ref_text <- function(r) paste0(sprintf("%s %s. %s. %s;%s:%s. doi:%s.", r$authors, r$title, r$journal, r$year, r$vol, r$pages, r$doi), if (.has_pmid(r)) sprintf(" PMID %s.", r$pmid) else "")
refs_bibtex <- function(ids) paste(vapply(ids, function(i) { r <- REFS[[i]]; sprintf("@article{%s,\n  author = {%s},\n  title = {%s},\n  journal = {%s},\n  year = {%s},\n  volume = {%s},\n  pages = {%s},\n  doi = {%s}%s\n}\n", r$key, gsub(", et al\\.", " and others", gsub("\\.$", "", gsub(", ", " and ", r$authors))), r$title, r$journal, r$year, r$vol, gsub("-", "--", r$pages), r$doi, if (.has_pmid(r)) sprintf(",\n  note = {PMID %s}", r$pmid) else "") }, ""), collapse = "\n")
refs_ris <- function(ids) paste(vapply(ids, function(i) { r <- REFS[[i]]; paste(c("TY  - JOUR", sprintf("AU  - %s", trimws(strsplit(gsub("\\.$", "", gsub(", et al", "", r$authors)), ", ")[[1]])), sprintf("TI  - %s", r$title), sprintf("JO  - %s", r$journal), sprintf("PY  - %s", r$year), sprintf("VL  - %s", r$vol), sprintf("SP  - %s", r$pages), sprintf("DO  - %s", r$doi), if (.has_pmid(r)) sprintf("AN  - %s", r$pmid), "ER  - ", ""), collapse = "\n") }, ""), collapse = "\n")

# ---------- figures: one drawing function, written as PDF, SVG and PNG so they cannot disagree ----------
.draw_forest <- function(est, title = NULL) {
  op <- graphics::par(mar = c(3.4, 8.6, 1.4, 1.2), mgp = c(1.9, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  mi <- graphics::par("mai"); lw <- max(graphics::strwidth(est$label, units = "inches")) + 0.15   # Round 98: room for the longest label (up to half the width)
  if (lw > mi[2]) { mi[2] <- min(lw, graphics::par("din")[1] * 0.5); graphics::par(mai = mi) }
  k <- nrow(est); rt <- isTRUE(est$ratio[1]); xl <- range(c(est$ci_lower, est$ci_upper, if (rt) 1 else 0), na.rm = TRUE); pad <- diff(xl) * 0.06   # Round 109: an OR axis is logarithmic, with no effect at 1
  if (rt) { xl <- xl * c(0.97, 1.03); pad <- 0 }
  graphics::plot(NA, xlim = xl + c(-pad, pad), ylim = c(k + 0.6, 0.4), yaxt = "n", xlab = est$unit[1], ylab = "", bty = "n", log = if (rt) "x" else "")
  graphics::abline(v = if (isTRUE(est$ratio[1])) 1 else 0, lty = 2, col = "grey55"); graphics::axis(2, at = seq_len(k), labels = est$label, tick = FALSE)
  graphics::segments(est$ci_lower, seq_len(k), est$ci_upper, seq_len(k), lwd = 1.6, col = est$col); graphics::points(est$estimate, seq_len(k), pch = 16, cex = 1.1, col = est$col)
  if (!is.null(title)) graphics::mtext(title, 3, 0.2, adj = 0, font = 2, cex = 0.8)
}
.draw_scatter <- function(v, lines, xlab, ylab) {
  op <- graphics::par(mar = c(4.4, 4.9, 0.8, 0.8), mgp = c(3.2, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  s <- ifelse(v$bx < 0, -1, 1); x <- v$bx * s; y <- v$by * s
  graphics::plot(x, y, pch = 16, col = "grey25", xlim = c(0, max(x) * 1.05), ylim = range(c(y, 0)) * 1.1, xlab = xlab, ylab = ylab, bty = "l")
  graphics::abline(h = 0, col = "grey75"); graphics::segments(x, y - 1.96 * v$sy, x, y + 1.96 * v$sy, col = "grey75", lwd = 0.6)
  for (l in lines) graphics::abline(a = l$a, b = l$b, col = l$col, lty = l$lty, lwd = 2)
  graphics::legend("topleft", legend = vapply(lines, `[[`, "", "name"), col = vapply(lines, `[[`, "", "col"), lty = vapply(lines, `[[`, 1, "lty"), lwd = 2, bty = "n", cex = 0.85)
}
# Round 39: leave-one-out plot (IVW without each variant; above 60 variants, the 30 whose removal changes it most) and funnel plot.
LOO_PLOT_MAX <- 60L; LOO_PLOT_SHOWN <- 30L
.loo_rows <- function(lf, full_b) { d <- lf; if (nrow(d) > LOO_PLOT_MAX) d <- d[order(-abs(d$loo_b - full_b))[seq_len(LOO_PLOT_SHOWN)], , drop = FALSE]; d[order(d$loo_b), , drop = FALSE] }
.draw_loo <- function(lf, full, z, unit) {
  d <- .loo_rows(lf, full$b); k <- nrow(d) + 1L
  b <- c(d$loo_b, full$b); lo <- c(d$loo_b - z * d$loo_se, full$lo); hi <- c(d$loo_b + z * d$loo_se, full$hi)
  op <- graphics::par(mar = c(3.4, 8.6, 0.8, 1.2), mgp = c(1.9, 0.6, 0), cex = 0.72, las = 1); on.exit(graphics::par(op))
  xl <- range(c(lo, hi, 0), na.rm = TRUE); pad <- diff(xl) * 0.06; cols <- c(rep("grey25", k - 1L), "#2b3fc7")
  graphics::plot(NA, xlim = xl + c(-pad, pad), ylim = c(k + 0.6, 0.4), yaxt = "n", xlab = unit, ylab = "", bty = "n")
  graphics::abline(v = 0, lty = 2, col = "grey55"); graphics::abline(v = full$b, lty = 3, col = "#2b3fc7")
  graphics::axis(2, at = seq_len(k), labels = c(paste("without", d$id), "All variants"), tick = FALSE, cex.axis = 0.85)
  graphics::segments(lo, seq_len(k), hi, seq_len(k), lwd = 1.2, col = cols); graphics::points(b, seq_len(k), pch = c(rep(16, k - 1L), 18), cex = c(rep(0.8, k - 1L), 1.4), col = cols)
}
.draw_funnel <- function(lf, lines) {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.9, 0.8, 0.8), mgp = c(3.2, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  graphics::plot(lf$wr, 1 / lf$wr_se, pch = 16, col = "grey25", xlab = "Single-variant estimate (Wald ratio)", ylab = "Precision (1 / standard error)", bty = "l")
  graphics::abline(v = 0, col = "grey75"); for (l in lines) graphics::abline(v = l$b, col = l$col, lty = l$lty, lwd = 2)
  if (length(lines)) graphics::legend("topright", legend = vapply(lines, `[[`, "", "name"), col = vapply(lines, `[[`, "", "col"), lty = vapply(lines, `[[`, 1, "lty"), lwd = 2, bty = "n", cex = 0.85)
  # Round 59: density of the single-variant estimates weighted by their precision (as TwoSampleMR::mr_density_plot)
  ok <- is.finite(lf$wr) & is.finite(lf$wr_se) & lf$wr_se > 0
  if (sum(ok) >= 3) { w <- 1 / lf$wr_se[ok]^2; dn <- stats::density(lf$wr[ok], weights = w / sum(w))
    graphics::plot(dn, main = "", xlab = "Single-variant estimate (Wald ratio)", ylab = "Density (weighted by precision)", bty = "l"); graphics::rug(lf$wr[ok], col = "grey45")
    graphics::abline(v = 0, col = "grey75"); for (l in lines) graphics::abline(v = l$b, col = l$col, lty = l$lty, lwd = 2) } else graphics::plot.new()
}
#' Leave-one-out and Wald ratio results of one pair, read from the loo diagnostics (NULL if the analysis did not run).
loo_facts <- function(cx) {
  dg <- Filter(function(d) identical(d$method, "loo"), cx$pay$diagnostics); if (!length(dg)) return(NULL)
  g <- function(n) { x <- Filter(function(d) identical(d$name, n), dg)
    data.frame(id = vapply(x, function(d) as.character(d$note), ""), v = vapply(x, function(d) as.numeric(d$value), 0), p = vapply(x, function(d) suppressWarnings(as.numeric(d$p %||% NA)), 0), stringsAsFactors = FALSE) }
  b <- g("loo_estimate"); if (!nrow(b)) return(NULL); se <- g("loo_se"); w <- g("wald_ratio"); ws <- g("wald_ratio_se")
  data.frame(id = b$id, loo_b = b$v, loo_se = se$v[match(b$id, se$id)], loo_p = b$p, wr = w$v[match(b$id, w$id)], wr_se = ws$v[match(b$id, ws$id)], wr_p = w$p[match(b$id, w$id)], stringsAsFactors = FALSE)
}
#' Round 43 (FR-105): the dominant-variant and leave-subset-out rows of the leave-one-out task, or NULL.
subset_facts <- function(pay) {
  dg <- Filter(function(d) identical(d$method, "loo"), pay$diagnostics); g <- function(n) { x <- Filter(function(d) identical(d$name, n), dg); if (length(x)) x[[1]] else NULL }
  dm <- g("dominant_share"); if (is.null(dm)) return(NULL)
  lb <- Filter(function(d) identical(d$name, "loo_estimate") && identical(as.character(d$note), as.character(dm$note)), dg)
  r <- g("lso_repeats"); v <- function(n) { x <- g(n); if (is.null(x)) NA_real_ else as.numeric(x$value) }
  list(dom_id = as.character(dm$note), dom_share = as.numeric(dm$value), dom_without = if (length(lb)) as.numeric(lb[[1]]$value) else NA_real_,
       lso = !is.null(r), lso_n = if (is.null(r)) NA_integer_ else as.integer(r$value), lso_note = if (is.null(r)) "" else as.character(r$note),
       lso_lo = v("lso_q025"), lso_med = v("lso_median"), lso_hi = v("lso_q975"), lso_min = v("lso_min"), lso_max = v("lso_max"))
}

.MCOL <- c(ivw_steiger = "#6b7f3a", bwmr = "#5a4a8a", grip = "#9a4f2a", mrmix = "#3a7a6b", egger_boot = "#c2703a", pwm = "#5f86c4", smode_nome = "#d98fb0", wmode_nome = "#b85c8a", rucker = "#6f7d3a", ivw_mre = "#2b3fc7", divw = "#5b6b8c", maxlik = "#3f6fb5", mbe = "#a0527a", cml = "#2f8a7a", robust_ivw = "#6b5b2e", pivw = "#4a5aa0", gsmr = "#2f7d6d", grapple = "#6a5a2f", ppmr = "#2f5a6a", egger_rob = "#b8743a", hetpen = "#7a5a9a", clr = "#4f7f8f", mvivwme = "#3b6fa0", mvgmm = "#8a6a3a", raps = "#7a5c2e", conmix = "#4f7a3a", lasso = "#7a3a6b", egger_simex = "#d4874f", egger = "#b5541c", wm = "#2f7fd6", smed = "#8c4fd6", mode = "#d6407a", smode = "#2a9d8f", presso = "#b8860b", radial = "#3d9a38", radegger = "#a04a3d")
.fig_cairo <- function() isTRUE(unname(capabilities("cairo")))
#' Clears R's cached height of "M" (used to centre text), keyed by the device address: a new device at the address of one just closed
#' could otherwise centre its text with the old device's units (labels about 11 px off at 300 dpi). A throwaway device measures "M" at an unusual size.
.fig_metric_reset <- function() tryCatch({ f <- tempfile(fileext = ".png"); grDevices::png(f, width = 60, height = 60, type = if (isTRUE(unname(capabilities("cairo")))) "cairo" else getOption("bitmapType")); on.exit({ grDevices::dev.off(); unlink(f) }); graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new(); graphics::strheight("M", cex = 0.0137); invisible(NULL) }, error = function(e) invisible(NULL))
.fig_write <- function(base, w_mm, h_mm, draw) {
  if (exists(".fig_record")) .fig_record(base, w_mm, h_mm, draw)   # Round 109: recorded for the figures script
  lc <- Sys.getlocale("LC_CTYPE"); if (!isTRUE(l10n_info()[["UTF-8"]])) { for (l in c("C.UTF-8", "en_US.UTF-8", "English_United States.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", l)))) break; on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", lc)), add = TRUE) }   # Round 109: symbols (phi, beta hat) drawn the same in every session
  w <- w_mm / 25.4; h <- h_mm / 25.4; out <- character()
  ok <- function(f, expr) { r <- tryCatch({ expr; TRUE }, error = function(e) FALSE); if (!r) { try(grDevices::dev.off(), silent = TRUE); if (file.exists(f)) unlink(f) }; if (r) out <<- c(out, f) }
  f <- paste0(base, ".pdf"); .fig_metric_reset(); ok(f, { if (.fig_cairo()) grDevices::cairo_pdf(f, width = w, height = h, family = "sans") else grDevices::pdf(f, width = w, height = h, useDingbats = FALSE); draw(); grDevices::dev.off() })   # Round 91: Cairo embeds the fonts, as journals ask
  if (isTRUE(unname(capabilities("cairo")))) {
    f <- paste0(base, ".svg"); .fig_metric_reset(); ok(f, { grDevices::svg(f, width = w, height = h); draw(); grDevices::dev.off() })
    f <- paste0(base, ".png"); .fig_metric_reset(); ok(f, { grDevices::png(f, width = w, height = h, units = "in", res = 300, type = "cairo"); draw(); grDevices::dev.off() })
  }
  out
}

# ---------- gathering what the report needs ----------
.est_df <- function(pay) {
  if (!length(pay$estimates)) return(NULL)
  d <- do.call(rbind, lapply(pay$estimates, function(e) data.frame(method = e$method, name = e$method_name, estimand = e$estimand, estimate = e$estimate, se = e$se, ci_lower = e$ci_lower, ci_upper = e$ci_upper, ci_level = e$ci_level, p = e$p, p_resolution = as.numeric(e$p_resolution %||% NA_real_), n_ivs = e$n_ivs_used, scale = e$scale %||% "beta", exposure_id = as.character(e$exposure_id %||% NA_character_), stringsAsFactors = FALSE)))
  d
}
.short <- c(loo = "Leave-one-out", steiger = "Steiger", ivw_steiger = "IVW (Steiger-filtered)", bwmr = "BWMR", grip = "MR-GRIP", mrmix = "MRMix", egger_boot = "MR-Egger (bootstrap)", pwm = "Pen. weighted median", smode_nome = "Simple mode (NOME)", wmode_nome = "Weighted mode (NOME)", rucker = "Rucker", sign = "Sign test", ivw_mre = "IVW", divw = "dIVW", maxlik = "Max. likelihood", mbe = "MBE", cml = "MR-cML", robust_ivw = "MR-Robust", pivw = "pIVW", cook = "Cook's D", mrclust = "MR-Clust", gsmr = "GSMR", grapple = "GRAPPLE", ppmr = "PPMR", mvgrapple = "MVMR GRAPPLE", studres = "Stud. residuals", egger_rob = "MR-Egger (robust)", hetpen = "HetPen", clr = "CLR", mvivwme = "MVMR-IVW (ME)", mvgmm = "MVMR GMM", raps = "MR-RAPS", conmix = "ConMix", lasso = "MR-Lasso", egger_simex = "MR-Egger (SIMEX)", egger = "MR-Egger", wm = "Weighted median", smed = "Simple median", mode = "Weighted mode", smode = "Simple mode", presso = "MR-PRESSO", radial = "Radial IVW", radegger = "Radial Egger")
pair_context <- function(h, p) {
  pay <- pair_payload(h$dir, p$folder); est <- .est_df(pay); ce <- est[est$estimand == "causal_effect", , drop = FALSE]
  v <- do.call(rbind, lapply(pay$variants, function(x) data.frame(id = x$id, bx = x$bx, by = x$by, sx = x$sx, sy = x$sy, keep = x$keep, stringsAsFactors = FALSE)))
  kept <- v[v$keep, , drop = FALSE]; Fs <- (kept$bx / kept$sx)^2
  fstats <- tryCatch({ f <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "01_data", "analysis_data_*.csv"))[1]; dd <- utils::read.csv(f, stringsAsFactors = FALSE)   # Round 109: from the analysis data (EAF, N)
    if (!"beta_x" %in% names(dd)) NULL else instrument_f_stats(dd[as.logical(dd$mr_keep), , drop = FALSE], p$ledger$x_scale %||% "unknown") }, error = function(e) NULL)
  st <- pay$settings %||% list(); dt <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "01_data", "analysis_data_*.csv"))
  dd <- if (length(dt)) utils::read.csv(dt[1], stringsAsFactors = FALSE) else NULL
  nx <- if (!is.null(dd$n_x)) stats::median(dd$n_x, na.rm = TRUE) else NA_real_; ny <- if (!is.null(dd$n_y)) stats::median(dd$n_y, na.rm = TRUE) else NA_real_
  rc <- file.path(h$dir, "03_pairs", p$folder, "02_scripts", "00_run_all.R"); sc <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "02_scripts", "[0-9]*_*.R"))
  srcname <- function(id, role) if (startsWith(id, "lib:")) tryCatch({ r <- library_recipe(sub("^lib:", "", id)); list(file = r$original$file, pmid = r$facts$pmid %||% r$source$pmid %||% "", author = r$source$author %||% r$facts$author %||% "", year = r$source$year %||% r$facts$year %||% "", ancestry = .gwas_pop(r) %||% r$facts$ancestry %||% "unknown",   # Round 80: the GWAS population OpenGWAS gives, not the clumping panel
       n = r$facts$n %||% "", ncase = r$facts$ncase %||% "", ncontrol = r$facts$ncontrol %||% "", build = r$facts$build %||% "", og = if (identical(r$source$type, "opengwas")) r$source else NULL, vcf = if (identical(r$source$type, "opengwas_vcf")) r$source else NULL,   # Round 100: complete GWAS downloaded from OpenGWAS
      
       mv = if (identical(r$source$type, "opengwas_mv")) r$source else NULL) }, error = function(e) list(file = "(file no longer in the library)", ancestry = "unknown", n = "")) else list(file = NA_character_, ancestry = "EUR", n = "")
  list(p = p, pay = pay, est = est, ce = ce, vars = v, kept = kept, mean_f = mean(Fs), min_f = min(Fs), n_weak = sum(Fs < 10), fstats = fstats, settings = st, nx = nx, ny = ny,
       simulated = isTRUE(p$simulated) || !startsWith(p$ex %||% "", "lib:"), ledger = p$ledger, pthr = p$pthr %||% 5e-8, src_x = srcname(p$ex, "exposure"), src_y = srcname(if (identical(p$med$role, "xm") && !is.null(p$med$m_out)) p$med$m_out else p$out, "outcome"), methods = unlist(p$methods), scripts = sc,
       overlap = tryCatch(overlap_get(h)[[overlap_key(p$ex %||% "", p$out %||% "")]], error = function(e) NULL),
       ties_sens = tryCatch({ ts <- Filter(function(x) identical(x$kind, "ties") && identical(x$parent, p$folder), ui_state_get(h$con, "sens") %||% list())
         if (length(ts) && dir.exists(file.path(h$dir, "03_pairs", ts[[1]]$folder))) list(pay = pair_payload(h$dir, ts[[1]]$folder), regions = ts[[1]]$regions %||% p$ledger$ties$regions, n_out = length(unlist(p$ledger$ties$out)), n_in = length(unlist(p$ledger$ties[["in"]])), folder = ts[[1]]$folder) else NULL }, error = function(e) NULL),
       outl = tryCatch({ o <- Filter(function(x) identical(x$parent, p$folder), ui_state_get(h$con, "outliers") %||% list())   # Round 42
         if (length(o)) { o <- o[[1]]; ch <- Filter(function(x) identical(x$set, o$chosen), o$scenarios); o$chosen_row <- if (length(ch)) ch[[1]] else NULL
           rm <- sort(unlist(o$chosen_row$ids)); sn <- Filter(function(x) identical(x$parent, p$folder) && length(rm) && identical(sort(unlist(x$remove)), rm), ui_state_get(h$con, "sens") %||% list())
           if (length(sn) && dir.exists(file.path(h$dir, "03_pairs", sn[[1]]$folder))) { o$chosen_pay <- pair_payload(h$dir, sn[[1]]$folder); e <- .est_df(o$chosen_pay); o$chosen_ivw <- e[e$estimand == "causal_effect" & e$method == "ivw_mre", , drop = FALSE] }
           o } else NULL }, error = function(e) NULL),
       eaf_sens = tryCatch({ es <- Filter(function(x) identical(x$kind, "eaf") && identical(x$parent, p$folder), ui_state_get(h$con, "sens") %||% list())   # Round 34
         if (length(es) && dir.exists(file.path(h$dir, "03_pairs", es[[length(es)]]$folder))) list(pay = pair_payload(h$dir, es[[length(es)]]$folder), remove = unlist(es[[length(es)]]$remove)) else NULL }, error = function(e) NULL),
       harm_sens = tryCatch({ hs <- Filter(function(x) grepl("^harm[123]$", x$kind %||% "") && identical(x$parent, p$folder) && dir.exists(file.path(h$dir, "03_pairs", x$folder)), ui_state_get(h$con, "sens") %||% list())   # Round 109
         lapply(hs, function(x) list(action = x$action, counts = x$counts, pay = pair_payload(h$dir, x$folder))) }, error = function(e) NULL))
}

#' Round 28: how equal P values were handled in local clumping, for the Methods sentence (inside its parenthesis).
# "6 of the primary instruments", "1 other": counts in the tie sentences (Round 29). A stronger lead can also remove its correlated
# neighbours, so the numbers left out and added can differ; they are not regions.
.n_inst <- function(n, plural, singular = sub("instruments", "instrument", plural)) sprintf("%d %s", as.integer(n), if (n == 1) singular else plural)
.tie_phrase <- function(C, T) {
  if (identical(C$tie_break, "equal P ordered by |beta/SE|")) return("; variants with equal P, such as P values too small to store, were ordered by the absolute value of beta/SE")
  if (is.null(T) || !isTRUE(T$tied > 0)) return("")
  sprintf("; %d candidates shared their P value with another (P values below the smallest value the source stores are equal), and among variants with equal P, PLINK kept the one that comes first in the reference panel, as in ieugwasr::ld_clump_local and OpenGWAS server clumping. %s",
    as.integer(T$tied), if (isTRUE(T$regions > 0)) sprintf("A sensitivity analysis that kept the variant with the larger |beta/SE| instead left out %s and added %s", .n_inst(length(unlist(T$out)), "of the primary instruments", "of the primary instruments"), .n_inst(length(unlist(T[["in"]])), "others", "other")) else "This changed no instrument")
}

# ---------- text ----------
.ci_phrase <- function(r, scale) { ratio <- scale %in% c("log_or"); f <- function(x) if (ratio) exp(x) else x   # Round 109: an OR keeps its log OR alongside
  ci <- if (is.finite(suppressWarnings(as.numeric(r$ci_lower))) && is.finite(suppressWarnings(as.numeric(r$ci_upper)))) sprintf("%s CI %s to %s", .pct(r$ci_level), .fnum(f(r$ci_lower)), .fnum(f(r$ci_upper))) else "no confidence interval"   # Round 110a: PPMR gives none
  sprintf("%s = %s (%s%s; P = %s)", if (ratio) "OR" else "beta", .fnum(f(r$estimate)), ci, if (ratio) sprintf("; log OR %s", .fnum(r$estimate)) else "", .fp_res(r$p, r$p_resolution)) }
#' Round 103: the packages behind the methods that were run, from each method's adapter (a version number or a pinned commit),
#' and the methods MR Studio implements itself. Before, the text named TwoSampleMR for every estimate.
# Round 104: methods whose estimator MR Studio computes itself although their script loads a package (Q-minimisation uses
# MendelianRandomization only to read the data and boot only for bootstrap intervals).
OWN_IMPL <- c(mvqhet = "MVMR Q-minimisation")
.own_name <- function(m) if (m %in% names(OWN_IMPL)) OWN_IMPL[[m]] else if (m %in% names(.short)) .short[[m]] else if (exists(".short_mv") && !is.na(.short_mv(m))) .short_mv(m) else m
software_pkgs <- function(used) {
  ad <- param_adapters(); ms <- intersect(used, names(ad)); pk <- list(); own <- character(0)
  for (m in ms) { pp <- ad[[m]]$packages %||% list(); if (!length(pp) || m %in% names(OWN_IMPL)) { own <- c(own, m); next }; for (n in names(pp)) if (is.null(pk[[n]])) pk[[n]] <- as.character(pp[[n]]) }
  list(ver = unname(vapply(names(pk), function(n) { v <- pk[[n]]; if (grepl("^[0-9a-f]{40}$", v)) sprintf("%s (commit %s)", n, substr(v, 1, 10)) else sprintf("%s (version %s)", n, v) }, "")), own = own)
}
software_sentence <- function(used) {
  sp <- software_pkgs(used); ver <- sp$ver; own <- sp$own
  if (!length(ver) && !length(own)) return("The estimates were computed")
  lst <- function(x) if (length(x) <= 1L) paste(x, collapse = "") else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
  nm <- .own_name
  paste0("The estimates were computed with ", lst(ver), if (length(own)) sprintf("%s MR Studio's own implementation of %s", if (length(ver) > 1L) ", and with" else if (length(ver)) " and with" else "", lst(vapply(own, nm, ""))) else "")
}
#' Round 103: results of the methods that give no causal estimate (the sign concordance test) or that choose a model (Rucker),
#' which were run and described in the Methods but reported nowhere.
method_result_text <- function(cx) {
  dg <- cx$pay$diagnostics %||% list(); g <- function(m, n) { x <- Filter(function(d) identical(d$method, m) && identical(d$name, n), dg); if (length(x)) x[[1]] else NULL }
  out <- character(0); s <- g("sign", "sign_concordance"); v <- if (is.null(s)) NA_real_ else suppressWarnings(as.numeric(s$value)); k <- if (is.null(s)) NA_real_ else suppressWarnings(as.numeric(s$df))
  if (is.finite(v) && is.finite(k) && k > 0) { same <- as.integer(round((v / 2 + 0.5) * k))
    out <- c(out, sprintf("Sign concordance test: %d of %d variants (%s%%) had associations with %s and %s in the same direction (binomial test P = %s).", same, as.integer(k), format(round(100 * same / k)), cx$p$exName, cx$p$outName, .fp(suppressWarnings(as.numeric(s$p))))) }
  r <- g("rucker", "rucker_model"); qd <- g("rucker", "q_diff")
  lab <- c(A = "the fixed-effect IVW model", B = "the random-effects IVW model", C = "the fixed-effect MR-Egger model", D = "the random-effects MR-Egger model")
  if (!is.null(r) && !is.null(r$note) && !is.na(r$note) && r$note %in% names(lab)) out <- c(out, sprintf("Rucker's model selection chose %s (model %s)%s.", lab[[r$note]], r$note,
    if (!is.null(qd) && is.finite(suppressWarnings(as.numeric(qd$p)))) sprintf("; the fall in Q from adding the MR-Egger intercept was %.2f (P = %s)", as.numeric(qd$value), .fp(as.numeric(qd$p))) else ""))
  out
}
#' Round 102: the limitation sentence about heterogeneity, for pairs whose IVW Cochran's Q has P below 0.05.
heterogeneity_limitation <- function(ctxs) {
  het <- vapply(ctxs, function(cx) { d <- Filter(function(d) identical(d$method, "ivw_mre") && identical(d$name, "cochran_q"), cx$pay$diagnostics %||% list())
    length(d) > 0 && is.numeric(d[[1]]$p) && isTRUE(d[[1]]$p < 0.05) }, TRUE)
  if (!any(het)) return(NULL)
  who <- if (length(ctxs) == 1L) "Cochran's Q showed" else sprintf("In %d of %d pairs, Cochran's Q showed", sum(het), length(ctxs))
  paste(who, "more heterogeneity between the variants than chance allows (P < 0.05), which usually reflects horizontal pleiotropy. The IVW estimate uses multiplicative random effects to allow for it, and the robust methods estimate the effect under weaker assumptions; heterogeneity was reported, not removed [Burgess 2023].")
}
#' The limitation sentences about sample overlap. Overlap is never calculated by MR Studio; the sentence says what the analyst declared.
overlap_limitations <- function(ctxs) {
  real <- Filter(function(cx) !isTRUE(cx$simulated), ctxs)
  bias <- "Overlap biases estimates toward the observational association, most strongly with weak instruments."
  untested <- "Sample overlap between the exposure and outcome studies was not tested. "
  if (!length(real)) return("The simulated exposure and outcome were generated independently, so their samples do not overlap. With real data, overlap biases estimates toward the observational association, most strongly with weak instruments.")
  nm <- function(cx) paste0(cx$p$exName %||% cx$p$ex %||% "exposure", " and ", cx$p$outName %||% cx$p$out %||% "outcome")
  ans <- vapply(real, function(cx) cx$overlap$answer %||% "unanswered", "")
  out <- character(0)
  none <- real[ans == "none"]; some <- real[ans == "some"]; rest <- real[!(ans %in% c("none", "some"))]
  det <- function(cx) paste0(nm(cx), if (!is.null(cx$overlap$pct)) sprintf(" (about %s%% of the people)", format(cx$overlap$pct)) else "", if (nzchar(cx$overlap$note %||% "")) paste0(" [", cx$overlap$note, "]") else "")
  if (length(none)) out <- c(out, paste0(untested, "The analyst stated that the exposure and outcome studies for ", paste(vapply(none, function(cx) paste0(nm(cx), if (nzchar(cx$overlap$note %||% "")) paste0(" [", cx$overlap$note, "]") else ""), ""), collapse = "; "), " do not share participants. This is a declaration and MR Studio did not verify it."))
  if (length(some)) out <- c(out, paste0("The analyst stated that the exposure and outcome studies share participants for: ", paste(vapply(some, det, ""), collapse = "; "), ". MR Studio did not correct for this. ", bias, " Consider an outcome from a different cohort or a method designed for overlapping samples."))
  if (length(rest)) out <- c(out, paste0(untested, "The cohorts of ", paste(vapply(rest, nm, ""), collapse = "; "), " were not checked by the analyst. ", bias))
  out
}
#' Steiger test numbers of one pair, read from the saved diagnostics (NULL when the test was not run).
# ---- Round 45: bidirectional MR. Two ordinary pairs with the roles swapped; this only reads them side by side. ----
BIDIR_METHODS_TEXT <- "To assess the direction of effect, each analysis was repeated with the roles of exposure and outcome reversed (bidirectional MR), using genome-wide significant instruments for the original outcome and the same methods and settings [Burgess 2023]."
BIDIR_LIMITS_TEXT <- "Bidirectional MR compares two separate analyses. It does not estimate the genetic correlation between the traits, and methods that model effects in both directions at once (bidirectional MR-cML, LHC-MR) were not run. Variants associated with both traits, or pleiotropy shared by the two sets of instruments, can produce apparent effects in both directions."
bidir_pairs <- function(ctxs) { out <- list(); n <- length(ctxs); if (n < 2) return(out)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) if (identical(ctxs[[i]]$p$ex, ctxs[[j]]$p$out) && identical(ctxs[[i]]$p$out, ctxs[[j]]$p$ex)) out[[length(out) + 1L]] <- c(i, j)
  out }
.kept_ids <- function(cx) vapply(Filter(function(v) isTRUE(v$keep), cx$pay$variants %||% list()), function(v) as.character(v$id), "")
bidir_verdict <- function(pf, pr) {
  sf <- isTRUE(pf < 0.05); sr <- isTRUE(pr < 0.05)
  if (sf && !sr) "forward" else if (!sf && sr) "reverse" else if (sf && sr) "both" else "neither" }
bidir_results_text <- function(a, b) {
  iv <- function(cx) cx$ce[cx$ce$method == "ivw_mre", , drop = FALSE]; fa <- iv(a); fb <- iv(b)
  X <- a$p$exName; Y <- a$p$outName
  if (!nrow(fa) || !nrow(fb)) return(sprintf("Bidirectional analysis (%s and %s): the IVW estimate is missing in one direction, so the two directions cannot be compared.", X, Y))
  v <- bidir_verdict(fa$p[1], fb$p[1]); sh <- length(intersect(.kept_ids(a), .kept_ids(b)))
  vt <- switch(v, forward = sprintf("The evidence supports an effect of %s on %s but not of %s on %s.", X, Y, Y, X),
                  reverse = sprintf("The evidence supports an effect of %s on %s but not of %s on %s, so an association in the first direction may reflect reverse causation.", Y, X, X, Y),
                  both = "Both directions gave P below 0.05, which can reflect a feedback loop between the traits or pleiotropy shared by the two sets of instruments; neither direction can be taken as causal on this evidence alone.",
                  neither = "Neither direction gave P below 0.05.")
  st <- function(cx) { f <- steiger_facts(cx); if (is.null(f)) NULL else if (f$ok) "supported" else "not supported" }
  sa <- st(a); sb <- st(b)
  sprintf("Bidirectional analysis (%s and %s): IVW %s -> %s %s; %s -> %s %s. %s %s%s", X, Y, X, Y, .ci_phrase(fa[1, ], fa$scale[1]), Y, X, .ci_phrase(fb[1, ], fb$scale[1]), vt,
    if (sh) sprintf("%d variant%s %s instruments in both directions.", sh, if (sh == 1) "" else "s", if (sh == 1) "was an" else "were") else "No variant was an instrument in both directions.",
    if (!is.null(sa) && !is.null(sb)) sprintf(" The Steiger test %s the direction %s -> %s and %s %s -> %s.", sa, X, Y, sb, Y, X) else "")
}

steiger_facts <- function(cx) {
  dg <- Filter(function(d) identical(d$method, "steiger"), cx$pay$diagnostics); if (!length(dg)) return(NULL)
  g <- function(n) { x <- Filter(function(d) identical(d$name, n), dg); if (length(x)) x[[1]] else NULL }
  d <- g("steiger_direction"); if (is.null(d)) return(NULL)
  bx <- g("steiger_exposure_binary"); by <- g("steiger_outcome_binary")   # absent in studies saved before Round 16 (continuous only then)
  list(ok = isTRUE(as.numeric(d$value) == 1), p = as.numeric(d$p), r2x = as.numeric(g("steiger_r2_exposure")$value), r2y = as.numeric(g("steiger_r2_outcome")$value),
       wrong = as.numeric(g("steiger_variants_wrong_way")$value), n = as.numeric(g("steiger_variants")$value),
       bin_x = !is.null(bx) && isTRUE(as.numeric(bx$value) == 1), bin_y = !is.null(by) && isTRUE(as.numeric(by$value) == 1),
       note_x = if (is.null(bx)) "" else as.character(bx$note), note_y = if (is.null(by)) "" else as.character(by$note))
}
#' TRUE when any pair ran the Steiger test with a binary trait (liability-scale variance explained, Lee 2012).
steiger_any_binary <- function(ctxs) any(vapply(ctxs, function(cx) { f <- steiger_facts(cx); !is.null(f) && (f$bin_x || f$bin_y) }, TRUE))
.prev_of <- function(note) sub(".*prevalence ", "", note)
# Round 95: a bidirectional study, and pairs whose exposure is a binary trait (read per unit of its log odds)
.ctx_bidir <- function(ctxs) any(vapply(ctxs, function(cx) isTRUE(cx$p$rev), TRUE))
.ctx_bin_x <- function(ctxs) Filter(function(cx) { id <- cx$p$ex; if (!is.character(id) || length(id) != 1L) return(FALSE); startsWith(id, "lib:") && identical(tryCatch(library_recipe(sub("^lib:", "", id))$facts$kind, error = function(e) NULL), "binary") }, ctxs)
binary_exposure_limitations <- function(ctxs) { b <- .ctx_bin_x(ctxs); if (!length(b)) return(NULL)
  sprintf("%s %s binary: %s estimate is the change in the outcome per unit increase in the log odds of the exposure (multiply by 0.693 for a doubling of its odds), and it reflects the effect of the liability to the trait rather than of having it, so it should be read with care [Burgess and Labrecque 2018].",
    paste(unique(vapply(b, function(cx) cx$p$exName %||% cx$p$ex, "")), collapse = " and "), if (length(unique(vapply(b, function(cx) cx$p$exName %||% cx$p$ex, ""))) == 1L) "is" else "are", if (length(b) == 1L) "its" else "each") }
steiger_limitations <- function(ctxs) {
  f <- lapply(ctxs, steiger_facts); ran <- !vapply(f, is.null, TRUE)
  if (!any(ran) && .ctx_bidir(ctxs)) return("The Steiger directionality test was not run. The direction of the effect was examined instead by estimating each pair in both directions (reverse MR), and jointly with bidirectional MR-cML where it ran.")
  if (!any(ran)) return("Directionality was not tested (the Steiger directionality test was not run), so reverse causation is not excluded by this analysis.")
  bad <- vapply(seq_along(f), function(i) ran[i] && !(f[[i]]$ok && is.finite(f[[i]]$p) && f[[i]]$p < 0.05), TRUE)
  c("The Steiger directionality test [Hemani 2017] compares the variance the variants explain in the exposure and in the outcome. It assumes the variants act on the outcome only through the exposure and is sensitive to measurement error and horizontal pleiotropy. It does not exclude confounding or reverse causation.",
    {
      pv <- unlist(lapply(ctxs[ran], function(cx) { fx <- steiger_facts(cx); c(if (fx$bin_x) sprintf("%s %s", cx$p$exName, .prev_of(fx$note_x)), if (fx$bin_y) sprintf("%s %s", cx$p$outName, .prev_of(fx$note_y))) }))
      if (length(pv)) sprintf("For binary traits the variance explained was computed on the liability scale [Lee 2012] and depends on the population prevalence entered by the analyst (%s); a different prevalence can change the result.", paste(pv, collapse = "; ")) else NULL
    },
    if (any(bad)) sprintf("For %s the test did not clearly support the assumed direction.", paste(vapply(ctxs[bad], function(cx) sprintf("%s -> %s", cx$p$exName, cx$p$outName), ""), collapse = "; ")),
    if (any(!ran)) sprintf("The test was not run for %s.", paste(vapply(ctxs[!ran], function(cx) sprintf("%s -> %s", cx$p$exName, cx$p$outName), ""), collapse = "; ")))
}
AVAIL_SIM <- "No external data were used: the traits are simulated by MR Studio."
AVAIL_SCRIPTS <- "The analysis scripts, the harmonised data tables used in every calculation, all result tables and the replay check are provided in the reproducibility package of this study; every figure can be redrawn from the figures script of its pair (05_figures/make_figures.R)."
#' Round 99: the Data availability sentence; OpenGWAS dataset IDs and imported files are named separately when a study has both.
.ogw_both <- function(sx, sy) (!is.null(sx$og) || !is.null(sx$vcf)) && (!is.null(sy$og) || !is.null(sy$vcf))   # Round 100
.avail_sentence <- function(src) { og <- grepl("^OpenGWAS:", src)
  sprintf("Summary statistics were %s.", paste(c(if (any(og)) sprintf("obtained from the public OpenGWAS database (api.opengwas.io), dataset IDs: %s", paste(sub("^OpenGWAS:", "", src[og]), collapse = "; ")),
    if (any(!og)) sprintf("read from the files %s [Add the public source, accession number or DOI of each file.]", paste(src[!og], collapse = "; "))), collapse = "; and ")) }
#' Round 50: univariable pairs go through the text below; multivariable ones get their own paragraphs (mv_add_text).
report_text <- function(ctxs, study_name) {
  mvx <- Filter(function(cx) isTRUE(cx$mv), ctxs); uni <- Filter(function(cx) !isTRUE(cx$mv), ctxs)
  txt <- if (length(uni)) .report_text_uni(uni, study_name) else list(methods = sub("run in MR Studio,", sprintf("run in MR Studio %s,", .mrs_ver()), MV_INTRO_TEXT, fixed = TRUE), results = character(), limitations = "The three core MR assumptions (relevance, independence from confounders, and no effect of the variants on the outcome except through the exposures) cannot be proven from summary data.",
    availability = character(), legends = character(), supp = list(), supp_out = list())
  if (length(mvx)) txt <- mv_add_text(txt, mvx)
  for (cx in uni) for (hs in cx$harm_sens %||% list()) { iv <- Filter(function(e) identical(e$method, "ivw_mre") && identical(e$estimand, "causal_effect"), hs$pay$estimates %||% list())   # Round 109: harmonisation action as a sensitivity analysis
    txt$methods <- c(txt$methods, sprintf("As a sensitivity analysis, %s -> %s was harmonised again with TwoSampleMR action %d (%s; primary analysis action %s): %s flipped, %s palindromic dropped, %s instruments kept%s.", cx$p$exName, cx$p$outName, as.integer(hs$action), HARM_ACTIONS[[as.character(hs$action)]], cx$ledger$harmonise_action %||% "2", hs$counts$flipped %||% "?", hs$counts$dropped_palindromic %||% "?", hs$counts$kept %||% "?",
      if (length(iv)) paste0("; IVW ", .ci_phrase(as.data.frame(iv[[1]][c("estimate", "ci_lower", "ci_upper", "ci_level", "p")]), iv[[1]]$scale)) else "")) }
  sa <- Filter(function(cx) "egger_simex" %in% (cx$methods %||% character()) && isTRUE(simex_i2gx(cx) < SIMEX_I2GX), uni)   # Round 109: SIMEX added by itself
  if (length(sa)) txt$methods <- c(txt$methods, sprintf("MR-Egger with SIMEX correction was added because I2GX was below 0.9 (%s), which means the uncertainty in the variant-exposure associations dilutes the MR-Egger estimate toward the null (violation of the NOME assumption) [Bowden 2016, Int J Epidemiol].",
    paste(vapply(sa, function(cx) sprintf("%s -> %s: %.2f", cx$p$exName, cx$p$outName, simex_i2gx(cx)), ""), collapse = "; ")))
  sci_thr_text(txt)   # Round 109: P cut points in scientific form (5e-6, not 0.000005 or 5e-06)
}
.report_text_uni <- function(ctxs, study_name) {
  txt <- list(methods = character(), results = character(), limitations = character(), availability = character(), legends = character())
  used <- unique(unlist(lapply(ctxs, function(x) x$methods))); anysim <- any(vapply(ctxs, `[[`, TRUE, "simulated")); anyimp <- any(!vapply(ctxs, `[[`, TRUE, "simulated"))
  mname <- c(ivw_steiger = "IVW after Steiger filtering", bwmr = "Bayesian weighted MR (BWMR)", grip = "MR-GRIP", mrmix = "MRMix", egger_boot = "MR-Egger with bootstrap standard errors", pwm = "penalised weighted median", smode_nome = "simple mode (NOME)", wmode_nome = "weighted mode (NOME)", rucker = "Rucker model selection", sign = "sign concordance test", ivw_mre = "inverse-variance weighted (IVW)", divw = "debiased inverse-variance weighted (dIVW)", maxlik = "two-sample maximum likelihood", mbe = "mode-based estimate (MBE)", cml = "constrained maximum likelihood (MR-cML)", robust_ivw = "robust IVW with penalised weights (MR-Robust)", pivw = "penalised inverse-variance weighted (pIVW)", mrclust = "MR-Clust clustering of the variants", gsmr = "generalised summary-data-based MR (GSMR)", grapple = "GRAPPLE", ppmr = "probabilistic MR (PPMR)", cook = "Cook's distance from the IVW regression", studres = "studentized residuals from the IVW regression", egger_rob = "MR-Egger with robust regression and penalised weights", hetpen = "heterogeneity-penalised model averaging", clr = "conditional likelihood ratio method", mvivwme = "MVMR-IVW allowing for measurement error in the exposure associations", mvgmm = "multivariable generalised method of moments (GMM)", raps = "MR-RAPS", conmix = "contamination mixture", lasso = "MR-Lasso", egger_simex = "MR-Egger with SIMEX correction", egger = "MR-Egger regression", wm = "weighted median", smed = "simple median", mode = "weighted mode", smode = "simple mode", presso = "MR-PRESSO", radial = "radial IVW", radegger = "radial MR-Egger")
  txt$methods <- c(
    paste0("We performed a two-sample Mendelian randomization (MR) analysis using summary-level genetic associations, run in MR Studio ", .mrs_ver(), ", an application that writes every analysis step and setting into R scripts that can be re-run. The analysis plan follows the structure recommended in current MR guidelines [Burgess 2023]; reporting against the STROBE-MR checklist [Skrivankova 2021] has not been checked automatically and must be completed by the authors."))
  for (cx in ctxs) {
    p <- cx$p; st <- cx$settings
    if (cx$simulated) {
      txt$methods <- c(txt$methods, sprintf("Exposure %s and outcome %s. The exposure and the outcome are simulated demonstration traits generated by MR Studio. They do not describe any real population, and no result from them may be reported as a finding.", p$exName, p$outName))
    } else {
      sx <- cx$src_x; sy <- cx$src_y
      nstr <- function(s) if (nzchar(as.character(s$n))) sprintf(", sample size %s", format(as.numeric(s$n), big.mark = ",", scientific = FALSE)) else ""
      side <- function(s, what) if (!is.null(s$vcf)) sprintf("were taken from the complete summary statistics of OpenGWAS dataset %s (file %s, downloaded %s; ancestry %s%s%s%s)", s$vcf$id, s$file, substr(s$vcf$downloaded %||% "", 1, 10), s$ancestry, nstr(s),   # Round 100
                                  if (nzchar(as.character(s$ncase %||% "")) && nzchar(as.character(s$ncontrol %||% ""))) sprintf(", %s cases and %s controls", format(as.numeric(s$ncase), big.mark = ",", scientific = FALSE), format(as.numeric(s$ncontrol), big.mark = ",", scientific = FALSE)) else "",
                                  if (as.character(s$build %||% "") %in% c("37", "38")) sprintf(", genome build GRCh%s", s$build) else "")
                                else if (!is.null(s$og)) sprintf("were obtained from the OpenGWAS database (dataset %s; ancestry %s%s; accessed %s)", s$og$id, s$ancestry, nstr(s), substr(s$og$downloaded_at %||% "", 1, 10))
                                else sprintf("were read from the file \"%s\" (ancestry declared as %s%s%s%s)", s$file, { an <- unname(LDREF_POP_NAMES[s$ancestry]); if (length(an) && !is.na(an)) an else s$ancestry }, nstr(s),
                                  if (nzchar(as.character(s$ncase %||% "")) && nzchar(as.character(s$ncontrol %||% ""))) sprintf(", %s cases and %s controls", format(as.numeric(s$ncase), big.mark = ",", scientific = FALSE), format(as.numeric(s$ncontrol), big.mark = ",", scientific = FALSE)) else "",
                                  if (as.character(s$build %||% "") %in% c("37", "38")) sprintf(", genome build GRCh%s", s$build) else "")   # Round 99
      txt$methods <- c(txt$methods, sprintf("Exposure %s and outcome %s. Genetic associations with the exposure %s, and associations with the outcome %s. %s",
        p$exName, p$outName, side(sx, "were"), side(sy, "were"),
        if ((!is.null(sx$og) || !is.null(sx$vcf)) && (!is.null(sy$og) || !is.null(sy$vcf))) "Ancestry and the units of the effect sizes are as reported by OpenGWAS; the genome build was not verified automatically." else "Ancestry, genome build and the units of the effect sizes are as declared by the analyst when files were imported; they were not verified automatically."))
      L <- cx$ledger
      if (!is.null(L) && !is.null(sx$og) && is.null(L$clump)) {
        pr <- sx$og$params; po <- if (!is.null(sy$og)) sy$og$params else NULL
        txt$methods <- c(txt$methods, sprintf("Instruments were the genome-wide significant variants of the exposure dataset (P < %s), %s (r\u00b2 < %s within %s kb, 1000 Genomes %s reference panel; ieugwasr %s). %s Of %d instruments, %d %s absent from the outcome%s. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure, %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used.",
          format(pr$pval), if (identical(pr$where, "server_live")) "clumped live on the OpenGWAS server with ieugwasr::tophits (force_server = TRUE)" else if (ogw_stored_clumps(pr)) "as returned by ieugwasr::tophits from the stored, pre-computed clumps of OpenGWAS (default settings, preclumped = 1)" else "clumped on the OpenGWAS server with ieugwasr::tophits",
          format(pr$r2), format(pr$kb, big.mark = ","), pr$pop, sx$og$packages$ieugwasr %||% "unknown",
          if (!is.null(po)) sprintf("The same variants were extracted from the outcome dataset with TwoSampleMR::extract_outcome_data (%s; TwoSampleMR %s).", if (isTRUE(po$proxies)) sprintf("proxies allowed at r\u00b2 of at least %s, palindromic variants aligned by allele frequency only when the minor allele frequency was below %s", format(po$rsq), format(po$maf_threshold)) else "no proxies", sy$og$packages$TwoSampleMR %||% "unknown") else "The outcome associations were read from an imported file; instruments missing from it were not replaced by proxy variants.",
          L$instruments, L$missing_in_outcome, if (L$missing_in_outcome == 1) "was" else "were", if (!is.null(L$by_proxy) && L$by_proxy > 0) sprintf(" (a further %d %s found only through proxy variants)", L$by_proxy, if (L$by_proxy == 1) "was" else "were") else "", L$flipped, L$dropped_palindromic, L$dropped_mismatch, L$kept))
        if (isTRUE(L$ties$server) && isTRUE(L$ties$tied > 0)) txt$methods <- c(txt$methods, sprintf("%d instruments shared their P value with another instrument (%s). Because clumping was done on the OpenGWAS server, which of several variants with equal P was kept in a region could not be checked.", as.integer(L$ties$tied), p_store_phrase(L$ties)))
      } else       if (!is.null(L) && !is.null(L$clump)) {
        C <- L$clump
        og_local <- !is.null(sx$og); po <- if (og_local && !is.null(sy$og)) sy$og$params else NULL
        if (og_local) txt$methods <- c(txt$methods, sprintf("Every variant of the exposure dataset with P < %s was downloaded from OpenGWAS without clumping (ieugwasr::tophits with clump = 0; ieugwasr %s) and clumped on the analyst's computer at download, as described next. Variants with alleles other than single bases, a missing or non-positive standard error, or a repeated ID were removed before clumping.%s",
          format(sx$og$params$pval), sx$og$packages$ieugwasr %||% "unknown",
          if (!is.null(po)) sprintf(" The independent instruments described next were looked up in the outcome dataset with TwoSampleMR::extract_outcome_data (%s; TwoSampleMR %s).", if (isTRUE(po$proxies)) sprintf("proxies allowed at r\u00b2 of at least %s, palindromic variants aligned by allele frequency only when the minor allele frequency was below %s", format(po$rsq), format(po$maf_threshold)) else "no proxies", sy$og$packages$TwoSampleMR %||% "unknown") else ""))
        txt$methods <- c(txt$methods, sprintf("Instruments were selected from the exposure summary statistics on the analyst's computer: the %d variants with P < %s were clumped with %s [Chang 2015] using the %s 1000 Genomes Phase 3 reference panel [1000 Genomes 2015] distributed by the MRC Integrative Epidemiology Unit%s (r\u00b2 < %s within %s kb, keeping the variant with the smallest P in each region%s). %d variants were removed as correlated with a stronger variant and %d were removed because they are not in the reference panel, leaving %d independent instruments, of which %d were absent from the outcome%s. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure (swapped or on the opposite strand), %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used. Every candidate variant, kept or removed, is listed with its index variant in the pair's instrument_selection.csv.",
          as.integer(L$below_threshold), format(L$threshold), C$plink %||% "PLINK 1.9", C$pop, if (isTRUE(C$panel_matches_pin)) sprintf(" (files identical to its archive 1kg.v3.tgz, SHA-256 %s)", LD_PANEL$sha256) else if (isFALSE(C$panel_matches_pin)) " (a copy whose files differ from the published archive)" else "", format(C$r2), format(C$kb, big.mark = ","), .tie_phrase(C, L$ties), as.integer(C$clumped), as.integer(C$not_in_panel), as.integer(C$kept), as.integer(L$missing_in_outcome), if (!is.null(L$by_proxy) && L$by_proxy > 0) sprintf(" (a further %d %s found only through proxy variants)", L$by_proxy, if (L$by_proxy == 1) "was" else "were") else "", as.integer(L$flipped), as.integer(L$dropped_palindromic), as.integer(L$dropped_mismatch), as.integer(L$kept)))
      } else if (!is.null(L)) txt$methods <- c(txt$methods, sprintf("Instruments were the variants in the exposure file with P < %s. The exposure file was declared by the analyst to contain only independent (clumped) variants; MR Studio did not perform clumping. Of %d candidate variants, %d were absent from the outcome file. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure (swapped or on the opposite strand), %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used.", format(L$threshold), L$candidates, L$missing_in_outcome, L$flipped, L$dropped_palindromic, L$dropped_mismatch, L$kept))
      if (!is.null(L) && identical(L$strand, "forward")) txt$methods <- .strand_forward_text(txt$methods, L, .ogw_both(sx, sy))   # Round 26
      if (!is.null(L$eaf_check)) txt$methods <- c(txt$methods, eaf_check_text(L$eaf_check, list(x = sx$ancestry, y = sy$ancestry)))   # Round 33
    }
    fl <- if (isTRUE(st$min_f > 0)) sprintf(" Variants with a single-variant F statistic below %s were left out.", st$min_f) else ""
    txt$methods <- c(txt$methods, sprintf("The strength of the instruments was summarised by the single-variant F statistic, (beta_x / se_x)^2: mean %.1f, minimum %.1f, with %d of %d variants below 10.%s", cx$mean_f, cx$min_f, cx$n_weak, nrow(cx$kept), fl))
    fx <- tryCatch(f_stats_text(cx$ledger$fstats %||% cx$fstats), error = function(e) "")   # Round 109: the three F statistics, and the MAF filter
    if (nzchar(fx)) txt$methods <- c(txt$methods, paste("Instrument strength, three F statistics:", fx))
    if (!is.null(cx$ledger$maf)) txt$methods <- c(txt$methods, maf_text(cx$ledger$maf))
  }
  ivw_model <- c(mre_floor = "multiplicative random effects, with the residual standard error not allowed below 1 so that the standard error is never smaller than that of a fixed-effect model (TwoSampleMR mr_ivw), the model the guidelines recommend as the default",
                 mre = "multiplicative random effects without a lower limit on the residual standard error (TwoSampleMR mr_ivw_mre)", fixed = "a fixed-effect model")
  st1 <- ctxs[[1]]$settings
  desc <- vapply(used, function(m) switch(m,
    ivw_mre = sprintf("The primary estimate was the %s estimator with %s [%s].", mname[["ivw_mre"]], ivw_model[[st1$ivw_model %||% "mre"]], if (identical(st1$ivw_model, "mre_floor")) "Burgess 2013, Burgess 2023, Hemani 2018" else "Burgess 2013, Hemani 2018"),
    egger = "MR-Egger regression was used to test for directional pleiotropy through its intercept and to estimate the effect allowing for it; variants were oriented to a positive exposure effect and confidence intervals use a t distribution with k - 2 degrees of freedom [Bowden 2015].",
    wm = sprintf("The weighted median estimator, which is consistent when at least half of the weight comes from valid instruments, was computed with a bootstrapped standard error (%s replications) [Bowden 2016].", format(st1$nboot %||% 1000)),
    bwmr = "Bayesian weighted MR (BWMR), which allows for weak polygenic pleiotropic effects and down-weights outlying variants, was computed by variational Bayes [Zhao 2020], with MR Studio's implementation of the published algorithm run to convergence.",
    ivw_steiger = "As a sensitivity analysis, variants that explained more variance in the outcome than in the exposure (Steiger filtering) were removed and the IVW estimate was recomputed on the remaining variants [Hemani 2017]; because filtering uses the outcome data to select instruments, this estimate is reported alongside, not instead of, the primary analysis.",
    grip = sprintf("MR-GRIP, a version of MR-Egger whose estimate and intercept do not depend on which allele is coded as the effect allele%s, was computed [Dudbridge 2025].", if (identical(st1$grip_variant, "weak")) ", adjusted for weak instruments," else ""),
    mrmix = "MRMix, a mixture-model estimator robust to a majority of invalid variants, was computed on standardised associations (using allele frequencies and sample sizes) and converted to the per-unit scale of the other methods with the median implied standard deviations [Qi 2019].",
    egger_boot = sprintf("MR-Egger regression with bootstrap standard errors (%s iterations) was computed [Bowden 2015].", format(st1$eb_nboot %||% 1000)),
    pwm = sprintf("The penalised weighted median, which down-weights variants with heterogeneous ratio estimates (penalty constant %s, %s bootstrap iterations), was computed [Burgess 2017].", format(st1$pwm_penk %||% 20), format(st1$pwm_nboot %||% 1000)),
    smode_nome = sprintf("The simple mode-based estimate with ratio standard errors assuming no measurement error in the exposure associations (NOME; bandwidth multiplier %s) was computed [Hartwig 2017].", format(st1$smn_phi %||% 1)),
    wmode_nome = sprintf("The weighted mode-based estimate with ratio standard errors assuming no measurement error in the exposure associations (NOME; bandwidth multiplier %s) was computed [Hartwig 2017].", format(st1$wmn_phi %||% 1)),
    rucker = sprintf("Heterogeneity was compared between the IVW model (Cochran's Q) and the MR-Egger model (Rucker's Q'), and the model was selected within Rucker's framework%s [Bowden 2018].", switch(st1$rk_variant %||% "point", point = "", bootstrap = sprintf(", with the %s of bootstrap estimates", st1$rk_stat %||% "median"), jackknife = sprintf(", with the %s of jackknife estimates", st1$rk_stat %||% "median"))),
    sign = "A sign concordance test assessed whether the variant-exposure and variant-outcome associations had the same direction more often than expected by chance (binomial test) [Hemani 2018].",
    maxlik = sprintf("The two-sample maximum likelihood estimator, which allows for uncertainty in the variant-exposure associations, was computed (model: %s%s) [Burgess 2013].", st1$ml_model %||% "default", if (isTRUE((st1$ml_psi %||% 0) != 0)) sprintf("; correlation between exposure and outcome associations %s", format(st1$ml_psi)) else ""),
    mbe = sprintf("The mode-based estimate (MendelianRandomization, %s mode, bandwidth multiplier %s, %s bootstrap iterations), which is consistent when the largest group of variants with the same ratio estimate is valid, was computed [Hartwig 2017].", st1$mbe_weighting %||% "weighted", format(st1$mbe_phi %||% 1), format(st1$mbe_iter %||% 10000)),
    cml = sprintf("Constrained maximum likelihood with model averaging by BIC%s (MR-cML), robust to correlated and uncorrelated pleiotropy, was computed with the smaller of the exposure and outcome GWAS sample sizes [Xue 2021].", if (identical(st1$cml_dp, "ma")) "" else " and data perturbation"),
    robust_ivw = sprintf("MR-Robust, the inverse-variance weighted estimate %s, was used as a sensitivity analysis robust to outlying variants [Rees 2019b].", switch(st1$rob_variant %||% "both", both = "with robust regression (MM-estimation) and penalised weights", robust = "with robust regression (MM-estimation)", penalized = "with penalised weights")),
    pivw = sprintf("The penalised inverse-variance weighted (pIVW) estimator, robust to many weak instruments, was computed with penalty %s and %s [Xu 2023].", format(st1$piv_lambda %||% 1), if (identical(st1$piv_ci, "wald")) "normal intervals" else "bootstrap Fieller intervals"),
    grapple = sprintf("GRAPPLE, a profile-likelihood estimator robust to weak instruments, balanced horizontal pleiotropy and sample overlap, was computed with %s loss%s [Wang 2021].", st1$gr_loss %||% "tukey", if (isTRUE((st1$gr_rho_xy %||% 0) != 0)) sprintf(" and an exposure-outcome error correlation of %s", format(st1$gr_rho_xy)) else ""),
    ppmr = "Probabilistic MR (PPMR), a likelihood method that uses correlated instruments with their LD from the 1000 Genomes reference panel and tests and controls for horizontal pleiotropy, was computed; it gives likelihood-ratio P values but no confidence interval [Yuan 2020].",
    gsmr = sprintf("Generalised summary-data-based MR (GSMR; gsmr2 %s), which uses every instrument together with their LD correlation from the 1000 Genomes reference panel (PLINK) and removes pleiotropic variants with the HEIDI-outlier test%s, was computed [Zhu 2018].", substr(PIN_GSMR2, 1, 7), if (identical(st1$gsmr_heidi, "off")) " switched off" else ""),
    mrclust = sprintf("Variants were clustered by their ratio estimates with MR-Clust (mixture model with null and junk clusters, number of clusters chosen by BIC; a variant was assigned when its allocation probability was at least %s and clusters with fewer than %s variants were not reported), to detect groups of variants with different causal estimates [Foley 2021].", format(st1$mc_prob %||% 0.8), format(st1$mc_min %||% 4)),
    cook = sprintf("Influential variants were identified with Cook's distance from the IVW regression (flagged above %s / k, k the number of variants) [Cook 1977].", format(st1$cook_mult %||% 4)),
    studres = sprintf("Outlying variants were identified with the externally studentized residuals of the IVW regression (flagged when |r| > %s).", st1$stud_cut %||% "2"),
    egger_rob = sprintf("MR-Egger %s was used as a sensitivity analysis less affected by outlying variants [Rees 2019b].", switch(st1$egr_variant %||% "both", both = "with robust regression (MM-estimation) and penalised weights", robust = "with robust regression (MM-estimation)", penalized = "with penalised weights")),
    hetpen = sprintf("Heterogeneity-penalised model averaging, which is consistent when a plurality of the variants are valid instruments, was computed with a prior probability of validity of %s; its confidence interval was found on a grid and can consist of several ranges [Burgess 2018].", format(st1$hp_prior %||% 0.5)),
    clr = "Weak-instrument robust tests were computed with the conditional likelihood ratio method, giving Anderson-Rubin, Kleibergen and conditional likelihood ratio confidence intervals (no point estimate), using the median sample sizes of the exposure and outcome GWAS [Wang 2022].",
    divw = sprintf("The debiased inverse-variance weighted (dIVW) estimator, which removes the bias of IVW toward the null from weak instruments, was computed with %s [Ye 2021].", if (identical(st1$dispersion, "none")) "a standard error that assumes no pleiotropy" else "a standard error that allows for balanced pleiotropy"),
    raps = sprintf("The robust adjusted profile score (MR-RAPS) estimator, which accounts for weak instruments and limits the influence of outlying variants, was computed with %s loss and %s [Zhao 2020].", c(huber = "Huber", tukey = "Tukey biweight", l2 = "squared-error")[[st1$loss %||% "huber"]], if (identical(st1$raps_dispersion, "none")) "no allowance for pleiotropy" else "a variance component for balanced pleiotropy"),
    egger_simex = sprintf("As a sensitivity analysis for regression dilution of MR-Egger (NOME violation, measured by I2GX), the weighted MR-Egger regression was corrected by simulation extrapolation (SIMEX; simex package, quadratic extrapolation, %s simulations per step, jackknife standard errors, t interval with k - 2 degrees of freedom) [Bowden 2016, Int J Epidemiol].", format(st1$simex_b %||% 100)),
    conmix = sprintf("The contamination mixture method, which is consistent when the largest group of variants with similar causal estimates consists of valid instruments, was computed by profile likelihood on a grid of step %s with %s for the spread of invalid estimates; it gives an interval but no standard error, and the interval can consist of more than one range [Burgess 2020].", format(st1$ci_step %||% 0.01), if (isTRUE((st1$psi %||% 0) > 0)) sprintf("psi = %s", format(st1$psi)) else "the package default (1.5 times the standard deviation of the ratio estimates)"),
    lasso = sprintf("MR-Lasso, which penalises variants with outlying estimates and is consistent when a majority of the variants are valid, was used to identify and exclude outlying variants; the post-lasso IVW estimate on the remaining variants is reported, with a %s interval [Rees 2019].", if (identical(st1$lasso_dist, "t-dist")) "t-distribution" else "normal"),
    smed = sprintf("The simple median estimator, which is consistent when at least half of the variants are valid instruments, was computed with a bootstrapped standard error (%s replications) [Bowden 2016].", format(st1$nboot %||% 1000)),
    mode = sprintf("The weighted mode estimator, which is consistent when the largest group of variants with similar ratio estimates is valid, was computed with a bandwidth multiplier phi = %s and a bootstrapped standard error (%s replications) [Hartwig 2017].", format(st1$phi %||% 1), format(st1$nboot %||% 1000)),
    smode = sprintf("The simple mode estimator (the unweighted counterpart of the weighted mode) was computed with a bandwidth multiplier phi = %s and a bootstrapped standard error (%s replications) [Hartwig 2017].", format(st1$phi %||% 1), format(st1$nboot %||% 1000)),
    presso = sprintf("MR-PRESSO (global test of horizontal pleiotropy, outlier test and distortion test) was run with %s simulated datasets and a significance threshold of %s. Variants it flags are reported but were not removed from any analysis [Verbanck 2018].", format(st1$nbdist %||% 1000), format(st1$signif %||% 0.05)),
    radial = sprintf("Radial IVW was used to flag outlying variants through their contribution to Cochran's Q statistic (outlier level %s, %s weights); flagged variants are reported but were not removed [Bowden 2018].", .radial_level_txt(st1), c("1" = "first-order", "2" = "second-order", "3" = "modified second-order")[[as.character(st1$weights %||% "1")]]),
    loo = "The influence of single variants was assessed with a leave-one-out analysis (the IVW estimate recomputed with each variant omitted in turn, with the same effects model), and the single-variant (Wald ratio) estimates were shown in a funnel plot, where asymmetry can indicate directional pleiotropy [Burgess 2023].",
    steiger = paste0("The Steiger directionality test, run together with the IVW analysis, compared the variance the variants explain in the exposure with the variance they explain in the outcome, summed over variants (TwoSampleMR::mr_steiger; measurement reliability set to 1 unless stated otherwise). For a continuous trait the variance explained was computed from each variant's effect, standard error and sample size",
      if (steiger_any_binary(ctxs)) "; for a binary trait it was computed on the liability scale from the log odds ratio, the effect allele frequency, the numbers of cases and controls and the population prevalence (TwoSampleMR::get_r_from_lor, logit model) [Lee 2012]" else "",
      ". It gives no causal estimate [Hemani 2017]."),
    radegger = sprintf("Radial MR-Egger was used to flag outlying variants and to estimate the effect allowing for directional pleiotropy (outlier level %s, %s weights); flagged variants are reported but were not removed [Bowden 2018].", .radial_level_txt(st1), c("1" = "first-order", "2" = "second-order", "3" = "modified second-order")[[as.character(st1$weights %||% "1")]]), ""), "")
  if (any(vapply(ctxs, function(cx) !is.null(cx$outl), TRUE))) desc <- c(desc, outl = OUTLIER_METHODS_TEXT, if (any(vapply(ctxs, function(cx) length(cx$outl$lookup$rows) > 0, TRUE))) c(outl_lookup = OUTLIER_LOOKUP_METHODS))   # Round 42
  txt$methods <- c(txt$methods, unname(desc[nzchar(desc)]), sprintf("Confidence intervals are %s intervals. %s in %s. Bootstrap-based methods use a fixed seed recorded in each script, and every pair was re-run from its saved scripts in a fresh R session; the replay check was %s.",
    .pct(ctxs[[1]]$ce$ci_level[1] %||% 0.95), software_sentence(used), R.version.string, paste(unique(vapply(ctxs, function(x) x$pay$replay$status %||% "not run", "")), collapse = " / ")))
  for (cx in ctxs) {
    ce <- cx$ce; if (!nrow(ce)) next; sc <- ce$scale[1]; ratio <- sc == "log_or"
    unit <- if (ratio) sprintf("odds of %s per unit of %s", cx$p$outName, cx$p$exName) else sprintf("units of %s per unit of %s (beta; often per SD, as in the GWAS)", cx$p$outName, cx$p$exName)
    lines <- vapply(seq_len(nrow(ce)), function(i) sprintf("%s: %s", .short[[ce$method[i]]], .ci_phrase(ce[i, ], sc)), "")
    ic <- utils::head(cx$est[cx$est$estimand == "pleiotropy_intercept" & cx$est$method %in% "egger", , drop = FALSE], 1)   # Round 86: MR-Egger's own (radial MR-Egger, SIMEX and others have one too; with several the text stopped)
    dg <- Filter(function(d) d$method == "ivw_mre" && d$name == "cochran_q", cx$pay$diagnostics)
    same <- sum(sign(ce$estimate) == sign(ce$estimate[ce$method == "ivw_mre"][1]), na.rm = TRUE)
    txt$results <- c(txt$results, sprintf("%s -> %s (%d instruments, effect expressed as %s). %s.%s%s %d of %d methods agree on the direction of the effect.", cx$p$exName, cx$p$outName, ce$n_ivs[1], unit, paste(lines, collapse = "; "),
      if (length(dg)) sprintf(" Cochran's Q for the IVW estimate was %.1f (P = %s).", dg[[1]]$value, .fp(dg[[1]]$p)) else "",
      if (nrow(ic)) sprintf(" The MR-Egger intercept was %s (P = %s)%s.", .fnum(ic$estimate, 4), .fp(ic$p), if (ic$p < 0.05) ", suggesting directional pleiotropy" else ", giving no evidence of directional pleiotropy") else "", same, nrow(ce)))
    ts <- cx$ties_sens   # Round 28: the tie sensitivity analysis, when it was run
    if (!is.null(ts)) { te <- .est_df(ts$pay); te <- te[te$estimand == "causal_effect", , drop = FALSE]; ti <- te[te$method == "ivw_mre", , drop = FALSE]
      if (nrow(ti) && any(ce$method == "ivw_mre")) txt$results <- c(txt$results, sprintf("Sensitivity analysis for tied P values: with the variant of larger |beta/SE| kept among variants with equal P (%s left out, %s added; %d instruments used), the IVW estimate was %s, compared with %s in the primary analysis.",
        .n_inst(ts$n_out %||% 0, "primary instruments"), .n_inst(ts$n_in %||% 0, "others", "other"), as.integer(ti$n_ivs[1]), .ci_phrase(ti[1, ], ti$scale[1]), .ci_phrase(ce[ce$method == "ivw_mre", , drop = FALSE][1, ], sc))) }
    es <- cx$eaf_sens   # Round 34: the analysis without the variants flagged by the allele-frequency check
    if (!is.null(es)) { ee <- .est_df(es$pay); ee <- ee[ee$estimand == "causal_effect" & ee$method == "ivw_mre", , drop = FALSE]
      if (nrow(ee) && any(ce$method == "ivw_mre")) txt$results <- c(txt$results, sprintf("Sensitivity analysis without the %d variant%s whose effect-allele frequency differed from the reference panel by more than %s (%s; %d instruments used): the IVW estimate was %s, compared with %s in the primary analysis.",
        length(es$remove), if (length(es$remove) == 1L) "" else "s", format(EAF_FLIP_THRESHOLD), paste(es$remove, collapse = ", "), as.integer(ee$n_ivs[1]), .ci_phrase(ee[1, ], ee$scale[1]), .ci_phrase(ce[ce$method == "ivw_mre", , drop = FALSE][1, ], sc))) }
    ol <- cx$outl   # Round 42: outlier module (FR-103)
    if (!is.null(ol)) txt$results <- c(txt$results, outlier_results_text(ol, ce, sc))
    sf <- steiger_facts(cx)
    if (!is.null(sf)) txt$results <- c(txt$results, sprintf("Steiger directionality test: the variants explained %s of the variance in %s and %s of the variance in %s (R2), %s (P = %s); %d of %d variants explained more variance in the outcome than in the exposure.",
      .fnum(sf$r2x, 3), cx$p$exName, .fnum(sf$r2y, 3), cx$p$outName, if (sf$ok) "so the exposure explains more, as expected if it acts on the outcome" else "so the outcome explains more, which does not match the assumed direction", .fp(sf$p), as.integer(sf$wrong), as.integer(sf$n)))
    txt$results <- c(txt$results, method_result_text(cx))
    for (dm in intersect(c("presso", "radial", "radegger"), cx$ce$method)) {
      dg <- Filter(function(d) d$method == dm, cx$pay$diagnostics); od <- Filter(function(d) grepl("_outliers$", d$name), dg)
      ids <- if (length(od) && !is.null(od[[1]]$note) && !is.na(od[[1]]$note) && nzchar(od[[1]]$note)) strsplit(od[[1]]$note, ";")[[1]] else character(0)
      gl <- Filter(function(d) d$name == "presso_global_test", dg); ql <- Filter(function(d) d$name == "radial_q", dg)
      txt$results <- c(txt$results, sprintf("%s flagged %s of %d variants%s%s. No variant was removed from the primary analysis.", .short[[dm]], if (length(ids)) length(ids) else "none", nrow(cx$kept), if (length(ids)) paste0(" (", paste(ids, collapse = ", "), ")") else "",
        if (length(gl)) sprintf("; the global test P was %s%s", if (!is.null(gl[[1]]$note) && !is.na(gl[[1]]$note)) "below " else "", .fp(gl[[1]]$p)) else if (length(ql)) sprintf("; the radial Q statistic was %.1f (P = %s)", ql[[1]]$value, .fp(ql[[1]]$p)) else ""))
    }
  }
  mm <- intersect(c("wm", "smed", "mode", "smode"), used)   # Round 31: name only the robust estimators that were run
  lim <- c(paste0("The three core MR assumptions (relevance, independence from confounders, and no effect of the variants on the outcome except through the exposure) cannot be proven from summary data.",
                  if (any(mm %in% c("wm", "smed")) && any(mm %in% c("mode", "smode"))) " The median and mode estimators relax the third assumption only under their own conditions (a majority or plurality of valid instruments)."
                  else if (length(mm) && all(mm %in% c("wm", "smed"))) " The median estimators relax the third assumption only if a majority of the instruments are valid."
                  else if (length(mm)) " The mode estimators relax the third assumption only if a plurality of the instruments are valid." else ""))
  if (anysim) lim <- c(lim, "SIMULATED DATA: at least one pair uses simulated demonstration traits. These results illustrate the software and must not be presented as scientific findings.")
  anyog <- any(vapply(ctxs, function(cx) !cx$simulated && !is.null(cx$src_x$og), TRUE)); anyfile <- any(vapply(ctxs, function(cx) !cx$simulated && ((is.null(cx$src_x$og) && is.null(cx$src_x$vcf)) || (is.null(cx$src_y$og) && is.null(cx$src_y$vcf))), TRUE))   # Round 100: OpenGWAS complete GWAS are not analyst-declared files
  xfile <- any(vapply(ctxs, function(cx) !cx$simulated && is.null(cx$src_x$og) && is.null(cx$ledger$clump), TRUE))   # Round 100: not when MR Studio clumped the file itself   # Round 99: the clumping sentence only when an exposure is an imported file
  if (anyfile) lim <- c(lim, paste0(if (xfile) "For imported exposure files, clumping was not performed by MR Studio: the independence of the instruments rests on the analyst's statement at import. " else "",
    "Ancestry and genome build of imported files are as declared by the analyst and were not verified; exposure and outcome must come from the same genome build when variants are matched by position."))
  fw <- vapply(ctxs, function(cx) identical(cx$ledger$strand, "forward"), TRUE); fog <- fw & vapply(ctxs, function(cx) .ogw_both(cx$src_x, cx$src_y), TRUE)   # Round 26; Round 100: OpenGWAS files apart
  if (any(fog)) lim <- c(lim, sprintf("In %s, palindromic variants (A/T or C/G) were kept and aligned by their alleles alone because both datasets are OpenGWAS files, which OpenGWAS harmonises to the forward strand of the reference genome. A sensitivity analysis without them is still advisable.", if (sum(fog) == 1L) "one pair" else sprintf("%d pairs", sum(fog))))
  nfw <- sum(fw & !fog)
  if (nfw) lim <- c(lim, sprintf("In %s, palindromic variants (A/T or C/G) were kept because the analyst declared that both datasets report alleles on the forward strand. This was not verified: if either dataset uses another strand convention, these variants enter with the wrong sign. A sensitivity analysis without them is advisable.", if (nfw == 1L) "one pair" else sprintf("%d pairs", nfw)))
  if (anyog) lim <- c(lim, "Instruments from OpenGWAS were clumped against 1000 Genomes reference panels (on its server or, where stated, on the analyst's computer with the same panel), so linkage disequilibrium was estimated in a reference sample and not in the study population. Proxy variants, where used, stand in for missing instruments and are listed in the harmonised data table.")
  { nf <- vapply(ctxs, function(cx) length(unique(unlist(lapply(cx$ledger$eaf_check %||% list(), function(e) if (isTRUE(e$run)) unlist(e$flagged_ids) else NULL)))), 0L)   # Round 33, distinct variants
    done <- all(vapply(ctxs[nf > 0], function(cx) !is.null(cx$eaf_sens), TRUE))   # Round 34: every affected pair has the analysis without them
    if (sum(nf)) lim <- c(lim, sprintf("In %s, %d variant%s had an effect-allele frequency that differed from the 1000 Genomes reference panel by more than %s, a possible allele or strand error. %s kept in the primary analysis; %s.",
      if (sum(nf > 0) == 1L) "one pair" else sprintf("%d pairs", sum(nf > 0)), sum(nf), if (sum(nf) == 1L) "" else "s", format(EAF_FLIP_THRESHOLD), if (sum(nf) == 1L) "It was" else "They were",
      if (done) sprintf("the Results report a sensitivity analysis without %s", if (sum(nf) == 1L) "it" else "them") else sprintf("a sensitivity analysis without %s is advisable", if (sum(nf) == 1L) "it" else "them"))) }
  if (any(vapply(ctxs, function(cx) !is.null(cx$outl), TRUE))) lim <- c(lim, "Outlier removal was a data-driven sensitivity analysis: the variants removed were chosen after seeing the data, so heterogeneity and pleiotropy P values after removal are optimistic and descriptive only. Heterogeneity can reflect real differences in the effects of valid variants and not only pleiotropy, and removing variants that act through the exposure pathway (vertical pleiotropy) can shift the estimate. The primary analysis used all instruments [Burgess 2023].")
  lim <- c(lim, overlap_limitations(ctxs),
           steiger_limitations(ctxs), binary_exposure_limitations(ctxs),
           { nfew <- sum(vapply(ctxs, function(x) "ivw_mre" %in% x$methods && !is.null(x$kept) && nrow(x$kept) < 4, TRUE))   # Round 30, rule R21b
             if (nfew) sprintf("In %s, fewer than 4 variants were used, so the heterogeneity between the variants could not be estimated reliably and the IVW confidence interval and P value are uncertain [Burgess 2023].", if (nfew == 1) "one pair" else paste(nfew, "pairs")) else NULL },
           heterogeneity_limitation(ctxs),
           if (any(vapply(ctxs, function(x) x$n_weak > 0, TRUE))) "Some variants have an F statistic below 10, which can bias estimates toward the null in two-sample MR." else NULL,
           if (any(c("presso", "radial", "radegger") %in% used)) "Variants flagged by MR-PRESSO or radial MR were reported but not removed from the primary analysis. Detectors flag variants by a statistical rule and cannot tell horizontal pleiotropy from valid heterogeneity; any analysis without flagged variants is a post-hoc sensitivity analysis and must be reported as one [Burgess 2023]." else "Outlier detection (MR-PRESSO, radial MR) was not run in this study. Pleiotropy-robust estimators were compared with the IVW estimate; agreement in direction does not prove the assumptions hold.",
           if (length(ctxs) > 1L) "No correction for multiple testing was applied across pairs." else NULL)   # Round 99: not with a single pair
  txt$limitations <- lim
  src <- unique(unlist(lapply(ctxs, function(cx) if (!cx$simulated) vapply(list(cx$src_x, cx$src_y), function(s) if (!is.null(s$vcf)) paste0("OpenGWAS:", s$vcf$id, gwas_cite(s$vcf, " "), " (complete summary statistics)") else if (!is.null(s$og)) paste0(s$file, gwas_cite(s$og, " ")) else paste0(s$file %||% NA_character_, gwas_cite(s, " ")), "")))); txt$src <- src   # Round 100: OpenGWAS complete GWAS by ID
  txt$availability <- c(if (length(src)) .avail_sentence(src) else AVAIL_SIM, AVAIL_SCRIPTS)
  fi <- 0L
  for (cx in ctxs) { ce <- cx$ce; if (!nrow(ce)) next; fi <- fi + 1L
    txt$legends <- c(txt$legends, sprintf("Figure %d. Forest plot of the Mendelian randomization estimates of %s on %s. Points are estimates and bars are %s confidence intervals for each method (%s). %d instruments were used.", 2L * fi - 1L, cx$p$exName, cx$p$outName, .pct(ce$ci_level[1]), paste(.short[ce$method], collapse = ", "), ce$n_ivs[1]),
                    sprintf("Figure %d. Scatter plot of variant effects. Each point is one variant: its effect on %s (x) and on %s (y), with %s intervals for the outcome effect. Variants are oriented so that the exposure effect is positive. Lines show the slope of each method; the MR-Egger line includes its intercept.", 2L * fi, cx$p$exName, cx$p$outName, "95%")) }
  # Round 39: leave-one-out results sentence per pair and the supplementary figure legends.
  sn <- 0L; txt$supp <- vector("list", length(ctxs))
  # Round 96: a radial plot for each pair that ran radial MR without an outlier analysis (which draws its own); figures are numbered in order
  .rad_ran <- function(cx) any(c("radial", "radegger") %in% cx$ce$method) && is.null(cx$outl)
  for (i in seq_along(ctxs)) { cx <- ctxs[[i]]; lf <- loo_facts(cx); iv <- cx$ce[cx$ce$method == "ivw_mre", , drop = FALSE]; if (is.null(lf) || !nrow(iv)) next
    s1 <- sn + 1L; s2 <- sn + 2L; s3 <- if (.rad_ran(cx)) sn + 3L else NA_integer_; sn <- if (is.na(s3)) s2 else s3; full <- iv$estimate[1]; ch <- lf$loo_b - full; j <- which.max(abs(ch)); flip <- lf$id[sign(lf$loo_b) != sign(full)]
    lost <- if (iv$p[1] < 0.05) lf$id[lf$loo_p >= 0.05] else character(0); nm <- function(x) if (length(x) <= 5) paste(x, collapse = ", ") else sprintf("%s and %d other variants", paste(x[1:5], collapse = ", "), length(x) - 5L)
    txt$results <- c(txt$results, sprintf("Leave-one-out analysis (%s -> %s): with each of the %d variants left out in turn, the IVW estimate ranged from %s to %s (all variants: %s); the largest change, %s%% of the estimate, came from leaving out %s. %s%s The single-variant estimates are shown in a funnel plot (Supplementary Figures S%d and S%d).",
      cx$p$exName, cx$p$outName, nrow(lf), .fnum(min(lf$loo_b)), .fnum(max(lf$loo_b)), .fnum(full), .fnum(100 * abs(ch[j]) / abs(full), 1), lf$id[j],
      if (length(flip)) sprintf("Leaving out %s reversed the direction of the estimate.", nm(flip)) else "No single variant reversed the direction of the estimate.",
      if (iv$p[1] < 0.05) (if (length(lost)) sprintf(" Leaving out %s raised its P value to 0.05 or above.", nm(lost)) else " Its P value stayed below 0.05 whichever variant was left out.") else "",
      s1, s2))
    sf2 <- subset_facts(cx$pay)   # Round 43 (FR-105)
    if (!is.null(sf2)) txt$results <- c(txt$results, paste0(sprintf("The strongest variant, %s, carried %s%% of the total IVW weight; without it the estimate was %s.", sf2$dom_id, .fnum(100 * sf2$dom_share, 1), .fnum(sf2$dom_without)),
      if (sf2$lso) sprintf(" Leaving out random subsets (%s, %d repeats with a fixed seed), 95%% of the IVW estimates lay between %s and %s (median %s; range %s to %s) [Burgess 2023].", sf2$lso_note, sf2$lso_n, .fnum(sf2$lso_lo), .fnum(sf2$lso_hi), .fnum(sf2$lso_med), .fnum(sf2$lso_min), .fnum(sf2$lso_max)) else ""))
    ln <- intersect(c("ivw_mre", "egger"), cx$ce$method)
    txt$supp[[i]] <- c(sprintf("Supplementary Figure S%d. Leave-one-out analysis of %s on %s: the IVW estimate with each variant left out in turn, with %s confidence intervals; the bottom row uses all %d variants.%s", s1, cx$p$exName, cx$p$outName, .pct(iv$ci_level[1]), nrow(lf),
                         if (nrow(lf) > LOO_PLOT_MAX) sprintf(" Only the %d variants whose removal changed the estimate most are shown; all %d are in the supplementary table.", LOO_PLOT_SHOWN, nrow(lf)) else ""),
                       sprintf("Supplementary Figure S%d. Funnel plot of the single-variant (Wald ratio) estimates of %s on %s against their precision (1 / standard error). Vertical lines show the %s estimate%s using all variants. Asymmetry can indicate directional pleiotropy. Right: density of the same estimates weighted by their precision, with a tick for each variant; several peaks can indicate groups of variants with different effects.", s2, cx$p$exName, cx$p$outName, paste(.short[ln], collapse = " and "), if (length(ln) > 1) "s" else ""),
                       if (!is.na(s3)) sprintf("Supplementary Figure S%d. Radial plot of %s on %s: each point is one variant, its ratio estimate times the square root of its weight against the square root of its weight; the line is the IVW estimate. Points far from the line contribute most to Cochran's Q [Bowden 2018].", s3, cx$p$exName, cx$p$outName) else NULL)
  }
  bd <- bidir_pairs(ctxs)   # Round 45: bidirectional MR
  if (length(bd)) { txt$methods <- c(txt$methods, BIDIR_METHODS_TEXT); txt$limitations <- c(txt$limitations, BIDIR_LIMITS_TEXT)
    for (ij in bd) txt$results <- c(txt$results, bidir_results_text(ctxs[[ij[1]]], ctxs[[ij[2]]])) }
  txt$supp_out <- vector("list", length(ctxs)); so <- sn   # Round 44: outlier figures follow the leave-one-out ones
  for (i in seq_along(ctxs)) { ol <- ctxs[[i]]$outl; if (is.null(ol)) next
    ag <- outlier_agreement(ol); txt$supp_out[[i]] <- outlier_legends(ol, ctxs[[i]], so + 1L, ag, loo = outlier_has_loo(ctxs[[i]])); so <- so + length(txt$supp_out[[i]]) }
  txt
}

# ---------- the package ----------
`%||%` <- function(a, b) if (is.null(a)) b else a
package_dir <- function(h) { root <- file.path(h$dir, "05_publication"); dir.create(root, showWarnings = FALSE); n <- length(list.dirs(root, recursive = FALSE)) + 1L; d <- file.path(root, sprintf("package_%03d", n)); dir.create(d); d }

#' Round 75b: the key under the summary table, naming each number: the estimate, its confidence level, the P value, the instruments.
.est_key <- function(ctxs) {
  lev <- unique(unlist(lapply(ctxs, function(x) x$ce$ci_level))); lev <- lev[is.finite(lev)]
  lv <- if (length(lev) == 1L) sprintf("%g%%", 100 * lev) else "the stated"
  sprintf("Estimate: the effect of the exposure on the outcome per unit of the exposure (β; for a binary outcome the log odds ratio, whose exponential is the odds ratio); for multivariable MR, the direct effect given the other exposures. Lower and upper confidence limits: the %s confidence interval. P value: two-sided. Instruments: the variants the method used.", lv)
}
.and_join <- function(x) if (length(x) < 3L) paste(x, collapse = " and ") else paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])   # Round 101: "A, B and C"
build_package <- function(h, items, width_mm = 180, pub = NULL, pairs = NULL) {
  pairs <- pairs %||% (ui_state_get(h$con, "pairs") %||% list()); pairs <- Filter(function(p) !is.null(p$folder) && !is.null(pair_results(h$dir, p$folder)), pairs)
  if (!length(pairs)) stop("There are no saved results to put in a package. Run the analysis first.")
  pairs <- c(Filter(function(p) !isTRUE(p$mv), pairs), Filter(function(p) isTRUE(p$mv), pairs))   # Round 50: univariable pairs first, so text and figures keep their numbering
  ctxs <- lapply(pairs, function(p) if (isTRUE(p$mv)) mv_context(h, p) else pair_context(h, p)); pk <- package_dir(h); made <- character(); note <- character()
  add <- function(f) made <<- c(made, substring(f, nchar(h$dir) + 2L))
  txt <- report_text(ctxs, basename(h$dir))
  pwr <- tryCatch(ui_state_get(h$con, "power"), error = function(e) NULL); if (!is.null(pwr$text)) txt$methods <- c(txt$methods, pwr$text)   # Round 109: the power calculation
  med <- if (any(vapply(pairs, function(p) !is.null(p$med), TRUE))) tryCatch(mediation_run(h, pairs), error = function(e) NULL) else NULL   # Round 57: mediation
  if (!is.null(med)) { txt$methods <- c(txt$methods, MEDIATION_METHODS_TEXT); txt$results <- c(txt$results, mediation_results_text(med, attr(med, "spec"))); txt$limitations <- c(setdiff(txt$limitations, "No correction for multiple testing was applied across pairs."), MEDIATION_LIMITS_TEXT) }   # Round 98: the steps of one mediation analysis are not separate tests
  bcm <- if (exists("bidir_cml_run") && any(vapply(pairs, function(p) isTRUE(p$rev), TRUE)) && requireNamespace("BiDirectCausal", quietly = TRUE)) tryCatch(Filter(is.data.frame, bidir_cml_run(h, pairs, reuse = TRUE)), error = function(e) list()) else list()   # Round 63; Round 95: the saved result is reused
  if (any(vapply(pairs, function(p) isTRUE(p$rev), TRUE))) txt$methods <- c(txt$methods, "The study was bidirectional: each pair was also analysed in the reverse direction as a separate two-sample MR, with the instruments of the outcome trait and their effects looked up in the exposure trait.")   # Round 95
  if (length(bcm)) { txt$methods <- c(txt$methods, BIDIR_CML_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(bcm, bidir_cml_text))) }
  lds <- if (exists("ldsc_pair_read")) Filter(Negate(is.null), lapply(Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs), function(p) tryCatch(ldsc_pair_read(h, p), error = function(e) NULL))) else list()   # Round 64
  if (length(lds)) { txt$methods <- c(txt$methods, LDSC_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(lds, ldsc_text))) }
  lhc <- if (exists("lhcmr_pair_read")) Filter(Negate(is.null), lapply(Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs), function(p) tryCatch(lhcmr_pair_read(h, p), error = function(e) NULL))) else list()   # Round 65
  if (length(lhc)) { txt$methods <- c(txt$methods, LHCMR_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(lhc, lhcmr_text))) }
  mlp <- if (exists("mrlap_pair_read")) Filter(Negate(is.null), lapply(Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs), function(p) tryCatch(mrlap_pair_read(h, p), error = function(e) NULL))) else list()   # Round 68
  if (length(mlp)) { txt$methods <- c(txt$methods, MRLAP_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(mlp, mrlap_text))) }
  aps <- if (exists("mrapss_pair_read")) Filter(Negate(is.null), lapply(Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs), function(p) tryCatch(mrapss_pair_read(h, p), error = function(e) NULL))) else list()   # Round 69
  if (length(aps)) { txt$methods <- c(txt$methods, MRAPSS_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(aps, mrapss_text))) }
  cau <- if (exists("cause_pair_read")) Filter(Negate(is.null), lapply(Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs), function(p) tryCatch(cause_pair_read(h, p), error = function(e) NULL))) else list()   # Round 70
  if (length(cau)) { txt$methods <- c(txt$methods, CAUSE_METHODS_TEXT); txt$results <- c(txt$results, unlist(lapply(cau, cause_text))); txt$limitations <- c(txt$limitations, unlist(lapply(cau, cause_ties_text))) }
  txt$methods <- gsub("of which 1 were absent", "of which 1 was absent", txt$methods, fixed = TRUE)   # Round 100: grammar
  adj <- c(if (length(mlp)) "MRlap", if (length(aps)) "MR-APSS", if (length(cau)) "CAUSE")   # Round 100: the overlap limitation names the methods that allow for it
  if (length(adj)) { txt$limitations <- sub("MR Studio did not correct for this. ", sprintf("The IVW and the other two-sample estimates were not corrected for it; %s, reported in the Results, allow%s for sample overlap. ", .and_join(adj), if (length(adj) == 1L) "s" else ""), txt$limitations, fixed = TRUE)
    txt$limitations <- sub(" Consider an outcome from a different cohort or a method designed for overlapping samples.", " Compare their estimates with the IVW estimate.", txt$limitations, fixed = TRUE) }   # Round 87: tied lead variants
  # Round 81: the LD reference of every genome-wide method that gave a result; non-European panels and the analyst's accepted choices for Mixed data
  if (exists("ldref_read")) {
    one <- Filter(function(p) !isTRUE(p$mv) && !is.null(p$folder), pairs); recs <- list()
    for (m in c("ldsc", "mrlap", "mrapss", "cause")) for (p in one) { rd <- get0(paste0(m, "_pair_read")); if (is.null(rd) || is.null(tryCatch(rd(h, p), error = function(e) NULL))) next
      r <- tryCatch(ldref_read(h, p, m), error = function(e) NULL); if (!is.null(r)) recs[[length(recs) + 1L]] <- c(r, list(labels = c(p$exName, p$outName))) }
    for (m in c("cause", "mrapss")) { pp <- unique(unlist(lapply(Filter(function(r) identical(r$method, m), recs), `[[`, "pop")))
      if (length(pp) && !identical(pp, "EUR")) { hit <- which(txt$methods == get(if (m == "cause") "CAUSE_METHODS_TEXT" else "MRAPSS_METHODS_TEXT"))
        if (length(hit)) txt$methods[hit] <- sub("the 1000 Genomes European panel", sprintf("the 1000 Genomes %s panel%s", paste(LDREF_POP_NAMES[pp], collapse = " or "), if (length(pp) > 1) ", as recorded for each pair in its LD reference file (04_across_pairs)" else ""), txt$methods[hit], fixed = TRUE) } }
    txt$limitations <- c(txt$limitations, ldref_limitations(recs))   # Round 101: one sentence per pair and panel choice
  }
  txt <- .pfix(txt)   # Round 91: "P = < 1e-300" reads "P < 1e-300"
  used <- unique(c("hemani2018", "burgess2023", "skrivankova2021", if (length(bcm)) "xue2022", if (length(lds)) c("bulik2015a", "bulik2015b", "grotzinger2019"), if (length(lhc)) "darrous2021", if (length(mlp)) "mounier2023", if (length(aps)) "hu2022", if (length(cau)) "morrison2020", if (!is.null(med)) c("carter2021", "carter2019", "sanderson2019"), unlist(METHOD_REFS[unlist(lapply(ctxs, function(x) x$methods))]), if (steiger_any_binary(ctxs)) "lee2012", if (length(.ctx_bin_x(ctxs))) "burgess2018", if (any(vapply(ctxs, function(cx) !is.null(cx$ledger$clump), TRUE))) c("chang2015", "kgp2015")))
  figs <- list(); wmm <- as.numeric(width_mm); wmm <- if (is.finite(wmm) && wmm >= 60) wmm else 180
  if (any(c("figs", "report", "text") %in% items)) {
    .say("Drawing the figures: forest and scatter plots of every pair, at the column width you chose", "Figures are written as PDF, SVG and PNG, so they can go straight into a manuscript.")
    dir.create(file.path(pk, "figures"), showWarnings = FALSE)
    for (ns in FIG_NAMESPACES) suppressWarnings(requireNamespace(ns, quietly = TRUE))   # Round 109: the same loaded packages for every figure, whatever ran before
    fig_warmup()
    FIG_REC$on <- TRUE; FIG_REC$items <- list()   # Round 109: record the figures for the figures script of each pair
    for (i in seq_along(ctxs)) { cx <- ctxs[[i]]; ce <- cx$ce; if (!nrow(ce)) next; tag <- sprintf("P%d", i)
      est <- if (ce$scale[1] == "log_or") transform(ce, label = .short[method], col = .MCOL[method], estimate = exp(estimate), ci_lower = exp(ci_lower), ci_upper = exp(ci_upper), unit = "odds ratio (log scale)", ratio = TRUE) else transform(ce, label = .short[method], col = .MCOL[method], unit = "beta: effect on outcome per unit of exposure (often 1 SD)", ratio = FALSE)   # Round 109: OR for a binary outcome
      if (isTRUE(cx$mv)) {   # Round 50: one forest of direct effects, grouped by exposure; no scatter or variant-level plots
        est$label <- paste0(ce$exposure_label, ": ", .short_mv(ce$method)); est$col <- "#2b3fc7"; est <- est[order(match(est$exposure_label, unique(est$exposure_label))), , drop = FALSE]
        fh <- max(40, 22 + 8 * nrow(est)); ff <- .fig_write(file.path(pk, "figures", paste0(tag, "_mvmr_forest")), wmm, fh, function() .draw_forest(est))
        figs[[i]] <- list(forest = ff, fh = fh); for (f in ff) add(f)
        xf <- tryCatch(mv_extra_figures(cx, pk, tag, wmm, (ce$ci_level %||% 0.95)[1]), error = function(e) list()); figs[[i]]$xtra <- xf; for (x in xf) for (f in x$files) add(f)   # Round 110: item 12
        next }
      fh <- max(40, 22 + 8 * nrow(est)); ff <- .fig_write(file.path(pk, "figures", paste0(tag, "_forest")), wmm, fh, function() .draw_forest(est))
      ceS <- ce[ce$method != "radegger", , drop = FALSE]   # the radial Egger intercept is on a weighted scale, so it cannot be drawn as a line on this plot
      ln <- lapply(seq_len(nrow(ceS)), function(j) { m <- ceS$method[j]; ic <- cx$est$estimate[cx$est$estimand == "pleiotropy_intercept" & cx$est$method == "egger"]; list(name = .short[[m]], col = .MCOL[[m]], lty = if (m == "ivw_mre") 1 else if (m == "egger") 2 else 3, a = if (m == "egger" && length(ic)) ic[1] else 0, b = ceS$estimate[j]) })
      sf <- .fig_write(file.path(pk, "figures", paste0(tag, "_scatter")), wmm, round(wmm * 0.75), function() .draw_scatter(cx$kept, ln, paste("Effect on", cx$p$exName), paste("Effect on", cx$p$outName)))
      figs[[i]] <- list(forest = ff, scatter = sf, fh = fh, sh = round(wmm * 0.75)); for (f in c(ff, sf)) add(f)
      xf <- tryCatch(extra_figures(cx, pk, tag, wmm, (ce$ci_level %||% 0.95)[1]), error = function(e) list()); figs[[i]]$xtra <- xf; for (x in xf) for (f in x$files) add(f)   # Round 109
      lf <- loo_facts(cx); iv <- ce[ce$method == "ivw_mre", , drop = FALSE]   # Round 39: supplementary leave-one-out and funnel plots
      if (!is.null(lf) && nrow(iv)) {
        z <- stats::qnorm(1 - (1 - iv$ci_level[1]) / 2); full <- list(b = iv$estimate[1], lo = iv$ci_lower[1], hi = iv$ci_upper[1])
        lh <- max(50, round(18 + 4.2 * (nrow(.loo_rows(lf, full$b)) + 1L)))
        lo <- .fig_write(file.path(pk, "figures", paste0(tag, "_leave_one_out")), wmm, lh, function() .draw_loo(lf, full, z, est$unit[1]))
        fl <- lapply(intersect(c("ivw_mre", "egger"), ce$method), function(m) list(name = .short[[m]], col = .MCOL[[m]], lty = if (m == "ivw_mre") 1 else 2, b = ce$estimate[ce$method == m][1]))
        fu <- .fig_write(file.path(pk, "figures", paste0(tag, "_funnel")), wmm, round(wmm * 0.75), function() .draw_funnel(lf, fl))
        figs[[i]]$loo <- lo; figs[[i]]$funnel <- fu; figs[[i]]$lh <- lh; for (f in c(lo, fu)) add(f)
        if (any(c("radial", "radegger") %in% ce$method) && is.null(cx$outl) && exists("draw_outlier_radial")) { rd <- .fig_write(file.path(pk, "figures", paste0(tag, "_radial")), wmm, round(wmm * 0.75), function() draw_outlier_radial(cx$pay, NULL, est$unit[1])); figs[[i]]$rad <- rd; for (f in rd) add(f) }   # Round 96
      }
      ol <- cx$outl   # Round 44: outlier comparison and agreement figures (FR-98)
      if (!is.null(ol) && exists("draw_outlier_comparison")) {
        oh <- max(45, round(20 + 7 * length(.outl_rows(ol))))
        oc <- .fig_write(file.path(pk, "figures", paste0(tag, "_outlier_comparison")), wmm, oh, function() draw_outlier_comparison(ol, est$unit[1]))
        ag <- outlier_agreement(ol); ah <- if (is.null(ag)) 0 else max(45, round(30 + 4.2 * nrow(ag$m)))
        oa <- if (is.null(ag)) character() else .fig_write(file.path(pk, "figures", paste0(tag, "_outlier_agreement")), min(wmm, 120), ah, function() draw_outlier_agreement(ag))
        rr <- .fig_write(file.path(pk, "figures", paste0(tag, "_outlier_radial")), wmm, round(wmm * 0.75), function() draw_outlier_radial(cx$pay, ol, est$unit[1]))   # Round 47
        lb <- loo_facts(cx); la <- if (outlier_has_loo(cx)) loo_facts(list(pay = ol$chosen_pay)) else NULL
        ob <- if (!is.null(la)) .fig_write(file.path(pk, "figures", paste0(tag, "_outlier_leave_one_out")), wmm, 55, function() draw_outlier_loo(lb, la, ce$estimate[ce$method == "ivw_mre"][1], ol$chosen_ivw$estimate[1], est$unit[1])) else character()
        figs[[i]]$oc <- oc; figs[[i]]$oa <- oa; figs[[i]]$oh <- oh; figs[[i]]$ah <- ah; figs[[i]]$or <- rr; figs[[i]]$ol <- ob; for (f in c(oc, oa, rr, ob)) add(f)
      } }
    FIG_REC$on <- FALSE; for (i in seq_along(ctxs)) if (!is.null(figs[[i]])) tryCatch(write_figures_script(h, ctxs[[i]], i, wmm, FIG_REC$items, vapply(figs[[i]]$xtra %||% list(), `[[`, "", "key")), error = function(e) note <<- c(note, paste("The figures script could not be written:", conditionMessage(e))))   # Round 109
    if (!is.null(med)) { md <- mediation_forest_df(med); if (nrow(md)) {   # Round 98: the mediation estimates as a figure
      lg <- txt$legends; length(lg) <- 2L * length(ctxs); lg[is.na(lg)] <- ""; k <- length(ctxs) + 1L; fh <- max(40, 22 + 8 * nrow(md))
      mf <- .fig_write(file.path(pk, "figures", "mediation_forest"), wmm, fh, function() .draw_forest(md)); figs[[k]] <- list(forest = mf, fh = fh); for (f in mf) add(f)
      txt$legends <- c(lg, mediation_legend(attr(med, "spec"), sum(grepl("^Figure [0-9]+\\.", lg)) + 1L), "") } }
    if (!length(Filter(Negate(is.null), figs))) note <- c(note, "No figures could be drawn on this computer.")
    if (!isTRUE(unname(capabilities("cairo")))) note <- c(note, "SVG and PNG figures need the cairo graphics library, which this R does not have. The PDF figures were written.")
  }
  if ("scripts" %in% items) { dst <- file.path(pk, "scripts_and_tables"); dir.create(dst); for (d in list.dirs(file.path(h$dir, "03_pairs"), recursive = FALSE)) { file.copy(d, dst, recursive = TRUE) }; if (dir.exists(file.path(h$dir, "04_across_pairs"))) file.copy(file.path(h$dir, "04_across_pairs"), dst, recursive = TRUE); .rm_empty_dirs(dst); add(dst) }   # Round 96: every pair folder, as in the capsule (outlier reruns too)
  if ("tables" %in% items) {
    .say("Building the supplement workbook: estimates, settings and the instruments of every pair", "A supplement with every instrument and setting lets reviewers check the analysis.")
    est <- do.call(rbind, lapply(seq_along(ctxs), function(i) { e <- ctxs[[i]]$est; data.frame(pair = i, exposure = ctxs[[i]]$p$exName, outcome = ctxs[[i]]$p$outName, e, check.names = FALSE, stringsAsFactors = FALSE) }))
    if (nrow(est) && "scale" %in% names(est)) { r <- est$scale %in% "log_or" & est$estimand %in% c("causal_effect", "direct_effect", "indirect_effect", "total_effect", "reverse_effect")   # Round 109: binary outcome as OR (95% CI) beside the log OR
      est$or <- ifelse(r, exp(est$estimate), NA_real_); est$or_ci_lower <- ifelse(r, exp(est$ci_lower), NA_real_); est$or_ci_upper <- ifelse(r, exp(est$ci_upper), NA_real_) }
    dgs <- do.call(rbind, lapply(seq_along(ctxs), function(i) { d <- Filter(function(x) !identical(x$method, "loo"), ctxs[[i]]$pay$diagnostics %||% list()); if (!length(d)) return(NULL)   # Round 103: every method's diagnostics in one sheet
      f <- function(k, num = TRUE) vapply(d, function(x) { v <- x[[k]]; if (is.null(v) || !length(v)) NA_character_ else as.character(v[[1]]) }, ""); nu <- function(x) suppressWarnings(as.numeric(x))
      data.frame(pair = i, method = vapply(d, function(x) if (x$method %in% names(.short)) .short[[x$method]] else x$method, ""), method_id = f("method"), name = f("name"), value = nu(f("value")), df = nu(f("df")), p = nu(f("p")), note = f("note"), stringsAsFactors = FALSE) }))
    mtab <- tryCatch(do.call(rbind, lapply(seq_along(ctxs), function(i) { cx <- ctxs[[i]]; if (isTRUE(cx$mv) || !any(c("mode", "smode", "smode_nome", "wmode_nome", "mbe") %in% (cx$ce$method %||% character()))) return(NULL)   # Round 109: item 9
      m <- mode_table(cx$kept, (cx$ce$ci_level %||% 0.95)[1]); if (is.null(m)) NULL else data.frame(pair = i, exposure = cx$p$exName, outcome = cx$p$outName, m, stringsAsFactors = FALSE) })), error = function(e) NULL)
    sheets <- list(Estimates = est, Settings = do.call(rbind, lapply(seq_along(ctxs), function(i) { s <- ctxs[[i]]$settings; data.frame(pair = i, setting = names(s), value = vapply(s, function(v) paste(v, collapse = ","), ""), stringsAsFactors = FALSE) })))
    if (!is.null(dgs) && nrow(dgs)) sheets <- c(sheets[1], list(Diagnostics = dgs), sheets[-1])
    for (i in seq_along(ctxs)) sheets[[sprintf("Pair %d variants", i)]] <- if (isTRUE(ctxs[[i]]$mv)) { f <- Sys.glob(file.path(h$dir, "03_pairs", ctxs[[i]]$p$folder, "01_data", "analysis_data_*.csv")); if (length(f)) utils::read.csv(f[1], stringsAsFactors = FALSE) else data.frame() } else transform(ctxs[[i]]$vars, ratio = by / bx, F = (bx / sx)^2)
    for (i in seq_along(ctxs)) { lf <- loo_facts(ctxs[[i]]); if (!is.null(lf)) sheets[[sprintf("Pair %d leave-one-out", i)]] <- stats::setNames(lf, c("variant_left_out", "ivw_estimate_without", "ivw_se_without", "ivw_p_without", "wald_ratio", "wald_ratio_se", "wald_ratio_p")) }
    if (!is.null(med)) { ms <- as.data.frame(med); attr(ms, "spec") <- NULL; sheets[["Mediation"]] <- ms }   # Round 98
    for (m in c(ldsc = "LDSC", mrlap = "MRlap", mrapss = "MR-APSS", cause = "CAUSE", lhcmr = "LHC-MR") |> names()) {   # Round 100: genome-wide methods in the supplement
      fs <- sort(Sys.glob(file.path(h$dir, "04_across_pairs", m, "*_results.csv"))); if (!length(fs)) next
      tb <- tryCatch(do.call(rbind, lapply(fs, function(f) { d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE); data.frame(pair = sub("_results\\.csv$", "", basename(f)), d, check.names = FALSE, stringsAsFactors = FALSE) })), error = function(e) NULL)
      if (is.data.frame(tb) && nrow(tb)) sheets[[c(ldsc = "LDSC", mrlap = "MRlap", mrapss = "MR-APSS", cause = "CAUSE", lhcmr = "LHC-MR")[[m]]]] <- tb }
    if (!is.null(mtab) && nrow(mtab)) sheets$`Mode estimates` <- mtab   # Round 109: mode estimates at phi 0.25, 0.5 and 1, simple, weighted and NOME
    f <- file.path(pk, "supplement_tables.xlsx"); xlsx_write(f, sheets); add(f)
  }
  if ("refs" %in% items) { .say("Collecting the citations of every method and tool named in the text", "Each reference was checked against PubMed when MR Studio was built."); f1 <- file.path(pk, "references.bib"); writeLines(refs_bibtex(used), f1, useBytes = TRUE); f2 <- file.path(pk, "references.ris"); writeLines(refs_ris(used), f2, useBytes = TRUE); add(f1); add(f2) }
  blocks_text <- function(with_figs = FALSE) {
    b <- list(); H <- function(t, l = 1) b[[length(b) + 1L]] <<- list(type = if (l == 1) "h1" else "h2", text = t); P <- function(t) for (x in t) b[[length(b) + 1L]] <<- list(type = "p", text = x)
    if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) b[[1]] <- list(type = "note", text = "SIMULATED DATA: not for publication. At least one pair in this document uses simulated demonstration traits.")
    H("Methods"); P(txt$methods); H("Results"); P(txt$results)
    if (with_figs) for (i in seq_along(figs)) { f <- figs[[i]]; if (is.null(f)) next; png <- grep("forest\\.png$", f$forest, value = TRUE); sp <- grep("scatter\\.png$", f$scatter, value = TRUE)
      if (length(png)) b[[length(b) + 1L]] <- list(type = "image", file = png, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * f$fh / wmm, alt = txt$legends[2L * i - 1L]); P(txt$legends[2L * i - 1L])
      if (length(sp)) { b[[length(b) + 1L]] <- list(type = "image", file = sp, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * f$sh / wmm, alt = txt$legends[2L * i]); P(txt$legends[2L * i]) } }
    if (with_figs && any(vapply(figs, function(f) !is.null(f$loo) || !is.null(f$oc), TRUE))) { H("Supplementary figures", 2)   # Round 39
      for (i in seq_along(figs)) { f <- figs[[i]]; if (is.null(f$loo) || is.null(.nth(txt$supp, i))) next
        lp <- grep("\\.png$", f$loo, value = TRUE); fp <- grep("\\.png$", f$funnel, value = TRUE)
        if (length(lp)) b[[length(b) + 1L]] <- list(type = "image", file = lp, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * f$lh / wmm, alt = txt$supp[[i]][1]); P(txt$supp[[i]][1])
        if (length(fp)) b[[length(b) + 1L]] <- list(type = "image", file = fp, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * 0.75, alt = txt$supp[[i]][2]); P(txt$supp[[i]][2])
        rp <- grep("\\.png$", f$rad, value = TRUE); if (length(rp) && length(txt$supp[[i]]) >= 3L) { b[[length(b) + 1L]] <- list(type = "image", file = rp, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * 0.75, alt = txt$supp[[i]][3]); P(txt$supp[[i]][3]) } }
      for (i in seq_along(figs)) { f <- figs[[i]]; so <- .nth(txt$supp_out, i); if (is.null(f$oc) || is.null(so)) next   # Round 44
        im <- list(comparison = list(f$oc, f$oh / wmm), agreement = list(f$oa, NA), radial = list(f$or, 0.75), loo = list(f$ol, 55 / wmm))
        for (k in names(so)) { pn <- grep("\\.png$", im[[k]][[1]], value = TRUE); if (!length(pn)) next
          w <- if (k == "agreement") min(wmm, 120) else min(wmm, 160); hh <- if (k == "agreement") f$ah else w * im[[k]][[2]]
          b[[length(b) + 1L]] <- list(type = "image", file = pn, width_mm = w, height_mm = hh, alt = so[[k]]); P(so[[k]]) } } }
    H("Limitations"); P(txt$limitations); H("Data availability"); P(txt$availability); H("References"); P(vapply(used, function(i) ref_text(REFS[[i]]), "")); b }
  if ("text" %in% items) { .say("Writing the draft Methods, Results and Limitations with your numbers and settings", "The draft names and cites every method; the limitations include every item the Check screen asked you to review."); f <- file.path(pk, "draft_methods_results_limitations.docx"); docx_write(f, c(list(list(type = "title", text = "Draft text: Methods, Results, Limitations")), list(list(type = "note", text = "Draft generated from the saved results and settings. Check every statement and edit it before use.")), blocks_text(FALSE))); add(f) }
  if ("report" %in% items) {
    .say("Writing the full report (Word and PDF): summary table, methods, results, figures, limitations and references")
    est <- do.call(rbind, lapply(seq_along(ctxs), function(i) { ce <- ctxs[[i]]$ce; if (!nrow(ce)) return(NULL); data.frame(Pair = i, Exposure = if (!is.null(ce$exposure_label)) ce$exposure_label else ctxs[[i]]$p$exName, Outcome = ctxs[[i]]$p$outName, Method = ce$name, `Estimate (β or log OR)` = .fnum(ce$estimate), `Lower confidence limit` = .fnum(ce$ci_lower), `Upper confidence limit` = .fnum(ce$ci_upper), `P value` = .fp_res(ce$p, ce$p_resolution), `Instruments (variants used)` = ce$n_ivs, check.names = FALSE, stringsAsFactors = FALSE) }))
    bl <- c(list(list(type = "title", text = sprintf("Mendelian randomization report: %s", basename(h$dir))), list(type = "note", text = sprintf("Generated by MR Studio %s on %s. A draft for the authors to check.", .mrs_ver_build(), format(Sys.Date())))), list(list(type = "h1", text = "Summary of estimates"), list(type = "table", data = est), list(type = "note", text = .est_key(ctxs))), blocks_text(TRUE))
    f <- file.path(pk, "report.docx"); docx_write(f, bl); add(f)
    fp <- file.path(pk, "report.pdf"); r46 <- tryCatch({ pdf_write(fp, bl); TRUE }, error = function(e) { if (file.exists(fp)) unlink(fp); note <<- c(note, paste("The PDF report could not be written:", conditionMessage(e))); FALSE }); if (r46) add(fp)   # Round 46
    esc <- function(x) .xe(x); svg_k <- 0L; svgs <- function(files) { s <- grep("\\.svg$", files, value = TRUE); if (!length(s)) return("<p><em>Figure not available on this computer.</em></p>"); svg_k <<- svg_k + 1L; p <- sprintf("f%d-", svg_k)
      x <- paste(sub("^<\\?xml[^>]*\\?>\\s*", "", readLines(s[1], warn = FALSE, encoding = "UTF-8")), collapse = "\n")
      gsub("url(#", paste0("url(#", p), gsub("href=\"#", paste0("href=\"#", p), gsub(" id=\"", paste0(" id=\"", p), x, fixed = TRUE), fixed = TRUE), fixed = TRUE) }
    sect <- function(t, ps) sprintf("<h2>%s</h2>%s", esc(t), paste(sprintf("<p>%s</p>", esc(ps)), collapse = "\n"))
    figh <- paste(vapply(seq_along(figs), function(i) { f <- figs[[i]]; if (is.null(f)) "" else if (is.null(f$scatter)) sprintf("<figure>%s<figcaption>%s</figcaption></figure>", svgs(f$forest), esc(txt$legends[2L * i - 1L])) else sprintf("<figure>%s<figcaption>%s</figcaption></figure><figure>%s<figcaption>%s</figcaption></figure>", svgs(f$forest), esc(txt$legends[2L * i - 1L]), svgs(f$scatter), esc(txt$legends[2L * i])) }, ""), collapse = "\n")
    supph <- paste(vapply(seq_along(figs), function(i) { f <- figs[[i]]; if (is.null(f$loo) || is.null(.nth(txt$supp, i))) "" else sprintf("<figure>%s<figcaption>%s</figcaption></figure><figure>%s<figcaption>%s</figcaption></figure>", svgs(f$loo), esc(txt$supp[[i]][1]), svgs(f$funnel), esc(txt$supp[[i]][2])) |> paste0(if (length(f$rad) && length(txt$supp[[i]]) >= 3L) sprintf("<figure>%s<figcaption>%s</figcaption></figure>", svgs(f$rad), esc(txt$supp[[i]][3])) else "") }, ""), collapse = "\n")   # Round 39; Round 96: radial
    supph <- paste0(supph, paste(vapply(seq_along(figs), function(i) { f <- figs[[i]]; so <- .nth(txt$supp_out, i); if (is.null(f$oc) || is.null(so)) "" else   # Round 44
      paste(vapply(names(so), function(k) { fl <- switch(k, comparison = f$oc, agreement = f$oa, radial = f$or, loo = f$ol); if (!length(fl)) "" else sprintf("<figure>%s<figcaption>%s</figcaption></figure>", svgs(fl), esc(so[[k]])) }, ""), collapse = "") }, ""), collapse = ""))
    supph <- paste0(supph, paste(vapply(seq_along(figs), function(i) paste(vapply(figs[[i]]$xtra %||% list(), function(x) sprintf("<figure>%s<figcaption>%s</figcaption></figure>", svgs(x$files), esc(x$caption)), ""), collapse = ""), ""), collapse = ""))   # Round 109
    if (nzchar(supph)) supph <- paste0("<h2>Supplementary figures</h2>", supph)
    tab <- paste0("<table><thead><tr>", paste(sprintf("<th>%s</th>", esc(names(est))), collapse = ""), "</tr></thead><tbody>", paste(vapply(seq_len(nrow(est)), function(i) paste0("<tr>", paste(sprintf("<td>%s</td>", esc(unlist(est[i, ]))), collapse = ""), "</tr>"), ""), collapse = ""), "</tbody></table>")
    sim <- if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) "<div class=\"warn\">SIMULATED DATA: not for publication. At least one pair uses simulated demonstration traits.</div>" else ""
    html <- paste0("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>MR report: ", esc(basename(h$dir)), "</title><style>body{font:16px/1.55 Georgia,'Times New Roman',serif;max-width:46rem;margin:2.5rem auto;padding:0 1.2rem;color:#1b1b1b}h1{font-size:1.9rem}h2{margin-top:2.2rem;font-size:1.3rem;border-bottom:1px solid #ccc;padding-bottom:.2rem}table{border-collapse:collapse;width:100%;font:14px/1.4 system-ui,sans-serif}th,td{border:1px solid #bbb;padding:.35rem .5rem;text-align:left}th{background:#f1f1f1}figure{margin:1.5rem 0}figure svg{max-width:100%;height:auto}figcaption{font-size:.9rem;color:#444}.warn{border:2px solid #b5541c;background:#fff4ec;padding:.7rem 1rem;margin:1rem 0;font:600 15px system-ui,sans-serif}.note{color:#555;font-style:italic}@media print{body{margin:0;max-width:none}h2{break-after:avoid}figure,table{break-inside:avoid}}</style></head><body><h1>Mendelian randomization report: ", esc(basename(h$dir)), "</h1><p class=\"note\">Generated by MR Studio ", esc(.mrs_ver_build()), " on ", format(Sys.Date()), ". A draft for the authors to check. To make a PDF, open this file in a browser and print it.</p>", sim, "<h2>Summary of estimates</h2>", tab, sect("Methods", txt$methods), sect("Results", txt$results), "<h2>Figures</h2>", figh, supph, sect("Limitations", txt$limitations), sect("Data availability", txt$availability), sprintf("<h2>References</h2><ol>%s</ol>", paste(sprintf("<li>%s</li>", esc(vapply(used, function(i) ref_text(REFS[[i]]), ""))), collapse = "")), "</body></html>")
    f <- file.path(pk, "report.html"); writeLines(html, f, useBytes = TRUE); add(f)
  }
  if ("capsule" %in% items) {
    .say("Packing the reproducibility capsule: every pair folder with its scripts, data, results and replay check, and the R version", "Anyone can rerun the whole analysis from the capsule and compare the numbers.")
    fs <- list.files(file.path(h$dir, c("03_pairs", "04_across_pairs")), recursive = TRUE, full.names = TRUE); rel <- substring(fs, nchar(h$dir) + 2L)   # Round 98: with the analyses across pairs
    ent <- setNames(lapply(fs, .file_raw), rel)
    ent[["SESSION.txt"]] <- .txt(c(sprintf("MR Studio reproducibility capsule, written %s", format(Sys.time())), sprintf("MR Studio version: %s", .mrs_ver_build()), R.version.string, sprintf("Platform: %s", R.version$platform), sprintf("TwoSampleMR commit pinned by the adapters: %s", ADAPTER_IVW_MRE_V11$packages$TwoSampleMR), "", "Each pair folder holds its scripts (02_scripts), data (01_data), results (04_tables), console output and replay check. Open a pair folder and run 02_scripts/00_run_all.R from it."))
    f <- file.path(pk, "reproducibility_capsule.zip"); zip_write(f, ent); add(f)
  }
  strobe <- NULL
  if ("strobe" %in% items) {
    .say("Filling the STROBE-MR checklist: each of the 20 items marked as filled from your study, partly filled, or for you to write")
    rows <- strobe_rows(ctxs, txt, tryCatch(ui_state_get(h$con, "plan") %||% list(), error = function(e) list()), items)
    f <- file.path(pk, "strobe_mr_checklist.docx"); n <- strobe_write(f, rows, basename(h$dir)); add(f)
    strobe <- list(filled = n[["FILLED"]], partly = n[["PARTLY"]], author = n[["AUTHOR"]])
  }
  readme <- c("MR Studio publication package", sprintf("Study: %s", basename(h$dir)), sprintf("Built: %s", format(Sys.time())), "", "Everything here is a draft for the authors to check. Numbers come from the saved result tables; texts are assembled from the recorded settings.",
              if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) "At least one pair uses SIMULATED demonstration traits: do not report these results as findings." else NULL, "", "Files:", made, if (length(note)) c("", "Notes:", note) else NULL)
  rf <- file.path(pk, "README.txt"); writeLines(readme, rf, useBytes = TRUE)
  list(dir = pk, folder = basename(pk), files = made, notes = note, simulated = any(vapply(ctxs, `[[`, TRUE, "simulated")), strobe = strobe)
}

#' Round 33: the Methods sentence of the allele-frequency check of both datasets of a pair.
eaf_check_text <- function(E, anc = list()) {
  one <- function(e, a = NULL) {
    if (is.null(e)) return(NULL)
    if (!isTRUE(e$run)) { an <- .eaf_anc(e, a); if (!is.na(an)) return(sprintf("The effect-allele frequencies of %s were not compared with a reference panel, because its ancestry (%s) is not one of the five 1000 Genomes ancestries.", e$trait, an))   # Round 101: one plain sentence
      w <- sub("\\.?$", ".", e$why %||% "no reason was recorded."); return(sprintf("The effect-allele frequencies of %s were not compared with a reference panel (%s", e$trait, sub("\\.$", ").", paste0(tolower(substr(w, 1, 1)), substring(w, 2))))) }
    j <- as.integer(e$judged %||% 0L); f <- as.integer(e$flagged %||% 0L); mx <- suppressWarnings(as.numeric(e$max_difference %||% NA))
    sprintf("The effect-allele frequencies of %s were compared with the 1000 Genomes %s reference panel: of the %d variant%s that could be judged (not palindromic, panel minor allele frequency below %s), %s (largest difference %s).%s",
            e$trait, e$pop, j, if (j == 1L) "" else "s", format(e$maf_limit %||% EAF_FLIP_MAF_LIMIT),
            if (f) sprintf("%d differed by more than %s, which suggests an allele or strand error", f, format(e$threshold %||% EAF_FLIP_THRESHOLD)) else sprintf("none differed by more than %s", format(e$threshold %||% EAF_FLIP_THRESHOLD)),
            if (is.finite(mx)) sprintf("%.2f", mx) else "none",
            if (f) " Flagged variants were kept in the analysis and are listed in allele_frequency_check.csv." else "")
  }
  ax <- if (!is.null(E$x) && !isTRUE(E$x$run)) .eaf_anc(E$x, anc$x) else NA; ay <- if (!is.null(E$y) && !isTRUE(E$y$run)) .eaf_anc(E$y, anc$y) else NA
  if (!is.na(ax) && identical(ax, ay)) return(sprintf("The effect-allele frequencies of %s and %s were not compared with a reference panel, because their ancestry (%s) is not one of the five 1000 Genomes ancestries.", E$x$trait, E$y$trait, ax))   # Round 101: both in one sentence
  paste(c(one(E$x, anc$x), one(E$y, anc$y)), collapse = " ")
}
# Round 101: the ancestry named in an "is not one of the five 1000 Genomes ancestries" reason (filled in when the ledger left it empty), else NA
.eaf_anc <- function(e, a = NULL) { w <- e$why %||% ""; if (!grepl("is not one of the five 1000 Genomes ancestries", w, fixed = TRUE)) return(NA_character_)
  m <- regmatches(w, regexpr("\\(([^()]*)\\) is not one", w)); an <- if (length(m)) sub("^\\(([^()]*)\\).*$", "\\1", m) else ""; if (!nzchar(an)) an <- a %||% "not reported"; an }

#' Round 26: the harmonisation sentence when both files were declared to report alleles on the forward strand (action 1).
.strand_forward_text <- function(m, L, og = FALSE) {   # Round 100: og, both datasets are OpenGWAS files
  m <- gsub("(action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous)",
            if (og) "(action 1: both datasets are OpenGWAS files, which OpenGWAS harmonises to the forward strand of the reference genome, so palindromic variants were kept and aligned by their alleles alone, without a check of allele frequency)"
            else "(action 1: the analyst declared that both datasets report alleles on the forward strand, so palindromic variants were kept and aligned by their alleles alone, without a check of allele frequency)", m, fixed = TRUE)
  k <- as.integer(L$kept_palindromic %||% 0L)
  gsub(sprintf("%d palindromic variants were dropped", as.integer(L$dropped_palindromic %||% 0L)),
       sprintf("%d palindromic variant%s %s kept", k, if (k == 1L) "" else "s", if (k == 1L) "was" else "were"), m, fixed = TRUE)
}
#' Round 109: I2GX of a pair's kept instruments (Bowden 2016), from the report context.
simex_i2gx <- function(cx) tryCatch(feas_i2gx(cx$kept$bx, cx$kept$sx), error = function(e) NA_real_)
#' Round 109: the publication of a GWAS for the report: "(Willer, 2013; PMID 24097068)", from the OpenGWAS record or the PMID given at import.
gwas_cite <- function(s, lead = "") {
  if (is.null(s)) return(""); g <- function(k) { v <- s[[k]]; if (is.null(v) || !length(v) || is.na(v[[1]])) "" else trimws(as.character(v[[1]])) }
  a <- sub("[ ,].*$", "", g("author")); y <- sub("\\.0+$", "", g("year")); p <- sub("^PMID:? *", "", g("pmid"), ignore.case = TRUE)
  parts <- c(paste(c(if (nzchar(a)) a, if (nzchar(y)) y), collapse = ", "), if (nzchar(p) && grepl("^[0-9]+$", p)) paste("PMID", p))
  parts <- parts[nzchar(parts)]; if (!length(parts)) "" else paste0(lead, "(", paste(parts, collapse = "; "), ")")
}

# ---- Round 109: further figures of every univariable pair ----
# MendelianRandomization's own mr_forest, mr_funnel and mr_plot (so the layout is the package's), the mode estimators' density
# (as mrrobust mrmodalplot), the MR-Egger plot with its confidence band (as mrrobust mreggerplot) and the SIMEX extrapolation (as
# mrrobust mreggersimex). Each is a list(key, files, caption); a figure that cannot be drawn is left out.
.MRF_METHODS <- c(ivw_mre = "ivw", smed = "median", wm = "wmedian", egger = "egger", maxlik = "maxlik", mbe = "mbe", mode = "mbe", conmix = "conmix")
#' Round 110: MendelianRandomization's mr_forest prints its plot itself (print(forest)), which opened R's default device and wrote Rplots.pdf
#' into the working folder. The package's plot calls run with a device that writes nothing; the device that was current stays current.
.fig_quiet <- function(expr) {
  prev <- grDevices::dev.cur(); grDevices::pdf(NULL); d <- grDevices::dev.cur()
  on.exit({ if (d %in% grDevices::dev.list()) grDevices::dev.off(d); if (prev > 1L && prev %in% grDevices::dev.list()) grDevices::dev.set(prev) }, add = TRUE)
  expr
}
extra_figures <- function(cx, pk, tag, wmm, ci_level = 0.95) {
  k <- cx$kept; if (is.null(k) || nrow(k) < 3L) return(list())
  if (!requireNamespace("MendelianRandomization", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) return(list())
  mri <- MendelianRandomization::mr_input(bx = k$bx, bxse = k$sx, by = k$by, byse = k$sy, snps = k$id, exposure = cx$p$exName, outcome = cx$p$outName)
  al <- 1 - ci_level; out <- list(); ran <- unique(cx$ce$method %||% character())
  put <- function(key, h, g, cap) { f <- .fig_write(file.path(pk, "figures", paste0(tag, "_", key)), wmm, h, function() print(g)); if (length(f)) out[[length(out) + 1L]] <<- list(key = key, files = f, caption = cap) }
  pc <- .pct(ci_level); nm <- sprintf("%s on %s", cx$p$exName, cx$p$outName)
  g <- tryCatch(.fig_quiet(mrf_restyle(MendelianRandomization::mr_forest(mri, alpha = al, ordered = TRUE), wmm)), error = function(e) NULL)
  if (!is.null(g)) put("mr_forest_variants", max(60, 20 + 3.2 * nrow(k)), g, sprintf("Variant forest plot (MendelianRandomization mr_forest) of %s: the ratio estimate of each variant with its %s confidence interval, ordered by estimate, and the IVW estimate (diamond).", nm, pc))
  mm <- unique(unname(.MRF_METHODS[intersect(names(.MRF_METHODS), ran)]))
  g <- if (length(mm)) tryCatch(.fig_quiet(.fig_seeded(mrf_methods_plot(mri, al, mm, wmm))), error = function(e) NULL) else NULL
  if (!is.null(g)) put("mr_forest_methods", round(wmm / 1.617), g, sprintf("Methods forest plot (MendelianRandomization mr_forest) of %s: the estimate and %s confidence interval of each method that MendelianRandomization can draw among those run.", nm, pc))
  g <- tryCatch(.fig_quiet(mrf_restyle(MendelianRandomization::mr_funnel(mri, CI = TRUE), wmm)), error = function(e) NULL)
  if (!is.null(g)) put("mr_funnel", round(wmm / 1.617), g, sprintf("Funnel plot (MendelianRandomization mr_funnel) of %s: each variant's ratio estimate (x) against its precision, |beta_x| / se_y (y), with 95%% intervals; solid line at zero, dashed line at the IVW estimate. It differs from the funnel plot above, which has the standard error on a reversed axis.", nm))
  g <- tryCatch(.fig_quiet(mrf_restyle(MendelianRandomization::mr_plot(mri, interactive = FALSE, labels = FALSE, line = "ivw", orientate = TRUE), wmm)), error = function(e) NULL)
  if (!is.null(g)) put("mr_plot", round(wmm / 1.617), g, sprintf("Scatter plot (MendelianRandomization mr_plot) of %s, variants oriented to a positive effect on the exposure, with error bars on both axes and the IVW line (variant labels are in the supplement tables).", nm))
  g <- tryCatch(.fig_quiet(mrf_restyle(MendelianRandomization::mr_plot(.fig_seeded(MendelianRandomization::mr_allmethods(mri, method = "all", alpha = al)), interactive = FALSE), wmm)), error = function(e) NULL)   # fixed random state: the robust fits sample at random
  if (!is.null(g)) put("mr_allmethods", round(wmm / 1.617), g, sprintf("All-methods plot (MendelianRandomization mr_allmethods, method = \"all\", and mr_plot) of %s: one line per estimator (medians, IVW and MR-Egger with their penalised and robust versions).", nm))
  if (any(c("mode", "smode", "smode_nome", "wmode_nome", "mbe") %in% ran)) { g <- mode_density_plot(k); if (!is.null(g)) put("mode_density", round(wmm * 0.6), g, sprintf("Density of the single-variant ratio estimates of %s (as mrrobust mrmodalplot), smoothed with bandwidth multipliers phi = 0.25, 0.5 and 1; the mode-based estimates are the peaks [Hartwig 2017].", nm)) }
  if (any(c("egger", "radegger", "egger_simex", "egger_boot", "egger_rob") %in% ran)) { g <- egger_plot(k, ci_level); if (!is.null(g)) put("mregger", round(wmm * 0.6), g, sprintf("MR-Egger plot of %s (as mrrobust mreggerplot): variants oriented to a positive effect on the exposure with %s confidence intervals, the MR-Egger line with its own intercept and the shaded %s confidence band.", nm, pc, pc)) }
  kn <- tryCatch(knigge_figures(k, sprintf("%s ~ %s", cx$p$outName, cx$p$exName)), error = function(e) NULL)   # Round 109: item 11
  if (!is.null(kn)) { qn <- sprintf("Flagged in red: the %d variant(s) whose contribution to Cochran's Q exceeds the chi-squared cut-off at 0.05 / k (MR Studio's Q step).", kn$n_flag)
    put("knigge_mr_plot", round(wmm * 0.85), kn$fig5, sprintf("MR plot of %s (layout of mr.plot, matthijsknigge/mendelianrandomization, redrawn): variant effects with 95%% error bars on both axes, coloured by the ratio z, with IVW and MR-Egger lines before and after Q filtering, marginal densities and rug marks. %s", nm, qn))
    put("knigge_forest", max(60, round(25 + 4 * (nrow(k) + 4))), kn$fig6, sprintf("Forest plot of %s (layout of mr.forest.plot, redrawn): each variant's ratio estimate with its 95%% CI, IVW and MR-Egger before and after Q filtering at the top, the IVW interval shaded and a dotted line at zero. %s", nm, qn))
    put("knigge_funnel_ivw", round(wmm * 0.75), kn$fig7, sprintf("Funnel plot for IVW of %s (layout of mr.funnel.plot, redrawn): ratio estimate against its standard error (reversed axis) with funnel limits around the IVW estimate before and after Q filtering. %s", nm, qn))
    put("knigge_funnel_egger", round(wmm * 0.75), kn$fig8, sprintf("Funnel plot for MR-Egger of %s (layout of mr.funnel.plot, redrawn), as the IVW funnel plot around the MR-Egger estimate. %s", nm, qn)) }
  if ("egger_simex" %in% ran) { g <- simex_plot(k); if (!is.null(g)) put("mregger_simex", round(wmm * 0.6), g, sprintf("SIMEX for MR-Egger of %s (as mrrobust mreggersimex): the MR-Egger slope (left) and intercept (right) against lambda; open circle the original estimate, dots the simulated estimates, solid line the quadratic fit, dashed line its extrapolation to lambda = -1 and cross the SIMEX estimate.", nm)) }
  out
}
# ---- Round 110: figures of every multivariable pair (item 12) ----
#' The variants of a multivariable pair as its methods used them: every exposure and the outcome present, kept by harmonisation (as .mv_core).
.mv_fig_data <- function(h, p) {
  f <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "01_data", "analysis_data_*.csv")); if (!length(f)) return(NULL)
  d <- utils::read.csv(f[1], stringsAsFactors = FALSE); K <- sum(grepl("^beta_x_[0-9]+$", names(d))); if (K < 2L) return(NULL)
  bxn <- sprintf("beta_x_%d", seq_len(K)); sxn <- sprintf("se_x_%d", seq_len(K)); keep <- if ("mr_keep" %in% names(d)) as.logical(d$mr_keep) %in% TRUE else rep(TRUE, nrow(d))
  use <- keep & is.finite(d$beta_y) & is.finite(d$se_y) & d$se_y > 0 & stats::complete.cases(d[, c(bxn, sxn)])
  d <- d[use, c("variant_id", bxn, sxn, "beta_y", "se_y"), drop = FALSE]; rownames(d) <- NULL; d
}
#' The MendelianRandomization MRMVInput of a multivariable pair (cx$mvdat), with the exposure labels of the pair.
.mv_mri <- function(cx) {
  d <- cx$mvdat; K <- sum(grepl("^beta_x_[0-9]+$", names(d))); labs <- as.character(cx$labels %||% character()); if (length(labs) != K) labs <- sprintf("Exposure %d", seq_len(K))
  MendelianRandomization::mr_mvinput(bx = as.matrix(d[, sprintf("beta_x_%d", seq_len(K)), drop = FALSE]), bxse = as.matrix(d[, sprintf("se_x_%d", seq_len(K)), drop = FALSE]),
    by = d$beta_y, byse = d$se_y, snps = d$variant_id, exposure = labs, outcome = cx$p$outName %||% "the outcome")
}
#' Round 110a: per-variant contributions to Cochran's Q of MVMR-IVW from the variants of a multivariable pair (cx$mvdat), largest first, with
#' the total, its degrees of freedom and P, and the Bonferroni cut-off (chi-squared 1 df at 0.05 / n).
mv_q_contrib <- function(d) {
  K <- sum(grepl("^beta_x_[0-9]+$", names(d))); bx <- as.matrix(d[, sprintf("beta_x_%d", seq_len(K)), drop = FALSE]); n <- nrow(d); if (n <= K) return(NULL)
  b <- stats::coef(stats::lm(d$beta_y ~ bx - 1, weights = d$se_y^-2)); q <- as.numeric((d$beta_y - bx %*% b)^2 / d$se_y^2); o <- order(-q)
  list(rows = data.frame(variant_id = as.character(d$variant_id[o]), q = q[o], share = q[o] / sum(q), stringsAsFactors = FALSE), total = sum(q), df = n - K, n = n,
       p = stats::pchisq(sum(q), n - K, lower.tail = FALSE), cut = stats::qchisq(1 - 0.05 / n, 1))
}
#' MendelianRandomization's mr_plot of an MRMVInput (static): each variant's association with the outcome against its fitted value from
#' the MVMR-IVW model (by regressed on all bx, through the origin, weights se_y^-2), 95% error bars on both axes (the horizontal ones are the
#' confidence interval of the fitted value) and the line through the origin with gradient 1. Same list(key, files, caption) as extra_figures.
mv_extra_figures <- function(cx, pk, tag, wmm, ci_level = 0.95) {
  d <- cx$mvdat; if (is.null(d) || nrow(d) < 3L) return(list())
  if (!requireNamespace("MendelianRandomization", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) return(list())
  mvi <- .mv_mri(cx); out <- list()
  put <- function(key, h, g, cap) { f <- .fig_write(file.path(pk, "figures", paste0(tag, "_", key)), wmm, h, function() print(g)); if (length(f)) out[[length(out) + 1L]] <<- list(key = key, files = f, caption = cap) }
  g <- tryCatch(.fig_quiet(mrf_restyle(MendelianRandomization::mr_plot(mvi, interactive = FALSE, error = TRUE, line = TRUE, orientate = FALSE, labels = FALSE), wmm)), error = function(e) NULL)
  if (!is.null(g)) { wr <- function(x, n) if (is.character(x)) paste(strwrap(x, n), collapse = "\n") else x   # the package's axis titles are long: wrapped, at 11 pt, so they fit
    g$labels$x <- wr(g$labels$x, 70); g$labels$y <- wr(g$labels$y, 40); g <- g + ggplot2::theme(axis.title = ggplot2::element_text(size = 11)) }
  if (!is.null(g)) put("mvmr_plot", round(wmm * 0.75), g, sprintf("Multivariable scatter plot (MendelianRandomization mr_plot of the multivariable input) of %s on %s, %d variants: each variant's association with the outcome (y) against its fitted value from the MVMR-IVW model of all exposures (x), with 95%% error bars on both axes (horizontal: the confidence interval of the fitted value) and the line through the origin with gradient 1. Variants far from the line fit the model poorly. Variant labels are in the supplementary data.",
    paste(mvi@exposure, collapse = ", "), mvi@outcome, nrow(d)))
  # Round 110a: each variant's contribution to Cochran's Q of MVMR-IVW (as the Results table "Which variants drive the heterogeneity?"):
  # q_j = (beta_y - sum_k b_k beta_x_k)^2 / se_y^2 with b the MVMR-IVW estimate (first-order weights, as mr_mvivw); about chi-squared 1 df
  # each without heterogeneity, so the Bonferroni cut-off at 0.05 / n marks the outlying variants. At most the 30 largest are drawn.
  qd <- tryCatch(mv_q_contrib(d), error = function(e) NULL)
  if (!is.null(qd)) { top <- utils::head(qd$rows, 30L); top$variant_id <- factor(top$variant_id, rev(top$variant_id)); top$over <- top$q > qd$cut
    g <- ggplot2::ggplot(top, ggplot2::aes(q, variant_id, fill = over)) + ggplot2::geom_col(width = 0.7) + ggplot2::geom_vline(xintercept = qd$cut, linetype = 2, linewidth = 0.4) +
      ggplot2::scale_fill_manual(NULL, values = c(`FALSE` = "#4682b4", `TRUE` = "#c0392b"), labels = c(`FALSE` = "Below the cut-off", `TRUE` = "Above the cut-off"), drop = FALSE) +
      ggplot2::labs(x = "Contribution to Cochran's Q (MVMR-IVW)", y = NULL) + ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom", panel.grid.major.y = ggplot2::element_blank())
    put("mvmr_q", max(60, round(22 + 3.6 * nrow(top))), g, sprintf("Contribution of each variant to Cochran's Q of MVMR-IVW for %s on %s (the %d largest of %d variants): (beta_y minus the fitted value)^2 / se_y^2. Total Q = %s on %d degrees of freedom (P = %s). Dashed line: the Bonferroni cut-off, the chi-squared (1 df) quantile at 0.05 / %d = %s; variants beyond it (red, %d) disagree with the model more than chance allows. Heterogeneity is reported, not removed.",
      paste(mvi@exposure, collapse = ", "), mvi@outcome, nrow(top), qd$n, formatC(qd$total, digits = 1, format = "f"), qd$df, format(signif(qd$p, 3)), qd$n, formatC(qd$cut, digits = 2, format = "f"), sum(qd$rows$q > qd$cut))) }
  # Round 110a: radial MVMR (RMVMR, Spiller, Bowden and Sanderson; optional package): each exposure's variants on the radial scale, without
  # and with correction for the other exposures, coloured by exposure. Drawn when RMVMR is installed.
  if (requireNamespace("RMVMR", quietly = TRUE)) {
    rp <- tryCatch(.fig_quiet(suppressMessages(suppressWarnings({ fit <- RMVMR::ivw_rmvmr(mvi, summary = FALSE); RMVMR::plot_rmvmr(mvi, fit) }))), error = function(e) NULL)
    nm <- function(g) { if (is.null(g)) return(NULL); K <- length(mvi@exposure); pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#0072B2", "#CC79A7")[seq_len(K)]   # RMVMR's palette
      suppressMessages(g + ggplot2::scale_colour_manual(values = stats::setNames(pal, paste0("Exposure_", seq_len(K))), labels = stats::setNames(mvi@exposure, paste0("Exposure_", seq_len(K))), name = "Exposure")) }   # the exposures by name (RMVMR writes Exposure_1, ...)
    rp$p1 <- nm(rp$p1); rp$p2 <- nm(rp$p2)
    if (!is.null(rp$p1)) put("mvmr_radial", round(wmm * 0.75), mrf_restyle(rp$p1, wmm), sprintf("Radial MVMR plot without correction (RMVMR plot_rmvmr) of %s on %s: for each exposure, its variants (F > 10) on the radial scale, sqrt(W_j) against beta_j sqrt(W_j), with the exposure's line; variants far from their line contribute most to heterogeneity.", paste(mvi@exposure, collapse = ", "), mvi@outcome))
    if (!is.null(rp$p2)) put("mvmr_radial_corrected", round(wmm * 0.75), mrf_restyle(rp$p2, wmm), sprintf("Radial MVMR plot with correction (RMVMR plot_rmvmr) of %s on %s: each variant's estimate corrected for its associations with the other exposures, on the radial scale and coloured by the exposure it was chosen for.", paste(mvi@exposure, collapse = ", "), mvi@outcome)) }
  out
}
#' Round 110: the interactive versions of MendelianRandomization's mr_plot (the package's default, interactive = TRUE: plotly, the variant
#' named on hover) for the Results screen; the package and the report keep the static figures. Each is list(key, caption, json), json the
#' plotly figure (plotly::plotly_json). A univariable pair gets the scatter plot (as the static mr_plot figure: oriented, IVW line), a
#' multivariable pair the plot of the outcome associations against the fitted values. mr_allmethods has no interactive version.
interactive_figures <- function(cx) {
  if (!requireNamespace("MendelianRandomization", quietly = TRUE) || !requireNamespace("plotly", quietly = TRUE)) return(list())
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)   # ggplotly measures text on a device; none is left open or written
  one <- function(key, p, cap) { j <- tryCatch(as.character(plotly::plotly_json(p, jsonedit = FALSE)), error = function(e) NULL); if (is.null(j)) NULL else list(key = key, caption = cap, json = j) }
  out <- if (isTRUE(cx$mv)) { d <- cx$mvdat; if (is.null(d) || nrow(d) < 3L) list() else { mvi <- .mv_mri(cx)
      list(one("mvmr_plot_interactive", MendelianRandomization::mr_plot(mvi, interactive = TRUE, error = TRUE, line = TRUE, orientate = FALSE),
        sprintf("Interactive multivariable scatter plot (MendelianRandomization mr_plot, interactive) of %s on %s: hover a point for the variant; outcome association against its fitted value from MVMR-IVW, 95%% error bars on both axes, line through the origin with gradient 1. Pinch or drag a box to zoom; when zoomed in, scroll with two fingers to move around; double-click to reset.", paste(mvi@exposure, collapse = ", "), mvi@outcome))) } }
    else { k <- cx$kept; if (is.null(k) || nrow(k) < 3L) list() else { mri <- MendelianRandomization::mr_input(bx = k$bx, bxse = k$sx, by = k$by, byse = k$sy, snps = k$id, exposure = cx$p$exName, outcome = cx$p$outName)
      list(one("mr_plot_interactive", MendelianRandomization::mr_plot(mri, interactive = TRUE, labels = FALSE, line = "ivw", orientate = TRUE),
        sprintf("Interactive scatter plot (MendelianRandomization mr_plot, interactive) of %s on %s: hover a point for its variant (rs number); variants oriented to a positive effect on the exposure, 95%% error bars on both axes and the IVW line. Pinch or drag a box to zoom; when zoomed in, scroll with two fingers to move around; double-click to reset.", cx$p$exName, cx$p$outName))) } }
  Filter(Negate(is.null), out)
}
#' Round 110: plotly.js as shipped with the R package plotly (offline; sent to the screen once, when the first interactive figure is shown).
plotly_js_text <- function() { f <- system.file("htmlwidgets", "lib", "plotlyjs", "plotly-latest.min.js", package = "plotly"); if (!nzchar(f)) return(NULL); paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n") }
#' Round 109: density of the ratio estimates at three bandwidths (Hartwig 2017 bandwidth h = phi * 0.9 * min(sd, mad) * n^(-1/5)), weighted.
mode_density_plot <- function(k) {
  r <- k$by / k$bx; se <- abs(k$sy / k$bx); w <- (1 / se^2) / sum(1 / se^2); s <- min(stats::sd(r), stats::mad(r)); if (!is.finite(s) || s <= 0) return(NULL)
  xs <- seq(min(r) - 3 * s, max(r) + 3 * s, length.out = 512)
  d <- do.call(rbind, lapply(c(0.25, 0.5, 1), function(phi) { h <- phi * 0.9 * s * length(r)^(-1 / 5); data.frame(x = xs, y = vapply(xs, function(x) sum(w * stats::dnorm(x, r, h)), 0), phi = format(phi)) }))
  ggplot2::ggplot(d, ggplot2::aes(x, y, colour = phi)) + ggplot2::geom_line(linewidth = 0.7) + ggplot2::labs(x = "IV estimates", y = "Density", colour = "phi") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
#' Round 109: MR-Egger plot with the confidence band of the fitted line.
# Round 109: styled as mrrobust mreggerplot (black squares, capped CIs, black line, grey band).
egger_plot <- function(k, ci_level = 0.95) {
  s <- sign(k$bx); d <- data.frame(bx = abs(k$bx), by = k$by * s, sx = k$sx, sy = k$sy)
  f <- stats::lm(by ~ bx, data = d, weights = 1 / sy^2); z <- stats::qnorm(1 - (1 - ci_level) / 2); rse <- max(1, summary(f)$sigma)
  nx <- data.frame(bx = seq(0, max(d$bx + z * d$sx), length.out = 100)); pr <- stats::predict(f, nx, se.fit = TRUE); sf <- pr$se.fit / summary(f)$sigma * rse
  b <- data.frame(bx = nx$bx, fit = pr$fit, lo = pr$fit - z * sf, hi = pr$fit + z * sf)
  ggplot2::ggplot() + ggplot2::geom_ribbon(data = b, ggplot2::aes(bx, ymin = lo, ymax = hi, fill = "MR-Egger 95% CI"), alpha = 0.35) + ggplot2::geom_line(data = b, ggplot2::aes(bx, fit, colour = "MR-Egger"), linewidth = 0.6) +
    ggplot2::geom_errorbar(data = d, ggplot2::aes(bx, ymin = by - z * sy, ymax = by + z * sy, linetype = "95% CIs"), width = max(d$bx) * 0.012, colour = "black") + ggplot2::geom_point(data = d, ggplot2::aes(bx, by, shape = "Instruments"), size = 1.6) +
    ggplot2::scale_fill_manual(NULL, values = "grey70") + ggplot2::scale_colour_manual(NULL, values = "black") + ggplot2::scale_linetype_manual(NULL, values = 1) + ggplot2::scale_shape_manual(NULL, values = 15) + ggplot2::guides(shape = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 2), colour = ggplot2::guide_legend(order = 3), fill = ggplot2::guide_legend(order = 4)) +
    ggplot2::labs(x = "Instrument-exposure associations", y = "Instrument-outcome associations") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
#' Round 109: SIMEX of the MR-Egger slope and intercept (simex package, quadratic extrapolation, as the egger_simex method).
simex_plot <- function(k, B = 200L, seed = 1L) {
  if (!requireNamespace("simex", quietly = TRUE)) return(NULL)
  s <- sign(k$bx); d <- data.frame(bx = abs(k$bx), by = k$by * s, w = 1 / k$sy^2, sx = k$sx / k$sy)
  d$byw <- d$by / k$sy; d$bxw <- d$bx / k$sy; d$iw <- 1 / k$sy
  fit <- stats::lm(byw ~ iw + bxw - 1, data = d, x = TRUE, y = TRUE); set.seed(seed)
  sx <- tryCatch(simex::simex(fit, B = B, measurement.error = d$sx, SIMEXvariable = "bxw", fitting.method = "quad", asymptotic = FALSE, lambda = c(0.5, 1, 1.5, 2)), error = function(e) NULL)
  if (is.null(sx)) return(NULL)
  th <- sx$SIMEX.estimates; lam <- th[, 1]; mk <- function(j, lab) { y <- th[, j]; q <- stats::lm(y ~ lam + I(lam^2)); gx <- seq(-1, 2, length.out = 100)
    list(pts = data.frame(lambda = lam, est = y, what = lab, orig = lam == 0), fit = data.frame(lambda = gx, est = stats::predict(q, data.frame(lam = gx)), what = lab, extra = gx < 0),
         sim = data.frame(lambda = -1, est = stats::predict(q, data.frame(lam = -1)), what = lab)) }
  a <- mk(which(colnames(th) == "bxw"), "MR-Egger slope"); b <- mk(which(colnames(th) == "iw"), "MR-Egger intercept")
  P <- rbind(a$pts, b$pts); Fq <- rbind(a$fit, b$fit); Sm <- rbind(a$sim, b$sim); P$what <- factor(P$what, c("MR-Egger slope", "MR-Egger intercept")); Fq$what <- factor(Fq$what, levels(P$what)); Sm$what <- factor(Sm$what, levels(P$what))
  ggplot2::ggplot() + ggplot2::geom_line(data = Fq[!Fq$extra, ], ggplot2::aes(lambda, est, linetype = "Quadratic fit")) + ggplot2::geom_line(data = Fq[Fq$extra, ], ggplot2::aes(lambda, est, linetype = "Extrapolation")) +
    ggplot2::geom_point(data = P[!P$orig, ], ggplot2::aes(lambda, est, shape = "Simulated estimates")) + ggplot2::geom_point(data = P[P$orig, ], ggplot2::aes(lambda, est, shape = "Original estimate"), size = 2.5) +
    ggplot2::geom_point(data = Sm, ggplot2::aes(lambda, est, shape = "SIMEX estimate"), size = 3) + ggplot2::facet_wrap(~ what, scales = "free_y", strip.position = "left") + ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank()) +
    ggplot2::scale_shape_manual(NULL, values = c("Original estimate" = 1, "Simulated estimates" = 16, "SIMEX estimate" = 4)) + ggplot2::scale_linetype_manual(NULL, values = c("Quadratic fit" = 1, "Extrapolation" = 2), breaks = c("Quadratic fit", "Extrapolation")) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) + ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2), linetype = ggplot2::guide_legend(nrow = 2)) + ggplot2::theme(legend.box = "vertical") +
    ggplot2::labs(x = expression(lambda), y = NULL) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", strip.placement = "outside", strip.background = ggplot2::element_blank(), text = ggplot2::element_text(size = 9))   # Round 109: as mrrobust mreggersimex (per-panel y titles, lambda symbol)
}

# ---- Round 109: Figures 5 to 8 of matthijsknigge/mendelianrandomization (mr.plot, mr.forest.plot, mr.funnel.plot), redrawn ----
# Redrawn in MR Studio's own code with the same layout (no code reused). The "after Q" lines use MR Studio's Cochran's Q step: a variant is
# flagged when its contribution to Q, (by - b bx)^2 / se_y^2 with b the IVW estimate, exceeds the chi-squared (1 df) cut-off at 0.05 / k.
knigge_q_flags <- function(k) { w <- 1 / k$sy^2; b <- sum(w * k$bx * k$by) / sum(w * k$bx^2); q <- (k$by - b * k$bx)^2 / k$sy^2; q > stats::qchisq(1 - 0.05 / nrow(k), 1) }
.kn_fits <- function(k, keep) { f <- function(d, egger) { w <- 1 / d$sy^2; if (egger) { m <- stats::lm(by ~ bx, data = d, weights = w); c(a = unname(stats::coef(m)[1]), b = unname(stats::coef(m)[2]), se = summary(m)$coefficients[2, 2] / min(1, summary(m)$sigma)) } else { b <- sum(w * d$bx * d$by) / sum(w * d$bx^2); c(a = 0, b = b, se = sqrt(1 / sum(w * d$bx^2))) } }
  s <- sign(k$bx); d <- data.frame(bx = abs(k$bx), by = k$by * s, sy = k$sy); list(ivw = f(d, FALSE), ivw_q = f(d[keep, ], FALSE), egger = f(d, TRUE), egger_q = f(d[keep, ], TRUE)) }
knigge_figures <- function(k, title) {
  if (nrow(k) < 4L) return(NULL); G <- ggplot2::ggplot; A <- ggplot2::aes
  fl <- knigge_q_flags(k); keep <- !fl; if (sum(keep) < 3L) keep[] <- TRUE; fits <- .kn_fits(k, keep)
  s <- sign(k$bx); d <- data.frame(id = k$id, bx = abs(k$bx), by = k$by * s, sx = k$sx, sy = k$sy, flag = fl); d$ratio <- d$by / d$bx; d$rse <- d$sy / d$bx; d$z <- d$ratio / d$rse
  ln <- data.frame(name = c("IVW", "IVW (Q-filtered)", "MR-Egger", "MR-Egger (Q-filtered)"), a = c(fits$ivw[["a"]], fits$ivw_q[["a"]], fits$egger[["a"]], fits$egger_q[["a"]]), b = c(fits$ivw[["b"]], fits$ivw_q[["b"]], fits$egger[["b"]], fits$egger_q[["b"]]),
                   se = c(fits$ivw[["se"]], fits$ivw_q[["se"]], fits$egger[["se"]], fits$egger_q[["se"]]))
  th <- ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0))
  f5 <- G(d, A(bx, by)) + ggplot2::geom_errorbar(A(ymin = by - 1.96 * sy, ymax = by + 1.96 * sy), width = 0, colour = "#c9c4e8") + ggplot2::geom_errorbarh(A(xmin = bx - 1.96 * sx, xmax = bx + 1.96 * sx), height = 0, colour = "#c9c4e8") +
    ggplot2::geom_point(A(colour = z), data = d[!d$flag, ]) + ggplot2::geom_point(data = d[d$flag, ], colour = "red") + ggplot2::scale_colour_gradient(name = "z", low = "#c6dbef", high = "#08306b") +
    ggplot2::geom_abline(A(intercept = a, slope = b, linetype = name), data = ln, colour = "#3f3fa0") + ggplot2::scale_linetype_manual(NULL, values = c(1, 2, 4, 3)) + ggplot2::guides(linetype = ggplot2::guide_legend(nrow = 2)) + ggplot2::theme(legend.box = "vertical") + ggplot2::geom_rug(sides = "bl", colour = "steelblue", alpha = 0.6) +
    ggplot2::labs(title = title, x = "SNP-exposure effect", y = "SNP-outcome effect") + th
  f5 <- tryCatch({ px <- G(d, A(bx)) + ggplot2::geom_density(fill = "steelblue", colour = NA) + ggplot2::theme_void(); py <- G(d, A(by)) + ggplot2::geom_density(fill = "steelblue", colour = NA) + ggplot2::coord_flip() + ggplot2::theme_void()
    g <- cowplot::insert_xaxis_grob(f5, cowplot::get_panel(px), grid::unit(0.15, "null"), position = "top"); cowplot::ggdraw(cowplot::insert_yaxis_grob(g, cowplot::get_panel(py), grid::unit(0.15, "null"), position = "right")) }, error = function(e) f5)
  fo <- rbind(data.frame(label = ln$name, est = ln$b, lo = ln$b - 1.96 * ln$se, hi = ln$b + 1.96 * ln$se, kind = "method", flag = FALSE), data.frame(label = d$id, est = d$ratio, lo = d$ratio - 1.96 * d$rse, hi = d$ratio + 1.96 * d$rse, kind = "snp", flag = d$flag))
  fo$label <- factor(fo$label, rev(fo$label))
  f6 <- G(fo, A(est, label)) + ggplot2::annotate("rect", xmin = ln$b[1] - 1.96 * ln$se[1], xmax = ln$b[1] + 1.96 * ln$se[1], ymin = -Inf, ymax = Inf, alpha = 0.35, fill = "steelblue") + ggplot2::annotate("rect", xmin = ln$b[3] - 1.96 * ln$se[3], xmax = ln$b[3] + 1.96 * ln$se[3], ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#7a3fc7") +   # Round 109: IVW (blue) and MR-Egger (purple) bands, as mr.forest.plot
    ggplot2::geom_vline(xintercept = 0, linetype = 3) +
    ggplot2::geom_errorbarh(A(xmin = lo, xmax = hi, colour = flag), height = 0) + ggplot2::geom_point(A(colour = flag, shape = kind)) + ggplot2::geom_point(data = fo[fo$kind == "method", ], colour = c("steelblue", "steelblue", "#3f3fa0", "#3f3fa0"), shape = 18, size = 3.5) + ggplot2::guides(shape = "none") + ggplot2::scale_colour_manual(NULL, values = c(`FALSE` = "black", `TRUE` = "red"), labels = c("kept", "flagged by Q")) +
    ggplot2::scale_shape_manual(NULL, values = c(method = 18, snp = 16)) + ggplot2::labs(title = title, x = "Causal estimate (95% CI)", y = NULL) + th
  fun <- function(m, lab) { b0 <- ln$b[ln$name == m]; b1 <- ln$b[ln$name == paste0(m, " (Q-filtered)")]; ys <- seq(0, max(d$rse) * 1.05, length.out = 50)
    fn <- rbind(data.frame(y = ys, lo = b0 - 1.96 * ys, hi = b0 + 1.96 * ys, which = "all variants"), data.frame(y = ys, lo = b1 - 1.96 * ys, hi = b1 + 1.96 * ys, which = "after Q filtering"))
    G(d, A(ratio, rse)) + ggplot2::geom_ribbon(A(y = y, xmin = -Inf, xmax = lo), data = fn[fn$which == "all variants", ], inherit.aes = FALSE, orientation = "y", fill = "grey88") + ggplot2::geom_ribbon(A(y = y, xmin = hi, xmax = Inf), data = fn[fn$which == "all variants", ], inherit.aes = FALSE, orientation = "y", fill = "grey88") + ggplot2::geom_line(A(lo, y, linetype = which), data = fn) + ggplot2::geom_line(A(hi, y, linetype = which), data = fn) + ggplot2::geom_vline(xintercept = c(b0, b1), linetype = c(1, 2)) +
      ggplot2::geom_point(A(colour = z), data = d[!d$flag, ]) + ggplot2::geom_point(data = d[d$flag, ], colour = "red") + ggplot2::scale_colour_gradient(name = "z", low = "#c6dbef", high = "#08306b") +
      ggplot2::scale_y_reverse() + ggplot2::geom_rug(sides = "bl", colour = "steelblue", alpha = 0.6) + ggplot2::scale_linetype_manual(NULL, values = c(1, 2)) + ggplot2::labs(title = paste0("Funnel, ", lab), x = "Causal estimate (beta_IV)", y = "Standard error") + th }
  list(fig5 = f5, fig6 = f6, fig7 = fun("IVW", "IVW"), fig8 = fun("MR-Egger", "MR-Egger"), n_flag = sum(fl))
}

#' Round 109 (item 9): the mode-based estimates at bandwidth multipliers phi = 0.25, 0.5 and 1, simple and weighted, each with
#' second-order (NOME not assumed) and first-order (NOME) ratio standard errors, as mrrobust mrmodal reports them
#' (MendelianRandomization::mr_mbe; bootstrap SEs with a fixed seed, so SEs differ by software).
mode_table <- function(k, ci_level = 0.95, iterations = 1000L, seed = 12345L) {
  if (nrow(k) < 3L || !requireNamespace("MendelianRandomization", quietly = TRUE)) return(NULL)
  mri <- MendelianRandomization::mr_input(bx = k$bx, bxse = k$sx, by = k$by, byse = k$sy)
  g <- expand.grid(phi = c(0.25, 0.5, 1), weighting = c("unweighted", "weighted"), stderror = c("delta", "simple"), stringsAsFactors = FALSE)   # mr_mbe calls the simple mode "unweighted"
  rows <- lapply(seq_len(nrow(g)), function(i) { f <- MendelianRandomization::mr_mbe(mri, weighting = g$weighting[i], stderror = g$stderror[i], phi = g$phi[i], seed = seed, iterations = iterations, alpha = 1 - ci_level)
    data.frame(estimator = paste(if (g$weighting[i] == "unweighted") "Simple mode" else "Weighted mode", if (g$stderror[i] == "simple") "(NOME)" else ""), phi = g$phi[i],
               estimate = f@Estimate, se = f@StdError, ci_lower = f@CILower, ci_upper = f@CIUpper, p = f@Pvalue, stringsAsFactors = FALSE) })
  out <- do.call(rbind, rows); out$estimator <- trimws(out$estimator); out[order(out$estimator, out$phi), ]
}

# ---- Round 109: a figures script in every pair folder (05_figures/make_figures.R), so the shipped figures can be redrawn in R ----
# While the package draws its figures, .fig_write records the data of each drawing (FIG_REC). Then every pair gets make_figures.R:
# the exact plotting functions of this MR Studio (dumped from the running engine, so the app and the script cannot differ), the
# recorded inputs (figure_inputs.rds) and the pair's instruments and estimates (figure_data.csv, figure_estimates.csv).
FIG_REC <- new.env()
FIG_NAMESPACES <- c("ggplot2", "MendelianRandomization", "cowplot", "simex", "RMVMR")   # Round 110a: RMVMR (radial MVMR, optional); loaded before any figure is drawn (loading them changes how text is drawn)
#' Draws a Greek letter once off-screen, so the symbol font is loaded before the first figure (its first use could differ).
fig_warmup <- function() tryCatch({ f <- tempfile(fileext = ".png"); grDevices::png(f, type = "cairo"); grid::grid.text(expression(phi == lambda)); grDevices::dev.off(); unlink(f) }, error = function(e) NULL)
#' Evaluates expr with a fixed random state, restored afterwards: robust fits (robustbase::lmrob in MendelianRandomization) draw random subsamples from the session state,
#' so without this the package and make_figures.R could draw them differently.
.fig_seeded <- function(expr, seed = 314159265L) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE); if (had) old <- get(".Random.seed", envir = globalenv()); ok <- RNGkind()
  on.exit({ RNGkind(ok[1], ok[2], ok[3]); if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()) })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); expr
}
.fig_record <- function(base, w_mm, h_mm, draw) {
  if (!isTRUE(FIG_REC$on)) return(invisible())
  env <- environment(draw); vars <- intersect(codetools::findGlobals(draw, merge = FALSE)$variables, ls(env, all.names = TRUE))
  FIG_REC$items[[length(FIG_REC$items) + 1L]] <- list(base = basename(base), w_mm = w_mm, h_mm = h_mm, body = body(draw), data = mget(vars, envir = env))
}
#' Every function and constant the given functions use, found by following their global names in env (the engine's environment).
.fig_deps <- function(start, env = environment(extra_figures)) {
  todo <- start; seen <- character()
  while (length(todo)) { n <- todo[1]; todo <- todo[-1]; if (n %in% seen || n %in% c(".fig_record", "FIG_REC") || !exists(n, envir = env, inherits = FALSE)) next
    seen <- c(seen, n); obj <- get(n, envir = env)
    if (is.function(obj)) { g <- codetools::findGlobals(obj); todo <- c(todo, setdiff(g, seen)) } }
  seen
}
write_figures_script <- function(h, cx, i, wmm, items, xkeys) {
  pd <- file.path(h$dir, "03_pairs", cx$p$folder); if (!dir.exists(pd)) return(invisible(NULL)); fd <- file.path(pd, "05_figures"); dir.create(fd, showWarnings = FALSE)
  tag <- sprintf("P%d", i); its <- Filter(function(x) startsWith(x$base, paste0(tag, "_")) && !(sub(paste0("^", tag, "_"), "", x$base) %in% xkeys), items)
  saveRDS(its, file.path(fd, "figure_inputs.rds"), version = 2)
  mv <- isTRUE(cx$mv); xfun <- if (mv) "mv_extra_figures" else "extra_figures"   # Round 110: a multivariable pair gets its script too
  utils::write.csv(if (mv) cx$mvdat else cx$kept, file.path(fd, "figure_data.csv"), row.names = FALSE); utils::write.csv(cx$ce, file.path(fd, "figure_estimates.csv"), row.names = FALSE)
  saveRDS(if (mv) list(mvdat = cx$mvdat, ce = cx$ce, labels = cx$labels) else list(kept = cx$kept, ce = cx$ce), file.path(fd, "figure_pair.rds"), version = 2)   # exact values (the CSV files round to 15 digits)
  calls <- unique(unlist(lapply(its, function(x) all.names(x$body)))); deps <- .fig_deps(unique(c(".fig_write", xfun, "fig_warmup", calls)))
  df <- tempfile(fileext = ".R"); dump(deps, file = df, envir = environment(extra_figures)); defs <- readLines(df, warn = FALSE); unlink(df)
  defs <- vapply(enc2utf8(defs), function(s) { cp <- utf8ToInt(s); if (all(cp < 128L)) return(s); paste(vapply(cp, function(c) if (c < 128L) intToUtf8(c) else sprintf("\\u%04x", c), ""), collapse = "") }, "", USE.NAMES = FALSE)   # ASCII only: other characters as \u escapes, so the script reads the same in any locale
  pk <- c("ggplot2", "MendelianRandomization", "cowplot", "simex", "RMVMR", "codetools"); ver <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), "")
  txt <- c("# make_figures.R: redraws every figure of this pair exactly as in the publication package.",
    sprintf("# Written by MR Studio %s with %s. Packages: %s.", .mrs_ver_build(), R.version.string, paste(sprintf("%s %s", pk, ver), collapse = ", ")),
    "# Run from this folder:  Rscript make_figures.R   (the figures are written to ./figures as PDF, SVG and PNG at 300 dpi)",
    "# The plotting functions below are the ones MR Studio used; figure_inputs.rds holds the data of each figure as drawn,",
    if (mv) "# figure_data.csv the variants (beta_x_k and se_x_k of each exposure k, beta_y, se_y) and figure_estimates.csv the estimates of the pair." else "# figure_data.csv the instruments (bx, by, sx, sy) and figure_estimates.csv the estimates of the pair.", "",
    "if (!isTRUE(l10n_info()[[\"UTF-8\"]])) for (l in c(\"C.UTF-8\", \"en_US.UTF-8\", \"English_United States.utf8\")) if (nzchar(suppressWarnings(Sys.setlocale(\"LC_CTYPE\", l)))) break   # Greek letters (phi) need a UTF-8 locale", "",
    "# ---- plotting functions of MR Studio ----", defs, "",
    sprintf("for (ns in %s) suppressWarnings(requireNamespace(ns, quietly = TRUE))   # as MR Studio, before any figure", paste(deparse(FIG_NAMESPACES), collapse = "")),
    "fig_warmup()",
    "# ---- the figures ----",
    "dir.create(\"figures\", showWarnings = FALSE)",
    "for (f in readRDS(\"figure_inputs.rds\")) { e <- list2env(f$data, parent = globalenv()); d <- eval(call(\"function\", NULL, f$body), e)",
    "  .fig_write(file.path(\"figures\", f$base), f$w_mm, f$h_mm, d) }",
    sprintf("cx <- c(readRDS(\"figure_pair.rds\"), list(p = list(exName = %s, outName = %s)))   # the instruments and estimates, exact (also in the CSV files)", deparse(cx$p$exName), deparse(cx$p$outName)),
    sprintf("invisible(%s(cx, \".\", %s, %d, %s))", xfun, deparse(tag), as.integer(wmm), format((cx$ce$ci_level %||% 0.95)[1])))
  writeLines(txt, file.path(fd, "make_figures.R"), useBytes = TRUE)
  invisible(file.path(fd, "make_figures.R"))
}

# ---- Round 109: figure styling matched to the originals (audited element by element against mrrobust, mendelianrandomization and MendelianRandomization; private_notes/r109_figcheck/AUDIT.md) ----
# Proposed styling for the Round 109 figures: drop-in replacements for egger_plot, simex_plot, mode_density_plot and knigge_figures in
# prototype/engine/report_core.R (signatures unchanged), plus wrappers for the MendelianRandomization figures (mrf_restyle, mrf_methods_plot).
# ASCII only: Greek letters only through plotmath (a Unicode phi prints as <U+03C6> in non-UTF-8 sessions).
# Reproducibility: every plot with more than one legend sets each guide's order explicitly (ggplot otherwise orders them
# differently between sessions), and simex_plot fixes RNGkind before set.seed, so a figures script redraws byte-identical output. Source report_core.R first, then this file.
#
# Scale: MR Studio writes these figures w_mm = 140 mm wide (vector PDF and SVG plus a 300 dpi PNG, by .fig_write). Every size is matched to
# the original at equal figure width ("same zoom"):
#  - mrrobust (Stata SVG, 5400 units across, from mreggerplot.ado / mreggersimex.ado / mrmodalplot.ado as exported on the mrrobust site):
#      1 unit = 140 / 5400 mm, ggplot linewidth = mm / 0.7528. 9.72 u -> 0.335, 7.29 u -> 0.251, 6.48 u -> 0.223, 4.86 u -> 0.167;
#      text 112.56 u -> 8.3 pt, 89.94 u -> 6.6 pt, 84.38 u -> 6.2 pt.
#  - Knigge (matthijsknigge/mendelianrandomization R/mr.plot.R, mr.forest.plot.R, mr.funnel.plot.R; README PNGs drawn at 72 dpi:
#      Fig 5 916 px, Fig 6 806 px, Figs 7-8 1140 px wide): every size of the package code times f = 140 mm / original width
#      (0.433, 0.492, 0.348). Text is kept at 6.5 pt or more so the report stays readable.
#  - MendelianRandomization 0.10.0 (mr_plot, mr_forest, mr_funnel; reference figures drawn at 7.29 x 4.51 in by pkgdown):
#      mrf_restyle scales the package's own ggplot by s = 140 / 185.2.

.MRR_LW <- c(thick = 0.335, fit = 0.251, axis = 0.223, thin = 0.167)
.MRR_GRID <- "#F0F0F0"
# Stata number format: no leading zero (".1", "-.05").
.stata_lab <- function(x) vapply(x, function(v) if (is.na(v)) NA_character_ else sub("^(-?)0\\.", "\\1.", format(signif(v, 6), trim = TRUE, drop0trailing = TRUE, scientific = FALSE)), "")
# Theme of mrrobust's Stata graphs: white, black L axes with outside ticks, dashed #F0F0F0 major grid on both axes, legend below.
.mrr_theme <- function(base = 8.3, legend = 6.6, grid = .MRR_LW[["thick"]]) {
  ggplot2::theme_classic(base_size = base, base_line_size = .MRR_LW[["axis"]], base_rect_size = .MRR_LW[["axis"]]) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = .MRR_GRID, linewidth = grid, linetype = "73"), panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black", linewidth = .MRR_LW[["axis"]]), axis.ticks = ggplot2::element_line(colour = "black", linewidth = .MRR_LW[["axis"]]),
                   axis.text = ggplot2::element_text(colour = "black", size = base), axis.title = ggplot2::element_text(colour = "black", size = base),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA), legend.position = "bottom", legend.key = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = legend), legend.title = ggplot2::element_blank())
}
# RColorBrewer "Blues"[2:9] as Knigge's scale_colour_gradientn, and brewer "Blues"[7] for the "p < 0.05" key.
.KN_BLUES <- c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
.KN_DARKBLUE <- "#00008B"; .KN_PURPLE <- "#A020F0"; .KN_STEEL <- "#4682B4"; .KN_GREY <- "#BEBEBE"; .KN_DARKGREY <- "#A9A9A9"; .KN_KEYBLUE <- "#2171B5"
.KN_PT <- c(fig5 = 2.35, fig6 = 1.06, funnel = 1.89)
.KN_ZLAB <- expression(frac(hat(beta)[IV], sigma[hat(beta)[IV]]))
.KN_CHILAB <- expression(chi[homogeneity]^2)
# theme_minimal as in Knigge's code, line sizes scaled by f, text floored at 5.5 pt; bold = TRUE for the 14 pt bold axis text of mr.plot and mr.forest.plot.
# Point sizes: the originals are 72 dpi PNGs whose antialiasing widens points, so sizes are set from the measured diameters
# (Fig 5 13 of 916 px, Figs 7-8 13 of 1140 px, Fig 6 6 of 806 px) rather than from size x f: .KN_PT.
# Text factor: the originals' 72 dpi hinting draws text about 1.38 times larger than its point size at the same zoom (tick glyph
# heights measured: Fig 5 15.3 vs 10.9, Fig 6 17.4 vs 12.7 per 1000 px of width), as with the points (.KN_PT). Legend keys and
# spacing follow the base size, as in theme_minimal; the base text has a 5 pt floor.
.KN_TXT <- 1.38
.kn_theme <- function(f, bold = FALSE) {
  b <- max(5, 11 * f * .KN_TXT)
  t <- ggplot2::theme_minimal(base_size = b, base_line_size = 0.5 * f, base_rect_size = 0.5 * f) + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0),
    legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.margin = ggplot2::margin(0, 0, 0, 0), legend.spacing.x = grid::unit(0.6 * b, "pt"), legend.key.spacing.x = grid::unit(0.15 * b, "pt"), legend.box.spacing = grid::unit(0.3 * b, "pt"),
    legend.key.width = grid::unit(1.6 * b, "pt"), legend.key.height = grid::unit(1.6 * b, "pt"), legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 0.3 * b, r = 0.3 * b)))
  if (bold) t <- t + ggplot2::theme(axis.text = ggplot2::element_text(size = 14 * f * .KN_TXT, face = "bold", colour = "grey30"), axis.title = ggplot2::element_text(size = 14 * f * .KN_TXT, face = "bold"))
  t
}

mode_density_plot <- function(k) {
  r <- k$by / k$bx; se <- abs(k$sy / k$bx); w <- (1 / se^2) / sum(1 / se^2); s <- min(stats::sd(r), stats::mad(r)); if (!is.finite(s) || s <= 0) return(NULL)
  xs <- seq(min(r) - 3 * s, max(r) + 3 * s, length.out = 512)
  d <- do.call(rbind, lapply(c(0.25, 0.5, 1), function(phi) { h <- phi * 0.9 * s * length(r)^(-1 / 5); data.frame(x = xs, y = vapply(xs, function(x) sum(w * stats::dnorm(x, r, h)), 0), phi = format(phi)) }))
  ggplot2::ggplot(d, ggplot2::aes(x, y, colour = phi)) + ggplot2::geom_line(linewidth = .MRR_LW[["thick"]]) +
    ggplot2::scale_colour_manual(NULL, values = c("0.25" = "#1A85FF", "0.5" = "#D41159", "1" = "#00BF7F"), breaks = c("0.25", "0.5", "1"), labels = parse(text = c("phi == \".25\"", "phi == \".5\"", "phi == 1"))) +
    ggplot2::scale_x_continuous(labels = .stata_lab) + ggplot2::scale_y_continuous(labels = .stata_lab) +
    ggplot2::labs(x = "IV estimates", y = "Density", colour = NULL) + .mrr_theme(legend = 8.3)
}

egger_plot <- function(k, ci_level = 0.95) {
  s <- sign(k$bx); d <- data.frame(bx = abs(k$bx), by = k$by * s, sx = k$sx, sy = k$sy)
  f <- stats::lm(by ~ bx, data = d, weights = 1 / sy^2); z <- stats::qnorm(1 - (1 - ci_level) / 2); rse <- max(1, summary(f)$sigma)
  nx <- data.frame(bx = seq(0, max(d$bx + z * d$sx), length.out = 100)); pr <- stats::predict(f, nx, se.fit = TRUE); sf <- pr$se.fit / summary(f)$sigma * rse
  b <- data.frame(bx = nx$bx, fit = pr$fit, lo = pr$fit - z * sf, hi = pr$fit + z * sf)
  cap <- 0.0036 * 1.1 * max(nx$bx)   # mrrobust rcap: 17 of 5400 units across the graph, about 0.36 % of the x span
  ggplot2::ggplot() + ggplot2::geom_ribbon(data = b, ggplot2::aes(bx, ymin = lo, ymax = hi, fill = "MR-Egger 95% CI"), alpha = 0.5) +
    ggplot2::geom_errorbar(data = d, ggplot2::aes(bx, ymin = by - z * sy, ymax = by + z * sy, linetype = "95% CIs"), width = cap, colour = "black", linewidth = .MRR_LW[["thin"]]) +
    ggplot2::geom_point(data = d, ggplot2::aes(bx, by, shape = "Instruments"), size = 0.5, stroke = 0.5, colour = "black") + ggplot2::geom_line(data = b, ggplot2::aes(bx, fit, colour = "MR-Egger"), linewidth = .MRR_LW[["thick"]]) +
    ggplot2::scale_fill_manual(NULL, values = "#CCCCCC") + ggplot2::scale_colour_manual(NULL, values = "black") + ggplot2::scale_linetype_manual(NULL, values = 1) + ggplot2::scale_shape_manual(NULL, values = 15) +
    ggplot2::scale_x_continuous(labels = .stata_lab) + ggplot2::scale_y_continuous(labels = .stata_lab) +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1, override.aes = list(size = 1.2)), linetype = ggplot2::guide_legend(order = 2), colour = ggplot2::guide_legend(order = 3), fill = ggplot2::guide_legend(order = 4)) +
    ggplot2::labs(x = "Instrument-exposure associations", y = "Instrument-outcome associations") + .mrr_theme()
}

# SIMEX as mrreggersimex: two panels side by side, each with its own y title and lambda x title, one shared legend below
# (SIMEX / Original / Simulated over Quadratic fit / Extrapolation). Falls back to a faceted plot when cowplot is missing.
simex_plot <- function(k, B = 200L, seed = 1L) {
  if (!requireNamespace("simex", quietly = TRUE)) return(NULL)
  s <- sign(k$bx); d <- data.frame(bx = abs(k$bx), by = k$by * s, w = 1 / k$sy^2, sx = k$sx / k$sy)
  d$byw <- d$by / k$sy; d$bxw <- d$bx / k$sy; d$iw <- 1 / k$sy
  fit <- stats::lm(byw ~ iw + bxw - 1, data = d, x = TRUE, y = TRUE); RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  sx <- tryCatch(simex::simex(fit, B = B, measurement.error = d$sx, SIMEXvariable = "bxw", fitting.method = "quad", asymptotic = FALSE, lambda = c(0.5, 1, 1.5, 2)), error = function(e) NULL)
  if (is.null(sx)) return(NULL)
  th <- sx$SIMEX.estimates; lam <- th[, 1]; mk <- function(j, lab) { y <- th[, j]; q <- stats::lm(y ~ lam + I(lam^2)); gx <- seq(-1, 2, length.out = 100)
    list(pts = data.frame(lambda = lam, est = y, what = lab, orig = lam == 0), fit = data.frame(lambda = gx, est = stats::predict(q, data.frame(lam = gx)), what = lab, extra = gx < 0),
         sim = data.frame(lambda = -1, est = stats::predict(q, data.frame(lam = -1)), what = lab)) }
  a <- mk(which(colnames(th) == "bxw"), "MR-Egger slope"); b <- mk(which(colnames(th) == "iw"), "MR-Egger intercept")
  P <- rbind(a$pts, b$pts); Fq <- rbind(a$fit, b$fit); Sm <- rbind(a$sim, b$sim); P$what <- factor(P$what, c("MR-Egger slope", "MR-Egger intercept")); Fq$what <- factor(Fq$what, levels(P$what)); Sm$what <- factor(Sm$what, levels(P$what))
  marks <- function(P, Fq, Sm) list(ggplot2::geom_line(data = Fq[!Fq$extra, ], ggplot2::aes(lambda, est, linetype = "Quadratic fit"), linewidth = .MRR_LW[["fit"]]), ggplot2::geom_line(data = Fq[Fq$extra, ], ggplot2::aes(lambda, est, linetype = "Extrapolation"), linewidth = .MRR_LW[["fit"]]),
    ggplot2::geom_point(data = P[!P$orig, ], ggplot2::aes(lambda, est, shape = "Simulated"), size = 0.8, stroke = 0.5), ggplot2::geom_point(data = P[P$orig, ], ggplot2::aes(lambda, est, shape = "Original"), size = 0.62, stroke = 0.38),
    ggplot2::geom_point(data = Sm, ggplot2::aes(lambda, est, shape = "SIMEX"), size = 0.96, stroke = 0.25),
    ggplot2::scale_shape_manual(NULL, values = c("Original" = 1, "Simulated" = 16, "SIMEX" = 4), breaks = c("SIMEX", "Original", "Simulated")), ggplot2::scale_linetype_manual(NULL, values = c("Quadratic fit" = 1, "Extrapolation" = 2), breaks = c("Quadratic fit", "Extrapolation")),
    ggplot2::scale_x_continuous(breaks = -1:2, labels = .stata_lab, expand = ggplot2::expansion(mult = 0.05)), ggplot2::scale_y_continuous(labels = .stata_lab),
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.5, stroke = 0.5)), linetype = ggplot2::guide_legend(order = 2, nrow = 1)),
    .mrr_theme(base = 6.2, legend = 6.6, grid = .MRR_LW[["fit"]]), ggplot2::theme(legend.box = "vertical", legend.spacing.y = grid::unit(0, "pt"), legend.margin = ggplot2::margin(0, 0, 0, 0)))
  pane <- function(w) ggplot2::ggplot() + marks(P[P$what == w, ], Fq[Fq$what == w, ], Sm[Sm$what == w, ]) + ggplot2::labs(x = expression(lambda), y = w)
  tryCatch({ p1 <- pane("MR-Egger slope"); p2 <- pane("MR-Egger intercept"); lg <- cowplot::get_plot_component(p1, "guide-box-bottom", return_all = FALSE)
    nl <- ggplot2::theme(legend.position = "none")
    cowplot::plot_grid(cowplot::plot_grid(p1 + nl, p2 + nl, nrow = 1), lg, ncol = 1, rel_heights = c(1, 0.24)) },
    error = function(e) ggplot2::ggplot() + marks(P, Fq, Sm) + ggplot2::facet_wrap(~ what, scales = "free_y", strip.position = "left") + ggplot2::labs(x = expression(lambda), y = NULL) +
      ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 6.2, colour = "black")))
}

# Knigge Figures 5 to 8 (mr.plot, mr.forest.plot, mr.funnel.plot), redrawn with the package's settings. Variants flagged by MR Studio's
# Q step are drawn as Knigge draws chochran.Q == 0 (red, legend "p > 0.05"); with no flagged variant the package's no-Q branch is used.
knigge_figures <- function(k, title) {
  if (nrow(k) < 4L) return(NULL); G <- ggplot2::ggplot; A <- ggplot2::aes
  fl <- knigge_q_flags(k); keep <- !fl; if (sum(keep) < 3L) keep[] <- TRUE; fits <- .kn_fits(k, keep); hq <- any(!keep)
  s <- sign(k$bx); d <- data.frame(id = k$id, bx = abs(k$bx), by = k$by * s, sx = k$sx, sy = k$sy, flag = fl & !keep); d$ratio <- d$by / d$bx; d$rse <- d$sy / d$bx; d$z <- d$ratio / d$rse
  d$q <- factor(ifelse(d$flag, "flag", "keep"), c("flag", "keep"))
  ln <- data.frame(name = c("IVW", "IVW.Q", "MRegger", "MRegger.Q"), a = c(fits$ivw[["a"]], fits$ivw_q[["a"]], fits$egger[["a"]], fits$egger_q[["a"]]), b = c(fits$ivw[["b"]], fits$ivw_q[["b"]], fits$egger[["b"]], fits$egger_q[["b"]]),
                   se = c(fits$ivw[["se"]], fits$ivw_q[["se"]], fits$egger[["se"]], fits$egger_q[["se"]]))
  zcol <- function(o) ggplot2::scale_colour_gradientn(name = .KN_ZLAB, colours = .KN_BLUES, guide = ggplot2::guide_colourbar(order = o))
  qshape <- function(f, o) list(ggplot2::scale_shape_manual(name = .KN_CHILAB, values = c(flag = 19, keep = 19), limits = c("flag", "keep"), labels = c("p > 0.05", "p < 0.05"), drop = FALSE),
    ggplot2::guides(shape = ggplot2::guide_legend(order = o, override.aes = list(shape = 21, colour = "white", fill = c("red", .KN_KEYBLUE), size = 5 * f * 1.35, alpha = 1))))
  pts <- function(f, x, y, ps) list(ggplot2::geom_point(A(.data[[x]], .data[[y]], colour = z, shape = q), data = d, size = ps, stroke = 0.5 * f, alpha = 0.8),
    ggplot2::geom_point(A(.data[[x]], .data[[y]]), data = d[d$flag, ], colour = "red", size = ps * 4.1 / 4, stroke = 0.5 * f),
    ggplot2::geom_rug(A(.data[[x]], .data[[y]], colour = z), data = d[!d$flag, ], sides = "bl", linewidth = 0.5 * f), ggplot2::geom_rug(A(.data[[x]], .data[[y]]), data = d[d$flag, ], colour = "red", sides = "bl", linewidth = 0.5 * f))
  # Fig 5 (mr.plot): f = 0.433. Error bars darkblue alpha 0.3; IVW dashed, MR-Egger solid, darkblue, before-Q lines alpha 0.2 when Q removed variants;
  # zero lines darkgrey; marginal densities steelblue (flagged variants red, alpha 0.3) on the axes, 0.2 of the panel.
  f <- 0.433; ln$alpha <- if (hq) c(0.2, 1, 0.2, 1) else 1
  p5 <- G(d, A(bx, by)) + ggplot2::geom_hline(yintercept = 0, colour = .KN_DARKGREY, linewidth = 0.5 * f) + ggplot2::geom_vline(xintercept = 0, colour = .KN_DARKGREY, linewidth = 0.5 * f) +
    ggplot2::geom_errorbar(A(ymin = by - 1.96 * sy, ymax = by + 1.96 * sy), width = 0, colour = .KN_DARKBLUE, alpha = 0.3, linewidth = 0.5 * f) + ggplot2::geom_errorbarh(A(xmin = bx - 1.96 * sx, xmax = bx + 1.96 * sx), height = 0, colour = .KN_DARKBLUE, alpha = 0.3, linewidth = 0.5 * f) +
    ggplot2::geom_abline(A(intercept = a, slope = b, linetype = name, alpha = I(alpha)), data = ln, colour = .KN_DARKBLUE, linewidth = 0.5 * f) + pts(f, "bx", "by", .KN_PT[["fig5"]]) + zcol(3) + qshape(f, 1) +
    ggplot2::scale_linetype_manual("Method", values = c(IVW = "dashed", IVW.Q = "dashed", MRegger = "solid", MRegger.Q = "solid")) +
    ggplot2::guides(linetype = ggplot2::guide_legend(order = 2, nrow = 1, override.aes = list(alpha = if (hq) c(0.2, 1, 0.2, 1) else 1, colour = .KN_DARKBLUE))) +
    ggplot2::labs(title = title, x = expression(hat(beta)[exposure]), y = expression(hat(beta)[outcome])) + .kn_theme(f, bold = TRUE) + ggplot2::theme(legend.box = "horizontal")
  # Marginals: steelblue densities of the kept variants. Knigge's code adds a red (alpha 0.3) geom_density of the p > 0.05 variants,
  # which the ggplot2 of the original figure drew as a low band over just their range; it is drawn here as that band (7 % of the blue peak).
  f5 <- tryCatch({ dens <- function(v) { g <- list(ggplot2::geom_density(A(x = val), data = data.frame(val = d[[v]][!d$flag]), fill = .KN_STEEL, colour = "black", alpha = 1, linewidth = 0.2 * f))
      if (any(d$flag)) { r <- range(d[[v]][d$flag]); if (diff(r) == 0) r <- r + c(-1, 1) * 0.005 * diff(range(d[[v]]))
        g <- c(g, list(ggplot2::annotate("rect", xmin = r[1], xmax = r[2], ymin = 0, ymax = 0.07 * max(stats::density(d[[v]][!d$flag])$y), fill = "red", colour = "black", alpha = 0.3, linewidth = 0.2 * f))) }
      g }
    xd <- cowplot::axis_canvas(p5, axis = "x") + dens("bx") + ggplot2::scale_y_continuous(expand = c(0, 0))
    yd <- cowplot::axis_canvas(p5, axis = "y", coord_flip = TRUE) + dens("by") + ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::coord_flip()
    g <- cowplot::insert_xaxis_grob(p5, xd, grid::unit(0.2, "null"), position = "top"); cowplot::ggdraw(cowplot::insert_yaxis_grob(g, yd, grid::unit(0.2, "null"), position = "right")) }, error = function(e) p5)
  # Fig 6 (mr.forest.plot): f = 0.492. Rows bottom to top: variants in data order, IVW, Egger and, with Q, IVW / Egger Cochran's Q. Points size 1.5,
  # error bars alpha 0.5 in the row colour, dotted zero line, no legend. Bands to the top variant row: with Q Egger grey 0.6, IVW grey 0.4,
  # Egger Q steelblue 0.6, IVW Q purple 0.4 (as the package draws them); without Q Egger purple 0.6, IVW darkblue 0.4.
  f <- 0.492; mo <- if (hq) c(1, 3, 2, 4) else c(1, 3)   # ln rows bottom to top: IVW, Egger, IVW Q, Egger Q
  mlab <- c("IVW", "IVW Cochran's Q", "Egger", "Egger Cochran's Q"); mcol <- if (hq) c(.KN_GREY, .KN_DARKBLUE, .KN_GREY, .KN_PURPLE) else c(.KN_DARKBLUE, NA, .KN_PURPLE, NA)
  fo <- rbind(data.frame(label = d$id, est = d$ratio, lo = d$ratio - 1.96 * d$rse, hi = d$ratio + 1.96 * d$rse, col = ifelse(d$flag, "red", "black")),
              data.frame(label = mlab[mo], est = ln$b[mo], lo = ln$b[mo] - 1.96 * ln$se[mo], hi = ln$b[mo] + 1.96 * ln$se[mo], col = mcol[mo]))
  fo$label <- factor(fo$label, unique(fo$label))
  band <- function(i, fill, alpha) ggplot2::annotate("rect", xmin = ln$b[i] - 1.96 * ln$se[i], xmax = ln$b[i] + 1.96 * ln$se[i], ymin = -Inf, ymax = nrow(d), fill = fill, alpha = alpha, colour = NA)
  bands <- if (hq) list(band(3, .KN_GREY, 0.6), band(1, .KN_GREY, 0.4), band(4, .KN_STEEL, 0.6), band(2, .KN_PURPLE, 0.4)) else list(band(3, .KN_PURPLE, 0.6), band(1, .KN_DARKBLUE, 0.4))
  f6 <- G(fo, A(est, label)) + ggplot2::geom_point(A(colour = I(col)), size = .KN_PT[["fig6"]], stroke = 0.5 * f) + ggplot2::geom_vline(xintercept = 0, linetype = 3, linewidth = 0.5 * f) +
    ggplot2::geom_errorbarh(A(xmin = lo, xmax = hi, colour = I(col)), height = 0.5, alpha = 0.5, linewidth = 0.5 * f) + bands +
    ggplot2::labs(title = title, x = expression(beta[IV]), y = "") + .kn_theme(f, bold = TRUE) + ggplot2::theme(legend.position = "none")
  # Figs 7-8 (mr.funnel.plot): f = 0.348. 95 % limits (1.96 se) dotted and 3.29 se limits twodash; before Q grey, after Q darkblue alpha 0.6;
  # estimate lines darkblue (alpha 0.3 before Q), dashed for IVW and solid for MR-Egger; grey alpha 0.3 shading outside the outer limits within the data range.
  f <- 0.348
  fun <- function(m, mname) { i0 <- which(ln$name == m); i1 <- which(ln$name == paste0(m, ".Q")); ys <- seq(0, max(d$rse), length.out = 60); lty <- if (m == "IVW") "dashed" else "solid"
    lim <- function(b) data.frame(y = ys, l95 = b - 1.96 * ys, u95 = b + 1.96 * ys, l99 = b - 3.29 * ys, u99 = b + 3.29 * ys)
    f0 <- lim(ln$b[i0]); f1 <- lim(ln$b[i1]); sh <- if (hq) f1 else f0; xr <- range(c(d$ratio, f0$l99, f0$u99, if (hq) c(f1$l99, f1$u99)))
    l2 <- function(fr, v, col, alpha, lt) ggplot2::geom_line(A(x, y), data = data.frame(x = fr[[v]], y = fr$y), inherit.aes = FALSE, colour = col, alpha = alpha, linetype = lt, linewidth = 0.5 * f)
    g <- G(d, A(ratio, rse)) + ggplot2::geom_ribbon(A(y = y, xmin = xr[1], xmax = pmax(l99, xr[1])), data = sh, inherit.aes = FALSE, orientation = "y", fill = .KN_GREY, alpha = 0.3) +
      ggplot2::geom_ribbon(A(y = y, xmin = pmin(u99, xr[2]), xmax = xr[2]), data = sh, inherit.aes = FALSE, orientation = "y", fill = .KN_GREY, alpha = 0.3) +
      l2(f0, "l95", .KN_GREY, 1, "dotted") + l2(f0, "u95", .KN_GREY, 1, "dotted") + l2(f0, "l99", .KN_GREY, 1, "twodash") + l2(f0, "u99", .KN_GREY, 1, "twodash")
    if (hq) g <- g + l2(f1, "l95", .KN_DARKBLUE, 0.6, "dotted") + l2(f1, "u95", .KN_DARKBLUE, 0.6, "dotted") + l2(f1, "l99", .KN_DARKBLUE, 0.6, "twodash") + l2(f1, "u99", .KN_DARKBLUE, 0.6, "twodash")
    nm <- c(mname, paste0(mname, ".Q")); vl <- data.frame(x = c(ln$b[i0], ln$b[i1])[if (hq) 1:2 else 1], which = nm[if (hq) 1:2 else 1], alpha = c(0.3, 1)[if (hq) 1:2 else 1])
    g + ggplot2::geom_vline(A(xintercept = x, linetype = which, alpha = I(alpha)), data = vl, colour = .KN_DARKBLUE, linewidth = 0.5 * f) +
      ggplot2::scale_linetype_manual("Method", values = stats::setNames(c(lty, lty), nm), breaks = vl$which) + ggplot2::guides(linetype = ggplot2::guide_legend(order = 1, override.aes = list(alpha = c(0.1, 1)[seq_len(nrow(vl))]))) +
      pts(f, "ratio", "rse", .KN_PT[["funnel"]]) + zcol(2) + qshape(f, 3) + ggplot2::scale_y_reverse() +
      ggplot2::labs(title = title, x = expression(beta[IV]), y = expression(sigma[beta[IV]])) + .kn_theme(f) }
  list(fig5 = f5, fig6 = f6, fig7 = fun("IVW", "Inverse-Variance Weighted method"), fig8 = fun("MRegger", "MR-egger method"), n_flag = sum(fl))
}

# ---- MendelianRandomization figures (mr_plot, mr_forest, mr_funnel) ----
# The package draws for a 7.29 x 4.51 in device (its pkgdown reference figures). mrf_restyle scales every layer's line width, point size and
# stroke, and the theme's text and lines, by s = w_mm / 185.2, so a figure written w_mm wide looks as the reference does at the same zoom.
# Colours, shapes and line types are the package's own and are left unchanged.
# Where to call it, in extra_figures (report_core.R): wrap each package plot before put(), for example
#   g <- tryCatch(mrf_restyle(MendelianRandomization::mr_forest(mri, alpha = al, ordered = TRUE), wmm), error = function(e) NULL)
#   g <- tryCatch(mrf_methods_plot(mri, al, mm, wmm), error = function(e) NULL)
#   g <- tryCatch(mrf_restyle(MendelianRandomization::mr_funnel(mri, CI = TRUE), wmm), error = function(e) NULL)
#   g <- tryCatch(mrf_restyle(MendelianRandomization::mr_plot(mri, interactive = FALSE, labels = FALSE, line = "ivw", orientate = TRUE), wmm), error = function(e) NULL)
#   g <- tryCatch(mrf_restyle(MendelianRandomization::mr_plot(MendelianRandomization::mr_allmethods(mri, method = "all", alpha = al), interactive = FALSE), wmm), error = function(e) NULL)
# and give mr_plot, mr_allmethods, mr_funnel and the methods forest the reference aspect: height round(wmm / 1.617) (the variant forest keeps its height per variant).
MRF_ASPECT <- 7.29 / 4.51
mrf_restyle <- function(g, w_mm = 140, ref_mm = 185.2) {
  if (!inherits(g, "ggplot")) return(g); s <- w_mm / ref_mm
  for (L in g$layers) { p <- L$aes_params
    if (inherits(L$geom, "GeomPoint")) { L$aes_params$size <- (if (is.null(p$size)) 1.5 else p$size) * s; L$aes_params$stroke <- (if (is.null(p$stroke)) 0.5 else p$stroke) * s }
    else if (inherits(L$geom, c("GeomText", "GeomLabel"))) L$aes_params$size <- (if (is.null(p$size)) 3.88 else p$size) * s
    else { L$aes_params$linewidth <- (if (!is.null(p$linewidth)) p$linewidth else if (!is.null(p$size)) p$size else 0.5) * s; L$aes_params$size <- NULL } }
  th <- g$theme; if (is.null(th$text)) g <- g + ggplot2::theme(text = ggplot2::element_text(size = 11), line = ggplot2::element_line(linewidth = 0.5), rect = ggplot2::element_rect(linewidth = 0.5))
  mods <- list(legend.key.size = grid::unit(1.2 * s, "lines"))
  for (nm in names(g$theme)) { el <- g$theme[[nm]]; sz <- tryCatch(el$size, error = function(e) NULL); lw <- tryCatch(el$linewidth, error = function(e) NULL)
    if (inherits(el, "element_text") && is.numeric(sz) && !inherits(sz, "rel")) mods[[nm]] <- ggplot2::element_text(size = sz * s)
    if (inherits(el, "element_line") && is.numeric(lw) && !inherits(lw, "rel")) mods[[nm]] <- ggplot2::element_line(linewidth = lw * s)
    if (inherits(el, "element_rect") && is.numeric(lw) && !inherits(lw, "rel")) mods[[nm]] <- ggplot2::element_rect(linewidth = lw * s) }
  g + do.call(ggplot2::theme, mods)
}
# Methods forest plot (mr_forest with snp_estimates = FALSE), restyled; call in place of the second mr_forest in extra_figures.
mrf_methods_plot <- function(mri, alpha = 0.05, methods = c("ivw", "median", "egger"), w_mm = 140) mrf_restyle(MendelianRandomization::mr_forest(mri, alpha = alpha, snp_estimates = FALSE, methods = methods), w_mm)
