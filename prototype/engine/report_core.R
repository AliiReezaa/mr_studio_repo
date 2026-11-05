# Publication package for MR Studio (Round 9). Base R only: no packages beyond the ones the engine already needs.
# Everything written here is a DRAFT for the author to check: numbers come from the saved result tables of the study,
# texts are assembled from the settings recorded in the generated scripts, and every limitation that applies is stated.
# Nothing here is an estimator.

# ---------- small helpers ----------
.xe <- function(s) { s <- as.character(s); s[is.na(s)] <- ""; s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE); s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("\"", "&quot;", s, fixed = TRUE) }
.fnum <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
.fp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, digits = 2, format = "e"), formatC(p, digits = 3, format = "f")))
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
METHOD_REFS <- list(ivw_mre = "burgess2013", egger = "bowden2015", wm = "bowden2016", smed = "bowden2016", mode = "hartwig2017", smode = "hartwig2017", presso = "verbanck2018", radial = "bowden2018", radegger = "bowden2018", steiger = "hemani2017")
REFS$hemani2017 <- list(pmid = "29149188", key = "Hemani2017", authors = "Hemani G, Tilling K, Davey Smith G.", title = "Orienting the causal relationship between imprecisely measured traits using GWAS summary data", journal = "PLoS Genet", year = "2017", vol = "13", pages = "e1007081", doi = "10.1371/journal.pgen.1007081")
REFS$skrivankova2021ee <- list(pmid = "34702754", key = "Skrivankova2021EE", authors = "Skrivankova VW, Richmond RC, Woolf BAR, et al.", title = "Strengthening the reporting of observational studies in epidemiology using mendelian randomisation (STROBE-MR): explanation and elaboration", journal = "BMJ", year = "2021", vol = "375", pages = "n2233", doi = "10.1136/bmj.n2233")
REFS$chang2015 <- list(pmid = "25722852", key = "Chang2015", authors = "Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ.", title = "Second-generation PLINK: rising to the challenge of larger and richer datasets", journal = "Gigascience", year = "2015", vol = "4", pages = "7", doi = "10.1186/s13742-015-0047-8")
REFS$kgp2015 <- list(pmid = "26432245", key = "1000Genomes2015", authors = "1000 Genomes Project Consortium, Auton A, Brooks LD, et al.", title = "A global reference for human genetic variation", journal = "Nature", year = "2015", vol = "526", pages = "68-74", doi = "10.1038/nature15393")
REFS$lee2012 <- list(pmid = "22714935", key = "Lee2012", authors = "Lee SH, Goddard ME, Wray NR, Visscher PM.", title = "A better coefficient of determination for genetic profile analysis", journal = "Genet Epidemiol", year = "2012", vol = "36", pages = "214-224", doi = "10.1002/gepi.21614")
REFS$verbanck2018 <- list(pmid = "29686387", key = "Verbanck2018", authors = "Verbanck M, Chen CY, Neale B, Do R.", title = "Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases", journal = "Nat Genet", year = "2018", vol = "50", pages = "693-698", doi = "10.1038/s41588-018-0099-7")
REFS$bowden2018 <- list(pmid = "29961852", key = "Bowden2018", authors = "Bowden J, Spiller W, Del Greco M F, et al.", title = "Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression", journal = "Int J Epidemiol", year = "2018", vol = "47", pages = "1264-1278", doi = "10.1093/ije/dyy101")
ref_text <- function(r) sprintf("%s %s. %s. %s;%s:%s. doi:%s. PMID %s.", r$authors, r$title, r$journal, r$year, r$vol, r$pages, r$doi, r$pmid)
refs_bibtex <- function(ids) paste(vapply(ids, function(i) { r <- REFS[[i]]; sprintf("@article{%s,\n  author = {%s},\n  title = {%s},\n  journal = {%s},\n  year = {%s},\n  volume = {%s},\n  pages = {%s},\n  doi = {%s},\n  note = {PMID %s}\n}\n", r$key, gsub(", et al\\.", " and others", gsub("\\.$", "", gsub(", ", " and ", r$authors))), r$title, r$journal, r$year, r$vol, gsub("-", "--", r$pages), r$doi, r$pmid) }, ""), collapse = "\n")
refs_ris <- function(ids) paste(vapply(ids, function(i) { r <- REFS[[i]]; paste(c("TY  - JOUR", sprintf("AU  - %s", trimws(strsplit(gsub("\\.$", "", gsub(", et al", "", r$authors)), ", ")[[1]])), sprintf("TI  - %s", r$title), sprintf("JO  - %s", r$journal), sprintf("PY  - %s", r$year), sprintf("VL  - %s", r$vol), sprintf("SP  - %s", r$pages), sprintf("DO  - %s", r$doi), sprintf("AN  - %s", r$pmid), "ER  - ", ""), collapse = "\n") }, ""), collapse = "\n")

# ---------- figures: one drawing function, written as PDF, SVG and PNG so they cannot disagree ----------
.draw_forest <- function(est, title = NULL) {
  op <- graphics::par(mar = c(3.4, 8.6, 1.4, 1.2), mgp = c(1.9, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  k <- nrow(est); xl <- range(c(est$ci_lower, est$ci_upper, 0), na.rm = TRUE); pad <- diff(xl) * 0.06
  graphics::plot(NA, xlim = xl + c(-pad, pad), ylim = c(k + 0.6, 0.4), yaxt = "n", xlab = est$unit[1], ylab = "", bty = "n")
  graphics::abline(v = if (isTRUE(est$ratio[1])) 1 else 0, lty = 2, col = "grey55"); graphics::axis(2, at = seq_len(k), labels = est$label, tick = FALSE)
  graphics::segments(est$ci_lower, seq_len(k), est$ci_upper, seq_len(k), lwd = 1.6, col = est$col); graphics::points(est$estimate, seq_len(k), pch = 16, cex = 1.1, col = est$col)
  if (!is.null(title)) graphics::mtext(title, 3, 0.2, adj = 0, font = 2, cex = 0.8)
}
.draw_scatter <- function(v, lines, xlab, ylab) {
  op <- graphics::par(mar = c(3.4, 3.8, 0.8, 0.8), mgp = c(2, 0.6, 0), cex = 0.78, las = 1); on.exit(graphics::par(op))
  s <- ifelse(v$bx < 0, -1, 1); x <- v$bx * s; y <- v$by * s
  graphics::plot(x, y, pch = 16, col = "grey25", xlim = c(0, max(x) * 1.05), ylim = range(c(y, 0)) * 1.1, xlab = xlab, ylab = ylab, bty = "l")
  graphics::abline(h = 0, col = "grey75"); graphics::segments(x, y - 1.96 * v$sy, x, y + 1.96 * v$sy, col = "grey75", lwd = 0.6)
  for (l in lines) graphics::abline(a = l$a, b = l$b, col = l$col, lty = l$lty, lwd = 2)
  graphics::legend("topleft", legend = vapply(lines, `[[`, "", "name"), col = vapply(lines, `[[`, "", "col"), lty = vapply(lines, `[[`, 1, "lty"), lwd = 2, bty = "n", cex = 0.85)
}
.MCOL <- c(ivw_mre = "#2b3fc7", egger = "#b5541c", wm = "#2f7fd6", smed = "#8c4fd6", mode = "#d6407a", smode = "#2a9d8f", presso = "#b8860b", radial = "#3d9a38", radegger = "#a04a3d")
.fig_write <- function(base, w_mm, h_mm, draw) {
  w <- w_mm / 25.4; h <- h_mm / 25.4; out <- character()
  ok <- function(f, expr) { r <- tryCatch({ expr; TRUE }, error = function(e) FALSE); if (!r) { try(grDevices::dev.off(), silent = TRUE); if (file.exists(f)) unlink(f) }; if (r) out <<- c(out, f) }
  f <- paste0(base, ".pdf"); ok(f, { grDevices::pdf(f, width = w, height = h, useDingbats = FALSE); draw(); grDevices::dev.off() })
  if (isTRUE(unname(capabilities("cairo")))) {
    f <- paste0(base, ".svg"); ok(f, { grDevices::svg(f, width = w, height = h); draw(); grDevices::dev.off() })
    f <- paste0(base, ".png"); ok(f, { grDevices::png(f, width = w, height = h, units = "in", res = 300, type = "cairo"); draw(); grDevices::dev.off() })
  }
  out
}

# ---------- gathering what the report needs ----------
.est_df <- function(pay) {
  if (!length(pay$estimates)) return(NULL)
  d <- do.call(rbind, lapply(pay$estimates, function(e) data.frame(method = e$method, name = e$method_name, estimand = e$estimand, estimate = e$estimate, se = e$se, ci_lower = e$ci_lower, ci_upper = e$ci_upper, ci_level = e$ci_level, p = e$p, n_ivs = e$n_ivs_used, scale = e$scale %||% "beta", stringsAsFactors = FALSE)))
  d
}
.short <- c(steiger = "Steiger", ivw_mre = "IVW", egger = "MR-Egger", wm = "Weighted median", smed = "Simple median", mode = "Weighted mode", smode = "Simple mode", presso = "MR-PRESSO", radial = "Radial IVW", radegger = "Radial Egger")
pair_context <- function(h, p) {
  pay <- pair_payload(h$dir, p$folder); est <- .est_df(pay); ce <- est[est$estimand == "causal_effect", , drop = FALSE]
  v <- do.call(rbind, lapply(pay$variants, function(x) data.frame(id = x$id, bx = x$bx, by = x$by, sx = x$sx, sy = x$sy, keep = x$keep, stringsAsFactors = FALSE)))
  kept <- v[v$keep, , drop = FALSE]; Fs <- (kept$bx / kept$sx)^2
  st <- pay$settings %||% list(); dt <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "01_data", "analysis_data_*.csv"))
  dd <- if (length(dt)) utils::read.csv(dt[1], stringsAsFactors = FALSE) else NULL
  nx <- if (!is.null(dd$n_x)) stats::median(dd$n_x, na.rm = TRUE) else NA_real_; ny <- if (!is.null(dd$n_y)) stats::median(dd$n_y, na.rm = TRUE) else NA_real_
  rc <- file.path(h$dir, "03_pairs", p$folder, "02_scripts", "00_run_all.R"); sc <- Sys.glob(file.path(h$dir, "03_pairs", p$folder, "02_scripts", "[0-9]*_*.R"))
  srcname <- function(id, role) if (startsWith(id, "lib:")) tryCatch({ r <- library_dataset(sub("^lib:", "", id))$recipe; list(file = r$original$file, ancestry = r$facts$ancestry %||% "unknown", n = r$facts$n %||% "", og = if (identical(r$source$type, "opengwas")) r$source else NULL) }, error = function(e) list(file = "(file no longer in the library)", ancestry = "unknown", n = "")) else list(file = NA_character_, ancestry = "EUR", n = "")
  list(p = p, pay = pay, est = est, ce = ce, vars = v, kept = kept, mean_f = mean(Fs), min_f = min(Fs), n_weak = sum(Fs < 10), settings = st, nx = nx, ny = ny,
       simulated = isTRUE(p$simulated) || !startsWith(p$ex %||% "", "lib:"), ledger = p$ledger, pthr = p$pthr %||% 5e-8, src_x = srcname(p$ex, "exposure"), src_y = srcname(p$out, "outcome"), methods = unlist(p$methods), scripts = sc,
       overlap = tryCatch(overlap_get(h)[[overlap_key(p$ex %||% "", p$out %||% "")]], error = function(e) NULL))
}

# ---------- text ----------
.ci_phrase <- function(r, scale) { ratio <- scale %in% c("log_or"); f <- function(x) if (ratio) exp(x) else x
  sprintf("%s = %s (%s CI %s to %s; P = %s)", if (ratio) "OR" else "beta", .fnum(f(r$estimate)), .pct(r$ci_level), .fnum(f(r$ci_lower)), .fnum(f(r$ci_upper)), .fp(r$p)) }
#' The limitation sentences about sample overlap. Overlap is never calculated by MR Studio; the sentence says what the analyst declared.
overlap_limitations <- function(ctxs) {
  real <- Filter(function(cx) !isTRUE(cx$simulated), ctxs)
  bias <- "Overlap biases estimates toward the observational association, most strongly with weak instruments."
  untested <- "Sample overlap between the exposure and outcome studies was not tested. "
  if (!length(real)) return(paste0(untested, bias))
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
steiger_limitations <- function(ctxs) {
  f <- lapply(ctxs, steiger_facts); ran <- !vapply(f, is.null, TRUE)
  if (!any(ran)) return("Directionality was not tested (the Steiger directionality test was not run), so reverse causation is not excluded by this analysis.")
  bad <- vapply(seq_along(f), function(i) ran[i] && !(f[[i]]$ok && is.finite(f[[i]]$p) && f[[i]]$p < 0.05), TRUE)
  c("The Steiger directionality test [Hemani 2017] compares the variance the variants explain in the exposure and in the outcome. It assumes the variants act on the outcome only through the exposure, it is sensitive to measurement error and horizontal pleiotropy. It does not exclude confounding or reverse causation.",
    {
      pv <- unlist(lapply(ctxs[ran], function(cx) { fx <- steiger_facts(cx); c(if (fx$bin_x) sprintf("%s %s", cx$p$exName, .prev_of(fx$note_x)), if (fx$bin_y) sprintf("%s %s", cx$p$outName, .prev_of(fx$note_y))) }))
      if (length(pv)) sprintf("For binary traits the variance explained was computed on the liability scale [Lee 2012] and depends on the population prevalence entered by the analyst (%s); a different prevalence can change the result.", paste(pv, collapse = "; ")) else NULL
    },
    if (any(bad)) sprintf("For %s the test did not clearly support the assumed direction.", paste(vapply(ctxs[bad], function(cx) sprintf("%s -> %s", cx$p$exName, cx$p$outName), ""), collapse = "; ")),
    if (any(!ran)) sprintf("The test was not run for %s.", paste(vapply(ctxs[!ran], function(cx) sprintf("%s -> %s", cx$p$exName, cx$p$outName), ""), collapse = "; ")))
}
report_text <- function(ctxs, study_name) {
  txt <- list(methods = character(), results = character(), limitations = character(), availability = character(), legends = character())
  used <- unique(unlist(lapply(ctxs, function(x) x$methods))); anysim <- any(vapply(ctxs, `[[`, TRUE, "simulated")); anyimp <- any(!vapply(ctxs, `[[`, TRUE, "simulated"))
  mname <- c(ivw_mre = "inverse-variance weighted (IVW)", egger = "MR-Egger regression", wm = "weighted median", smed = "simple median", mode = "weighted mode", smode = "simple mode", presso = "MR-PRESSO", radial = "radial IVW", radegger = "radial MR-Egger")
  txt$methods <- c(
    "We performed a two-sample Mendelian randomization (MR) analysis using summary-level genetic associations, run in MR Studio, an application that writes every analysis step and setting into R scripts that can be re-run. The analysis plan follows the structure recommended in current MR guidelines [Burgess 2023]; reporting against the STROBE-MR checklist [Skrivankova 2021] has not been checked automatically and must be completed by the authors.")
  for (cx in ctxs) {
    p <- cx$p; st <- cx$settings
    if (cx$simulated) {
      txt$methods <- c(txt$methods, sprintf("Exposure %s and outcome %s. The exposure and the outcome are simulated demonstration traits generated by MR Studio. They do not describe any real population, and no result from them may be reported as a finding.", p$exName, p$outName))
    } else {
      sx <- cx$src_x; sy <- cx$src_y
      nstr <- function(s) if (nzchar(as.character(s$n))) sprintf(", sample size %s", format(as.numeric(s$n), big.mark = ",", scientific = FALSE)) else ""
      side <- function(s, what) if (!is.null(s$og)) sprintf("were obtained from the OpenGWAS database (dataset %s; ancestry %s%s; accessed %s)", s$og$id, s$ancestry, nstr(s), substr(s$og$downloaded_at %||% "", 1, 10))
                                else sprintf("were read from the file \"%s\" (ancestry declared as %s%s)", s$file, s$ancestry, nstr(s))
      txt$methods <- c(txt$methods, sprintf("Exposure %s and outcome %s. Genetic associations with the exposure %s, and associations with the outcome %s. %s",
        p$exName, p$outName, side(sx, "were"), side(sy, "were"),
        if (!is.null(sx$og) && !is.null(sy$og)) "Ancestry and the units of the effect sizes are as reported by OpenGWAS; the genome build was not verified automatically." else "Ancestry, genome build and the units of the effect sizes are as declared by the analyst when files were imported; they were not verified automatically."))
      L <- cx$ledger
      if (!is.null(L) && !is.null(sx$og) && is.null(L$clump)) {
        pr <- sx$og$params; po <- if (!is.null(sy$og)) sy$og$params else NULL
        txt$methods <- c(txt$methods, sprintf("Instruments were the genome-wide significant variants of the exposure dataset (P < %s), %s (r\u00b2 < %s within %s kb, 1000 Genomes %s reference panel; ieugwasr %s). %s Of %d instruments, %d %s absent from the outcome%s. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure, %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used.",
          format(pr$pval), if (identical(pr$where, "server_live")) "clumped live on the OpenGWAS server with ieugwasr::tophits (force_server = TRUE)" else if (ogw_stored_clumps(pr)) "as returned by ieugwasr::tophits from the stored, pre-computed clumps of OpenGWAS (default settings, preclumped = 1)" else "clumped on the OpenGWAS server with ieugwasr::tophits",
          format(pr$r2), format(pr$kb, big.mark = ","), pr$pop, sx$og$packages$ieugwasr %||% "unknown",
          if (!is.null(po)) sprintf("The same variants were extracted from the outcome dataset with TwoSampleMR::extract_outcome_data (%s; TwoSampleMR %s).", if (isTRUE(po$proxies)) sprintf("proxies allowed at r\u00b2 of at least %s, palindromic variants aligned by allele frequency only when the minor allele frequency was below %s", format(po$rsq), format(po$maf_threshold)) else "no proxies", sy$og$packages$TwoSampleMR %||% "unknown") else "The outcome associations were read from an imported file.",
          L$instruments, L$missing_in_outcome, if (L$missing_in_outcome == 1) "was" else "were", if (!is.null(L$by_proxy) && L$by_proxy > 0) sprintf(" (a further %d %s found only through proxy variants)", L$by_proxy, if (L$by_proxy == 1) "was" else "were") else "", L$flipped, L$dropped_palindromic, L$dropped_mismatch, L$kept))
      } else       if (!is.null(L) && !is.null(L$clump)) {
        C <- L$clump
        og_local <- !is.null(sx$og); po <- if (og_local && !is.null(sy$og)) sy$og$params else NULL
        if (og_local) txt$methods <- c(txt$methods, sprintf("Every variant of the exposure dataset with P < %s was downloaded from OpenGWAS without clumping (ieugwasr::tophits with clump = 0; ieugwasr %s) and clumped on the analyst's computer at download, as described next. Variants with alleles other than single bases, a missing or non-positive standard error, or a repeated ID were removed before clumping.%s",
          format(sx$og$params$pval), sx$og$packages$ieugwasr %||% "unknown",
          if (!is.null(po)) sprintf(" The independent instruments described next were looked up in the outcome dataset with TwoSampleMR::extract_outcome_data (%s; TwoSampleMR %s).", if (isTRUE(po$proxies)) sprintf("proxies allowed at r\u00b2 of at least %s, palindromic variants aligned by allele frequency only when the minor allele frequency was below %s", format(po$rsq), format(po$maf_threshold)) else "no proxies", sy$og$packages$TwoSampleMR %||% "unknown") else ""))
        txt$methods <- c(txt$methods, sprintf("Instruments were selected from the exposure summary statistics on the analyst's computer: the %d variants with P < %s were clumped with %s [Chang 2015] using the %s 1000 Genomes Phase 3 reference panel [1000 Genomes 2015] distributed by the MRC Integrative Epidemiology Unit (r\u00b2 < %s within %s kb, keeping the variant with the smallest P in each region%s). %d variants were removed as correlated with a stronger variant and %d were removed because they are not in the reference panel, leaving %d independent instruments, of which %d were absent from the outcome%s. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure (swapped or on the opposite strand), %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used. Every candidate variant, kept or removed, is listed with its index variant in the pair's instrument_selection.csv.",
          as.integer(L$below_threshold), format(L$threshold), C$plink %||% "PLINK 1.9", C$pop, format(C$r2), format(C$kb, big.mark = ","), if (identical(C$tie_break, "equal P ordered by |beta/SE|")) "; variants with equal P, such as P values too small to store, were ordered by the absolute value of beta/SE" else "", as.integer(C$clumped), as.integer(C$not_in_panel), as.integer(C$kept), as.integer(L$missing_in_outcome), if (!is.null(L$by_proxy) && L$by_proxy > 0) sprintf(" (a further %d %s found only through proxy variants)", L$by_proxy, if (L$by_proxy == 1) "was" else "were") else "", as.integer(L$flipped), as.integer(L$dropped_palindromic), as.integer(L$dropped_mismatch), as.integer(L$kept)))
      } else if (!is.null(L)) txt$methods <- c(txt$methods, sprintf("Instruments were the variants in the exposure file with P < %s. The exposure file was declared by the analyst to contain only independent (clumped) variants; MR Studio did not perform clumping. Of %d candidate variants, %d were absent from the outcome file. Exposure and outcome alleles were aligned with TwoSampleMR::harmonise_data (action 2: palindromic variants are aligned using allele frequency and dropped when the strand is ambiguous): %d variants had their alleles aligned to the exposure (swapped or on the opposite strand), %d palindromic variants were dropped and %d variants were dropped because the alleles did not match. %d instruments were used.", format(L$threshold), L$candidates, L$missing_in_outcome, L$flipped, L$dropped_palindromic, L$dropped_mismatch, L$kept))
    }
    fl <- if (isTRUE(st$min_f > 0)) sprintf(" Variants with a single-variant F statistic below %s were left out.", st$min_f) else ""
    txt$methods <- c(txt$methods, sprintf("The strength of the instruments was summarised by the single-variant F statistic, (beta_x / se_x)^2: mean %.1f, minimum %.1f, with %d of %d variants below 10.%s", cx$mean_f, cx$min_f, cx$n_weak, nrow(cx$kept), fl))
  }
  ivw_model <- c(mre = "multiplicative random effects", mre_floor = "random effects when the variants are over-dispersed and fixed effect otherwise", fixed = "a fixed-effect model")
  st1 <- ctxs[[1]]$settings
  desc <- vapply(used, function(m) switch(m,
    ivw_mre = sprintf("The primary estimate was the %s estimator with %s [Burgess 2013, Hemani 2018].", mname[["ivw_mre"]], ivw_model[[st1$ivw_model %||% "mre"]]),
    egger = "MR-Egger regression was used to test for directional pleiotropy through its intercept and to estimate the effect allowing for it; variants were oriented to a positive exposure effect and confidence intervals use a t distribution with k - 2 degrees of freedom [Bowden 2015].",
    wm = sprintf("The weighted median estimator, which is consistent when at least half of the weight comes from valid instruments, was computed with a bootstrapped standard error (%s replications) [Bowden 2016].", format(st1$nboot %||% 1000)),
    smed = sprintf("The simple median estimator, which is consistent when at least half of the variants are valid instruments, was computed with a bootstrapped standard error (%s replications) [Bowden 2016].", format(st1$nboot %||% 1000)),
    mode = sprintf("The weighted mode estimator, which is consistent when the largest group of variants with similar ratio estimates is valid, was computed with a bandwidth multiplier phi = %s and a bootstrapped standard error (%s replications) [Hartwig 2017].", format(st1$phi %||% 1), format(st1$nboot %||% 1000)),
    smode = sprintf("The simple mode estimator (the unweighted counterpart of the weighted mode) was computed with a bandwidth multiplier phi = %s and a bootstrapped standard error (%s replications) [Hartwig 2017].", format(st1$phi %||% 1), format(st1$nboot %||% 1000)),
    presso = sprintf("MR-PRESSO (global test of horizontal pleiotropy, outlier test and distortion test) was run with %s simulated datasets and a significance threshold of %s. Variants it flags are reported but were not removed from any analysis [Verbanck 2018].", format(st1$nbdist %||% 1000), format(st1$signif %||% 0.05)),
    radial = sprintf("Radial IVW was used to flag outlying variants through their contribution to Cochran's Q statistic (outlier level alpha = %s, %s weights); flagged variants are reported but were not removed [Bowden 2018].", format(st1$alpha %||% 0.05), c("1" = "first-order", "2" = "second-order", "3" = "modified second-order")[[as.character(st1$weights %||% "1")]]),
    steiger = paste0("The Steiger directionality test, run together with the IVW analysis, compared the variance the variants explain in the exposure with the variance they explain in the outcome, summed over variants (TwoSampleMR::mr_steiger; measurement reliability set to 1 unless stated otherwise). For a continuous trait the variance explained was computed from each variant's effect, standard error and sample size",
      if (steiger_any_binary(ctxs)) "; for a binary trait it was computed on the liability scale from the log odds ratio, the effect allele frequency, the numbers of cases and controls and the population prevalence (TwoSampleMR::get_r_from_lor, logit model) [Lee 2012]" else "",
      ". It gives no causal estimate [Hemani 2017]."),
    radegger = sprintf("Radial MR-Egger was used to flag outlying variants and to estimate the effect allowing for directional pleiotropy (outlier level alpha = %s, %s weights); flagged variants are reported but were not removed [Bowden 2018].", format(st1$alpha %||% 0.05), c("1" = "first-order", "2" = "second-order", "3" = "modified second-order")[[as.character(st1$weights %||% "1")]]), ""), "")
  txt$methods <- c(txt$methods, unname(desc[nzchar(desc)]), sprintf("Confidence intervals are %s intervals. All estimates come from TwoSampleMR (commit %s)%s in %s. Bootstrap-based methods use a fixed seed recorded in each script, and every pair was re-run from its saved scripts in a fresh R session; the replay check was %s.",
    .pct(ctxs[[1]]$ce$ci_level[1] %||% 0.95), substr(ADAPTER_IVW_MRE_V11$packages$TwoSampleMR, 1, 10), if (any(c("presso", "radial", "radegger") %in% used)) sprintf(", with %s%s%s", if ("presso" %in% used) sprintf("MRPRESSO (commit %s)", substr(PIN_MRPRESSO, 1, 10)) else "", if ("presso" %in% used && any(c("radial", "radegger") %in% used)) " and " else "", if (any(c("radial", "radegger") %in% used)) sprintf("RadialMR (commit %s)", substr(PIN_RADIALMR, 1, 10)) else "") else "", R.version.string, paste(unique(vapply(ctxs, function(x) x$pay$replay$status %||% "not run", "")), collapse = " / ")))
  for (cx in ctxs) {
    ce <- cx$ce; if (!nrow(ce)) next; sc <- ce$scale[1]; ratio <- sc == "log_or"
    unit <- if (ratio) sprintf("odds of %s per unit of %s", cx$p$outName, cx$p$exName) else sprintf("units of %s per unit of %s", cx$p$outName, cx$p$exName)
    lines <- vapply(seq_len(nrow(ce)), function(i) sprintf("%s: %s", .short[[ce$method[i]]], .ci_phrase(ce[i, ], sc)), "")
    ic <- cx$est[cx$est$estimand == "pleiotropy_intercept", , drop = FALSE]; dg <- Filter(function(d) d$method == "ivw_mre" && d$name == "cochran_q", cx$pay$diagnostics)
    same <- sum(sign(ce$estimate) == sign(ce$estimate[ce$method == "ivw_mre"][1]), na.rm = TRUE)
    txt$results <- c(txt$results, sprintf("%s -> %s (%d instruments, effect expressed as %s). %s.%s%s %d of %d methods agree on the direction of the effect.", cx$p$exName, cx$p$outName, ce$n_ivs[1], unit, paste(lines, collapse = "; "),
      if (length(dg)) sprintf(" Cochran's Q for the IVW estimate was %.1f (P = %s).", dg[[1]]$value, .fp(dg[[1]]$p)) else "",
      if (nrow(ic)) sprintf(" The MR-Egger intercept was %s (P = %s)%s.", .fnum(ic$estimate, 4), .fp(ic$p), if (ic$p < 0.05) ", suggesting directional pleiotropy" else ", giving no evidence of directional pleiotropy") else "", same, nrow(ce)))
    sf <- steiger_facts(cx)
    if (!is.null(sf)) txt$results <- c(txt$results, sprintf("Steiger directionality test: the variants explained %s of the variance in %s and %s of the variance in %s (R2), %s (P = %s); %d of %d variants explained more variance in the outcome than in the exposure.",
      .fnum(sf$r2x, 3), cx$p$exName, .fnum(sf$r2y, 3), cx$p$outName, if (sf$ok) "so the exposure explains more, as expected if it acts on the outcome" else "so the outcome explains more, which does not match the assumed direction", .fp(sf$p), as.integer(sf$wrong), as.integer(sf$n)))
    for (dm in intersect(c("presso", "radial", "radegger"), cx$ce$method)) {
      dg <- Filter(function(d) d$method == dm, cx$pay$diagnostics); od <- Filter(function(d) grepl("_outliers$", d$name), dg)
      ids <- if (length(od) && !is.null(od[[1]]$note) && !is.na(od[[1]]$note) && nzchar(od[[1]]$note)) strsplit(od[[1]]$note, ";")[[1]] else character(0)
      gl <- Filter(function(d) d$name == "presso_global_test", dg); ql <- Filter(function(d) d$name == "radial_q", dg)
      txt$results <- c(txt$results, sprintf("%s flagged %s of %d variants%s%s. No variant was removed from the primary analysis.", .short[[dm]], if (length(ids)) length(ids) else "none", nrow(cx$kept), if (length(ids)) paste0(" (", paste(ids, collapse = ", "), ")") else "",
        if (length(gl)) sprintf("; the global test P was %s%s", if (!is.null(gl[[1]]$note) && !is.na(gl[[1]]$note)) "below " else "", .fp(gl[[1]]$p)) else if (length(ql)) sprintf("; the radial Q statistic was %.1f (P = %s)", ql[[1]]$value, .fp(ql[[1]]$p)) else ""))
    }
  }
  lim <- c("The three core MR assumptions (relevance, independence from confounders, and no effect of the variants on the outcome except through the exposure) cannot be proven from summary data. The median and mode estimators relax the third assumption only under their own conditions (a majority or plurality of valid instruments).")
  if (anysim) lim <- c(lim, "SIMULATED DATA: at least one pair uses simulated demonstration traits. These results illustrate the software and must not be presented as scientific findings.")
  anyog <- any(vapply(ctxs, function(cx) !cx$simulated && !is.null(cx$src_x$og), TRUE)); anyfile <- any(vapply(ctxs, function(cx) !cx$simulated && (is.null(cx$src_x$og) || is.null(cx$src_y$og)), TRUE))
  if (anyfile) lim <- c(lim, "For imported files, clumping was not performed by MR Studio: the independence of the instruments rests on the analyst's statement at import. Ancestry and genome build are declared, not verified, and the exposure and outcome files must come from the same genome build.")
  if (anyog) lim <- c(lim, "Instruments from OpenGWAS were clumped against 1000 Genomes reference panels (on its server or, where stated, on the analyst's computer with the same panel), so linkage disequilibrium was estimated in a reference sample and not in the study population. Proxy variants, where used, stand in for missing instruments and are listed in the harmonised data table.")
  lim <- c(lim, overlap_limitations(ctxs),
           steiger_limitations(ctxs),
           if (any(vapply(ctxs, function(x) x$n_weak > 0, TRUE))) "Some variants have an F statistic below 10, which can bias estimates toward the null in two-sample MR." else NULL,
           if (any(c("presso", "radial", "radegger") %in% used)) "Variants flagged by MR-PRESSO or radial MR were reported but not removed from the primary analysis. Detectors flag variants by a statistical rule and cannot tell horizontal pleiotropy from valid heterogeneity; any analysis without flagged variants is a post-hoc sensitivity analysis and must be reported as one [Burgess 2023]." else "Outlier detection (MR-PRESSO, radial MR) was not run in this study. Pleiotropy-robust estimators were compared with the IVW estimate; agreement in direction does not prove the assumptions hold.",
           "No correction for multiple testing was applied across pairs.")
  txt$limitations <- lim
  src <- unique(unlist(lapply(ctxs, function(cx) c(if (!cx$simulated) cx$src_x$file, if (!cx$simulated) cx$src_y$file))))
  txt$availability <- c(if (length(src)) sprintf("Summary statistics were %s: %s.%s", if (all(grepl("^OpenGWAS:", src))) "obtained from the public OpenGWAS database (api.opengwas.io), dataset IDs" else "read from the following files", paste(sub("^OpenGWAS:", "", src), collapse = "; "), if (all(grepl("^OpenGWAS:", src))) "" else " [Add the public source, accession number or DOI of each file.]") else "No external data were used: the traits are simulated by MR Studio.",
                        "The analysis scripts, the harmonised data tables used in every calculation, all result tables and the replay check are provided in the reproducibility package of this study.")
  fi <- 0L
  for (cx in ctxs) { ce <- cx$ce; if (!nrow(ce)) next; fi <- fi + 1L
    txt$legends <- c(txt$legends, sprintf("Figure %d. Forest plot of the Mendelian randomization estimates of %s on %s. Points are estimates and bars are %s confidence intervals for each method (%s). %d instruments were used.", 2L * fi - 1L, cx$p$exName, cx$p$outName, .pct(ce$ci_level[1]), paste(.short[ce$method], collapse = ", "), ce$n_ivs[1]),
                    sprintf("Figure %d. Scatter plot of variant effects. Each point is one variant: its effect on %s (x) and on %s (y), with %s intervals for the outcome effect. Variants are oriented so that the exposure effect is positive. Lines show the slope of each method; the MR-Egger line includes its intercept.", 2L * fi, cx$p$exName, cx$p$outName, "95%")) }
  txt
}

# ---------- the package ----------
`%||%` <- function(a, b) if (is.null(a)) b else a
package_dir <- function(h) { root <- file.path(h$dir, "05_publication"); dir.create(root, showWarnings = FALSE); n <- length(list.dirs(root, recursive = FALSE)) + 1L; d <- file.path(root, sprintf("package_%03d", n)); dir.create(d); d }

build_package <- function(h, items, width_mm = 180, pub = NULL, pairs = NULL) {
  pairs <- pairs %||% (ui_state_get(h$con, "pairs") %||% list()); pairs <- Filter(function(p) !is.null(p$folder) && !is.null(pair_results(h$dir, p$folder)), pairs)
  if (!length(pairs)) stop("There are no saved results to put in a package. Run the analysis first.")
  ctxs <- lapply(pairs, function(p) pair_context(h, p)); pk <- package_dir(h); made <- character(); note <- character()
  add <- function(f) made <<- c(made, substring(f, nchar(h$dir) + 2L))
  txt <- report_text(ctxs, basename(h$dir)); used <- unique(c("hemani2018", "burgess2023", "skrivankova2021", unlist(METHOD_REFS[unlist(lapply(ctxs, function(x) x$methods))]), if (steiger_any_binary(ctxs)) "lee2012", if (any(vapply(ctxs, function(cx) !is.null(cx$ledger$clump), TRUE))) c("chang2015", "kgp2015")))
  figs <- list(); wmm <- as.numeric(width_mm); wmm <- if (is.finite(wmm) && wmm >= 60) wmm else 180
  if (any(c("figs", "report", "text") %in% items)) {
    dir.create(file.path(pk, "figures"), showWarnings = FALSE)
    for (i in seq_along(ctxs)) { cx <- ctxs[[i]]; ce <- cx$ce; if (!nrow(ce)) next; tag <- sprintf("P%d", i)
      est <- transform(ce, label = .short[method], col = .MCOL[method], unit = if (ce$scale[1] == "log_or") "log odds ratio" else "effect on outcome per unit of exposure", ratio = FALSE)
      fh <- max(40, 22 + 8 * nrow(est)); ff <- .fig_write(file.path(pk, "figures", paste0(tag, "_forest")), wmm, fh, function() .draw_forest(est))
      ceS <- ce[ce$method != "radegger", , drop = FALSE]   # the radial Egger intercept is on a weighted scale, so it cannot be drawn as a line on this plot
      ln <- lapply(seq_len(nrow(ceS)), function(j) { m <- ceS$method[j]; ic <- cx$est$estimate[cx$est$estimand == "pleiotropy_intercept"]; list(name = .short[[m]], col = .MCOL[[m]], lty = if (m == "ivw_mre") 1 else if (m == "egger") 2 else 3, a = if (m == "egger" && length(ic)) ic[1] else 0, b = ceS$estimate[j]) })
      sf <- .fig_write(file.path(pk, "figures", paste0(tag, "_scatter")), wmm, round(wmm * 0.75), function() .draw_scatter(cx$kept, ln, paste("Effect on", cx$p$exName), paste("Effect on", cx$p$outName)))
      figs[[i]] <- list(forest = ff, scatter = sf, fh = fh, sh = round(wmm * 0.75)); for (f in c(ff, sf)) add(f) }
    if (!length(Filter(Negate(is.null), figs))) note <- c(note, "No figures could be drawn on this computer.")
    if (!isTRUE(unname(capabilities("cairo")))) note <- c(note, "SVG and PNG figures need the cairo graphics library, which this R does not have. The PDF figures were written.")
  }
  if ("scripts" %in% items) { dst <- file.path(pk, "scripts_and_tables"); dir.create(dst); for (p in pairs) { file.copy(file.path(h$dir, "03_pairs", p$folder), dst, recursive = TRUE) }; add(dst) }
  if ("tables" %in% items) {
    est <- do.call(rbind, lapply(seq_along(ctxs), function(i) { e <- ctxs[[i]]$est; data.frame(pair = i, exposure = ctxs[[i]]$p$exName, outcome = ctxs[[i]]$p$outName, e, check.names = FALSE, stringsAsFactors = FALSE) }))
    sheets <- list(Estimates = est, Settings = do.call(rbind, lapply(seq_along(ctxs), function(i) { s <- ctxs[[i]]$settings; data.frame(pair = i, setting = names(s), value = vapply(s, function(v) paste(v, collapse = ","), ""), stringsAsFactors = FALSE) })))
    for (i in seq_along(ctxs)) sheets[[sprintf("Pair %d variants", i)]] <- transform(ctxs[[i]]$vars, ratio = by / bx, F = (bx / sx)^2)
    f <- file.path(pk, "supplement_tables.xlsx"); xlsx_write(f, sheets); add(f)
  }
  if ("refs" %in% items) { f1 <- file.path(pk, "references.bib"); writeLines(refs_bibtex(used), f1, useBytes = TRUE); f2 <- file.path(pk, "references.ris"); writeLines(refs_ris(used), f2, useBytes = TRUE); add(f1); add(f2) }
  blocks_text <- function(with_figs = FALSE) {
    b <- list(); H <- function(t, l = 1) b[[length(b) + 1L]] <<- list(type = if (l == 1) "h1" else "h2", text = t); P <- function(t) for (x in t) b[[length(b) + 1L]] <<- list(type = "p", text = x)
    if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) b[[1]] <- list(type = "note", text = "SIMULATED DATA: not for publication. At least one pair in this document uses simulated demonstration traits.")
    H("Methods"); P(txt$methods); H("Results"); P(txt$results)
    if (with_figs) for (i in seq_along(figs)) { f <- figs[[i]]; if (is.null(f)) next; png <- grep("forest\\.png$", f$forest, value = TRUE); sp <- grep("scatter\\.png$", f$scatter, value = TRUE)
      if (length(png)) b[[length(b) + 1L]] <- list(type = "image", file = png, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * f$fh / wmm, alt = txt$legends[2L * i - 1L]); P(txt$legends[2L * i - 1L])
      if (length(sp)) b[[length(b) + 1L]] <- list(type = "image", file = sp, width_mm = min(wmm, 160), height_mm = min(wmm, 160) * f$sh / wmm, alt = txt$legends[2L * i]); P(txt$legends[2L * i]) }
    H("Limitations"); P(txt$limitations); H("Data availability"); P(txt$availability); H("References"); P(vapply(used, function(i) ref_text(REFS[[i]]), "")); b }
  if ("text" %in% items) { f <- file.path(pk, "draft_methods_results_limitations.docx"); docx_write(f, c(list(list(type = "title", text = "Draft text: Methods, Results, Limitations")), list(list(type = "note", text = "Draft generated from the saved results and settings. Check every statement and edit it before use.")), blocks_text(FALSE))); add(f) }
  if ("report" %in% items) {
    est <- do.call(rbind, lapply(seq_along(ctxs), function(i) { ce <- ctxs[[i]]$ce; if (!nrow(ce)) return(NULL); data.frame(Pair = i, Exposure = ctxs[[i]]$p$exName, Outcome = ctxs[[i]]$p$outName, Method = ce$name, Estimate = .fnum(ce$estimate), `CI lower` = .fnum(ce$ci_lower), `CI upper` = .fnum(ce$ci_upper), P = .fp(ce$p), Instruments = ce$n_ivs, check.names = FALSE, stringsAsFactors = FALSE) }))
    bl <- c(list(list(type = "title", text = sprintf("Mendelian randomization report: %s", basename(h$dir))), list(type = "note", text = sprintf("Generated by MR Studio on %s. A draft for the authors to check.", format(Sys.Date())))), list(list(type = "h1", text = "Summary of estimates"), list(type = "table", data = est)), blocks_text(TRUE))
    f <- file.path(pk, "report.docx"); docx_write(f, bl); add(f)
    esc <- function(x) .xe(x); svgs <- function(files) { s <- grep("\\.svg$", files, value = TRUE); if (length(s)) paste(sub("^<\\?xml[^>]*\\?>\\s*", "", readLines(s[1], warn = FALSE, encoding = "UTF-8")), collapse = "\n") else "<p><em>Figure not available on this computer.</em></p>" }
    sect <- function(t, ps) sprintf("<h2>%s</h2>%s", esc(t), paste(sprintf("<p>%s</p>", esc(ps)), collapse = "\n"))
    figh <- paste(vapply(seq_along(figs), function(i) { f <- figs[[i]]; if (is.null(f)) "" else sprintf("<figure>%s<figcaption>%s</figcaption></figure><figure>%s<figcaption>%s</figcaption></figure>", svgs(f$forest), esc(txt$legends[2L * i - 1L]), svgs(f$scatter), esc(txt$legends[2L * i])) }, ""), collapse = "\n")
    tab <- paste0("<table><thead><tr>", paste(sprintf("<th>%s</th>", esc(names(est))), collapse = ""), "</tr></thead><tbody>", paste(vapply(seq_len(nrow(est)), function(i) paste0("<tr>", paste(sprintf("<td>%s</td>", esc(unlist(est[i, ]))), collapse = ""), "</tr>"), ""), collapse = ""), "</tbody></table>")
    sim <- if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) "<div class=\"warn\">SIMULATED DATA: not for publication. At least one pair uses simulated demonstration traits.</div>" else ""
    html <- paste0("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>MR report: ", esc(basename(h$dir)), "</title><style>body{font:16px/1.55 Georgia,'Times New Roman',serif;max-width:46rem;margin:2.5rem auto;padding:0 1.2rem;color:#1b1b1b}h1{font-size:1.9rem}h2{margin-top:2.2rem;font-size:1.3rem;border-bottom:1px solid #ccc;padding-bottom:.2rem}table{border-collapse:collapse;width:100%;font:14px/1.4 system-ui,sans-serif}th,td{border:1px solid #bbb;padding:.35rem .5rem;text-align:left}th{background:#f1f1f1}figure{margin:1.5rem 0}figure svg{max-width:100%;height:auto}figcaption{font-size:.9rem;color:#444}.warn{border:2px solid #b5541c;background:#fff4ec;padding:.7rem 1rem;margin:1rem 0;font:600 15px system-ui,sans-serif}.note{color:#555;font-style:italic}@media print{body{margin:0;max-width:none}h2{break-after:avoid}figure,table{break-inside:avoid}}</style></head><body><h1>Mendelian randomization report: ", esc(basename(h$dir)), "</h1><p class=\"note\">Generated by MR Studio on ", format(Sys.Date()), ". A draft for the authors to check. To make a PDF, open this file in a browser and print it.</p>", sim, "<h2>Summary of estimates</h2>", tab, sect("Methods", txt$methods), sect("Results", txt$results), "<h2>Figures</h2>", figh, sect("Limitations", txt$limitations), sect("Data availability", txt$availability), sprintf("<h2>References</h2><ol>%s</ol>", paste(sprintf("<li>%s</li>", esc(vapply(used, function(i) ref_text(REFS[[i]]), ""))), collapse = "")), "</body></html>")
    f <- file.path(pk, "report.html"); writeLines(html, f, useBytes = TRUE); add(f)
  }
  if ("capsule" %in% items) {
    fs <- list.files(file.path(h$dir, "03_pairs"), recursive = TRUE, full.names = TRUE); rel <- substring(fs, nchar(h$dir) + 2L)
    ent <- setNames(lapply(fs, .file_raw), rel)
    ent[["SESSION.txt"]] <- .txt(c(sprintf("MR Studio reproducibility capsule, written %s", format(Sys.time())), R.version.string, sprintf("Platform: %s", R.version$platform), sprintf("TwoSampleMR commit pinned by the adapters: %s", ADAPTER_IVW_MRE_V11$packages$TwoSampleMR), "", "Each pair folder holds its scripts (02_scripts), data (01_data), results (04_tables), console output and replay check. Open a pair folder and run 02_scripts/00_run_all.R from it."))
    f <- file.path(pk, "reproducibility_capsule.zip"); zip_write(f, ent); add(f)
  }
  strobe <- NULL
  if ("strobe" %in% items) {
    rows <- strobe_rows(ctxs, txt, tryCatch(ui_state_get(h$con, "plan") %||% list(), error = function(e) list()), items)
    f <- file.path(pk, "strobe_mr_checklist.docx"); n <- strobe_write(f, rows, basename(h$dir)); add(f)
    strobe <- list(filled = n[["FILLED"]], partly = n[["PARTLY"]], author = n[["AUTHOR"]])
  }
  readme <- c("MR Studio publication package", sprintf("Study: %s", basename(h$dir)), sprintf("Built: %s", format(Sys.time())), "", "Everything here is a draft for the authors to check. Numbers come from the saved result tables; texts are assembled from the recorded settings.",
              if (any(vapply(ctxs, `[[`, TRUE, "simulated"))) "At least one pair uses SIMULATED demonstration traits: do not report these results as findings." else NULL, "", "Files:", made, if (length(note)) c("", "Notes:", note) else NULL)
  rf <- file.path(pk, "README.txt"); writeLines(readme, rf, useBytes = TRUE)
  list(dir = pk, folder = basename(pk), files = made, notes = note, simulated = any(vapply(ctxs, `[[`, TRUE, "simulated")), strobe = strobe)
}
