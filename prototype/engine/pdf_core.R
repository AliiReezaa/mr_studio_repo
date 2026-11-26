# Round 46: the report as PDF, written by R itself from the same blocks as the Word report (no LaTeX, browser or extra program).
# A4 portrait; text is word-wrapped with the device's own font metrics; figures are the PNG files of the package; tables are
# drawn with wrapped cells. cairo_pdf is used where available so that Unicode (beta, squared, arrows) prints correctly.
.pdf_page <- list(w = 8.27, h = 11.69, ml = 0.9, mr = 0.9, mt = 0.85, mb = 0.85)
.pdf_style <- list(title = list(fs = 17, face = 2, after = 10), h1 = list(fs = 13.5, face = 2, before = 10, after = 5), h2 = list(fs = 11.5, face = 2, before = 8, after = 4),
                   p = list(fs = 10, face = 1, after = 6), note = list(fs = 9, face = 3, after = 6, col = "grey30"))

pdf_write <- function(path, blocks) {
  P <- .pdf_page; W <- P$w - P$ml - P$mr; y <- 0; page <- 0L
  open <- if (isTRUE(unname(capabilities("cairo")))) function() grDevices::cairo_pdf(path, width = P$w, height = P$h, onefile = TRUE, family = "sans")
          else function() grDevices::pdf(path, width = P$w, height = P$h, family = "Helvetica", onefile = TRUE)
  if (exists(".fig_metric_reset")) .fig_metric_reset()   # Round 109a: no cached text height from a closed device (report_core.R)
  open(); on.exit(grDevices::dev.off())
  newpage <- function() { grid::grid.newpage(); page <<- page + 1L; y <<- P$mt
    grid::grid.text(as.character(page), x = grid::unit(P$w / 2, "inches"), y = grid::unit(P$mb / 2, "inches"), gp = grid::gpar(fontsize = 8, col = "grey45")) }
  room <- function(h) { if (page == 0L || y + h > P$h - P$mb) newpage() }
  gp <- function(st) grid::gpar(fontsize = st$fs, fontface = st$face, col = st$col %||% "black")
  wrap <- function(txt, st, width) {
    words <- strsplit(gsub("\\s+", " ", trimws(txt)), " ", fixed = TRUE)[[1]]; if (!length(words)) return(character())
    grid::pushViewport(grid::viewport(gp = gp(st))); on.exit(grid::popViewport())
    ww <- grid::convertWidth(grid::stringWidth(words), "inches", valueOnly = TRUE); sp <- grid::convertWidth(grid::stringWidth(" "), "inches", valueOnly = TRUE)
    lines <- character(); cur <- character(); cw <- 0
    for (i in seq_along(words)) { add <- if (length(cur)) sp + ww[i] else ww[i]
      if (length(cur) && cw + add > width) { lines <- c(lines, paste(cur, collapse = " ")); cur <- words[i]; cw <- ww[i] } else { cur <- c(cur, words[i]); cw <- cw + add } }
    c(lines, paste(cur, collapse = " "))
  }
  lh <- function(st) st$fs * 1.38 / 72
  text_block <- function(txt, st) {
    if (!is.null(st$before) && y > P$mt + 0.01) y <<- y + st$before / 72
    for (ln in wrap(txt, st, W)) { room(lh(st)); y <<- y + lh(st)
      grid::grid.text(ln, x = grid::unit(P$ml, "inches"), y = grid::unit(P$h - y + lh(st) * 0.28, "inches"), just = c("left", "bottom"), gp = gp(st)) }
    y <<- y + (st$after %||% 0) / 72
  }
  table_block <- function(d) {
    st <- list(fs = 7.5, face = 1); sh <- list(fs = 7.5, face = 2); d <- as.data.frame(lapply(d, as.character), stringsAsFactors = FALSE, check.names = FALSE)
    wts <- vapply(names(d), function(n) max(nchar(n) + 2, min(22, max(nchar(d[[n]]), na.rm = TRUE))), 1); cw <- W * wts / sum(wts); pad <- 0.04
    row <- function(vals, s, line_above) {
      cells <- lapply(seq_along(vals), function(j) wrap(ifelse(is.na(vals[j]), "", vals[j]), s, cw[j] - 2 * pad)); h <- max(lengths(cells), 1) * lh(s) + 2 * pad
      room(h); top <- y; x0 <- P$ml
      if (line_above) grid::grid.lines(grid::unit(c(P$ml, P$ml + W), "inches"), grid::unit(P$h - c(top, top), "inches"), gp = grid::gpar(col = "grey60", lwd = 0.5))
      for (j in seq_along(cells)) { for (k in seq_along(cells[[j]])) grid::grid.text(cells[[j]][k], x = grid::unit(x0 + pad, "inches"), y = grid::unit(P$h - top - pad - k * lh(s) + lh(s) * 0.28, "inches"), just = c("left", "bottom"), gp = gp(s)); x0 <- x0 + cw[j] }
      y <<- top + h }
    row(names(d), sh, TRUE); for (i in seq_len(nrow(d))) row(unlist(d[i, ]), st, TRUE)
    grid::grid.lines(grid::unit(c(P$ml, P$ml + W), "inches"), grid::unit(P$h - c(y, y), "inches"), gp = grid::gpar(col = "grey60", lwd = 0.5)); y <<- y + 8 / 72
  }
  image_block <- function(b) {
    if (!file.exists(b$file) || !requireNamespace("png", quietly = TRUE)) { text_block("[Figure not available in this PDF; see the figures folder.]", .pdf_style$note); return(invisible()) }
    img <- png::readPNG(b$file); wi <- min(W, (b$width_mm %||% 160) / 25.4); hi <- wi * dim(img)[1] / dim(img)[2]
    if (hi > P$h - P$mt - P$mb) { hi <- P$h - P$mt - P$mb; wi <- hi * dim(img)[2] / dim(img)[1] }
    room(hi + 0.05); grid::grid.raster(img, x = grid::unit(P$ml + (W - wi) / 2, "inches"), y = grid::unit(P$h - y, "inches"), width = grid::unit(wi, "inches"), height = grid::unit(hi, "inches"), just = c("left", "top"))
    y <<- y + hi + 6 / 72
  }
  newpage()   # start the first page before any text is measured, so no blank page comes first
  for (b in blocks) switch(b$type,
    title = text_block(b$text, .pdf_style$title), h1 = text_block(b$text, .pdf_style$h1), h2 = text_block(b$text, .pdf_style$h2),
    p = text_block(b$text, .pdf_style$p), note = text_block(b$text, .pdf_style$note), table = table_block(b$data), image = image_block(b), NULL)
  invisible(page)   # number of pages written
}
