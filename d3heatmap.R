
d3heatmap <<- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                       "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                                     w), k_row, k_col, symm = FALSE, revC, scale = c("none", 
                                                                                                                                                                                     "row", "column"), na.rm = TRUE, labRow = rownames(x), 
                    labCol = colnames(x), cexRow, cexCol, digits = 3L, cellnote, 
                    cellnote_scale = FALSE, theme = NULL, colors = "RdYlBu", 
                    width = NULL, height = NULL, xaxis_height = 80, yaxis_width = 120, 
                    xaxis_font_size = NULL, yaxis_font_size = NULL, brush_color = "#0000FF", 
                    show_grid = TRUE, anim_duration = 500, ...) 
{
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) 
    stop("x must be a matrix")
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  if(!is.null(labRow)) rownames(x) <- labRow else rownames(x) <- paste(1:nrow(x))
  if(!is.null(labCol)) colnames(x) <- labCol else colnames(x) <- paste(1:ncol(x))
  
  if (!missing(cexRow)) {
    if (is.numeric(cexRow)) {
      xaxis_font_size <- cexRow * 14
    }
    else {
      warning("cexRow is not numeric. It is ignored")
    }
  }
  if (!missing(cexCol)) {
    if (is.numeric(cexCol)) {
      yaxis_font_size <- cexCol * 14
    }
    else {
      warning("cexCol is not numeric. It is ignored")
    }
  }
  dendrogram <- match.arg(dendrogram)
  if (missing(Rowv)) {
    Rowv <- dendrogram %in% c("both", "row")
  }
  if (missing(Colv)) {
    Colv <- dendrogram %in% c("both", "column")
  }
  if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
  }
  if (is.numeric(Rowv)) {
    Rowv <- reorderfun(as.dendrogram(hclustfun(distfun(x))), 
                       Rowv)
  }
  if (d3heatmap:::is.dendrogram(Rowv)) {
    Rowv <- rev(Rowv)
    rowInd <- order.dendrogram(Rowv)
    if (nr != length(rowInd)) 
      stop("Row dendrogram is the wrong size")
  }
  else {
    if (!is.null(Rowv) && !is.na(Rowv) && !identical(Rowv, 
                                                     FALSE)) 
      warning("Invalid value for Rowv, ignoring")
    Rowv <- NULL
    rowInd <- 1:nr
  }
  if (identical(Colv, "Rowv")) {
    Colv <- Rowv
  }
  if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
  }
  if (is.numeric(Colv)) {
    Colv <- reorderfun(as.dendrogram(hclustfun(distfun(t(x)))), 
                       Colv)
  }
  if (d3heatmap:::is.dendrogram(Colv)) {
    colInd <- order.dendrogram(Colv)
    if (nc != length(colInd)) 
      stop("Col dendrogram is the wrong size")
  }
  else {
    if (!is.null(Colv) && !is.na(Colv) && !identical(Colv, 
                                                     FALSE)) 
      warning("Invalid value for Colv, ignoring")
    Colv <- NULL
    colInd <- 1:nc
  }
  if (missing(revC)) {
    if (symm) {
      revC <- TRUE
    }
    else if (d3heatmap:::is.dendrogram(Colv) & d3heatmap:::is.dendrogram(Rowv) & 
             identical(Rowv, rev(Colv))) {
      revC <- TRUE
    }
    else {
      revC <- FALSE
    }
  }
  if (revC) {
    Colv <- rev(Colv)
    colInd <- rev(colInd)
  }
  x <- x[rowInd, colInd]
  if (!missing(cellnote)) 
    cellnote <- cellnote[rowInd, colInd]
  if (!missing(k_row) | !missing(k_col)) 
    dendextend::assign_dendextend_options()
  if (d3heatmap:::is.dendrogram(Rowv) & !missing(k_row)) {
    Rowv <- dendextend::color_branches(Rowv, k = k_row)
  }
  if (d3heatmap:::is.dendrogram(Colv) & !missing(k_col)) {
    Colv <- dendextend::color_branches(Colv, k = k_col)
  }
  rowDend <- if (d3heatmap:::is.dendrogram(Rowv)) 
    d3heatmap:::dendToTree(Rowv)
  else NULL
  colDend <- if (d3heatmap:::is.dendrogram(Colv)) 
    d3heatmap:::dendToTree(Colv)
  else NULL
  scale <- match.arg(scale)
  if (!cellnote_scale) 
    x_unscaled <- x
  if (scale == "row") {
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
    x <- sweep(x, 1, apply(x, 1, sd, na.rm = na.rm), "/")
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
    x <- sweep(x, 2, apply(x, 2, sd, na.rm = na.rm), "/")
  }
  if (missing(cellnote)) {
    if (cellnote_scale) {
      cellnote <- round(x, digits = digits)
    }
    else {
      cellnote <- round(x_unscaled, digits = digits)
    }
  }
  if (is.null(dim(cellnote))) {
    if (length(cellnote) != nr * nc) {
      stop("Incorrect number of cellnote values")
    }
    dim(cellnote) <- dim(x)
  }
  if (!identical(dim(x), dim(cellnote))) {
    stop("cellnote matrix must have same dimensions as x")
  }
  mtx <- list(data = as.character(t(cellnote)), dim = dim(x), 
              rows = rownames(x), cols = colnames(x))
  if (is.factor(x)) {
    colors <- scales::col_factor(colors, x, na.color = "gray")
  }
  else {
    rng <- range(x, na.rm = TRUE)
    if (scale %in% c("row", "column")) {
      rng <- c(max(abs(rng)), -max(abs(rng)))
    }
    colors <- scales::col_numeric(colors, rng, na.color = "gray")
  }
  imgUri <- d3heatmap:::encodeAsPNG(t(x), colors)
  options <- NULL
  options <- c(options, list(xaxis_height = xaxis_height, 
                             yaxis_width = yaxis_width, xaxis_font_size = xaxis_font_size, 
                             yaxis_font_size = yaxis_font_size, brush_color = brush_color, 
                             show_grid = show_grid, anim_duration = anim_duration))
  if (is.null(rowDend)) {
    c(options, list(yclust_width = 0))
  }
  if (is.null(colDend)) {
    c(options, list(xclust_height = 0))
  }
  payload <- list(rows = rowDend, cols = colDend, matrix = mtx, 
                  image = imgUri, theme = theme, options = options)
  htmlwidgets::createWidget(name = "d3heatmap", payload, width = width, 
                            height = height, package = "d3heatmap", sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE))
}


