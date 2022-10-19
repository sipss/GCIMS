mat_to_nativeRaster <- function(x, colormap, rangex = NULL)  {
  require_pkgs("farver")
  if (length(x) == 0) {
    return(
      structure(
        numeric(0L),
        dim = c(0L, 0L),
        class = "nativeRaster",
        channels = 4L
      )
    )
  }
  if (is.numeric(colormap)) {
    cols_to_ints <- colormap
  } else {
    cols_to_ints <- farver::encode_native(colormap)
  }
  if (is.null(rangex)) {
    rangex <- range(x)
  }
  breaks <- seq(from = rangex[1], to = rangex[2], length.out = length(colormap))
  xdim <- dim(x)
  rev_cols <- seq.int(ncol(x), 1L, by = -1L)
  x <- x[, rev_cols]
  x <- findInterval(x, breaks, rightmost.closed = TRUE)
  x <- cols_to_ints[x]
  x <- matrix(x, nrow = xdim[1], ncol = xdim[2], byrow = FALSE)

  structure(
    x,
    dim = c(xdim[2], xdim[1]),
    class = "nativeRaster",
    channels = 4L
  )
}
