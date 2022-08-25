colour_to_integer <- function(x) {
  # this function vs nara::colour_to_integer(x)
  # for x of length 256 (typical colormap):
  #   this takes 400ms, nara takes 200ms
  #   this takes 36KB RAM, nara takes 1KB
  # for x of length 1E7 (a whole matrix)
  #   this function takes 1.17s, nara takes 400ms
  #   this function uses 1.3GB RAM, nara uses 38MB (!!!)
  # Gives the same output, but since this is R and not C it is less optimized
  #
  # However typical colormaps can be of length 256, so we don't really care

  if (is.numeric(x)) {
    return(x)
  }

  # grDevices::col2rgb takes the RAM
  rgba <- grDevices::col2rgb(x, alpha = TRUE)
  col_as_int <- rgba[1,] + 256L*rgba[2,] + 256L*256L*rgba[3,] + 256*256*256*rgba[4,]

  my_ints <- as.integer(
    ifelse(
      rgba[4,] > 127L,
      -2^32 + col_as_int,
      col_as_int
    )
  )
  my_ints
}


mat_to_nativeRaster <- function(x, colormap, rangex = NULL)  {
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
  cols_to_ints <- colour_to_integer(colormap)
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

# Precompute the default colormap to speed up plotting
COLORMAP_VIRIDIS_256_A_m1 <- colour_to_integer(viridisLite::viridis(256L, direction = -1, option = "A"))
