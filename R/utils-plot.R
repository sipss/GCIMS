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


#' Resolve an `intensity_range` argument (shared by `plot,GCIMSSample-method`
#' and `plot,GCIMSDataset-method`) into `c(min, max)`
#'
#' @param intensity_range One of `"global"`, `"ranged"`, a numeric vector of
#' length 2, or a length-2 vector/list whose elements are independently a
#' number, `"global"` or `"ranged"`.
#' @param global_range A zero-argument function returning the full/uncropped
#' `c(min, max)`. Only called if actually referenced.
#' @param ranged_range A zero-argument function returning the `c(min, max)`
#' of exactly what will be plotted. Only called if actually referenced.
#' @return A numeric vector `c(min, max)`
#' @noRd
resolve_intensity_range <- function(intensity_range, global_range, ranged_range) {
  resolve_endpoint <- function(value, position) {
    if (is.numeric(value) && length(value) == 1) {
      return(value)
    }
    if (identical(value, "global")) {
      return(global_range()[position])
    }
    if (identical(value, "ranged")) {
      return(ranged_range()[position])
    }
    cli_abort('Each element of intensity_range should be a number, "global" or "ranged"')
  }

  if (is.character(intensity_range) && length(intensity_range) == 1) {
    if (intensity_range == "global") {
      return(global_range())
    } else if (intensity_range == "ranged") {
      return(ranged_range())
    }
    cli_abort('intensity_range should be "global", "ranged", or a length-2 vector/list of numbers, "global" or "ranged"')
  }

  if (length(intensity_range) != 2) {
    cli_abort('intensity_range should be "global", "ranged", or a length-2 vector/list of numbers, "global" or "ranged"')
  }

  c(
    resolve_endpoint(intensity_range[[1]], 1L),
    resolve_endpoint(intensity_range[[2]], 2L)
  )
}

#' Make a plot interactive
#'
#' Wraps the `plt` with [plotly::ggplotly()] and sets the `xaxis` and `yaxis`
#' ticks to `"auto"`, so the axis labels are updated when zooming.
#'
#' @param plt A ggplot plot
#'
#' @return A plotly plot
#' @export
#'
#' @examples
#' d <- data.frame(x = c(1,2), y=c(1,2))
#' plt <- ggplot2::ggplot(d) + 
#'   ggplot2::geom_point(ggplot2::aes(x = x, y = y))
#' plot_interactive(plt)
plot_interactive <- function(plt) {
  require_pkgs("plotly")
  plotly::layout(
    plotly::ggplotly(plt),
    xaxis = list(tickmode = "auto"),
    yaxis = list(tickmode = "auto")
  )
}
