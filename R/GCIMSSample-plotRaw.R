#' Topographical plot of a GCIMSSample object
#'
#' @param object A [GCIMSSample] object
#' @param dt_range The drift time range to plot (in milliseconds)
#' @param rt_range The retention time range to plot (in seconds)
#' @return A plot of the GCIMSSample
#' @examples
#' dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3),
#'   gc_column = "Optional column name",
#'   drift_gas = "nitrogen",
#'   drift_tube_length = 98.0 # in mm
#' )
#' plotRaw(dummy_obj)
#' @export
setMethod(
  "plotRaw",
  "GCIMSSample",
  function(object, dt_range = NULL, rt_range = NULL) {
    intmat <- intensity(object, dt_range = dt_range, rt_range = rt_range)
    intens_long <- reshape2::melt(intmat, value.name = "Intensity")

    cubic_root_trans <- scales::trans_new(
      name = "cubic_root",
      transform = function(x) sign(x)*abs(x)^(1/3),
      inverse = function(x) sign(x)*abs(x)^3,
      breaks = function(x, n = 5) {
        x <- x[is.finite(x)]
        if (length(x) == 0) {
          return(numeric())
        }
        rng <- range(x)
        rng <- sign(rng)*abs(rng)^(1/3)
        out <- labeling::extended(rng[1], rng[2], n)
        out <- sign(out)*abs(out)^3
        out
      }
    )



    p <- ggplot2::ggplot(
      intens_long,
      mapping = ggplot2::aes(x = .data$dt_ms, y = .data$rt_s, fill = .data$Intensity)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      viridis::scale_fill_viridis(discrete = FALSE, option = "A", direction = -1, trans = cubic_root_trans) +
      ggplot2::labs(
        x = "Drift Time (ms)",
        y = "Retention Time (s)",
        fill = "Intensity"
      ) +
      ggplot2::theme_minimal()
    p
  }
)
