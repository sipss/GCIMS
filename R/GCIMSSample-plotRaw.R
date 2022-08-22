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
    intens_long <- tidy(object, dt_range = dt_range, rt_range = rt_range)

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



    p <- ggplot2::ggplot(intens_long) +
      ggplot2::geom_raster(
        mapping = ggplot2::aes(x = .data$dt_ms, y = .data$rt_s, fill = .data$Intensity),
        interpolate = FALSE
      ) +
      viridis::scale_fill_viridis(discrete = FALSE, option = "A", direction = -1, trans = cubic_root_trans) +
      ggplot2::labs(
        x = "Drift Time (ms)",
        y = "Retention Time (s)",
        fill = "Intensity",
        caption = object@description
      ) +
      ggplot2::theme_minimal()
    p
  }
)

#' Tidy the raw data into a long data frame
#'
#' @param x A [GCIMSSample] object
#' @inheritParams intensity,GCIMSSample-method
#' @param ... unused
#' @return A data frame with `dt_ms`, `rt_s` and `Intensity` columns
#' @export
tidy.GCIMSSample <- function(x, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, ...) {
  intmat <- intensity(x, dt_range = dt_range, rt_range = rt_range, rt_idx = rt_idx, dt_idx = dt_idx)
  intens_long <- reshape2::melt(intmat, value.name = "Intensity")
  intens_long
}


#' Add peak list rectangles to a raw plot
#'
#' @param plt The output of [plotRaw()] when applied to a [GCIMSSample]
#' @param peaklist The output of [peaks()] of a GCIMSSample
#'
#' @return The given `plt` with rectangles showing the ROIs and crosses showing the apexes
#' @export
#'
add_peaklist_rect <- function(plt, peaklist) {
  dt_range <- ggplot2::layer_scales(plt)$x$range$range
  rt_range <- ggplot2::layer_scales(plt)$y$range$range
  if (is.null(dt_range)) {
    dt_range <- c(-Inf, Inf)
  }
  if (is.null(rt_range)) {
    rt_range <- c(-Inf, Inf)
  }

  peaklist <- as.data.frame(peaklist)

  peaklist_to_plot_rect <- dplyr::filter(
    peaklist,
    .data$dt_min_ms <= max(dt_range),
    .data$dt_max_ms >= min(dt_range),
    .data$rt_min_s <= max(rt_range),
    .data$rt_max_s >= min(rt_range)
  )

  peaklist_to_plot_apex <- dplyr::filter(
    peaklist,
    .data$dt_apex_ms >= min(dt_range),
    .data$dt_apex_ms <= max(dt_range),
    .data$rt_apex_s >= min(rt_range),
    .data$rt_apex_s <= max(rt_range)
  )

  if (length(unique(peaklist$SampleID)) == 1) {
    # Only one sample, use green
    plt <- plt +
      ggplot2::geom_rect(
        data = peaklist_to_plot_rect,
        mapping = ggplot2::aes(
          xmin = .data$dt_min_ms,
          xmax = .data$dt_max_ms,
          ymin = .data$rt_min_s,
          ymax = .data$rt_max_s
        ),
        color = "green",
        alpha = 0.3
      ) +
      ggplot2::geom_point(
        data = peaklist_to_plot_apex,
        mapping = ggplot2::aes(
          x = .data$dt_apex_ms,
          y = .data$rt_apex_s
        ),
        color = "green",
        shape = "x"
      )

  } else {
    plt <- plt +
      ggplot2::geom_rect(
        data = peaklist_to_plot_rect,
        mapping = ggplot2::aes(
          xmin = .data$dt_min_ms,
          xmax = .data$dt_max_ms,
          ymin = .data$rt_min_s,
          ymax = .data$rt_max_s,
          color = .data$SampleID
        ),
        alpha = 0.3
      ) +
      ggplot2::geom_point(
        data = peaklist_to_plot_apex,
        mapping = ggplot2::aes(
          x = .data$dt_apex_ms,
          y = .data$rt_apex_s,
          color = .data$SampleID
        ),
        shape = "x"
      )

  }
  plt
}

