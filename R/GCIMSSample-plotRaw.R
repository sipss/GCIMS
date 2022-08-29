#' Topographical plot of a GCIMSSample object
#'
#' @param object A [GCIMSSample] object
#' @param dt_range A numeric vector of length 2 with the drift time range to plot (in milliseconds)
#' @param rt_range A numeric vector of length 2 with the retention time range to plot (in seconds)
#' @param ... Ignored
#' @param remove_baseline Set to `TRUE` to subtract the estimated baseline first
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
  function(object, dt_range = NULL, rt_range = NULL, ..., remove_baseline = FALSE) {
    dt <- dtime(object)
    rt <- rtime(object)
    idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range)

    intmat <- intensity(object, idx)
    if (isTRUE(remove_baseline)) {
      basel <- baseline(object)[idx$dt_logical, idx$rt_logical]
      intmat <- intmat - basel
    }
    minmax <- range(intmat)

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

    intmat_trans <- cubic_root_trans$transform(intmat)
    nr <- mat_to_nativeRaster(intmat_trans, COLORMAP_VIRIDIS_256_A_m1)

    # The geom_rect is fake and it is only used to force the fill legend to appear
    gplt <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        xmin = idx$dt_ms_min, xmax = idx$dt_ms_min,
        ymin = idx$rt_s_min, ymax = idx$rt_s_min,
        ggplot2::aes(fill = .data$x),
        data = data.frame(x = NA_real_)
      ) +
      ggplot2::annotation_raster(
        nr,
        xmin = idx$dt_ms_min, xmax = idx$dt_ms_max,
        ymin = idx$rt_s_min, ymax = idx$rt_s_max
      ) +
      ggplot2::scale_fill_viridis_c( # This has to match with the COLORMAP above
        direction = -1,
        option = "A",
        limits = minmax,
        na.value = "#00000000",
        trans = cubic_root_trans
      ) +
      ggplot2::lims(
        x = c(idx$dt_ms_min, idx$dt_ms_max),
        y = c(idx$rt_s_min, idx$rt_s_max)
      ) +
      ggplot2::labs(
        x = "Drift time (ms)",
        y = "Retention time (s)",
        fill = "Intensity (a.u.)"
      ) +
      ggplot2::theme_minimal()
    gplt
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
#'
#' @param plt The output of [plotRaw()] when applied to a [GCIMSSample]
#' @param peaklist A data frame with at least the columns: `dt_min_ms`, `dt_max_ms`, `rt_min_s`, `rt_max_s`
#' and optionally additional columns (e.g. the column given to `color_by`)
#' @param color_by A character with a column name of `peaklist`. Used to color the border of
#' the added rectangles
#'
#' @details
#' If `peaklist` includes `dt_apex_ms` and `rt_apex_s` a cross will be plotted on the peak apex.
#'
#' @return The given `plt` with rectangles showing the ROIs and crosses showing the apexes
#' @export
#'
add_peaklist_rect <- function(plt, peaklist, color_by = NULL, col_prefix = "") {
  dt_range <- ggplot2::layer_scales(plt)$x$range$range
  rt_range <- ggplot2::layer_scales(plt)$y$range$range
  if (is.null(dt_range)) {
    dt_range <- c(-Inf, Inf)
  }
  if (is.null(rt_range)) {
    rt_range <- c(-Inf, Inf)
  }

  if (is.null(color_by)) {
    color_by_sym <- NULL
    how_many_colors <- 1L
  } else {
    color_by_sym <- rlang::sym(color_by)
    how_many_colors <- length(unique(peaklist[[color_by]]))
  }

  peaklist <- as.data.frame(peaklist)
  if (col_prefix != "") {
    datacols <- c("dt_min_ms", "dt_max_ms", "rt_min_s", "rt_max_s")
    datasyms <- rlang::data_syms(paste0(col_prefix, datacols))
    names(datasyms) <- datacols
    peaklist <- peaklist |>
      dplyr::select(-dplyr::any_of(datacols)) |>
      dplyr::rename(!!!datasyms)
  }

  peaklist_to_plot_rect <- peaklist |>
    dplyr::filter(
      .data$dt_min_ms <= max(dt_range),
      .data$dt_max_ms >= min(dt_range),
      .data$rt_min_s <= max(rt_range),
      .data$rt_max_s >= min(rt_range)
    )

  has_apex <- all(c("dt_apex_ms", "rt_apex_s") %in% colnames(peaklist))

  if (has_apex) {
    peaklist_to_plot_apex <- peaklist |>
      dplyr::filter(
        .data$dt_apex_ms >= min(dt_range),
        .data$dt_apex_ms <= max(dt_range),
        .data$rt_apex_s >= min(rt_range),
        .data$rt_apex_s <= max(rt_range)
      )
  }

  if (how_many_colors == 1) {
    # Only one color, use green
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
      )
    if (has_apex) {
      plt <- plt +
        ggplot2::geom_point(
          data = peaklist_to_plot_apex,
          mapping = ggplot2::aes(
            x = .data$dt_apex_ms,
            y = .data$rt_apex_s
          ),
          color = "green",
          shape = "x"
        )
    }
  } else {
    plt <- plt +
      ggplot2::geom_rect(
        data = peaklist_to_plot_rect,
        mapping = ggplot2::aes(
          xmin = .data$dt_min_ms,
          xmax = .data$dt_max_ms,
          ymin = .data$rt_min_s,
          ymax = .data$rt_max_s,
          color = !!color_by_sym
        ),
        alpha = 0.3
      )
    if (has_apex) {
      plt <- plt +
        ggplot2::geom_point(
          data = peaklist_to_plot_apex,
          mapping = ggplot2::aes(
            x = .data$dt_apex_ms,
            y = .data$rt_apex_s,
            color = !!color_by_sym
          ),
          shape = "x"
        )
    }
  }
  plt
}

