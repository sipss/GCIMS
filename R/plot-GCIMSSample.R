#' Topographical plot of a GCIMSSample object
#'
#' @param x A [GCIMSSample] object
#' @param dt_range A numeric vector of length 2 with the drift time range to plot (in milliseconds)
#' @param rt_range A numeric vector of length 2 with the retention time range to plot (in seconds)
#' @param ... Ignored
#' @param remove_baseline Set to `TRUE` to subtract the estimated baseline first
#' @param trans The transformation to the intensity values. "cubic_root" is the default. "intensity" is also valid.
#' See the `trans` argument in [ggplot2::continuous_scale()] for other possibilities.
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
#' plot(dummy_obj)
#' @export
setMethod(
  "plot",
  "GCIMSSample",
  function(x, dt_range = NULL, rt_range = NULL, ..., remove_baseline = FALSE, trans = "cubic_root") {
  dt <- dtime(x)
  rt <- rtime(x)
  idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range)

  intmat <- intensity(x, idx)
  if (isTRUE(remove_baseline)) {
    basel <- baseline(x)[idx$dt_logical, idx$rt_logical]
    intmat <- intmat - basel
  }
  mat_to_gplot(
    intmat,
    dt_min = idx$dt_ms_min,
    dt_max = idx$dt_ms_max,
    rt_min = idx$rt_s_min,
    rt_max = idx$rt_s_max,
    trans = trans
  )
})


mat_to_gplot <- function(intmat, dt_min = NULL, dt_max = NULL, rt_min = NULL, rt_max = NULL, trans = "cubic_root") {
  require_pkgs(c("farver", "viridisLite"))
  if (is.null(dt_min)) {
    dt_min <- as.numeric(rownames(intmat)[1L])
  }
  if (is.null(dt_max)) {
    dt_max <- as.numeric(rownames(intmat)[nrow(intmat)])
  }
  if (is.null(rt_min)) {
    rt_min <- as.numeric(colnames(intmat)[1L])
  }
  if (is.null(rt_max)) {
    rt_max <- as.numeric(colnames(intmat)[ncol(intmat)])
  }
  minmax <- range(intmat)

  if (is.character(trans)) {
    trans_func <- paste0(trans, "_trans")
    scales_envir <- getNamespace("scales")
    if (trans_func == "cubic_root_trans") {
      trans <- cubic_root_trans()
    } else if (trans_func %in% ls(scales_envir)) {
      trans <- get(trans_func, envir = scales_envir)()
    } else {
      rlang::abort("unknown trans value")
    }
  } else if (!inherits(trans, "trans")) {
    rlang::abort("unknown trans value")
  }
  intmat_trans <- trans$transform(intmat)
  colormap <- farver::encode_native(
    viridisLite::viridis(256L, direction = -1, option = "A")
  )
  nr <- mat_to_nativeRaster(intmat_trans, colormap)

  # The geom_rect is fake and it is only used to force the fill legend to appear
  # The geom_rect  limits are used to help set the plot limits
  # The geom_rect data, that contains the limits as well, is there because ggplotly
  # raises a javascript error otherwise: "Uncaught Error: Something went wrong with axis scaling"
  # in setScale (bug not yet reported to plotly)
  gplt <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = dt_min, xmax = dt_min,
      ymin = rt_min, ymax = rt_min,
      ggplot2::aes(fill = .data$x),
      data = data.frame(
        x = NA_real_,
        dt_ms_min = dt_min, dt_ms_max = dt_max,
        rt_s_min = rt_min, rt_s_max = rt_max
      )
    ) +
    ggplot2::annotation_raster(
      nr,
      xmin = dt_min, xmax = dt_max,
      ymin = rt_min, ymax = rt_max
    ) +
    ggplot2::scale_fill_viridis_c( # This has to match with the COLORMAP above
      direction = -1,
      option = "A",
      limits = minmax,
      na.value = "#00000000",
      trans = trans
    ) +
    ggplot2::lims(
      x = c(dt_min, dt_max),
      y = c(rt_min, rt_max)
    ) +
    ggplot2::labs(
      x = "Drift time (ms)",
      y = "Retention time (s)",
      fill = "Intensity (a.u.)"
    ) +
    ggplot2::theme_minimal()
  gplt
}

#' Cubic root transformation
#'
#' A scales transformation to be used with ggplot2.
#'
#' This function is exported because we are using it in vignettes, but it may
#' become unavailable in future versions
#'
#' @return A scale transformation object of name "cubic_root"
#'
#' @examples
#' library(ggplot2)
#' x <- 1:10
#' y <- x^3
#' df <- data.frame(x = x, y = y)
#' ggplot(data.frame(x=x, y=y)) +
#'  geom_point(aes(x = x, y = y)) +
#'   scale_y_continuous(trans=cubic_root_trans())
#' @export
cubic_root_trans <- function() {
  require_pkgs(c("scales", "labeling"))
  scales::trans_new(
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
}

#' Turn the intensity matrix into a data frame
#'
#' @param x A [GCIMSSample] object
#' @inheritParams intensity,GCIMSSample-method
#' @inheritParams base::as.data.frame
#' @param ... unused
#' @return A data frame with `dt_ms`, `rt_s` and `Intensity` columns
#' @export
as.data.frame.GCIMSSample <- function(x, row.names = NULL, optional = FALSE, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, ...) {
  intmat <- intensity(x, dt_range = dt_range, rt_range = rt_range, rt_idx = rt_idx, dt_idx = dt_idx)
  intens_long <- reshape2::melt(intmat, value.name = "Intensity")
  if (!is.null(row.names)) {
    rownames(intens_long) <- row.names
  }
  intens_long
}

# P40 <- unname(Polychrome::createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80)))
# Polychrome::swatch(P40)
P40 <- c("#FE1C2E", "#1CE50D", "#1C0DFE", "#E6B8C2", "#FD00E3", "#0DD1FC",
         "#FC9900", "#226516", "#A80D47", "#70228D", "#CACD00", "#386376",
         "#3BDDC0", "#C70DFD", "#7F5A00", "#EE76C4", "#859BFE", "#F900A7",
         "#E4B4FE", "#733556", "#C1CB97", "#F49386", "#D94B00", "#7CDA78",
         "#F92A75", "#F27CFF", "#EFBB5F", "#A3CAC8", "#005698", "#7A0DC8",
         "#685F51", "#A2BFEE", "#91380D", "#957AA1", "#7A9A00", "#950068",
         "#C76D89", "#00B1B9", "#C187FF", "#0042BD")

#' Add peak list rectangles to a raw plot
#'
#'
#' @param plt The output of [plot()] when applied to a [GCIMSSample]
#' @param peaklist A data frame with at least the columns: `dt_min_ms`, `dt_max_ms`, `rt_min_s`, `rt_max_s`
#' and optionally additional columns (e.g. the column given to `color_by`)
#' @param color_by A character with a column name of `peaklist`. Used to color the border of
#' the added rectangles
#' @param col_prefix After clustering, besides `dt_min_ms`, we also have
#' @param pdata A phenotype data data frame, with a SampleID column to be merged into peaklist so color_by can specify
#' a phenotype
#' `freesize_dt_min_ms`. Use `col_prefix = "freesize_"` to plot the `freesize` version
#' @param palette A character vector with color names to use drawing the rectangles. Use `NULL` to let `ggplot2` set the defaults.
#' @details
#' If `peaklist` includes `dt_apex_ms` and `rt_apex_s` a cross will be plotted on the peak apex.
#'
#' @return The given `plt` with rectangles showing the ROIs and crosses showing the apexes
#' @examples
#' dt <- 1:10
#' rt <- 1:10
#' int <- matrix(0.0, nrow = length(dt), ncol = length(rt))
#'
#' int[2, 4:8] <- c(.5, .5, 1, .5, 0.5)
#' int[3, 4:8] <- c(0.5, 2, 2, 2, 0.5)
#' int[4, 4:8] <- c(1, 2, 5, 2, 1)
#' int[5, 4:8] <- c(0.5, 2, 2, 2, 0.5)
#' int[6, 4:8] <- c(.5, .5, 1, .5, 0.5)
#'
#' dummy_obj <-GCIMSSample(
#'   drift_time = dt,
#'   retention_time = rt,
#'   data = int
#' )
#' plt <- plot(dummy_obj)
#'
#' # Add a rectangle on top of the plot
#' rect <- data.frame(
#'   dt_min_ms = 2.75,
#'   dt_max_ms = 5.6,
#'   rt_min_s = 4.6,
#'   rt_max_s = 7.4
#' )
#'
#' add_peaklist_rect(
#'   plt = plt,
#'   peaklist = rect
#' )
#' @export
add_peaklist_rect <- function(plt, peaklist, color_by = NULL, col_prefix = "", pdata = NULL,
                              palette = P40) {
  cli_warn(
    c(
      "{.code add_peaklist_rect} is deprecated",
      "i" = "Replace {.code add_peaklist_rect(plt, peaklist)} with {.code plt + overlay_peaklist(peaklist)} instead"
    )
  )
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
    if (!is.null(pdata)) {
      peaklist <- dplyr::left_join(peaklist, pdata, by = "SampleID")
    }
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
    show_legend <- nrow(peaklist_to_plot_rect) <= 10
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
        alpha = 0.3,
        show.legend = show_legend
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
          shape = "x",
          show.legend = show_legend
        )
    }
    if (!is.null(palette)) {
      palette <- unname(palette)
      colours_in_pal <- length(palette)
      if (colours_in_pal >= how_many_colors) {
        palette <- palette[seq_len(how_many_colors)]
      } else {
        # recycle palette:
        recycle_times <- floor(how_many_colors / colours_in_pal)
        palette <- rep(palette, times = recycle_times)
        palette <- c(palette, palette[seq_len(how_many_colors - length(palette))])
      }
      plt <- plt +
        ggplot2::scale_color_manual(values = palette)
    }
  }
  plt
}



#' Overlay a peak list to a plot
#'
#' @param peaklist A data frame with at least the columns: `dt_min_ms`, `dt_max_ms`, `rt_min_s`, `rt_max_s`
#' and optionally additional columns (e.g. the column given to `color_by`)
#' @param ... Ignored.
#' @inheritParams dt_rt_range_normalization
#' @param pdata A phenotype data data frame, with a `SampleID` column. This column will be used to
#' merge `pdata` with `peaklist`, so `color_by` can specify a phenotype.
#' @param color_by A character with a column name of `peaklist` or `pdata`. Used to color the border of
#' the added rectangles and apices. A string with a color name is also acceptable.
#' @param mapping_roi A 4-elements named character vector with the names of the columns from `peaklist` that will
#' be used as the rectangle coordinates.
#' @param palette A character vector with color names to use drawing the rectangles. Use `NULL` to let `ggplot2` set the defaults.
#' @details
#' If `peaklist` includes `dt_apex_ms` and `rt_apex_s` a cross will be plotted on the peak apex.
#'
#' @return A list with the ggplot layers to overlay, including `geom_rect` and possibly `geom_point` and `scale_fill_manual`.
#' @examples
#' dt <- 1:10
#' rt <- 1:10
#' int <- matrix(0.0, nrow = length(dt), ncol = length(rt))
#'
#' int[2, 4:8] <- c(.5, .5, 1, .5, 0.5)
#' int[3, 4:8] <- c(0.5, 2, 2, 2, 0.5)
#' int[4, 4:8] <- c(1, 2, 5, 2, 1)
#' int[5, 4:8] <- c(0.5, 2, 2, 2, 0.5)
#' int[6, 4:8] <- c(.5, .5, 1, .5, 0.5)
#'
#' dummy_obj <-GCIMSSample(
#'   drift_time = dt,
#'   retention_time = rt,
#'   data = int
#' )
#' plt <- plot(dummy_obj)
#'
#' # Add a rectangle on top of the plot
#' rect <- data.frame(
#'   dt_min_ms = 2.75,
#'   dt_max_ms = 5.6,
#'   rt_min_s = 4.6,
#'   rt_max_s = 7.4
#' )
#'
#' plt + overlay_peaklist(
#'   peaklist = rect
#' )
#' @export
#'
overlay_peaklist <- function(
    peaklist = NULL,
    ...,
    dt_range = NULL,
    rt_range = NULL,
    pdata = NULL,
    color_by = NULL,
    mapping_roi = c(
      "dt_min_ms" = "dt_min_ms",
      "dt_max_ms" = "dt_max_ms",
      "rt_min_s" = "rt_min_s",
      "rt_max_s" = "rt_max_s"
    ),
    palette = P40
) {
  if (is.null(dt_range)) {
    dt_range <- c(-Inf, Inf)
  }
  if (is.null(rt_range)) {
    rt_range <- c(-Inf, Inf)
  }

  if (!inherits(peaklist, "data.frame")) {
    # coerces if it's a S4Vectors::DataFrame.
    # does not coerce if it's a tibble
    peaklist <- as.data.frame(peaklist)
  }
  if (!is.null(pdata)) {
    fulldata <- dplyr::left_join(peaklist, pdata, by = "SampleID")
  } else {
    fulldata <- peaklist
  }

  if (is.null(color_by)) {
    color_by <- "green"
    how_many_colors <- 1L
  } else if (purrr::is_scalar_character(color_by)) {
    if (color_by %in% colnames(fulldata)) {
      color_by_sym <- rlang::sym(color_by)
      how_many_colors <- length(unique(fulldata[[color_by]]))
    } else {
      is_color <- tryCatch({
        farver::encode_native(color_by)
        TRUE
      }, error = function(e) FALSE)
      if (is_color) {
        how_many_colors <- 1L
      } else {
        cli_abort(
          c(
            "Invalid {.arg color_by}",
            "x" = "{.arg color_by} {color_by} is neither a valid color name nor a column of peaklist or pdata"
          )
        )
      }
    }
  } else {
    cli_abort("{.arg color_by} should be {.code NULL} or a string")
  }

  if (!identical(names(mapping_roi), c("dt_min_ms", "dt_max_ms", "rt_min_s", "rt_max_s"))) {
    cli_abort(
      c("i" = 'mapping_roi should be a named vector with names c("dt_min_ms", "dt_max_ms", "rt_min_s", "rt_max_s")')
    )
  }

  fulldata <- fulldata |>
    dplyr::select(-dplyr::any_of(names(mapping_roi)[!mapping_roi %in% names(mapping_roi)])) |>
    dplyr::rename(dplyr::all_of(mapping_roi))


  peaklist_to_plot_rect <- fulldata |>
    dplyr::filter(
      .data$dt_min_ms <= max(dt_range),
      .data$dt_max_ms >= min(dt_range),
      .data$rt_min_s <= max(rt_range),
      .data$rt_max_s >= min(rt_range)
    )

  has_apex <- all(c("dt_apex_ms", "rt_apex_s") %in% colnames(fulldata))

  if (has_apex) {
    peaklist_to_plot_apex <- fulldata |>
      dplyr::filter(
        .data$dt_apex_ms >= min(dt_range),
        .data$dt_apex_ms <= max(dt_range),
        .data$rt_apex_s >= min(rt_range),
        .data$rt_apex_s <= max(rt_range)
      )
  }

  out <- list()
  if (how_many_colors == 1) {
    # Only one color, use green
    out <- c(
      out,
      ggplot2::geom_rect(
        data = peaklist_to_plot_rect,
        mapping = ggplot2::aes(
          xmin = .data$dt_min_ms,
          xmax = .data$dt_max_ms,
          ymin = .data$rt_min_s,
          ymax = .data$rt_max_s
        ),
        color = color_by,
        alpha = 0.3
      )
    )
    if (has_apex) {
      out <- c(
        out,
        ggplot2::geom_point(
          data = peaklist_to_plot_apex,
          mapping = ggplot2::aes(
            x = .data$dt_apex_ms,
            y = .data$rt_apex_s
          ),
          color = color_by,
          shape = "x"
        )
      )
    } else {
      # This else block is needed just because of plotly:
      # plotly can't place an annotation_raster between geoms, it must be either below or
      # above the whole plot.
      # See https://github.com/sipss/GCIMS/issues/16
      #
      # When plotly finds a raster between geoms, chooses to place it "above" or "below"
      # all of them by checking how many layers it has above and how many below, and picking the
      # one with the majority (e.g. if the raster is the 3rd out of 4 layers it would be placed above)
      # so to force it to appear below, we add a geom_blank on top, that will be invisible, but will
      # make the annotation_raster to be 2/4 (First half) instead of 2/3 (second out of 3)
      out <- c(
        out,
        ggplot2::geom_blank()
      )
    }
  } else {
    show_legend <- how_many_colors <= 10
    out <- c(
      out,
      ggplot2::geom_rect(
        data = peaklist_to_plot_rect,
        mapping = ggplot2::aes(
          xmin = .data$dt_min_ms,
          xmax = .data$dt_max_ms,
          ymin = .data$rt_min_s,
          ymax = .data$rt_max_s,
          color = !!color_by_sym
        ),
        alpha = 0.3,
        show.legend = show_legend
      ),
      ggplot2::geom_blank()
    )
    if (has_apex) {
      out <- c(
        out,
        ggplot2::geom_point(
          data = peaklist_to_plot_apex,
          mapping = ggplot2::aes(
            x = .data$dt_apex_ms,
            y = .data$rt_apex_s,
            color = !!color_by_sym
          ),
          shape = "x",
          show.legend = show_legend
        )
      )
    } else {
      # This else block is needed just because of plotly:
      # plotly can't place an annotation_raster between geoms, it must be either below or
      # above the whole plot.
      # See https://github.com/sipss/GCIMS/issues/16
      #
      # When plotly finds a raster between geoms, chooses to place it "above" or "below"
      # all of them by checking how many layers it has above and how many below, and picking the
      # one with the majority (e.g. if the raster is the 3rd out of 4 layers it would be placed above)
      # so to force it to appear below, we add a geom_blank on top, that will be invisible, but will
      # make the annotation_raster to be 2/4 (First half) instead of 2/3 (second out of 3)
      out <- c(
        out,
        ggplot2::geom_blank()
      )
    }
    if (!is.null(palette)) {
      palette <- unname(palette)
      colours_in_pal <- length(palette)
      if (colours_in_pal >= how_many_colors) {
        palette <- palette[seq_len(how_many_colors)]
      } else {
        # recycle palette:
        recycle_times <- floor(how_many_colors / colours_in_pal)
        palette <- rep(palette, times = recycle_times)
        palette <- c(palette, palette[seq_len(how_many_colors - length(palette))])
      }
      out <- c(
        out,
        ggplot2::scale_color_manual(values = palette)
      )
    }
    # Even if we used show.legend = FALSE, plotly (4.10.1.9000) needs
    # guides(colour = "none").
    # So we force it here anyway:
    if (!show_legend) {
      out <- c(
        out,
        list(ggplot2::guides(colour = "none"))
      )
    }
  }
  out
}
