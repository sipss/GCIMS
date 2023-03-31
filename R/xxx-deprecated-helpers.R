#' Get a sample from a GCIMSDataset object
#'
#' @param object a [GCIMSDataset] object
#' @param sample Either an integer (sample index) or a string (sample name)
#' @return The corresponding [GCIMSSample]
#' @keywords internal
#' @export
getSample <- function(object, sample) {
  cli_warn(
    "getSample(object, sample) is deprecated. Use object$getSample(sample) instead.",
    frequency = "once",
    frequency_id = "getSample-method-deprecated"
  )
  object$getSample(sample)
}

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


#' Runs all delayed operations on the object
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @param keep_intermediate A logical, whether to keep the intermediate files of
#' the previous realization once this one finishes. If `NA`, keeping will depend
#' on the `object`.
#' @return The same [GCIMSDataset] object, without pending operations
#' @export
#' @examples
#' base_dir <- system.file("extdata", "sample_formats", package = "GCIMS")
#' annot <- data.frame(SampleID = "Sample1", FileName = "small.mea.gz")
#' dataset <- GCIMSDataset$new(annot, base_dir)
#' print(dataset)
#' realize(dataset)
#' print(dataset)
#'
realize <- function(object, keep_intermediate = NA) {
  cli_warn(
    "realize(object) is deprecated, use object$realize() instead",
    .frequency = "once",
    .frequency_id = "realize-deprecated",
  )
  object$realize(keep_intermediate = keep_intermediate)
}


setMethod(
  "show",
  "GCIMSDataset",
  function(object) {
    cli_warn(
      "show(gcimsdataset) is deprecated. Use gcimsdataset$print() instead.",
      frequency = "once",
      frequency_id = "show-gcimsdataset-deprecated"
    )
    object$print()
  }
)
