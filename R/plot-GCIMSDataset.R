#' Topographical plot of several samples of a GCIMSDataset, on a shared color scale
#'
#' Several samples are plotted side by side, on the same intensity color
#' scale, so they are visually comparable.
#'
#' @param x A [GCIMSDataset] object
#' @param sample A number, a string, or a vector of numbers/strings with the
#' sample index(es) or name(s) to plot. If `NULL` (the default), all samples
#' are plotted.
#' @inheritParams dt_rt_range_normalization
#' @param ... Ignored
#' @param remove_baseline Set to `TRUE` to subtract each sample's estimated
#' baseline first
#' @param trans The transformation to the intensity values, see `plot,GCIMSSample-method`
#' @param intensity_range Controls the shared color scale across all plotted
#' samples. One of:
#' - `"global"` (the default): each sample's full, uncropped intensity range,
#'   cached on the dataset (computed once, alongside TIS/RIC, regardless of
#'   `sample`/`dt_range`/`rt_range`). Free, and keeps the scale stable across
#'   different calls with different `sample`/`dt_range`/`rt_range` selections.
#'   Not valid together with `remove_baseline = TRUE`, since the cache holds
#'   raw intensities.
#' - `"ranged"`: the range of exactly what's plotted (respecting `dt_range`,
#'   `rt_range` and `remove_baseline`). Requires loading every selected
#'   sample's data twice (once to compute the range, once to plot).
#' - A numeric vector of length 2, `c(min, max)`: fixed limits.
#' - A length-2 vector/list whose elements are independently a number,
#'   `"global"` or `"ranged"`, e.g. `list(min = 0, max = "global")`.
#' @param ncol Number of columns for [cowplot::plot_grid()]. `NULL` lets `cowplot` decide.
#' @return A combined plot with one panel per sample
#' @export
setMethod(
  "plot",
  "GCIMSDataset",
  function(x, sample = NULL, dt_range = NULL, rt_range = NULL, ...,
           remove_baseline = FALSE, trans = "cubic_root",
           intensity_range = "global", ncol = NULL) {
    require_pkgs("cowplot")
    sample_names_all <- sampleNames(x)
    if (is.null(sample)) {
      sample <- sample_names_all
    }
    selected_names <- sample_name_or_number_to_both(sample, sample_names_all)$name

    samples <- stats::setNames(
      purrr::map(selected_names, function(nm) x$getSample(nm)),
      selected_names
    )

    cropped_intensity <- function(s) {
      dt <- dtime(s)
      rt <- rtime(s)
      idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range)
      intmat <- intensity(s, idx)
      if (isTRUE(remove_baseline)) {
        basel <- baseline(s)[idx[["dt_logical"]], idx[["rt_logical"]]]
        intmat <- intmat - basel
      }
      intmat
    }

    ranged_range_cache <- NULL
    get_ranged_range <- function() {
      if (is.null(ranged_range_cache)) {
        ranges <- purrr::map(samples, function(s) range(cropped_intensity(s)))
        ranged_range_cache <<- range(unlist(ranges))
      }
      ranged_range_cache
    }
    get_global_range <- function() {
      if (isTRUE(remove_baseline)) {
        cli_abort(
          c(
            "intensity_range can't use {.val global} together with {.code remove_baseline = TRUE}",
            "i" = "The cached global range holds raw intensities. Use {.code intensity_range = \"ranged\"} instead"
          )
        )
      }
      dataset_intensity_range(x)
    }

    limits <- resolve_intensity_range(intensity_range, get_global_range, get_ranged_range)

    panels <- purrr::map(selected_names, function(nm) {
      plot(
        samples[[nm]],
        dt_range = dt_range, rt_range = rt_range,
        remove_baseline = remove_baseline, trans = trans,
        intensity_range = limits
      ) + ggplot2::labs(title = nm)
    })

    cowplot::plot_grid(plotlist = panels, ncol = ncol)
  }
)
