#' Topographical plot of several samples of a GCIMSDataset, on a shared color scale
#'
#' Several samples are plotted side by side, on the same intensity color
#' scale, so they are visually comparable. Paginated: only one page's worth
#' of samples is ever loaded into memory at a time for rendering.
#'
#' @param x A [GCIMSDataset] object
#' @param sample A number, a string, or a vector of numbers/strings with the
#' sample index(es) or name(s) to plot. If `NULL` (the default), all samples
#' are plotted (across as many pages as needed).
#' @inheritParams dt_rt_range_normalization
#' @param ... Ignored
#' @param remove_baseline Set to `TRUE` to subtract each sample's estimated
#' baseline first
#' @param trans The transformation to the intensity values, see `plot,GCIMSSample-method`
#' @param intensity_range Controls the shared color scale across all plotted
#' samples (all pages, not just the current one). One of:
#' - `"global"` (the default): each sample's full, uncropped intensity range,
#'   regardless of `sample`/`dt_range`/`rt_range`. With `remove_baseline =
#'   FALSE`, this comes from a cache computed once, alongside TIS/RIC, so
#'   it's free. With `remove_baseline = TRUE`, the cache (raw intensities)
#'   doesn't apply, so each selected sample's full, uncropped, baseline-removed
#'   range is computed instead (loading every selected sample's data an extra
#'   time).
#' - `"ranged"`: the range of exactly what's plotted (respecting `dt_range`,
#'   `rt_range` and `remove_baseline`), across *every* selected sample, not
#'   just the current page -- so the scale is still comparable page to page.
#'   This is deliberately not page-scoped: a cheaper, page-local range would
#'   make pages incomparable to each other, defeating the purpose of paging
#'   through a dataset to compare samples. Costs loading every selected
#'   sample's data twice (once to compute the range, once to plot).
#' - A numeric vector of length 2, `c(min, max)`: fixed limits.
#' - A length-2 vector/list whose elements are independently a number,
#'   `"global"` or `"ranged"`, e.g. `list(min = 0, max = "global")`.
#' @param nrow,ncol Page grid shape, forwarded to [cowplot::plot_grid()]. Page
#' capacity is `nrow * ncol`. If both are `NULL` (the default), they're picked
#' from the number of selected samples: an exact fit for 1-6 samples (e.g. 5
#' or 6 samples -> 2x3), otherwise 3x3 (so page capacity maxes out at 9, and
#' more than 9 selected samples span multiple pages). If only one of
#' `nrow`/`ncol` is given, the other becomes `3` when there are more than 9
#' selected samples, or just enough to fit them all on one page otherwise.
#' @param page Which page to plot, 1-based. Errors if out of bounds.
#' @return A combined plot with one panel per sample, for the requested page
#' @export
setMethod(
  "plot",
  "GCIMSDataset",
  function(x, sample = NULL, dt_range = NULL, rt_range = NULL, ...,
           remove_baseline = FALSE, trans = "cubic_root",
           intensity_range = "global", nrow = NULL, ncol = NULL, page = 1) {
    require_pkgs("cowplot")
    sample_names_all <- sampleNames(x)
    if (is.null(sample)) {
      sample <- sample_names_all
    }
    selected_names <- sample_name_or_number_to_both(sample, sample_names_all)$name
    num_samples <- length(selected_names)
    if (num_samples == 0) {
      cli_abort("No samples to plot")
    }

    grid_dims <- resolve_page_grid(nrow, ncol, num_samples)
    nrow <- grid_dims[["nrow"]]
    ncol <- grid_dims[["ncol"]]
    page_size <- nrow * ncol
    n_pages <- ceiling(num_samples / page_size)
    if (page < 1 || page > n_pages) {
      cli_abort(
        c(
          "{.arg page} = {page} is out of bounds",
          "i" = "With {num_samples} sample(s) laid out on a {nrow}x{ncol} grid, there are {n_pages} page(s)"
        )
      )
    }
    page_start <- (page - 1) * page_size + 1
    page_end <- min(page * page_size, num_samples)
    page_names <- selected_names[page_start:page_end]

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
        # Scoped to every selected sample (all pages), not just the current
        # one, so the scale stays comparable across pages. Each sample is
        # loaded, ranged and discarded one at a time to avoid holding them
        # all in RAM at once.
        ranges <- purrr::map(selected_names, function(nm) {
          range(cropped_intensity(x$getSample(nm)))
        })
        ranged_range_cache <<- range(unlist(ranges))
      }
      ranged_range_cache
    }
    global_range_cache <- NULL
    get_global_range <- function() {
      if (is.null(global_range_cache)) {
        global_range_cache <<- if (isTRUE(remove_baseline)) {
          ranges <- purrr::map(selected_names, function(nm) {
            s <- x$getSample(nm)
            range(intensity(s) - baseline(s))
          })
          range(unlist(ranges))
        } else {
          dataset_intensity_range(x)
        }
      }
      global_range_cache
    }

    limits <- resolve_intensity_range(intensity_range, get_global_range, get_ranged_range)

    # Only the current page's samples are loaded for rendering, one at a
    # time, so peak memory never holds more than one raw sample plus the
    # (much smaller, native-raster-encoded) panels built so far.
    panels <- purrr::map(page_names, function(nm) {
      s <- x$getSample(nm)
      plot(
        s,
        dt_range = dt_range, rt_range = rt_range,
        remove_baseline = remove_baseline, trans = trans,
        intensity_range = limits
      ) + ggplot2::labs(title = nm)
    })

    cowplot::plot_grid(plotlist = panels, nrow = nrow, ncol = ncol)
  }
)
