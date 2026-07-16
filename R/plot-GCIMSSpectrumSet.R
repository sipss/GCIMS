#' @describeIn GCIMSSpectrumSet-class plot method
#' @param x A [GCIMSSpectrumSet] object to plot
#' @param color_by The name of a `pData(x)` column (or `"SampleID"`) used to
#' color the spectra
#' @param ... Ignored
#' @return A ggplot2 plot object
#' @export
setMethod(
  "plot",
  "GCIMSSpectrumSet",
  function(x, color_by = "SampleID", ...) {
    sample_names <- sampleNames(x)
    if (length(sample_names) == 0) {
      cli_abort("Can't plot an empty GCIMSSpectrumSet")
    }
    df <- dplyr::bind_rows(purrr::map(
      sample_names,
      function(sample_id) {
        spec <- x[[sample_id]]
        data.frame(
          SampleID = sample_id,
          drift_time_ms = dtime(spec),
          intensity = unname(intensity(spec))
        )
      }
    ))

    pd <- pData(x)
    if (!is.null(pd) && "SampleID" %in% colnames(pd)) {
      df <- dplyr::left_join(df, pd, by = "SampleID")
    }
    if (!color_by %in% colnames(df)) {
      cli_abort("{.val {color_by}} is not a column of {.code pData(x)} (or {.val SampleID})")
    }

    ggplot2::ggplot(df) +
      ggplot2::geom_line(
        mapping = ggplot2::aes(
          x = .data$drift_time_ms,
          y = .data$intensity,
          color = .data[[color_by]],
          group = .data$SampleID
        )
      ) +
      ggplot2::labs(
        x = "Drift time (ms)",
        y = "Intensity (a.u.)",
        color = color_by
      )
  }
)
