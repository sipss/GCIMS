#' @describeIn GCIMSChromatogramSet-class plot method
#' @param x A [GCIMSChromatogramSet] object to plot
#' @param color_by The name of a `pData(x)` column (or `"SampleID"`) used to
#' color the chromatograms
#' @param ... Ignored
#' @return A ggplot2 plot object
#' @export
setMethod(
  "plot",
  "GCIMSChromatogramSet",
  function(x, color_by = "SampleID", ...) {
    sample_names <- sampleNames(x)
    if (length(sample_names) == 0) {
      cli_abort("Can't plot an empty GCIMSChromatogramSet")
    }
    df <- dplyr::bind_rows(purrr::map(
      sample_names,
      function(sample_id) {
        chrom <- x[[sample_id]]
        data.frame(
          SampleID = sample_id,
          retention_time_s = rtime(chrom),
          intensity = unname(intensity(chrom))
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
          x = .data$retention_time_s,
          y = .data$intensity,
          color = .data[[color_by]],
          group = .data$SampleID
        )
      ) +
      ggplot2::labs(
        x = "Retention time (s)",
        y = "Intensity (a.u.)",
        color = color_by
      )
  }
)
