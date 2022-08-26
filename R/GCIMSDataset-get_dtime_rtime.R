#' Get A reference drift time vector for the dataset
#' @param object A GCIMSDataset
#' @return a drift time vector
#' @export
setMethod("dtime", "GCIMSDataset", function(object, sample = 0L) {
  if (hasDelayedOps(object) || is.null(object@envir$dt_ref)) {
    object <- extract_dtime_rtime(object)
    object <- realize(object)
  }
  if (is.null(sample)) {
    return(object@envir$dt_ref)
  } else {
    sample_idx_name <- sample_name_or_number_to_both(sample)
    time_metrics <- as.list(object@envir$dt_rt_metrics[sample_idx_name$idx,])
    return(
      seq(
        from = time_metrics$dt_min_ms,
        to = time_metrics$dt_max_ms,
        length.out = time_metrics$dt_length_points
      )
    )
  }
})

#' Get a reference retention time vector for the dataset
#'
#' @param object A GCIMSDataset
#' @return a retention time vector
#' @export
setMethod("rtime", "GCIMSDataset", function(object, sample = NULL) {
  if (hasDelayedOps(object) || is.null(object@envir$rt_ref)) {
    object <- extract_dtime_rtime(object)
    object <- realize(object)
  }
  if (is.null(sample)) {
    return(object@envir$rt_ref)
  } else {
    sample_idx_name <- sample_name_or_number_to_both(sample)
    time_metrics <- as.list(object@envir$dt_rt_metrics[sample_idx_name$idx,])
    return(
      seq(
        from = time_metrics$rt_min_s,
        to = time_metrics$rt_max_s,
        length.out = time_metrics$rt_length_points
      )
    )
  }
})

.extract_dtime_rtime_fun_extract <- function(gcimssample) {
  # Extract the drift time and ret time range and steps
  dt <- dtime(gcimssample)
  rt <- rtime(gcimssample)
  # A data frame with one row will be extracted from each sample:
  data.frame(
    dt_min_ms = dt[1],
    dt_max_ms = dt[length(dt)],
    dt_step_ms = stats::median(diff(dt)),
    dt_length_points = length(dt),
    rt_min_s = rt[1],
    rt_max_s = rt[length(rt)],
    rt_step_s = stats::median(diff(rt)),
    rt_length_points = length(rt)
  )
}

.extract_dtime_rtime_fun_aggregate <- function(gcimsdataset, extracted_objects) {
  # extracted_objects is a named list. The names are the SampleIDs
  # the values are the output of fun_extract for the corresponding sample.
  # In our case, the value is a data frame with the drift and retention time ranges
  # so we bind_rows for all of them:
  dt_rt_metrics <- dplyr::bind_rows(extracted_objects, .id = "SampleID")
  max_dt_min <- max(dt_rt_metrics$dt_min_ms)
  min_dt_max <- min(dt_rt_metrics$dt_max_ms)
  min_dt_step <- min(dt_rt_metrics$dt_step_ms)
  max_dt_length <- floor(round((min_dt_max - max_dt_min)/min_dt_step, digits = 8)) + 1L
  max_rt_min <- max(dt_rt_metrics$rt_min_s)
  min_rt_max <- min(dt_rt_metrics$rt_max_s)
  min_rt_step <- min(dt_rt_metrics$rt_step_s)
  max_rt_length <- floor(round((min_rt_max - max_rt_min)/min_rt_step, digits = 8)) + 1L
  dt_ref <- seq(from = max_dt_min, to = min_dt_max, length.out = max_dt_length)
  rt_ref <- seq(from = max_rt_min, to = min_rt_max, length.out = max_rt_length)
  gcimsdataset@envir$dt_rt_metrics <- dt_rt_metrics
  gcimsdataset@envir$dt_ref <- dt_ref
  gcimsdataset@envir$rt_ref <- rt_ref

  # Check all have the same dimensions:
  has_different_dt_start <- any(dt_rt_metrics$dt_min_ms != dt_rt_metrics$dt_min_ms[1])
  has_different_dt_end <- any(dt_rt_metrics$dt_max_ms != dt_rt_metrics$dt_max_ms[1])
  has_different_rt_start <- any(dt_rt_metrics$rt_min_s != dt_rt_metrics$rt_min_s[1])
  has_different_rt_end <- any(dt_rt_metrics$rt_max_s != dt_rt_metrics$rt_max_s[1])

  needs_cutting <- any(has_different_dt_start, has_different_dt_end, has_different_rt_start, has_different_rt_end)

  has_different_dt_step <- any(dt_rt_metrics$dt_step_ms != dt_rt_metrics$dt_step_ms[1])
  has_different_rt_step <- any(dt_rt_metrics$rt_step_ms != dt_rt_metrics$rt_step_ms[1])

  needs_interpolate <- any(has_different_dt_step, has_different_rt_step)

  if (needs_interpolate) {
    axesHeterogeneity(gcimsdataset) <- "needs_interpolate"
  } else if (needs_cutting) {
    axesHeterogeneity(gcimsdataset) <- "needs_cutting"
  } else {
    axesHeterogeneity(gcimsdataset) <- "all_equal"
  }
  # fun_aggregate must return the GCIMSDataset object
  gcimsdataset
}


extract_dtime_rtime <- function(object) {
  delayed_op <- GCIMSDelayedOp(
    name = "extract_dtime_rtime",
    fun_extract = .extract_dtime_rtime_fun_extract,
    fun_aggregate = .extract_dtime_rtime_fun_aggregate
  )
  object <- appendDelayedOp(object, delayed_op)
  object
}

"axesHeterogeneity<-" <- function(object, value) {
  valid_values <- c("needs_interpolate", "needs_cutting", "all_equal")
  if (!value %in% valid_values) {
    rlang::abort(
      message = c(
        "Invalid value for axesHeterogeneity",
        "x" = glue("You gave {value}"),
        "i" = glue("Expected one of: {paste0(valid_values, collapse = ', ')}")
        )
    )
  }
  object@envir$axes_heterogeneity <- value
  object
}

axesHeterogeneity <- function(object) {
  object@envir$axes_heterogeneity
}



