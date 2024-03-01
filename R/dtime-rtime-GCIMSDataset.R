#' Get A reference drift time vector for the dataset
#' @param object A GCIMSDataset
#' @param sample A number or a string with the sample index or name. If `NULL`, the reference drift time is returned
#' @return a drift time vector
#' @export
setMethod("dtime", "GCIMSDataset", function(object, sample = NULL) {
  if (object$hasDelayedOps() || is.null(object$dt_ref)) {
    object$extract_dtime_rtime()
    object$realize()
  }
  if (is.null(sample)) {
    return(object$dt_ref)
  } else {
    sample <- object$getSample(sample = sample)
    return(dtime(sample))
  }
})

#' Get a reference retention time vector for the dataset
#'
#' @param object A GCIMSDataset
#' @return a retention time vector
#' @param sample A number or a string with the sample index or name. If `NULL`, the reference drift time is returned
#' @export
setMethod("rtime", "GCIMSDataset", function(object, sample = NULL) {
  if (object$hasDelayedOps() || is.null(object$rt_ref)) {
    object$extract_dtime_rtime()
    object$realize()
  }
  if (is.null(sample)) {
    return(object$rt_ref)
  } else {
    sample <- object$getSample(sample = sample)
    return(rtime(sample))
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
  gcimsdataset$dt_ref <- dt_ref
  gcimsdataset$rt_ref <- rt_ref

  # fun_aggregate must return the GCIMSDataset object
  gcimsdataset
}

.extract_inv_k0 <- function(object){
  inv_k0 <- object@inverse_reduced_mobility
  data.frame(
      inv_k0_min = inv_k0[1],
      inv_k0_max = inv_k0[length(inv_k0)],
      inv_k0_step = stats::median(diff(inv_k0))
    )
}

.extract_inv_k0_fun_aggregate <- function(object, extracted_objects){
  inv_k0_samples <- dplyr::bind_rows(extracted_objects, .id = "SampleID")

  max_inv_k0_min <- max(inv_k0_samples$inv_k0_min)
  min_inv_k0_max <- min(inv_k0_samples$inv_k0_max)
  min_inv_k0_step <- min(inv_k0_samples$inv_k0_step)
  max_inv_k0_length <- floor(round((min_inv_k0_max - max_inv_k0_min)/min_inv_k0_step, digits = 8)) + 1L
  inv_k0_ref <- seq(from = max_inv_k0_min, to = min_inv_k0_max, length.out = max_inv_k0_length)

  object$inv_k0_ref <- inv_k0_ref

}

extract_dtime_rtime <- function(object) {
  object$extract_dtime_rtime()
  object
}




