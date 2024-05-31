#' Align a GCIMS dataset
#'
#' The alignment uses a multiplicative correction in drift time and a
#' Parametric Time Warping correction in retention time
#' @param object A [GCIMSDataset] object, modified in-place
#' @param method_rt Method for alignment, should be "ptw" or "pow"
#' if pow is selected the package "pow" must be installed, to do so visit:
#' https://github.com/sipss/pow
#' @param align_dt if `TRUE` the drift time axis will be aligned using a multiplicative correction
#' @param align_ip if `TRUE` a multiplicative correction will be done in retention time before applying the other algorithm
#' @param reference_sample_idx One number, the index of the sample to use as reference for the alignment in retention time, if NULL the reference will be calculated automatically depending on the method
#' @param ... additional parameters for POW alignment
#' @return The modified [GCIMSDataset]
#' @export
setMethod(
  "align",
  "GCIMSDataset",
  function(object,
           method_rt = "ptw",
           align_dt = TRUE,
           align_ip = TRUE,
           reference_sample_idx = NULL,
           ...) {
    if (align_dt | align_ip){
      tis_matrix <- getTIS(object)
      ric_matrix <- getRIC(object)
      dt <- dtime(object)
      rt <- rtime(object)

      # Optimize drift time alignment parameters:
      rip_position <- apply(tis_matrix, 1L, which.max)
      rip_ref_idx <- round(stats::median(rip_position, na.rm = TRUE))
      rip_ref_ms <- dt[rip_ref_idx]

      # IP alignment parameters
      mins <- apply(ric_matrix, 1L ,which.min)
      rt_ref <- rt[1 : (length(rt) - (max(mins) - min(mins)))]
      min_start <- min(mins) - 1

      delayed_op <- DelayedOperation(
        name = "prealign",
        fun = prealign,
        params = list(align_dt = align_dt,
                      align_ip = align_ip,
                      rip_ref_ms = rip_ref_ms,
                      min_start = min_start,
                      rt_ref = rt_ref),
        fun_extract = .align_fun_extract,
        fun_aggregate = .align_fun_aggregate
      )
      object$appendDelayedOp(delayed_op)
      object$extract_dtime_rtime()
      object$extract_RIC_and_TIS()
    }

    if (method_rt == "ptw" | method_rt == "pow"){
      rt <- rtime(object)
      ric_matrix <- getRIC(object)
      # Optimize ret time alignment parameters:
      if (is.null(reference_sample_idx)){
        if (method_rt == "pow"){
          reference_sample_idx <- pow::select_reference(ric_matrix)
        } else {
          reference_sample_idx <- ptw::bestref(ric_matrix)$best.ref
        }
      }
      # Select reference RIC
      ric_ref <- as.numeric(ric_matrix[reference_sample_idx, ])

      delayed_op <- DelayedOperation(
        name = "align",
        fun = align,
        params = list(method_rt = method_rt,
                      ric_ref = ric_ref,
                      ric_ref_rt = rt,
                      ...),
        fun_extract = .align_fun_extract,
        fun_aggregate = .align_fun_aggregate
      )
      object$appendDelayedOp(delayed_op)
      object$extract_dtime_rtime()
      object$extract_RIC_and_TIS()
    }
    invisible(object)
  }
)

.align_fun_extract <- function(x) {
  list(
    dt_kcorr = x@proc_params$align$dt_kcorr,
    w = x@proc_params$align$w
  )
}

.align_fun_aggregate <- function(ds, extracted_obj) {
  if (is.null(ds$align)) {
    ds$align <- list()
  }
  ds$align[["dt_kcorr"]] <- purrr::map_dbl(extracted_obj, "dt_kcorr")
  ds$align[["w"]] <- purrr::map(extracted_obj, "w")
  ds
}

#' Plots to interpret alignment results
#'
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @return A list with plots created with ggplot2.
#' @export
alignPlots <- function(object) {
  eval_poly <- function(x, coefs) {
    if (is.matrix(coefs)) {
      degree <- ncol(coefs) - 1L
    } else {
      degree <- length(coefs) - 1L
    }
    if (degree == 0L) {
      x_pows <- matrix(1., nrow = 1L, ncol = length(x))
    } else {
      x_pows <- rbind(rep(1, length(x)), t(stats::poly(x, degree = degree,  raw = TRUE, simple = TRUE)))
    }
    coefs %*% x_pows
  }

  get_corrected_dt <- function(dt_before, dt_kcorr) {
    dt_poly_coefs <- matrix(0, nrow = length(dt_kcorr), ncol = 2)
    dt_poly_coefs[,2] <- dt_kcorr
    dt_after <- eval_poly(dt_before, dt_poly_coefs)
    dt_before_mat <- matrix(dt_before, nrow = length(dt_kcorr), ncol = length(dt_before), byrow = TRUE)
    dt_diff <- dt_after - dt_before_mat
    dimnames(dt_diff) <- list(
      SampleID = names(dt_kcorr),
      drift_time_ms = dt_before
    )
    tibble::as_tibble(reshape2::melt(dt_diff, value.name = "correction_ms"))
  }

  get_corrected_rt <- function(rt_before, w) {
    rt_diff <- t(as.data.frame(lapply(w, function(w, rt){(w / seq_along(w)) * rt - rt}, rt_before)))
    colnames(rt_diff) <- rt
    rt_diff <- reshape2::melt(rt_diff, value.name = "correction_s")
    colnames(rt_diff) <- c("SampleID", "ret_time_s", "correction_s")
    rt_diff
  }

  object$realize()
  rt <- rtime(object)
  dt <- dtime(object)

  rt_diff <- get_corrected_rt(rt, object$align$w)
  rt_diff_plot <- ggplot2::ggplot(rt_diff) +
    ggplot2::geom_line(ggplot2::aes(x = .data$ret_time_s, y = .data$correction_s, group = .data$SampleID, color = .data$SampleID)) +
    ggplot2::labs(x = "Retention time (s)", y = "Ret. time correction (s)", color = "SampleID")


  dt_diff <- get_corrected_dt(dt, object$align$dt_kcorr)

  dt_diff_plot <- ggplot2::ggplot(dt_diff) +
    ggplot2::geom_line(ggplot2::aes(x = .data$drift_time_ms, y = .data$correction_ms, group = .data$SampleID, color = .data$SampleID)) +
    ggplot2::labs(x = "Drift time (ms)", y = "Drift time correction (ms)", color = "SampleID")

  dt_kcorr_plot <- ggplot2::ggplot(
    data.frame(
      x = names(object$align$dt_kcorr),
      y = object$align$dt_kcorr
    )
  ) + ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::labs(x = "SampleID", y = "Multiplicative factor (drift time correction, unitless)")

  list(
    rt_diff_plot = rt_diff_plot,
    dt_diff_plot = dt_diff_plot,
    dt_kcorr_plot = dt_kcorr_plot
  )
}
