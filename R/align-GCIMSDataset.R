#' Align a GCIMS dataset
#'
#' The alignment uses a multiplicative correction in drift time and a
#' Parametric Time Warping correction in retention time
#' @param object A [GCIMSDataset] object, modified in-place
#' @param method_rt Method for alignment, should be "ptw" or "pow"
#' if pow is selected the package "pow must be installed, to do so visit:
#' https://github.com/sipss/pow
#' @param align_dt if `TRUE` the drift time axis will be aligned using a multiplicative correction
#' @param align_ip if TRUE a multiplicative correction will be done in retention time before applying the other algorithm
#' @param ... additional parameters for POW alignment
#' @return The modified [GCIMSDataset]
#' @export
setMethod(
  "align",
  "GCIMSDataset",
  function(object, method_rt = "ptw", align_dt = TRUE, align_ip = TRUE, ...) {
    tis_matrix <- getTIS(object)
    ric_matrix <- getRIC(object)
    dt <- dtime(object)
    rt <- rtime(object)

    align_params <- alignParams(
      dt = dt,
      rt = rt,
      tis_matrix = tis_matrix,
      ric_matrix = ric_matrix,
      method_rt = method_rt,
      align_dt = align_dt,
      align_ip = align_ip,
      ...)

    delayed_op <- DelayedOperation(
      name = "align",
      fun = align,
      params = align_params,
      fun_extract = .align_fun_extract,
      fun_aggregate = .align_fun_aggregate
    )
    object$appendDelayedOp(delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object$extract_dtime_rtime()
    object$extract_RIC_and_TIS()
    invisible(object)
  }
)

.align_fun_extract <- function(x) {
  list(
    dt_kcorr = x@proc_params$align$dt_kcorr,
    rt_poly_coefs = x@proc_params$align$rt_poly_coefs
  )
}

.align_fun_aggregate <- function(ds, extracted_obj) {
  if (is.null(ds$align)) {
    ds$align <- list()
  }
  ds$align[["dt_kcorr"]] <- purrr::map_dbl(extracted_obj, "dt_kcorr")
  poly_order <- purrr::map_int(
    extracted_obj,
    function(obj) {
      length(obj$rt_poly_coefs) - 1L
    }
  )
  ds$align[["rt_poly_order"]] <- poly_order
  poly_coefs <- matrix(0.0, nrow = length(extracted_obj), ncol = max(poly_order) + 1L)
  poly_coefs[,2] <- 1
  for (i in  seq_along(extracted_obj)) {
    coefs <- extracted_obj[[i]]$rt_poly_coefs
    poly_coefs[i, seq_along(coefs)] <- coefs
  }
  dimnames(poly_coefs) <- list(
    SampleID = sampleNames(ds),
    PolyOrder = paste0("Order_", seq_len(ncol(poly_coefs)) - 1L)
  )
  ds$align[["rt_poly_coefs"]] <- poly_coefs
  ds
}


alignParams <- function(dt, rt, tis_matrix, ric_matrix, method_rt, align_dt, align_ip, ...) {
  # Optimize ret time alignment parameters:
  if (method_rt == "ptw"){
    ref_ric_sample_idx <- ptw::bestref(ric_matrix)$best.ref
  } else {
    ref_ric_sample_idx <- pow::select_reference(ric_matrix)
  }
  # Select reference RIC
  ric_ref <- as.numeric(ric_matrix[ref_ric_sample_idx, ])

  # Optimize drift time alignment parameters:
  rip_position <- apply(tis_matrix, 1L, which.max)
  rip_ref_idx <- round(stats::median(rip_position, na.rm = TRUE))
  rip_ref_ms <- dt[rip_ref_idx]

  # IP alignment parameters

  ip_position <- apply(ric_matrix, 1L ,which.min)
  ip_ref_idx <- round(stats::median(ip_position, na.rm = TRUE))
  ip_ref_s <- rt[ip_ref_idx]

  list(rip_ref_ms = rip_ref_ms,
       ric_ref = ric_ref,
       ric_ref_rt = rt,
       ip_ref_s = ip_ref_s,
       method_rt = method_rt,
       align_dt = align_dt,
       align_ip = align_ip,
       ...)
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

  get_corrected_rt <- function(rt_before, rt_poly_coefs) {
    rt_after <- eval_poly(rt_before, rt_poly_coefs)
    rt_before_mat <- matrix(rt_before, nrow = nrow(rt_poly_coefs), ncol = length(rt_before), byrow = TRUE)
    rt_diff <- rt_after - rt_before_mat
    dimnames(rt_diff) <- list(
      SampleID = rownames(rt_poly_coefs),
      ret_time_s = rt_before
    )
    tibble::as_tibble(reshape2::melt(rt_diff, value.name = "correction_s"))
  }

  object$realize()
  rt <- rtime(object)
  dt <- dtime(object)

  rt_diff <- get_corrected_rt(rt, object$align$rt_poly_coefs)

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
  ) + ggplot2::geom_col(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::labs(x = "SampleID", y = "Multiplicative factor (drift time correction, unitless)") +
    ggplot2::coord_flip()


  list(
    rt_diff_plot = rt_diff_plot,
    dt_diff_plot = dt_diff_plot,
    dt_kcorr_plot = dt_kcorr_plot
  )
}
