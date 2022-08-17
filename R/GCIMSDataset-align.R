#' Align a GCIMS dataset
#'
#' The alignment uses a multiplicative correction in drift time and a
#' Parametric Time Warping correction in retention time
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @return The modified [GCIMSDataset]
#' @export
setMethod(
  "align",
  "GCIMSDataset",
  function(object) {
    object <- extract_dtime_rtime(object)
    tis_matrix <- getTIS(object)
    ric_matrix <- getRIC(object)


    # Optimize ret time alignment parameters:
    ref_ric_sample_idx <- find_reference_ric(ric_matrix)
    # Select reference RIC
    ric_ref_rt <- rtime(object)
    ric_ref <- as.numeric(ric_matrix[ref_ric_sample_idx, ])

    # Optimize drift time alignment parameters:
    rip_position <- apply(tis_matrix, 1L, which.max)
    rip_ref_idx <- round(stats::median(rip_position, na.rm = TRUE))
    dt <- dtime(object)
    rip_ref_ms <- dt[rip_ref_idx]

    delayed_op <- GCIMSDelayedOp(
      name = "align",
      fun = align,
      params = list(rip_ref_ms = rip_ref_ms, ric_ref = ric_ref, ric_ref_rt = ric_ref_rt),
      fun_extract = function(x) {
        list(
          dt_kcorr = x@proc_params$align$dt_kcorr,
          rt_poly_coefs = x@proc_params$align$rt_poly_coefs
        )
      },
      fun_aggregate = function(ds, extracted_obj) {
        if (!"align" %in% names(ds@envir)) {
          ds@envir$align <- list()
        }
        ds@envir$align[["dt_kcorr"]] <- purrr::map_dbl(extracted_obj, "dt_kcorr")
        poly_order <- purrr::map_int(
          extracted_obj,
          function(obj) {
            length(obj$rt_poly_coefs) - 1L
          }
        )
        ds@envir$align[["rt_poly_order"]] <- poly_order
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
        ds@envir$align[["rt_poly_coefs"]] <- poly_coefs
        ds
      }
    )
    object <- appendDelayedOp(object, delayed_op)

    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

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
      x_pows <- rbind(rep(1, length(x)), t(poly(x, degree = degree,  raw = TRUE, simple = TRUE)))
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

  object <- realize(object)
  rt <- rtime(object)
  dt <- dtime(object)

  rt_diff <- get_corrected_rt(rt, object@envir$align$rt_poly_coefs)

  rt_diff_plot <- ggplot2::ggplot(rt_diff) +
    ggplot2::geom_line(ggplot2::aes(x = .data$ret_time_s, y = .data$correction_s, group = .data$SampleID, color = .data$SampleID)) +
    ggplot2::labs(x = "Retention time (s)", y = "Ret. time correction (s)", color = "SampleID")


  dt_diff <- get_corrected_dt(dt, object@envir$align$dt_kcorr)

  dt_diff_plot <- ggplot2::ggplot(dt_diff) +
    ggplot2::geom_line(ggplot2::aes(x = .data$drift_time_ms, y = .data$correction_ms, group = .data$SampleID, color = .data$SampleID)) +
    ggplot2::labs(x = "Drift time (ms)", y = "Drift time correction (ms)", color = "SampleID")

  dt_kcorr_plot <- ggplot2::qplot(x = names(object@envir$align$dt_kcorr), y = object@envir$align$dt_kcorr, geom = "col") +
    labs(x = "SampleID", y = "Multiplicative factor (drift time correction, unitless)") +
    ggplot2::coord_flip()


  list(
    rt_diff_plot = rt_diff_plot,
    dt_diff_plot = dt_diff_plot,
    dt_kcorr_plot = dt_kcorr_plot
  )
}
