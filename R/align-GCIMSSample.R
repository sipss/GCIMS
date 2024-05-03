#' Align a GCIMSSample object, in drift and retention time
#' @param object A [GCIMSSample] object
#' @param rip_ref_ms The reference position of the Reactant Ion Peak in the dataset (in ms)
#' @param ric_ref The reference Reverse Ion Chromatogram
#' @param ric_ref_rt The retention times corresponding to `ric_ref`
#' @param min_start minimun injection point, to calculate where to begin the spectrums and cut as few points as posible, to be used in injection point alignment
#' @param rt_ref retention time reference for alignment to injection point
#' @param method_rt Method for alignment, should be "ptw" or "pow"
#' @param align_dt if `TRUE`, align the drift time axis using a multiplicative correction
#' @param align_ip if TRUE a multiplicative correction will be done in retention time before applying the other algorithm
#' @param ploynomial_order_ptw One number, by default 2, the maximum order of the polynomial for the parametric time warping alignment.
#' @param ... Additional arguments passed on to the alignment method.
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "align", "GCIMSSample",
  function(object, rip_ref_ms, ric_ref, ric_ref_rt, min_start, rt_ref, method_rt, align_dt, align_ip, ploynomial_order_ptw, ...){
    if (all(is.na(object@data))) {
      cli_abort("All the data matrix of {description(object)} are missing values. Align is impossible")
    }

    # Align in drift time
    if (align_dt) {
      object <- alignDt(object, rip_ref_ms = rip_ref_ms)
    } else{
      object@proc_params$align$dt_kcorr <- 1
    }

    if (all(is.na(object@data))) {
      cli_abort("After aligning drift times, all the data matrix of {description(object)} are missing values. This should not happen")
    }

    # Align in retention time
    if (align_ip){
      object <- alignRt_ip(object, min_start = min_start, rt_ref = rt_ref)
    }
    if (method_rt == "ptw") {
      object <- alignRt_ptw(object, ric_ref = ric_ref, ric_ref_rt = ric_ref_rt, ploynomial_order_ptw)
    } else if(method_rt == "pow"){
      object <- alignRt_pow(object, ric_ref = ric_ref, ric_ref_rt = ric_ref_rt, ...)
    } else {
      if (!"align" %in% names(object@proc_params)) {
        object@proc_params$align <- list()
      }
      object@proc_params$align$rt_min_s <- ric_ref_rt[1]
      object@proc_params$align$rt_max_s <- ric_ref_rt[length(ric_ref_rt)]
      object@proc_params$align$rt_poly_coefs <- c(0, 1)
    }


    if (all(is.na(object@data))) {
      cli_abort("After aligning drift and retention times, all the data matrix of {description(object)} are missing values. This should not happen")
    }

    dt_range <- c(
      object@proc_params$align$dt_min_ms,
      object@proc_params$align$dt_max_ms
    )
    rt_range <- c(
      object@proc_params$align$rt_min_s,
      object@proc_params$align$rt_max_s
    )
    object <- subset(object, dt_range = dt_range, rt_range = rt_range)
    if (all(is.na(object@data))) {
      cli_abort("After aligning and subsetting all the data matrix of {description(object)} are missing values. This should not happen")
    }
    object
  })

#' Align a GCIMSSample in drift time with a multiplicative correction
#' @param x A [GCIMSSample] object
#' @param rip_ref_ms The position of the RIP in ms
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "alignDt", "GCIMSSample",
  function(x, rip_ref_ms) {
    tis <- getTIS(x)
    dt <- dtime(x)
    rip_pos_ms <- dt[which.max(tis)]
    Kcorr <- rip_ref_ms/rip_pos_ms

    int_mat <- intensity(x)
    dt_corr <- Kcorr * dt
    for (j in seq_len(ncol(int_mat))) {
      int_mat[, j] <- signal::interp1(dt_corr, int_mat[, j], dt, extrap = TRUE)
    }
    intensity(x) <- int_mat

    # Save kcorr and non-extrapolated dt range:
    dt_beg <- dt[1]
    dt_end <- dt[length(dt)]
    dt_corr_beg <- dt_corr[1]
    dt_corr_end <- dt_corr[length(dt_corr)]

    if (dt_corr_beg > dt_beg) {
      dt_min_ms <- dt[dt >= dt_corr_beg][1]
    } else {
      dt_min_ms <- dt_beg
    }

    if (dt_corr_end < dt_end) {
      dt_max_ms <- utils::tail(dt[dt <= dt_corr_end], 1L)
    } else {
      dt_max_ms <- dt_end
    }


    if (!"align" %in% names(x@proc_params)) {
      x@proc_params$align <- list()
    }
    x@proc_params$align$dt_kcorr <- Kcorr
    x@proc_params$align$dt_min_ms <- dt_min_ms
    x@proc_params$align$dt_max_ms <- dt_max_ms
    # Return the sample
    x
  })


#' Align a GCIMSSample in retention time with a multiplicative correction
#' @param x A [GCIMSSample] object
#' @param min_start minimun injection point, to calculate where to begin the spectrums and cut as few points as posible
#' @param rt_ref retention time reference
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "alignRt_ip","GCIMSSample",
  function(x, min_start, rt_ref) {
    ric <- getRIC(x)
    injection_point <- which.min(ric)
    x@retention_time <- rt_ref
    x@data <- x@data[, (injection_point - min_start):((injection_point - min_start)+length(rt_ref)-1)]
    x
  })

#' Align a GCIMSSample in retention time with parametric optimized warping
#' @param x A [GCIMSSample] object
#' @param ric_ref The reference Reverse Ion Chromatogram
#' @param ric_ref_rt The retention times corresponding to `ric_ref`
#' @param lambdas a vector with the penalties to test the POW
#' @param p By default `10`, meaning to use one every `10` points to validate.
#' @param max_it Maximum number of iterations
#' @param lambda1 Regularization parameter for second derivative of warp
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "alignRt_pow", "GCIMSSample",
  function(x,
           ric_ref,
           ric_ref_rt,
           lambdas = pracma::logspace(-2, 4, 31),
           p = 10,
           max_it = 5000,
           lambda1 = 10^6) {
    ric <- getRIC(x)
    v <- rep(1, length(ric_ref))
    iv <- seq(2, length(ric_ref) - 1, by = p)
    v[iv] <- 0L
    W <- Matrix::Diagonal(x = v)
    result_val <- pow:::val(ric, ric_ref, W, iv, lambdas, fom ="rms", lambda1 = lambda1)
    e_ix <- result_val$e
    ti_ix <- result_val$ti
    e_ix[ti_ix == 1] <- NA
    lambdas[ti_ix == 1] <- NA
    best_lambda <- lambdas[which.min(e_ix)]
    w <- pow::pow(ric, best_lambda, ric_ref, max_it = max_it, lambda1= lambda1)

    int <- intensity(x)
    inter <- t(apply(int, 1, pow::interpolation, w = w, return = FALSE))
    sel <- !is.na(inter[1,])
    x@retention_time <- ric_ref_rt[sel]
    x@data <- inter[,sel]


    if (!"align" %in% names(x@proc_params)) {
      x@proc_params$align <- list()
    }
    x@proc_params$align$warp <- w
    #x@proc_params$align$rt_min_s <- x@retention_time[1]
    #x@proc_params$align$rt_max_s <- x@retention_time[length(x@retention_time)]
    x@proc_params$align$rt_poly_coefs <- c(0, 1)

    return(x)
  })

#' Align a GCIMSSample in retention time using parametric time warping
#' @param x A [GCIMSSample] object
#' @param ric_ref The reference Reverse Ion Chromatogram
#' @param ric_ref_rt The retention times corresponding to `ric_ref`
#' @param ploynomial_order_ptw maximum order of the polynomial to be used
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "alignRt_ptw", "GCIMSSample",
  function(x, ric_ref, ric_ref_rt, ploynomial_order_ptw) {
    optimize_polynomial_order <- function(ric_sample, ric_ref) {
      correction_type_options <- seq.int(0, ploynomial_order_ptw)
      poly_orders <- seq.int(1, ploynomial_order_ptw)
      xi <- seq_len(length(ric_sample))
      corr <- numeric(length(poly_orders) + 1L)
      corr[1] <- stats::cor(ric_ref, ric_sample, use = "complete.obs")
      for (j in seq_along(poly_orders)) {
        poly <- rep(0, poly_orders[j])
        poly[2] <- 1
        corr[j + 1L] <- stats::cor(
          ric_ref,
          as.numeric(ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = poly)$warped.sample[1,]),
          use = "complete.obs"
        )
      }

      ### RIGHT NOW IF THIS IS UNCOMMENTED 1ST ORDER POLYNOMIAL IS GOING TO BE SELECTED
      # # Initialize index:
      # idx_max <- idx_sel <- idx_zero <- idx_sign <- length(poly_orders) + 1L
      # # Check when correlation decreases for the first time or does not change when increasing degree of the polynomial.
      # diff_corr_i <- diff(corr)
      # if (all(sign(diff_corr_i) == 1)) {
      #   # Correlation is always increasing (don't do anything)
      # } else if (any(diff_corr_i == 0)) {
      #   # Correlation is equal for at least for  different degrees of the polynomial.
      #   idx_zero <- which(diff_corr_i == 0)[1]
      # } else if (any(sign(diff_corr_i) == -1)) {
      #   # Correlation decreases at least for two consecutive degrees of the polynomial.
      #   idx_sign <- which(sign(diff_corr_i) == -1)[1]
      # }
      # # Combine indexes and look for the minimum
      # idx_combine <- c(idx_sign, idx_zero)
      # if (any(idx_combine < idx_max)) {
      #   # Select the minimum index in which the correlation is still increasing
      #   idx_sel <- min(idx_combine)
      # }
      correction_type_options[which.max(corr)]
    }

    # if (missing(ric_ref) && !missing(y)) {
    #   cli_warn("Please provide ric_ref as a named argument")
    #   ric_ref <- y
    # }

    ric <- getRIC(x)
    rt <- rtime(x)

    if (identical(rt, ric_ref_rt)) {
      ric_ref_interp <- ric_ref
    } else {
      ric_ref_interp <- signal::interp1(x = ric_ref_rt, y = ric_ref, xi = rt, extrap = FALSE)
    }

    poly_order <- optimize_polynomial_order(ric, ric_ref_interp)

    if (!"align" %in% names(x@proc_params)) {
      x@proc_params$align <- list()
    }


    if (poly_order == 0L) {
      x@proc_params$align$rt_poly_order <- 2
      x@proc_params$align$rt_poly_coefs <- c(0, 1)
      x@proc_params$align$rt_min_s <- rt[1]
      x@proc_params$align$rt_max_s <- rt[length(rt)]
      return(x)
    }
    # Correct retention time axis using the reference RIC and sample RIC.
    poly <- rep(0, poly_order + 1L)
    poly[2] <- 1
    align_result <- ptw::ptw(ref = ric_ref_interp, samp = ric, init.coef = poly)
    rt_corr_idx <- align_result$warp.fun[1,]
    # Interpolate data using the correction
    int_mat <- intensity(x)
    rt_idx <- seq_len(ncol(int_mat))
    for (i in seq_len(nrow(int_mat))) {
      int_mat[i,] <- signal::interp1(x = rt_corr_idx, y = int_mat[i,], xi = rt_idx, extrap = FALSE)
    }

    intensity(x) <- int_mat

    rt_beg <- rt_idx[1]
    rt_end <- rt_idx[length(rt_idx)]
    rt_corr_beg <- rt_corr_idx[1]
    rt_corr_end <- rt_corr_idx[length(rt_corr_idx)]

    if (rt_corr_beg > rt_beg) {
      rt_min_idx <- rt_idx[rt_idx >= rt_corr_beg][1]
    } else {
      rt_min_idx <- rt_beg
    }
    rt_min_s <- rt[rt_min_idx]


    if (rt_corr_end < rt_end) {
      rt_max_idx <- utils::tail(rt_idx[rt_idx <= rt_corr_end], 1L)
    } else {
      rt_max_idx <- rt_end
    }
    rt_max_s <- rt[rt_max_idx]

    x@proc_params$align$rt_poly_order <- poly_order
    x@proc_params$align$rt_poly_coefs <- as.numeric(align_result$warp.coef[1L,])
    x@proc_params$align$rt_min_s <- rt_min_s
    x@proc_params$align$rt_max_s <- rt_max_s

    x
  })
