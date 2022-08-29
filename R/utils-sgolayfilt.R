fft_plan <- function(n, effort = 0) {
  if (requireNamespace("fftw", quietly = TRUE)) {
    fftw::planFFT(n, effort = effort)
  } else {
    NULL
  }
}

fft <- function(x, inverse = FALSE, plan = NULL) {
  if (requireNamespace("fftw", quietly = TRUE)) {
    if (!is.null(plan)) {
      fftw::FFT(x, inverse = inverse, plan = plan)
    } else {
      fftw::FFT(x, inverse = inverse)
    }
  } else {
    stats::fft(x, inverse = inverse)
  }
}

convolve_prep <- function(fft_x, conj_fft_y, plan) {
  n <- length(fft_x)
  x <- fft(fft_x * conj_fft_y, inverse = TRUE, plan = plan)
  Re(x)/n
}

#' Savitzky-Golay filtering
#'
#' @param x A numeric matrix or vector
#' @inheritParams signal::sgolayfilt
#' @param rowwise If `TRUE`, Apply the filter by rows instead of by columns
#' @param engine "fft" Uses the Fast Fourier Transform to apply the filter. "signal" uses `signal::sgolayfilt`. Both give
#' the same results, fft is usually faster on matrices.
#'
#' @return A matrix or vector of the same dimensions or length as `x`, with the result of the filter
#' @export
#'
#' @examples
#' x <- runif(300)
#' y <- sgolayfilt(x, p=2, n = 21)
sgolayfilt <- function(x, p = 3, n = p + 3 - p %% 2, m = 0, ts = 1, rowwise = FALSE, engine = c("auto", "fft", "signal")) {
  engine <- match.arg(engine)
  if (inherits(p, "sgolayFilter") || (!is.null(dim(p)) && dim(p) > 1)) {
    filt <- p
  } else {
    filt <- signal::sgolay(p, n, m, ts)
  }
  return_matrix <- TRUE

  orig_engine <- engine
  if (engine == "auto") {
    if (is.matrix(x) && requireNamespace("fftw", quietly = TRUE)) {
      engine <- "fft"
    } else {
      engine <- "signal"
    }
  } else if (engine == "fft" && anyNA(x)) {
    if (orig_engine == "fft") {
      rlang::warn("The fft engine does not handle missing values. Using engine = 'signal' instead")
    }
    engine <- "signal"
  }
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    return_matrix <- FALSE
    rowwise <- FALSE
  }
  if (engine == "signal") {
    if (rowwise) {
      out <- t(apply(x, 1, function(xvec) signal::sgolayfilt(xvec, filt)))
    } else {
      out <- apply(x, 2L, function(xvec) signal::sgolayfilt(xvec, filt))
    }
    if (!return_matrix) {
      attr(out, "dim") <- NULL
    }
    return(out)
  } else if (engine != "fft") {
    rlang::abort(glue("Unknown engine {engine}"))
  }
  # fft engine:

  if (rowwise) {
    num_ser <- nrow(x)
    len <- ncol(x)
  } else {
    num_ser <- ncol(x)
    len <- nrow(x)
  }

  n <- nrow(filt)
  k <- floor(n/2)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  conv_coefs <- filt[k + 1L, n:1L]
  #len_pow2 <- stats::nextn(len, 2)
  fft_length <- length(conv_coefs) + len - 1L
  conv_coefs_padded <- c(rev(conv_coefs), rep(0, len - 1L))
  # Roughly at commit e357073d727e811396fd84269cbb64c71058c8c1
  # The vignette took 5.1 minutes to build with effort = 3 and
  # 3.2 minutes to build with effort = 0. So we just hardcode effort=0.
  # FIXME: Maybe effort=1 or effort=2 is better, but I'll check later
  plan <- fft_plan(fft_length, effort = 0)
  conj_fft_y_prep <- Conj(fft(conv_coefs_padded, plan = plan))
  x_padded <- numeric(fft_length)
  for (i in seq_len(num_ser)) {
    if (rowwise) {
      xvec <- x[i,]
    } else {
      xvec <- x[,i]
    }
    first_points <- filt[1:k, ] %*% xvec[1:n]
    x_padded[length(conv_coefs):fft_length] <- xvec
    fft_x_prep <- fft(x_padded, plan = plan)
    center_points <- convolve_prep(
      fft_x_prep,
      conj_fft_y_prep,
      plan = plan
    )[n:len]
    last_points <- filt[(k + 2L):n, ] %*% xvec[(len - n + 1L):len]
    if (rowwise) {
      out[i,] <- c(first_points, center_points, last_points)
    } else {
      out[,i] <- c(first_points, center_points, last_points)
    }
  }
  if (!return_matrix) {
    attr(out, "dim") <- NULL
  }
  return(out)
  out
}

# This code can be used to benchmark the performance.
# We could add some plots to this, to easily show the time and memory requirements
# of each engine as a function of n, p and the matrix dimensions.
# engine = "auto" outperforms the alternatives:
# filt <- signal::sgolay(p = 2, n = 5, m = 0, ts = 1)
# x <- matrix(runif(4000*4000), ncol = 4000)
# bm_cols <- bench::mark(
#   signal = {mysgolayfilt(x, filt, engine = "signal")},
#   me = {mysgolayfilt(x, filt, engine = "auto")}
# )
# print(bm_cols)
# filt <- signal::sgolay(p = 2, n = 25, m = 0, ts = 1)
# x <- matrix(runif(4000*2), ncol = 2)
# bm_cols <- bench::mark(
#   signal = {mysgolayfilt(x, filt, engine = "signal")},
#   me = {mysgolayfilt(x, filt, engine = "auto")}
# )
# print(bm_cols)
#
#
#
# filt <- signal::sgolay(p = 2, n = 25, m = 0, ts = 1)
# x <- matrix(runif(4000*4000), ncol = 4000)
# bm_cols <- bench::mark(
#   signal = {mysgolayfilt(x, filt, engine = "signal")},
#   me = {mysgolayfilt(x, filt, engine = "auto")}
# )
# print(bm_cols)
#
#
# bm_rows <- bench::mark(
#   signal = {mysgolayfilt(x, filt, engine = "signal", rowwise = TRUE)},
#   me = {mysgolayfilt(x, filt, engine = "auto", rowwise = TRUE)}
# )
# print(bm_rows)
#
# x <- runif(1E6)
# bm_num <- bench::mark(
#   signal = {mysgolayfilt(x, filt, engine = "signal")},
#   me = {mysgolayfilt(x, filt, engine = "auto")}
# )
# print(bm_num)
