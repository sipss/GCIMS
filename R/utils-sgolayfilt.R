convolve_prepare <- function(x, conj = FALSE, plan = NULL, impl = "auto") {
  if (conj) {
    do_conj <- Conj
  } else {
    do_conj <- identity
  }
  if (!is.matrix(x)) {
    do_conj(stats::fft(x))
  } else {
    do_conj(stats::mvfft(x))
  }
}

convolve_do <- function(fft_x, conj_fft_y) {
  if (!is.matrix(fft_x)) {
    Re(stats::fft(fft_x * conj_fft_y, inverse = TRUE))/length(fft_x)
  } else {
    Re(stats::mvfft(fft_x * conj_fft_y, inverse = TRUE))/nrow(fft_x)
  }
}


sgolayfilt_impl <- function(x, filt, rowwise, return_matrix, engine = "fft") {
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

  coefs_for_first_points <- filt[1:k, , drop = FALSE]
  coefs_for_last_points <- filt[(k + 2L):n, , drop = FALSE]

  if (rowwise) {
    out[, 1L:k] <- x[, 1L:n, drop = FALSE] %*% t(coefs_for_first_points)
    out[, (len - k + 1L):len] <- x[, (len - n + 1L):len, drop = FALSE] %*% t(coefs_for_last_points)
  } else {
    out[1L:k, ] <- coefs_for_first_points %*% x[1L:n, , drop = FALSE]
    out[(len - k + 1L):len,] <- coefs_for_last_points %*% x[(len - n + 1L):len, , drop = FALSE]
  }

  if (engine == "fft") {
    conv_coefs <- filt[k + 1L, n:1L]
    fft_length <- length(conv_coefs) + len - 1L
    conv_coefs_padded <- c(rev(conv_coefs), rep(0, len - 1L))
    conj_fft_y_prep <- convolve_prepare(conv_coefs_padded, conj = TRUE)
    x_padded <- numeric(fft_length)
    center_points_idx <- (k + 1L):(len - k)
    for (i in seq_len(num_ser)) {
      if (rowwise) {
        x_padded[length(conv_coefs):fft_length] <- x[i,]
      } else {
        x_padded[length(conv_coefs):fft_length] <- x[,i]
      }
      fft_x_prep <- convolve_prepare(x_padded, conj = FALSE)
      center_points <- convolve_do(
        fft_x_prep,
        conj_fft_y_prep
      )[n:len]
      if (rowwise) {
        out[i, center_points_idx] <- center_points
      } else {
        out[center_points_idx, i] <- center_points
      }
    }
  } else if (engine == "fft2") {
    conv_coefs <- filt[k + 1L, n:1L]
    fft_length <- length(conv_coefs) + len - 1L
    conv_coefs_padded <- c(rev(conv_coefs), rep(0, len - 1L))
    conj_fft_y_prep <- convolve_prepare(conv_coefs_padded, conj = TRUE)
    xmat_padded <- matrix(0, nrow = fft_length, ncol = num_ser)
    center_points_idx <- (k + 1L):(len - k)
    if (rowwise) {
      xmat_padded[length(conv_coefs):fft_length,] <- t(x)
    } else {
      xmat_padded[length(conv_coefs):fft_length,] <- x
    }
    fft_x_prep <- convolve_prepare(xmat_padded, conj = FALSE)
    if (rowwise) {
      out[, center_points_idx] <- t(convolve_do(
        fft_x_prep,
        conj_fft_y_prep
      )[n:len,,drop = FALSE])
    } else {
      out[center_points_idx, ] <- convolve_do(
        fft_x_prep,
        conj_fft_y_prep
      )[n:len,,drop = FALSE]
    }


  } else if (engine == "filter") {
    center_points_idx <- (k + 1L):(len - k)
    filt_cent <- filt[k + 1L, n:1L]
    embedded_signal <- matrix(0, nrow = len - 2*k, ncol = n)
    for (i in seq_len(num_ser)) {
      if (rowwise) {
        for (j in n:1L) {
          first_idx <- (n - j + 1L)
          last_idx <- first_idx + len - 2*k - 1
          embedded_signal[,j] <- x[i,first_idx:last_idx]
        }
        out[i, center_points_idx] <- embedded_signal %*% filt_cent
      } else {
        for (j in n:1L) {
          first_idx <- (n - j + 1L)
          last_idx <- first_idx + len - 2*k - 1
          embedded_signal[,j] <- x[first_idx:last_idx, i]
        }
        out[center_points_idx, i] <- embedded_signal %*% filt_cent
      }
    }
  } else {
    stop("Wrong engine. Use fft or filter")
  }
  if (!return_matrix) {
    attr(out, "dim") <- NULL
  }
  out
}


choose_engine <- function(x, filter_length, orig_engine) {
  engine <- orig_engine
  if (engine == "filter") {
    return(engine)
  }
  if (engine == "auto") {
    if (is.matrix(x)) {
      if (filter_length > 100) {
        engine <- "fft"
      } else {
        engine <- "filter"
      }
    }
  }
  if (engine == "fft" && anyNA(x)) {
    if (orig_engine == "fft") {
      warn(
        message = c(
          'Switching sgolayfilt engine from "fft" to "filter"',
          "!" = "The fft engine does not handle missing values.",
          "i" = "Using engine = 'filter' instead."
        )
      )
    }
    engine <- "filter"
  }
  engine
}

#' Savitzky-Golay filtering
#'
#' @noRd
#' @param x A numeric matrix or vector
#' @inheritParams signal::sgolayfilt
#' @param rowwise If `TRUE`, Apply the filter by rows instead of by columns
#' @param engine "fft" Uses the Fast Fourier Transform to apply the filter. "filter" uses a convolution. Both give
#' the same results, fft is usually faster on larger filter lengths, and larger matrices.
#'
#' @return A matrix or vector of the same dimensions or length as `x`, with the result of the filter
#'
#' @examples
#' x <- runif(300)
#' y <- sgolayfilt(x, p=2, n = 21)
sgolayfilt <- function(x, p = 3, n = p + 3 - p %% 2, m = 0, ts = 1, rowwise = FALSE,
                       engine = c("auto", "fft", "fft2", "filter")) {
  engine <- match.arg(engine)
  if (inherits(p, "sgolayFilter") || (!is.null(dim(p)) && dim(p) > 1)) {
    filt <- p
  } else {
    filt <- signal::sgolay(p, n, m, ts)
  }
  return_matrix <- TRUE
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    return_matrix <- FALSE
    rowwise <- FALSE
  }
  engine <- choose_engine(x = x, filter_length = nrow(filt), orig_engine = engine)
  out <- sgolayfilt_impl(x, filt, rowwise, return_matrix, engine = engine)
  out
}
