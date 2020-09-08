#' Estimate spectra baseline using psalsa
#'
#' @param spectra  matrix with one spectrum per row
#' @param lambda smoothing parameter (generally 1e5 - 1e8)
#' @param p      asymmetry parameter
#' @param k      peak height parameter (usually 5\% of maximum intensity)
#' @param maxit  max number of iterations
#'
#' @export
psalsa <- function(spectra, lambda = 1E7, p = 0.001, k = -1, maxit = 25) {
  if (is.matrix(spectra)) {
    estbaseline <- 0*spectra
    for (i in 1:nrow(spectra)) {
      estbaseline[i, ] <- psalsa(spectra[i,], lambda, p, k, maxit)
    }
  } else {
    estbaseline <- psalsa_one(spectra, lambda, p, k, maxit)
  }
  return (list(baseline = estbaseline, corrected = spectra - estbaseline))
}


psalsa_one <- function(y, lambda = 1e+07, p = 0.001, k = -1, maxit = 25) {
  z <- 0 * y
  w <- z + 1
  if (k == -1) {
    k <- max(y)/20
  }

  dimensions <- dim(y)
  if (is.null(dimensions)) {
    num_samples <- 1
    num_points <- length(y)
  } else {
    num_samples <- dimensions[1]
    num_points <- dimensions[2]
  }
  d <- 0*numeric(num_points) - 1

  for (it in 1:maxit) {
    z <- ptw::whit2(y, lambda, w)

    d_geq_old = (d >= 0)
    d <- y - z
    w[d >= 0] = p*exp(-d[d >= 0]/k)
    w[d < 0]  = 1 - p;
    if (all((d >= 0) == d_geq_old)) {
      break
    }

  }
  return(z)
}

