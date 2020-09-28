#' Baseline Removal using Psalsa algorithm


#' @param dir_in           The input directory.
#' @param dir_out          The output directory.
#' @param samples          The set of samples to be processed.
#' @param by_rows          Logical. Direction to apply the function. If TRUE it by rows (drift time direction).
#'                         If FALSE, applied by columns (that is the retention time direction).
#' @param lambda           Smoothing parameter (generally 1e5 - 1e8)
#' @param p                Asymmetry parameter
#' @param k                Peak height parameter (usually 5\% of maximum intensity)
#' @return A baseline removed gcims dataset.
#' @family Baseline removal functions
#' @export
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }




gcims_baseline_removal <- function(dir_in, dir_out, samples, by_rows,lambda, p, k){


  print(" ")
  print("  /////////////////////////")
  print(" /    Baseline removal   /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0;
  for (i in (samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    if (by_rows == FALSE){
      aux <- t(aux)
    }
    aux2 <- psalsa(aux,lambda, p, k)
    M <- aux2$corrected

    if (by_rows == FALSE){
      M = t(M)
    }

    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

}
}



#' Estimate spectra baseline using psalsa
#'
#' @param spectra  matrix with one spectrum per row
#' @param lambda smoothing parameter (generally 1e5 - 1e8)
#' @param p      asymmetry parameter
#' @param k      peak height parameter (usually 5\% of maximum intensity)
#' @param maxit  max number of iterations
#' @family Baseline removal functions
#' @export
psalsa <- function(spectra, lambda = 1E7, p = 0.001, k = -1, maxit = 25) {
  if (is.matrix(spectra)) {
    estbaseline <- 0*spectra
    for (i in 1:nrow(spectra)) {
      #estbaseline[i, ] <- psalsa(spectra[i,], lambda, p, k, maxit)
      estbaseline[i, ] <- psalsa_one(spectra[i,], lambda, p, k, maxit)
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

