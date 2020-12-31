#' Remove baseline of a dataset using psalsa
#'
#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param time            It indicates if the correction is going to be in the
#'                        drift time or in the retention time. It should be
#'                        introduce "Retention" for correcting the retention
#'                        time; or "Drift" for the drift time.
#' @param lambda          Smoothing parameter (generally 1e5 - 1e8)
#' @param p               Asymmetry parameter
#' @param k               Peak height parameter (usually 5\% of maximum intensity)
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
#'
gcims_baseline_removal <- function(dir_in, dir_out, samples,
                                   by_rows,lambda, p, k){

  print(" ")
  print("  /////////////////////////")
  print(" /    Baseline removal   /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  m = -1;
  for (i in c(0, samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)

    if (time == "Drift"){
      aux <- t(aux)
    } else if (time == "Retention"){
    }

    psalsa_results <- psalsa(aux,lambda, p, k)
    aux  <- psalsa_results$corrected

    if (time == "Drift"){
      aux <- t(aux)
    } else if (time == "Retention"){
    }

    aux_list$data$data_df <- aux
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

}
}






#' Estimate the baseline for all rows in a matrix using psalsa
#'
#' @param data            Matrix with one spectrum (or chromatogram) per row.
#' @param lambda          Smoothing parameter (generally 1e5 - 1e8).
#' @param p               Asymmetry parameter.
#' @param k               Peak height parameter (usually 5\% of maximum
#'                        intensity).
#' @param maxit           Max number of iterations.
#' @return a list with:
#'     - baseline:  A matrix with all the estimated baselines.
#'     - corrected: A matrix of corrected data.
#'
#' @family Baseline removal functions.
#' @export
#'
#' @references {
#'
#' Oller-Moreno, S., Pardo, A., Jiménez-Soto, J. M., Samitier, J., & Marco, S. (2014, February).
#' Adaptive Asymmetric Least Squares baseline estimation for analytical instruments.
#' In 2014 IEEE 11th International Multi-Conference on Systems, Signals & Devices (SSD14) (pp. 1-5). IEEE.
#'
#'  }
#'
#' @author Sergio Oller-Moreno,  \email{soller@@.ibecbarcelona.eu}
#'
#'
#' Institute for Bioengieering of Catalonia. Signal and Information Processing for Sensing Systems group. Barcelona, Spain.
#'
psalsa <- function(data, lambda = 1E7, p = 0.001, k = -1, maxit = 25) {
  if (is.matrix(data)) {
    estbaseline <- 0*data
    for (i in 1:nrow(data)) {
      #estbaseline[i, ] <- psalsa(spectra[i,], lambda, p, k, maxit)
      estbaseline[i, ] <- psalsa_one(data[i,], lambda, p, k, maxit)
    }
  } else {
    estbaseline <- psalsa_one(data, lambda, p, k, maxit)
  }
  return (list(baseline = estbaseline, corrected = data - estbaseline))

}








#' Estimate baseline of a curve using psalsa
#'
#' @param y               Either a chromatogram or a spectrum
#' @param lambda          Smoothing parameter (generally 1e5 - 1e8)
#' @param p               Asymmetry parameter
#' @param k               Peak height parameter (usually 5\% of maximum
#'                        intensity)
#' @param maxit           Max number of iterations
#' @return The estimated baseline
#' @family Baseline removal functions
#' @export
#'
#' @references {
#'
#' Oller-Moreno, S., Pardo, A., Jiménez-Soto, J. M., Samitier, J., & Marco, S. (2014, February).
#' Adaptive Asymmetric Least Squares baseline estimation for analytical instruments.
#' In 2014 IEEE 11th International Multi-Conference on Systems, Signals & Devices (SSD14) (pp. 1-5). IEEE.
#'
#'  }
#'
#' @author Sergio Oller-Moreno,  \email{soller@@.ibecbarcelona.eu}
#'
#'
#' Institute for Bioengieering of Catalonia. Signal and Information Processing for Sensing Systems group. Barcelona, Spain.
#'
#'
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
      if (any(is.nan(w))) {
        break
      }
      if (any((d >= 0) == d_geq_old)) {
        break
      }

   # }
  }
  return(z)
}

