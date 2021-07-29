#' Remove baseline of a dataset using psalsa
#'
#' @param dir_in          Input directory. Where input data files are
#'   loaded from.
#' @param dir_out         Output directory. Where data files after baseline
#'   correction are stored.
#' @param samples         A vector. Set of samples to which remove the baseline
#'   (e.g.: c(1, 2, 3)).
#' @param time             Sets the dimension to be corrected: drift time or
#'   retention time. Introduce "Retention" for retention time; or "Drift" for
#'   drift time.
#' @param lambda          Smoothing parameter (generally 1e3 - 1e9).
#' @param p               Asymmetry parameter (generally 1e-3 - 1e-6).
#' @param k               Peak height parameter (usually 5\% of maximum intensity).
#' @details \code{gcims_baseline_removal} performs baseline correction on sample
#'   chromatograms (in retention time) or spectra (in drift time).This method
#'   relies in is Psalsa algorithm. Psalsa is based on the Asymmetric Least
#'   Squares (ALS) algorithm for baseline estimation, although slightly modified
#'   to reject the effects of large peaks above the baseline more easily. Psalsa
#'   algorithm depends on three parameters: a penalty on the second derivative
#'   (\code{lambda}), a penalty on the value of the baseline with respect the
#'   value of the signal (\code{p}), and third parameter (\code{k}), that
#'   modifies substantially the value of (\code{p}) in the event of large
#'   intensity peaks.
#' @family Baseline removal functions
#' @return A set of S3 objects.
#' @export
#' @importFrom ptw whit2
#' @references {
#'
#'   Eilers, P.H.C. (2003) A perfect smoother. Anal. Chem.  p. 3631–6.
#'   https://doi.org/10.1021/ac034173t.
#'
#'   Oller-Moreno, S., Pardo, A., Jiménez-Soto, J. M., Samitier, J., & Marco, S.
#'   (2014, February). Adaptive Asymmetric Least Squares baseline estimation for
#'   analytical instruments. In 2014 IEEE 11th International Multi-Conference on
#'   Systems, Signals & Devices (SSD14) (pp. 1-5). IEEE.
#'   https://doi.org/10.1109/SSD.2014.6808837.
#'
#'   }
#'
#' @author Sergio Oller-Moreno,  \email{soller@@.ibecbarcelona.eu}
#'
#'
#' Institute for Bioengieering of Catalonia. Signal and Information Processing for Sensing Systems group. Barcelona, Spain.
#'
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3, 7)
#'
#' # Example of Drift time Baseline Removal:
#' # Before:
#' gcims_plot_spec(dir_in, samples, rt_value = 60, dt_range = NULL, colorby = "Class")
#'
#' # After:
#' time <- "Drift"
#' lambda <- 1E3
#' p <- 0.001
#' k <- -1
#' gcims_baseline_removal(dir_in, dir_out, samples, time, lambda, p, k)
#' gcims_plot_spec(dir_out, samples, rt_value = 60, dt_range = NULL, colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'
gcims_baseline_removal <- function(dir_in, dir_out, samples,
                                   time,lambda = 1E7, p = 0.001, k = -1){


  #-------------#
  #     MAIN    #
  #-------------#


  print(" ")
  print("  /////////////////////////")
  print(" /    Baseline removal   /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  m = -1;
  for (i in c(0, samples)){
    m = m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)

    if (time == "Drift"){
      aux <- t(aux)
      lambda <- 1E7
      p <- 0.001
      k <- -1
    } else if (time == "Retention"){
      lambda <- 1E4
      p <- 0.1
      k <- -1
    }

    psalsa_results <- psalsa(aux, lambda, p, k)
    aux  <- psalsa_results$corrected

    if (time == "Drift"){
      aux <- t(aux)
    } else if (time == "Retention"){
    }

    aux_list$data$data_df <- round(aux)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }
}






#---------------#
#   FUNCTIONS   #
#---------------#


#---------------#
#   psalsa_one  #
#---------------#

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
    if (all((d >= 0) == d_geq_old)) {#any
      break
    }

    # }
  }
  return(z)
}


#------------#
#   psalsa   #
#------------#

psalsa <- function(data, lambda = 1E7, p = 0.001, k = -1, maxit = 25) {
  if (is.matrix(data)) {
    estbaseline <- 0*data
    for (i in 1:nrow(data)) {
      estbaseline[i, ] <- psalsa_one(data[i,], lambda, p, k, maxit)
    }
  } else {
    estbaseline <- psalsa_one(data, lambda, p, k, maxit)
  }
  return (list(baseline = estbaseline, corrected = data - estbaseline))

}






