#' Data Smoothing using Savitzky-Golay filters


#' @param dir_in           Input directory. Where input data files are loaded
#'   from.
#' @param dir_out          Output directory. Where smoothed data files are
#'   stored.
#' @param samples          A vector. Set of samples to be filtered (e.g.: c(1,
#'   2, 3)).
#' @param time             Sets the dimension to be corrected: drift time or
#'   retention time. Introduce "Retention" for retention time; or "Drift" for
#'   drift time.
#' @param filter_length    Numerical. Length of the filter.
#' @param polynomial_order Numerical. Order of the polynomial.
#'
#' @details `gcims_smoothing` performs digital smoothing on sample
#'   chromatograms (in retention time) or spectra (in drift time). This
#'   smoothing process is required to reduce the noise level of data samples.
#'   `gcims_smoothing` relies in Savitzky-Golay filters to reduce noise
#'   contribution in data. In a nutshell, Savitzky-Golay filters reconstruct
#'   data by fitting successive sets of data points of given length to
#'   low-degree polynomials. Short *filter_length* values don't reduce
#'   noise significantly, while large ones may distort signal shape. Similarly,
#'   polynomials with low *polynomial_order* may not be flexible enough to
#'   properly reconstruct the signal, while polynomials with high
#'   *polynomial_order* can be too flexible and model noise on the
#'   reconstructed signal. A trade-off for these parameter values should be
#'   found.
#' @return A set of S3 objects.
#' @family Smoothing functions
#' @export
#' @references { Savitzky, A. and Golay, M.J.E. (1964) Smoothing and
#'  Differentiation of Data by Simplified Least Squares Procedures. Analytical
#'  Chemistry, 36, 1627–39. https://doi.org/10.1021/ac60214a047.
#'
#'  Schafer, R.W. (2011) What is a savitzky-golay filter? IEEE Signal Processing
#'  Magazine, Institute of Electrical and Electronics Engineers Inc. 28, 111–7.
#'  https://doi.org/10.1109/MSP.2011.941097. }
#'
#' @importFrom signal sgolayfilt
#'
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Retention time smoothing:
#' # Before:
#' gcims_plot_chrom(dir_in, samples, dt_value = 8.5,  rt_range = NULL, colorby = "Class")
#'
#' # After:
#' gcims_smoothing(dir_in, dir_out, samples, time = "Retention")
#' gcims_plot_chrom(dir_out, samples, dt_value = 8.5,  rt_range = NULL, colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'
gcims_smoothing <- function (dir_in, dir_out, samples, time,
                             filter_length = 19, polynomial_order = 2){

  print(" ")
  print("  /////////////////////////")
  print(" /   Filtering the data /")
  print("/////////////////////////")
  print(" ")

  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  m = -1
  for (i in c(0, samples)){
    m = m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }

    aux_list <- readRDS(file.path(dir_in, paste0("M", i, ".rds")))
    aux <- as.matrix(aux_list$data$data_df)

    if (time == "Drift"){
      aux <- t(aux)
    }

    for (j in seq_len(nrow(aux))) {
      aux[j, ] <- signal::sgolayfilt(aux[j, ], p = polynomial_order, n = filter_length)
    }

    if (time == "Drift"){
      aux <- t(aux)
    }

    aux_list$data$data_df <- round(aux)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}
