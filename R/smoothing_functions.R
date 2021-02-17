#' Data Smoothing using Savitzky-Golay filters


#' @param dir_in           The input directory.
#' @param dir_out          The output directory.
#' @param samples          The set of samples to be processed.
#' @param time             It indicates if the correction is going to be in the
#'                         drift time or in the retention time. It should be
#'                         introduce "Retention" for correcting the retention
#'                         time; or "Drift" for the drift time.
#' @param filter_length    Numerical. Length of the filter.
#' @param polynomial_order Numerical. Order of the polynomial.
#' @return A filtered  gcims dataset.
#' @family Smoothing functions
#' @export
#' @importFrom signal sgolayfilt
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

  setwd(dir_in)
  m = -1
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
    } else if (time == "Retention"){
    }

    n <- dim(aux)[1]

    for (j in (1:n)){
      aux[j, ] <- sgolayfilt(aux[j, ], p = polynomial_order, n = filter_length)
    }

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


