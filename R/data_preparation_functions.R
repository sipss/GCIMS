
#' Data preparation


#' @param dir_in           Input directory. Where input data files are loaded
#'   from.
#' @param dir_out          Output directory. Where prepared data files are
#'   stored.
#' @param samples          A vector. Set of samples to which remove the baseline
#'   (e.g.: c(1, 2, 3)).
#' @param params           A list of lists with parameter values to perform
#'                         filtering and decimation in retention and drift time axes.
#' @description Interpolates, filters and decimates data.
#' @details  `gcims_prepare_data()` prepares samples for the subsequent data
#'  pre-processing stage. Each of the samples is interpolated, filtered (optional),
#'  and decimated (optional). To select which pre-processing steps apply to data use the input
#'  variable `params`, described below:
#' \describe{
#'   \item{`filter$do, filter$order_rt, filter$length_rt, filter$order_dt`, `filter$length_dt`}{
#'   Set of variables that control data filtering using Savitzky-Golay filters.  `do` is a Boolean variable.
#'   If TRUE, digital filters are applied. `order_rt` and `length_rt` are positive integers that control
#'   the polynomial order and the filter length in retention time, respectively. The analogous control
#'   variables for drift time filtering are `order_dt` and `length_dt`.}
#'   \item{`decimate$do, decimate$dfactor_rt, decimate$factor_dt`}{
#'   Set of variables that control data decimation. `do` is a Boolean variable is a Boolean variable.
#'   If TRUE, decimation is performed on data. `factor_rt` and `factor_dt` are positive integers that
#'   control the decimation factors in retention and drift times, respectively.}
#'   }
#' @note Note that filter length must be an odd number bigger than the
#'  polynomial order of the filter.
#' @return A set of S3 objects. Additionally, it returns a list with two variables needed to perform sample alignment in drift
#'  and retention time axes: `tis` and `rics`. First variable (a matrix) corresponds to the Total Ion Spectra of
#'  samples, while second one (a matrix) to their Reactive Ion Chromatograms, respectively. These outputs are input variables for
#'  the function `gcims_align_data()`.
#' @family Data preparation functions
#' @export
#' @references {García, S., Luengo, J. and Herrera, F., 2015. Data preprocessing
#'  in data mining (Vol. 72, pp. 59-139). Cham, Switzerland: Springer International
#'  Publishing.
#'  \doi{10.1007/978-3-319-10247-4}
#'   }
#'
#' @importFrom signal sgolayfilt
#' @importFrom signal interp1
#'
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3,7)
#'
#' # Example of digital smoothing and decimation, in both axes.
#' # To show the effect data preparation two EIC are compared:
#'
#' # Before:
#' gcims_plot_chrom(dir_in, samples, dt_value = 8.5,  rt_range = NULL, colorby = "Class")
#'
#' # Data preparation:
#' alignment_data <- gcims_prepare_data(
#'   dir_in = dir_in,
#'   dir_out = dir_out,
#'   samples = samples,
#'   params = list(
#'     filter = list(
#'       do = TRUE,
#'       order_rt = 2,
#'       length_rt = 19,
#'       order_dt = 2,
#'       length_dt = 19
#'     ),
#'     decimate = list(
#'       do = TRUE,
#'       factor_rt = 2,
#'       factor_dt = 2
#'     )
#'   )
#' )
#'
#' # After:
#' gcims_plot_chrom(dir_out, samples, dt_value = 8.5,  rt_range = NULL, colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' }
#'
gcims_prepare_data <- function (dir_in, dir_out, samples, params){

  # Function to know if a number is integer
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }

  # SOME GENERAL CHECKS
  do_filter <- params$filter$do
  if(!is.logical(do_filter)){
    stop("the parameter params$filter$do must be logical")
  }

  do_decimate <- params$decimate$do
  if(!is.logical(do_decimate)){
    stop("the parameter params$decimate$do must be logical")
  }

  #CHECK filtering values
  if (do_filter == TRUE){
    order_rt <- params$filter$order_rt
    length_rt <-params$filter$length_rt
    order_dt <-params$filter$order_dt
    length_dt <-params$filter$length_dt

    # Give some reasonable numbers in case the user doesn't provide any´.
    if(is.null(order_rt)){
      order_rt = 2
    }
    if(is.null(length_rt)){
      length_rt = 19
    }
    if(is.null(order_dt)){
      order_dt = 2
    }
    if(is.null(length_dt)){
      length_dt = 19
    }

    # TEST that the introduced numbers are positive integers and that the selection does not produces an error

    if(any((c(order_rt, length_rt, order_dt, length_dt)) < 0)){
      stop("Make sure that filter orders and lengths are positive")
    }
    if(all(is.wholenumber(c(order_rt, length_rt, order_dt, length_dt)))){
    } else {
      stop("Make sure that filter orders and lengths are whole numbers")
    }
    if(order_rt > length_rt){
      stop("The length of the filter must be higher than the order of the polynomial (Retention time).")
    }
    if(order_dt > length_dt){
      stop("The length of the filter must be higher than the order of the polynomial (Drift time).")
    }
  }

  # CHECK decimation values
  if(do_decimate == TRUE){
    factor_rt <- params$decimate$factor_rt
    factor_dt <- params$decimate$factor_dt
  }

  # TEST that the introduced numbers are positive integers and that the selection does not produces an error
  if(any((c(factor_rt,factor_dt)) < 0)){
    stop("Make sure that decimation factors are positive")
  }
  if(all(is.wholenumber(c(factor_rt,factor_dt)))){
  } else {
    stop("ake sure that decimation factors are whole numbers")
  }


  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  m <- 0
  for (i in  samples){
    m <- m + 1
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux_list <- interpolate(aux_list,m)
    if(do_filter == TRUE){
      aux_list <- smoothing(aux_list, m,order_rt, length_rt, order_dt, length_dt)
    }
    if(do_decimate == TRUE){
      aux_list <- decimate_impl(aux_list, m, factor_rt, factor_dt)
    }
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))

    # Create reference information for alignment
    if (m == 1){
      # Compute the dimensions of aux after the first pre-processing stage
      aux <- as.matrix(aux_list$data$data_df)
      # Create the variables rics and tis
      dimensions <- dim(aux)
      rics <- matrix(0, nrow = length(samples) , ncol = dimensions[2])
      tis <- matrix(0, nrow = length(samples) , ncol = dimensions[1])
    }
    rics[m, ] <- compute_ric(aux_list)
    tis[m, ]  <- compute_tis(aux_list)
  }

    alignment_data <- list(rics = rics, tis = tis)
    return(alignment_data)

}

interpolate <- function(aux_list, sample_index){
  # Interpolation function
  interpolate <- function(aux_list, time){
    aux <- as.matrix(aux_list$data$data_df)
    if(time == "Retention"){
      x <- aux_list$data$retention_time
    }else if(time == "Drift"){
      aux <- t(aux)
      x <- aux_list$data$drift_time
    }
    step_x <- (x[length(x)]- x[1]) / (length(x) - 1)
    xi <- seq(from = x[1],
              by = step_x,
              length.out = length(x))

    n <- dim(aux)[1]
    for (j in (1:n)){
      aux[j, ] <- signal::interp1(x, aux[j, ], xi, method = "linear", extrap = TRUE)
    }
    if(time == "Retention"){
      aux_list$data$retention_time <- xi
    } else if (time == "Drift"){
      aux <- t(aux)
      aux_list$data$drift_time <- xi
    }
    aux_list$data$data_df <- aux
    return(aux_list)
  }
  aux_list <- interpolate(aux_list,"Retention")
  aux_list <- interpolate(aux_list,"Drift")
  aux_list$data$data_df <- round(aux_list$data$data_df)
  return(aux_list)
}

smoothing <- function (aux_list, sample_index,
                                  order_rt, length_rt,
                                  order_dt, length_dt){
  # Filtering function
  digital_filter <- function(aux_list, time, order, length){
    aux <- as.matrix(aux_list$data$data_df)
    if (time == "Drift"){
      aux <- t(aux)
    }
    for (j in seq_len(nrow(aux))) {
      aux[j, ] <- signal::sgolayfilt(aux[j, ], p = order, n = length)
    }
    if (time == "Drift"){
      aux <- t(aux)
    }
    aux_list$data$data_df <- aux
    return(aux_list)
  }
  aux_list <-  digital_filter(aux_list,"Retention",order_rt, length_rt)
  aux_list <-  digital_filter(aux_list,"Drift",order_dt, length_dt)
  aux_list$data$data_df <- round(aux_list$data$data_df)
  return(aux_list)
}

decimate_impl <- function(aux_list, sample_index, factor_rt, factor_dt){
  aux <- as.matrix(aux_list$data$data_df)
  rt_index <- seq(from = 1, to = dim(aux)[2], by = factor_rt)
  dt_index <- seq(from = 1, to = dim(aux)[1], by = factor_dt)
  aux_list$data$retention_time <-aux_list$data$retention_time[rt_index]
  aux_list$data$drift_time <-aux_list$data$drift_time[dt_index]
  aux_list$data$data_df<- aux[dt_index,rt_index]
  return(aux_list)

}

compute_ric <- function(aux_list){
  aux <-aux_list$data$data_df
  ric_pos <- which.max(rowSums(aux))
  ric <- aux[ric_pos, ]
  ric <- max(ric) - ric
  ric <- ric/sum(ric)
  return(ric)
}

compute_tis <- function(aux_list){
  aux <- as.matrix(aux_list$data$data_df)
  tis <- rowSums(aux)
}




