

#' Removes the Reactant Ion Peak (RIP) from samples

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where RIP removed data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their RIP has to be removed.
#' @details `gcims_remove_rip` substitutes the RIP by its corresponding
#'   linear approximation to the RIP baseline, for every spectrum in a sample.
#'   This process is repeated for all samples in `samples`. Use this
#'   function if you are interested in enhancing the contrast of peaks of sample
#'   images / chromatograms / spectra to be obtained from
#'   `gcims_view_sample` / `gcims_plot_chrom` /
#'   `gcims_plot_spec`.
#' @return A Set of S3 objects.
#' @family Utility functions
#' @export
#' @importFrom pracma findpeaks
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Reactant Ion Peak removal
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' gcims_remove_rip(dir_in, dir_out, samples)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#'
gcims_remove_rip <- function(dir_in, dir_out, samples){

  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  for (i in samples){
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- (as.matrix(aux_list$data$data_df))

    # Compute the total ion spectra
    aux_2 <- rowSums(aux)
    peaks_info <- findpeaks(aux_2)

    # Look for rip position
    rip_pos_ind <- which.max(peaks_info[ , 1])
    rip_pos <- peaks_info[rip_pos_ind, 2]

    # Look for the beginning and ending of the RIP (searching the closest minima to it)
    valleys_info <- findpeaks(-aux_2)
    valleys_pos <- valleys_info[ , 2]
    closest_valley_ind <- which.min(abs(valleys_pos - rip_pos))


    # Select the RIP region
    if(valleys_pos[closest_valley_ind] < rip_pos){
      rip_bounds <- valleys_pos[c(closest_valley_ind,closest_valley_ind + 1)]
    } else if (valleys_pos[closest_valley_ind] > rip_pos){
      if(closest_valley_ind > 1){
        rip_bounds <- valleys_pos[c(closest_valley_ind - 1,closest_valley_ind)]
      } else {
        rip_bounds <- valleys_pos[c(1, closest_valley_ind)]
      }

    }


    # Erase the peak (with style...)
    for (j in (1:dim(aux)[2])){
      aux[rip_bounds[1]: rip_bounds[2], j] <- seq(from = aux[rip_bounds[1], j],
                                                  to = aux[rip_bounds[2], j],
                                                  length.out = length(rip_bounds[1]: rip_bounds[2]))
    }

    aux_list$data$data_df <- round(aux)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}

#' Reshape samples
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where reshaped samples are stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to be reshaped.
#' @description Force that all samples have the same minimal data dimensions.
#' @details It may happen that not all samples in a dataset have the same size (dimensions).
#'   When this happens, an adjustment of the dimensions of samples is nedeed.`gcims_reshape_samples()`
#'   ensures that all samples in a dataset have the same dimensions (number of data points)
#'   in retention and drift times. `gcims_reshape_samples()` checks what are the minimum retention /
#'   drift time ranges a cuts all samples according to these ranges.
#'
#' @return A set of S3 objects.
#' @family Utility functions
#' @export
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3, 7)
#'
#' # Example of reshaping samples
#' # (all samples must have the same
#' # retention and drift time dimensions).
#'
#' # In this particular example, data dimensions is samples
#' # 3 and 7 is the same, so the expected behavior of
#' # the function gcims_reshaphe_samples() is keeping the
#' # same dimensionality.
#'
#' # Before:
#' nrow_before <- matrix(0, nrow = length(samples), ncol = 1)
#' ncol_before <- matrix(0, nrow = length(samples), ncol = 1)
#' for (i in seq_along(samples)){
#'   aux_string <- paste0("M", samples[i], ".rds")
#'   aux_list <- readRDS(file.path(dir_in, aux_string))
#'   nrow_before[i, 1] <- nrow(as.matrix(aux_list$data$data_df))
#'   ncol_before[i, 1] <- ncol(as.matrix(aux_list$data$data_df))
#' }
#'
#' # After:
#' gcims_reshape_samples(dir_in, dir_out, samples)
#' nrow_after <- matrix(0, nrow = length(samples), ncol = 1)
#' ncol_after <- matrix(0, nrow = length(samples), ncol = 1)
#' for (i in seq_along(samples)){
#'   aux_string <- paste0("M", samples[i], ".rds")
#'   aux_list <- readRDS(file.path(dir_out, aux_string))
#'   nrow_after[i, 1] <- nrow(as.matrix(aux_list$data$data_df))
#'   ncol_after[i, 1] <- ncol(as.matrix(aux_list$data$data_df))
#' }
#'
#' reshaping_info <- cbind(nrow_before, ncol_before, nrow_after, ncol_after)
#' colnames(reshaping_info) <- c("rows_before", "columns_before", "rows_after", "columns_after")
#' rownames(reshaping_info) <- c("M3", "M7")
#'
#' print(reshaping_info)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
gcims_reshape_samples <- function(dir_in, dir_out, samples) {

  dimensions <- list(NULL)
  for (i in seq_along(samples)) {
    aux_string <- paste0("M", samples[i], ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)
    dimensions[[i]] <- dim(aux)
  }

  rts <- NULL
  dts <- NULL
  for (i in seq_along(samples)) {
    dts <- c(dts, dimensions[[i]][1])
    rts <- c(rts, dimensions[[i]][2])
  }

  rts <- min(rts)
  dts <- min(dts)

  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }
  for (i in seq_along(samples)) {
    aux_string <- paste0("M", samples[i], ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)
    aux <- aux[1:dts, 1:rts]
    aux_list$data$data_df <- aux
    aux_list$data$retention_time <- aux_list$data$retention_time[1:rts]
    aux_list$data$drift_time <- aux_list$data$drift_time[1:dts]
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", samples[i], ".rds")))
  }
}


#' Cut samples


#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where cut data files are
#'   stored.
#' @param samples         A vector. Set of samples to be cut (e.g.: c(1, 2,
#'   3)).
#' @param rt_range        A  vector of two components. Beginning and end of the
#'   retention time cut.If NULL the complete retention time range is used.
#' @param dt_range        A  vector of two components. Beginning and end of the
#'   drift time cut. If NULL the complete drift time range is used.
#' @return A set of S3 objects.
#' @description Cut samples in a retention time - drift time rectangle.
#' @details `gcims_cut_samples()` cuts a sample in a retention time - drift time
#'   rectangle according to the retention time / drift time ranges given by
#'   function arguments `rt_range` / `dt_range`. Use this function to
#'   focus on the retention time - drift time region where chemical information
#'   is more abundant, that is, where you can find a high peak densities by
#'   visual inspection.
#' @note By reducing the size of data, the computational time of signal
#'   pre-processing stage can reduced substantially.
#' @family Utility functions
#' @export
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example Sample data cutting:
#' # Before:
#' gcims_view_sample(
#'   dir_in, sample_num = samples,
#'   rt_range = NULL, dt_range = NULL, transform = FALSE
#' )
#'
#' # After:
#' rt_range <-c(70, 125)
#' dt_range <- c(8, 9.25)
#' gcims_cut_samples(dir_in, dir_out, samples, rt_range, dt_range)
#' gcims_view_sample(
#'   dir_out, sample_num = samples,
#'   rt_range = NULL, dt_range = NULL, transform = FALSE
#' )
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#
gcims_cut_samples <- function(dir_in, dir_out, samples, rt_range, dt_range){

  #---------------#
  #   FUNCTIONS   #
  #---------------#


  #-----------------#
  #   cut_samples   #
  #-----------------#


  cut_samples <- function(loaded_sample, rt_range = NULL, dt_range = NULL) {

    aux_list <- loaded_sample
    aux <- (as.matrix(aux_list$data$data_df))

    #SOME CHECKS
    retention_time <- aux_list$data$retention_time
    drift_time <- aux_list$data$drift_time
    cond_1_rt <- (rt_range[1] - retention_time[1]) < 0
    cond_2_rt <- (rt_range[2] - retention_time[length(retention_time)]) > 0
    cond_1_dt <-(dt_range[1] - drift_time[1]) < 0
    cond_2_dt <-(dt_range[2] - drift_time[length(drift_time)]) > 0


    if(is.null(rt_range)){# old
      rt_ind <- c(1, dim(aux)[2])

    } else{
      if(cond_1_rt | cond_2_rt){
        stop("Retention time range out of bounds.")
      }
      rt_ind  <- c(which.min(abs(retention_time - rt_range[1])), which.min(abs(retention_time - rt_range[2])))
      if( rt_ind[1] == rt_ind[2]){
        stop("Initial and Final retention time values can't be equal in the variable rt_range.")
      }
    }



    if(is.null(dt_range)){# old
      dt_ind <- c(1, dim(aux)[1])
    } else{
      if(cond_1_dt | cond_2_dt){
        stop("Drift time range out of bounds.")
      }
      dt_ind  <- c(which.min(abs(drift_time - dt_range[1])), which.min(abs(drift_time - dt_range[2])))
      if( dt_ind[1] == dt_ind[2]){
        stop("Initial and Final drift time values can't be equal in the variable dt_range.")
      }
    }

    sel_index_rt <- rt_ind[1]: rt_ind[2]
    sel_index_dt <- dt_ind[1]: dt_ind[2]

    if(is.null(rt_range)){

    } else if (methods::is(sel_index_rt, "integer") & (sel_index_rt[2] > sel_index_rt[1])){
    } else {
      stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
    }

    if(is.null(dt_range)){

    } else if (methods::is(sel_index_dt, "integer") & (sel_index_dt[2] > sel_index_dt[1])){
    } else {
      stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
    }

    aux_list$data$retention_time <- retention_time[sel_index_rt]
    aux_list$data$drift_time  <- drift_time[sel_index_dt]
    aux_list$data$data_df <- round(aux_list$data$data_df[sel_index_dt, sel_index_rt]) #transposed when matrix!

    return(aux_list)
  }


  #-------------#
  #     MAIN    #
  #-------------#

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(samples)) {
    aux_string <- paste0("M", samples[i], ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    M <- cut_samples(aux_list, rt_range, dt_range)
    saveRDS(M, file = file.path(dir_out, paste0("M", samples[i], ".rds")))
  }

}




#' Shift samples in retention time.
#'
#' Perform a shift in retention time so to the injection point of all samples
#' coincide.
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where decimated data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to be shifted in retention time.
#' @return A set of S3 objets.
#' @family Utility functions
#' @export
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#'
#' # Before:
#' samples <- c(3, 7)
#' gcims_plot_chrom(dir_in, samples, dt_value = NULL,  rt_range = NULL, colorby = "Class")
#'
#' # After:
#' gcims_shift_rt(dir_in, dir_out, samples)
#' gcims_plot_chrom(dir_out, samples, dt_value = NULL,  rt_range = NULL, colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#'

gcims_shift_rt <- function(dir_in, dir_out, samples){

  tics <- NULL
  for (i in samples){
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string)) #new
    aux <- as.matrix(aux_list$data$data_df)
    tic <- colSums(aux)
    tics <- rbind(tics, tic)
  }

  mins <- apply(tics, 1, which.min)
  # plot(x = c(1:64), y = mins)
  # intensities <- apply(tics, 1, min)
  # plot(x = c(1:64), y = intensities)
  referenceindex <- min(mins)
  referencetime <- samples[which.max(mins)]

  for (i in samples){
    if (i == referencetime) {
      file.copy(
        from = file.path(dir_in, paste0("M", referencetime, ".rds")),
        to = file.path(dir_out, paste0("M", referencetime, ".rds"))
      )
    } else {
      reference <- readRDS(file.path(dir_in, paste0("M", referencetime, ".rds")))
      reference <- reference$data$data_df
      referencetic <- colSums(reference)
      reference <- reference[,-c(1:abs(referenceindex - which.min(referencetic)))]
      dimension <- dim(reference)[2]
      aux_string <- paste0("M", i, ".rds")
      aux_list <- readRDS(file.path(dir_in, aux_string)) #new
      aux <- as.matrix(aux_list$data$data_df)
      auxtic <- colSums(aux)
      if (abs(referenceindex - which.min(auxtic)) > 0){
        aux <- aux[,-c(1:abs(referenceindex - which.min(auxtic)))]
      } else {
        aux <- aux
      }
      aux_list$data$retention_time <- aux_list$data$retention_time[1:dimension]
      aux_list$data$data_df <- aux[,c(1:dimension)]
      M <- aux_list
      saveRDS(M, file = file.path(dir_out, paste0("M", i, ".rds")))
      }
  }
}


new_progress_bar <- function(...) {
  if (!requireNamespace("progress", quietly = TRUE)) {
    rlang::inform(
      message = c("i" = 'Use install.packages("progress") to get a progress bar'),
      class = "GCIMS_suggest_install_progress",
      .frequency = "once",
      .frequeny_id = "suggest_progress_installation"
    )
    # This code path returns a dummy object exposing a compatible API of
    # progress:progress_bar$new(...).
    dummy_pbar <- list(
      tick = function() {}
    )
    return(dummy_pbar)
  }
  progress::progress_bar$new(...)
}

require_pkgs <- function(pkg, msgs = NULL, ...) {
  have_pkgs <- purrr::map_lgl(pkg, function(p) {requireNamespace(p, quietly = TRUE)})
  names(have_pkgs) <- pkg
  if (!all(have_pkgs)) {
    missing_pkgs <- names(have_pkgs)[!have_pkgs]
    aval_pkgs <- rownames(utils::available.packages())
    missing_cran_pkgs <- intersect(missing_pkgs, aval_pkgs)
    missing_bioc_pkgs <- setdiff(missing_pkgs, missing_cran_pkgs)
    if (length(missing_bioc_pkgs) > 0) {
      if (!"BiocManager" %in% rownames(utils::installed.packages())) {
        missing_cran_pkgs <- c(missing_cran_pkgs, "BiocManager")
      }
    }
    if (length(missing_cran_pkgs) > 0) {
      missing_cran_pkgs <- deparse(missing_cran_pkgs)
    }
    if (length(missing_bioc_pkgs) > 0) {
      missing_bioc_pkgs <- deparse(missing_bioc_pkgs)
    }
    parent_call <- format(rlang::caller_call())
    rlang::abort(
      message = c(
        glue::glue("{parent_call} requires additional packages. Please install them. You may want to use:", parent_call = parent_call),
        glue::glue("    install.packages({missing_cran_pkgs}) and", missing_cran_pkgs = missing_cran_pkgs),
        glue::glue("    BiocManager::install({missing_bioc_pkgs})", missing_bioc_pkgs = missing_bioc_pkgs),
        msgs
      ),
      ...
    )
  }
}


units_to_points <- function(length_phys, step_phys, must_odd = FALSE) {
  length_pts <- round(length_phys/step_phys)
  if (must_odd) {
    length_pts <- length_pts + (length_pts %% 2 == 0) # the filter length in points
  }
  length_pts
}
