
#' Cut samples in a retention time - drift time rectangle


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param rt_range        A  vector of two components. Beginning and end
#'                        of the retention time cut.
#' @param dt_range        A  vector of two components. Beginning and end
#'                        of the drift time cut.
#' @return An cut gcims dataset.
#' @family Cutting functions
#' @export
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example Sample data cutting:
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' rt_range <-c(70, 125)
#' dt_range <- c(8, 9.25)
#' gcims_cut_samples(dir_in, dir_out, samples, rt_range, dt_range)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#
gcims_cut_samples <- function(dir_in, dir_out, samples, rt_range, dt_range){


  print(" ")
  print("  /////////////////////")
  print(" /    Cut Samples    /")
  print("/////////////////////")
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
    M <- cut_samples(aux_list, rt_range, dt_range)
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }

}



#' Cuts retention time and drift time
#'
#' @param loaded_sample   Current sample to be cut.
#' @param rt_range        A  vector of two components. Beginning and end
#'                        of the retention time cut.
#' @param dt_range        A  vector of two components. Beginning and end
#'                        of the drift time cut.
#' @return An cut matrix.
#' @family Cutting functions
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
#'
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

  } else if((class(sel_index_rt) == "integer") & (sel_index_rt[2] > sel_index_rt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
  }

  if(is.null(dt_range)){

  } else if((class(sel_index_dt) == "integer") & (sel_index_dt[2] > sel_index_dt[1])){
  } else {
    stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
  }

  aux_list$data$retention_time <- retention_time[sel_index_rt]
  aux_list$data$drift_time  <- drift_time[sel_index_dt]
  aux_list$data$data_df <- round(aux_list$data$data_df[sel_index_dt, sel_index_rt]) #transposed when matrix!

  return(aux_list)
}
