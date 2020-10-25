
#' Select or create a reference sample


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param ref_number      Reference sample. If NULL the reference
#'                        is the mean sample.
#' @return A reference sample
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

gcims_create_reference <- function(dir_in, dir_out, samples, ref_number){


  print(" ")
  print("  ////////////////////////////////////")
  print(" /     Create a reference sample    /")
  print("////////////////////////////////////")
  print(" ")


  setwd(dir_in)
  if (length(samples) == 0){
    stop("No vector with sample numbers was provided")
  } else if(!is.numeric(samples)){
    stop("samples must be a vector of integer values")
  } else if(!all(samples == round(samples))){
    stop("samples must be a vector of integer values")
  } else if(length(unique(samples)) != length(samples)){
    stop("the vector of sample numbers contains repeated samples")
  } else if(any(samples == 0)){
    stop("samples contains sample 0. A reference sample can't be used to create another reference sample")
  } else if(any(samples < 0)){
      stop("no sample in samples can have a negative index")
  }


  if (length(ref_number) > 1){
    stop("ref_number must be a scalar or NULL")
  }else if(!is.numeric(ref_number) & !is.null(ref_number)){
    stop("ref_number must be a scalar or NULL")
  }else if ((!is.null(ref_number)) & (sum(samples == ref_number) == 0)){
    stop("ref_number is not containend in samples")
  }


  if (sum(samples == ref_number) == 1){
    print(paste0("Reference Sample is: M", ref_number))
    aux_string <- paste0("M", ref_number, ".rds")
    M <- readRDS(aux_string)
    setwd(dir_out)
    saveRDS(M, file = paste0("M", 0, ".rds"))
    setwd(dir_in)
  }

  if (is.null(ref_number)){
      print(paste0("Reference Sample is the average sample"))
      aux_string <- paste0("M", samples[1], ".rds")

      aux_list <- readRDS(aux_string) #new
      aux <- t(as.matrix(aux_list$data$data_df)) #new

      #aux <- readRDS(aux_string)
      acc_sample <- matrix(0, dim(aux)[1], dim(aux)[2])
      rm(aux_string, aux)

      m <- 0
      for (i in samples){
        m <- m + 1
        aux_string <- paste0("M", i, ".rds")
        aux_list <- readRDS(aux_string) #new
        aux <- t(as.matrix(aux_list$data$data_df)) #new
        acc_sample <- acc_sample + aux
      }

      aux_list$data$data_df <- t(acc_sample / length(samples))
      M <- aux_list

    setwd(dir_out)
    saveRDS(M, file = paste0("M", 0, ".rds"))
    setwd(dir_in)
  }
}


