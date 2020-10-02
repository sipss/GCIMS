
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
  if ((!is.null(ref_number)) & (sum(samples == ref_number) == 0)){
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
      aux <- readRDS(aux_string)
      acc_sample <- matrix(0, dim(aux)[1], dim(aux)[2])
      rm(aux_string, aux)

      m <- 0
      for (i in c(samples)){
        m <- m + 1
        aux_string <- paste0("M", i, ".rds")
        aux <- readRDS(aux_string)
        acc_sample <- acc_sample + aux
      }

      M <- acc_sample / length(samples)

    setwd(dir_out)
    saveRDS(M, file = paste0("M", 0, ".rds"))
    setwd(dir_in)
  }
}


