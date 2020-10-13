#' Read GC-IMS Samples from .csv files


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param skip           Number of lines to skip before reading gcims data.
#' @return An R object that contains the matrix and the retention times
#' @family Reading functions
#' @export
#' @importFrom readr read_csv cols
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_read_samples <- function(dir_in, dir_out, skip = 0) {
  setwd(dir_in)
  print(" ")
  print("  /////////////////////////")
  print(" /    Reading Samples    /")
  print("/////////////////////////")
  print(" ")

  metadata <- list(condition = NULL, experiment = NULL, group = NULL,
                   index = NULL, raw_path = NULL, replicate = NULL)
  data <- list(retention_time = NULL, drift_time = NULL, data_df = NULL)

  files <- list.files(pattern = ".csv")
  if (length(files) == 0){
    stop("This folder does not contains any .csv file")
  }

  m <- 0
  for (i in 1:length(files)){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(files)))
    aux_string <- paste0("M", i, ".rds")
    single_file <- read_csv(files[i], skip = skip,
                     progress = FALSE, col_names = FALSE,
                     col_types = cols(.default = "d"))
    dd <- single_file[-1, -c(1:2)]
    data$data_df = dd
    dd_list <- list(metadata = metadata, data = data)
    setwd(dir_out)
    saveRDS(dd_list, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


#' Read GC-IMS Samples from matlab

#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @return An R object that contains the matrix and the retention times
#' @family Reading functions
#' @export
#' @importFrom R.matlab readMat
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_read_mat <- function(dir_in, dir_out) {
  setwd(dir_in)

  print(" ")
  print("  /////////////////////////")
  print(" /    Reading Samples    /")
  print("/////////////////////////")
  print(" ")

  metadata <- NULL
  files <- list.files(pattern = ".mat")
  if (length(files) == 0){
    stop("This folder does not contains any .mat file")
  }
  m <- 0
  for (i in 1:length(files)){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(files)))
    aux_string <- paste0("M", i, ".rds")
    dd <- readMat(files[i])[[1]]
    dd <- list(metadata = metadata, data = dd)
    setwd(dir_out)
    saveRDS(dd, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}
