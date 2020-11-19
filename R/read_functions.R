#' Read GC-IMS Samples from .csv files


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param skip           Number of lines to skip before reading gcims data.
#' @return An R object that contains the matrix and the retention times
#' @family Reading functions
#' @export
#' @importFrom readr read_csv cols type_convert
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

  sftwr <- menu(c("Yes", "No"), title="Were your samples converter to .csv using the VOCal software?")

  if (sftwr == 1){
    files <- list.files(pattern = ".csv")
    if (length(files) == 0){
      stop("This folder does not contains any .csv file")
    }

    m <- 0
    for (i in 1:length(files)){
      m <- m + 1
      print(paste0("Sample ", m, " of ", length(files)))
      aux_string <- paste0("M", i, ".rds")

      # METADATA
      #metadata$experiment <- files[i]
      #metadata$raw_path
      # DATA
      single_file <- read_csv(files[i], skip = skip,
                              progress = FALSE, col_names = FALSE,
                              col_types = cols(.default = "c"))

      dd <- t(single_file[-c(1:2), -c(1:2)])
      # Depende del tipo de instrumento OJO
      data$data_df = - type_convert(dd, col_types = cols(.default = "d"))                    # Signo Menos solo para gasdormund programa nuevo
      data$retention_time <- as.numeric(single_file[1, -c(1:2)])
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <-type_convert(single_file[-c(1:2), 2] , col_types = cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
      data$drift_time <-  data$drift_time[[1]]                                                # Depende del tipo de instrumento OJO

      # Join
      dd_list <- list(metadata = metadata, data = data)
      setwd(dir_out)
      saveRDS(dd_list, file = paste0("M", i, ".rds"))
      setwd(dir_in)
    }
  }
  else {
    files <- list.files(pattern = ".csv")
    if (length(files) == 0){
      stop("This folder does not contains any .csv file")
    }

    m <- 0
    for (i in 1:length(files)){
      m <- m + 1
      print(paste0("Sample ", m, " of ", length(files)))
      aux_string <- paste0("M", i, ".rds")

      # METADATA
      #metadata$experiment <- files[i]
      #metadata$raw_path
      # DATA
      single_file <- read_csv(files[i], skip = skip,
                              progress = FALSE, col_names = FALSE,
                              col_types = cols(.default = "c"))

      dd <- single_file[-c(1:2), -c(1:2)]
      # Depende del tipo de instrumento OJO
      data$data_df = - type_convert(dd, col_types = cols(.default = "d"))                    # Signo Menos solo para gasdormund programa nuevo
      data$retention_time <- as.numeric(single_file[1, -c(1:2)])
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <-type_convert(single_file[-c(1:2), 2] , col_types = cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
      data$drift_time <-  data$drift_time[[1]]                                                # Depende del tipo de instrumento OJO

      # Join
      dd_list <- list(metadata = metadata, data = data)
      setwd(dir_out)
      saveRDS(dd_list, file = paste0("M", i, ".rds"))
      setwd(dir_in)
    }
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
