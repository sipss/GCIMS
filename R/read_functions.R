#' Read GC-IMS Samples from .csv files


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param file           Name of the metadata file.
#' @param sftwr          Use 1 if the software that has converted the original files into
#' csv format is the new one. Use 2 if the software is the old one.
#' @return An R object that contains the matrix and the retention times
#' @family Reading functions
#' @export
#' @importFrom readr read_csv cols type_convert
#' @importFrom utils menu
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_read_samples <- function(dir_in, dir_out, file, sftwr) {
  setwd(dir_in)
  print(" ")
  print("  /////////////////////////")
  print(" /    Reading Samples    /")
  print("/////////////////////////")
  print(" ")

  metadata <- list(NULL)
  data <- list(retention_time = NULL, drift_time = NULL, data_df = NULL)

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
      metadata <- i
      # DATA
      single_file <- read_csv(files[i],
                              progress = FALSE, skip = 1, col_names = FALSE,
                              col_types = cols(.default = "c"))

      dd <- as.data.frame(t(single_file[-1, -1]))
      # Depende del tipo de instrumento OJO
      data$data_df <- type_convert(dd, col_types = cols(.default = "d"))                    # Signo Menos solo para gasdormund programa viejo
      data$retention_time <- type_convert(single_file[-1, 1] , col_types = cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
      data$retention_time <-  data$retention_time[[1]]
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <- as.numeric(single_file[1, -(dim(single_file)[2])])#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.

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
      metadata <- i
      # DATA
      single_file <- read_csv(files[i], skip = 130,
                              progress = FALSE, col_names = FALSE,
                              col_types = cols(.default = "c"))

      dd <- single_file[-c(1:2), -c(1:2)]
      # Depende del tipo de instrumento OJO
      data$data_df = - type_convert(dd, col_types = cols(.default = "d"))                    # Signo Menos solo para gasdormund programa nuevo
      data$retention_time <- as.numeric(single_file[1, -c(1:2)])
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <- type_convert(single_file[-c(1:2), 2] , col_types = cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
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


#' Read metadata form a csv

#' @param dir_in          The input directory.
#' @param samples         The set of samples to be processed.
#' @param file            Name of the file that contains the metadata. It must be an excel file
#' @return An R object that contains the metadata information
#' @family Reading functions
#' @export
#' @importFrom readxl read_excel
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_read_metadata <- function(dir_in, samples, file) {
  setwd(dir_in)

  print(" ")
  print("  //////////////////////////")
  print(" /    Reading Metadata    /")
  print("//////////////////////////")
  print(" ")

  Metadatafile <- read_excel(file)
  metadata <- NULL
  m <- 0
  for (i in 1:length(samples)){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string)
    metadata <- Metadatafile[which(Metadatafile$Name == aux_list$metadata),]
    aux_list$metadata <- metadata
    saveRDS(aux_list, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }
}

