#' Read GC-IMS samples from .csv files
#' @param dir_in          Input directory. Where the raw data files are stored.
#' @param dir_out         Output directory. Where data files are stored as R
#'   objects.
#' @param sftwr           Selection of the software version. Use 1 if the
#'   software that has converted the original files into csv format is VOCal.
#'   Use 2 if the software is LAV.
#'
#' @details `gcims_read_samples` stores a set of S3 objects in the
#'   directory `dir_out` (one per sample sample stored in `dir_in`).
#'   Each object is a list containing the following variables: `metadata`
#'   and `data`.In `metadata` you can find information related to the
#'   sample that can be added by the user. On the other hand in `data` you
#'   can find information of the measurement.By default, when the
#'   `metadata` field is created only includes one variable, called
#'   `Name`, although the user can include more fields using the function
#'   [gcims_read_metadata()]. Regarding `data`, it is consist in
#'   three fields:  \describe{
#'      \item{`retention_time`}{ The vector of retention times associated to a sample measurement.}
#'      \item{`drift_time`}{ The vector of drift times associated to a sample
#'   measurement.}
#'      \item{`data_df`}{ A dataframe containing the intensities
#'        of a GCIMS measurement of a sample. Each column in the dataset corresponds
#'        to a particular drift time, while each row to a different retention time.}
#'   }
#'
#' @return A set of S3 objects.
#' @note In the current version of the package, only .csv files converted from
#'   LAV or VOCal software of G.A.S GmbH can be used by the function
#'   `gcims_read_samples`.
#' @family Reading functions.
#' @export
#' @examples
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' gcims_read_samples(dir_in, dir_out, sftwr = 2)
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' print(files)
#' file.remove(files)
#'
gcims_read_samples <- function(dir_in, dir_out, sftwr) {
  stopifnot(is.character(dir_out), length(dir_out) == 1)

  metadata <- list(Name = NULL)
  data <- list(retention_time = NULL, drift_time = NULL, data_df = NULL)

  if (sftwr == 1){
    files <- list.files(path = dir_in, pattern = ".csv")
    if (length(files) == 0){
      stop("This folder does not contains any .csv file")
    }

    if (!dir.exists(dir_out)) {
      dir.create(dir_out, recursive = TRUE)
    }

    for (i in seq_along(files)){
      aux_string <- paste0("M", i, ".rds")
      # METADATA
      metadata$Name <- i
      # DATA
      single_file <- readr::read_csv(
        file.path(dir_in, files[i]),
        progress = FALSE, skip = 1, col_names = FALSE,
        col_types = readr::cols(.default = "c")
      )

      dd <- as.data.frame(t(single_file[-1, -1]))
      # Depende del tipo de instrumento OJO
      data$data_df <- readr::type_convert(dd, col_types = readr::cols(.default = "d"))                    # Signo Menos solo para gasdormund programa viejo
      data$retention_time <- readr::type_convert(single_file[-1, 1] , col_types = readr::cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
      data$retention_time <-  data$retention_time[[1]]
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <- as.numeric(single_file[1, -(dim(single_file)[2])])#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.

      # Join
      dd_list <- list(metadata = metadata, data = data)
      saveRDS(dd_list, file = file.path(dir_out, paste0("M", i, ".rds")))
    }
  }
  else {
    files <- list.files(path = dir_in, pattern = ".csv")
    if (length(files) == 0){
      stop("This folder does not contains any .csv file")
    }

    for (i in seq_along(files)){
      aux_string <- paste0("M", i, ".rds")

      # METADATA
      metadata$Name <- i
      # DATA
      single_file <- readr::read_csv(file.path(dir_in, files[i]), skip = 130,
                              progress = FALSE, col_names = FALSE,
                              col_types = readr::cols(.default = "c"))

      dd <- single_file[-c(1:2), -c(1:2)]
      # Depende del tipo de instrumento OJO
      data$data_df = - readr::type_convert(dd, col_types = readr::cols(.default = "d"))                    # Signo Menos solo para gasdormund programa nuevo
      data$retention_time <- as.numeric(single_file[1, -c(1:2)])
      # Depende del tipo de instrumento OJO
      #print(str(single_file))
      data$drift_time <- readr::type_convert(single_file[-c(1:2), 2] , col_types = readr::cols(.default = "d"))#[-c(1:2), 2]       # el primer indexado es para sacar de la lista. El segundo para los datos.
      data$drift_time <-  data$drift_time[[1]]                                                # Depende del tipo de instrumento OJO

      # Join
      dd_list <- list(metadata = metadata, data = data)
      saveRDS(dd_list, file = file.path(dir_out, paste0("M", i, ".rds")))
    }
  }
}


#' Read GC-IMS samples from .mea or .mea.gz files and save them as RDS list objects
#' @param dir_in          Input directory. Where the raw data files are stored.
#' @param dir_out         Output directory. Where data files are stored as R
#'   objects.
#'
#' @details `gcims_read_samples()` stores a set of S3 objects in the
#'   directory `dir_out` (one per sample sample stored in `dir_in`).
#'   Each object is a list containing the following variables: `metadata`
#'   and `data`.In `metadata` you can find information related to the
#'   sample that can be added by the user. On the other hand in `data` you
#'   can find information of the measurement.By default, when the
#'   `metadata` field is created only includes one variable, called
#'   `Name`, although the user can include more fields using the function
#'   [gcims_read_metadata()]. Regarding `data`, it is consist in
#'   three fields:  \describe{
#'      \item{`retention_time`}{ The vector of retention times associated to a sample measurement.}
#'      \item{`drift_time`}{ The vector of drift times associated to a sample
#'   measurement.}
#'      \item{`data_df`}{ A dataframe containing the intensities
#'        of a GCIMS measurement of a sample. Each column in the dataset corresponds
#'        to a particular drift time, while each row to a different retention time.}
#'   }
#'
#' @return The filenames that have been created in `dir_out`
#' @family Reading functions.
#' @export
#' @examples
#' dir_in <- system.file("extdata", "sample_formats", package = "GCIMS")
#' dir_out <- tempfile("dir_out")
#' dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
#' gcims_read_mea(dir_in, dir_out)
#' list.files(dir_out)
#'
gcims_read_mea <- function(dir_in, dir_out) {
  files <- list.files(path = dir_in, pattern = "(\\.mea(\\.gz)?)$", full.names = TRUE)
  outfiles <- c()
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(files)) {
    metadata <- list(Name = NULL)
    data <- list(retention_time = NULL, drift_time = NULL, data_df = NULL)

    aux_string <- paste0("M", i, ".rds")
    single_mea <- read_mea(files[i])

    # Metadata
    metadata$Name <- i
    # Data
    data$drift_time <- dtime(single_mea)
    data$retention_time <- rtime(single_mea)
    data$data_df <- intensity(single_mea)

    # Join
    dd_list <- list(metadata = metadata, data = data)
    outfile <- file.path(dir_out, aux_string)
    saveRDS(dd_list, file = outfile)
    outfiles <- c(outfiles, outfile)
  }
  names(outfiles) <- files
  outfiles
}


#' Read GC-IMS Samples from MATLAB files
#'
#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @return An invisible NULL. This function is called to convert mat files to rds files.
#' @family Reading functions
#' @export
#' @examples
#' # We can create a mat file:
#' mat_samples_dir <- tempfile("mat_samples")
#' dir.create(mat_samples_dir)
#' matfile <- file.path(mat_samples_dir, "sample.mat")
#'
#' obj <- list(
#'   metadata = list(Name = 1),
#'   data = list(
#'     drift_time = 1:3,
#'     retention_time = 1:4,
#'     data_df = matrix(1:12, nrow = 3, ncol = 4)
#'   )
#' )
#' R.matlab::writeMat(matfile, dd = obj)
#' dir_out <- tempfile("rds_samples")
#' dir.create(dir_out)
#' gcims_read_mat(mat_samples_dir, dir_out)
#' list.files(dir_out, full.names = TRUE)
gcims_read_mat <- function(dir_in, dir_out) {
  require_pkgs("R.matlab")

  files <- list.files(path = dir_in, pattern = "\\.mat$")
  if (length(files) == 0){
    stop("This folder does not contains any .mat file")
  }

  for (i in seq_along(files)) {

    dd <- R.matlab::readMat(file.path(dir_in, files[i]))[[1]]
    dd <- list(
      metadata = list(Name = i),
      data = dd
    )
    if (!dir.exists(dir_out)) {
      dir.create(dir_out)
    }
    saveRDS(dd, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
  invisible(NULL)
}


#' Read metadata from a csv

#' @param dir_in          Input directory. Where data files are stored.
#' @param samples         Vector of positive integers. Set of sample
#'   identifiers.
#' @param file            Name of the file that contains the metadata (an Excel
#'   file).
#' @details `gcims_read_metadata` reads a metadata Excel file. This file,
#' contains:
#' \itemize{
#'   \item A number of named columns. First column is always named 'Name'. The rest of the columns names
#'   correspond to new fields that must be included as metadata in the S3 object.
#'   \item A number of rows, each of them associated to a sample in the dataset. Thus in a row, you can find
#'   all metadata information of a sample in the dataset.
#' }
#' After that `gcims_read_metadata` looks for the match between the S3
#' object metadata field 'Name' value, and its corresponding value
#' found in column 'Name' of the metadata file. Then metadata information is
#' copied from the metadata file to the S3 object (column names in the metadata
#' file define the new fields of metadata in the S3 object). This is done for
#' all the samples of the dataset.
#' @return A set of S3 objects.
#' @family Reading functions
#' @note The metadata file has to be in the same directory as the data samples.
#' @export
#' @examples
#' # Prepare a directory:
#' dir_in <- tempfile("dir_in")
#' dir.create(dir_in, recursive = TRUE)
#' file.copy(
#'   system.file("extdata","add_meta", "Metadata.xlsx", package = "GCIMS"),
#'   dir_in
#' )
#' file.copy(
#'   system.file("extdata","add_meta", "M1.rds", package = "GCIMS"),
#'   dir_in
#' )
#' M1 <- readRDS(file.path(dir_in, "M1.rds"))
#' print(M1$metadata)
#' samples <- 1
#' file <- "Metadata.xlsx"
#' gcims_read_metadata(dir_in, samples, file)
#' M1 <- readRDS(file.path(dir_in, "M1.rds"))
#' print(M1$metadata)
#' M1$metadata$Class <- NULL
#' M1$metadata$Bottle <- NULL
#' M1$metadata$Brand <- NULL
#' print(M1$metadata)
#'
gcims_read_metadata <- function(dir_in, samples, file) {

  if (is.character(file)) {
    Metadatafile <- readxl::read_excel(file.path(dir_in, file))
  } else if (is.data.frame(file)) {
    Metadatafile <- file
  } else {
    stop("file must be an Excel file name")
  }

  metadata <- NULL
  for (i in seq_along(samples)){
    aux_string <- file.path(dir_in, paste0("M", samples[i], ".rds"))
    aux_list <- readRDS(aux_string)
    metadata <- Metadatafile[which(Metadatafile$Name == aux_list$metadata$Name),]
    aux_list$metadata <- metadata
    saveRDS(aux_list, file = aux_string)
  }
}
