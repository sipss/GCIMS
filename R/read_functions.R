#' Read GC-IMS samples from .csv files


#' @param dir_in          Input directory. Where the raw data files are stored.
#' @param dir_out         Output directory. Where data files are stored as R
#'   objects.
#' @param sftwr           Selection of the software version. Use 1 if the
#'   software that has converted the original files into csv format is VOCal.
#'   Use 2 if the software is LAV.
#'
#' @details \code{gcims_read_samples} stores a set of S3 objects in the
#'   directory \code{dir_out} (one per sample sample stored in \code{dir_in}).
#'   Each object is a list containing the following variables: \code{metadata}
#'   and \code{data}.In \code{metadata} you can find information related to the
#'   sample that can be added by the user. On the other hand in \code{data} you
#'   can find information of the measurement.By default, when the
#'   \code{metadata} field is created only includes one variable, called
#'   \code{Name}, although the user can include more fields using the function
#'   \code{\link{gcims_read_metadata}}. Regarding \code{data}, it is consist in
#'   three fields:  \describe{ \item{\code{retention_time}}{ The vector of
#'   retention times associated to a sample measurement.}
#'   \item{\code{drift_time}}{ The vector of drift times associated to a sample
#'   measurement.} \item{\code{data_df}}{ A dataframe containing the intensities
#'   of a GCIM measurement of a sample. Each column in the dataset corresponds
#'   to a particular drift time, while each row to a different retention time.}
#'   }
#'
#' @return A set of S3 objects.
#' @note In the current version of the package, only .csv files converted from
#'   LAV or VOCal software of G.A.S GmbH can be used by the function
#'   \code{gcims_read_samples}.
#' @family Reading functions.
#' @export
#' @importFrom readr read_csv cols type_convert
#' @importFrom utils menu
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' gcims_read_samples(dir_in, dir_out, sftwr = 2)
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' print(files)
#' file.remove(files)
#' setwd(current_dir)
#'
gcims_read_samples <- function(dir_in, dir_out, sftwr) {
  setwd(dir_in)
  print(" ")
  print("  /////////////////////////")
  print(" /    Reading Samples    /")
  print("/////////////////////////")
  print(" ")

  metadata <- list(Name = NULL)
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
      metadata$Name <- i
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
      metadata$Name <- i
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


#' Read .mea files (from GAS Dortmund)
#'
#' This function reads a .mea file (supporting gzip compressed .mea.gz files as well)
#' and returns a GCIMS object
#'
#' @param filename A .mea or a .mea.gz path to a file
#'
#' @return A GC-IMS sample
#' @export
#'
#' @examples
#' mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
#' sample <- read_mea(mea_file)
#'
read_mea <- function(filename) {
  string_keys <- c(
    "ADIO gpident no", "ADIO name", "ADIO serial", "ADIO version",
    "Class", "Firmware date", "Firmware version", "Machine name",
    "Machine serial", "Machine type", "Recognized substances", "Sample",
    "Sample number", "Program", "Status", "Timestamp", "Drift Gas",
    "Filter", "GC Column", "Sensor data")
  unitless_integer_keys <- c('Chunk averages', 'Chunk sample count', "Chunks count")
  value_unit_keys <- c(
    "Board temperature", 'Chunk sample rate',
    'Chunk trigger duration', 'Chunk trigger repetition',
    "Chunk voltrange", "Flow1 setpoint", "Flow2 setpoint",
    "Start Flow1", "Start flow2", "Start temp1",
    "Start temp2", "Start temp3", "Temp1 setpoint",
    "Temp2 setpoint", "Temp3 setpoint",
    "EPC ambient pressure", "EPC1 end-pressure", "EPC1 pressure", "EPC2 end-pressure",
    "EPC2 pressure", "Flow Epc 1", "Flow Epc 2", "Flow record interval",
    "nom Drift Potential Difference", "nom Drift Tube Length", "Pressure Ambient",
    "Pressure Epc 1", "Pressure Epc 2", "Sensor block", "Sensor drift", "Sensor inject",
    "Start flow1", "Start flow2", sprintf("Start temp %d", 1:6), sprintf("Temp %d setpoint", 1:6),
    "Pressure record interval")
  last_key <- "Timestamp"

  params <- list()

  if (endsWith(filename, ".gz")) {
    con <- gzfile(filename, "rb")
  } else {
    con <- file(filename, "rb")
  }
  on.exit(close(con))

  key <- ""
  while (key != last_key) {
    tline <- readLines(con, n=1, encoding = "windows-1252")
    tline <- iconv(tline, from = "windows-1252", to="UTF-8")
    key_values <- stringr::str_trim(stringr::str_split_fixed(tline, "=", n=2)[1,], side = "both")
    key <- key_values[1]
    value <- key_values[2]

    if (key %in% string_keys) {
      if (nchar(value) > 1 && startsWith(value, '"') && endsWith(value, '"')) {
        value <- substr(value, start = 2, stop = nchar(value)-1)
      }
      params[[key]] <- value
    } else if (key %in% unitless_integer_keys) {
      params[[key]] <- as.integer(value)
    } else if (key %in% value_unit_keys) {
      pat <- stringr::regex("(.*)\\[(.*)\\]")
      val_unit <- stringr::str_trim(stringr::str_match(value, pat)[1,], side = "both")
      value <- val_unit[2]
      value_as_num <- stringr::str_split(value, stringr::fixed(" "))[[1]]
      suppressWarnings({
        value_as_num <- as.numeric(value_as_num)
      })
      if (!anyNA(value_as_num)) {
        value <- value_as_num
      }
      unit <- val_unit[3]
      params[[key]] <- list(value=value, unit=unit)
    } else {
      warning(sprintf("Unknown key: %s (please open an issue to implement this)", key))
      params[[key]] <- value
    }
  }
  # drift time in ms, we want chunk sample rate in kHz
  if (params[["Chunk sample rate"]][["unit"]] == "kHz") {
    drift_time_sample_rate_khz <- params[["Chunk sample rate"]][["value"]]
  } else {
    stop(sprintf("Expected Chunk sample rate to be in kHz, found %s instead. Please open an issue to implement this", params[["Chunk sample rate"]][["unit"]]))
  }
  drift_time <- seq(from=0.0, by=1./drift_time_sample_rate_khz, length.out=params[["Chunk sample count"]])
  # create retention time vector. Example:
  # time between two spectra measurements (in sec)
  # we are averaging 32 spectra, and one spectra is taken every 21 ms:
  # (32+1)*21/1000 = 0.692 s.
  # FIXME: Ask GAS?: To match the retention time from LAV 2.0.0, we must
  # multiply the (number of averages +1) by the trigger repetition time.
  # I do not understand where the +1 comes from, but they are the ground truth.
  # Check what happens with VOCal.
  if (params[["Chunk trigger repetition"]][["unit"]] == "ms") {
    retention_time_step_s <- params[["Chunk trigger repetition"]][["value"]] / 1000
  } else if (params[["Chunk trigger repetition"]][["unit"]] == "s") {
    retention_time_step_s <- params[["Chunk trigger repetition"]][["value"]]
  } else {
    stop(sprintf("Expected Chunk trigger repetition to be in ms, found %s instead. Please open an issue to implement this", params[["Chunk sample rate"]][["unit"]]))
  }
  ret_time_step <- (params[['Chunk averages']]+1)*retention_time_step_s
  ret_time <- seq(from = 0.0, by = ret_time_step, length.out = params[["Chunks count"]])
  # read the actual data:
  npoints <- params[["Chunks count"]]*params[["Chunk sample count"]]
  data <- readBin(con, what = "int", n = npoints, signed = TRUE, size = 2L, endian = "big")
  if (length(data) != npoints) {
    stop(sprintf("binary data from mea file had %d points, and %d were expected", length(data), npoints))
  }
  data <- matrix(data, nrow = params[["Chunk sample count"]], ncol = params[["Chunks count"]])
  last_byte_is_a_zero <- readBin(con, what = "integer", n = 1, size = 1L)
  if (last_byte_is_a_zero != 0L) {
    stop("read_mea expected one zero byte at the end of the data stream.")
  }
  # FIXME: Filter data so we get the same values as LAV / VOCal
  # The data in the mea file has some values that differ by multiples of 256
  # from the corresponding values found in the csv files
  list(drift_time = drift_time, ret_time = ret_time, data = data, params = params)
}


#' Read GC-IMS Samples from MATLAB files
#'
#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @return NULL
#' @family Reading functions
#' @export
#' @examples
#' \dontrun{
#' gcims_read_mat("/path/to/dir/with/mat/files", "where_to_save_samples_in_rds/")
#' }

gcims_read_mat <- function(dir_in, dir_out) {
  print(" ")
  print("  /////////////////////////")
  print(" /    Reading Samples    /")
  print("/////////////////////////")
  print(" ")

  files <- list.files(path = dir_in, pattern = ".mat$")
  if (length(files) == 0){
    stop("This folder does not contains any .mat file")
  }
  for (i in seq_along(files)) {
    print(paste0("Sample ", i, " of ", length(files)))
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
}


#' Read metadata from a csv

#' @param dir_in          Input directory. Where data files are stored.
#' @param samples         Vector of positive integers. Set of sample
#'   identifiers.
#' @param file            Name of the file that contains the metadata (an Excel
#'   file).
#' @details \code{gcims_read_metadata} reads a metadata Excel file. This file,
#' contains:
#' \itemize{
#'   \item A number of named columns. First column is always named 'Name'. The rest of the columns names
#'   correspond to new fields that must be included as metadata in the S3 object.
#'   \item A number of rows, each of them associated to a sample in the dataset. Thus in a row, you can find
#'   all metadata information of a sample in the dataset.
#' }
#' After that \code{gcims_read_metadata} looks for the match between the S3
#' object metadata field 'Name' value, and its corresponding value
#' found in column 'Name' of the metadata file. Then metadata information is
#' copied from the metadata file to the S3 object (column names in the metadata
#' file define the new fields of metadata in the S3 object). This is done for
#' all the samples of the dataset.
#' @return A set of S3 objects.
#' @family Reading functions
#' @note The metadata file has to be in the same directory as the data samples.
#' @export
#' @importFrom readxl read_excel
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata","add_meta", package = "GCIMS")
#' setwd(dir_in)
#' M1 <- readRDS("M1.rds")
#' samples <- 1
#' print(M1$metadata)
#' file <- "Metadata.xlsx"
#' gcims_read_metadata(dir_in, samples, file)
#' M1 <- readRDS("M1.rds")
#' print(M1$metadata)
#' M1$metadata$Class <- NULL
#' M1$metadata$Bottle <- NULL
#' M1$metadata$Brand <- NULL
#' print(M1$metadata)
#' saveRDS(M1, "M1.rds")
#' setwd(current_dir)
#'
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
    metadata <- Metadatafile[which(Metadatafile$Name == aux_list$metadata$Name),]
    aux_list$metadata <- metadata
    saveRDS(aux_list, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }
}

