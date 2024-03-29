#' Read .mea files (from GAS Dortmund)
#'
#' This function reads a .mea file (supporting gzip compressed .mea.gz files as well)
#' and returns a GCIMS object
#'
#' @param filename A .mea or a .mea.gz path to a file
#'
#' @return The GC-IMS sample in a [GCIMSSample-class] object
#'
#' @details
#'
#' Thanks to Daniel Sanders and Thomas Wortelmann from  [GAS Dortmund](https://www.gas-dortmund.de/) for providing the
#' necessary support to implement this function.
#'
#' @examples
#' mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
#' sample <- read_mea(mea_file)
#'
#' @export
read_mea <- function(filename) {
  string_keys <- c(
    "ADIO gpident no", "ADIO name", "ADIO serial", "ADIO version",
    "Class", "Firmware date", "Firmware version", "Machine name",
    "Machine serial", "Machine type", "Recognized substances", "Sample",
    "Sample number", "Program", "Status", "Timestamp", "Drift Gas",
    "Filter", "GC Column", "Sensor data"
  )
  unitless_integer_keys <- c(
    'Chunk averages', 'Chunk sample count',
    "Chunks count", "Pump power"
  )
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
    "Start flow1", "Start flow2",
    sprintf("Start temp %d", seq.int(1,6)),
    sprintf("Temp %d setpoint", seq.int(1,6)),
    "Pressure record interval"
  )

  filename_l <- tolower(filename)
  if (!endsWith(filename_l, ".mea.gz") && !endsWith(filename_l, ".mea")) {
    cli_abort("read_mea expects either a .mea or .mea.gz file")
  }
  params <- list()


  tryCatch({
    if (endsWith(filename, ".gz")) {
      con <- gzfile(filename, "rb")
    } else {
      con <- file(filename, "rb")
    }
    # Read all file in memory (chunks of 64MB, so hopefully we'll have enough with one)
    sixtyfour_mb <- 64*1024*1024
    read_file_in_memory <- list(readBin(con = con, what = "raw", n = sixtyfour_mb))
    while (length(read_file_in_memory[[length(read_file_in_memory)]]) == sixtyfour_mb) {
      read_file_in_memory <- append(read_file_in_memory, list(readBin(con = con, what = "raw", n = sixtyfour_mb)))
    }
    if (length(read_file_in_memory) == 1) {
      read_file_in_memory <- read_file_in_memory[[1]]
    } else {
      read_file_in_memory <- purrr::flatten_raw(read_file_in_memory)
    }
  }, finally = {
    close(con)
  })

  # Returns the index of the first time "a" appears in "x".
  # efficient for raw vectors
  match_raw <- function(a, x, slice_size = 8192L) {
    offset <- 0
    x_len <- length(x)
    start_idx <- 1
    end_idx <- slice_size
    while (TRUE) {
      comp <- which(x[start_idx:end_idx] == a)
      if (length(comp) != 0) {
        return(offset + comp[1])
      } else {
        offset <- offset + slice_size
        start_idx <- start_idx + slice_size
        end_idx <- end_idx + slice_size
        if (start_idx > length(x)) {
          return(NA_integer_)
        } else if (end_idx > length(x)) {
          end_idx <- length(x)
        }
      }
    }
  }

  # Find where the header ends
  # match(as.raw(0L), read_file_in_memory) is painfully slow due to internal character coercions.
  head_data_sep <- match_raw(as.raw(0L), read_file_in_memory, slice_size = 8192L)
  header_as_text <- iconv(
    rawToChar(utils::head(read_file_in_memory, n = head_data_sep - 1)),
    from = "windows-1252",
    to = "UTF-8"
  )

  key_values_list <- strsplit(header_as_text, split = "\n", fixed = TRUE)[[1]]
  key_values_matrix <- stringr::str_split_fixed(key_values_list, "=", n = 2)
  key_values_matrix <- matrix(stringr::str_trim(key_values_matrix, side = "both"), ncol = 2)
  # Parse header:
  for (i in seq_len(nrow(key_values_matrix))) {
    key <- key_values_matrix[i, 1]
    value <- key_values_matrix[i, 2]

    if (key %in% string_keys) {
      if (nchar(value) > 1 && startsWith(value, '"') && endsWith(value, '"')) {
        value <- substr(value, start = 2, stop = nchar(value) - 1)
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
      params[[key]] <- list(value = value, unit = unit)
    } else {
      cli_warn(
        message = "Unknown key: {key} (please open an issue to implement this)",
        .frequency = "once",
        .frequency_id = glue("GCIMS-warn-unknown-key-{key}")
      )
      params[[key]] <- value
    }
  }

  # Convert and shape the data
  npoints <- params[["Chunks count"]]*params[["Chunk sample count"]]
  data <- readBin(
    read_file_in_memory[(head_data_sep + 1):length(read_file_in_memory)],
    what = "integer",
    n = npoints,
    signed = TRUE,
    size = 2,
    endian = "little"
  )
  data <- matrix(data, nrow = params[["Chunk sample count"]], ncol = params[["Chunks count"]])

  # Build drift time, retention time, parse gc_column...

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
  # After a great discussion with GAS, the +1 exists for backwards compatibility and it is not going to change.
  if (params[["Chunk trigger repetition"]][["unit"]] == "ms") {
    retention_time_step_s <- params[["Chunk trigger repetition"]][["value"]] / 1000
  } else if (params[["Chunk trigger repetition"]][["unit"]] == "s") {
    retention_time_step_s <- params[["Chunk trigger repetition"]][["value"]]
  } else {
    stop(sprintf("Expected Chunk trigger repetition to be in ms, found %s instead. Please open an issue to implement this", params[["Chunk sample rate"]][["unit"]]))
  }
  ret_time_step <- (params[['Chunk averages']] + 1)*retention_time_step_s
  ret_time <- seq(from = 0.0, by = ret_time_step, length.out = params[["Chunks count"]])

  drift_tube_length <- numeric(0)
  if (!is.null(params[["nom Drift Tube Length"]])) {
    drift_tube_length <- params[["nom Drift Tube Length"]][["value"]]
    if (params[["nom Drift Tube Length"]][["unit"]] == "\u00b5m") {
      drift_tube_length <- drift_tube_length / 1000
    } else if (params[["nom Drift Tube Length"]][["unit"]] == "mm") {
      # do nothing
    } else {
      stop(sprintf("Report this: Unit conversion %s->mm not implemented in drift tube length", params[["nom Drift Tube Length"]][["unit"]]))
    }
  }

  gc_column <- ifelse(!is.null(params[["GC Column"]]), params[["GC Column"]], character(0L))
  drift_gas <- ifelse(!is.null(params[["Drift Gas"]]), params[["Drift Gas"]], character(0L))

  if (!is.null(names(filename))) {
    description <- names(filename)
  } else {
    description <- tools::file_path_sans_ext(basename(filename), compression = TRUE)
  }

  GCIMSSample(
    drift_time = drift_time,
    retention_time = ret_time,
    data = data,
    gc_column = gc_column,
    drift_tube_length = drift_tube_length,
    drift_gas = drift_gas,
    history = sprintf("Sample loaded from %s", filename),
    filepath = filename,
    description = description,
    params = params
  )
}


#' write a mea file (experimental)
#'
#' This function writes a GCIMSSample into a .mea file. The implementation
#' misses a lot of metadata, and may not be fully compliant with LAV or VOCal,
#' so it is not recommended in production environments.
#'
#' @param object A `GCIMSSample` object
#' @param filename The filename to write it to (will have .mea appended if no extension is given)
#'
#' @return Nothing
#' @noRd
#' @examples
#' obj <- GCIMSSample(
#'  drift_time = 1:2,
#'  retention_time = 1:3,
#'  data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' mea_file <- tempfile(fileext = ".mea")
#' write_mea(obj, mea_file)
#' file.remove(mea_file)
write_mea <- function(object, filename) {
  if (!inherits(filename, "character")) {
    stop("filename should be a string")
  }
  if (length(filename) != 1) {
    stop("filename should be of length one")
  }
  if (endsWith(filename, ".mea.gz")) {
    con <- gzfile(filename, "wb")
  } else {
    if (!endsWith(filename, ".mea")) {
      filename <- paste0(filename, ".mea")
    }
    con <- file(filename, "wb")
  }
  on.exit(close(con))

  dt <- dtime(object)
  rt <- rtime(object)
  dtime_rate_khz <- 1./(dt[2] - dt[1])
  rtime_rate_hz <-  1./(rt[2] - rt[1])
  fmt <- c(
    "key_value_unit" = "%-31s= %f [%s]",
    "key_integer" = "%-31s= %d",
    "key_string" = '%-31s= "%s"'
  )
  FAKE_CHUNK_AVGS <- 1
  FAKE_TIMESTAMP <- "1970-01-01T00:00:00"
  header <- c(
    sprintf(fmt["key_string"], "ADIO gpident no", "5551D615"), #fake
    sprintf(fmt["key_string"], "ADIO name", "ADIO TYP02"), #fake
    sprintf(fmt["key_string"], "ADIO serial", "ADIO-10010"), #fake
    sprintf(fmt["key_string"], "ADIO version", "V1.20"), #fake
    sprintf(fmt["key_value_unit"], "Board temperature", 34, "\u00B0C"), #fake
    sprintf(fmt["key_integer"], "Chunk averages", FAKE_CHUNK_AVGS), #fake
    sprintf(fmt["key_integer"], "Chunk sample count", length(dt)),
    sprintf(fmt["key_value_unit"], "Chunk sample rate", dtime_rate_khz, "kHz"),
    sprintf(fmt["key_value_unit"], "Chunk trigger duration", 100, "\u00B5s"), #fake
    sprintf(fmt["key_value_unit"], "Chunk trigger repetition", 1000/rtime_rate_hz/(FAKE_CHUNK_AVGS + 1), "ms"),
    sprintf(fmt["key_value_unit"], "Chunk voltrange", 10.000, "V"),
    sprintf(fmt["key_integer"], "Chunks count", length(rt)),
    sprintf(fmt["key_string"], "Class", "program/pos"), # required
    sprintf(fmt["key_string"], "Drift Gas", "nitrogen"),
    sprintf(fmt["key_value_unit"], "EPC ambient pressure", 100.516, "kPa"), # fake
    sprintf(fmt["key_value_unit"], "EPC1 end-pressure", 0.503, "kPa"), # fake
    sprintf(fmt["key_value_unit"], "EPC1 pressure",0.500, "kPa"), # fake
    sprintf(fmt["key_value_unit"], "EPC2 end-pressure", 164.430, "kPa"), # fake
    sprintf(fmt["key_value_unit"], "EPC2 pressure", 141.011, "kPa"), # fake
    sprintf(fmt["key_string"], "Filter", "SG8"),
    sprintf(fmt["key_string"], "Firmware date", "2019-08-06"),
    sprintf(fmt["key_string"], "Firmware version", "2.54"),
    sprintf(fmt["key_string"], "Machine name", "GAScontrol"),
    sprintf(fmt["key_string"], "Machine serial", "1H1-00071"),
    sprintf(fmt["key_string"], "Machine type", "FlavourSpec\u00AE"),
    #sprintf(fmt["key_string"], "ADIO gpident no", "5551D615"),
    #sprintf(fmt["key_string"], "ADIO gpident no", "5551D615"),
    sprintf(fmt["key_string"], "Timestamp", FAKE_TIMESTAMP),
    ""
  )
  binheader <- iconv(paste0(header, collapse = "\n"), from = "UTF-8", to = "windows-1252", toRaw = TRUE)[[1]]
  writeBin(binheader, con = con, useBytes = TRUE)
  writeBin(object = as.raw(0), con = con, size = 1)
  writeBin(as.integer(intensity(object)),
           con, size = 2L, endian = "little", useBytes = FALSE)
  invisible(NULL)
}
