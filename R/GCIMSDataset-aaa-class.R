#' GCIMSDataset class
#'
#' GCIMSDataset is an S4 class to store a dataset
#'
#' The actual data is not stored in memory, but read/saved from/to files as
#' needed, so the dataset object scales with large number of samples.
#'
#' @export
#' @examples
#' # Create a new GCIMSDataset with methods::new()
#' dummy_obj <- methods::new(
#'   "GCIMSDataset",
#'   pData = data.frame(SampleID = character(), filename = character(0)),
#'   base_dir = tempdir()
#' )
#' @importClassesFrom S4Vectors DataFrame
# All the slots are in an environment because I want the dataset to be mutable.
# @slot pData A data frame with at least the SampleID and filename columns.
# @slot base_dir A directory containing the file names described in `pData`
# @slot scratch_dir A directory to save intermediate results.
# @slot delayed_ops Delayed operations
# @slot TIS A matrix of n_samples vs drift time, with the Total Ion Spectrum of each sample
# @slot RIC A matrix of n_samples vs retention time, with the Reverse Ion Chromatogram of each sample
# @slot dt_ref A numeric drift time of reference
# @slot rt_ref A numeric retention time of reference
methods::setClass(
  Class = "GCIMSDataset",
  slots = c(
    # pData = "DataFrame",
    # base_dir = "character",
    # scratch_dir = "character",
    # delayed_ops = "list",
    # TIS = "matrixOrNULL",
    # RIC = "matrixOrNULL",
    # dt_ref = "numericOrNULL",
    # rt_ref = "numericOrNULL",
    envir = "environment"
  )
)

check_files <- function(filenames, base_dir) {
  files_exist <- purrr::map_lgl(
    filenames,
    function(filename, base_dir) {
      full_path <- file.path(base_dir, filename)
      file.exists(full_path)
    },
    base_dir = base_dir
  )
  missing_files <- filenames[!files_exist]
  if (length(missing_files) > 0) {
    if (length(missing_files) > 10) {
      to_show <- 5
      and_n_more <- paste0(" (and ", length(missing_files) - to_show, " more)")
    } else {
      to_show <- length(missing_files)
      and_n_more <- ""
    }
    rlang::abort(
      message = c(
        "Files not found",
        "x" = paste0("The following files", and_n_more, " were not found"),
        "i" = utils::head(missing_files, to_show)
      )
    )
  }

}

abort_if_errors <- function(errors, title = "Errors found") {
  if (length(errors) > 0) {
    names(errors) <- rep("x", length(errors))
    rlang::abort(message = c(title, errors))
  }
}

validate_pData <- function(pData) {
  pData <- S4Vectors::DataFrame(pData)
  pheno_names <- tolower(colnames(pData))
  sampleid_col <- which(pheno_names == "sampleid")
  filename_col <- which(pheno_names == "filename")
  errors <- character(0L)
  if (length(sampleid_col) == 0) {
    errors <- c(errors, "pData should have a SampleID column with a sample name")
  }
  if (length(sampleid_col) > 1) {
    errors <- c(errors, "pData should only have one SampleID column (don't use sampleid or other casing)")
  }
  if (length(filename_col) == 0) {
    errors <- c(errors, "pData should have a FileName column with the paths to the sample files")
  }
  if (length(filename_col) > 1) {
    errors <- c(errors, "pData should only have one FileName column (don't use filename or other casing)")
  }
  abort_if_errors(errors, title = "pData is not valid")

  # CaSiNg MaY bE wRoNg, sO We sEt iT juSt in CaSe.
  colnames(pData)[filename_col] <- "FileName"
  colnames(pData)[sampleid_col] <- "SampleID"
  # Return corrected pData, with SampleID on the first column and FileName on the second column
  pData <- pData[, c("SampleID", "FileName", setdiff(colnames(pData), c("SampleID", "FileName")))]
  pData
}

validate_base_dir <- function(base_dir) {
  errors <- character(0L)
  if (!rlang::is_string(base_dir)) {
    errors <- c(errors, "base_dir should be a string")
  }
  if (!dir.exists(base_dir)) {
    errors <- c(errors, "base_dir does not exist")
  }
  abort_if_errors(errors, title = "base_dir is not valid")
  normalizePath(base_dir, mustWork = TRUE)
}

validate_scratch_dir <- function(scratch_dir) {
  errors <- character(0L)
  if (!rlang::is_string(scratch_dir)) {
    errors <- c(errors, "scratch_dir should be a string")
  }
  dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(scratch_dir)) {
    errors <- c(errors, "scratch_dir does not exist and could not be created")
  }
  abort_if_errors(errors, title = "scratch_dir is not valid")
  normalizePath(scratch_dir, mustWork = TRUE)
}

validate_keep_intermediate <- function(keep_intermediate) {
  keep_intermediate <- as.logical(keep_intermediate)
  if (is.na(keep_intermediate)) {
    rlang::abort("keep_intermediate must be either TRUE or FALSE")
  }
  keep_intermediate
}

methods::setMethod(
  "initialize", "GCIMSDataset",
  function(.Object, pData, base_dir, scratch_dir = tempfile("GCIMSDataset_tempdir_"), keep_intermediate = FALSE) {
    pData <- validate_pData(pData)
    base_dir <- validate_base_dir(base_dir)
    check_files(pData$FileName, base_dir)
    scratch_dir <- validate_scratch_dir(scratch_dir)
    keep_intermediate <- validate_keep_intermediate(keep_intermediate)
    .Object@envir <- rlang::new_environment()
    .Object@envir$pData <- pData
    .Object@envir$base_dir <- base_dir
    .Object@envir$scratch_dir <- scratch_dir
    .Object@envir$dt_ref <- NULL
    .Object@envir$rt_ref <- NULL
    .Object@envir$TIS <- NULL
    .Object@envir$RIC <- NULL
    .Object@envir$delayed_ops <- NULL
    .Object@envir$previous_ops <- list()
    .Object@envir$hasheddir <- ""
    .Object@envir$keep_intermediate <- keep_intermediate
    canRealize(.Object) <- TRUE
    .Object <- appendDelayedOp(.Object, GCIMSDelayedOp(name = "read_sample", fun = read_sample))
    # Some sample stats:
    .Object <- extract_dtime_rtime(.Object)
    .Object <- extract_RIC_and_TIS(.Object)
    .Object
  }
)

CurrentHashedDir <- function(object) {
  hasheddir <- object@envir$hasheddir
  if (hasheddir == "") {
    return(NULL)
  }
  file.path(
    object@envir$scratch_dir,
    hasheddir
  )
}

"CurrentHashedDir<-" <- function(object, value) {
  object@envir$hasheddir <- basename(value)
  object
}


NextHashedDir <- function(object) {
  if (!hasDelayedOps(object)) {
    return(CurrentHashedDir(object))
  }
  # A hash chain:
  current_hash <- object@envir$hasheddir
  if (current_hash == "") {
    # The first hash comes from the base_dir, sampleID and FileNames
    pd <- pData(object)
    if (!all(c("SampleID", "FileName") %in% colnames(pd))) {
      rlang::abort(c("Unexpected Error", "x" = "Expected pData with SampleID and FileName columns"))
    }
    current_hash <- digest::digest(
      list(
        object@envir$base_dir,
        pd[["SampleID"]],
        pd[["FileName"]]
      )
    )
  }
  for (op in purrr::keep(object@envir$delayed_ops, modifiesSample)) {
    current_hash <- paste0(name(op), "_", digest::digest(list(current_hash, hashableDelayedOp(op))))
  }
  file.path(
    object@envir$scratch_dir,
    current_hash
  )
}


read_sample <- function(filename) {
  filename_l <- tolower(filename)
  if (endsWith(filename_l, ".mea.gz") || endsWith(filename_l, ".mea")) {
    return(read_mea(filename))
  } else if (endsWith(filename_l, ".rds")) {
    obj <- readRDS(filename)
    if (methods::is(obj, "GCIMSSample")) {
      return(obj)
    } else if (methods::is(obj, "list")) {
      # FIXME: Remove this if we drop the old API
      #Convert old list to GCIMSSample:
      sample <- GCIMSSample(
        drift_time = obj$data$drift_time,
        retention_time = obj$data$retention_time,
        data = as.matrix(obj$data$data_df),
      )
      return(sample)
    }
  } else {
    rlang::abort(glue("Support for reading {filename} not yet implemented"))
  }
}


#' @describeIn GCIMSDataset-class Constructor method
#'
#' @param pData A data frame with at least the `SampleID` and `FileName` columns.
#' @param base_dir A directory containing the file names described in `pData`
#' @param scratch_dir A directory to save intermediate results.
#' @param keep_intermediate A logical. Whether to keep or not intermediate files when realizing the GCIMSDataset object
#' @export
#' @return A GCIMSDataset object
#'
#' @examples
#' # Create a new GCIMSDataset with the convenient constructor function:
#' dummy_obj <-GCIMSDataset(
#'   pData = data.frame(SampleID = character(), filename = character(0)),
#'   base_dir = tempdir()
#' )
GCIMSDataset <- function(pData, base_dir, scratch_dir = tempfile("GCIMSDataset_tempdir_"), keep_intermediate = FALSE) {
  methods::new("GCIMSDataset", pData, base_dir, scratch_dir, keep_intermediate)
}


