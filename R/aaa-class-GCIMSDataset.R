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
    abort(
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
    abort(message = c(title, errors))
  }
}

validate_pData <- function(pData) {
  # Convert to S4Vectors::DataFrame so it fits in a slot
  # FIXME: We can review and possibly avoid this conversion:
  pData <- S4Vectors::DataFrame(pData)
  # Find columns where SampleID and FileName are, tolerating other case conventions
  pheno_names <- tolower(colnames(pData))
  sampleid_col <- which(pheno_names == "sampleid")
  colnames(pData)[sampleid_col] <- "SampleID"

  # A place to store errors for reporting afterwards:
  errors <- character(0L)

  # One and only one SampleID is always mandatory
  if (length(sampleid_col) == 0) {
    errors <- c(errors, "pData should have a SampleID column with a sample name")
  }
  if (length(sampleid_col) > 1) {
    errors <- c(errors, "pData should only have one SampleID column (don't use sampleid or other casing)")
  }

  # SampleID should be a character column. If it is integer we warn.
  sampleids <- pData[[sampleid_col]]
  if (is.factor(sampleids)) {
    sampleids <- as.character(sampleids)
  }
  if (is.integer(sampleids)) {
    cli_warn("pData$SampleID is coerced to a character column (integers were found instead).")
    sampleids <- as.character(sampleids)
  }
  pData[[sampleid_col]] <- sampleids

  # If we are creating a dataset reading from files, then samples==NULL and we need a FileName column
  filename_col <- which(pheno_names == "filename")
  if (length(filename_col) > 1) {
    errors <- c(errors, "pData should only have one FileName column (don't use filename or other casing)")
  }
  if (length(filename_col) > 0) {
    colnames(pData)[filename_col] <- "FileName"
  }
  if (length(filename_col) == 0) {
    pData[["FileName"]] <- NA_character_
  }

  abort_if_errors(errors, title = "pData is not valid")

  # Return corrected pData, with SampleID on the first column and FileName on the second column
  pData <- pData[, c("SampleID", "FileName", setdiff(colnames(pData), c("SampleID", "FileName")))]

  pData
}

validate_pData_samples <- function(pData, samples) {
  samples <- validate_samples(samples)
  if (!is.null(samples) && is.null(pData)) {
    pData <- data.frame(SampleID = names(samples))
  }
  if (is.null(samples) && is.null(pData)) {
      pData <- data.frame(
        SampleID = character(0L),
        FileName = character(0L)
      )
      samples <- structure(list(), names = character(0))
      return(list(pData = pData, samples = samples))
  }

  pData <- validate_pData(pData)
  # A place to store errors for reporting afterwards:
  errors <- character(0L)


  if (is.null(samples)) {
    # Must be loaded from disk
    if (!"FileName" %in% colnames(pData)) {
      errors <- c(errors, "pData should have a FileName column with the paths to the sample files")
    }
  } else {
    # Samples are on RAM.
    if (any(!is.na(pData[["FileName"]]))) {
      cli_inform("GCIMSDataset: pData$FileName will be ignored since samples are already given")
    }
    # ensure no sample is missing
    missing_samples <- which(!pData[["SampleID"]] %in% names(samples))
    missing_annotations <- which(!names(samples) %in% pData[["SampleID"]])
    if (length(missing_annotations) > 0) {
      cli_warn(
        c(
          "{length(missing_annotations)} samples are not present in the pData and will not be included in the GCIMSDataset",
          "i" = "For instance samples {head(names(samples)[missing_annotations])}"
        )
      )
    }
    if (length(missing_samples) > 0) {
      errors <- c(
        errors,
        "Samples given as a list, but are missing in pData.",
        "x" = "Samples for {length(missing_samples)} SampleIDs were not found",
        "i" = "Either remove them from pData or provide the samples",
        "i" = "For instance: {head(sampleids[missing_samples])}"
      )
    }
  }

  abort_if_errors(errors, title = "pData is not valid")

  if (!is.null(samples)) {
    # pData is valid, samples is valid, but we ensure they match:
    samples <- samples[pData$SampleID]
  }
  list(
    pData = pData,
    samples = samples
  )
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

validate_scratch_dir <- function(scratch_dir, on_ram) {
  errors <- character(0L)
  if (on_ram) {
    if (is.null(scratch_dir) ||
        identical(scratch_dir, NA) ||
        identical(scratch_dir, NA_character_)) {
      return(NA_character_)
    }
    if (length(scratch_dir) != 1 || !is.character(scratch_dir)) {
      errors <- c(errors, "scratch_dir must be a string of length 1 or NA")
    }
    abort_if_errors(errors, title = "scratch_dir is not valid")
    return(scratch_dir)
  }
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
  errors <- c()
  if (!rlang::is_bool(keep_intermediate)) {
    errors <- c(errors, "keep_intermediate must be either TRUE or FALSE")
  }
  abort_if_errors(errors, title = "keep_intermediate is not valid")
  keep_intermediate
}

validate_on_ram <- function(on_ram) {
  errors <- c()
  if (length(on_ram) != 1 || is.na(on_ram)) {
    errors <- c(errors, "on_ram must be either TRUE or FALSE")
  }
  on_ram <- as.logical(on_ram)
  abort_if_errors(errors, title = "on_ram is not valid")
  on_ram
}

validate_samples <- function(samples) {
  if (is.null(samples)) {
    # Expected to read them from disk
    return(samples)
  }
  if (!is.list(samples)) {
    cli_abort(
      c(
        "samples should be either NULL or a named list of GCIMSSample objects.",
        "i" = "Found {class(samples)} instead"
      )
    )
  }

  if (is.null(names(samples))) {
    cli_abort("samples should be a named list of GCIMSSample objects. Sample names should correspond to SampleIDs")
  }
  empty_names_idx <- which(names(samples) == "")
  if (length(empty_names_idx) > 0) {
    cli_abort("{length(empty_names_idx)} samples did not have names. Please name them (e.g. {head(empty_names_idx)})")
  }
  # check for duplicated names:
  repeated_names <- unique(names(samples)[duplicated(names(samples))])
  if (length(repeated_names) > 0) {
    cli_abort("sample names should be unique. Samples {repeated_names} appeared more than once")
  }
  are_gcimssample <- purrr::map_lgl(samples, function(x) methods::is(x, "GCIMSSample"))
  not_gcims_samples <- which(!are_gcimssample)
  if (length(not_gcims_samples) > 0) {
    cli_abort("{length(not_gcims_samples)} samples were not of type `GCIMSSample`. (e.g. at {head(not_gcims_samples)}))")
  }
  samples
}

validate_parser <- function(parser) {
  if (!identical(parser, "default") && !is.function(parser)) {
    cli_abort("Parser must be either 'default' or a function")
  }
  parser
}


methods::setMethod(
  "initialize", "GCIMSDataset",
  function(
    .Object,
    pData = NULL,
    ...,
    samples = NULL,
    base_dir = NULL,
    parser = "default",
    scratch_dir = tempfile("GCIMSDataset_tempdir_"),
    keep_intermediate = FALSE,
    on_ram = FALSE
  ) {
    # The GCIMSDataset object always uses the pData dataframe to store phenotype data,
    # (often also called metadata or sample annotations).
    #
    # Besides, there are several possibilities for initializing the GCIMSDataset object:
    #  (a) Initialize reading samples from "disk", at file.path(base_dir, pData$FileName)
    #    (a.1) If on_ram, the samples slot will be populated with GCIMSSample objects. This is suitable for small datasets
    #    (a.2) If not on_ram, the samples slot will not be populated, scratch_dir will be used instead to store sample objects
    #
    #  (b) Initialize using a named "list" of GCIMSSample objects, with names according to pData$SampleID
    #    (b.1) If on_ram, the samples slot will be populated with that list
    #    (b.2) If not on_ram, the samples will be dumped to scratch_dir
    #
    #
    pds <- validate_pData_samples(pData, samples)
    pData <- pds$pData
    samples <- pds$samples
    on_ram <- validate_on_ram(on_ram)
    scratch_dir <- validate_scratch_dir(scratch_dir, on_ram)
    keep_intermediate <- validate_keep_intermediate(keep_intermediate)
    parser <- validate_parser(parser)

    # We want the GCIMSDataset object to be mutable, so any pending delayed operation
    # can be applied in-place if needed.
    # Therefore, instead of using immutable slots for our attributes, we will
    # use an environment where we will place.
    .Object@envir <- rlang::new_environment()
    .Object@envir$pData <- pData
    .Object@envir$scratch_dir <- scratch_dir
    .Object@envir$dt_ref <- NULL
    .Object@envir$rt_ref <- NULL
    .Object@envir$TIS <- NULL
    .Object@envir$RIC <- NULL
    .Object@envir$delayed_ops <- NULL
    .Object@envir$previous_ops <- list()
    .Object@envir$hasheddir <- ""
    .Object@envir$keep_intermediate <- keep_intermediate
    .Object@envir$on_ram <- on_ram
    .Object@envir$samples <- NULL # Only used if on_ram is TRUE
    canRealize(.Object) <- TRUE
    # How samples will be loaded:
    if (is.null(samples)) {
      # Check base_dir and files:
      base_dir <- validate_base_dir(base_dir)
      check_files(pData$FileName, base_dir)
      # From disk
      .Object <- appendDelayedOp(
        .Object,
        GCIMSDelayedOp(
          name = "read_sample",
          fun = read_sample,
          params = list(
            base_dir = base_dir,
            parser = parser
          )
        )
      )
      # Some sample stats:
      .Object <- extract_dtime_rtime(.Object)
      .Object <- extract_RIC_and_TIS(.Object)
    } else {
      # From list, to RAM or to disk?
      if (on_ram) {
        .Object@envir$samples <- unname(samples)
      } else {
        # Dump samples to disk and prepare...
        CurrentHashedDir(.Object) <- "from_list"
        save_to <- CurrentHashedDir(.Object)
        purrr::walk2(
          pData$SampleID,
          samples,
          function(sample_name, sample, save_to) {
            next_filename <- file.path(save_to, paste0(sample_name, ".rds"))
            saveRDS(sample, next_filename)
          },
          save_to = save_to,
          .progress = "Saving list of samples to disk..."
        )
      }
    }
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
    # The first hash comes from the sampleID and FileNames
    pd <- pData(object)
    if (!all(c("SampleID", "FileName") %in% colnames(pd))) {
      abort(c("Unexpected Error", "x" = "Expected pData with SampleID and FileName columns"))
    }
    current_hash <- digest::digest(
      list(
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


read_sample <- function(filename, base_dir, parser = "default") {
  filename <- file.path(base_dir, filename)
  filename_l <- tolower(filename)
  sample <- NULL
  if (identical(parser, "default")) {
    if (endsWith(filename_l, ".mea.gz") || endsWith(filename_l, ".mea")) {
      sample <- read_mea(filename)
    }
    else if (endsWith(filename_l, ".rds")) {
      sample <- readRDS(filename)
    } else {
      cli::cli_abort("Support for reading {filename} not yet implemented")
    }
  } else {
    sample <- parser(filename)
  }
  if (!methods::is(sample, "GCIMSSample")) {
    cli::cli_abort("R Object in {filename} is not of type GCIMSSample")
  }
  sample
}


#' @describeIn GCIMSDataset-class Constructor method
#'
#' @param pData A data frame with at least the `SampleID` and `FileName` columns.
#' @param base_dir A directory containing the file names described in `pData`
#' @param parser Either a string `"default"` or a function that takes a filename and returns a [GCIMSSample-class] object
#' @param scratch_dir A directory to save intermediate results.
#' @param keep_intermediate A logical. Whether to keep or not intermediate files when realizing the GCIMSDataset object
#' @param on_ram A logical. If `TRUE`, samples are kept on RAM. This is faster
#' as long as you have enough memory to keep all samples. If this is `TRUE`, then
#' `scratch_dir` and `keep_intermediate` are ignored.
#' @export
#' @return A GCIMSDataset object
#'
#' @examples
#' # Create a new GCIMSDataset with the convenient constructor function:
#' dummy_obj <-GCIMSDataset(
#'   pData = data.frame(SampleID = character(), filename = character(0)),
#'   base_dir = tempdir()
#' )
GCIMSDataset <- function(pData, base_dir, parser = "default", scratch_dir = tempfile("GCIMSDataset_tempdir_"), keep_intermediate = FALSE, on_ram = FALSE) {
  methods::new(
    "GCIMSDataset", pData = pData, base_dir = base_dir, parser = parser,
    scratch_dir = scratch_dir, keep_intermediate = keep_intermediate, on_ram = on_ram
  )
}

#' @describeIn GCIMSDataset-class Constructor method
#'
#' @inheritParams GCIMSDataset
#' @param samples A named list of [GCIMSSample] objects. names should match `pData$SampleID`
#' @export
#' @return A GCIMSDataset object
#'
#' @examples
#' # Create a new GCIMSDataset with the convenient constructor function:
#' sample1 <- GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' dummy_obj <- GCIMSDataset_fromList(
#'   pData = data.frame(SampleID = "Sample1", Sex = "female"),
#'   samples = list(Sample1 = sample1)
#' )
GCIMSDataset_fromList <- function(samples, pData=NULL, scratch_dir = tempfile("GCIMSDataset_tempdir_"), keep_intermediate = FALSE, on_ram = TRUE) {
  methods::new(
    "GCIMSDataset", pData = pData, samples = samples,
    scratch_dir = scratch_dir, keep_intermediate = keep_intermediate, on_ram = on_ram
  )
}

