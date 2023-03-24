#' GCIMSDataset
#'
#' @description GCIMSDataset is an R6 class to store a dataset.
#'
#' When the dataset is created, the `on_ram` option controls whether the actual
#' data is stored not in memory or it is read/saved from/to files as
#' needed, so the dataset object scales with large number of samples.
#'
#'
#' @importClassesFrom S4Vectors DataFrame
#' @export GCIMSDataset
#' @exportClass GCIMSDataset
GCIMSDataset <- R6::R6Class("GCIMSDataset",
  public = list(
    # FIXME Make some private, accessors, validation...
    # Fields:
    #' @field pData A data frame with at least the SampleID and filename columns.
    pData = NULL,
    #' @field scratch_dir A directory to save intermediate results.
    scratch_dir = NA_character_, # string
    #' @field delayed_ops Delayed operations
    delayed_ops = list(), # FIXME: This should be private
    #' @field previous_ops History of previous operations
    previous_ops = list(), # FIXME: This should be private
    #' @field hasheddir The basename of the current hashed directory in scracth_dir
    hasheddir = NA_character_, # string
    #' @field keep_intermediate Whether to keep intermediate files
    keep_intermediate = FALSE, # logical
    #' @field on_ram Whether samples should be kept on RAM or stored on disk.
    on_ram = FALSE, # logical
    #' @field samples List of samples, if on_ram or NULL otherwise
    samples = NULL, # list or NULL
    #' @field can_realize Whether we can realize the object, internal
    can_realize = FALSE, # logical
    #' @field align To store alignment results
    align = NULL, # list or NULL
    #' @field peaks To store the peak list
    peaks = NULL, # data frame or NULL
    #' @field TIS A matrix of n_samples vs drift time, with the Total Ion Spectrum of each sample
    TIS = NULL, # matrix or NULL
    #' @field RIC A matrix of n_samples vs retention time, with the Reverse Ion Chromatogram of each sample
    RIC = NULL, # matrix or NULL
    #' @field dt_ref A numeric drift time of reference
    dt_ref = NULL, # numeric or NULL
    #' @field rt_ref A numeric retention time of reference
    rt_ref = NULL, # numeric or NULL
    #' @field userData A list to store arbitrary data in the dataset
    userData = list(), # list
    # Methods:
    #' @description
    #' Create a new GCIMSDataset object
    #' @param pData A data frame holding phenotype data for the samples (or `NULL`). The data frame
    #'  should at least have a `SampleID` column, and a `filename` column if samples are stored in files.
    #' @param ... Unused
    #' @param samples A named list of `GCIMSSample` objects to be included in the dataset (or `NULL`). Names
    #' should correspond to the `SampleID` column in the `pData` data frame.
    #' @param base_dir The base directory. Sample `i` is found on `file.path(base_dir, pData$filename[i])`.
    #' @param parser Function that takes a file path and returns a [GCIMSSample] object. Use `"default"` to use the
    #' default parser in the GCIMS package, that supports `.mea` files (from GAS). Check
    #' out `vignette("importing-custom-data-formats", package = "GCIMS")` for more information
    #' @param scratch_dir A directory where intermediate and processed samples will be stored
    #' @param keep_intermediate If `TRUE`, intermediate results will not be deleted (ignored if `on_ram` is `TRUE`).
    #' @param on_ram If `TRUE`, samples are not stored on disk, but rather kept on RAM. Set it to `TRUE` only with
    #' small datasets.
    #' @examples
    #' dummy_dataset <- GCIMSDataset$new(
    #'   pData = data.frame(SampleID = character(), filename = character(0)),
    #'   base_dir = tempdir()
    #' )
    initialize = function(pData = NULL,
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

        self$pData <- pData
        self$scratch_dir <- scratch_dir
        self$dt_ref <- NULL
        self$rt_ref <- NULL
        self$TIS <- NULL
        self$RIC <- NULL
        self$delayed_ops <- list()
        self$previous_ops <- list()
        self$hasheddir <- ""
        self$keep_intermediate <- keep_intermediate
        self$on_ram <- on_ram
        self$samples <- NULL # Only used if on_ram is TRUE
        self$align <- NULL
        self$peaks <- NULL
        self$can_realize <- TRUE
        self$userData <- list()

        # How samples will be loaded:
        if (is.null(samples)) {
          # Check base_dir and files:
          base_dir <- validate_base_dir(base_dir)
          check_files(pData$FileName, base_dir)
          # From disk

          self$appendDelayedOp(
            operation = GCIMSDelayedOp(
              name = "read_sample",
              fun = read_sample,
              params = list(
                base_dir = base_dir,
                parser = parser
              )
            )
          )
          # Some sample stats:
          self$extract_dtime_rtime()
          self$extract_RIC_and_TIS()
        } else {
          # From list, to RAM or to disk?
          if (on_ram) {
            self$samples <- unname(samples)
          } else {
            # Dump samples to disk and prepare...
            self$hasheddir <- "from_list"
            save_to <- self$currentHashedDir()
            dir.create(save_to, showWarnings = FALSE, recursive = TRUE)
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
        self
      },
    #' @description prints the dataset to the screen
    print = function() {
      # FIXME: Refactor
      show(self)
    },
    #' @description Executes any pending action
    #' @param keep_intermediate Keep intermediate results. If `NA`, the setting from the dataset is used
    realize = function(keep_intermediate = NA) {
      if (!self$hasDelayedOps()) {
        return(self)
      }

      if (!canRealize(self)) {
        cli_abort(
          message = c(
            "UnexpectedError",
            paste0(
              "Tried to realize a dataset that could not be realized. This is ",
              "a bug in the package. Please report it to https://github.com/sipss/GCIMS/issues"
            )
          )
        )
      }
      canRealize(self) <- FALSE
      on.exit({canRealize(self) <- TRUE})

      self <- optimize_delayed_operations(self)
      if (isTRUE(self$on_ram)) {
        self <- realize_ram(self)
      } else {
        self <- realize_disk(self, keep_intermediate = keep_intermediate)
      }
      invisible(self)
    },
    #' @description
    #' Appends a delayed operation to the dataset so it will run afterwards
    #' @param operation A [GCIMSDelayedOp] object
    #' @return The modified GCIMSDataset object
    appendDelayedOp = function(operation) {
      self$delayed_ops <- c(self$delayed_ops, list(operation))
      self
    },
    #' @description
    #' Find out if the dataset has pending operations
    #' @return Returns `TRUE` if the dataset has pending operations, `FALSE` otherwise
    hasDelayedOps = function() {
      length(self$delayed_ops) > 0
    },
    # FIXME: Private?
    #' @description
    #' The directory where processed samples are saved
    #' @return A path to the directory where processed samples are kept, or `NULL`
    currentHashedDir = function() {
      hasheddir <- self$hasheddir
      if (hasheddir == "") {
        return(NULL)
      }
      file.path(
        self$scratch_dir,
        hasheddir
      )
    },
    # FIXME: Private
    #' @description
    #' The directory where samples will be saved after processing
    #' @return A path
    nextHashedDir = function() {
      if (!self$hasDelayedOps()) {
        return(self$currentHashedDir())
      }

      # A hash chain:
      current_hash <- self$hasheddir
      if (current_hash == "") {
        # The first hash comes from the SampleID and FileNames
        pd <- self$pData
        if (!all(c("SampleID", "FileName") %in% colnames(pd))) {
          cli_abort(c("Unexpected Error", "x" = "Expected pData with SampleID and FileName columns"))
        }
        current_hash <- digest::digest(
          list(
            pd[["SampleID"]],
            pd[["FileName"]]
          )
        )
      }
      for (op in purrr::keep(self$delayed_ops, modifiesSample)) {
        current_hash <- paste0(name(op), "_", digest::digest(list(current_hash, hashableDelayedOp(op))))
      }
      file.path(
        self$scratch_dir,
        current_hash
      )
    },
    #' @description
    #' Get a sample from a GCIMSDataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @return The GCIMSSample object
    getSample = function(sample) {
      self$realize()
      sample_id_num <- sample_name_or_number_to_both(sample, sampleNames(self))
      if (self$on_ram) {
        gcimssample <- self$samples[[sample_id_num$idx]]
      } else {
        filename <- paste0(sample_id_num$name, ".rds")
        current_intermediate_dir <- self$currentHashedDir()
        sample_file <- file.path(current_intermediate_dir, filename)
        if (!file.exists(sample_file)) {
          cli_abort("File not found: {sample_file} should have been created")
        }
        gcimssample <- readRDS(sample_file)
      }
      if (!methods::is(gcimssample, "GCIMSSample")) {
        cli_abort("Expected a GCIMSSample object, but it was not found")
      }
      updateObject(gcimssample)
    },
    #' @description
    #' Sets an action to extract the reference retention and drift times
    extract_dtime_rtime = function() {
      delayed_op <- GCIMSDelayedOp(
        name = "extract_dtime_rtime",
        fun_extract = .extract_dtime_rtime_fun_extract,
        fun_aggregate = .extract_dtime_rtime_fun_aggregate
      )
      self$appendDelayedOp(delayed_op)
      self
    },
    #' @description
    #' Get the Reverse Ion Chromatogram
    #' @return A matrix with the reverse ion chromatograms for all samples
    getRIC = function() {
      if (self$hasDelayedOps() || is.null(self$RIC)) {
        self$extract_RIC_and_TIS()
        self$realize()
      }
      out <- self$RIC
      dimnames(out) <- list(
        SampleID = self$sampleNames,
        retention_time_s = self$rt_ref
      )
      out
    },
    #' @description
    #' Extracts the RIC and the TIS
    #' @return The GCIMSDataset
    extract_RIC_and_TIS = function() {
      self$extract_dtime_rtime()
      delayed_op <- GCIMSDelayedOp(
        name = "extract_RIC_and_TIS",
        fun = NULL,
        fun_extract = .extract_RIC_and_TIS_fun_extract,
        fun_aggregate = .extract_RIC_and_TIS_fun_aggregate
      )
      self$appendDelayedOp(delayed_op)
      invisible(self)
    }
  ),
  active = list(
    #' @field sampleNames The sample names of the GCIMSDataset samples
    sampleNames = function(value) {
      # Getter
      if (missing(value)) return(self$pData$SampleID)
      # Setter
      if (nrow(self$pData) != length(value)) {
        cli_abort(
          c(
            "Invalid sample names",
            "x" = "The number of sample names given ({length(value)}) != Number of samples ({nrow(self$pData)})"
          )
        )
      }
      if (anyNA(value) || anyDuplicated(value)) {
        cli_abort("Sample names must be unique and not missing")
      }
      private$rename_intermediate_files(value)
      self$pData$SampleID <- value
      self
    }
  ),
  private = list(
    rename_intermediate_files = function(new_names) {
      current_dir <- self$currentHashedDir()
      if (is.null(current_dir)) {
        return()
      }
      current_dir_old <- paste0(current_dir, "_old")
      if (!dir.exists(current_dir)) {
        return()
      }
      file.rename(current_dir, current_dir_old)
      dir.create(current_dir, showWarnings = FALSE, recursive = TRUE)
      old_names <- names(self)
      for (i in seq_along(new_names)) {
        old_name <- paste0(old_names[i], ".rds")
        new_name <- paste0(new_names[i], ".rds")
        if (file.exists(file.path(current_dir_old, old_name))) {
          file.rename(
            file.path(current_dir_old, old_name),
            file.path(current_dir, new_name)
          )
        }
      }
      unlink(current_dir_old, recursive = TRUE)
      invisible(NULL)
    }
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
    cli_abort(
      message = c(
        "Files not found",
        "x" = "The following files {and_n_more} were not found",
        "i" = "{utils::head(missing_files, to_show)}"
      )
    )
  }

}
# So we can use S4 methods:
setOldClass("GCIMSDataset")

abort_if_errors <- function(errors, title = "Errors found") {
  if (length(errors) > 0) {
    names(errors) <- rep("x", length(errors))
    cli_abort(message = c(title, errors))
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
    cli_warn(
      c("samples should be a named list of GCIMSSample objects.",
        "i" = "Sample names should correspond to SampleIDs",
        "!" = "Using Sample1...Sample{length(samples)} as SampleID by default"
      )
    )
    names(samples) <- paste0("Sample", seq_len(length(samples)))
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



CurrentHashedDir <- function(object) {
  cli_warn("CurrentHashedDir(obj) is deprecated. Use obj$currentHashedDir() instead")
  object$currentHashedDir()
}

"CurrentHashedDir<-" <- function(object, value) {
  cli_warn("{.code CurrentHashedDir(obj) <- value} is deprecated. Use {.code obj$hasheddir <- basename(value)} instead")
  object$hasheddir <- basename(value)
  object
}

NextHashedDir <- function(object) {
  cli_warn("Deprecated: NextHashedDir(obj) is deprecated. Use obj$nextHashedDir() instead.")
  object$nextHashedDir()
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
      cli_abort("Support for reading {filename} not yet implemented")
    }
  } else {
    sample <- parser(filename)
  }
  if (!methods::is(sample, "GCIMSSample")) {
    cli_abort("R Object in {filename} is not of type GCIMSSample")
  }
  sample
}


#' @name GCIMSDataset-old-constructor
#' @title GCIMSDataset old constructor
#'
#' @param pData A data frame with at least the `SampleID` and `FileName` columns.
#' @param base_dir A directory containing the file names described in `pData`
#' @param parser Either a string `"default"` or a function that takes a filename and returns a [GCIMSSample-class] object
#' @param scratch_dir A directory to save intermediate results.
#' @param keep_intermediate A logical. Whether to keep or not intermediate files when realizing the GCIMSDataset object
#' @param on_ram A logical. If `TRUE`, samples are kept on RAM. This is faster
#' as long as you have enough memory to keep all samples. If this is `TRUE`, then
#' `scratch_dir` and `keep_intermediate` are ignored.
#' @return A GCIMSDataset object
#'
#' @examples
#' # Create a new GCIMSDataset with the convenient constructor function:
#' dummy_obj <-GCIMSDataset$new(
#'   pData = data.frame(SampleID = character(), filename = character(0)),
#'   base_dir = tempdir()
#' )
NULL

#' GCIMSDataset_fromList
#'
#' @param samples A named list of [GCIMSSample] objects. names should match `pData$SampleID`
#' @param pData A data frame with at least the SampleID and filename columns.
#' @param scratch_dir A directory to save intermediate results.
#' @param keep_intermediate Whether to keep sample files for intermediate results. Only used if `on_ram=FALSE`
#' @param on_ram logical. Whether the dataset should be kept stored on RAM or on disk.
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
#' @export
GCIMSDataset_fromList <- function(samples, pData=NULL,
                                  scratch_dir = tempfile("GCIMSDataset_tempdir_"),
                                  keep_intermediate = FALSE, on_ram = TRUE) {
  GCIMSDataset$new(
    pData = pData,
    samples = samples,
    scratch_dir = scratch_dir,
    keep_intermediate = keep_intermediate,
    on_ram = on_ram
  )
}

