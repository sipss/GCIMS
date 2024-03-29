#' GCIMSDataset
#' @aliases GCIMSDataset-class
#'
#' @description GCIMSDataset is an R6 class to store a dataset.
#'
#' When the dataset is created, the `on_ram` option controls whether the actual
#' data is stored not in memory or it is read/saved from/to files as
#' needed, so the dataset object scales with large number of samples.
#'
#' # Constructors:
#'
#' - [`GCIMSDataset$new()`](#method-GCIMSDataset-new)
#' - [`GCIMSDataset$new_from_list()`](#constructor-GCIMSDataset-new-from-list): Create a new GCIMSDataset from a list of samples
#' - [`GCIMSDataset$new_from_saved_dir()`](#constructor-GCIMSDataset-new-from-saved-dir): Create a new on disk GCIMSDataset from a directory
#'
#'
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="constructor-GCIMSDataset-new-from-list"></a>}}
#' \if{latex}{\out{\hypertarget{constructor-GCIMSDataset-new-from-list}{}}}
#' ## Constructor `new_from_list()`
#'
#' Create a new GCIMSDataset object from a list of samples. Note that with this
#' constructor `on_ram` is `TRUE` by default
#'
#' ### Usage
#'
#' ```
#' GCIMSDataset$new_from_list(
#'   samples,
#'   pData=NULL,
#'   scratch_dir = NULL,
#'   keep_intermediate = FALSE,
#'   on_ram = TRUE
#' )
#' ```
#'
#' ### Arguments
#'
#' See [`GCIMSDataset$new()`](#method-GCIMSDataset-new)
#'
#' ## Constructor `new_from_saved_dir()`
#'
#' Creates a new GCIMSDataset object from a directory where a GCIMSDataset with
#' `on_ram=FALSE` was saved.
#'
#' ### Usage
#'
#' ```
#' GCIMSDataset$new_from_saved_dir(
#'   input_dir,
#'   scratch_dir = dirname(input_dir)
#' )
#' ```
#'
#' ### Arguments
#'
#' - `input_dir`: The path to the directory where the `dataset.rds` is saved and all the corresponding `sample_*.rds` files are.
#'    Typically a subdirectory of `scratch_dir`.
#' - `scratch_dir`: The new scratch directory where further processing samples will be saved. By default it is the parent of `input_dir`.
#'
#' @export GCIMSDataset
#' @exportClass GCIMSDataset
GCIMSDataset <- R6::R6Class("GCIMSDataset",
  public = list(
    # Fields
    # IMPORTANT: If you add a field, make sure to update the subset method to ensure consistency

    #' @field pData A data frame with at least the SampleID and filename columns.
    pData = NULL,
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
                          base_dir = NULL,
                          ...,
                          samples = NULL,
                          parser = "default",
                          scratch_dir = NULL,
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
        keep_intermediate <- validate_keep_intermediate(keep_intermediate)
        parser <- validate_parser(parser)

        self$pData <- pData
        if (is.null(samples)) {
          # Check base_dir and files:
          base_dir <- validate_base_dir(base_dir)
          check_files(pData$FileName, base_dir)
          samples <- pData$FileName
          names(samples) <- pData$SampleID
        }
        if (on_ram) {
          private$delayed_dataset <- DelayedDatasetRAM$new(
            samples = samples,
            sample_class = "GCIMSSample"
          )
        } else {
          scratch_dir <- validate_scratch_dir(scratch_dir)
          private$delayed_dataset <- DelayedDatasetDisk$new(
            samples = samples,
            scratch_dir = scratch_dir,
            keep_intermediate = keep_intermediate,
            sample_class = "GCIMSSample"
          )
        }
        if (is.character(samples)) {
          # So this is going to be our first action: read the samples
          self$appendDelayedOp(
            operation = DelayedOperation(
              name = "read_sample",
              fun = read_sample,
              params = list(
                base_dir = base_dir,
                parser = parser
              )
            )
          )
        }
        private$updateSampleDescriptions(names(samples))
        # Optimize dtime and RIC and TIS extraction on every realize action:
        private$delayed_dataset$registerOptimization(optimize_RIC_TIS)
        # Extract RIC and TIS when we first realize:
        self$extract_dtime_rtime()
        self$extract_RIC_and_TIS()
        self
      },
    #' @description prints the dataset to the screen
    print = function() {
      outstring <- yaml::as.yaml(private$describe_as_list())
      cat(outstring)
    },
    #' @description
    #' Create a new dataset containing a subset of the samples
    #' @param samples A numeric vector (sample indices), a character vector (sample names)
    #' or a logical vector of the length equal to the number of samples in the
    #' dataset (`TRUE` elements will be subset)
    #' @param inplace if `TRUE` subset happens in-place, otherwise subset will return a copy.
    #' @param new_scratch_dir A new scratch directory, only used if `inplace=FALSE` and the dataset is on-disk.
    #' @return A GCIMSDataset (new or the current one depending on `inplace`), with the requested sample subset
    subset = function(samples, inplace = FALSE, new_scratch_dir = NA) {
      self$realize()

      if (inplace) {
        new <- self
      } else {
        # FIXME: Efficiency. Copying all samples and then removing some is not smart
        new <- self$copy(scratch_dir = new_scratch_dir)
      }
      new$.impl__subset__(samples)
      new
    },
    #' @description
    #' Do not call this method. It does an inplace subset. Use
    #' `obj$subset(samples, inplace = TRUE)` instead
    #' @param samples A numeric vector (sample indices), a character vector (sample names)
    #' or a logical vector of the length equal to the number of samples in the
    #' dataset (`TRUE` elements will be subset)
    #' @return The given GCIMSDataset object, with a subset of the samples
    .impl__subset__ = function(samples) {
      # FIXME: .impl__subset__ could become a private method if R6 allowed for a bit of more magic...
      sample_info <- sample_name_or_number_to_both(samples, self$sampleNames)
      if (!is.null(self$pData)) {
        self$pData <- self$pData[sample_info$logical, , drop = FALSE]
      }
      self$align <- NULL
      if (!is.null(self$peaks)) {
        self$peaks <- self$peaks[self$peaks$SampleID %in% sample_info$name, , drop = FALSE]
      }

      self$TIS <- NULL
      self$RIC <- NULL
      self$dt_ref <- NULL
      self$rt_ref <- NULL
      private$delayed_dataset$subset(samples)
      self
    },
    #'
    #' @description
    #' Appends a delayed operation to the dataset so it will run afterwards
    #' @param operation A [DelayedOperation-class] object
    #' @return The modified GCIMSDataset object
    appendDelayedOp = function(operation) {
      private$delayed_dataset$appendDelayedOp(operation)
      self
    },
    #' @description
    #' Find out if the dataset has pending operations
    #' @return Returns `TRUE` if the dataset has pending operations, `FALSE` otherwise
    hasDelayedOps = function() {
      private$delayed_dataset$hasDelayedOps()
    },
    #' @description
    #' Execute all pending operations on the dataset
    #' @param keep_intermediate logical or `NA`. Only when the analysis is on disk, keep intermediate result files.
    #' If `NA`, the `keep_intermediate` option given at the dataset initialization takes precedence.
    #' @return The dataset object, invisibly
    realize = function(keep_intermediate = NA) {
      private$delayed_dataset$realize(keep_intermediate = keep_intermediate, dataset = self)
      invisible(self)
    },
    #' @description
    #' Get a sample from a GCIMSDataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @return The GCIMSSample object
    getSample = function(sample) {
      private$delayed_dataset$getSample(sample, dataset = self)
    },
    #' @description
    #' Sets an action to extract the reference retention and drift times
    extract_dtime_rtime = function() {
      delayed_op <- DelayedOperation(
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
      delayed_op <- DelayedOperation(
        name = "extract_RIC_and_TIS",
        fun = NULL,
        fun_extract = .extract_RIC_and_TIS_fun_extract,
        fun_aggregate = .extract_RIC_and_TIS_fun_aggregate
      )
      self$appendDelayedOp(delayed_op)
      invisible(self)
    },
    #' @description
    #' Whether the dataset is saved on disk or stored in RAM
    #' @return `TRUE` if on disk, `FALSE` otherwise
    is_on_disk = function() {
      inherits(private$delayed_dataset, "DelayedDatasetDisk")
    },
    #' @description
    #' Creates a copy of the dataset. If the dataset is stored on disk, then
    #' a new `scratch_dir` must be used.
    #' @param scratch_dir The scratch directory where samples being processed will be stored, if the copy is on disk.
    #' @return A new GCIMSDataset object
    copy = function(scratch_dir = NA) {
      # I should just override clone(), but that's not possible. See:
      # https://github.com/r-lib/R6/pull/273
      # Meanwhile I'll work with this copy() method.
      # No one should use clone() with a GCIMSDataset object, unless they are aware of what they are doing.
      # A shallow clone is useless in our case (we should enforce deep=TRUE) because we need to copy the DelayedDataset
      # Unfortunately the on-disk dataset can't be just cloned, because the scratch_dir must be changed.

      # On RAM datasets must only enforce deep=TRUE when cloning.
      # Besides that limitation, the clone method is fine.
      if (!self$is_on_disk()) {
        return(self$clone(deep = TRUE))
      }

      scratch_dir <- validate_scratch_dir(scratch_dir)
      # On disk datasets must use deep=TRUE AND update the scratch directory:
      # The scratch directory must be different:
      if (private$delayed_dataset$scratchDir == scratch_dir) {
        cli_abort(
          c(
            "Could not create a copy of the dataset",
            "x" = "The given scratch_dir {.path scratch_dir} is equal to current one and should have been different",
            "i" = "Please provide a new with {.code obj$copy(scratch_dir = )}"
          )
        )
      }

      new_dataset <- self$clone(deep = TRUE)
      # This could be a private method?
      new_dataset$updateScratchDir(scratch_dir)
      new_dataset
    },
    #' @description
    #' For on-disk datasets, copy all samples to a new scratch dir.
    #' This is useful when creating copies of the dataset, using the `dataset$copy()` method.
    #' @param scratch_dir The new `scratch_dir`, must be different from the current one
    #' @param override_current_dir Typically used only internally, overrides the location of the samples. Useful when we are
    #' loading a dataset from a directory and the directory was moved since it was saved.
    updateScratchDir = function(scratch_dir, override_current_dir = NULL) {
      if (!self$is_on_disk()) {
        cli_warn("obj$updateScratchDir() on a on-RAM dataset has no effect")
        return(invisible(NULL))
      }
      scratch_dir <- validate_scratch_dir(scratch_dir)
      # On disk datasets must use deep=TRUE AND update the scratch directory:
      # The scratch directory must be different so further analysis steps do not get mixed up,
      # unless we are loading a dataset saved on disk instead of copying a dataset
      if (is.null(override_current_dir) && private$delayed_dataset$scratchDir == scratch_dir) {
        cli_abort(
          c(
            "Could not create a copy of the dataset",
            "x" = "The given scratch_dir {.path scratch_dir} is equal to current one and should have been different",
            "i" = "Please provide a new with {.code obj$updateScratchDir(scratch_dir = )}"
          )
        )
      }
      private$delayed_dataset$updateScratchDir(
        new_scratch_dir = scratch_dir,
        dataset = self,
        override_current_dir = override_current_dir
      )
      invisible(NULL)
    },
    #' @description Get the directory where processed samples are being saved, on on-disk datasets.
    #' @return Either a path or `NULL`. `NULL` is returned if samples have not been saved (either
    #' because have not been loaded or because the dataset is stored on RAM)
    getCurrentDir = function() {
      if (!self$is_on_disk()) {
        cli_warn("{.code GCIMSDataset$getCurrentDir()} only makes sense for on-disk datasets")
        return(NULL)
      }
      private$delayed_dataset$getCurrentDir()
    }
  ),
  active = list(
    #' @field sampleNames The sample names of the GCIMSDataset samples
    sampleNames = function(value) {
      # Getter
      if (missing(value)) return(private$delayed_dataset$sampleNames)
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
      # Rename the delayed_dataset
      private$delayed_dataset$sampleNames <- value
      # Update pData accordingly
      self$pData$SampleID <- value
      # And the sample descriptions:
      self$updateSampleDescriptions(value)
      self
    }
  ),
  private = list(
    # @field The DelayedDataset object
    # IMPORTANT: If you add a field here, make sure the subset method
    # accounts for it
    delayed_dataset = NULL,
    updateSampleDescriptions = function(new_names) {
      op = DelayedOperation(
        name = "setSampleNamesAsDescription",
        fun = function(sample, sample_name) {
          description(sample) <- sample_name
          sample
        },
        params_iter = list(
          sample_name = setNames(new_names, new_names)
        )
      )
      self$appendDelayedOp(op)
    },
    describe_as_list = function() {
      out <- list()
      if (self$is_on_disk()) {
        disk_path <- self$getCurrentDir()
        if (is.null(disk_path)) {
          stored_info <- glue::glue("Stored on disk (not loaded yet)")
        } else {
          stored_info <- glue::glue("Stored on disk (at {disk_path})")
        }
      } else {
        stored_info <- "Stored on RAM"
      }
      root_txt <- "A GCIMSDataset"
      sample_info <- glue::glue("With {length(self$sampleNames)} samples")
      pheno_info <- phenos_to_string(self$pData)
      out[[root_txt]] <- list(
        sample_info,
        stored_info,
        pheno_info,
        private$delayed_dataset$history_as_list(),
        private$delayed_dataset$pending_as_list()
      )
      out
    }
  )
)





optimize_RIC_TIS <- function(delayed_ops) {
  # Extra operations that extract the rtime and dtime or TIS and RIC from the object can be delayed
  # Find them:
  where_extract_times <- purrr::map_lgl(delayed_ops, function(op) {name(op) == "extract_dtime_rtime"})
  where_extract_RIC_TIS <- purrr::map_lgl(delayed_ops, function(op) {name(op) == "extract_RIC_and_TIS"})
  where_extract <- where_extract_times | where_extract_RIC_TIS
  # Not found, return:
  if (!any(where_extract)) {
    return(delayed_ops)
  }
  # Get one action of each, if present:
  if (any(where_extract_times)) {
    extract_dtime_rtime <- delayed_ops[[which(where_extract_times)[1]]]
  } else {
    extract_dtime_rtime <- NULL
  }

  if (any(where_extract_RIC_TIS)) {
    extract_RIC_TIS <- delayed_ops[[which(where_extract_RIC_TIS)[1]]]
  } else {
    extract_RIC_TIS <- NULL
  }
  # Remove all those actions
  delayed_ops[where_extract] <- NULL
  # Return them appended at the end
  return(c(delayed_ops, extract_dtime_rtime, extract_RIC_TIS))
}


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

validate_scratch_dir <- function(scratch_dir) {
  errors <- character(0L)
  if (is.null(scratch_dir) || is.na(scratch_dir)) {
    inform_temp_dir <- !is.null(scratch_dir)
    scratch_dir <- tempfile("GCIMSDataset_tempdir_")
    dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
    if (inform_temp_dir) {
      cli_inform(
        c(
          "Creating a temporary directory for processing data",
          "i" = "The directory {.path {scratch_dir}} has been created to store the dataset and process it",
          "i" = "This directory will be deleted when the R session closes.",
          "i" = "If you want to save your results, please consider providing a {.code scratch_dir} when creating the dataset.",
          "i" = "You can omit this message if you set {.code scratch_dir = NULL}."
        )
      )
    }
  } else {
    if (!rlang::is_string(scratch_dir)) {
      errors <- c(errors, "scratch_dir should be a string")
    }
    dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
  }
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
      cli_abort(
        message = c(
          "Support for reading {.path filename} not yet implemented",
          "i" = "Please use a custom {.code parser}.",
          "i" = "Checkout the vignette at {.url https://sipss.github.io/GCIMS/articles/importing-custom-data-formats.html}"
        )
      )
    }
  } else {
    sample <- parser(filename)
  }
  sample
}

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
                                  scratch_dir = NULL,
                                  keep_intermediate = FALSE, on_ram = TRUE) {
  GCIMSDataset$new(
    pData = pData,
    samples = samples,
    scratch_dir = scratch_dir,
    keep_intermediate = keep_intermediate,
    on_ram = on_ram
  )
}

# FIXME: Fix consistency in method capitalization?

GCIMSDataset$new_from_list <- function(samples, pData=NULL,
                                       scratch_dir = NULL,
                                       keep_intermediate = FALSE, on_ram = TRUE) {
  GCIMSDataset$new(
    pData = pData,
    samples = samples,
    scratch_dir = scratch_dir,
    keep_intermediate = keep_intermediate,
    on_ram = on_ram
  )
}

GCIMSDataset$new_from_saved_dir <- function(input_dir, scratch_dir = dirname(input_dir)) {
  rds_file <- file.path(input_dir, "dataset.rds")
  if (!file.exists(rds_file)) {
    cli_abort(c(
      "Can't load GCIMSDataset",
      "x" = "File {.path {rds_file}} does not exist"
    ))
  }
  dataset <- readRDS(rds_file)
  if (!dataset$is_on_disk()) {
    return(dataset)
  }
  scratch_dir <- validate_scratch_dir(scratch_dir)
  dataset$updateScratchDir(
    scratch_dir,
    override_current_dir = input_dir
  )
  dataset
}
