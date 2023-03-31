#' DelayedDatasetDisk
#'
#' @description
#' A class that contains a dataset where its samples are backed on disk.
#' Each sample is stored on a file, and when queued operations are executed files
#' are loaded and saved on a new folder.
#'
#' @details
#'
#' This class is not exported, but if you want to use it reach us at
#' <https://github.com/sipss/GCIMS/issues/> and we will export it.
#'
#' @keywords internal
DelayedDatasetDisk <- R6::R6Class(
  "DelayedDatasetDisk",
  inherit = DelayedDatasetBase,
  public = list(
    #' @description
    #' Create a delayed dataset on disk
    #' @param samples A named vector. The names are sample ids, the values are
    #'  either filenames or sample objects. If they are not filenames then the objects are dumped to disk.
    #'  If they are filenames, the filenames are relative to `base_dir`.
    #' @param scratch_dir The directory where samples being processed will be saved
    #' @param keep_intermediate A logical value, whether intermediate realization steps should be saved.
    #' @param dataset The dataset R6 object where aggregated results from queued operations are stored
    #' @param dataset_class The class of the given dataset, used just to validate the contract between
    #' the delayed actions and the dataset. If `NULL` action return values are not checked
    #' @param sample_class The class of the samples in the dataset, used just to validate the contract between
    #' the delayed actions and the samples. If `NULL` action return values are not checked
    #' @return The DelayedDatasetDisk object
    initialize = function(samples, scratch_dir, keep_intermediate = FALSE, dataset = NULL, dataset_class = NULL, sample_class = NULL) {
      super$initialize(dataset = dataset, dataset_class = dataset_class, sample_class = sample_class)
      private$scratch_dir <- scratch_dir
      private$keep_intermediate <- keep_intermediate
      if (is.character(samples)) {
        # Samples are filenames, that we will read
        private$samples <- samples
      } else {
        # Samples are objects, we dump them to disk:
        # Dump samples to disk and prepare...
        private$hasheddir <- "from_list"
        save_to <- private$currentHashedDir()
        dir.create(save_to, showWarnings = FALSE, recursive = TRUE)
        purrr::iwalk(
          samples,
          function(sample, sample_name, save_to) {
            next_filename <- sample_rds_basenames(sample_name, save_to)
            saveRDS(sample, next_filename)
          },
          save_to = save_to,
          .progress = "Saving list of samples to disk..."
        )
        private$samples <- setNames(nm = names(samples))
      }
      self
    },
    #' @description
    #' Get a sample from the dataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @return The sample object
    getSample = function(sample) {
      self$realize()
      sample_id_num <- sample_name_or_number_to_both(sample, names(private$samples))

      filename <- paste0("sample_", sample_id_num$name, ".rds")
      current_intermediate_dir <- private$currentHashedDir()
      sample_file <- file.path(current_intermediate_dir, filename)
      if (!file.exists(sample_file)) {
        cli_abort("File not found: {.path sample_file} should exist")
      }
      out <- readRDS(sample_file)
      out <- updateObject(out)
      validObject(out)
      out
    },
    #' @description
    #' Get the path to the location of the processed samples
    #'
    #' @return The full path to the directory where samples are saved
    getCurrentDir = function() {
      private$currentHashedDir()
    }
  ),
  active = list(
    #' @field sampleNames The character vector with unique sample names. Renaming samples renames files on `obj$getCurrentDir()` as well
    sampleNames = function(value) {
      # Getter
      if (missing(value)) return(names(private$samples))
      # Setter
      if (anyNA(value) || anyDuplicated(value)) {
        cli_abort("Sample names must be unique and not missing")
      }
      if (length(value) != length(private$samples)) {
        cli_abort("sampleNames should be a vector of length {length(private$samples)} and not of length {length(value)}.")
      }
      private$rename_samples_and_files(value)
      value
    }
  ),
  private = list(
    # @field samples Named list of samples. Names are sample ids. Values are either filenames or sample objects
    samples = structure(list(), names = character(0)), # named list
    # @field scratch_dir The full directory where the samples get processed
    scratch_dir = "",
    # @field hasheddir The base of the folder name for the current location of the samples
    hasheddir = "",
    # @field keep_intermediate logical. Whether to keep intermediate results
    keep_intermediate = FALSE,
    # @description
    # The directory where processed samples are saved
    # @return A path to the directory where processed samples are kept, or `NULL`
    currentHashedDir = function() {
      hasheddir <- private$hasheddir
      if (hasheddir == "") {
        return(NULL)
      }
      file.path(
        private$scratch_dir,
        hasheddir
      )
    },
    # @description
    # The directory where samples will be saved after processing
    # @return A path
    nextHashedDir = function() {
      if (!self$hasDelayedOps()) {
        return(private$currentHashedDir())
      }

      # A hash chain:
      current_hash <- private$hasheddir
      if (current_hash == "") {
        # The first hash comes from the sample filenames and sample ids
        current_hash <- digest::digest(private$samples)
      }
      for (op in purrr::keep(private$delayed_ops, modifiesSample)) {
        current_hash <- paste0(name(op), "_", digest::digest(list(current_hash, hashableDelayedOp(op))))
      }
      file.path(
        private$scratch_dir,
        current_hash
      )
    },
    # @description
    # Implement the realize action
    # @param ... ignored
    realize_impl = function(keep_intermediate = NA, ...) {
      delayed_ops <- private$delayed_ops
      current_dir <- private$currentHashedDir()
      if (is.null(current_dir)) {
        if (name(delayed_ops[[1]]) != "read_sample") {
          cli_abort(
            message = c(
              "UnexpectedError",
              "x" = "The first operation should have been named '{.code read_sample}' instead of '{.code {name(delayed_ops[[1]])}}'",
              "i" = "This is an unexpected problem. You can try deleting the dataset and restart again",
              "i" = "If this happens again, please open an issue at {.url https://github.com/sipss/GCIMS/issues}"
            )
          )
        }
      }
      next_dir <- private$nextHashedDir()
      if (!rlang::is_string(next_dir)) {
        cli_abort(
          message = c(
            "UnexpectedError",
            "x" = "next_dir should not have been {.code {next_dir}}",
            "i" = "This is an unexpected problem. You can try deleting the dataset and restart again",
            "i" = "If this happens again, please open an issue at {.url https://github.com/sipss/GCIMS/issues}"
          )
        )
      }
      dir.create(next_dir, showWarnings = FALSE, recursive = TRUE)

      sample_names <- names(private$samples)

      if (is.null(current_dir)) {
        current_filenames <- unname(private$samples)
      } else {
        current_filenames <- sample_rds_basenames(sample_names, current_dir)
      }
      next_filenames <- sample_rds_basenames(sample_names, next_dir)
      extracted_results <- mymapply(
        FUN = realize_one_sample_disk,
        sample_name = sample_names,
        current_filename = current_filenames,
        next_filename = next_filenames,
        MoreArgs = list(
          delayed_ops = delayed_ops,
          is_first_step = is.null(current_dir),
          sample_class = private$sample_class
        ),
        SIMPLIFY = FALSE
      )

      private$aggregate_all_results(extracted_results)
      # Make delayed_ops history and remove them from pending
      private$move_delayed_to_history()

      private$hasheddir <- basename(next_dir)
      if (is.na(keep_intermediate)) {
        keep_intermediate <- private$keep_intermediate
      }
      if (!keep_intermediate && !is.null(current_dir) && current_dir != next_dir) {
        unlink(current_dir, recursive = TRUE)
      }
      return()
    },
    rename_samples_and_files = function(new_names) {
      current_dir <- private$currentHashedDir()
      if (is.null(current_dir)) {
        return()
      }
      if (!dir.exists(current_dir)) {
        return()
      }
      current_dir_old <- paste0(current_dir, "_old")
      file.rename(current_dir, current_dir_old)
      dir.create(current_dir, showWarnings = FALSE, recursive = TRUE)
      old_names <- self$sampleNames
      for (i in seq_along(new_names)) {
        old_name <- sample_rds_basenames(old_names[i])
        new_name <- sample_rds_basenames(new_names[i])
        if (file.exists(file.path(current_dir_old, old_name))) {
          file.rename(
            file.path(current_dir_old, old_name),
            file.path(current_dir, new_name)
          )
        }
      }
      # update
      names(private$samples) <- new_names
      unlink(current_dir_old, recursive = TRUE)
      invisible(NULL)
    }
  )
)

sample_rds_basenames <- function(sample_names, prefix_dir = NULL) {
  base_file <- sprintf("sample_%s.rds", sample_names)
  if (is.null(prefix_dir)) {
    base_file
  } else {
    file.path(prefix_dir, base_file)
  }
}

realize_one_sample_disk <- function(sample_name, current_filename, next_filename, delayed_ops, is_first_step, sample_class) {
  if (is_first_step) {
    # we are going to read_sample, we are starting
    sample_obj <- current_filename
  } else {
    # We load the file
    if (!file.exists(current_filename)) {
      cli_abort(
        message = c(
          "UnexpectedError",
          "x" = "The file {.path current_filename} should exist",
          "i" = "Try re running the pipeline from scratch or report the error to {.url https://github.com/sipss/GCIMS/issues}"
        )
      )
    }
    sample_obj <- readRDS(current_filename)
  }
  res <- realize_one_sample_ram(sample_obj, sample_name, delayed_ops, sample_class = sample_class)

  # saveRDS, or just copy
  if (res$needs_re_saving || is.null(current_filename)) {
    saveRDS(res$sample_obj, next_filename)
  } else if (current_filename != next_filename) {
    if (!file.copy(
      current_filename,
      next_filename
    )) {
      cli_abort("Could not copy {.path {current_filename}} to {.path {next_filename}}.")
    }
  }
  res$extracted_objects
}
