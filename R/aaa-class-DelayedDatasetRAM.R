#' DelayedDatasetRAM
#'
#' @description
#' A class that contains a dataset where its samples are stored in memory (RAM).
#' This is usually faster than storing them on files, but the size of the
#' dataset becomes limited by the amount of RAM available.
#'
#' @details
#'
#' This class is not exported, but if you want to use it reach us at
#' <https://github.com/sipss/GCIMS/issues/> and we will export it.
#'
#' @keywords internal
DelayedDatasetRAM <- R6::R6Class(
  "DelayedDatasetRAM",
  inherit = DelayedDatasetBase,
  public = list(
    #' @description
    #' Create a delayed dataset on RAM
    #' @param samples A named list of samples.
    #' @param sample_class The class of the samples in the dataset, used just to validate the contract between
    #' the delayed actions and the samples. If `NULL` action return values are not checked
    initialize = function(samples, sample_class = NULL) {
      sample_names <- names(samples)
      if (length(samples) > 0 && is.null(sample_names)) {
        cli_abort("samples should be a named list")
      }
      if (any(names(samples) == "") || any(duplicated(names(samples)))) {
        cli_abort("sample names should be unique and non-empty")
      }
      private$samples <- samples
      super$initialize(sample_class = sample_class)
    },
    #' @description
    #' Get a sample from the dataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @param dataset The dataset so we can realize if there are enqueued actions
    #' @return The sample object
    getSample = function(sample, dataset) {
      self$realize(dataset = dataset)
      sample_id_num <- sample_name_or_number_to_both(sample, names(private$samples))
      out <- private$samples[[sample_id_num$idx]]
      out <- updateObject(out)
      validObject(out)
      out
    },
    #' @description
    #' Subsets some samples
    #' @param sample A numeric vector with indices, a character vector with names or a logical vector
    #' @return the delayed dataset modified in-place
    subset = function(sample) {
      if (self$hasDelayedOps()) {
        cli_abort("Can't subset a delayed dataset with pending operations")
      }
      sample_info <- sample_name_or_number_to_both(sample, self$sampleNames)
      private$samples <- private$samples[sample_info$name]
      self
    }
  ),
  active = list(
    #' @field sampleNames The character vector with unique sample names.
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
      names(private$samples) <- value
      value
    }
  ),
  private = list(
    # @field samples Named list of samples. Names are sample ids. Values are either filenames or sample objects
    samples = structure(list(), names = character(0)), # named list
    # @description
    # Implement the realize action
    # @param ... ignored
    # @param dataset A dataset object to be modified in place by the delayed operations
    realize_impl = function(..., dataset) {
      sample_names <- names(private$samples)
      results <- BiocParallel::bpmapply(
        FUN = realize_one_sample_ram,
        sample_name = names(private$samples),
        sample_obj = private$samples,
        MoreArgs = list(
          delayed_ops = private$delayed_ops,
          sample_class = private$sample_class
        ),
        SIMPLIFY = FALSE
      )
      private$samples <- purrr::map(results, "sample_obj")
      names(private$samples) <- sample_names

      extracted_results <- purrr::map(results, "extracted_objects")

      private$aggregate_all_results(extracted_results, dataset = dataset)
      private$move_delayed_to_history()
      return()
    }
  )
)


realize_one_sample_ram <- function(sample_obj, sample_name, delayed_ops, sample_class) {
  out <- vector("list", length = length(delayed_ops))
  needs_re_saving <- FALSE
  # Execute:
  for (i in seq_along(delayed_ops)) {
    result <- apply_op_to_sample(delayed_ops[[i]], sample_obj, sample_name)
    sample_obj <- result[["sample"]]
    if (!is.null(sample_class) && !inherits(sample_obj, sample_class)) {
      cli_abort(
        c("Invalid action {delayed_ops[[i]]@name}",
          "i" = "The action should have returned a sample of class {sample_name}, but the object was of class {class(sample_obj)}",
          "i" = "If you did not write the action, please report this to {.url https://github.com/sipss/GCIMS/issues}"
          )
      )
    }
    out[i] <- result["extracted_obj"]
    if (result$needs_resaving) {
      needs_re_saving <- TRUE
    }
  }
  list(
    extracted_objects = out,
    sample_obj = sample_obj,
    needs_re_saving = needs_re_saving
  )
}
