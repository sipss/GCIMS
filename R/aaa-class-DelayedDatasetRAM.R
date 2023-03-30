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
    #' @param dataset The dataset that stores aggregated results from the processing
    initialize = function(samples, dataset) {
      sample_names <- names(samples)
      if (length(samples) > 0 && is.null(sample_names)) {
        cli_abort("samples should be a named list")
      }
      if (any(names(samples) == "") || any(duplicated(names(samples)))) {
        cli_abort("sample names should be unique and non-empty")
      }
      private$samples <- samples
      super$initialize(dataset = dataset)
    },
    #' @description
    #' Get a sample from the dataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @return The sample object
    getSample = function(sample) {
      self$realize()
      sample_id_num <- sample_name_or_number_to_both(sample, names(private$samples))
      out <- private$samples[[sample_id_num$idx]]
      out <- updateObject(out)
      validObject(out)
      out
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
      # FIXME: Set sample description
      value
    }
  ),
  private = list(
    # @field samples Named list of samples. Names are sample ids. Values are either filenames or sample objects
    samples = structure(list(), names = character(0)), # named list
    # @description
    # Implement the realize action
    # @param ... ignored
    realize_impl = function(...) {
      sample_names <- names(private$samples)
      results <- mymapply(
        FUN = realize_one_sample_ram,
        sample_name = names(private$samples),
        sample_obj = private$samples,
        MoreArgs = list(
          delayed_ops = private$delayed_ops
        ),
        SIMPLIFY = FALSE
      )
      private$samples <- purrr::map(results, "sample_obj")
      names(private$samples) <- sample_names

      extracted_results <- purrr::map(results, "extracted_objects")

      private$aggregate_all_results(extracted_results)
      private$move_delayed_to_history()
      return()
    }
  )
)


realize_one_sample_ram <- function(sample_obj, sample_name, delayed_ops) {
  out <- vector("list", length = length(delayed_ops))
  needs_re_saving <- FALSE
  # Execute:
  for (i in seq_along(delayed_ops)) {
    result <- apply_op_to_sample(delayed_ops[[i]], sample_obj, sample_name)
    sample_obj <- result[["sample"]]
    # FIXME: Validate sample_obj
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
