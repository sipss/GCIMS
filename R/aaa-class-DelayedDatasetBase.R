#' DelayedDatasetBase
#'
#' @description
#' Base R6 class to handle a dataset that can queue execute operations on it
#'
#' This is an abstract class. See [DelayedDatasetRAM] and [DelayedDatasetDisk]
#' for specific implementations.
#'
#' This class is not exported, but if you want to use it reach us at
#' <https://github.com/sipss/GCIMS/issues/> and we will export it.
#'
#' @keywords internal
DelayedDatasetBase <- R6::R6Class(
  "DelayedDatasetBase",
  public = list(
    #' @description
    #' Execute all the queued operations
    #' @param dataset The dataset where extracted results will be aggregated on
    #' @param ... Arguments used by children classes, passed on to realize_impl(...)
    realize = function(..., dataset) {
      if (!self$hasDelayedOps()) {
        return()
      }

      if (!private$can_realize) {
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
      private$can_realize <- FALSE
      on.exit({private$can_realize <- TRUE})
      private$optimize_delayed_operations()
      private$realize_impl(..., dataset = dataset)
    },
    #' @description
    #' Create a new object, initialize the slots
    #' @param sample_class The class of the samples in the dataset, used just to validate the contract between
    #' the delayed actions and the samples. If `NULL` action return values are not checked
    initialize = function(sample_class = NULL) {
      private$can_realize <- TRUE
      private$sample_class <- sample_class
    },
    #' @description
    #' The list of queued operations are usually processed sequentially over each sample when
    #' the delayed dataset is realized.
    #'
    #' In some cases we want to operate on the list of queued operations before processing them,
    #' for instance to remove redundant operations or change their order, to improve performance.
    #'
    #' Write a function that takes and returns a list of queued operations when used
    #' as follows:
    #'
    #'     queued_ops <- your_function(queued_ops)
    #'
    #' And call `obj$registerOptimization(your_function)` to register it, so every time
    #' `obj$realize()` is called, the list of queued operations will be passed to your function
    #' and optimized accordingly.
    #'
    #' You can register as many optimizations as you like.
    #'
    #' @param optimization A function that takes a list of queued operations and
    #' returns a more optimized version of it.
    #' @return The object
    #'
    registerOptimization = function(optimization) {
      private$registered_optimizations <- c(private$registered_optimizations, optimization)
      self
    },
    #' @description
    #' Queues an operation to the dataset so it will run afterwards
    #' @param operation A [DelayedOperation] object
    #' @return The modified dataset object
    appendDelayedOp = function(operation) {
      private$delayed_ops <- c(private$delayed_ops, list(operation))
      self
    },
    #' @description
    #' Find out if the dataset has pending operations
    #' @return Returns `TRUE` if the dataset has pending operations, `FALSE` otherwise
    hasDelayedOps = function() {
      length(private$delayed_ops) > 0
    },
    #' @description
    #' Get a sample from the dataset
    #'
    #' @param sample Either an integer (sample index) or a string (sample name)
    #' @param dataset The dataset object so we can realize it if needed
    #' @return The sample object
    getSample = function(sample, dataset) {
      cli_abort(
        message = c(
          "Not implemented error",
          "{.code getSample} is a virtual function and should be implemented in classes inheriting from us"
        )
      )
    },
    #' @description
    #' Get a list that describes the already executed operations
    #' @return A list, with the operations
    history_as_list = function() {
      # history info:
      # Previous operations
      pops <- private$previous_ops
      pops <- purrr::keep(pops, modifiesSample)
      pops <- purrr::map(pops, describeAsList)
      if (length(pops) > 0) {
        history_info <- list("History" = pops)
      } else {
        history_info <- "No previous history"
      }
      history_info
    },
    #' @description
    #' Get a list that describes the queued operations
    #' @return A list, with the operations
    pending_as_list = function() {
      # Pending operations
      pops <- private$delayed_ops
      pops <- purrr::keep(pops, modifiesSample)
      pops <- purrr::map(pops, describeAsList)
      if (length(pops) > 0) {
        pending_info <- list("Queued operations" = pops)
      } else {
        pending_info <- "No operations enqueued"
      }
      pending_info
    }
  ),
  active = list(
    #' @field sampleNames The character vector with unique sample names (virtual method)
    sampleNames = function(value) {
      cli_abort(
        message = c(
          "Not implemented error",
          "{.code sampleNames} is a virtual active binding and should be implemented in classes inheriting from us"
        )
      )
    }
  ),
  private = list(
    # @field previous_ops History of previous operations
    previous_ops = list(),
    # @field delayed_ops Delayed operations
    delayed_ops = list(),
    # @field can_realize Whether we can realize or not
    can_realize = FALSE,
    # @field registered_optimizations List of functions that take a list of delayed_ops and return it after optimizing it
    registered_optimizations = list(),
    # @field sample_class A string with the class of the sample objects, for validation of the functions that modify samples.
    sample_class = NULL,
    # @description
    # Implement the realize action
    # @param ... ignored
    realize_impl = function(..., dataset) {
      cli_abort(
        message = c(
          "Not implemented error",
          "{.code realize_impl(..., dataset)} is a virtual function and should be implemented in classes inheriting from us"
        )
      )
    },
    optimize_delayed_operations = function() {
      if (!self$hasDelayedOps() || length(private$registered_optimizations) == 0) {
        return()
      }
      delayed_ops <- private$delayed_ops
      for (optim in private$registered_optimizations) {
        delayed_ops <- optim(delayed_ops)
      }
      private$delayed_ops <- delayed_ops
      return()
    },
    move_delayed_to_history = function() {
      # Make delayed_ops history and remove them from pending
      private$previous_ops <- c(private$previous_ops, private$delayed_ops)
      private$delayed_ops <- list()
      invisible(NULL)
    },
    aggregate_all_results = function(extracted_results, dataset) {
      for (i in seq_along(private$delayed_ops)) {
        delayed_op <- private$delayed_ops[[i]]
        # Extract i-th result for all samples:
        extracted_result <- purrr::map(extracted_results, i)
        # Let the delayed operation aggregate the extracted results and save them in dataset
        private$aggregate_result(delayed_op, extracted_result, dataset)
      }
      return()
    },
    aggregate_result = function(delayed_op, extracted_result, dataset) {
      if (is.null(delayed_op@fun_aggregate)) {
        return()
      }
      f <- delayed_op@fun_aggregate
      f(dataset, extracted_result)
      return()
    }
  )
)
