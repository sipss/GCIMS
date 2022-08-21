setClassUnion("functionOrNULL", c("function", "NULL"))

#' GCIMSDelayedOp class
#'
#' @description
#' GCIMSDelayedOp is an S4 class to store a delayed operation
#'
#' Delayed operations are not applied to the dataset immediately, but rather
#' when some data from the dataset is required. Delaying the calculation has
#' the advantage that some operations that are applied to each sample can be applied
#' sequentially to each sample, reducing the amount of reads and writes to
#' the data files, and increasing the performance.
#'
#' @slot name A named for de delayed operation, only used for printing.
#' @slot fun A function that takes a [GCIMSSample] and returns a [GCIMSSample] (modified)
#' @slot params A named list with additional arguments to be passed to function
#' @slot fun_extract A function that takes a modified [GCIMSSample] and returns an extracted object.
#' @slot fun_aggregate A function that takes a [GCIMSDataset] and a list of extracted objects and returns a modified [GCIMSDataset].
#'
#' @export
methods::setClass(
  Class = "GCIMSDelayedOp",
  slots = c(
    name = "character",
    fun = "functionOrNULL",
    params = "list",
    fun_extract = "functionOrNULL",
    fun_aggregate = "functionOrNULL"
  )
)

#' @describeIn GCIMSDelayedOp class
#'
#' @param ... See the slots section in this page
#' @return A GCIMSDelayedOp object
GCIMSDelayedOp <- function(...) {
  methods::new("GCIMSDelayedOp", ...)
}

methods::setMethod(
  "initialize", "GCIMSDelayedOp",
  function(.Object, name, fun = NULL, params = list(), fun_extract = NULL, fun_aggregate = NULL) {
    if (!grepl(pattern = "^[a-zA-Z0-9][-a-zA-Z_0-9]*$", x = name)) {
      rlang::abort(
        message = c(
          "Invalid delayed operation name",
          "x" = glue("{name} is not a valid operation name"),
          "i" = "The operation name must start with an alphanumeric character",
          "i" = "The operation name may only contain alphanumeric characters hyphens and underscores"
          )
      )
    }
    .Object@name <- name
    .Object@fun <- fun
    .Object@params <- params
    .Object@fun_extract <- fun_extract
    .Object@fun_aggregate <- fun_aggregate
    .Object
  }
)

aggregate_result <- function(delayed_op, extracted_result, dataset) {
  if (is.null(delayed_op@fun_aggregate)) {
    return(dataset)
  }
  dataset_class <- class(dataset)
  f <- delayed_op@fun_aggregate
  dataset <- f(dataset, extracted_result)
  if (!methods::is(dataset, dataset_class)) {
    rlang::abort(
      message = c(
        "Delayed operation contract was broken",
        "x" = glue("The delayed action {name(delayed_op)} has a `fun_aggregate` slot that does not return a {dataset_class} object"),
        "i" = "If you did not write the delayed action, this is not your fault. Please report this error at https://github.com/sipss/GCIMS."
      )
    )
  }
  dataset
}


name <- function(x) {
  x@name
}

modifiesSample <- function(delayed_op) {
  # If there is a fun, we assume it modifies the sample
  !is.null(delayed_op@fun)
}

apply_op_to_sample <- function(delayed_op, sample) {
  fun <- delayed_op@fun
  extracted_obj <- NULL
  if (!is.null(fun)) {
    params <- delayed_op@params
    sample <- do.call(
      fun,
      c(list(sample), params)
    )
  }
  if (!is.null(delayed_op@fun_extract)) {
    f <- delayed_op@fun_extract
    extracted_obj <- f(sample)
  }
  list(
    sample = sample,
    needs_resaving = modifiesSample(delayed_op),
    extracted_obj = extracted_obj
  )
}

methods::setMethod(
  "describeAsList", "GCIMSDelayedOp",
  function(object) {
    num_params <- length(object@params)
    txt <- name(object)
    if (num_params == 0) {
      return(txt)
    }
    p <- object@params
    for (i in seq_along(p)) {
      if (is.atomic(p[[i]]) && length(p[[i]]) > 10) {
        p[[i]] <- glue("< A {mode(p[[i]])} vector of {length(p[[i]])} elements >")
      }
    }
    out <- list()
    out[[txt]] <- p
    return(out)
  }
)

methods::setMethod(
  "show", "GCIMSDelayedOp",
  function(object) {
    cat(yaml::as.yaml(describeAsList(object)))
    invisible(NULL)
  }
)

