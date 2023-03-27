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

#' Create a [GCIMSDelayedOp] object
#'
#'
#' Delayed operations enable GCIMS to process our samples faster on big datasets.
#' See the details section for details on how they work.
#'
#' @details
#'
#' Let's say we have a pipeline with two actions (e.g. smooth() and detectPeaks()).
#' and we want to apply it to a dataset with two samples (e.g s1, s2).
#'
#' This is a simple pseudocode to execute all actions in all samples:
#'
#' ```
#' dataset = list(s1, s2)
#' actions = list(smooth, detectPeaks)
#' for (action in actions) {
#'   for (i in seq_along(dataset)) {
#'       dataset[[i]] <- action(dataset[[i]])
#'   }
#' }
#' ```
#'
#' When the dataset is big, samples are stored in disk, and loaded/saved when used:
#'
#' ```
#' dataset = list(s1, s2)
#' actions = list(smooth, detectPeaks)
#' for (action in actions) {
#'   for (i in seq_along(dataset)) {
#'       sample <- read_from_disk(i)
#'       sample <- action(sample)
#'       save_to_disk(sample)
#'   }
#' }
#' ```
#'
#' So actually, we can avoid "saving and loading" by changing the loop order:
#'
#' ```
#' dataset = list(s1, s2)
#' actions = list(smooth, detectPeaks)
#' for (i in seq_along(dataset)) {
#'   sample <- read_from_disk(i)
#'   for (action in actions) {
#'       sample <- action(sample)
#'   }
#'   save_to_disk(sample)
#' }
#' ```
#'
#' This requires that when we apply an operation to the dataset, the operation
#' is delayed, so we can stack many delayed operations and run them all at once.
#'
#' The GCIMSDelayedOp class allows us to store all pending actions and run them
#' afterwards when the data is needed.
#'
#' Besides, samples can be processed in parallel if enough cores and RAM are
#' available.
#'
#' The GCIMSDelayedOp class also considers that sometimes we want to extract
#' some information from each sample (e.g. the Reverse Ion Chromatogram)
#' and build some matrix with the Reverse Ion Chromatograms of all samples. It changes
#' the loops above, so after each action modifies each sample, we can extract something
#' out of the sample and save it. After all actions have executed, we can aggregate
#' the results we have extracted and save them into the dataset. This is used for instance
#' in the [getRIC()] implementation, to extract the RIC from each sample and afterwards
#' aggregate it into a matrix. This is implemented here with the `fun_extract` and
#' `fun_aggregate` functions.
#'
#'
#' @param name A named for de delayed operation, only used for printing.
#' @param fun A function that takes a [GCIMSSample] and returns a [GCIMSSample] (modified)
#' @param params A named list with additional arguments to be passed to function
#' @param fun_extract A function that takes a modified [GCIMSSample] and returns an extracted object.
#' @param fun_aggregate A function that takes a [GCIMSDataset] and a list of extracted objects and returns a modified [GCIMSDataset].
#' @return A GCIMSDelayedOp object
GCIMSDelayedOp <- function(name, fun = NULL, params = list(), fun_extract = NULL, fun_aggregate = NULL) {
  methods::new("GCIMSDelayedOp", name = name, fun = fun, params = params, fun_extract = fun_extract, fun_aggregate = fun_aggregate)
}

methods::setMethod(
  "initialize", "GCIMSDelayedOp",
  function(.Object, name, fun = NULL, params = list(), fun_extract = NULL, fun_aggregate = NULL) {
    if (!grepl(pattern = "^[a-zA-Z0-9][-a-zA-Z_0-9]*$", x = name)) {
      cli_abort(
        message = c(
          "Invalid delayed operation name",
          "x" = "{name} is not a valid operation name",
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

# FIXME: Create an R6 parent class of GCIMSDataset named DelayedDataset with
# all the realize and delayed stuff
# FIXME: This should be a private method of GCIMSDataset
aggregate_result <- function(delayed_op, extracted_result, dataset) {
  if (is.null(delayed_op@fun_aggregate)) {
    return(dataset)
  }
  dataset_class <- class(dataset)
  f <- delayed_op@fun_aggregate
  dataset <- f(dataset, extracted_result)
  if (!inherits(dataset, dataset_class)) {
    cli_abort(
      message = c(
        "Delayed operation contract was broken",
        "x" = "The delayed action {name(delayed_op)} has a `fun_aggregate` slot that does not return a {dataset_class} object",
        "i" = "If you did not write the delayed action, this is not your fault. Please report this error at https://github.com/sipss/GCIMS."
      )
    )
  }
  dataset
}

hashableDelayedOp <- function(object) {
  # Functions have associated environments that change from session to session.
  # Since we don't really care about those, we just return format() of the functions
  # so the digest/hash is reproducible.
  list(
    name = object@name,
    fun = format(object@fun),
    params = object@params,
    fun_extract = format(object@fun_extract),
    fun_aggregate = format(object@fun_aggregate)
  )
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

#' Format a delayed operation as a list
#' @noRd
#' @return A list with a brief description/representation of the object
#'
describeAsList <- function(object) {
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
    if (inherits(p[[i]], "data.frame")) {
      p[[i]] <- glue("< A data.frame with {nrow(p[[i]])} rows and {ncol(p[[i]])} columns >")
    }
    if (is.function(p[[i]])) {
      p[[i]] <- "< function >"
    }
  }
  out <- list()
  out[[txt]] <- p
  return(out)
}

methods::setMethod(
  "show", "GCIMSDelayedOp",
  function(object) {
    cat(yaml::as.yaml(describeAsList(object)))
    invisible(NULL)
  }
)

