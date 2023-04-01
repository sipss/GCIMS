#' Delayed Operation class
#'
#' @description
#' DelayedOperation is an S4 class to store a delayed operation
#'
#' Delayed operations are not applied to the dataset immediately, but rather
#' when some data from the dataset is required. When working on large datasets,
#' keeping all samples in RAM may be impossible, and the [DelayedDatasetDisk]
#' architecture becomes convenient, where samples are stored in a directory, loaded
#' processed and saved individually.
#'
#' Under such arquitecture, it is more efficient to load a sample, run as many operations
#' as possible on it and save the sample, instead of loading a sample, running one
#' operation, saving the sample.
#'
#' See how to create such delayed operations and more details at
#' `vignette("creating-a-workflow-step", package = "GCIMS")`.
#'
#' @slot name A named for de delayed operation, only used for printing.
#' @slot fun A function that takes a sample object and returns a sample object, usually with some change (filtered,...)
#' @slot params A named list with additional arguments to be passed to `fun`
#' @slot params_iter A named list with additional arguments to be passed to
#' `fun`. Compared to `params`, each argument must be a named list of length the number of samples, so
#' each sample will receive its corresponding parameter according to its name
#' @slot fun_extract A function that takes the modified sample object returned
#'  by `fun` and extracts some component out of it. This component will be stored in the dataset for faster access.
#' @slot fun_aggregate A function that takes a dataset object and a list of extracted results (the output of all `fun_extract` calls)
#' and modifies the dataset.
#'
#'
#' @export
methods::setClass(
  Class = "DelayedOperation",
  slots = c(
    name = "character",
    fun = "functionOrNULL",
    params = "list",
    params_iter = "list",
    fun_extract = "functionOrNULL",
    fun_aggregate = "functionOrNULL"
  )
)

#' Create a [DelayedOperation] object
#'
#'
#' Delayed operations enables us to process our samples faster on big datasets.
#' See the details section for details on how they work.
#'
#' @details
#'
#' Let's say we have a pipeline with two actions (e.g. smooth() and detectPeaks()).
#' and we want to apply it to a dataset with two samples (e.g s1, s2).
#'
#' This is a simple pseudocode to execute all actions in all samples. The code
#' is written so you can get an idea of how :
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
#' The DelayedOperation class allows us to store all pending actions and run them
#' afterwards when the data is needed.
#'
#' Besides, samples can be processed in parallel if enough cores and RAM are
#' available.
#'
#' The DelayedOperation class also considers that sometimes we want to extract
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
#' @param fun A function that takes a sample and returns a modified sample
#' @param params A named list with additional arguments to be passed to function
#' @param params_iter A named list with additional arguments to be passed to
#' function. Compared to `params`, each argument must be a named list of length the number of samples, so
#' each sample will receive its corresponding parameter according to its name
#' @param fun_extract A function that takes a modified sample and returns an extracted object.
#' @param fun_aggregate A function that takes a dataset and a list of extracted objects and returns a modified dataset.
#' @export
#' @return A [DelayedOperation] object
DelayedOperation <- function(name, fun = NULL, params = list(), params_iter = list(), fun_extract = NULL, fun_aggregate = NULL) {
  methods::new("DelayedOperation", name = name, fun = fun, params = params, params_iter = params_iter, fun_extract = fun_extract, fun_aggregate = fun_aggregate)
}

methods::setMethod(
  "initialize", "DelayedOperation",
  function(.Object, name, fun = NULL, params = list(), params_iter = list(), fun_extract = NULL, fun_aggregate = NULL) {
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
    .Object@params_iter <- params_iter
    .Object@fun_extract <- fun_extract
    .Object@fun_aggregate <- fun_aggregate
    .Object
  }
)

hashableDelayedOp <- function(object) {
  # Functions have associated environments that change from session to session.
  # Since we don't really care about those, we just return format() of the functions
  # so the digest/hash is reproducible.
  list(
    name = object@name,
    fun = format(object@fun),
    params = object@params,
    params_iter = object@params_iter,
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

apply_op_to_sample <- function(delayed_op, sample, sample_name) {
  fun <- delayed_op@fun
  extracted_obj <- NULL
  if (!is.null(fun)) {
    params <- delayed_op@params
    params_iter <- purrr::map(delayed_op@params_iter, sample_name)
    sample <- do.call(
      fun,
      c(list(sample), params_iter, params)
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
  "show", "DelayedOperation",
  function(object) {
    cat(yaml::as.yaml(describeAsList(object)))
    invisible(NULL)
  }
)

