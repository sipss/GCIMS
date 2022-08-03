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
#' @slot changes_sample A logical value. If `fun` doesn't change sample, set it to `FALSE` so we avoid re-saving it.
#' @slot fun_extract A function that takes a modified [GCIMSSample] and returns an extracted object.
#' @slot fun_aggregate A function that takes a [GCIMSDataset] and a list of extracted objects and returns a modified [GCIMSDataset].
#'
#' @export
methods::setClass(
  Class = "GCIMSDelayedOp",
  slots = c(
    name = "character",
    fun = "function",
    params = "list",
    changes_sample = "logical",
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
  function(.Object, name, fun, params, changes_sample = TRUE, fun_extract = NULL, fun_aggregate = NULL) {
    .Object@name <- name
    .Object@fun <- fun
    .Object@params <- params
    .Object@changes_sample <- changes_sample
    .Object@fun_extract <- fun_extract
    .Object@fun_aggregate <- fun_aggregate
    .Object
  }
)

methods::setMethod(
  "show", "GCIMSDelayedOp",
  function(object) {
    outstring <- c(
      glue("- Delayed op {object@name}{ifelse(length(object@params) > 0, ' with parameters', '')}:"),
      paste0(purrr::imap_chr(object@params, function(v, n) glue(" - {n} = {v}")), collapse = "\n"),
    )
    cat(paste0(outstring, collapse = "\n"))
    invisible(NULL)
  }
)

