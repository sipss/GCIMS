#' GCIMSSample class
#'
#' GCIMS Sample is an S4 class to store one sample
#' with the drift and retention time ranges and other relevant attributes
#' (GC column, drift tube length...) if available
#'
#' The actual spectra is stored in the `data` slot, in a matrix,
#' where the first index (rows) corresponds to drift times and the
#' second to retention times (columns).
#'
#'
#' @slot drift_time numeric. (required)
#' @slot retention_time numeric. (required)
#' @slot data matrix A matrix with drift time in the rows and retention time
#'  in columns. (required)
#'
#' @return A GCIMSSample object
#' @export
#'
methods::setClass(
  Class = "GCIMSSample",
  slots = c(
    drift_time = "numeric",
    retention_time = "numeric",
    data = "matrix", # [dt, rt]
    gc_column = "character",
    drift_tube_length = "numeric",
    params = "list" # arbitrary parameters from the instrument
  )
)


#' @name GCIMSSample-methods
#' @title Methods for the GCIMSSample class
NULL

setGeneric("driftTime", function(x) standardGeneric("driftTime"))
setMethod("driftTime", "GCIMSSample", function(x) x@drift_time)

setGeneric("retentionTime", function(x) standardGeneric("retentionTime"))
setMethod("retentionTime", "GCIMSSample", function(x) x@retention_time)

setGeneric("intensity", function(x, ...) standardGeneric("intensity"))
setMethod("intensity", "GCIMSSample", function(x, dtrange = NULL, rtrange = NULL) {
  dt <- driftTime(x)
  if (!is.null(dtrange)) {
    dtmin <- min(dtrange)
    dtmax <- max(dtrange)
    dt_idx <- dt >= dtmin & dt < dtmax
  } else {
    dt_idx <- seq_along(dt)
  }

  rt <- retentionTime(x)
  if (!is.null(rtrange)) {
    rtmin <- min(rtrange)
    rtmax <- max(rtrange)
    rt_idx <- rt >= rtmin & rt < rtmax
  } else {
    rt_idx <- seq_along(rt)
  }
  x@data[dt_idx, rt_idx]
})

setMethod("show",
          "GCIMSSample",
          function(object) {
            cat(sprintf("A GCIMS Sample\n"))
            # FIXME: Give more details about the sample
          }
)

setValidity("GCIMSSample", function(object) {
  success <- TRUE
  issues <- c()
  dt <- driftTime(object)
  rt <- retentionTime(object)
  x <- object@data
  if (!inherits(x, "matrix")) {
    issues <- c(issues, "data slot should have a matrix")
  }
  if (nrow(x) != length(dt)) {
    issues <- c(issues, "Length of drift_time does not match the number of rows in data")
    success <- FALSE
  }
  if (ncol(x) != length(rt)) {
    issues <- c(issues, "Length of retention_time does not match the number of columns in data")
    success <- FALSE
  }
  if (!success) {
    return(issues)
  }
  success
})


#' @describeIn GCIMSSample-methods Simple subsetter for \code{GCIMSSample} objects
#'
#' @return \code{[}: object \code{x} with features \code{i} and cells \code{j}
#'
#' @param i index for drift time to subset
#' @param j index for retention time to subset
#' @param ... ignored
#' @export
#' @method [ GCIMS
#'
#' @examples
#' # `[' examples
#'
"[.GCIMS" <- function(x, i, j, ...) {
  if (missing(x = i) && missing(x = j)) {
    return(x)
  }
  dt <- driftTime(x)
  rt <- retentionTime(x)
  if (missing(x = i)) {
    i <- seq_along(dt)
  }
  if (missing(x = j)) {
    j <- seq_along(rt)
  }

  if (is.logical(x = i)) {
    if (length(i) != length(x = dt)) {
      stop("Incorrect number of logical values provided to subset drift time")
    }
  }
  if (is.logical(x = j)) {
    if (length(j) != length(x = rt)) {
      stop("Incorrect number of logical values provided to subset retention time")
    }
  }
  return(subset.GCIMSSample(x = x, dt_idx = i, rt_idx = j, ...))
}


dim.GCIMSSample <- function(x) {
  return(dim(x@data))
}


dimnames.GCIMSSample <- function(x) {
  return(dimnames(x@data))
}


#' @describeIn GCIMSSample-methods Subset a \code{\link{GCIMSSample-class}} object
#'
#' @param dt_idx A vector of drift time indices to keep
#' @param rt_idx A vector of retention time indices to keep
#'
#' @return \code{subset}: A subsetted \code{GCIMSSample} object
#'
#' @importFrom rlang enquo
#
#' @aliases subset
#' @seealso \code{\link[base]{subset}}
#'
#' @export
#' @method subset GCIMSSample
#'
subset.GCIMSSample <- function(
  x,
  dt_idx = NULL,
  rt_idx = NULL,
  ...
) {

  dt <- driftTime(x)
  rt <- retentionTime(x)

  if (is.null(dt_idx)) {
    dt_idx <- seq_along(dt)
  }
  if (is.null(rt_idx)) {
    rt_idx <- seq_along(rt)
  }

  if (all(dt_idx == seq_along(dt)) &&
      all(rt_idx == seq_along(rt))) {
    return(x)
  }

  new_obj <- x
  new_obj@drift_time <- dt[dt_idx]
  new_obj@retention_time <- rt[rt_idx]
  new_obj@data <- x@data[dt_idx, rt_idx, drop = FALSE]
  new_obj
}


#' @describeIn GCIMSSample-methods Add metadata
#'
#' @param value A vector of length number of samples with metadata to add; \strong{note}:
#' can pass \code{NULL} to remove metadata or an associated object
#'
#' @return \code{[[<-}: \code{x} with the metadata or associated objects added
#' as \code{i}; if \code{value} is \code{NULL}, removes metadata \code{i}
#' from object \code{x}
#'
#' @export
#'
setMethod( # because R doesn't allow S3-style [[<- for S4 classes
  f = '[[<-',
  signature = c('x' = 'GCIMSSample'),
  definition = function(x, i, ..., value) {
    if (!is.character(x = i)) {
      stop("'i' must be a character", call. = FALSE)
    }
    x@metadata[[i]] <- value
    x
  }
)
