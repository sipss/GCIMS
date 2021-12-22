#' GCIMS dataset class
#'
#' GCIMS is an S4 class to store one or more samples
#' with the same drift and retention time ranges.
#'
#' The actual spectra is stored in the `data` slot, in a 3-D array,
#' where the first index corresponds to drift time, the second to retention time
#' and the third one to the sample.
#'
#' @slot drift_time numeric.
#' @slot retention_time numeric.
#' @slot metadata data.frame.
#' @slot data array. A cubic array with drift time in the rows, retention time
#'  in columns and the sample index in the third dimension
#'
#' @return
#' @export
#'
methods::setClass(
  Class = "GCIMS",
  slots = c(
    drift_time = "numeric",
    retention_time = "numeric",
    metadata = "data.frame",
    data = "array" # [dt, rt, sample],
    # long tubo deriva, ... drift gas...
  )
)


#' @name GCIMS-methods
#' @title Methods for the GC-IMS class
NULL


setMethod("length", "GCIMS", function(x) dim(x)[3])

setGeneric("driftTime", function(x) standardGeneric("driftTime"))
setMethod("driftTime", "GCIMS", function(x) x@drift_time)

setGeneric("retentionTime", function(x) standardGeneric("retentionTime"))
setMethod("retentionTime", "GCIMS", function(x) x@retention_time)

setGeneric("intensity", function(x, ...) standardGeneric("intensity"))
setMethod("intensity", "GCIMS", function(x, dtrange = NULL, rtrange = NULL, sample_idx = NULL) {
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

  if (is.null(sample_idx)) {
    sample_idx <- seq_len(length(x))
  }
  x@data[dt_idx, rt_idx, sample_idx]
})

setMethod("show",
          "GCIMS",
          function(object) {
            cat(sprintf("A GCIMS dataset with %d samples\n", dim(object)[3]))
          }
)

setValidity("GCIMS", function(object) {
  success <- TRUE
  issues <- c()
  dt <- driftTime(object)
  rt <- retentionTime(object)
  metadata <- object@metadata
  x <- object@data
  if (dim(x)[1] != length(dt)) {
    issues <- c(issues, "Length of drift_time does not match first data dimension")
    success <- FALSE
  }
  if (dim(x)[2] != length(rt)) {
    issues <- c(issues, "Length of retention_time does not match second data dimension")
    success <- FALSE
  }
  if (dim(x)[3] != nrow(metadata)) {
    issues <- c(issues, "The number of rows of metadata does not match the third data dimension")
    success <- FALSE
  }
  if (!success) {
    return(issues)
  }
  success
})



#' @describeIn GCIMS-methods Autocompletion for \code{$} access on a
#' \code{GCIMS} object
#'
#' @inheritParams utils::.DollarNames
#'
#' @importFrom utils .DollarNames
#' @export
#' @method .DollarNames GCIMS
#'
".DollarNames.GCIMS" <- function(x, pattern = '') {
  meta.data <- as.list(x = colnames(x = x@metadata))
  names(x = meta.data) <- unlist(x = meta.data)
  return(.DollarNames(x = meta.data, pattern = pattern))
}

#' @describeIn GCIMS-methods Metadata access for \code{GCIMS} objects
#'
#' @return \code{$}: metadata column \code{i} for object \code{x};
#' \strong{note}: unlike \code{[[}, \code{$} drops the shape of the metadata
#' to return a vector instead of a data frame
#'
#' @export
#' @method $ GCIMS
#'
#' @examples
#' \dontrun{
#' # Get metadata using `$'
#' head(mydataset$groups)
#' }
#'
"$.GCIMS" <- function(x, i, ...) {
  return(x[[i, drop = TRUE]])
}


#' @describeIn GCIMS-methods Metadata setter for \code{GCIMS} objects
#'
#' @return \code{$<-}: object \code{x} with metadata \code{value} saved as
#' \code{i}
#'
#' @export
#' @method $<- GCIMS
#'
#' @examples
#' # Add metadata using the `$' operator
#' set.seed(42)
#' mydataset$value <- sample(1:3, size = length(mydataset), replace = TRUE)
#' head(mydataset[["value"]])
#'
"$<-.GCIMS" <- function(x, i, ..., value) {
  x[[i]] <- value
  return(x)
}


#' @describeIn GCIMS-methods Simple subsetter for \code{GCIMS} objects
#'
#' @return \code{[}: object \code{x} with features \code{i} and cells \code{j}
#'
#' @export
#' @method [ GCIMS
#'
#' @examples
#' # `[' examples
#'
"[.GCIMS" <- function(x, i, j, k, ...) {
  if (missing(x = i) && missing(x = j) && missing(x = k)) {
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
  if (missing(x = k)) {
    k <- seq_len(length(x))
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
  if (is.logical(x = k)) {
    if (length(k) != length(x)) {
      stop("Incorrect number of logical values provided to subset samples")
    }
  }
  return(subset.GCIMS(x = x, dt_idx = i, rt_idx = j, samples_idx = k, ...))
}


dim.GCIMS <- function(x) {
  return(dim(x@data))
}


dimnames.GCIMS <- function(x) {
  return(dimnames(x@data))
}



#' @describeIn GCIMS-methods Subset a \code{\link{GCIMS}} object
#'
#' @param subset Logical expression metadata filters to keep
#' @param dt_idx A vector of drift time indices to keep
#' @param rt_idx A vector of retention time indices to keep
#' @param samples_idx A vector of sample indices to keep
#'
#' @return \code{subset}: A subsetted \code{GCIMS} object
#'
#' @importFrom rlang enquo
#
#' @aliases subset
#' @seealso \code{\link[base]{subset}} \code{\link{WhichSamples}}
#'
#' @export
#' @method subset GCIMS
#'
#' @examples
#' # `subset' examples
#' subset(mydataset, subset = Sex > 'Male')
#' subset(mydataset, subset = Sex > 'Male', rt = 1:100, dt = 10:20)
#'
subset.GCIMS <- function(
  x,
  subset,
  dt_idx = NULL,
  rt_idx = NULL,
  samples_idx = NULL,
  ...
) {
  if (!missing(x = subset)) {
    subset <- rlang::enquo(arg = subset)
    metadata <- x@metadata
    metadata$tmp_row_idx <- seq_len(nrow(metadata))
    samples_idx_subset <- dplyr::filter(metadata, !!subset)$tmp_row_idx
    if (is.null(samples_idx)) {
      samples_idx <- samples_idx_subset
    } else {
      samples_idx <- intersect(samples_idx, samples_idx_subset)
    }
  }

  dt <- driftTime(x)
  rt <- retentionTime(x)

  if (is.null(dt_idx)) {
    dt_idx <- seq_along(dt)
  }
  if (is.null(rt_idx)) {
    rt_idx <- seq_along(rt)
  }
  if (is.null(samples_idx)) {
    samples_idx <- seq_len(length(x))
  }

  if (all(dt_idx == seq_along(dt)) &&
      all(rt_idx == seq_along(rt)) &&
      all(samples_idx == seq_len(length(x)))) {
    return(x)
  }

  new(
    "GCIMS",
    drift_time = dt[dt_idx],
    retention_time = rt[rt_idx],
    metadata = x@metadata[samples_idx, , drop = FALSE],
    data = x@data[dt_idx, rt_idx, samples_idx, drop = FALSE]
  )
}
