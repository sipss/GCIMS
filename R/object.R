#' GCIMSSample class
#'
#' @description
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
#' @slot gc_column character. (optional) The type of chromatographic column used
#' @slot drift_tube_length numeric (optional) The length of the drift tube, in mm
#' @slot drift_gas character. (optional) The drift gas used (e.g "nitrogen")
#' @slot params list (optional) Arbitrary list of parameters and annotations
#'
#' @export
#' @seealso [GCIMSSample-methods]
#' @examples
#' # Create a new GCIMSSample with methods::new()
#' dummy_obj <-methods::new(
#'   "GCIMSSample",
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3),
#'   gc_column = "Optional column name",
#'   drift_gas = "nitrogen",
#'   drift_tube_length = 98.0 # in mm
#' )
methods::setClass(
  Class = "GCIMSSample",
  slots = c(
    drift_time = "numeric",
    retention_time = "numeric",
    data = "matrix", # [dt, rt]
    gc_column = "character",
    drift_gas = "character",
    drift_tube_length = "numeric",
    params = "list" # arbitrary parameters from the instrument
  )
)



#' @describeIn GCIMSSample Create a GCIMSSample object
#'
#' @param ... See the slots section in this page
#' @export
#'
#' @examples
#' # Create a new GCIMSSample with the convenient constructor function:
#' dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3),
#'   gc_column = "Optional column name",
#'   drift_gas = "nitrogen",
#'   drift_tube_length = 98.0 # in mm
#' )
GCIMSSample <- function(
  ...
) {
  methods::new("GCIMSSample", ...)
}

#' @name GCIMSSample-methods
#' @title Methods for the GCIMSSample class
NULL


#' @describeIn GCIMSSample Get drift time numeric vector
#'
#' @return A numeric vector with the drift time
#'
#' @param object A GCIMSSample object
#' @export
#'
setGeneric("driftTime", function(object) standardGeneric("driftTime"))
setMethod("driftTime", "GCIMSSample", function(object) object@drift_time)

#' @describeIn GCIMSSample Get retention time numeric vector
#'
#' @return A numeric vector with the retention time
#'
#' @param object A GCIMSSample object
#' @export
#'
setGeneric("retentionTime", function(object) standardGeneric("retentionTime"))
setMethod("retentionTime", "GCIMSSample", function(object) object@retention_time)


#' @describeIn GCIMSSample Get the intensity matrix
#'
#' @return A numeric matrix with drift time in the rows and retention time in columns
#'
#' @param object A GCIMSSample object
#' @param dtrange The minimum and maximum drift times to extract (length 2 vector)
#' @param rtrange The minimum and maximum retention times to extract (length 2 vector)
#' @export
#' @examples
#' mea_file <- system.file("extdata",  "sample_formats", "small.mea.gz", package = "GCIMS")
#' gcims_sample <- read_mea(mea_file)
#' my_matrix <- intensityMatrix(gcims_sample, dtrange = c(7, 8), rtrange = c(1,30))
#'
setGeneric("intensityMatrix", function(object, ...) standardGeneric("intensityMatrix"))
setMethod("intensityMatrix", "GCIMSSample", function(object, dtrange = NULL, rtrange = NULL) {
  dt <- driftTime(object)
  if (!is.null(dtrange)) {
    dtmin <- min(dtrange)
    dtmax <- max(dtrange)
    dt_idx <- dt >= dtmin & dt < dtmax
  } else {
    dt_idx <- seq_along(dt)
  }

  rt <- retentionTime(object)
  if (!is.null(rtrange)) {
    rtmin <- min(rtrange)
    rtmax <- max(rtrange)
    rt_idx <- rt >= rtmin & rt < rtmax
  } else {
    rt_idx <- seq_along(rt)
  }
  out <- object@data[dt_idx, rt_idx]
  rownames(out) <- dt[dt_idx]
  colnames(out) <- rt[rt_idx]
  out
})

setMethod(
  "show",
  "GCIMSSample",
  function(object) {
    axes <- list(
      "drift time" = list(value=driftTime(object), unit="ms"),
      "retention time" = list(value=retentionTime(object), unit="s")
    )
    outstring <- "A GCIMS Sample"
    for (axis_name in  names(axes)) {
      axis <- axes[[axis_name]][["value"]]
      axis_unit <- axes[[axis_name]][["unit"]]
      if (length(axis) == 0) {
        first <- NaN
        last <- NaN
        res <- NaN
      } else if (length(axis) == 1) {
        first <- axis[1]
        last <- axis[1]
        res <- NaN
      } else {
        first <- axis[1]
        last <- axis[length(axis)]
        res <- axis[2]-axis[1]
      }
      outstring <- c(
        outstring,
        paste0(" with ", axis_name, " from ", first, " to ", last, " ", axis_unit,
               " (", "step: ", res, " ", axis_unit , ", ",
               "points: ", length(axis), ")")
      )
    }
    # FIXME: Give more details about the sample
    cat(paste0(outstring, collapse = "\n"))
  }
)

setValidity("GCIMSSample", function(object) {
  success <- TRUE
  issues <- c()
  dt <- driftTime(object)
  rt <- retentionTime(object)
  x <- object@data
  if (length(dt) == 0) {
    issues <- c(issues, "drift_time should be a numeric vector")
    success <- FALSE
  }
  if (length(rt) == 0) {
    issues <- c(issues, "retention_time should be a numeric vector")
    success <- FALSE
  }
  if (!inherits(x, "matrix")) {
    issues <- c(issues, "data slot should have a matrix")
    success <- FALSE
  }
  if (nrow(x) != length(dt) || nrow(x) == 0) {
    issues <- c(issues, "data should have as many rows as the drift_time length")
    success <- FALSE
  }
  if (ncol(x) != length(rt) || ncol(x) == 0) {
    issues <- c(issues, "data should have as many columns as the retention_time length")
    success <- FALSE
  }
  if (!success) {
    return(issues)
  }
  success
})

#' @describeIn GCIMSSample Topographical plot of a GC-IMS Sample
#'
#' @return
#'
#' @param ... passed to [graphics::filled.contour]
#' @export
#' @method image GCIMSSample
#'
#' @examples
#' # `[' examples
#'
"image.GCIMSSample" <- function(x, ...) {
  # filled.contour includes a legend
  graphics::filled.contour(
    x=x@drift_time,
    y = x@retention_time/60.0,
    z=x@data,
    xlab = "Drift time (ms)",
    ylab = "Retention time (min)",
    ...
  )
}


#' @describeIn GCIMSSample-methods Simple subsetter for [GCIMSSample-class] objects
#'
#' @return `[`: object `x` with features `i` and cells `j`
#'
#' @param i index for drift time to subset
#' @param j index for retention time to subset
#' @param ... ignored
#' @export
#' @method [ GCIMSSample
#'
#' @examples
#' # `[' examples
#'
"[.GCIMSSample" <- function(x, i, j, ...) {
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


#' @describeIn GCIMSSample-methods Subset a [GCIMSSample-class] object
#'
#' @param x A GCIMSSample object
#' @param dt_idx A vector of drift time indices to keep
#' @param rt_idx A vector of retention time indices to keep
#'
#' @return `subset`: A subsetted `GCIMSSample` object
#'
#' @aliases subset
#' @seealso [base::subset()]
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
