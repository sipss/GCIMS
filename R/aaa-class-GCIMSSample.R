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
#' @slot history character. A character vector with a summary of information of the
#' processing details the sample has gone through already
#' @slot filepath character. A string with the path to the raw data
#' @slot description A string (optional). A sample name or ID or description used in plots
#' @slot proc_params list (internal). Data processing parameters computed and used internally.
#' @slot peaks A data frame (internal). The peak list, typically set using [findPeaks()].
#'  Use [peaks()] to get/set this.
#' @slot peaks_debug_info A list with arbitrary debug information from [findPeaks()].
#' @slot baseline A matrix of the same dimensions as `data` with the baseline. Use [estimateBaseline()] to estimate it
#' and `baseline()` to get or set it.
#' @slot class_version "numeric_version" (internal) The GCIMSSample object defines
#' internally a class version, so if a GCIMSSample object is saved, the GCIMS
#' package is updated and the GCIMSSample class has changed during the upgrade
#' it may be possible to upgrade the previously saved object when it's loaded.
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
    history = "character",
    filepath = "character",
    description = "character",
    class_version = "numeric_version",
    params = "list", # arbitrary parameters from the instrument
    proc_params = "list",
    peaks = "DataFrameOrNULL",
    peaks_debug_info = "listOrNULL", # arbitrary debug info from findPeaks()
    baseline = "matrixOrNULL"
  )
)

.CURRENT_GCIMSSAMPLE_CLASS_VERSION <- numeric_version("0.0.4")

methods::setMethod(
  "initialize", "GCIMSSample",
  function(.Object, drift_time, retention_time, data, ...) {
    # class_version and proc_params are internal and should not be given
    # drift_time, retention_time, data are mandatory and already given outside of ...
    invalid_dot_names <-  c(
      "drift_time", "retention_time", "data",
      "class_version", "proc_params", "baseline"
    )
    valid_dot_names <- setdiff(methods::slotNames("GCIMSSample"), invalid_dot_names)
    names(valid_dot_names) <- rep("*", length(valid_dot_names))
    dots <- list(...)
    if (length(dots) > 0) {
      if (is.null(names(dots)) || any(nchar(names(dots)) == 0) ) {
        cli_abort(
          c(
            "Error creating GCIMSSample object",
            "x" = "All extra arguments should be named",
            "Valid extra arguments are:",
            valid_dot_names
          )
        )
      }

      if (!all(names(dots) %in% valid_dot_names)) {
        wrong_dot_names <- setdiff(names(dots), c(valid_dot_names))
        names(wrong_dot_names) <- rep("x", length(wrong_dot_names))
        cli_abort(
          c(
            "Invalid named arguments in GCIMSSample()",
            wrong_dot_names,
            "Valid options are",
            "*" = "drift_time",
            "*" = "retention_time",
            "*" = "data",
            valid_dot_names
          )
        )
      }
    }

    .Object@class_version <- .CURRENT_GCIMSSAMPLE_CLASS_VERSION
    .Object@drift_time <- drift_time
    .Object@retention_time <- retention_time
    .Object@data <- data
    .Object@peaks <- NULL
    .Object@peaks_debug_info <- list()
    if (length(dots) == 0) {
      return(.Object)
    }

    for (arg in names(dots)) {
      methods::slot(.Object, arg) <- dots[[arg]]
    }
    .Object
  })


#' Create a [GCIMSSample-class] object
#'
#' @param drift_time,retention_time,data,... See the Slots section in [GCIMSSample-class] page
#' @export
#' @return A GCIMSSample object
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
  drift_time, retention_time, data, ...
) {
  methods::new("GCIMSSample", drift_time, retention_time, data, ...)
}



#' @name GCIMSSample-methods
#' @title Methods for the GCIMSSample class
#' @param object A [GCIMSSample] object
NULL



setValidity("GCIMSSample", function(object) {
  success <- TRUE
  issues <- c()
  dt <- dtime(object)
  rt <- rtime(object)
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

  if (!is.null(object@baseline) && any(dim(object@baseline) != dim(object@data))) {
    issues <- c(issues, "baseline should be of the same dimensions as data")
    success <- FALSE
  }

  if (!success) {
    return(issues)
  }
  success
})


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
  dt <- dtime(x)
  rt <- rtime(x)
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



#' @describeIn GCIMSSample-methods Dimension of the data matrix
#'
#' @param x A GCIMSSample object
#'
#' @return An integer vector with the number of rows and columns of the matrix
#' @export
#'
#' @examples
#' obj <- GCIMSSample(drift_time=1:2, retention_time=1:3, data = matrix(1:6, nrow=2, ncol=3))
#' dim(obj)
dim.GCIMSSample <- function(x) {
  return(dim(x@data))
}


#' @describeIn GCIMSSample-methods Subset a [GCIMSSample-class] object
#'
#' @param x A GCIMSSample object
#' @inheritParams dt_rt_range_normalization
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
  dt_range = NULL,
  rt_range = NULL,
  ...
) {

  dt <- dtime(x)
  rt <- rtime(x)
  idx <- dt_rt_range_normalization(
    dt = dt, rt = rt,
    dt_range = dt_range, rt_range = rt_range,
    dt_idx = dt_idx, rt_idx = rt_idx
  )

  new_obj <- x
  new_obj@drift_time <- dt[idx$dt_logical]
  new_obj@retention_time <- rt[idx$rt_logical]
  new_obj@data <- x@data[idx$dt_logical, idx$rt_logical, drop = FALSE]
  if (!is.null(new_obj@baseline)) {
    new_obj@baseline <- new_obj@baseline[idx$dt_logical, idx$rt_logical, drop = FALSE]
  }
  new_obj
}


#' Get the extracted ion chromatogram
#' @param object A [GCIMSSample] object
#' @inheritParams dt_rt_range_normalization
#' @param aggregate Function that takes the subsetted intensity matrix according to the region of interest
#'  and aggregates the drift times, returning a vector representing the chromatogram intensity. `colSums` by default.
#' @export
#' @return A [GCIMSChromatogram] object
#' @examples
#' x <- GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' getChromatogram(x)
#' # Take the maximum intensity in the region for each retention time:
#' sp1 <- getChromatogram(x, aggregate = function(x) apply(x, 2, max))
getChromatogram <- function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, aggregate = colSums) {
  dt <- dtime(object)
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
  int_mat_region <- intensity(object, idx)
  intens <- aggregate(int_mat_region)
  if (length(intens) != ncol(int_mat_region)) {
    cli_abort("The aggregation outcome should be a vector of length {ncol(int_mat_region)} but {length(intens)} was returned.")
  }
  basel <- baseline(object, idx, .error_if_missing = FALSE)
  if (!is.null(basel)) {
    basel <- aggregate(basel)
  }
  GCIMSChromatogram(
    retention_time = rt[idx[["rt_logical"]]],
    intensity = intens,
    drift_time_idx = unique(c(idx[["dt_idx_min"]], idx[["dt_idx_max"]])),
    drift_time_ms = unique(c(idx[["dt_ms_min"]], idx[["dt_ms_max"]])),
    description = object@description,
    baseline = basel
  )
}

#' Get IMS spectrum from a sample
#'
#' @inheritParams dt_rt_range_normalization
#' @param object A GCIMSSample object
#' @param aggregate Function that takes the subsetted intensity matrix according to the region of interest
#'  and aggregates the retention times, returning a vector representing the spectrum intensity. `rowSums` by default.
#' @export
#' @examples
#' x <- GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' getSpectrum(x, rt_idx = 2)
#'
#' # Take the maximum intensity in the region for each drift time:
#' sp1 <- getSpectrum(x, aggregate = function(x) apply(x, 1, max))
getSpectrum <- function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, aggregate = rowSums) {
  dt <- dtime(object)
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
  int_mat_region <- intensity(object, idx)
  intens <- aggregate(int_mat_region)
  if (length(intens) != nrow(int_mat_region)) {
    cli_abort("The aggregation outcome should be a vector of length {nrow(int_mat_region)} but {length(intens)} was returned.")
  }
  basel <- baseline(object, idx, .error_if_missing = FALSE)
  if (!is.null(basel)) {
    basel <- aggregate(basel)
  }
  GCIMSSpectrum(
    drift_time = dt[idx[["dt_logical"]]],
    intensity = intens,
    retention_time_idx = unique(c(idx[["rt_idx_min"]], idx[["rt_idx_max"]])),
    retention_time_s = unique(c(idx[["rt_s_min"]], idx[["rt_s_max"]])),
    description = object@description,
    baseline = basel
  )
}

