# S3 classes can be used as S4 slots registering them as "old" classes
# (old because S4 is newer than S3, not because S3 is deprecated, not great naming)
methods::setOldClass("numeric_version")


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
    peaks = "DataFrameOrNULL"
  )
)

.CURRENT_GCIMSSAMPLE_CLASS_VERSION <- numeric_version("0.0.3")

methods::setMethod(
  "initialize", "GCIMSSample",
  function(.Object, drift_time, retention_time, data, ...) {
    .Object@class_version <- .CURRENT_GCIMSSAMPLE_CLASS_VERSION
    .Object@drift_time <- drift_time
    .Object@retention_time <- retention_time
    .Object@data <- data
    .Object@peaks <- NULL
    dots <- list(...)
    if (length(dots) == 0) {
      return(.Object)
    }
    if (is.null(names(dots)) || any(nchar(names(dots)) == 0) ) {
      rlang::abort(c("Error creating GCIMSSample object",
                     "x" = "All arguments should be named"))
    }
    # class_version and proc_params are internal and should not be given
    # drift_time, retention_time, data are mandatory and already given outside of ...
    invalid_dot_names <-  c("drift_time", "retention_time", "data", "class_version", "proc_params")
    valid_dot_names <- setdiff(methods::slotNames("GCIMSSample"), invalid_dot_names)
    if (!all(names(dots) %in% valid_dot_names)) {
      wrong_dot_names <- setdiff(names(dots), valid_dot_names)
      names(wrong_dot_names) <- rep("x", length(wrong_dot_names))
      rlang::abort(c("Invalid named arguments in GCIMSSample initialization", wrong_dot_names))
    }
    for (arg in names(dots)) {
      methods::slot(.Object, arg) <- dots[[arg]]
    }
    .Object
  })


#' @describeIn GCIMSSample Create a GCIMSSample object
#'
#' @param ... See the slots section in this page
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
  ...
) {
  methods::new("GCIMSSample", ...)
}



#' @name GCIMSSample-methods
#' @title Methods for the GCIMSSample class
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
#' @param dt_idx A vector of drift time indices to keep
#' @param rt_idx A vector of retention time indices to keep
#' @param dt_range The drift time range to keep (in milliseconds)
#' @param rt_range The retention time range to keep (in seconds)
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
  dt_range = NULL,
  rt_range = NULL,
  ...
) {

  dt <- dtime(x)
  rt <- rtime(x)

  if (!is.null(dt_range)) {
    dt_idx <- which(dt >= min(dt_range) & dt <= max(dt_range))
  }

  if (!is.null(rt_range)) {
    rt_idx <- which(rt >= min(rt_range) & rt <= max(rt_range))
  }

  if (is.null(dt_idx)) {
    dt_idx <- seq_along(dt)
  }
  if (is.null(rt_idx)) {
    rt_idx <- seq_along(rt)
  }

  if (identical(dt_idx, seq_along(dt)) &&
      identical(rt_idx, seq_along(rt))) {
    return(x)
  }

  new_obj <- x
  new_obj@drift_time <- dt[dt_idx]
  new_obj@retention_time <- rt[rt_idx]
  new_obj@data <- x@data[dt_idx, rt_idx, drop = FALSE]
  new_obj
}


#' @describeIn GCIMSSample-methods Get the extracted ion chromatogram
#' @export
setMethod("getEIC", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
  dt <- dtime(object)
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
  intens <- intensity(object, idx)
  GCIMSChromatogram(
    retention_time = rt[idx[["rt_logical"]]],
    intensity = colSums(intens),
    drift_time_idx = unique(c(idx[["dt_idx_min"]], idx[["dt_idx_max"]])),
    drift_time_ms = unique(c(idx[["dt_ms_min"]], idx[["dt_ms_max"]])),
    description = object@description
  )
})

#' @describeIn GCIMSSample-methods Get IMS spectrum
#' @param object A GCIMSSample object
#' @export
setMethod("getIMS", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
  dt <- dtime(object)
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
  intens <- intensity(object, idx)
  GCIMSSpectrum(
    drift_time = dt[idx[["dt_logical"]]],
    intensity = rowSums(intens),
    retention_time_idx = unique(c(idx[["rt_idx_min"]], idx[["rt_idx_max"]])),
    retention_time_s = unique(c(idx[["rt_s_min"]], idx[["rt_s_max"]])),
    description = object@description
  )
})


dt_rt_range_normalization <- function(dt = numeric(0L), rt = numeric(0L), dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
  out <- list(
    dt_ms_min = NA_real_,
    dt_ms_max = NA_real_,
    dt_idx_min = NA_integer_,
    dt_idx_max = NA_integer_,
    dt_logical = rep(NA, length(dt)),
    rt_s_min = NA_real_,
    rt_s_max = NA_real_,
    rt_idx_min = NA_integer_,
    rt_idx_max = NA_integer_,
    rt_logical = rep(NA, length(rt))
  )
  if (!is.null(dt_range) && !is.null(dt_idx)) {
    rlang::abort("Please provide either dt_range or dt_idx, but not both in the same call")
  }

  if (!is.null(dt_range)) {
    if (length(dt_range) == 2) {
      out[["dt_ms_min"]] <- min(dt_range)
      out[["dt_ms_max"]] <- max(dt_range)
      out[["dt_logical"]] <- dt >= out[["dt_ms_min"]] & dt <= out[["dt_ms_max"]]
    } else if (length(dt_range) == 1) {
      dt_distance <- abs(dt - dt_range)
      dt_idx_closest <- which.min(dt_distance)
      dt_ms_distance <- dt_distance[dt_idx_closest]
      if (dt_ms_distance > (dt[2] - dt[1])) {
        rlang::abort(glue("The given dt_range {dt_range} is not in [{dt[1]} - {dt[length(dt)]}]"))
      }
      out[["dt_ms_min"]] <- dt_range
      out[["dt_ms_max"]] <- dt_range
      out[["dt_logical"]] <- rep(FALSE, length(dt))
      out[["dt_logical"]][dt_idx_closest] <- TRUE
    } else {
      stop("dt_range should be of length 2 or a single number")
    }
    dt_idx <- which(out[["dt_logical"]])
    if (length(dt_idx) > 0) {
      out[["dt_idx_min"]] <- dt_idx[1]
      out[["dt_idx_max"]] <- dt_idx[length(dt_idx)]
    }
  } else if (!is.null(dt_idx)) {
    if (is.logical(dt_idx)) {
      if (length(dt_idx) != length(dt)) {
        rlang::abort(
          sprintf(
            "dt_idx is a logical vector of length %d and it should be of length %d",
            length(dt_idx),
            length(dt)
          )
        )
      }
      out[["dt_logical"]] <- dt_idx
      dt_idx <- which(out[["dt_logical"]])
      if (length(dt_idx) > 0) {
        out[["dt_idx_min"]] <- dt_idx[1]
        out[["dt_idx_max"]] <- dt_idx[length(dt_idx)]
      }
      out[["dt_ms_min"]] <- dt[out[["dt_idx_min"]]]
      out[["dt_ms_max"]] <- dt[out[["dt_idx_max"]]]
    } else if (is.numeric(dt_idx)) {
      out[["dt_idx_min"]] <- min(dt_idx)
      out[["dt_idx_max"]] <- max(dt_idx)
      if (out[["dt_idx_min"]] < 1 || out[["dt_idx_max"]] > length(dt)) {
        rlang::abort("dt_idx out of range")
      }
      out[["dt_ms_min"]] <- dt[out[["dt_idx_min"]]]
      out[["dt_ms_max"]] <- dt[out[["dt_idx_max"]]]
      out[["dt_logical"]] <- dt >= out[["dt_ms_min"]] & dt <= out[["dt_ms_max"]]
    }
  } else {
    if (length(dt) > 0) {
      out[["dt_ms_min"]] <- dt[1L]
      out[["dt_ms_max"]] <- dt[length(dt)]
      out[["dt_logical"]] <- rep(TRUE, length(dt))
      out[["dt_idx_min"]] <- 1L
      out[["dt_idx_max"]] <- length(dt)
    } else {
      out[["dt_ms_min"]] <- integer(0L)
      out[["dt_ms_max"]] <- integer(0L)
      out[["dt_logical"]] <- logical(0L)
      out[["dt_idx_min"]] <- integer(0L)
      out[["dt_idx_max"]] <- integer(0L)
    }
  }

  if (!is.null(rt_range) && !is.null(rt_idx)) {
    rlang::abort("Please provide either rt_range or rt_idx, but not both in the same call")
  }

  if (!is.null(rt_range)) {
    if (length(rt_range) == 2) {
      out[["rt_s_min"]] <- min(rt_range)
      out[["rt_s_max"]] <- max(rt_range)
      out[["rt_logical"]] <- rt >= out[["rt_s_min"]] & rt <= out[["rt_s_max"]]
    } else if (length(rt_range) == 1) {
      rt_distance <- abs(rt - rt_range)
      rt_idx_closest <- which.min(rt_distance)
      rt_s_distance <- rt_distance[rt_idx_closest]
      if (rt_s_distance > (rt[2] - rt[1])) {
        rlang::abort(glue("The given rt_range {rt_range} is not in [{rt[1L]} - {rt[length(rt)]}]"))
      }
      out[["rt_s_min"]] <- rt_range
      out[["rt_s_max"]] <- rt_range
      out[["rt_logical"]] <- rep(FALSE, length(rt))
      out[["rt_logical"]][rt_idx_closest] <- TRUE
    } else {
      stop("rt_range should be of length 2 or a single number")
    }

    rt_idx <- which(out[["rt_logical"]])
    if (length(rt_idx) > 0) {
      out[["rt_idx_min"]] <- rt_idx[1]
      out[["rt_idx_max"]] <- rt_idx[length(rt_idx)]
    }
  } else if (!is.null(rt_idx)) {
    if (is.logical(rt_idx)) {
      if (length(rt_idx) != length(rt)) {
        rlang::abort(
          sprintf(
            "rt_idx is a logical vector of length %d and it should be of length %d",
            length(rt_idx),
            length(rt)
          )
        )
      }
      out[["rt_logical"]] <- rt_idx
      rt_idx <- which(out[["rt_logical"]])
      if (length(rt_idx) > 0) {
        out[["rt_idx_min"]] <- rt_idx[1]
        out[["rt_idx_max"]] <- rt_idx[length(rt_idx)]
      }
      out[["rt_s_min"]] <- rt[out[["rt_idx_min"]]]
      out[["rt_s_max"]] <- rt[out[["rt_idx_max"]]]
    } else if (is.numeric(rt_idx)) {
      out[["rt_idx_min"]] <- min(rt_idx)
      out[["rt_idx_max"]] <- max(rt_idx)
      if (out[["rt_idx_min"]] < 1 || out[["rt_idx_max"]] > length(rt)) {
        rlang::abort("rt_idx out of range")
      }
      out[["rt_s_min"]] <- rt[out[["rt_idx_min"]]]
      out[["rt_s_max"]] <- rt[out[["rt_idx_max"]]]
      out[["rt_logical"]] <- rt >= out[["rt_s_min"]] & rt <= out[["rt_s_max"]]
    }
  } else {
    if (length(rt) > 0) {
      out[["rt_s_min"]] <- rt[1L]
      out[["rt_s_max"]] <- rt[length(rt)]
      out[["rt_logical"]] <- rep(TRUE, length(rt))
      out[["rt_idx_min"]] <- 1L
      out[["rt_idx_max"]] <- length(rt)
    } else {
      out[["rt_s_min"]] <- integer(0L)
      out[["rt_s_max"]] <- integer(0L)
      out[["rt_logical"]] <- logical(0L)
      out[["rt_idx_min"]] <- integer(0L)
      out[["rt_idx_max"]] <- integer(0L)
    }
  }
  class(out) <- "dt_rt_range_normalization"
  out
}

