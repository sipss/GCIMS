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
    class_version = "numeric_version",
    params = "list" # arbitrary parameters from the instrument
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
    dots <- list(...)
    if (length(dots) == 0) {
      return(.Object)
    }
    if (is.null(names(dots)) || any(nchar(names(dots)) == 0) ) {
      rlang::abort(c("Error creating GCIMSSample object",
                     "x" = "All arguments should be named"))
    }
    # class_version is internal and should not be given
    # drift_time, retention_time, data are mandatory and already given outside of ...
    invalid_dot_names <-  c("drift_time", "retention_time", "data", "class_version")
    valid_dot_names <- setdiff(methods::slotNames("GCIMSSample"), invalid_dot_names)
    if (!all(names(dots) %in% valid_dot_names)) {
      wrong_dot_names <- setdiff(names(dots), valid_dot_names)
      # This naming is part of rlang formatting (see ?rlang::abort or
      # ?rlang::format_error_bullets)
      names(wrong_dot_names) <- "x"
      rlang::abort(c("Invalid named arguments in GCIMSSample initialization", wrong_dot_names))
    }
    for (arg in names(dots)) {
      methods::slot(.Object, arg) <- dots[[arg]]
    }
    .Object
  })

methods::setMethod("smooth","GCIMSSample",function(arr, p, n, m){
  dx = sgolayfilt(arr, p, n, m)
  return(dx)
})

#' Updates old saved GCIMSSample object to the latest version
#'
#' This function is useful when you have saved a [GCIMSSample] object
#' with a previous version of the GCIMS package and you want to load it
#' using a new version of the package.
#'
#' The function allows you to update the old object, adding missing
#' slots, etc so it is fully compatible with the new class definition.
#'
#'
#' @param object A [GCIMSSample] object
#'
#' @return The updated [GCIMSSample] object
#' @export
#'
#' @examples
#' obj <- GCIMSSample(drift_time=1:2, retention_time=1:3, data = matrix(1:6, nrow=2, ncol=3))
#' # Here we modify the object as if it was older:
#' attr(obj, 'history') <- NULL
#' attr(obj, "class_version") <- numeric_version("0.0.1")
#' # Update the object:
#' newobj <- UpdateGCIMSSample(obj)
#'
UpdateGCIMSSample <- function(object) {
  if (!"class_version" %in% methods::slotNames(object)) {
    rlang::abort("The object is too old and can't be updated")
  }
  obj_version <- methods::slot(object, "class_version")
  if (obj_version <= "0.0.1") {
    # Migrate objects from 0.0.1 to the latest version
    object@history <- "Unknown initial history"
    object@filepath <- ""
    object@class_version <- numeric_version("0.0.3")
    return(object)
  }
  if (obj_version == "0.0.2") {
    object@filepath <- ""
    object@class_version <- numeric_version("0.0.3")
    return(object)

  }
  return(object)
}


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



#' @describeIn GCIMSSample-methods Get the drift time vector
#' @export
setMethod("dtime", "GCIMSSample", function(object) object@drift_time)

#' @describeIn GCIMSSample-methods Get the retention time vector
#' @import ProtGenerics
#' @export
setMethod("rtime", "GCIMSSample", function(object) object@retention_time)

#' @describeIn GCIMSSample-methods Get the intensity matrix
#' @import ProtGenerics
#' @param dt_range The minimum and maximum drift times to extract (length 2 vector)
#' @param rt_range The minimum and maximum retention times to extract (length 2 vector)
#' @param dt_idx A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)
#' @param rt_idx A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)
#' @examples
#' mea_file <- system.file("extdata",  "sample_formats", "small.mea.gz", package = "GCIMS")
#' gcims_sample <- read_mea(mea_file)
#' my_matrix <- intensity(gcims_sample, dt_range = c(7, 8), rt_range = c(1,30))
setMethod("intensity", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
  dt <- dtime(object)
  if (!is.null(dt_range) && !is.null(dt_idx)) {
    rlang::abort("Please provide either dt_range or dt_idx, but not both in the same call")
  }
  if (!is.null(dt_range)) {
    dtmin <- min(dt_range)
    dtmax <- max(dt_range)
    dt_idx <- dt >= dtmin & dt < dtmax
  } else if (!is.null(dt_idx)) {
    if (is.logical(dt_idx) && length(dt_idx) != length(dt)) {
      rlang::abort(
        sprintf(
          "dt_idx is a logical vector of length %d and it should be of length %d",
          length(dt_idx),
          length(dt)
        )
      )
    } else if (is.numeric(dt_idx) && (min(dt_idx) < 1 || max(dt_idx) > length(dt))) {
      rlang::abort("dt_idx out of range")
    }
  } else {
    dt_idx <- seq_along(dt)
  }

  rt <- rtime(object)
  if (!is.null(rt_range) && !is.null(rt_idx)) {
    rlang::abort("Please provide either rt_range or rt_idx, but not both in the same call")
  }
  if (!is.null(rt_range)) {
    rtmin <- min(rt_range)
    rtmax <- max(rt_range)
    rt_idx <- rt >= rtmin & rt < rtmax
  } else if (!is.null(rt_idx)) {
    if (is.logical(rt_idx) && length(rt_idx) != length(rt)) {
      rlang::abort(
        sprintf(
          "rt_idx is a logical vector of length %d and it should be of length %d",
          length(rt_idx),
          length(rt)
        )
      )
    } else if (is.numeric(rt_idx) && (min(rt_idx) < 1 || max(rt_idx) > length(rt))) {
      rlang::abort("rt_idx out of range")
    }
  } else {
    rt_idx <- seq_along(rt)
  }
  out <- object@data[dt_idx, rt_idx]
  dimnames(out) <- list(dt_ms = dt[dt_idx], rt_s = rt[rt_idx])
  out
})

setMethod(
  "show",
  "GCIMSSample",
  function(object) {
    axes <- list(
      "drift time" = list(value=dtime(object), unit="ms"),
      "retention time" = list(value=rtime(object), unit="s")
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
    if (length(object@history) > 0) {
      outstring <- c(outstring, "History:", sprintf("- %s", object@history))
    } else {
      outstring <- c(outstring, "History:", "No history available")
    }
    # FIXME: Give more details about the sample
    cat(paste0(outstring, collapse = "\n"))
  }
)

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

#' @describeIn GCIMSSample Topographical plot of a GC-IMS Sample
#'
#' @return
#'
#' @param x A GCIMSSample object
#' @param ... passed to [graphics::filled.contour]
#' @export
#' @method image GCIMSSample
#'
#' @examples
#'
#' dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3),
#'   gc_column = "Optional column name",
#'   drift_gas = "nitrogen",
#'   drift_tube_length = 98.0 # in mm
#' )
#' image(dummy_obj)
"image.GCIMSSample" <- function(x, ...) {
  # filled.contour includes a legend
  dots <- list(...)
  if (!"main" %in% names(dots)) {
    dots[["main"]] <- basename(methods::slot(x, "filepath"))
  }
  rlang::exec(
    graphics::filled.contour,
    x=x@drift_time,
    y = x@retention_time/60.0,
    z=x@data,
    xlab = "Drift time (ms)",
    ylab = "Retention time (min)",
    !!!dots
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
    dt_idx <- which(dt >= min(dt_range) & dt < max(dt_range))
  }

  if (!is.null(rt_range)) {
    rt_idx <- which(rt >= min(rt_range) & rt < max(rt_range))
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
setMethod("getEIC", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL) {
  rlang::warn(
    message = c(
      "Deprecation warning",
      "x" = "getEIC will be deprecated",
      "i" = "We will provide a chromatogram() method, aligned with xcms, MSnbase and ProtGenerics"
    ),
    .frequency = "once",
    .frequency_id = "GCIMS::getEIC deprecation"
  )
  intens <- intensity(object, rt_range = rt_range, dt_range = dt_range)
  colSums(intens)
})

#' @describeIn GCIMSSample-methods Get IMS spectrum
#' @param object A GCIMSSample object
#' @export
setMethod("getIMSScan", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL) {
  intens <- intensity(object, rt_range = rt_range, dt_range = dt_range)
  rowSums(intens)
})

