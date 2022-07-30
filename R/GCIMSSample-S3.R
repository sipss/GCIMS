#' Create a GCIMSSample (S3 version)
#'
#' This is mostly compatible with our previous list-based implementation.
#'
#' @param drift_time A numeric vector
#' @param retention_time A numeric vector
#' @param intensity A matrix, with size length(drift_time) x length(retention_time)
#' @param metadata A named list of name-value annotations
#'
#' @return A GCIMSSampleS3 object
#' @export
#'
#' @examples
#' obj <- new_GCIMSSampleS3(
#'   drift_time = 1:3,
#'   retention_time = 1:5,
#'   intensity = matrix(1:15, nrow = 3, ncol = 5),
#'   metadata = list(Name = "Sample1")
#' )
new_GCIMSSampleS3 <- function(drift_time, retention_time, intensity, metadata) {
  stopifnot(is.matrix(intensity))
  stopifnot(length(drift_time) == nrow(intensity))
  stopifnot(length(retention_time) == ncol(intensity))
  stopifnot(is.list(metadata))

  object <- list(
    metadata = metadata,
    data = list(
      drift_time = drift_time,
      retention_time = retention_time,
      data_df = intensity
    )
  )
  class(object) <- "GCIMSSample_S3"
  object
}

# So we can use S4 generic methods with this S3 class:
methods::setOldClass("GCIMSSample_S3")


#' Get the drift time vector of a GCIMSSample_S3 object
#' @param object The GCIMSSample_S3 object
#' @export
#' @examples
#' obj <- new_GCIMSSampleS3(
#'   drift_time = 1:3,
#'   retention_time = 1:5,
#'   intensity = matrix(1:15, nrow = 3, ncol = 5),
#'   metadata = list(Name = "Sample1")
#' )
#' all(dtime(obj) == 1:3)
setMethod("dtime", "GCIMSSample_S3", function(object) unclass(object)$data$drift_time)


#' Get the retention time vector of a GCIMSSample_S3 object
#' @importMethodsFrom ProtGenerics rtime
#' @param object The GCIMSSample_S3 object
#' @export
#' @examples
#' obj <- new_GCIMSSampleS3(
#'   drift_time = 1:3,
#'   retention_time = 1:5,
#'   intensity = matrix(1:15, nrow = 3, ncol = 5),
#'   metadata = list(Name = "Sample1")
#' )
#' all(rtime(obj) == 1:5)
setMethod("rtime", "GCIMSSample_S3", function(object) unclass(object)$data$retention_time)

#' Get the intensity matrix of a GCIMSSample_S3 object
#' @importMethodsFrom ProtGenerics intensity
#' @param object The GCIMSSample_S3 object
#' @param dt_range The minimum and maximum drift times to extract (length 2 vector)
#' @param rt_range The minimum and maximum retention times to extract (length 2 vector)
#' @param dt_idx A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)
#' @param rt_idx A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)
#' @export
#' @examples
#' m <- matrix(1:15, nrow = 3, ncol = 5)
#' obj <- new_GCIMSSampleS3(
#'   drift_time = 1:3,
#'   retention_time = 1:5,
#'   intensity = m,
#'   metadata = list(Name = "Sample1")
#' )
#' intensity(obj)
setMethod(
  "intensity",
  "GCIMSSample_S3",
  function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
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
    out <- object$data$data_df[dt_idx, rt_idx]
    dimnames(out) <- list(dt_ms = dt[dt_idx], rt_s = rt[rt_idx])
    out
  }
)
