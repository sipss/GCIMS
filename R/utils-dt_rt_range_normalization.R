#' Index normalization
#'
#' Users may want to specify:
#'
#'  - A single retention/drift time to extract, in physical units or in indices
#'  - A range of retention/drift time to extract, in physical units or in indices
#'
#' This function is internal to the [GCIMS-package] and it is meant to be used
#' to normalize those specifications, returning an object that computes the rest
#' so the caller function can easily get the indices or physical units it needs.
#'
#' The behaviour is intuitive in corner cases, where the conversion from physical
#' units to indices may have rounding issues:
#' - If a single value is requested in physical units the closest one will be returned. If the
#'   requested single value is outside the valid range, then an error is given. A tolerance of
#'   one time step is accepted to prevent rounding issues.
#' - If a range value is requested in physical units, all the points within the range
#'   are included. The interval is considered closed at both ends, for consistency with how
#'   indexes behave in R.
#'
#' @param dt A numeric vector with the drift time
#' @param rt A numeric vector with the retention time
#' @param dt_range The minimum and maximum drift times to extract (length 2 vector)
#' @param rt_range The minimum and maximum retention times to extract (length 2 vector)
#' @param dt_idx A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)
#' @param rt_idx A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)
#' @return A named list with class `dt_rt_range_normalization`. The elements of the list are:
#' - `dt_ms_min`, `dt_ms_max`: Drift time in milliseconds (min and max). This is the same value if the user did want a single value
#' - `dt_idx_min`, `dt_idx_max`: The indices corresponding to the requested range limits in drift time
#' - `dt_logical`: A logical vector of length equal to the given drift time vector, with `TRUE` where the range is included and `FALSE` elsewhere. Useful for indexing.
#' - `rt_s_min`, `rt_s_max`: Retention time in seconds (min and max). This is the same value if the user did want a single value
#' - `rt_idx_min`, `rt_idx_max`: The indices corresponding to the requested range limits in retention time.
#' - `rt_logical`: A logical vector of length equal to the given retention time vector, with `TRUE` where the range is included and `FALSE` elsewhere. Useful for indexing.
#' @keywords internal
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
  class(out) <- "dt_rt_range_normalization"

  # Shortcut in case we already have computed the indices and passed them here:
  if (!is.null(dt_range) && inherits(dt_range, "dt_rt_range_normalization")) {
    return(dt_range)
  }
  if (!is.null(dt_range) && !is.null(dt_idx)) {
    cli_abort("Please provide either dt_range or dt_idx, but not both in the same call")
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
        cli_abort("The given dt_range {dt_range} is not in [{dt[1]} - {dt[length(dt)]}]")
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
        cli_abort("dt_idx is a logical vector of length {length(dt_idx)} and it should be of length {length(dt)}")
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
        cli_abort("dt_idx out of range")
      }
      out[["dt_ms_min"]] <- dt[out[["dt_idx_min"]]]
      out[["dt_ms_max"]] <- dt[out[["dt_idx_max"]]]
      # If dt_idx is numeric, it can include some of the indices only:
      out[["dt_logical"]] <- logical(length(dt))
      out[["dt_logical"]][dt_idx] <- TRUE
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
    cli_abort("Please provide either rt_range or rt_idx, but not both in the same call")
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
        cli_abort("The given rt_range {rt_range} is not in [{rt[1L]} - {rt[length(rt)]}]")
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
        cli_abort("rt_idx is a logical vector of length {length(rt_idx)} and it should be of length {length(rt)}")
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
        cli_abort("rt_idx out of range")
      }
      out[["rt_s_min"]] <- rt[out[["rt_idx_min"]]]
      out[["rt_s_max"]] <- rt[out[["rt_idx_max"]]]
      # If rt_idx is numeric, it can include some of the indices only:
      out[["rt_logical"]] <- logical(length(rt))
      out[["rt_logical"]][rt_idx] <- TRUE
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
  out
}
