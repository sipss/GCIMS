#' Peak detection for a GCIMSSample
#' @param object A [GCIMSSample] object
#' @param noise_level     Scalar number. The number of times the standard deviation
#'   above the noise level needed to detect a peak.
#' @param verbose If `TRUE`, debug information will be printed
#' @return The modified [GCIMSSample], with a peak list
#' @importMethodsFrom xcms findPeaks
#' @export
setMethod(
  "findPeaks",
  "GCIMSSample",
  function(object, noise_level = 3, verbose = FALSE) {
    dt <- dtime(object)
    rt <- rtime(object)
    int_mat <- intensity(object)
    peak_list <- peak_detection(
      drift_time = dt,
      retention_time = rt,
      int_mat = int_mat,
      noise_level = noise_level,
      verbose = verbose
    )
    peaks(object) <- peak_list
    object
  }
)


#' @importMethodsFrom ProtGenerics peaks
#' @export
setMethod(
  "peaks",
  "GCIMSSample",
  function(object) {
    if (is.null(object@peaks)) {
      stop("Please run findPeaks() on the sample first")
    }
    object@peaks
  }
)

#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSSample",
  function(object, value) {
    value <- S4Vectors::as.data.frame(value)
    object@peaks <- value
    object
  }
)
