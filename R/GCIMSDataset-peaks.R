#' Peak detection on the GCIMS dataset
#' @param object A [GCIMSDataset] object, modified in-place
#' @param noise_level     Scalar number. The number of times the standard deviation
#'   above the noise level needed to detect a peak.
#' @param verbose If `TRUE`, debug information will be printed
#' @return The modified [GCIMSDataset], with a peak list
#' @importMethodsFrom xcms findPeaks
#' @export
setMethod(
  "findPeaks",
  "GCIMSDataset",
  function(object, noise_level = 3, verbose = FALSE) {
    delayed_op <- GCIMSDelayedOp(
      name = "findPeaks",
      fun = findPeaks,
      params = list(noise_level = noise_level, verbose = verbose),
      fun_extract = peaks,
      fun_aggregate = function(ds, objs) {
        p <- dplyr::bind_rows(purrr::map(objs, as.data.frame), .id = "SampleID")
        peaks(ds) <- p
        ds
      }
    )
    object <- appendDelayedOp(object, delayed_op)
    invisible(object)
  }
)


#' @importMethodsFrom ProtGenerics peaks
#' @export
setMethod(
  "peaks",
  "GCIMSDataset",
  function(object) {
    object <- realize(object)
    if (!"peaks" %in% names(object@envir)) {
      stop("Please run findPeaks() on your dataset first")
    }
    object@envir$peaks
  }
)

#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSDataset",
  function(object, value) {
    object@envir$peaks <- methods::as(value, "DataFrame")
    object
  }
)
