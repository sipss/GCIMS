#' @describeIn GCIMSDataset-class Get the peak list
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
    x <- tibble::as_tibble(object@envir$peaks)
    x$UniqueID <- sprintf("%s/%s", x$SampleID, x$PeakID)
    x
  }
)

#' @describeIn GCIMSDataset-class Set the peak list
#' @param value The data frame with a peak list
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
