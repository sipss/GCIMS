#' Get the peak list
#'
#' @importMethodsFrom ProtGenerics peaks
#' @param object A [GCIMSDataset] object
#' @return A data frame with the detected peaks
#' @export
setMethod(
  "peaks",
  "GCIMSDataset",
  function(object) {
    object$realize()
    if (is.null(object$peaks)) {
      cli_abort("Please run {.code findPeaks()} on your dataset first")
    }
    x <- tibble::as_tibble(object$peaks)
    x$UniqueID <- sprintf("%s/%s", x$SampleID, x$PeakID)
    x
  }
)

#' @describeIn peaks-GCIMSDataset-method Set the peak list
#'
#' @param object A [GCIMSDataset] object
#' @param value The data frame with a peak list
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSDataset",
  function(object, value) {
    object$peaks <- methods::as(value, "DataFrame")
    object
  }
)
