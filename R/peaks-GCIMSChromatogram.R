
#' @describeIn GCIMSChromatogram-class Get the peak list
#' @importMethodsFrom ProtGenerics peaks
#' @return A data frame with the peaks in the chromatogram
#' @export
setMethod(
  "peaks",
  "GCIMSChromatogram",
  function(object) {
    if (is.null(object@peaks)) {
      stop("Please run findPeaks() on the chromatogram first")
    }
    tibble::as_tibble(cbind(SampleID = object@description, object@peaks))
  }
)

#' @describeIn GCIMSChromatogram-class Set the peak list
#' @param value A data frame with the peak list
#' @return The GCIMSChromatogram object
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSChromatogram",
  function(object, value) {
    value <- methods::as(value, "DataFrame")
    # We don't want this column:
    value$SampleID <- NULL
    object@peaks <- value
    object
  }
)
