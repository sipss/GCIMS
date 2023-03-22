
#' @describeIn GCIMSChromatogram-class Get the peak list
#' @importMethodsFrom ProtGenerics peaks
#' @family GCIMSChromatogram
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
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @family GCIMSChromatogram
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
