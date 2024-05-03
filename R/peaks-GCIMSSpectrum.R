#' @describeIn GCIMSSpectrum-class Get the peak list
#' @importMethodsFrom ProtGenerics peaks
#' @return A data frame with the peaks in the spectrum
#' @export
setMethod(
  "peaks",
  "GCIMSSpectrum",
  function(object) {
    if (is.null(object@peaks)) {
      stop("Please run findPeaks() on the spectrum first")
    }
    if(nrow(object@peaks)==0){
      p <- NULL
    }else{
      p <- tibble::as_tibble(cbind(SampleID = object@description, object@peaks))
    }
    return(p)
  }
)

#' @describeIn GCIMSSpectrum-class Set the peak list
#' @param value A data frame with the peak list
#' @return The GCIMSSpectrum object
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSSpectrum",
  function(object, value) {
    value <- methods::as(value, "DataFrame")
    # We don't want this column:
    value$SampleID <- NULL
    object@peaks <- value
    object
  }
)
