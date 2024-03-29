
#' @describeIn GCIMSSample-methods Get the peak list
#' @importMethodsFrom ProtGenerics peaks
#' @export
setMethod(
  "peaks",
  "GCIMSSample",
  function(object) {
    if (is.null(object@peaks)) {
      stop("Please run findPeaks() on the sample first")
    }
    if(nrow(object@peaks)==0){
      p <- NULL
    }else{
      p <- tibble::as_tibble(cbind(SampleID = object@description, object@peaks))
    }
    return(p)
  }
)

#' @describeIn GCIMSSample-methods Set the peak list
#' @param value A data frame with the peak list
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSSample",
  function(object, value) {
    value <- methods::as(value, "DataFrame")
    # We don't want this column:
    value$SampleID <- NULL
    object@peaks <- value
    object
  }
)
