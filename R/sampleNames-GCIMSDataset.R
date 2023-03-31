#' Sample names
#'
#' @param object GCIMSDataset object
#' @param value A character vector of length the number of samples with the sample names
#' @return The [GCIMSDataset] object
#' @importMethodsFrom Biobase sampleNames
#' @importMethodsFrom Biobase "sampleNames<-"
setMethod("sampleNames", "GCIMSDataset", function(object) object$sampleNames)

#' @describeIn sampleNames-GCIMSDataset-method Sample names
setReplaceMethod("sampleNames", "GCIMSDataset", function(object, value) {
  object$sampleNames <- value
  object
})


