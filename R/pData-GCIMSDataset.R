#' @describeIn GCIMSDataset-class Get the phenotype data
#' @importMethodsFrom Biobase pData
setMethod(
  "pData",
  "GCIMSDataset",
  function(object) {
    object@envir$pData
  }
)

#' @describeIn GCIMSDataset-class Set the phenotype data
#' @param value The data frame with annotations, it should have a FileName column and a SampleID column.
#' @importMethodsFrom Biobase "pData<-"
setReplaceMethod(
  "pData",
  "GCIMSDataset",
  function(object, value) {
    value <- validate_pData(value)
    current_pData <- Biobase::pData(object)
    if (nrow(value) != nrow(current_pData)) {
      abort(
        message = c(
          "Invalid replacement for pData",
          "x" = glue("The number of rows of the given pData ({nrow(value)}) does not match the rows of current pData ({nrow(current_pData)})")
        )
      )
    }
    if (!all(current_pData$FileName == value$FileName)) {
      warn(
        message = c(
          "FileName changed in pData",
          "i" = "It is not advisable to replace the filenames, please create a new dataset instead if possible"
        )
      )
    }
    object@envir$pData <- value
    object
  }
)
