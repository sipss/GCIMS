#' Get/Set the phenotype data
#'
#' @importMethodsFrom Biobase pData
#' @param object A [GCIMSDataset] object
#' @return A data frame with the phenotype data
setMethod(
  "pData",
  "GCIMSDataset",
  function(object) {
    tibble::as_tibble(object$pData)
  }
)

#' @describeIn pData-GCIMSDataset-method Set pData
#' @param object A [GCIMSDataset] object
#' @param value The data frame with annotations, it should have a FileName column and a SampleID column.
#' @importMethodsFrom Biobase "pData<-"
setReplaceMethod(
  "pData",
  "GCIMSDataset",
  function(object, value) {
    value <- validate_pData(value)
    current_pData <- Biobase::pData(object)
    if (nrow(value) != nrow(current_pData)) {
      cli_abort(
        c(
          "Invalid replacement for pData",
          "x" = "The number of rows given ({nrow(value)}) != Number of rows expected ({nrow(current_pData)})"
        )
      )
    }
    if (any(!is.na(current_pData$FileName) & current_pData$FileName != value$FileName)) {
      cli_warn(
        c(
          "FileName changed in pData",
          "i" = "It is not advisable to replace the filenames, please create a new dataset instead if possible"
        )
      )
    }
    object$pData <- value
    object
  }
)
