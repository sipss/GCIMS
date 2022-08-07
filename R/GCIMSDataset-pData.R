#' @importMethodsFrom Biobase pData
setMethod(
  "pData",
  "GCIMSDataset",
  function(object) {
    object@envir$pData
  }
)

#' @importMethodsFrom Biobase "pData<-"
setReplaceMethod(
  "pData",
  "GCIMSDataset",
  function(object, value) {
    value <- validate_pData(value)
    current_pData <- Biobase::pData(object)
    if (nrow(value) != nrow(current_pData)) {
      rlang::abort(
        message = c(
          "Invalid replacement for pData",
          "x" = glue("The number of rows of the given pData ({nrow(value)}) does not match the rows of current pData ({nrow(current_pData)})")
        )
      )
    }
    if (!all(current_pData$FileName == value$FileName)) {
      rlang::warn(
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
