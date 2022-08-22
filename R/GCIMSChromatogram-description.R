#' @describeIn GCIMSChromatogram-class Get the description
#' @export
methods::setMethod(
  "description", "GCIMSChromatogram",
  function(object) {
    object@description
  }
)

#' @describeIn GCIMSChromatogram-class Set the description
#' @param value A string with the description
#' @export
methods::setMethod(
  "description<-", "GCIMSChromatogram",
  function(object, value) {
    if (!rlang::is_string(value)) {
      rlang::abort("description should be a string")
    }
    object@description <- value
  }
)
