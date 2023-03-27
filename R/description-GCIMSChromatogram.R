#' @describeIn GCIMSChromatogram-class Get the description
#' @family GCIMSChromatogram
#' @param object A [GCIMSChromatogram] object
#' @export
methods::setMethod(
  "description", "GCIMSChromatogram",
  function(object) {
    object@description
  }
)

#' @describeIn GCIMSChromatogram-class Set the description
#' @param object A GCIMSChromatogram object
#' @param value A string with the description
#' @family GCIMSChromatogram
#' @export
methods::setMethod(
  "description<-", "GCIMSChromatogram",
  function(object, value) {
    if (!rlang::is_string(value)) {
      cli_abort("description should be a string")
    }
    object@description <- value
    object
  }
)
