#' @describeIn GCIMSChromatogram-class Get the description
#' @param object A [GCIMSChromatogram] object
#' @return The description of the chromatogram
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
#' @return The chromatogram object
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
