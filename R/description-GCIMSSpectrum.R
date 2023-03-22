#' @describeIn GCIMSSpectrum-class Get the description
#' @export
methods::setMethod(
  "description", "GCIMSSpectrum",
  function(object) {
    object@description
  }
)

#' @describeIn GCIMSSpectrum-class Set the description
#' @param value A string with the description
#' @export
methods::setMethod(
  "description<-", "GCIMSSpectrum",
  function(object, value) {
    if (!rlang::is_string(value)) {
      cli_abort("description should be a string")
    }
    object@description <- value
    object
  }
)
