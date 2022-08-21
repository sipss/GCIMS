#' @export
methods::setMethod(
  "description", "GCIMSSpectrum",
  function(object) {
    object@description
  }
)

#' @export
methods::setMethod(
  "description<-", "GCIMSSpectrum",
  function(object, value) {
    if (!rlang::is_string(value)) {
      rlang::abort("description should be a string")
    }
    object@description <- value
  }
)
