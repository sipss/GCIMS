#' @describeIn GCIMSSample-methods Get the description
#' @export
methods::setMethod(
  "description", "GCIMSSample",
  function(object) {
    object@description
  }
)

#' @describeIn GCIMSSample-methods Set the description
#' @param value A string with the description
#' @export
methods::setMethod(
  "description<-", "GCIMSSample",
  function(object, value) {
    if (!rlang::is_string(value)) {
      rlang::abort("description should be a string")
    }
    object@description <- value
    object
  }
)
