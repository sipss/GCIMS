#' @importMethodsFrom ProtGenerics filterRt
#' @export
setMethod(
  "filterRt",
  "GCIMSSample",
  function(object, rt_range) {
    subset(object, rt_range = rt_range)
  }
)

#' @export
setMethod(
  "filterDt",
  "GCIMSSample",
  function(object, dt_range) {
    subset(object, dt_range = dt_range)
  }
)
