#' Filter GCIMSSample samples by retention time
#'
#' @param object A [GCIMSSample] object
#' @param rt_range A numeric vector of length 2 with the retention time range to keep, in seconds
#' @export
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' s <- filterRt(s, rt_range = c(5, 50))
#' @importMethodsFrom ProtGenerics filterRt
setMethod(
  "filterRt",
  "GCIMSSample",
  function(object, rt_range) {
    subset(object, rt_range = rt_range)
  }
)

#' Filter GCIMSSample samples by drift time
#'
#' @param object A [GCIMSSample] object
#' @param dt_range A numeric vector of length 2 with the drift time range to keep, in milliseconds
#' @export
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' s <- filterDt(s, dt_range = c(5, 9.5))
setMethod(
  "filterDt",
  "GCIMSSample",
  function(object, dt_range) {
    subset(object, dt_range = dt_range)
  }
)
