#' Get the total ion spectrum
#'
#' @param object A [GCIMSSample] object
#'
#' @return A numeric vector with the total ion spectrum
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' tis <- getTIS(s)
setMethod("getTIS", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  tis <- rowSums(intmat)
  tis
})

#' Get the reverse ion chromatogram
#'
#' @param object A [GCIMSSample] object
#'
#' @return A numeric vector with the reverse ion chromatogram
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' ric <- getRIC(s)
setMethod("getRIC", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  tis <- rowSums(intmat)
  ric_pos <- which.max(tis)
  ric <- intmat[ric_pos, ]
  ric <- max(ric) - ric
  ric <- ric/sum(ric)
  ric
})


