#' GCIMSChromatogramSet class
#'
#' @description
#' GCIMSChromatogramSet is an S4 class to store one [GCIMSChromatogram] per
#' sample of a [GCIMSDataset], together with a copy of `pData()` so plots can
#' use the dataset's annotations.
#'
#' Samples are not required to share a common retention time axis: each
#' chromatogram keeps its own, exactly as extracted from its sample. No
#' interpolation is performed.
#'
#' @slot chromatograms A named list of [GCIMSChromatogram] objects, one per
#' sample, named after their `SampleID`.
#' @slot pData A `DataFrame` with the phenotype data, or `NULL`.
#'
#' @export
#' @family GCIMSChromatogram
methods::setClass(
  Class = "GCIMSChromatogramSet",
  slots = c(
    chromatograms = "list",
    pData = "DataFrameOrNULL"
  )
)

methods::setMethod(
  "initialize", "GCIMSChromatogramSet",
  function(.Object, chromatograms = list(), pData = NULL) {
    if (!rlang::is_named(chromatograms) && length(chromatograms) > 0) {
      cli_abort("chromatograms should be a named list, with the SampleID of each chromatogram as its name")
    }
    if (!all(purrr::map_lgl(chromatograms, inherits, "GCIMSChromatogram"))) {
      cli_abort("All elements of chromatograms should be GCIMSChromatogram objects")
    }
    .Object@chromatograms <- chromatograms
    .Object@pData <- pData
    .Object
  }
)

#' Create a [GCIMSChromatogramSet-class] object
#'
#' @param chromatograms A named list of [GCIMSChromatogram] objects, one per
#' sample, named after their `SampleID`.
#' @param pData A `data.frame`/`DataFrame`/tibble with the phenotype data, or `NULL`.
#' @return A [GCIMSChromatogramSet-class] object
#' @export
#' @family GCIMSChromatogram
GCIMSChromatogramSet <- function(chromatograms = list(), pData = NULL) {
  if (!is.null(pData) && !inherits(pData, "DataFrame")) {
    pData <- S4Vectors::DataFrame(pData)
  }
  methods::new("GCIMSChromatogramSet", chromatograms = chromatograms, pData = pData)
}

#' @describeIn GCIMSChromatogramSet-class Get the sample names
#' @param object A [GCIMSChromatogramSet] object
#' @return A character vector with the sample names
#' @export
setMethod("sampleNames", "GCIMSChromatogramSet", function(object) {
  nms <- names(object@chromatograms)
  if (is.null(nms)) character(0) else nms
})

#' @describeIn GCIMSChromatogramSet-class Get the phenotype data
#' @param object A [GCIMSChromatogramSet] object
#' @return A tibble with the phenotype data, or `NULL` if not set
#' @export
setMethod("pData", "GCIMSChromatogramSet", function(object) {
  if (is.null(object@pData)) {
    return(NULL)
  }
  tibble::as_tibble(object@pData)
})

#' @describeIn GCIMSChromatogramSet-class Number of chromatograms (samples) in the set
#' @param x A [GCIMSChromatogramSet] object
#' @return An integer with the number of chromatograms
#' @export
setMethod("length", "GCIMSChromatogramSet", function(x) length(x@chromatograms))

#' @describeIn GCIMSChromatogramSet-class Extract the chromatogram of a single sample
#' @param x A [GCIMSChromatogramSet] object
#' @param i A number or a string with the sample index or name
#' @return The [GCIMSChromatogram] of the requested sample
#' @export
setMethod("[[", "GCIMSChromatogramSet", function(x, i) x@chromatograms[[i]])
