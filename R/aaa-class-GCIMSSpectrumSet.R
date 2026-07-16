#' GCIMSSpectrumSet class
#'
#' @description
#' GCIMSSpectrumSet is an S4 class to store one [GCIMSSpectrum] per sample of
#' a [GCIMSDataset], together with a copy of `pData()` so plots can use the
#' dataset's annotations.
#'
#' Samples are not required to share a common drift time axis: each spectrum
#' keeps its own, exactly as extracted from its sample. No interpolation is
#' performed.
#'
#' @slot spectra A named list of [GCIMSSpectrum] objects, one per sample,
#' named after their `SampleID`.
#' @slot pData A `DataFrame` with the phenotype data, or `NULL`.
#'
#' @export
#' @family GCIMSSpectrum
methods::setClass(
  Class = "GCIMSSpectrumSet",
  slots = c(
    spectra = "list",
    pData = "DataFrameOrNULL"
  )
)

methods::setMethod(
  "initialize", "GCIMSSpectrumSet",
  function(.Object, spectra = list(), pData = NULL) {
    if (!rlang::is_named(spectra) && length(spectra) > 0) {
      cli_abort("spectra should be a named list, with the SampleID of each spectrum as its name")
    }
    if (!all(purrr::map_lgl(spectra, inherits, "GCIMSSpectrum"))) {
      cli_abort("All elements of spectra should be GCIMSSpectrum objects")
    }
    if (!is.null(pData)) {
      if (!"SampleID" %in% colnames(pData)) {
        cli_abort("pData should have a SampleID column")
      }
      if (!setequal(as.character(pData[["SampleID"]]), names(spectra))) {
        cli_abort(
          c(
            "pData$SampleID does not match the names of spectra",
            "i" = "Both should refer to exactly the same set of samples"
          )
        )
      }
    }
    .Object@spectra <- spectra
    .Object@pData <- pData
    .Object
  }
)

#' Create a [GCIMSSpectrumSet-class] object
#'
#' @param spectra A named list of [GCIMSSpectrum] objects, one per sample,
#' named after their `SampleID`.
#' @param pData A `data.frame`/`DataFrame`/tibble with the phenotype data, or `NULL`.
#' @return A [GCIMSSpectrumSet-class] object
#' @export
#' @family GCIMSSpectrum
GCIMSSpectrumSet <- function(spectra = list(), pData = NULL) {
  if (!is.null(pData) && !inherits(pData, "DataFrame")) {
    pData <- S4Vectors::DataFrame(pData)
  }
  methods::new("GCIMSSpectrumSet", spectra = spectra, pData = pData)
}

#' @describeIn GCIMSSpectrumSet-class Get the sample names
#' @param object A [GCIMSSpectrumSet] object
#' @return A character vector with the sample names
#' @export
setMethod("sampleNames", "GCIMSSpectrumSet", function(object) {
  nms <- names(object@spectra)
  if (is.null(nms)) character(0) else nms
})

#' @describeIn GCIMSSpectrumSet-class Get the phenotype data
#' @param object A [GCIMSSpectrumSet] object
#' @return A tibble with the phenotype data, or `NULL` if not set
#' @export
setMethod("pData", "GCIMSSpectrumSet", function(object) {
  if (is.null(object@pData)) {
    return(NULL)
  }
  tibble::as_tibble(object@pData)
})

#' @describeIn GCIMSSpectrumSet-class Number of spectra (samples) in the set
#' @param x A [GCIMSSpectrumSet] object
#' @return An integer with the number of spectra
#' @export
setMethod("length", "GCIMSSpectrumSet", function(x) length(x@spectra))

#' @describeIn GCIMSSpectrumSet-class Extract the spectrum of a single sample
#' @param x A [GCIMSSpectrumSet] object
#' @param i A number or a string with the sample index or name
#' @return The [GCIMSSpectrum] of the requested sample
#' @export
setMethod("[[", "GCIMSSpectrumSet", function(x, i) x@spectra[[i]])
