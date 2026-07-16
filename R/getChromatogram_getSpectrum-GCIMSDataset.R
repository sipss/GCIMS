#' Get a chromatogram from each sample of a dataset
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @param aggregate Function that takes the subsetted intensity matrix of each
#' sample according to the region of interest and aggregates the drift times,
#' returning a vector representing the chromatogram intensity. `colSums` by
#' default.
#' @return A [GCIMSChromatogramSet], with one [GCIMSChromatogram] per sample
#' (each on its own retention time axis, no interpolation across samples) and
#' a copy of `pData(object)`
#' @export
setMethod(
  "getChromatogram",
  "GCIMSDataset",
  function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, aggregate = colSums) {
    object$realize()
    sample_names <- sampleNames(object)
    chromatograms <- stats::setNames(
      purrr::map(sample_names, function(sample_name) {
        sample <- object$getSample(sample_name)
        getChromatogram(
          sample,
          dt_range = dt_range, rt_range = rt_range,
          dt_idx = dt_idx, rt_idx = rt_idx,
          aggregate = aggregate
        )
      }),
      sample_names
    )
    GCIMSChromatogramSet(chromatograms = chromatograms, pData = pData(object))
  }
)

#' Get a spectrum from each sample of a dataset
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @param aggregate Function that takes the subsetted intensity matrix of each
#' sample according to the region of interest and aggregates the retention
#' times, returning a vector representing the spectrum intensity. `rowSums`
#' by default.
#' @return A [GCIMSSpectrumSet], with one [GCIMSSpectrum] per sample (each on
#' its own drift time axis, no interpolation across samples) and a copy of
#' `pData(object)`
#' @export
setMethod(
  "getSpectrum",
  "GCIMSDataset",
  function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, aggregate = rowSums) {
    object$realize()
    sample_names <- sampleNames(object)
    spectra <- stats::setNames(
      purrr::map(sample_names, function(sample_name) {
        sample <- object$getSample(sample_name)
        getSpectrum(
          sample,
          dt_range = dt_range, rt_range = rt_range,
          dt_idx = dt_idx, rt_idx = rt_idx,
          aggregate = aggregate
        )
      }),
      sample_names
    )
    GCIMSSpectrumSet(spectra = spectra, pData = pData(object))
  }
)
