#' Deprecated function. Use getChromatogram() instead
#' @param ... Arguments passed to [getChromatogram()]
#' @export
getEIC <- function(...) {
  rlang::warn(c(
    "Deprecation warning",
    c("i" = "getEIC() is being renamed to getChromatogram()."),
    c("i" = "Please use getChromatogram() instead")
  ),
  .frequency = "once",
  .frequency_id = "getEIC-deprecation"
  )
  getChromatogram(...)
}

#' Deprecated function. Use getSpectrum() instead
#' @param ... Arguments passed to [getSpectrum()]
#' @export
getIMS <- function(...) {
  rlang::warn(c(
    "Deprecation warning",
    c("i" = "getIMS() is being renamed to getSpectrum()."),
    c("i" = "Please use getSpectrum() instead")
  ),
  .frequency = "once",
  .frequency_id = "getSpectrum-deprecation"
  )
  getSpectrum(...)
}


#' Deprecated function. Use getSample() instead
#' @param ... Arguments passed to [getSample()]
#' @export
getGCIMSSample <- function(...) {
  rlang::warn(c(
    "Deprecation warning",
    c("i" = "getGCIMSSample() is being renamed to getSample()."),
    c("i" = "Please use getSample() instead")
  ),
  .frequency = "once",
  .frequency_id = "getGCIMSSample-deprecation"
  )
  getSample(...)
}
