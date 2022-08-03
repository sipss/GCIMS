# R packages, when built (using R CMD build) collate all R/*.R files either
# alphabetically (based on the C locale alphabet order) or in any arbitrary
# way defined in the DESCRIPTION "Collate:" field. The collation order defines
# the order of function definitions in the package.
#
# Generic methods for S4 classes must be declared with setGeneric before they
# are defined with setMethod.
#
# This file it is named AllGenerics.R so it is collated *first*, and it includes
# all the generic declarations in the package.

#' Get drift time vector
#'
#' @return A numeric vector with the drift time
#'
#' @param object An object with drift time
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#' @examples
#' x <- dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' dtime(x) # c(1,2)
#'
setGeneric("dtime", function(object, ...) standardGeneric("dtime"))

#' Get extracted ion chromatogram
#'
#' @return A numeric vector with the extracted ion chromatogram
#'
#' @param object An object that has ion chromatograms
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#' @examples
#' x <- dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' getEIC(x) # c(`1` = 3, `2` = 7, `3` = 11)
setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))


#' Get IMS scans
#'
#' @return An ion mobility spectrum
#' @param object An object that has IMS scans
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#' @examples
#' x <- dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' getIMSScan(x) # c(`1` = 9, `2` = 12)
setGeneric("getIMSScan", function(object, ...) standardGeneric("getIMSScan"))
#' Plot raw data
#'
#' @return A plot
#' @param object An object to plot raw data
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotRaw", function(object, ...) standardGeneric("plotRaw"))
