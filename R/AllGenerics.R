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
#' @export
#'
setGeneric("driftTime", function(object) standardGeneric("driftTime"))

#' Get retention time vector
#'
#' @return A numeric vector with the retention time
#'
#' @param object An object with retention time
#' @export
#'
setGeneric("retentionTime", function(object) standardGeneric("retentionTime"))

#' Get the intensity matrix
#'
#' @return A numeric matrix with drift time in the rows and retention time in columns
#'
#' @param object An object with an intensity matrix
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#'
setGeneric("intensityMatrix", function(object, ...) standardGeneric("intensityMatrix"))


#' Get extracted ion chromatogram
#'
#' @return A numeric vector with the extracted ion chromatogram
#'
#' @param object An object that has ion chromatograms
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#'
setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))


#' Get IMS scans
#'
#' @return An ion mobility spectrum
#' @param object An object that has IMS scans
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getIMSScan", function(object, ...) standardGeneric("getIMSScan"))
