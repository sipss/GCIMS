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

setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))


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
#' getIMS(x, rt_idx = 2)
setGeneric("getIMS", function(object, ...) standardGeneric("getIMS"))


#' Get a GCIMSSample
#'
#' @return A [GCIMSSample] object
#' @param object An object where we can extract GCIMSSamples from
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getGCIMSSample", function(object, ...) standardGeneric("getGCIMSSample"))

#' Get the Total Ion Spectrum
#'
#' @return The Total Ion Spectrum as a numeric vector or a matrix
#' (depending if the object is one sample or several)
#' @param object An object to extract its total ion spectrum
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getTIS", function(object, ...) standardGeneric("getTIS"))

#' Get the Reverse Ion Chromatogram
#'
#' @return The Reverse Ion Chromatogram, as a numeric vector or a
#' matrix (depending if the object is one sample or several)
#' @param object An object to extract its reverse ion chromatogram
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getRIC", function(object, ...) standardGeneric("getRIC"))


#' Plot raw data
#'
#' @return A plot
#' @param object An object to plot raw data
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotRaw", function(object, ...) standardGeneric("plotRaw"))

#' Plot total ion spectrum
#'
#' @return A plot
#' @param object An object to plot the total ion spectrum
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotTIS", function(object, ...) standardGeneric("plotTIS"))

#' Plot Reverse Ion Chromatogram
#'
#' @return A plot
#' @param object An object to plot its reverse ion chromatogram
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotRIC", function(object, ...) standardGeneric("plotRIC"))


#' Filter in Drift time
#' @return The object, modified
#' @param object An object to filter
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("filterDt", function(object, ...) standardGeneric("filterDt"))


#' This function returns a list of relevant object-related information to be
#' used when printing.
#'
#' @keywords internal
setGeneric("describeAsList", function(object, ...) standardGeneric("describeAsList"))

#' Decimate an object
#'
#' @param object An object to decimate
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("decimate", function(object, ...) standardGeneric("decimate"))


#' Align an object
#'
#' @param object An object to align
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("align", function(object, ...) standardGeneric("align"))


#' Align an object in drift time
#'
#' @param x An object to align
#' @param y Another object to use as reference
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("alignDt", function(x, y, ...) standardGeneric("alignDt"))


setGeneric("pData", getGeneric("pData", package = "Biobase"))
setGeneric("pData<-", getGeneric("pData<-", package = "Biobase"))
setGeneric("sampleNames", getGeneric("sampleNames", package = "Biobase"))
setGeneric("sampleNames<-", getGeneric("sampleNames<-", package = "Biobase"))

setGeneric("updateObject", getGeneric("updateObject", package = "BiocGenerics"))

setGeneric("alignRt", getGeneric("alignRt", package = "ProtGenerics"))
setGeneric("filterRt", getGeneric("filterRt", package = "ProtGenerics"))
setGeneric("intensity", getGeneric("intensity", package = "ProtGenerics"))
setGeneric("intensity<-", getGeneric("intensity<-", package = "ProtGenerics"))
setGeneric("peaks", getGeneric("peaks", package = "ProtGenerics"))
setGeneric("peaks<-", getGeneric("peaks<-", package = "ProtGenerics"))
setGeneric("rtime", getGeneric("rtime", package = "ProtGenerics"))
setGeneric("smooth", getGeneric("smooth", package = "ProtGenerics"))

# FIXME: Ask the xcms maintainer to move the findPeaks generic to ProtGenerics
setGeneric("findPeaks", getGeneric("findPeaks", package = "xcms"))
