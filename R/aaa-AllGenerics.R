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
setClassUnion("matrixOrNULL", c("array", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))

# S3 classes can be used as S4 slots registering them as "old" classes
# (old because S4 is newer than S3, not because S3 is deprecated, not great naming)
setOldClass("numeric_version")


#' GCIMS Generics
#'
#' Generics defined at the GCIMS package. We are open to moving them
#' to an existing generics-only package if you need so.
#'
#' @name GCIMS-generics
NULL

#' @describeIn GCIMS-generics Get drift time vector
#'
#' @return A numeric vector with the drift time
#'
#' @param object An object with drift time
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
#' @examples
#' x <- GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3)
#' )
#' dtime(x) # c(1,2)
#'
setGeneric("dtime", function(object, ...) standardGeneric("dtime"))

#' @describeIn GCIMS-generics Get chromatogram
#'
#' @return A chromatogram
#'
#' @param object An object with chromatograms
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getChromatogram", function(object, ...) standardGeneric("getChromatogram"))

#' @describeIn GCIMS-generics Get spectrums
#'
#' @return An ion mobility spectrum
#' @param object An object that has spectra
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getSpectrum", function(object, ...) standardGeneric("getSpectrum"))

#' @describeIn GCIMS-generics Get a Sample
#'
#' @param object An object where we can extract samples from
#' @param ... Further arguments, possibly used by downstream methods.
#' @return A sample object
#' @export
setGeneric("getSample", function(object, ...) standardGeneric("getSample"))

#' @describeIn GCIMS-generics Get the Total Ion Spectrum
#'
#' @return The Total Ion Spectrum as a numeric vector or a matrix
#' (depending if the object is one sample or several)
#' @param object An object to extract its total ion spectrum
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getTIS", function(object, ...) standardGeneric("getTIS"))

#' @describeIn GCIMS-generics Get the Reverse Ion Chromatogram
#'
#' @return The Reverse Ion Chromatogram, as a numeric vector or a
#' matrix (depending if the object is one sample or several)
#' @param object An object to extract its reverse ion chromatogram
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("getRIC", function(object, ...) standardGeneric("getRIC"))


#' @describeIn GCIMS-generics Plot total ion spectrum
#'
#' @return A plot
#' @param object An object to plot the total ion spectrum
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotTIS", function(object, ...) standardGeneric("plotTIS"))

#' @describeIn GCIMS-generics Plot Reverse Ion Chromatogram
#'
#' @return A plot
#' @param object An object to plot its reverse ion chromatogram
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("plotRIC", function(object, ...) standardGeneric("plotRIC"))


#' @describeIn GCIMS-generics Filter in Drift time
#' @return The object, modified
#' @param object An object to filter
#' @param ... Further arguments, possibly used by downstream methods.
#' @export
setGeneric("filterDt", function(object, ...) standardGeneric("filterDt"))


#' Format as a list
#'
#' @return A list with a brief description/representation of the object
#'
#' @keywords internal
setGeneric("describeAsList", function(object, ...) standardGeneric("describeAsList"))

#' @describeIn GCIMS-generics Decimate an object
#'
#' @param object An object to decimate
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("decimate", function(object, ...) standardGeneric("decimate"))


#' @describeIn GCIMS-generics Align an object
#'
#' @param object An object to align
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("align", function(object, ...) standardGeneric("align"))


#' @describeIn GCIMS-generics Align an object in drift time
#'
#' @param x An object to align
#' @param y Another object to use as reference
#' @param ... Further arguments, possibly used by downstream methods.
#' @return The object, modified
#' @export
setGeneric("alignDt", function(x, y, ...) standardGeneric("alignDt"))

#' @importMethodsFrom Biobase pData
#' @export
setGeneric("pData", getGeneric("pData", package = "Biobase"))

#' @importMethodsFrom Biobase "pData<-"
#' @export
setGeneric("pData<-", getGeneric("pData<-", package = "Biobase"))


#' @importMethodsFrom Biobase sampleNames
#' @export
setGeneric("sampleNames", getGeneric("sampleNames", package = "Biobase"))

#' @importMethodsFrom Biobase "sampleNames<-"
#' @export
setGeneric("sampleNames<-", getGeneric("sampleNames<-", package = "Biobase"))


#' @importMethodsFrom BiocGenerics updateObject
#' @export
setGeneric("updateObject", getGeneric("updateObject", package = "BiocGenerics"))


#' @importMethodsFrom ProtGenerics alignRt
#' @export
setGeneric("alignRt", getGeneric("alignRt", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics filterRt
#' @export
setGeneric("filterRt", getGeneric("filterRt", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics intensity
#' @export
setGeneric("intensity", getGeneric("intensity", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics "intensity<-"
#' @export
setGeneric("intensity<-", getGeneric("intensity<-", package = "ProtGenerics"))

#' @importMethodsFrom Biobase description
#' @export
setGeneric("description", getGeneric("description", package = "Biobase"))

#' @importMethodsFrom Biobase "description<-"
#' @export
setGeneric("description<-", getGeneric("description<-", package = "Biobase"))


#' @importMethodsFrom ProtGenerics peaks
#' @export
setGeneric("peaks", getGeneric("peaks", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setGeneric("peaks<-", getGeneric("peaks<-", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics rtime
#' @export
setGeneric("rtime", getGeneric("rtime", package = "ProtGenerics"))

#' @importMethodsFrom ProtGenerics smooth
#' @export
setGeneric("smooth", getGeneric("smooth", package = "ProtGenerics"))


# FIXME: Ask the xcms maintainer to move the findPeaks generic to ProtGenerics
#        Redefining creates a conflict between packages, but reusing the generic
#        requires a dependency I would rather avoid
#' Find Peaks in an object
#' @param object An object to find peaks on
#' @param ... Additional arguments for downstream methods
#' @return The object, with found peaks
#'
#' @export
setGeneric("findPeaks", function(object, ...) standardGeneric("findPeaks"))
#setGeneric("findPeaks", getGeneric("findPeaks", package = "xcms"))


#' @describeIn GCIMS-generics Estimate the baseline in an object
#' @param object An object to estimate the baseline
#' @param ... Additional arguments for downstream methods
#' @return The object, with a baseline estimated
#'
#' @export
setGeneric("estimateBaseline", function(object, ...) standardGeneric("estimateBaseline"))


#' @describeIn GCIMS-generics Get the baseline of an object
#' @param object An object to get the baseline
#' @param ... Additional arguments for downstream methods
#' @return The baseline of the object
#'
#' @export
setGeneric("baseline", function(object, ...) standardGeneric("baseline"))

#' @describeIn GCIMS-generics Set the baseline of an object
#' @param object An object to set the baseline
#' @param value baseline to set
#' @return The object
#'
#' @export
setGeneric("baseline<-", function(object, value) standardGeneric("baseline<-"))


#' @describeIn GCIMS-generics Integrate peaks of an object
#' @param object An object to get the baseline
#' @param ... Additional arguments for downstream methods
#' @return The object, with integrated peaks
#'
#' @export
setGeneric("integratePeaks", function(object, ...) standardGeneric("integratePeaks"))
