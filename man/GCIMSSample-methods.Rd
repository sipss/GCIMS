% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-class-GCIMSSample.R,
%   R/description-GCIMSSample.R, R/peaks-GCIMSSample.R
\name{GCIMSSample-methods}
\alias{GCIMSSample-methods}
\alias{[.GCIMSSample}
\alias{dim.GCIMSSample}
\alias{subset.GCIMSSample}
\alias{subset}
\alias{description,GCIMSSample-method}
\alias{description<-,GCIMSSample,ANY-method}
\alias{peaks,GCIMSSample-method}
\alias{peaks<-,GCIMSSample-method}
\title{Methods for the GCIMSSample class}
\usage{
\method{[}{GCIMSSample}(x, i, j, ...)

\method{dim}{GCIMSSample}(x)

\method{subset}{GCIMSSample}(x, dt_idx = NULL, rt_idx = NULL, dt_range = NULL, rt_range = NULL, ...)

\S4method{description}{GCIMSSample}(object)

\S4method{description}{GCIMSSample,ANY}(object) <- value

\S4method{peaks}{GCIMSSample}(object)

\S4method{peaks}{GCIMSSample}(object) <- value
}
\arguments{
\item{x}{A GCIMSSample object}

\item{i}{index for drift time to subset}

\item{j}{index for retention time to subset}

\item{...}{ignored}

\item{dt_idx}{A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)}

\item{rt_idx}{A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)}

\item{dt_range}{The minimum and maximum drift times to extract (length 2 vector)}

\item{rt_range}{The minimum and maximum retention times to extract (length 2 vector)}

\item{object}{A \link{GCIMSSample} object}

\item{value}{A data frame with the peak list}
}
\value{
\code{[}: object \code{x} with features \code{i} and cells \code{j}

An integer vector with the number of rows and columns of the matrix

\code{subset}: A subsetted \code{GCIMSSample} object
}
\description{
Methods for the GCIMSSample class
}
\section{Functions}{
\itemize{
\item \code{[}: Simple subsetter for \linkS4class{GCIMSSample} objects

\item \code{dim(GCIMSSample)}: Dimension of the data matrix

\item \code{subset(GCIMSSample)}: Subset a \linkS4class{GCIMSSample} object

\item \code{description(GCIMSSample)}: Get the description

\item \code{description(object = GCIMSSample) <- value}: Set the description

\item \code{peaks(GCIMSSample)}: Get the peak list

\item \code{peaks(GCIMSSample) <- value}: Set the peak list

}}
\examples{
# `[' examples

obj <- GCIMSSample(drift_time=1:2, retention_time=1:3, data = matrix(1:6, nrow=2, ncol=3))
dim(obj)
}
\seealso{
\code{\link[base:subset]{base::subset()}}
}
