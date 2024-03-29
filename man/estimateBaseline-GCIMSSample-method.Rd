% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline-GCIMSSample.R
\name{estimateBaseline,GCIMSSample-method}
\alias{estimateBaseline,GCIMSSample-method}
\alias{baseline,GCIMSSample-method}
\alias{baseline<-,GCIMSSample-method}
\title{Estimate the baseline of a GCIMS Sample using a connect local minima algorithm}
\usage{
\S4method{estimateBaseline}{GCIMSSample}(
  object,
  dt_peak_fwhm_ms,
  dt_region_multiplier,
  rt_length_s,
  remove = TRUE
)

\S4method{baseline}{GCIMSSample}(
  object,
  dt_range = NULL,
  rt_range = NULL,
  dt_idx = NULL,
  rt_idx = NULL,
  .error_if_missing = TRUE
)

\S4method{baseline}{GCIMSSample}(object) <- value
}
\arguments{
\item{object}{A \link{GCIMSSample} object}

\item{dt_peak_fwhm_ms}{Full Width at Half Maximum in milliseconds. Used to
determine the length of the regions where local minima are searched.}

\item{dt_region_multiplier}{A multiplier to calculate the region}

\item{rt_length_s}{The length of the baseline region. It should be
comparable or longer than the peak width}

\item{remove}{A boolean, if TRUE it removes the baseline from the intensity}

\item{dt_range}{The minimum and maximum drift times to extract (length 2 vector)}

\item{rt_range}{The minimum and maximum retention times to extract (length 2 vector)}

\item{dt_idx}{A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)}

\item{rt_idx}{A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)}

\item{.error_if_missing}{A logical. If \code{TRUE}, raise error if baseline has not been estimated. If \code{FALSE} returns \code{NULL} instead.}

\item{value}{A matrix with the sample baseline of the same dimensions as \code{dim(object)}}
}
\value{
The modified \link{GCIMSSample}
}
\description{
The baseline is estimated by connecting local minima and interpolating from those.
The local minima are identified as "the minima in each region of length x"
The length of the regions are estimated as \verb{fwhm * a multiplier / 2.3482}. This
assumes it's several times
}
\section{Functions}{
\itemize{
\item \code{baseline(GCIMSSample)}: Get the baseline

\item \code{baseline(GCIMSSample) <- value}: Set the baseline

}}
