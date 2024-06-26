% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align-GCIMSSample.R
\name{alignRt_ip}
\alias{alignRt_ip}
\title{Align a GCIMSSample in retention time with a multiplicative correction}
\usage{
alignRt_ip(object, min_start, rt_ref)
}
\arguments{
\item{object}{A \link{GCIMSSample} object}

\item{min_start}{minimun injection point, to calculate where to begin the spectrums and cut as few points as posible}

\item{rt_ref}{retention time reference}
}
\value{
The modified \link{GCIMSSample}
}
\description{
Align a GCIMSSample in retention time with a multiplicative correction
}
