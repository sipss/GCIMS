% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-GCIMSSample.R
\name{plot,GCIMSSample,ANY-method}
\alias{plot,GCIMSSample,ANY-method}
\title{Topographical plot of a GCIMSSample object}
\usage{
\S4method{plot}{GCIMSSample,ANY}(
  x,
  dt_range = NULL,
  rt_range = NULL,
  ...,
  remove_baseline = FALSE,
  trans = "cubic_root"
)
}
\arguments{
\item{x}{A \link{GCIMSSample} object}

\item{dt_range}{A numeric vector of length 2 with the drift time range to plot (in milliseconds)}

\item{rt_range}{A numeric vector of length 2 with the retention time range to plot (in seconds)}

\item{...}{Ignored}

\item{remove_baseline}{Set to \code{TRUE} to subtract the estimated baseline first}

\item{trans}{The transformation to the intensity values. "cubic_root" is the default. "intensity" is also valid.
See the \code{trans} argument in \code{\link[ggplot2:continuous_scale]{ggplot2::continuous_scale()}} for other possibilities.}
}
\value{
A plot of the GCIMSSample
}
\description{
Topographical plot of a GCIMSSample object
}
\examples{
dummy_obj <-GCIMSSample(
  drift_time = 1:2,
  retention_time = 1:3,
  data = matrix(1:6, nrow = 2, ncol = 3),
  gc_column = "Optional column name",
  drift_gas = "nitrogen",
  drift_tube_length = 98.0 # in mm
)
plot(dummy_obj)
}
