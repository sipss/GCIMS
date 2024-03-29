% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter-GCIMSSample.R
\name{filterDt,GCIMSSample-method}
\alias{filterDt,GCIMSSample-method}
\title{Filter GCIMSSample samples by drift time}
\usage{
\S4method{filterDt}{GCIMSSample}(object, dt_range)
}
\arguments{
\item{object}{A \link{GCIMSSample} object}

\item{dt_range}{The minimum and maximum drift times to extract (length 2 vector)}
}
\value{
A subset of the sample, only in the selected \code{dt_range}
}
\description{
Filter GCIMSSample samples by drift time
}
\examples{
sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
s <- read_mea(sample_file)
s <- filterDt(s, dt_range = c(5, 9.5))
}
