% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-io-mea.R
\name{read_mea}
\alias{read_mea}
\title{Read .mea files (from GAS Dortmund)}
\usage{
read_mea(filename)
}
\arguments{
\item{filename}{A .mea or a .mea.gz path to a file}
}
\value{
The GC-IMS sample in a \linkS4class{GCIMSSample} object
}
\description{
This function reads a .mea file (supporting gzip compressed .mea.gz files as well)
and returns a GCIMS object
}
\details{
Thanks to Daniel Sanders and Thomas Wortelmann from  \href{https://www.gas-dortmund.de/}{GAS Dortmund} for providing the
necessary support to implement this function.
}
\examples{
mea_file <- system.file("extdata/sample_formats/small.mea.gz", package = "GCIMS")
sample <- read_mea(mea_file)

}
