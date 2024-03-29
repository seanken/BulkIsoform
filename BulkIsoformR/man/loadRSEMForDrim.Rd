% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadRSEM.R
\name{loadRSEMForDrim}
\alias{loadRSEMForDrim}
\title{Load RSEM in testing format}
\usage{
loadRSEMForDrim(files, nams = c())
}
\arguments{
\item{files}{A name array of isoform files from RSEM}

\item{nams}{An array of names, if not given just uses files as names. Used as sample names.}
}
\value{
Ready from DRIM or other testing
}
\description{
A function to load a list of data from RSEM as count matrix ready for DTU testing
}
