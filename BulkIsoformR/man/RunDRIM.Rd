% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunDRIM.R
\name{RunDRIM}
\alias{RunDRIM}
\title{Run DRIM-Seq}
\usage{
RunDRIM(cts, meta, form, coef = 2, min_samps_feature_expr = 3)
}
\arguments{
\item{cts}{A dataframe of counts as produced by the load DRIM functions for RSEM and Salmon. First two columns should be feature_id and gene_id, the transcript and gene names.}

\item{meta}{A dataframe with meta data. Needs a column sample_id whose names match the names of samples in cts.}

\item{form}{The formula to use for testing, using columns in meta.}

\item{coef}{The coefficient to return}

\item{min_samps_feature_expr}{The number of samples that a given trascript must be expressed in.}
}
\value{
A dmTest object from DRIM-seq
}
\description{
Runs DRIM-seq to see to see differential transcript useage (DTU).
}
