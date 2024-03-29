% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cosmic_signatures_v3.2}
\alias{cosmic_signatures_v3.2}
\title{Single base substitution signatures from COSMIC v3.2}
\format{An object of class \code{matrix} with 78 rows and 96 columns.}
\source{
\url{https://cancer.sanger.ac.uk/cosmic/signatures/SBS/}

\url{https://cancer.sanger.ac.uk/signatures/downloads/}
}
\usage{
cosmic_signatures_v3.2
}
\description{
A data set of 78 SBS signatures from COSMIC v3.2 (March 2021).
With a few exceptions, these signatures were extracted
using SigProfiler (as described in Alexandrov et al., 2020) from the 2,780
whole-genome variant calls produced by the ICGC/TCGA Pan Cancer Analysis of
Whole Genomes (PCAWG) Network.
The following signatures are annotated as "Possible sequencing artefacts":
SBS27, SBS43, SBS45-60.
}
\keyword{datasets}
