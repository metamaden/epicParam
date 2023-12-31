% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epicParam-class.R
\name{epicParam}
\alias{epicParam}
\title{Make new object of class epicParam}
\usage{
epicParam(
  bulkExpression,
  referenceExpression,
  cellScaleFactors = NULL,
  z.var = NULL,
  returnInfo = FALSE
)
}
\arguments{
\item{bulkExpression}{Bulk mixed signals matrix of samples, which can be 
matched to single-cell samples.}

\item{referenceExpression}{Signature matrix of cell type-specific signals. 
If not provided, can be computed from a
provided ExpressionSet containing single-cell data.}

\item{cellScaleFactors}{Cell size factor transformations of length equal to 
the K cell types to deconvolve.}

\item{z.var}{Variance table for the reference matrix.}

\item{returnInfo}{Whether to return metadata and original method outputs 
with predicted proportions.}
}
\description{
Main constructor for class \linkS4class{epicParam}.
}
\references{
Racle, Julien, and David Gfeller. 2020. “EPIC: A Tool to Estimate the 
Proportions of Different Cell Types from Bulk Gene Expression Data.” In 
Bioinformatics for Cancer Immunotherapy: Methods and Protocols, edited by 
Sebastian Boegel, 233–48. Methods in Molecular Biology. New York, NY: 
Springer US. https://doi.org/10.1007/978-1-0716-0327-7_17.
}
