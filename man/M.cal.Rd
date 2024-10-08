% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M.cal.R
\name{M.cal}
\alias{M.cal}
\title{Calculation of projection matrix.}
\usage{
M.cal(Y_raw, Design_matrix, z.states.k, factor.to.test)
}
\arguments{
\item{Y_raw}{A list of matrix of pseudo-bulk data for each cell type.}

\item{Design_matrix}{the design matrix of all the covariates of the data
model for gene expression.}

\item{z.states.k}{DE states in each cell type(0: normal; 1: disease).}

\item{factor.to.test}{A phenotype name, e.g. "disease", or a vector of
contrast terms, e.g. c("disease", "case", "control").}
}
\value{
the projection matrix
}
\description{
It calculates the projection matrix of data linear model based on
design matrix.
}
