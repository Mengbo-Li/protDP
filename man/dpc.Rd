% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpc.R
\name{dpc}
\alias{dpc}
\title{Detection probability curve for label free shotgun proteomics data
assuming observed normal intensities}
\usage{
dpc(nuis, maxit = 100, eps = 1e-04, b1.upper = 1)
}
\arguments{
\item{nuis}{The list of nuisance parameters.}

\item{maxit}{Maximum number of iterations.}

\item{eps}{Convergence tolerance.}

\item{b1.upper}{Upper bound for beta_1. Typically 1.}
}
\value{
Fitted parameters and the fitting history.
}
\description{
Detection probability curve for label free shotgun proteomics data
assuming observed normal intensities
}
\examples{
## See the vignettes.
}
