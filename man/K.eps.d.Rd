% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Basic.R
\name{K.eps.d}
\alias{K.eps.d}
\title{Covariance Matrix for noise}
\usage{
K.eps.d(d, r, k, rho.eps)
}
\arguments{
\item{d}{an m*m matrix representing the distance.}

\item{r}{a numeric value between 0 and 1 representing the proportion of the independent error.}

\item{k}{a numerical value acting as the power term inside exponential.}

\item{rho.eps}{a numerical value representing the spatial range of noise.}
}
\value{
a matrix representing the covariance for noise
}
\description{
Function to calculate Covariance matrix for epsilon (Multivariate normal distribution)
}
\examples{
A <- matrix(runif(4),2)
K.eps.d(A,1,1,1)
}
