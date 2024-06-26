% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spatial_Detection.R
\name{Neigh_Partition}
\alias{Neigh_Partition}
\title{Neighbor Partition}
\usage{
Neigh_Partition(S.neigh)
}
\arguments{
\item{S.neigh}{a list consists of m vector, i-th vector is the indexes of the i-th location's neighbors.}
}
\value{
an integer vector representing the indexes for implementing nonparametric empirical Bayes;
}
\description{
To make the neighbors of different location do not intersect.
}
\examples{
m <- 900
point <- matrix(seq(0,1,length.out=m), ncol=1)
Dist.p <- as.matrix(dist(point))
hh.seq <- rep(2,m)
S.neigh <- lapply(1:m,
                 function(i){
                 ind1 <- order(Dist.p[i,])[2:(hh.seq[i]+1)]
                 ind1})
ind <- Neigh_Partition(S.neigh)
}
