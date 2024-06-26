% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spatial_Detection.R
\name{Spatial_Detect_exact_grp_BH_down}
\alias{Spatial_Detect_exact_grp_BH_down}
\title{Two dimension detection (BH based)}
\usage{
Spatial_Detect_exact_grp_BH_down(
  Tm,
  Ta,
  Va,
  VmVa.cov,
  ind,
  q = 0.1,
  max.rej = NULL,
  pis = NULL,
  cutoff = NULL,
  Tm.star = Inf,
  Ta.star = Inf,
  const = 0,
  seed = 0,
  dig = 7,
  tau.tm = 1,
  tau.ta = 1,
  mua = NULL,
  is.quick.stop = F,
  stop.step = 10,
  is.fullset = F,
  is.fixRectShape = F,
  sig_ratio = NULL,
  EmpMethod = c("NPEB", "DeCor")
)
}
\arguments{
\item{Tm}{a numerical vector representing the primary statistics.}

\item{Ta}{a numerical vector representing the auxiliary statistics.}

\item{Va}{a numerical vector representing the variance of auxiliary statistics.}

\item{VmVa.cov}{a numerical vector representing the covariance between primary and auxiliary statistics.}

\item{ind}{an integer vector representing the indexes for implementing nonparametric empirical Bayes.}

\item{q}{a numerical value representing the target fdr level.}

\item{max.rej}{an integer indicating the current maximum rejection number that controls fdr under q.}

\item{pis}{a numerical vector representing the probability of being null.}

\item{cutoff}{a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;}

\item{Tm.star}{a numerical value representing the threshold found by one dimensional detection based on the primary statistics.}

\item{Ta.star}{a numerical value representing the threshold found by one dimensional detection based on the auxiliary statistics.}

\item{const}{a numerical value ensuring the selection stability, usually be target fdr level.}

\item{seed}{a numerical value indicating the simulation seed (can be omit).}

\item{dig}{an integer representing the digit precision.}

\item{tau.tm}{a numerical value representing the censor level for p-values.}

\item{tau.ta}{a numerical value representing the censor level for p-values.}

\item{mua}{a numerical vector representing the population mean of auxiliary statistics (Usually unavailable).}

\item{is.quick.stop}{a bool value representing whether we should stop early.}

\item{stop.step}{an integer indicating the maximal step for no valuable finding(estimate fdr<target fdr level).}

\item{is.fullset}{a bool value representing whether we take all partitions of full set for implementing NPEB.}

\item{is.fixRectShape}{a bool value representing whether the shape of rectangle is fixed}

\item{sig_ratio}{a numerical vector representing the estimated ratio of magnitude.}

\item{EmpMethod}{a character indicating which method we choose to estimate the empirical distribution.}
}
\value{
list;
selected is an integer vector consisting of selected locations;
tm0 a numerical value representing the threshold for the primary statistics;
ta0 a numerical value representing the threshold for the auxiliary statistics;
final.fdr a numerical value representing the final estimated false discovery rate.
}
\description{
Searching for the optimal threshold for the two dimensional BH procedure based on the primary and auxiliary statistics.
}
\examples{
library(dplyr)
library(fda)
library(qvalue)
# Basic Setting
set.seed(1234)
dig <- 7 # digits precision
n <- 1 # observation number at each point
m <- 120 # point size
q <- 0.1 # fdr level
const <- q # const for stablizing
# Generate data
point <- matrix(seq(0,1,length.out=m), ncol=1)
Dist.p <- as.matrix(dist(point))
I_S <- Init_Setting_1D(mu_type = "Dense",Cov_type = "Weak",
                       magnitude=2,mu_gen_machine="mvnorm",point=point)
mu <- I_S$mu
Sigma.eps.p <- I_S$Sigma.eps.p
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma.eps.p)
X <- matrix(X, nrow = n)
sgm <- sqrt(diag(Sigma.eps.p))
Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm # Generate primary statistics
Tm <- round(Tm,dig)
p.value <- 1 - pnorm(Tm) # Generate one-side p.value
# Run BH algorithm
result.BH <- qvalue(p.value,pi0 = 1)
selected.BH <- which(result.BH$qvalues<=q)

# Run 1D algorithm (1D(BH))
res.1D <- OneD_Detect(Tm, q,  pis = NULL,const = const, tau.tm = 1)
selected.1D <- res.1D$selected
Tm.star <- res.1D$tm.min
max.rej <- res.1D$max.rej

# Run 2D algorith with 2 neighbors (2D(BH))
hh <- 2 # neighbor number
hh.seq <- rep(hh,m)
Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X,
                                 Dist = Dist.p, Sigma.eps = Sigma.eps.p,
                                 detect.m = "top.k")
Ta <- Neigh_Detect_res$Ta
Va <- Neigh_Detect_res$Va
VmVa.cov <- Neigh_Detect_res$VmVa.cov
ind <- Neigh_Detect_res$ind
mua <- Neigh_Detect_res$mua

res.2D <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, VmVa.cov, ind,
                                           q, max.rej,
                                           Tm.star = Tm.star,
                                           Ta.star = Inf,
                                           const = const)
selected.2D <- res.2D$selected
tm <- res.2D$tm0
ta <- res.2D$ta0

# fdp and power
fdp(selected.BH,mu);fdp(selected.1D,mu);fdp(selected.2D,mu)
Pow(selected.BH,mu);Pow(selected.1D,mu);Pow(selected.2D,mu)
}
