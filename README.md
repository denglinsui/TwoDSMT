
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TwodSpatial

<!-- badges: start -->
<!-- badges: end -->

The goal of TwoDSMT is to conduct multiple hypothesis testing on the
spatial domain using neighboring information.

## Installation

You can install the development version of TwoDSMT like so:

``` r
R CMD INSTALL TwoDSMT_0.0.0.9000.tar.gz
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TwoDSMT)
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
# Run SABHA algorithm
tau = 0.5;eps = 0.1; TV_bd = 2
alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
SABHA_method = function(pvals, qhat, alpha, tau){
pvals[pvals>tau] = Inf
khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
which(qhat*pvals<=alpha*khat/length(pvals))
}
sab.selected <- SABHA_method(p.value, qhat, q, tau)


# Run 2D algorith with 2 neighbors (2D(SA))
ws.sabha.fun <- function(x){1/x}
hh <- 2
hh.seq <- rep(hh,m)
res.2D.sabha <- TwoDwSMT(X, hh.seq=hh.seq,pis=qhat,ws.fun = ws.sabha.fun,Dist.p=Dist.p,
                        Sigma.eps.p=Sigma.eps.p, q=q, const = const,n.group.max=5)
selected.2D.sabha <- res.2D.sabha$selected
tm <- res.2D.sabha$tm0
ta <- res.2D.sabha$ta0
# fdp and power
fdp(sab.selected,mu);fdp(selected.2D.sabha,mu)
#> [1] 0
#> [1] 0
Pow(sab.selected,mu);Pow(selected.2D.sabha,mu)
#> [1] 0.6363636
#> [1] 0.8181818
```
