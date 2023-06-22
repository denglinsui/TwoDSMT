# Some basic functions
# Criterion for discovery: FDP, POWER
# Covariance matrix
#
#' FDP evaluation
#'
#' Function to evaluate FDP for upper tail hypothesis.
#' @param selected a vector consists of the indexes of selected signals.
#' @param mu a vector of the true population mean.
#'
#' @return FDP: a numerical value, the number of false discoveries/the number of discoveries.
#' @export
#'
#' @examples
#' mu <- rnorm(100);
#' selected <- sample(100,10);
#' fdp(selected,mu)
fdp = function(selected,mu){
  sum(mu[selected] <= 0) / max(1, length(selected))
}

#' POWER evaluation
#'
#' Function to evaluate POWER for upper tail hypothesis.
#' @param selected a vector consists of the indexes of selected signals.
#' @param mu a vector of the true population mean.
#'
#' @return POWER: a numeric value, the number of true discoveries/the number of true signals.
#' @export
#'
#' @examples
#'  mu <- rnorm(100);
#'  selected <- sample(100,10);
#'  Pow(selected,mu)
Pow = function(selected,mu){
  if(sum(mu>0)==0){
    return(0)
  }else{
    return(sum(mu[selected] > 0) / sum(mu > 0))
  }
}


#' Covariance Matrix for mu
#'
#' Function to calculate Covariance matrix for mu (Multivariate normal distribution).
#' @param d an m*m matrix representing the distance, m is the number of locations.
#' @param k a numerical value acting as the power term inside exponential.
#' @param rho.mu a numerical value representing the spatial range of mu.
#' @param sig2.mu a numerical value representing error variance.
#'
#' @return an m*m matrix representing the covariance for mu.
#' @export
#'
#' @examples
#' A <- matrix(runif(4),2)
#' K.mu.d(A,1,1,1)
K.mu.d <- function(d, k, rho.mu,sig2.mu){
  sig2.mu*exp(-(d/rho.mu)^k)
}


#' Covariance Matrix for noise
#'
#' Function to calculate Covariance matrix for epsilon (Multivariate normal distribution)
#' @param d an m*m matrix representing the distance.
#' @param r a numeric value between 0 and 1 representing the proportion of the independent error.
#' @param k a numerical value acting as the power term inside exponential.
#' @param rho.eps a numerical value representing the spatial range of noise.
#'
#' @return a matrix representing the covariance for noise
#' @export
#'
#' @examples
#' A <- matrix(runif(4),2)
#' K.eps.d(A,1,1,1)
K.eps.d <- function(d, r, k, rho.eps){
  (1-r)*(d==0) + r*exp(-(d/rho.eps)^k)
}
