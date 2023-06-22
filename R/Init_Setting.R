#====================================
## Functions for Generating Signals
#====================================

#=================================
## 1 dimension
##
#' One dimension setting
#'
#' Basic function for initializing signal mu (Fix magnitude): 1 Dimension
#' @param point an m*1 matrix or numeric vector indicating the point location on one dimensional space, m is the number of locations.
#' @param method a character indicating the generating mechanism of mu, taking one of uc.spline, uc.unif and mvnorm.
#' @param magnitude a numeric value indicating the magnitude of signals.
#' @param mu.mean a numeric value shifting the signal (mvnorm mechanism).
#' @param Sigma.mu a matrix representing the covariance matrix for mu (mvnorm mechanism).
#' @param num_bump a integer indicating the number of bumps.(uc.unif & uc.spline mechanism).
#' @param single_bump_prop a numeric value indicating how many space one bump occupies.
#'
#' @return a list; mu is a numeric vector representing the true population mean; pis is a numeric vector representing the probability of being null.
#' @importFrom dplyr filter mutate
#' @importFrom fda create.bspline.basis eval.basis
#' @importFrom stats pnorm rbinom
#' @importFrom MASS mvrnorm
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' point <- sort(runif(100))
#' res1 <- mu.fix.gen.1D(point = point, method = "uc.spline",
#'                       magnitude = 1)
#' res2 <- mu.fix.gen.1D(point = point, method = "uc.unif",
#'                       magnitude = 1)
#' Dist.p <- as.matrix(dist(point))
#' Sigma.mu.p <- K.mu.d(Dist.p, k =1, rho.mu = 0.3, sig2.mu=1)
#' res3 <- mu.fix.gen.1D(point = point, method = "mvnorm",
#'                       magnitude = 1, mu.mean=0,
#'                       Sigma.mu = Sigma.mu.p)
mu.fix.gen.1D <- function(point, method, magnitude,mu.mean,Sigma.mu,
                          num_bump=2,single_bump_prop=1/6 #For uc.spline
){
  if(is.unsorted(point)){
    stop("The input point should be ordered")
  }
  if(is.vector(point)){
    point <- matrix(point,ncol = 1)
  }
  if(!(method %in% c("uc.unif","uc.spline","mvnorm"))){
    stop("Please make sure your input method is one of uc.unif, uc.spline and mvnorm")
  }

  m <- nrow(point)

  if(method == "uc.unif"){
    if(num_bump*single_bump_prop>1) warning("The condition for Splines is not satisfied.")
    #=== Choose the boundary
    center.ind <- floor(seq(0,1,length.out=num_bump+2)[-c(1,num_bump+2)]*m)
    each_er <- floor(m*single_bump_prop/2)

    #=== Init point
    point.basis <- point %>%  as.data.frame() %>% mutate(mu =0)
    colnames(point.basis) <- c("x","mu")

    #=== Start Generating
    for(cen.ind in center.ind){
      ind.low <- cen.ind-each_er; ind.up <- cen.ind+each_er
      p.low <- point[ind.low,1];p.up <- point[ind.up,1]

      #=== Init Basis
      n.break <- 3
      n.order <- 4
      n.end <- n.order + n.break -2

      sp.basis <- create.bspline.basis(c(p.low, p.up), dropind=c(1,n.end),
                                       breaks = seq(p.low, p.up, length.out = n.break),
                                       norder = n.order)

      point.basis.tmp <- point.basis[point.basis$x>=p.low & point.basis$x<=p.up,]
      point.basis.tmp.mu <- point.basis.tmp$mu+ 1 * rowSums(eval.basis(point.basis.tmp$x,basisobj = sp.basis))
      point.basis$mu[point.basis$x>=p.low & point.basis$x<=p.up] <-point.basis.tmp.mu

    }

    pis <- 1 * point.basis$mu
    pis[pis<=0.01] <- 0.01 #normalize
    mu <- magnitude * rbinom(n=length(pis),size=1,prob=pis)
  }


  if(method == "uc.spline"){
    if(num_bump*single_bump_prop>1) warning("The condition for Splines is not satisfied.")
    #=== Choose the boundary
    center.ind <- floor(seq(0,1,length.out=num_bump+2)[-c(1,num_bump+2)]*m)
    each_er <- floor(m*single_bump_prop/2)

    #=== Init point
    point.basis <- point %>%  as.data.frame() %>% mutate(mu =0)
    colnames(point.basis) <- c("x","mu")

    #=== Start Generating
    for(cen.ind in center.ind){
      ind.low <- cen.ind-each_er; ind.up <- cen.ind+each_er
      p.low <- point[ind.low,1];p.up <- point[ind.up,1]

      #=== Init Basis
      n.break <- 3
      n.order <- 4
      n.end <- n.order + n.break -2

      sp.basis <- create.bspline.basis(c(p.low, p.up), dropind=c(1,n.end),
                                       breaks = seq(p.low, p.up, length.out = n.break),
                                       norder = n.order)
      point.basis.tmp <- point.basis[point.basis$x>=p.low & point.basis$x<=p.up,]
      point.basis.tmp.mu <- point.basis.tmp$mu+ 1 * rowSums(eval.basis(point.basis.tmp$x,basisobj = sp.basis))
      point.basis$mu[point.basis$x>=p.low & point.basis$x<=p.up] <- point.basis.tmp.mu

    }

    mu <- magnitude * point.basis$mu
    pis <- as.numeric(mu<=0)

  }

  if(method == "mvnorm"){
    mu <- magnitude * mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = Sigma.mu)
    pis <- as.numeric(mu<=0)
  }

#  if(method == "mixture"){
#    m <- nrow(point)
#    pis <- numeric(m)
#    location <- as.numeric(abs(point[,1]-5)<1) + as.numeric(abs(point[,1]-18)<2)
#    pis[location] <- 0.9
#    pis[-location] <- 0.1
#
#    location <- rbinom(m,1,pis)
#    mu <- location*magnitude
#    pis <- 1-pis
#  }
  return(list(mu=mu,
              pis=pis))
}

#=================================


#=================================
## 2 dimension
##
#' Two dimension setting
#'
#' Init mu (Fix magnitude): 2 Dimension
#' @param point an m*2 matrix or numeric vector indicating the point location on one, m is the number of locations.
#' @param method a character indicating the generating mechanism of mu, taking one of uc.square, uc.circle and uc.mix.
#' @param magnitude a numeric value indicating the magnitude of signals.
#' @param mu.mean a numeric value shifting the signal (mvnorm mechanism).
#' @param Sigma.mu a matrix representing the covariance matrix for mu (mvnorm mechanism).
#'
#' @return a list; mu is a matrix representing the true population mean; pis is a matrix representing the probability of being null.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' m <- 900 # point size
#' m_x <- 30
#' m_y <- 30
#' point_x <- seq(0,1,length.out=m_x)
#' point_y <- seq(0,1,length.out=m_x)
#' point <- cbind(rep(point_x,times=length(point_y)),
#'                rep(point_y,each=length(point_x)))
#' res1 <- mu.fix.gen.2D(point,method="uc.circle",magnitude=1)
#' res2 <- mu.fix.gen.2D(point,method="uc.square",magnitude=1)
#' res3 <- mu.fix.gen.2D(point,method="uc.mix",magnitude=1)
mu.fix.gen.2D <- function(point, method, magnitude,mu.mean,Sigma.mu){
  m <- nrow(point)

  if(method == "uc.circle"|method == "uc.mix"){
    location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
    pis <- numeric(m)
    pis[location] <- 0.9
    pis[!location] <- 0.01
    mu <- magnitude * rbinom(n=length(pis),size=1,prob=pis)
    pis <- as.numeric(mu<=0)

  }

  if(method == "uc.square"|method == "uc.mix"){
    #=== Init Basis
    dis.low <-3/4 +sqrt(2)/16-1/8
    dis.up <- 3/4+sqrt(2)/16+1/8
    n.break <- 3
    n.order <- 4
    n.end <- n.order + n.break -2
    sp.basis <- create.bspline.basis(c(dis.low,dis.up), dropind=c(1,n.end),
                                     breaks = seq(dis.low, dis.up, length.out = n.break),
                                     norder = n.order)

    point.basis <- point %>%  as.data.frame() %>% mutate(mu.x = 0) %>% mutate(mu.y = 0)
    colnames(point.basis) <- c("x","y","mu.x","mu.y")

    #=== Basis product
    point.basis.tmp <- point.basis[point.basis$x>=dis.low & point.basis$x<=dis.up,]
    point.basis.tmp.mu.x <- sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(point.basis.tmp$x,basisobj = sp.basis))
    point.basis$mu.x[point.basis$x>=dis.low & point.basis$x<=dis.up] <-  point.basis.tmp.mu.x

    point.basis.tmp <- point.basis[point.basis$y>=dis.low & point.basis$y<=dis.up,]
    point.basis.tmp.mu.y <- sqrt(n.end)/sqrt(n.end-2) * rowSums(eval.basis(point.basis.tmp$y,basisobj = sp.basis))
    point.basis$mu.y[point.basis$y>=dis.low & point.basis$y<dis.up] <-point.basis.tmp.mu.y
    #sqrt(n.end)/sqrt(n.end-2) is for normalization
    point.basis <- point.basis %>% mutate(mu = point.basis$mu.x*point.basis$mu.y) %>% mutate(mu = mu)
    if(method == "uc.mix"){
      mu <- mu + magnitude * point.basis$mu
    }else{
      mu <- magnitude * point.basis$mu
    }
    pis <- as.numeric(mu<=0)
  }

  if(method == "mvnorm"){
    mu <- magnitude*mvrnorm(n = 1, mu = rep(mu.mean,m), Sigma = Sigma.mu)
    pis <- as.numeric(mu<=0)
  }

  if(method == "mixture"){
    m <- nrow(point)
    pis <- numeric(m)
    location <- (point[,1]-1/2)^2+(point[,2]-1/2)^2<(1/4)^2
    pis[location] <- 0.9
    pis[-location] <- 0.1

    location <- rbinom(m,1,pis)
    mu <- location*magnitude
    pis <- 1-pis
  }
  return(list(mu=mu,
              pis=pis))
}
#=================================



#' One dimension setting
#'
#' Generate mu & covariance matrix (One dimensional case)
#' @param mu_type a character indicating the type of mu, taking one of the value Sparse, Medium and Dense.
#' @param Cov_type a character indicating the type of covariance, taking one of the value Weak, Medium and Strong.
#' @param mu_gen_machine a character indicating generating mechanism, taking one of the value uc.spline, uc.unif and mvnorm.
#' @param magnitude a numeric value indicating the magnitude of signals.
#' @param point an m*1 matrix or numeric vector indicating the point location on one dimensional space.
#' @param single_bump_prop a numeric value indicating how many space one bump occupies.
#' @param Dist.p a m*m matrix representing the distance between points.
#'
#' @return list; mu is a numeric vector representing the true population mean; Sigma.eps.p is a matrix representing the covariance among noises.
#' @importFrom stats dist
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' m <- 900 # point size
#' point <- matrix(seq(0,1,length.out=m), ncol=1)
#' res <- Init_Setting_1D(mu_type = "Medium",Cov_type = "Medium",
#'                        mu_gen_machine="uc.spline",point=point)
Init_Setting_1D <- function(mu_type = "Medium",
                            Cov_type = "Medium",
                            mu_gen_machine = "uc.spline",
                            magnitude = 1,
                            point = point,
                            single_bump_prop = 1/20,
                            Dist.p = NULL){
  m <- nrow(point)
  #=== Generate Dist.p
  if(is.null(Dist.p)){
    Dist.p <- as.matrix(dist(point))
  }
  #===========================================================================
  #=== Generate mu
  #single_bump_prop <- 1/10
  if(mu_gen_machine=="uc.spline" | mu_gen_machine == "uc.unif"){
    ind.mu <- switch(mu_type,
                     Sparse = 1,
                     Medium = 2,
                     Dense = 4)
    num_bump = 1:ind.mu
    mu.res <- mu.fix.gen.1D(point, mu_gen_machine, magnitude,
                            num_bump=num_bump[ind.mu],single_bump_prop=single_bump_prop)
    mu <- mu.res$mu
  }else if(mu_gen_machine=="mvnorm"){
    mu.mean <- switch(mu_type,
                      Sparse = -2.5,
                      Medium = -2,
                      Dense = -1)
    Sigma.mu <- K.mu.d(Dist.p, k =1, rho.mu = 0.3, sig2.mu=3)
    #set.seed(10)
    mu.res <- mu.fix.gen.1D(point, mu_gen_machine,
                            magnitude = magnitude, mu.mean = mu.mean,
                            Sigma.mu = Sigma.mu)
    mu <- mu.res$mu
  }
  #===========================================================================

  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0.1,0.3,0.5)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)

  # For testing
  r <- c(0.5,0.9,0.9)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)

  r <- c(0.5,0.8,0.6)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)

  ind.Cov <- switch(Cov_type,
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
#  if(is.null(Dist.p)){
#    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
#  }else{
    Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
#  }
  #===========================================================================

  return(list(mu = mu,
              Sigma.eps.p = Sigma.eps.p))
}

#' Two dimension setting
#'
#' Generate mu & covariance matrix (Two dimensional case)
#' @param mu_type a character indicating the type of mu, taking one of the value Sparse, Medium and Dense.
#' @param Cov_type a character indicating the type of covariance, taking one of the value Weak, Medium and Strong.
#' @param mu_gen_machine a character indicating generating mechanism, taking one of the value uc.square, uc.circle and uc.mix.
#' @param magnitude a numeric value indicating the magnitude of signals.
#' @param point an m*2 matrix or numeric vector indicating the point location on two dimensional space.
#' @param Dist.p The distance matrix between points.
#'
#' @return list; mu is a matrix representing the true population mean; Sigma.eps.p is a matrix representing the covariance among noises.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' m <- 900 # point size
#' m_x <- 30
#' m_y <- 30
#' point_x <- seq(0,1,length.out=m_x)
#' point_y <- seq(0,1,length.out=m_x)
#' point <- cbind(rep(point_x,times=length(point_y)),
#'                rep(point_y,each=length(point_x)))
#' res <- Init_Setting_2D(mu_type = "Medium",Cov_type = "Medium",
#'                        mu_gen_machine="uc.square",point=point)
Init_Setting_2D <- function(mu_type = "Medium",
                            Cov_type = "Medium",
                            mu_gen_machine = NULL,
                            magnitude = 1,
                            point = point,
                            #single_bump_prop = 1/20,
                            Dist.p = NULL){
  m <- nrow(point)

  #===========================================================================
  #=== Generate mu
  #single_bump_prop <- 1/10
  methods <- c("uc.square","uc.circle","uc.mix")
  ind.mu <- switch(mu_type,
                   Sparse = 1,
                   Medium = 2,
                   Dense = 3)
  num_bump = 1:ind.mu
  mu.res <- mu.fix.gen.2D(point, method=methods[ind.mu], magnitude)
  mu <- mu.res$mu
  #===========================================================================

  #===========================================================================
  #=== Generate Covariance Matrix
  r <- c(0.1,0.3,0.5)
  k <- c(1,1,2)
  rho.eps <- c(0.1,0.1,0.3)

  # For testing
  r <- c(0.5,0.9,0.9)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)

  r <- c(0.5,0.8,0.6)
  k <- c(1,1,2)
  rho.eps <- c(0.05,0.1,0.2)

  ind.Cov <- switch(Cov_type,
                    Weak = 1,
                    Medium = 2,
                    Strong = 3)
#  if(is.null(Dist.p)){
#    Sigma.eps.p <- Sigma.eps(m, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov],point)
#  }else{
  Sigma.eps.p <- K.eps.d(Dist.p, r[ind.Cov], k[ind.Cov], rho.eps[ind.Cov])
#  }
  #===========================================================================

  return(list(mu = mu,
              Sigma.eps.p = Sigma.eps.p))
}

