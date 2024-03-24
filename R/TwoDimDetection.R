#' 2d-SMT: Two dimensional spatial multiple testing procedure

#' @description
#' Searching for the optimal threshold for the two dimensional spatial multiple testing procedure (2d-SMT) based on the primary and auxiliary statistics.
#' @param X a n*m matrix representing the observations, m is the number of locations, n is the number of observations at each location.
#' @param Sigma.eps.p a m*m matrix representing the covariance matrix for noises.
#' @param hh.seq an integer vector indicating the number of neighbors for each location.
#' @param Dist.p a m*m matrix representing the distance between points.s
#' @param detect.m a character representing the detecting method, currently only the k-nearest neighbor, "top.k", is available.
#' @param q a numerical value representing the target fdr level.
#' @param pis  a numerical vector representing the probability of being null.
#' @param cutoff a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param seed a numerical value indicating the simulation seed (can be omit).
#' @param dig an integer representing the digit precision.
#' @param tau.tm a numerical value representing the censor level for p-values.
#' @param tau.ta a numerical value representing the censor level for p-values.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' Tm is a numerical vector consisting of primary statistics;
#' Ta is a numerical vector consisting of auxiliary statistics;
#' tm0 a numerical value representing the threshold for the primary statistics;
#' ta0 a numerical value representing the threshold for the auxiliary statistics;
#' final.fdr a numerical value representing the final estimated false discovery rate.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' library(qvalue)
#' # Basic Setting
#' set.seed(1234)
#' dig <- 7 # digits precision
#' n <- 1 # observation number at each point
#' m <- 120 # point size
#' q <- 0.1 # fdr level
#' const <- q # const for stablizing
#' # Generate data
#' point <- matrix(seq(0,1,length.out=m), ncol=1)
#' Dist.p <- as.matrix(dist(point))
#' I_S <- Init_Setting_1D(mu_type = "Dense",Cov_type = "Weak",
#'                        magnitude=2,mu_gen_machine="mvnorm",point=point)
#' mu <- I_S$mu
#' Sigma.eps.p <- I_S$Sigma.eps.p
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma.eps.p)
#' X <- matrix(X, nrow = n)
#' sgm <- sqrt(diag(Sigma.eps.p))
#' Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm # Generate primary statistics
#' Tm <- round(Tm,dig)
#' p.value <- 1 - pnorm(Tm) # Generate one-side p.value
#'# Run BH algorithm
#' result.BH <- qvalue(p.value,pi0 = 1)
#' selected.BH <- which(result.BH$qvalues<=q)
#'
#'# Run 1D algorithm (1D(BH))
#' res.1D <- OneD_Detect(Tm, q,  pis = NULL,const = const, tau.tm = 1)
#' selected.1D <- res.1D$selected
#' Tm.star <- res.1D$tm.min
#' max.rej <- res.1D$max.rej
#'
#'# Run 2D algorith with 2 neighbors (2D(BH))
#' hh <- 2
#' hh.seq <- rep(hh,m)
#' res.2D <- TwoDSMT(X, hh.seq=hh.seq, Dist.p=Dist.p, Sigma.eps.p=Sigma.eps.p, q=q, const = const)
#' selected.2D <- res.2D$selected
#' tm <- res.2D$tm0
#' ta <- res.2D$ta0
#'
#' # fdp and power
#' fdp(selected.BH,mu);fdp(selected.2D,mu)
#' Pow(selected.BH,mu);Pow(selected.2D,mu)
TwoDSMT <- function(X,
                   hh.seq,
                   Dist.p,
                   Sigma.eps.p,
                   detect.m = "top.k",
                   q=0.1,
                   pis = NULL, cutoff = NULL,
                   const = 0,
                   seed = 0,
                   dig=7,
                   tau.tm = 1,
                   tau.ta = 1){
  n <- dim(X)[1]
  m <- dim(X)[2]
  sgm <- sqrt(diag(Sigma.eps.p))
  Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm # Generate primary statistics
  Tm <- round(Tm,dig)

  # Draw necessary results from one dimensional detection
  res.1D <- OneD_Detect(Tm, q,  pis = pis,const = const, tau.tm = tau.tm)
  selected.1D <- res.1D$selected
  Tm.star <- res.1D$tm.min
  max.rej <- res.1D$max.rej

  # Generate prelimaries for two dimensional detection
  Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X, Dist = Dist.p, Sigma.eps = Sigma.eps.p, detect.m = detect.m)
  Ta <- Neigh_Detect_res$Ta
  Va <- Neigh_Detect_res$Va
  VmVa.cov <- Neigh_Detect_res$VmVa.cov
  ind <- Neigh_Detect_res$ind
  mua <- Neigh_Detect_res$mua

  # Run 2dBH
  res.2D <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, VmVa.cov, ind,
                                             q, max.rej,
                                             Tm.star = Tm.star,
                                             Ta.star = Inf,
                                             const = const,
                                             seed = seed,
                                             mua = mua)
  return(list(selected = res.2D$selected,
              Tm = Tm,
              Ta = Ta,
              tm0 = res.2D$tm0,
              ta0 = res.2D$ta0,
              final.fdr=res.2D$final.fdr))
}


#' 2d-SMT(w): Two dimensional weighted spatial multiple testing procedure
#' @description
#' Searching for the optimal threshold for the two dimensional spatial multiple testing procedure (2d-SMT) procedure based on the primary and auxiliary statistics.
#' @param X a n*m matrix representing the observations, m is the number of locations, n is the number of observations at each location.
#' @param Sigma.eps.p a m*m matrix representing the covariance matrix for noises.
#' @param hh.seq an integer vector indicating the number of neighbors for each location.
#' @param ws a numerical vector indicating the non-null weights.
#' @param pis  a numerical vector representing the probability of being null.
#' @param ws.fun a function for weighting, depending on the null proportion.
#' @param Dist.p a m*m matrix representing the distance between points.s
#' @param detect.m a character representing the detecting method, currently only the k-nearest neighbor, "top.k", is available.
#' @param q a numerical value representing the target fdr level.
#' @param cutoff a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param seed a numerical value indicating the simulation seed (can be omit).
#' @param dig an integer representing the digit precision.
#' @param tau.tm a numerical value representing the censor level for p-values.
#' @param tau.ta a numerical value representing the censor level for p-values.
#' @param n.group.max an integer determining how many groups do we use to accelerate.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' Tm is a numerical vector consisting of primary statistics;
#' Ta is a numerical vector consisting of auxiliary statistics;
#' tm0 a numerical value representing the threshold for the primary statistics;
#' ta0 a numerical value representing the threshold for the auxiliary statistics;
#' final.fdr a numerical value representing the final estimated false discovery rate.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(fda)
#' library(qvalue)
#' # Basic Setting
#' set.seed(1234)
#' dig <- 7 # digits precision
#' n <- 1 # observation number at each point
#' m <- 120 # point size
#' q <- 0.1 # fdr level
#' const <- q # const for stablizing
#' # Generate data
#' point <- matrix(seq(0,1,length.out=m), ncol=1)
#' Dist.p <- as.matrix(dist(point))
#' I_S <- Init_Setting_1D(mu_type = "Dense",Cov_type = "Weak",
#'                        magnitude=2,mu_gen_machine="mvnorm",point=point)
#' mu <- I_S$mu
#' Sigma.eps.p <- I_S$Sigma.eps.p
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma.eps.p)
#' X <- matrix(X, nrow = n)
#' sgm <- sqrt(diag(Sigma.eps.p))
#' Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm # Generate primary statistics
#' Tm <- round(Tm,dig)
#' p.value <- 1 - pnorm(Tm) # Generate one-side p.value
#'# Run SABHA algorithm
#' tau = 0.5;eps = 0.1; TV_bd = 2
#' alpha_ADMM = 10^2; beta = 10^3; eta = 5; max_iters = 5000; converge_thr = 1e-4 # parameters for ADMM
#' ADMM_params = c(alpha_ADMM,beta,eta,max_iters,converge_thr)
#' qhat = Solve_q_TV_1dim(p.value, tau, eps, TV_bd, ADMM_params)
#' SABHA_method = function(pvals, qhat, alpha, tau){
#' pvals[pvals>tau] = Inf
#' khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
#' which(qhat*pvals<=alpha*khat/length(pvals))
#' }
#' sab.selected <- SABHA_method(p.value, qhat, q, tau)
#'
#'
#'# Run 2D algorith with 2 neighbors (2D(SA))
#' ws.sabha.fun <- function(x){1/x}
#' hh <- 2
#' hh.seq <- rep(hh,m)
#' res.2D.sabha <- TwoDwSMT(X, hh.seq=hh.seq,pis=qhat,ws.fun = ws.sabha.fun,Dist.p=Dist.p,
#'                         Sigma.eps.p=Sigma.eps.p, q=q, const = const,n.group.max=5)
#' selected.2D.sabha <- res.2D.sabha$selected
#' tm <- res.2D.sabha$tm0
#' ta <- res.2D.sabha$ta0
#' # fdp and power
#' fdp(sab.selected,mu);fdp(selected.2D.sabha,mu)
#' Pow(sab.selected,mu);Pow(selected.2D.sabha,mu)
TwoDwSMT <- function(X,
                   hh.seq,
                   pis = NULL,
                   ws.fun=NULL,
                   ws = NULL,
                   Dist.p,
                   Sigma.eps.p,
                   detect.m = "top.k",
                   q=0.1,
                   cutoff = NULL,
                   const = 0,
                   seed = 0,
                   dig=7,
                   tau.tm = 1,
                   tau.ta = 1,
                   n.group.max=5){
  n <- dim(X)[1]
  m <- dim(X)[2]
  sgm <- sqrt(diag(Sigma.eps.p))
  Tm <- apply(X,2,function(x){sum(x)/sqrt(n)})/sgm # Generate primary statistics
  Tm <- round(Tm,dig)

  # Draw necessary results from one dimensional detection
  res.1D <- OneD_Detect_w(Tm, q, pis = pis,
                          ws = ws,
                          ws.fun = ws.fun,
                          const = const,
                          tau.tm = 1)
  tm.star <- res.1D$tm.min
  max.rej <- res.1D$max.rej
  selected.1D <- res.1D$selected


  # Generate prelimaries for two dimensional detection
  Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X, Dist = Dist.p, Sigma.eps = Sigma.eps.p, detect.m = detect.m)
  Ta <- Neigh_Detect_res$Ta
  Va <- Neigh_Detect_res$Va
  VmVa.cov <- Neigh_Detect_res$VmVa.cov
  ind <- Neigh_Detect_res$ind
  mua <- Neigh_Detect_res$mua

  # Run 2dBH
  res.2D <-  Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
                                                    VmVa.cov, ind,
                                                    q,
                                                    ws = ws,
                                                    max.rej = max.rej,
                                                    pis = pis,
                                                    pws.tm.star = tm.star,
                                                    pws.ta.star = 0,
                                                    const = const,
                                                    ws.fun = ws.fun,
                                                    n.group.max = n.group.max)

  return(list(selected = res.2D$selected,
              Tm = Tm,
              Ta = Ta,
              tm0 = res.2D$tm0,
              ta0 = res.2D$ta0,
              final.fdr=res.2D$final.fdr))
}

