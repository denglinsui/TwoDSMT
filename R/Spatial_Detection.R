## Adaptive Neighbor
##
# Adapt_Neigh <- function(Dist, X, Hh){
#   m <- dim(Dist)[1]
#   # Hh: the upper bound
#
#   nei.order <- matrix(apply(Dist,1,order)[2:(Hh+1),],nrow=Hh)
#
#   order.X <- matrix(ncol = m, nrow = Hh)
#   for(i in 1:m){
#     order.X[,i] <-  X[nei.order[,i]]
#   }
#
#   signif <- matrix(apply(order.X, 2,
#                          function(x){
#                            cumsum(x)/sqrt(1:Hh)
#                          }),nrow=Hh)
#   opt.h <- apply(signif, 2, which.max)
#
#   return(opt.h)
# }



#' Neighbor Detection
#'
#' Neighbor Detection: 1. Detect the neighors (K-nearest neighbors) of each location according to the Distance matrix; 2. Calculate auxiliary test statistics;
#' 3. Calculate the variance of primary and auxiliary test statistics and covariance between them
#' @param hh an integer vector indicating the number of neighbors for each location.
#' @param X a n*m matrix representing the observations, m is the number of locations, n is the number of observations at each location.
#' @param Dist a m*m matrixes representing the distance matrix between points.
#' @param Sigma.eps a m*m matrix representing the covariance matrix for noises.
#' @param detect.m a character representing the detecting method, currently only the k-nearest neighbor, "top.k", is available.
#' @param mu mu is a vector representing the true population mean (Usually unavailable).
#'
#' @return list; Ta is a numerical vector representing the auxiliary statistics;
#' Va is a numerical vector representing the variance of auxiliary statistics;
#' VmVa.cov is a numerical vector representing the covariance between primary and auxiliary statistics;
#' ind is an integer vector representing the indexes for implementing nonparametric empirical Bayes;
#' mua is a numerical vector representing the population mean of auxiliary statistics (Usually unavailable).
#' sig_ratio is a numerical vector representing the estimated ratio of magnitude.
#' @export
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
#' hh <- 2 # neighbor number
#' hh.seq <- rep(hh,m)
#' Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X,
#'                                  Dist = Dist.p, Sigma.eps = Sigma.eps.p,
#'                                  detect.m = "top.k")
#' Ta <- Neigh_Detect_res$Ta
#' Va <- Neigh_Detect_res$Va
#' VmVa.cov <- Neigh_Detect_res$VmVa.cov
#' ind <- Neigh_Detect_res$ind
#' mua <- Neigh_Detect_res$mua
#' sig_ratio <- Neigh_Detect_res$sig_ratio
#'
#' res.2D <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, VmVa.cov, ind,
#'                                            q, max.rej,
#'                                            Tm.star = Tm.star,
#'                                            Ta.star = Inf,
#'                                            sig_ratio=sig_ratio,
#'                                            const = const)
#' selected.2D <- res.2D$selected
#' tm <- res.2D$tm0
#' ta <- res.2D$ta0
#'
#' # fdp and power
#' fdp(selected.BH,mu);fdp(selected.1D,mu);fdp(selected.2D,mu)
#' Pow(selected.BH,mu);Pow(selected.1D,mu);Pow(selected.2D,mu)
Neigh_Detect <- function(hh, X, Dist, Sigma.eps, detect.m = "top.k",mu=NULL){
  m <- dim(X)[2]
  n <- dim(X)[1]
  if(detect.m == "top.k"){
    #The neighbor is determined by the top hh nearest point
    S.neigh <- lapply(1:m,
                      function(i){
                        ind1 <- order(Dist[i,])[2:(hh[i]+1)]
                        ind1})
  }

  # if(detect.m == "rad.h"){
  #   #The neighbor is chosen based on the cutoff distance
  #   S.neigh <- lapply(1:m,
  #                     function(i){
  #                       ind1 <- which(Dist[i,]<=hh[i])
  #                       ind <- ind1[ind1!=i]
  #                       ind})
  # }

  ind = Neigh_Partition(S.neigh)
  #ind = 1:m

  Va <- sapply(1:m,
               function(i){
                 ind.Va = S.neigh[[i]]
                 sqrt(sum(Sigma.eps[ind.Va,ind.Va]))})

  Ta <- sapply(1:m,
               function(i){
                 ind.Ta = S.neigh[[i]]
                 sum(X[,ind.Ta]) /(sqrt(n))})

  sig_ratio <- hh*sqrt(diag(Sigma.eps))/Va

  Ta <- Ta/Va
  if(!is.null(mu)){
    mua <- sapply(1:m,
                  function(i){
                    ind.Ta = S.neigh[[i]]
                    sum(mu[ind.Ta]) /(sqrt(n))})
    mua <- mua/Va
  }else{
    mua <-NULL
  }

  VmVa.cov <- sapply(1:m,
                     function(i){
                       ind.Cov = S.neigh[[i]]
                       sum(Sigma.eps[i,ind.Cov])})
  VmVa.cov <- VmVa.cov/Va/sqrt(diag(Sigma.eps))
  Va <- Va/Va
  return(list(Ta = Ta,
              Va = Va,
              VmVa.cov = VmVa.cov,
              ind = ind,
              mua = mua,
              sig_ratio=sig_ratio))
}

#' Neighbor Partition
#'
#' Neighbor Partition: To make the neighbors of different location do not intersect.
#' @param S.neigh a list consists of m vector, i-th vector is the indexes of the i-th location's neighbors.
#'
#' @return an integer vector representing the indexes for implementing nonparametric empirical Bayes;
#' @export
#'
#' @examples
#' m <- 900
#' point <- matrix(seq(0,1,length.out=m), ncol=1)
#' Dist.p <- as.matrix(dist(point))
#' hh.seq <- rep(2,m)
#' S.neigh <- lapply(1:m,
#'                  function(i){
#'                  ind1 <- order(Dist.p[i,])[2:(hh.seq[i]+1)]
#'                  ind1})
#' ind <- Neigh_Partition(S.neigh)
Neigh_Partition <- function(S.neigh){
  S.len <- lapply(S.neigh,length)
  #S.len.order <- order(unlist(S.len))
  S.len.order <- 1:length(unlist(S.len))
  #S.len.order <- sample(seq(1,length(unlist(S.len)),3))

  Union <- NULL
  ind <- NULL
  for(ind.order in S.len.order){
    S.neigh.cur <- S.neigh[[ind.order]]
    if(length(intersect(S.neigh.cur,Union))==0){
      ind <- c(ind,ind.order)
      Union <- union(Union,S.neigh.cur)
    }
  }

  return(ind)
}

#Cufoff for exact: not identical to BH like
# cutoff.gen <- function(Tm, Ta, Tm.star = Inf, Ta.star = Inf){
#   # No ties are allowed!
#   index <- order(Tm, decreasing = T)
#   cutoff <- matrix(c(Inf, Inf), ncol = 2)
#   Ta.tmp <- c()
#   for(i in 1:length(Tm)){
#     index.Tm <- index[i]
#     Ta.tmp <- c(Ta.tmp,Ta[index.Tm])
#     if(Tm[index.Tm]<=Tm.star){
#       Ta.choose <- sort(Ta.tmp[Ta.tmp<=min(Ta[index.Tm],Ta.star)])
#       #index.Ta.pre <- index[1:i]
#       #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]
#
#       #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
#       #                sort(Ta[index.Ta])),
#       #              ncol = 2)
#       tmp <- matrix(c(rep(Tm[index.Tm], length(Ta.choose)),
#                       Ta.choose),
#                     ncol = 2)
#       # OneStep Optimize
#       #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
#       cutoff <- rbind(cutoff, tmp)
#     }
#   }
#
#   return(cutoff)
# }

#' Cutoff Initialization (BH based)
#'
#' Initialize the candidate cutoff set based on primary and auxiliary statistics for searching the optimal threshold.
#' @param Tm a numerical vector representing the auxiliary statistics.
#' @param Ta a numerical vector representing the primary statistics.
#' @param Tm.star a numerical value representing the threshold found by one dimensional detection based on Tm.
#' @param Ta.star a numerical value representing the threshold found by one dimensional detection based on Ta.
#'
#' @return list;
#' cutoff is a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenious;
#' ind.Tm is an integer vector representing the indexes where the first coordinate of cutoffs, Tm, changes.
#' @export
#'
#' @examples
#' library(ggplot2)
#' set.seed(1234)
#' Tm <- rnorm(20)
#' Ta1 <- rnorm(20)
#' Ta2 <- rnorm(20,Tm,0.5) # depends on Tm
#' cutoff1 <- cutoff.gen.rec(Tm,Ta1)
#' cutoff2 <- cutoff.gen.rec(Tm,Ta2)
#' # Visualize finding cutoffs
#' Type <- c(rep("Cutoff",dim(cutoff1$cutoff)[1]),rep("Data",20))
#' qplot(x=c(cutoff1$cutoff[,1],Tm),y=c(cutoff1$cutoff[,2],Ta1),
#'       colour = Type, shape=Type, xlab="Tm", ylab="Ta")+
#'       scale_shape_manual(values = c(20,4))
#' dim(cutoff1$cutoff)[1] # The number of cutoffs
#' Type <- c(rep("Cutoff",dim(cutoff2$cutoff)[1]),rep("Data",20))
#' qplot(x=c(cutoff2$cutoff[,1],Tm),y=c(cutoff2$cutoff[,2],Ta2),
#'       colour = Type, shape=Type, xlab="Tm", ylab="Ta")+
#'       scale_shape_manual(values = c(20,4))
#' dim(cutoff2$cutoff)[1] # The number of cutoffs
cutoff.gen.rec <- function(Tm, Ta, Tm.star = Inf, Ta.star = Inf){
  #=== Ties are allowed
  #=== We generate the corresponding replicates
  index <- order(Tm, decreasing = T)
  cutoff <- matrix(c(Inf, Inf), ncol = 2)
  Ta.tmp <- c()
  index.Tm.change <- c(1)
  index.Tm.change.end <- 2
  for(i in 1:length(Tm)){
    index.Tm.change[i+1] <- index.Tm.change.end
    index.Tm <- index[i]
    Ta.tmp[i] <- Ta[index.Tm]
    if(Tm[index.Tm]<=Tm.star){
      Ta.choose <- sort(Ta.tmp[Ta.tmp<=min(Ta[index.Tm],Ta.star)],
                        decreasing = T)
      #index.Ta.pre <- index[1:i]
      #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]

      #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
      #                sort(Ta[index.Ta])),
      #              ncol = 2)
      Ta.choose.len <- length(Ta.choose)
      tmp <- matrix(c(rep(Tm[index.Tm], Ta.choose.len),
                      Ta.choose),
                    ncol = 2)
      #record the position of Tm
      index.Tm.change.end <- index.Tm.change.end + Ta.choose.len
      # OneStep Optimize
      #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
      cutoff <- rbind(cutoff, tmp)
    }
  }
  index.Tm.change[length(Tm)+2] <- index.Tm.change.end

  return(list(cutoff = cutoff,
              ind.Tm = index.Tm.change))
}



#' Cutoff Initialization (weighted BH based)
#'
#' Initialize the candidate cutoff set based on weighted p-values for primary and auxiliary statistics for searching the optimal threshold.
#' @param pws.tm a numerical vector representing the weighted p-values of primary statistics.
#' @param pws.ta a numerical vector representing the weighted p-values of representing the auxiliary statistics.
#' @param pws.tm.star a numerical value representing the threshold found by one dimensional detection based on the weighted p-values of primary statistics.
#' @param pws.ta.star a numerical value representing the threshold found by one dimensional detection based on the weighted p-values of auxiliary statistics.
#'
#' @return list;
#' cutoff is a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;
#' ind.Tm is an integer vector representing the indexes where the first coordinate of cutoffs, Tm, changes.
#' @export
#'
#' @examples
#' library(ggplot2)
#' set.seed(1234)
#' Tm <- rnorm(20)
#' Ta1 <- rnorm(20)
#' Ta2 <- rnorm(20,Tm,0.5) # depends on Tm
#' pval.tm <- 1-pnorm(Tm)
#' pval.ta1 <- 1-pnorm(Ta1)
#' pval.ta2 <- 1-pnorm(Ta2)
#' cutoff1 <- cutoff.pw.gen.rec(pval.tm,pval.ta1)
#' cutoff2 <- cutoff.pw.gen.rec(pval.tm,pval.ta2)
#' # Visualize finding cutoffs
#' Type <- c(rep("Cutoff",dim(cutoff1$cutoff)[1]),rep("Data",20))
#' qplot(x=c(cutoff1$cutoff[,1],pval.tm),y=c(cutoff1$cutoff[,2],pval.ta1),
#'       colour = Type, shape=Type, xlab="pval.tm", ylab="pval.ta")+
#'       scale_shape_manual(values = c(20,4))
#' dim(cutoff1$cutoff)[1] # The number of cutoffs
#' Type <- c(rep("Cutoff",dim(cutoff2$cutoff)[1]),rep("Data",20))
#' qplot(x=c(cutoff2$cutoff[,1],pval.tm),y=c(cutoff2$cutoff[,2],pval.ta2),
#'       colour = Type, shape=Type, xlab="pval.tm", ylab="pval.ta")+
#'       scale_shape_manual(values = c(20,4))
#' dim(cutoff2$cutoff)[1] # The number of cutoffs
cutoff.pw.gen.rec <- function(pws.tm, pws.ta, pws.tm.star = 0, pws.ta.star = 0){
  #=== Ties are allowed
  #=== The cutoffs are generated adjusting for the features of weighted p-values
  #=== A large amount of weighted p-value is 1, the replicates are useless in this case.

  pws.tm.rev <- -pws.tm
  pws.ta.rev <- -pws.ta

  index <- order(pws.tm.rev, decreasing = T)
  cutoff <- matrix(c(0, 0), ncol = 2) # rejects none.
  ta.tmp <- c()
  index.tm.change <- c(1)
  index.tm.change.end <- 2
  for(i in 1:length(pws.tm.rev)){
    index.tm.change[i+1] <- index.tm.change.end
    index.tm <- index[i]
    ta.tmp[i] <- pws.ta.rev[index.tm]
    if(pws.tm.rev[index.tm]<=-pws.ta.star){
      ta.choose <- sort(ta.tmp[ta.tmp<=min(pws.ta.rev[index.tm],pws.ta.star)],
                        decreasing = T)

      #--- Remove replicates with 1 with respect to pws.ta
      if(sum(ta.choose==-1)!=0){
        ta.choose <- c(ta.choose[ta.choose>-1],-1)
      }

      #index.Ta.pre <- index[1:i]
      #index.Ta <- index[which(Ta[index.Ta.pre]<=Ta[index.Tm])]

      #tmp <- matrix(c(rep(Tm[index.Tm], length(index.Ta)),
      #                sort(Ta[index.Ta])),
      #              ncol = 2)
      ta.choose.len <- length(ta.choose)
      tmp <- matrix(c(rep(pws.tm.rev[index.tm], ta.choose.len),
                      ta.choose),
                    ncol = 2)
      #record the position of Tm
      index.tm.change.end <- index.tm.change.end + ta.choose.len
      # OneStep Optimize
      #tmp <- tmp[tmp[,1]<=Tm.star & tmp[,2]<=Ta.star,]
      cutoff <- rbind(cutoff, tmp)
    }

    #=--- Remove replicates with 1 with respect to pws.tm
    if(pws.tm.rev[index.tm]==-1) break
  }
  index.tm.change[length(index.tm.change)+1] <- index.tm.change.end

  return(list(cutoff = -cutoff,
              ind.tm = index.tm.change))
}


#' False discovery probability
#'
#' Estimated false discovery probability for one location.
#' @param tm a numerical value representing the threshold for the primary statistics.
#' @param ta a numerical value representing the threshold for the auxiliary statistics
#' @param mm a numerical vector representing the estimated support of the prior distribution for the auxiliary statistics.
#' @param normalized.prob a numerical vector representing the corresponding estimated probability of the prior distribution for the auxiliary statistics.
#' @param vVa a numerical value representing the variance of the auxiliary statistic. (Usually be one, so we omit using it here.)
#' @param vVmVa.cov a numerical value representing the covariance between the primary and auxiliary statistic. (In most case, it is the correlation between them.)
#'
#' @return a numeric value representing the estimated false discovery probability.
#' @importFrom pbivnorm pbivnorm
#' @importFrom stats pnorm
#' @importFrom stats pnorm
#' @importFrom REBayes GLmix
#' @export
#'
#' @examples
#' m <- 900
#' set.seed(1234)
#' Tm <- rnorm(m)
#' Ta <- rnorm(m)
#' mm <- REBayes::GLmix(x = Ta)
#' normalized.prob <- mm$y / sum(mm$y)
#' tm0 <- 1
#' ta0 <- 1
#' L.val <- L.cal(tm0,ta0,mm,normalized.prob,1,0)
#' # estimated fdp
#' (L.val+0.1)/sum(Tm>=tm0&Ta>=ta0)
L.cal <- function(tm,ta,
                  mm,normalized.prob,
                  vVa,vVmVa.cov){
  #=== Estimating L
  x1 <- -tm
  x2 <- tm
  y1 <- -ta- mm$x #* vVa
  y2 <- ta- mm$x #* vVa
  A1 <- pbivnorm(x = x1, y = y1, rho = vVmVa.cov)
  A2 <- pbivnorm(x = x2, y = y1, rho = vVmVa.cov)
  A3 <- pbivnorm(x = x2, y = y2, rho = vVmVa.cov)
  A4 <- pbivnorm(x = x1, y = y2, rho = vVmVa.cov)

  B1 <- pnorm(x1)
  B2 <- pnorm(x2)
  C1 <- pnorm(y1)
  C2 <- pnorm(y2)

  if(F){
    L <- sum(normalized.prob * (1 + A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4))
  }
  #One side
  L <- sum(normalized.prob *(1-B2 -C2 + A3))

  return(L)
}


#' One dimension detection for primary statistics (BH based)
#'
#‘ One Dimension Detection based on the primary statistics.
#' This is an equivalent expression for BH procedure if we set const=0 and tau.tm=1.
#' @param Tm a numerical vector representing the auxiliary statistics.
#' @param q a numerical value representing the target fdr level.
#' @param pis a numerical vector representing the probability of being null.
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param tau.tm a numerical value representing the censor level for p-values.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' tm.min is the one dimensional threshold;
#' max.rej is an interger indicating the number of selected locations.
#'
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' fdr.level <- 0.1
#' Tm <- rnorm(100)
#' res.1D.Tm <- OneD_Detect(Tm,q=fdr.level)
OneD_Detect <- function(Tm,
                        q=0.1,
                        pis = NULL,
                        const = 0,
                        tau.tm = 1){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }
  tau.Tm <- qnorm(1-tau.tm)
  p <- length(Tm)
  Tm.ind <- rank(-Tm)

  FDP <- (sum(pis)*(1-pnorm(Tm))+const)/Tm.ind

  if(sum(FDP<=q)!=0){
    threshold <- max(tau.Tm,min(Tm[FDP<=q]))
    selected <- which(Tm>=threshold)
    tm.min <- threshold
    max.rej <- length(selected)
  }else{
    selected <- NULL
    tm.min <- Inf
    max.rej <- 0
  }
  return(list(selected = selected,
              tm.min = tm.min,
              max.rej = max.rej))
}



#' One dimension detection for auxiliary statistics (BH based)
#'
#' One Dimension Detection based on the auxiliary statistics.
#' @param Ta a numerical vector representing the auxiliary statistics.
#' @param q a numerical value representing the target fdr level.
#' @param ind an integer vector representing the indexes for implementing nonparametric empirical Bayes.
#' @param pis a numerical vector representing the probability of being null.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' ta.min is the one dimensional threshold for auxiliary statistics;
#' max.rej is an interger indicating the number of selected locations.
#' @export
#'
#' @examples
#' fdr.level <- 0.1
#' Ta <- rnorm(100)
#' res.1D.Ta <- OneD_Detect_Ta(Ta,q=fdr.level,ind=1:100)
OneD_Detect_Ta <- function(Ta,
                           q=0.1,
                           ind,
                           pis = NULL){
  if(is.null(pis)){
    pis <- rep(1,length(Ta))
  }

  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  mm <- GLmix(x = eta[ind])
  normalized.prob <- mm$y / sum(mm$y)

  p <- length(Ta)
  Ta.ind <- rank(-Ta)

  FD <- pnorm(Ta)

  FD.individual <- sapply(Ta, function(ta){1 - sum(pnorm(ta- mm$x)*normalized.prob)})
  FDP <- sum(pis)*FD.individual/Ta.ind

  if(sum(FDP<=q)!=0){
    threshold <- min(Ta[FDP<=q])
    selected <- which(Ta>threshold)
    ta.min <- threshold
    max.rej <- length(selected)
  }else{
    selected <- NULL
    ta.min <- Inf
    max.rej <- 0
  }
  return(list(selected = selected,
              ta.min = ta.min,
              max.rej = max.rej))
}


#' One dimension detection for primary statistics (weighted BH based)
#'
#' One Dimension Detection (weighted) based on the primary statistics.
#' This is an equivalent expression for weighted BH procedure if we set const=0 and tau.tm=1.
#' The weights and pis can be assigned separately.
#' @param Tm a numerical vector representing the auxiliary statistics.
#' @param q a numerical value representing the target fdr level.
#' @param pis a numerical vector representing the probability of being null.
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param ws.fun a function for weighting, depending on the null proportion.
#' @param ws a numerical vector indicating the non-null weights.
#' @param tol a numerical value to control the calculation precision.
#' @param tau.tm a numerical value representing the censor level for p-values.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' tm.min is the one dimensional threshold;
#' max.rej is an interger indicating the number of selected locations.
#' @export
#'
#' @examples
#' fdr.level <- 0.1
#' Tm <- rnorm(100)
#' res.1D.w.Tm <- OneD_Detect_w(Tm,q=fdr.level)
OneD_Detect_w <- function(Tm,
                          q=0.1,
                          pis = NULL,
                          const = 0,
                          #ws.fun = function(pis){(1-pis)/(pis)},
                          ws.fun = NULL,
                          ws = NULL,
                          tol = 1e-9,
                          tau.tm = 1){
  #Initialization
  m <- length(Tm)
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }

  # Automatically calculate the weight
  # pis is given previously or is assigned by its default value:1.
  if(is.null(ws)&is.null(ws.fun)){
    ws <- rep(1,m)
  }else if(is.null(ws)&(!is.null(ws.fun))){
    nu <- 10e-5
    pis.new <- pis
    pis.new[which(pis.new<nu)] <- nu # stabilization
    pis.new[which(pis.new>1-nu)] <- 1-nu # stabilization
    ws <- ws.fun(pis.new)
  }

  nu <- 1e-5
  ws[ws<nu] <- nu # Make sure no NaN generate in the 2D process

  pv.tm <- 1-pnorm(Tm)
  pws<-pmin(rep(1,m),pv.tm/ws)
  for(tm.search in sort(pws,decreasing = T)){
    FD <- sum(pis * (1-pnorm(qnorm(1-pmin(rep(tau.tm,m),ws*tm.search)))))+const
    Discover <- sum(Tm>=qnorm(1-pmin(rep(tau.tm,m),ws*tm.search))-tol)
    #print(Discover)
    fdp.hat <- FD/max(Discover,1)
    if(fdp.hat<q) break
  }


  if(sum(fdp.hat<=q)!=0){
    selected <- which(Tm>=qnorm(1-pmin(rep(tau.tm,m),ws*tm.search))-tol)
    max.rej <- length(selected)
    tm.min <- tm.search
  }else{
    selected <- NULL
    tm.min <- 0
    max.rej <- 0
    #Equivalent Expression for LAWS or SABHA

  }
  return(list(selected = selected,
              tm.min = tm.min,
              max.rej = max.rej))
}


#' Two dimension detection (BH based)
#'
#' Searching for the optimal threshold for the two dimensional BH procedure based on the primary and auxiliary statistics.
#' @param Tm a numerical vector representing the primary statistics.
#' @param Ta a numerical vector representing the auxiliary statistics.
#' @param Va a numerical vector representing the variance of auxiliary statistics.
#' @param VmVa.cov a numerical vector representing the covariance between primary and auxiliary statistics.
#' @param ind an integer vector representing the indexes for implementing nonparametric empirical Bayes.
#' @param q a numerical value representing the target fdr level.
#' @param max.rej an integer indicating the current maximum rejection number that controls fdr under q.
#' @param pis a numerical vector representing the probability of being null.
#' @param cutoff a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;
#' @param Tm.star a numerical value representing the threshold found by one dimensional detection based on the primary statistics.
#' @param Ta.star a numerical value representing the threshold found by one dimensional detection based on the auxiliary statistics.
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param seed a numerical value indicating the simulation seed (can be omit).
#' @param dig an integer representing the digit precision.
#' @param tau.tm a numerical value representing the censor level for p-values.
#' @param tau.ta a numerical value representing the censor level for p-values.
#' @param mua a numerical vector representing the population mean of auxiliary statistics (Usually unavailable).
#' @param is.quick.stop a bool value representing whether we should stop early.
#' @param stop.step an integer indicating the maximal step for no valuable finding(estimate fdr<target fdr level).
#' @param is.fullset a bool value representing whether we take all partitions of full set for implementing NPEB.
#' @param is.fixRectShape a bool value representing whether the shape of rectangle is fixed
#' @param EmpMethod a character indicating which method we choose to estimate the empirical distribution.
#' @param sig_ratio a numerical vector representing the estimated ratio of magnitude.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
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
#' hh <- 2 # neighbor number
#' hh.seq <- rep(hh,m)
#' Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X,
#'                                  Dist = Dist.p, Sigma.eps = Sigma.eps.p,
#'                                  detect.m = "top.k")
#' Ta <- Neigh_Detect_res$Ta
#' Va <- Neigh_Detect_res$Va
#' VmVa.cov <- Neigh_Detect_res$VmVa.cov
#' ind <- Neigh_Detect_res$ind
#' mua <- Neigh_Detect_res$mua
#'
#' res.2D <- Spatial_Detect_exact_grp_BH_down(Tm, Ta, Va, VmVa.cov, ind,
#'                                            q, max.rej,
#'                                            Tm.star = Tm.star,
#'                                            Ta.star = Inf,
#'                                            const = const)
#' selected.2D <- res.2D$selected
#' tm <- res.2D$tm0
#' ta <- res.2D$ta0
#'
#' # fdp and power
#' fdp(selected.BH,mu);fdp(selected.1D,mu);fdp(selected.2D,mu)
#' Pow(selected.BH,mu);Pow(selected.1D,mu);Pow(selected.2D,mu)
Spatial_Detect_exact_grp_BH_down <- function(Tm, Ta, Va, VmVa.cov, ind,
                                             q=0.1, max.rej = NULL,
                                             pis = NULL, cutoff = NULL,
                                             Tm.star = Inf, Ta.star = Inf,
                                             const = 0,
                                             seed = 0,
                                             dig=7,
                                             tau.tm = 1,
                                             tau.ta = 1,
                                             mua = NULL,
                                             is.quick.stop=F,
                                             stop.step=10,
                                             is.fullset = F,
                                             is.fixRectShape = F,
                                             sig_ratio = NULL,
                                             EmpMethod = c("NPEB","DeCor")){
  if(is.null(pis)){
    pis <- rep(1,length(Tm))
  }

  if(is.null(sig_ratio)){
    sig_ratio <- rep(1,length(Tm))
  }

  EmpMethod <- match.arg(EmpMethod)

  if(is.fixRectShape){
    if(EmpMethod=="NPEB"){
      avrg_c <- mean(sig_ratio)
    }else{
      avrg_c <- mean((sig_ratio-VmVa.cov)/(1-VmVa.cov^2))
    }
  }


  ind.step <- 0
  #Initialization
  m <- length(Tm)
  eta <- Ta
  #=== Perform Non-parametric Emprical Bayesian
  if(EmpMethod == "NPEB"){
    if(is.fullset){
      ## We just consider one dimensional case
      # mm <- GLmix(x = eta[ind])
      # normalized.prob <- mm$y / sum(mm$y)
      # num.par <- 1
      # for(ind.par in (ind[1]+1):(ind[2]-1)){
      #   ind.cur <- ind+ind.par-1
      #   ind.cur <- ind.cur[ind.cur<=m] # avoid out-of-subscript
      #
      #   mm.tmp <- GLmix(x = eta[ind.cur])
      #   mm$x <- c(mm$x,mm.tmp$x)
      #   normalized.prob <- c(normalized.prob,mm.tmp$y / sum(mm.tmp$y))
      #   num.par <- num.par+1
      # }
      # normalized.prob <- normalized.prob/num.par

      mm <- GLmix(x = eta)
      normalized.prob <- mm$y / sum(mm$y)
    }else{
      #=== Perform Non-parametric Emprical Bayesian
      mm <- GLmix(x = eta[ind])
      normalized.prob <- mm$y / sum(mm$y)
    }

    if(!is.null(mua)){
      mm.pre <- table(mua)
      mm$x <- as.numeric(names(mm.pre))
      mm$y <- mm.pre/sum(mm.pre)
      normalized.prob <- mm$y / sum(mm$y)
    }
    #mm$x <- rep(0,length(mm$x))
    #=== Estimate false discovery via NPEB
    fdr.est <- function (tm, ta, NP.max, etype){
      p <- length(Tm)
      L <- numeric(p)
      NP <- sum(Tm >= tm & Ta >= ta)
      NP <- ifelse(is.na(NP), 0, NP)

      if (NP == 0){
        FDP <- 0
        FD <- 0
      }else if(NP<NP.max){
        FDP <- NA
        FD <- NA
      }else{
        #group_by (Ta,VmVa.cov)
        grp.val <- unique(cbind(Va,VmVa.cov))
        for(j in 1:nrow(grp.val)){
          vVa <- grp.val[j,1]
          vVmVa.cov <- grp.val[j,2]
          row.ind <- which(Va==vVa & VmVa.cov==vVmVa.cov)
          L[row.ind] <- L.cal(tm,ta,
                              mm,normalized.prob,
                              vVa,vVmVa.cov)
        }

        if (etype == 'FDR') {
          # When pi0.est, this is simply summation
          FD <- sum(pis * L) + const
          FDP <- FD / NP
          #print(paste(FD,NP,FDP))
        }
      }
      return(list(FDP = FDP,
                  NP = NP,
                  FD = FD))
    }
  }else if (EmpMethod == "DeCor"){
    Ta <- Ta - VmVa.cov *Tm # Create Decorrelate Statistics
    #=== Estimate false discovery via decor statistics
    fdr.est <- function (tm, ta, NP.max, etype){
      p <- length(Tm)
      L <- numeric(p)
      NP <- sum(Tm >= tm & Ta >= ta)
      NP <- ifelse(is.na(NP), 0, NP)

      if (NP == 0){
        FDP <- 0
        FD <- 0
      }else if(NP<NP.max){
        FDP <- NA
        FD <- NA
      }else{
        L <- (1-pnorm(tm))*as.integer(Ta>=ta)

        if (etype == 'FDR') {
          # When pi0.est, this is simply summation
          FD <- sum(pis * L) + const
          FDP <- FD / NP
          #print(paste(FD,NP,FDP))
        }
      }
      return(list(FDP = FDP,
                  NP = NP,
                  FD = FD))
    }
  }
  pv.tm <- 1-pnorm(Tm)
  pv.ta <- 1-pnorm(Ta)

  #=== Initialize the searching grid
  if(is.null(cutoff)){
    Tm.cutoff <- Tm[pv.tm<=tau.tm & pv.ta<=tau.ta] # Add tau-censoring
    Ta.cutoff <- Ta[pv.tm<=tau.tm & pv.ta<=tau.ta]
    if(is.fixRectShape){
      # Tm.cutoff from 1dDetect can be insuitable for this method
      # Therefore, we need to adjust it
      Tm.star <- Inf
      max.rej <- 0
    }
    cutoff.rec <- cutoff.gen.rec(Tm.cutoff,Ta.cutoff,Tm.star,Ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.Tm <- cutoff.rec$ind.Tm

    if(is.fixRectShape){
      cutoff[,2] <- cutoff[,1] * avrg_c
    }
  }

  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej

  index.Tm <- max.rej+1
  m_prime <- length(ind.Tm)-1
  while((index.Tm <= m_prime)&((!is.quick.stop)|ind.step<=stop.step)){
    #print(index.Tm)
    i.up <- ind.Tm[index.Tm+1]-1
    cur.ind.Tm <- ind.Tm[index.Tm]
    i.down <- cur.ind.Tm

    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]

    NP.down <- sum(Tm >= tm.down & Ta >= ta.down)
    i.down <- i.down + max(0, NP.max - NP.down)

    if(is.fixRectShape){
      i.down <- i.up # For each tm, we only consider one ta
    }

    while(i.up >= i.down){
      #print(paste(i.up,i.down))

      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]

      obj <- fdr.est(tm.down, ta.down, NP.max, etype="FDR")


      #print(c(tm.down,ta.down,obj$FDP))

      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Init ind.step
        ind.step <- 0
        # Save to find minimum FDP
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #print(paste(tm.down,ta.down,obj$FD,obj$NP,obj$FDP))
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down

        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        #print(paste(tm.down,obj$FD,obj$NP,obj$FDP))

        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        i.down <- i.down + max(1, min.REJ-NP.down)
        #i.down <- i.down + min.REJ - NP.down
      }
    }
    ind.step <- ind.step+1
    index.Tm <- index.Tm + 1
  }

  if(is.null(FDP)){
    selected <- NULL
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  #print(paste0("seed",seed,"start"))
  pos <- Tm >= tm0 & Ta >= ta0
  #print(paste0("seed",seed,"end"))

  selected <- which(pos)

  final.obj <- fdr.est(tm0, ta0, NP.max, etype="FDR")
  final.fdr <- final.obj$FDP
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0,
              final.fdr = final.fdr))
}




#' Two dimension detection for primary statistics (weighted BH based)
#'
#' Searching for the optimal threshold for the two dimensional weighted BH procedure based on primary and auxiliary statistics.
#' Practically, we reorder p-values according to the non-null weights ws.
#' @param Tm a numerical vector representing the primary statistics.
#' @param Ta a numerical vector representing the auxiliary statistics.
#' @param Va a numerical vector representing the variance of auxiliary statistics.
#' @param VmVa.cov a numerical vector representing the covariance between primary and auxiliary statistics.
#' @param ind an integer vector representing the indexes for implementing nonparametric empirical Bayes.
#' @param q a numerical value representing the target fdr level.
#' @param max.rej an integer indicating the current maximum rejection number that controls fdr under q.
#' @param pis a numerical vector representing the probability of being null.
#' @param ws a numerical vector indicating the non-null weights.
#' @param cutoff a matrix with two columns representing the candidate cutoff set and it should be noted that the cutoffs are ordered for searching convenience;
#' @param pws.tm.star a numerical value representing the threshold found by one dimensional detection with non-null weights based on the primary statistics.
#' @param pws.ta.star a numerical value representing the threshold found by one dimensional detection with non-null weights based on the auxiliary statistics.
#' @param const a numerical value ensuring the selection stability, usually be target fdr level.
#' @param seed a numerical value indicating the simulation seed (can be omit).
#' @param dig an integer representing the digit precision.
#' @param ws.fun a function for weighting, depending on the null proportion.
#' @param n.group.max an integer determining how many groups do we use to accelerate.
#' @param is.exact.group a logical representing whether the group is exact.
#' @param tol a numerical value to control the calculation precision.
#' @param tau.tm a numerical value representing the censor level for p-values.
#' @param tau.ta a numerical value representing the censor level for p-values.
#' @param adj.up a numerical value to ensure the precision during calculation, it should be below 1;
#' @param adj.ud a numerical value to ensure the precision during calculation, it should be above 1;
#' @param mua a numerical vector representing the population mean of auxiliary statistics (Usually unavailable).
#' @param is.quick.stop a bool value representing whether we should stop early.
#' @param stop.step an integer indicating the maximal step for no valuable finding(estimate fdr<target fdr level).
#' @param is.fullset a bool value representing whether we take all partitons of full set for implementing NPEB
#' @param EmpMethod a character indicating which method we choose to estimate the empirical distribution.
#'
#' @return list;
#' selected is an integer vector consisting of selected locations;
#' tm0 a numerical value representing the threshold for the weighted p-value of primary statistics;
#' ta0 a numerical value representing the threshold for the weighted p-value of auxiliary statistics;
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
#'# Run 1D algorithm (SABHA)
#' ws.sabha.fun <- function(x){1/x}
#' res.1D.sabha <- OneD_Detect_w(Tm, q, pis = qhat,
#'                               ws.fun = ws.sabha.fun,
#'                               const = const,
#'                               tau.tm = 1)
#' tm.star.sabha <- res.1D.sabha$tm.min
#' max.rej.sabha <- res.1D.sabha$max.rej
#' selected.1D.sabha <- res.1D.sabha$selected
#'
#'# Run 2D algorith with 2 neighbors (2D(SA))
#' hh <- 2 # neighbor number
#' hh.seq <- rep(hh,m)
#' Neigh_Detect_res <- Neigh_Detect(hh = hh.seq, X = X,
#'                                  Dist = Dist.p, Sigma.eps = Sigma.eps.p,
#'                                  detect.m = "top.k")
#' Ta <- Neigh_Detect_res$Ta
#' Va <- Neigh_Detect_res$Va
#' VmVa.cov <- Neigh_Detect_res$VmVa.cov
#' ind <- Neigh_Detect_res$ind
#' mua <- Neigh_Detect_res$mua
#'
#' res.2D.sabha <- Spatial_Detect_exact_BH_down_reTm_reTa(Tm, Ta, Va,
#'                                                        VmVa.cov, ind,
#'                                                        q,
#'                                                        max.rej = max.rej.sabha,
#'                                                        pis = qhat,
#'                                                        pws.tm.star = tm.star.sabha,
#'                                                        pws.ta.star = 0,
#'                                                        const = const,
#'                                                        ws.fun = ws.sabha.fun,
#'                                                        n.group.max = 5)
#' selected.2D.sabha <- res.2D.sabha$selected
#' tm <- res.2D.sabha$tm0
#' ta <- res.2D.sabha$ta0
#' # fdp and power
#' fdp(sab.selected,mu);fdp(selected.1D.sabha,mu);fdp(selected.2D.sabha,mu)
#' Pow(sab.selected,mu);Pow(selected.1D.sabha,mu);Pow(selected.2D.sabha,mu)
Spatial_Detect_exact_BH_down_reTm_reTa <- function(Tm, Ta, Va, VmVa.cov, ind,
                                                   q=0.1, max.rej = NULL,
                                                   pis = NULL, # The prior probability of being null.
                                                   ws = NULL, # weights for reordering.
                                                   cutoff = NULL,
                                                   pws.tm.star = 0, pws.ta.star = 0,
                                                   const = 0,
                                                   seed = 0,
                                                   dig=7,
                                                   ws.fun = function(pis){(1-pis)/pis},
                                                   n.group.max = 10,
                                                   is.exact.group = F,
                                                   tol=1e-9,
                                                   tau.tm=1,
                                                   tau.ta=1,
                                                   adj.up = 0.9,adj.ud = 1.1,
                                                   mua = NULL,
                                                   is.quick.stop=F,
                                                   stop.step=10,
                                                   is.fullset=F,
                                                   EmpMethod = c("NPEB","DeCor")){
  if(tau.tm<0 | tau.tm>1|tau.ta<0 | tau.ta>1){
    stop("Please check the input: tau must between 0 and 1.")
  }

  eps.in <- 10^{-dig}

  #==== Initialization
  m <- length(Tm)
  eta <- Ta
  ind.step <- 0

  EmpMethod <- match.arg(EmpMethod)

  #---- Initialize pis and ws
  if(is.null(pis)){
    pis <- rep(1,m)
  }

  if(is.null(ws)&is.null(ws.fun)){
    ws <- rep(1,m) # no weight
  }else if(is.null(ws)&(!is.null(ws.fun))){
    nu <- 10e-5
    pis.new <- pis
    pis.new[which(pis.new<nu)] <- nu # stabilization
    pis.new[which(pis.new>1-nu)] <- 1-nu # stabilization
    ws <- ws.fun(pis.new)
  }

  nu <- 1e-5
  ws[ws<nu] <- nu # Make sure no NaN generate in the process

  if(EmpMethod == "NPEB"){
    #=== Perform Non-parametric Emprical Bayesian
    if(is.fullset){
      # ## We just consider one dimensional case
      # mm <- GLmix(x = eta[ind])
      # normalized.prob <- mm$y / sum(mm$y)
      # num.par <- 1
      # for(ind.par in (ind[1]+1):(ind[2]-1)){
      #   ind.cur <- ind+ind.par-1
      #   ind.cur <- ind.cur[ind.cur<=m] # avoid out-of-subscript
      #
      #   mm.tmp <- GLmix(x = eta[ind.cur])
      #   mm$x <- c(mm$x,mm.tmp$x)
      #   normalized.prob <- c(normalized.prob,mm.tmp$y / sum(mm.tmp$y))
      #   num.par <- num.par+1
      # }
      # normalized.prob <- normalized.prob/num.par
      mm <- GLmix(x = eta)
      normalized.prob <- mm$y / sum(mm$y)
    }else{
      #=== Perform Non-parametric Emprical Bayesian
      mm <- GLmix(x = eta[ind])
      normalized.prob <- mm$y / sum(mm$y)
    }

    if(!is.null(mua)){
      mm.pre <- table(mua)
      mm$x <- as.numeric(names(mm.pre))
      mm$y <- mm.pre/sum(mm.pre)
      normalized.prob <- mm$y / sum(mm$y)
    }
    #mm$x <- rep(0,length(mm$x))
    #=== Estimate false discovery
    #=== tm,ta is the threshold with respect to the weighted p-values! Not from tm,ta.
    #=== group.value and group.ind is to accelarate the calculation
    fdr.est <- function (tm, ta, ws, NP.max =0, etype="FDR",
                         q = 1,# Target fdr value
                         Grp_Val_s = NULL,
                         Grp_Ind_s = NULL,
                         Grp_Val_s_up = NULL,
                         #ws.up_est = NULL,
                         is.exact.group = F,
                         tau.tm = 1,
                         tau.ta = 1,
                         ud_est.FD.val = 0,
                         up_est.FD.val = 0,
                         FD.ctl = NULL,
                         adj.up = 0.9,
                         adj.ud = 1.1,
                         need.exact = FALSE){

      p <- length(Tm)
      L <- numeric(p)
      L.up_est <- numeric(p)
      L.ud_est<- numeric(p)

      #--- Corresponding Rejection
      NP <- sum((Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm))-tol) &
                  (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta))-tol))

      #--- The index is for accelarating the calculation
      #--- In this way, we could calculate some false discovery probability in one dimension
      ind.tmp.tm <- ws*tm>=tau.tm
      ind.tmp.ta <- ws*ta>=tau.ta
      ind.tmp.NaN <- (ws*tm<=1e-14)|(ws*ta<=1e-14) # Avoid NaN
      #NP <- sum(Tm >= tm & Ta >= ta)
      #NP <- ifelse(is.na(NP), 0, NP)

      if(is.null(FD.ctl)){
        FD.ctl <- NP*q
      }

      if (NP == 0){
        FDP <- 0
        FD <- 0
      }else if(NP<NP.max){
        FDP <- NA
        FD <- NA
      }else{
        ind.part0 <- which(ind.tmp.NaN) #This part will induce NaN
        ind.partm <- which(ind.tmp.tm&ind.tmp.ta&(!ind.tmp.NaN))
        ind.parta <- which(!ind.tmp.tm&(ind.tmp.ta)&(!ind.tmp.NaN))
        ind.part3 <- which((ind.tmp.tm)&(!ind.tmp.ta)&(!ind.tmp.NaN))
        ind.part4 <- which((!ind.tmp.tm)&(!ind.tmp.ta)&(!ind.tmp.NaN))

        # Common part
        L[ind.part0] <- ws[ind.part0]*tm # A more conservative estimation
        L[ind.partm] <- 1
        #L[ind.parta] <- tau.tm
        L[ind.parta] <- (1-pnorm(qnorm(1-pmin(rep(tau.tm,length(ind.parta)),ws[ind.parta]*tm))))

        L.part0 <- sum(pis[ind.part0]*L[ind.part0])
        L.partm <- sum(pis[ind.partm]*L[ind.partm])
        L.parta <- sum(pis[ind.parta]*L[ind.parta])
        L.part0.ud_est <- L.part0.up_est <- L.part0#
        L.partm.ud_est <- L.partm.up_est <- L.partm#
        L.parta.ud_est <- L.parta.up_est <- L.parta#
        #=== With Group
        # If the group is exact, we only need the first part!
        if(!(is.null(Grp_Val_s))&!(is.null(Grp_Ind_s))){
          for(i in (1:nrow(Grp_Val_s))){
            ind.tmp <- Grp_Ind_s[[i]]
            val.ud_est <- Grp_Val_s[i,]; val.up_est <- Grp_Val_s_up[i,]#Vm,VmVa,ws
            ws.ud_est <- val.ud_est[3]; ws.up_est <-val.up_est[3]
            interct.ind3 <- intersect(ind.tmp,ind.part3)
            if(length(interct.ind3)!=0){
              ind.part3.tmp <- intersect(ind.tmp,ind.part3)
              L.ud_est[interct.ind3] <- 1-sum(normalized.prob * pnorm(qnorm(1-ws.ud_est*ta),
                                                                      mean = mm$x))
              L.up_est[interct.ind3] <- 1-sum(normalized.prob * pnorm(qnorm(1-min(1,ws.up_est*ta)),
                                                                      mean = mm$x))
            }
          }
          L.part3.ud_est <- sum(pis[ind.part3]*L.ud_est[ind.part3])
          L.part3.up_est <- sum(pis[ind.part3]*L.up_est[ind.part3])

          for(i in (1:nrow(Grp_Val_s))){
            ind.tmp <- Grp_Ind_s[[i]]
            val.ud_est <- Grp_Val_s[i,]; val.up_est <- Grp_Val_s_up[i,]#Vm,VmVa,ws
            Va.ud_est <- val.ud_est[1]; VmVa.cov.ud_est <- val.ud_est[2]; ws.ud_est <- val.ud_est[3]
            Va.up_est <- val.up_est[1]; VmVa.cov.up_est <- val.up_est[2]; ws.up_est <- val.up_est[3]
            #print(ws.tmp*tm)
            #print(ws.tmp*ta)
            #print(ws.tmp*ta)
            interct.ind4 <- intersect(ind.tmp,ind.part4)
            if(length(interct.ind4)!=0){
              #print(paste0(ws.ud_est*tm,ws.ud_est*ta))
              if((ws.ud_est*tm<=1e-14)|(ws.ud_est*ta<=1e-14)){
                L.ud_est[interct.ind4] <- ws.ud_est*tm
              }else{
                L.ud_est[interct.ind4] <- L.cal(qnorm(1-ws.ud_est*tm),qnorm(1-ws.ud_est*ta),
                                                mm,normalized.prob,
                                                Va.ud_est,VmVa.cov.ud_est)
              }
              L.up_est[interct.ind4] <- L.cal(qnorm(1-min(1-1e-14,ws.up_est*tm)),
                                              qnorm(1-min(1-1e-14,ws.up_est*ta)),
                                              mm,normalized.prob,Va.up_est,VmVa.cov.up_est)
            }
          }
          L.part4.ud_est <- sum(pis[ind.part4]*L.ud_est[ind.part4])
          L.part4.up_est <- sum(pis[ind.part4]*L.up_est[ind.part4])

          #-- A "Conservative" Estimation for FDP
          FD.ud_est.org <- L.part0+L.partm+L.parta+L.part3.ud_est+L.part4.ud_est  + const
          FD.ud_est <- FD.ud_est.org + ud_est.FD.val
          FD.up_est.org <- L.part0+L.partm+L.parta+L.part3.up_est+L.part4.up_est + const
          FD.up_est <- FD.up_est.org - up_est.FD.val
          FDP.ud_est <- FD.ud_est/NP; FDP.up_est <- FD.up_est/NP
          #print(c(FDP.ud_est,FDP.up_est))
        }else{
          FDP.ud_est = 0 # To make sure the calculation goes to the second part.
          FDP.up_est = 1
        }

        #=== Without Group
        # If the conservative estimation indicate this threshold is possible to be active
        # We further explore it
        # adj.up = 0.9;adj.ud = 1.1
        if(is.exact.group){
          FDP <- FDP.ud_est
          FD <- FD.ud_est
        }else if((FDP.up_est>q*adj.up & FDP.ud_est<=q*adj.ud) | (need.exact == TRUE)){
          #if(FDP.ud_est<=q*adj.ud){
          #
          #if(FDP.ud_est<=q & !is.exact.group){
          #if(T){ # For Test
          #-- Va=1
          for(j in (1:p)[ind.part3]){
            L[j] <- 1-sum(normalized.prob * pnorm(qnorm(1-ws[j]*ta),
                                                  mean = mm$x))
          }
          L.part3 <- sum(pis[ind.part3]*L[ind.part3])
          #No longer to be grouped
          for(j in (1:p)[ind.part4]){
            L[j] <- L.cal(qnorm(1-ws[j]*tm),qnorm(1-ws[j]*ta),
                          mm,normalized.prob,
                          Va[j],VmVa.cov[j])
          }
          L.part4 <- sum(pis[ind.part4]*L[ind.part4])

          FD <- L.part0+L.partm+L.parta+L.part3+L.part4 + const
          FDP <- FD/NP

          ud_est.FD.val <- FD-FD.ud_est.org # Update the approximate diff.
          up_est.FD.val <- FD.up_est.org - FD # Update the approximate diff.
          #print(c(FD.ud_est.org/NP,FDP.ud_est,FDP,FDP.up_est,FD.up_est.org/NP))
          #print(c(FD.ud_est.org,FD.ud_est,FD,FD.up_est,FD.up_est.org))
        }else if(FDP.up_est<=q*adj.up){
          FDP <- FDP.up_est
          FD <- FD.up_est
        }else{
          FDP <- FDP.ud_est
          FD <- FD.ud_est
        }
        #if (etype == 'FDR') {
        # When pi0.est, this is simply summation
        #FD <- sum(pis * L) + const
        #FDP <- FD / NP
        #print(paste(FDP,p*(1-B2)/NP))
        #}
      }
      return(list(FDP = FDP,
                  NP = NP,
                  FD = FD,
                  ud_est.FD.val = ud_est.FD.val,
                  up_est.FD.val = up_est.FD.val))
    }

  }else if(EmpMethod == "DeCor"){
    Ta <- Ta - VmVa.cov* Tm
    #=== Estimate false discovery
    fdr.est <- function (tm, ta, ws, NP.max =0, etype="FDR",
                         q = 1,# Target fdr value
                         Grp_Val_s = NULL,
                         Grp_Ind_s = NULL,
                         Grp_Val_s_up = NULL,
                         #ws.up_est = NULL,
                         is.exact.group = F,
                         tau.tm = 1,
                         tau.ta = 1,
                         ud_est.FD.val = 0,
                         up_est.FD.val = 0,
                         FD.ctl = NULL,
                         adj.up = 0.9,
                         adj.ud = 1.1,
                         need.exact = FALSE){

      p <- length(Tm)
      L <- numeric(p)
      L.up_est <- numeric(p)
      L.ud_est<- numeric(p)

      #--- Corresponding Rejection
      NP <- sum((Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm))-tol) &
                  (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta))-tol))

      #--- The index is for accelarating the calculation
      #--- In this way, we could calculate some false discovery probability in one dimension

      if (NP == 0){
        FDP <- 0
        FD <- 0
      }else if(NP<NP.max){
        FDP <- NA
        FD <- NA
      }else{
        # Obtain the expected false discoveries
        L <- (1-pnorm(qnorm(1-pmin(rep(tau.tm,m),ws*tm))-tol)) * as.integer(Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta))-tol)

        if (etype == 'FDR') {
          # When pi0.est, this is simply summation
          FD <- sum(pis * L) + const
          FDP <- FD / NP
          ud_est.FD.val <- up_est.FD.val <- FD
          #print(paste(FD,NP,FDP))
        }
      }
      return(list(FDP = FDP,
                  NP = NP,
                  FD = FD,
                  ud_est.FD.val = ud_est.FD.val,
                  up_est.FD.val = up_est.FD.val))
    }

  }

  #=== Initialize reordered p-values
  #--- weighted p value for Tm
  pv.tm <- 1-pnorm(Tm)
  #pws.tm <- pmin(rep(1,m),pv.tm/ws) #pv.tm/ws
  pws.tm <- pv.tm/ws
  #--- weighted p value for Ta
  pv.ta <- 1-pnorm(Ta)
  #pws.ta<-  pmin(rep(1,m),pv.ta/ws)#pv.ta/ws
  pws.ta<-  pv.ta/ws

  #=== Initialize the searching grid
  #--- The cutoff is based on the weighted p-value
  if(is.null(cutoff)){
    pws.tm.cutoff <- pws.tm[pv.tm<=tau.tm & pv.ta<=tau.ta]
    pws.ta.cutoff <- pws.ta[pv.tm<=tau.tm & pv.ta<=tau.ta]
    cutoff.rec <- cutoff.pw.gen.rec(pws.tm.cutoff,pws.ta.cutoff,
                                    pws.tm.star,pws.ta.star)
    cutoff <- cutoff.rec$cutoff
    ind.pws.tm <- cutoff.rec$ind.tm
  }

  #=== Initilize the group
  #print("init")
  res.hg <- Hard_Group(Org_Value = cbind(Va,VmVa.cov,ws),n.group.max = n.group.max)
  Grp_Val_s <- res.hg$Grp_Val_summarize
  Grp_Val_s_up <- res.hg$Grp_Val_up_summarize
  Grp_Ind_s <- res.hg$Grp_Ind_summarize
  #print("end")
  #ws.up_est <- res.hg$Grp_Value
  if(EmpMethod == "NPEB"){
    is.exact.group <- res.hg$is.exact.group
  }else if(EmpMethod == "DeCor"){
    is.exact.group <- T
  }
  #print(is.exact.group)
  #res.hg_up <- Hard_Group_up(Org_Value = cbind(Va,VmVa.cov,ws),n.group.max = n.group.max)
  #Grp_Val_s_up <- res.hg_up$Grp_Val_summarize
  #Grp_Ind_s_up <- res.hg_up$Grp_Ind_summarize
  #is.exact.group_up <- res.hg_up$is.exact.group

  #=== Searching for the optimal threshold
  FDP <- NULL
  NP <- NULL
  #tmta <- NULL
  tm.cand.set <- NULL
  ta.cand.set <- NULL
  NP.max <- max.rej

  index.pws.tm<- max.rej+1

  #-- Adjust for the weighted p-value
  m_prime <- length(ind.pws.tm)-1
  FD.record <- matrix(numeric(0),ncol=2)
  while((index.pws.tm <= m_prime)&&((!is.quick.stop)|ind.step<=stop.step)){
    #print(index.pws.tm)
    i.up <- ind.pws.tm[index.pws.tm+1]-1
    cur.ind.pws.tm <- ind.pws.tm[index.pws.tm]
    i.down <- cur.ind.pws.tm

    tm.up <- cutoff[i.up,1]
    ta.up <- cutoff[i.up,2]

    # We only consider the rej > current max rej
    tm.down <- cutoff[i.down,1]
    ta.down <- cutoff[i.down,2]

    NP.down <- sum((Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm.down))-tol) &
                     (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta.down))-tol))
    i.down <- i.down + max(0, NP.max - NP.down)
    ud_est.FD.val <- 0 # Initialize the estimated error
    up_est.FD.val <- 0
    while(i.up >= i.down){

      #=== Search Down
      # If not stand, we change i.down to accelerate
      tm.down <- cutoff[i.down,1]
      ta.down <- cutoff[i.down,2]
      obj <- fdr.est(tm = tm.down,ta = ta.down,ws = ws,
                     NP.max = NP.max, etype="FDR",q = q,
                     Grp_Val_s = Grp_Val_s, Grp_Ind_s = Grp_Ind_s,
                     Grp_Val_s_up = Grp_Val_s_up,
                     is.exact.group = is.exact.group,
                     tau.tm = tau.tm,
                     tau.ta = tau.ta,
                     adj.up = adj.up,adj.ud = adj.ud,
                     ud_est.FD.val = ud_est.FD.val,
                     up_est.FD.val = up_est.FD.val)
      ud_est.FD.val <- obj$ud_est.FD.val
      up_est.FD.val <- obj$up_est.FD.val
      FD.record <- rbind(FD.record,
                         c(ud_est.FD.val, up_est.FD.val))
      #  print(paste(i.up,i.down,obj$FDP))
      if(is.na(obj$FDP)){
        # In case of ties
        i.down <- i.down + max(0, NP.max - obj$NP)
      } else if(obj$FDP<=q){
        # Save to find minimum FDP
        ind.step <- 0
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        tm.cand.set <- c(tm.cand.set, tm.down)
        ta.cand.set <- c(ta.cand.set, ta.down)
        #tmta <- c(tmta, paste(tm.down, ta.down))
        NP.max <- obj$NP
        tm0 <- tm.down
        ta0 <- ta.down

        i.down <- i.down +1
      } else{
        # obj$FDP>q
        # Actually, we could alse use this strategy on the x-axis to accelarate
        # We find a way to accelerate the calculation
        FD.down <- obj$FD
        NP.down <- obj$NP
        min.REJ <- ceiling(FD.down/q)
        # Actually, min.REJ>NP.down
        # Current rejecting number is not enough
        #i.down <- i.down + max(0, min.REJ-NP.down)
        i.down.tmp <- i.down + max(min.REJ - NP.down,1)
        i.down <- ifelse((ta.down!=1)&(ta.up==1)&(i.down.tmp>i.up),
                         i.up,i.down.tmp)
      }
    }
    ind.step <- ind.step+1
    index.pws.tm <- index.pws.tm + 1
  }

  if(sum(FD.record[,1]<0 |  FD.record[,2]<0)>0){
    warning_mess <- paste0("seed:",seed,
                           " The surrogate FDP(ud) and FDP(up) are not exact upper and under estimator of FDP with at least ",
                           "Misspecified Under:",sum(FD.record[,1]<0),
                           " and Misspecified Upper:",sum(FD.record[,2]<0)
    )
    warning(warning_mess)
  }

  if(is.null(FDP)){
    selected <- NULL
    tm0 <- pws.tm.star # Conservative, to make the program run temproally.
    # For most methods, guaranteed by the one dimensional approach, we don't reach this step;
    # For LAWS, things becomes different because of the numerical digits
    ta0 <- max(1/ws)+tol
  }else{
    NP.cand <- NP[NP==max(NP)]
    FDP.cand <- FDP[NP==max(NP)]
    #tmta.cand <- tmta[NP==max(NP)]
    tm.cand <- tm.cand.set[NP==max(NP)]
    ta.cand <- ta.cand.set[NP==max(NP)]
    tm0 <- tm.cand[which.min(FDP.cand)]
    ta0 <- ta.cand[which.min(FDP.cand)]
    #tmta <- unique(tmta.cand[which.min(FDP.cand)])
    #tm0 <- as.numeric(unlist(strsplit(tmta,split=' '))[1])
    #ta0 <- as.numeric(unlist(strsplit(tmta,split=' '))[2])
  }
  #print(paste0("seed",seed,"start"))
  ind.tmp <- which(ws*tm0>1)
  #indd <- 3;tm0 <- tm.cand[indd];tm0 <- tm.cand[indd]
  # Report Selected point
  pos <- (Tm>= qnorm(1-pmin(rep(tau.tm,m),ws*tm0))-tol) &
    (Ta>=qnorm(1-pmin(rep(tau.ta,m),ws*ta0))-tol)
  #print(paste0("seed",seed,"end"))
  selected <- which(pos)

  final.obj <- fdr.est(tm = tm0,ta = ta0,ws = ws, etype="FDR",q = q,
                       Grp_Val_s = Grp_Val_s, Grp_Ind_s = Grp_Ind_s,
                       Grp_Val_s_up = Grp_Val_s_up, need.exact=T)
  final.fdr <- final.obj$FDP
  # Report Selected
  #fdp(selected,mu);Pow(selected,mu)
  return(list(selected = selected,
              tm0 = tm0,
              ta0 = ta0,
              final.fdr=final.fdr))
}

