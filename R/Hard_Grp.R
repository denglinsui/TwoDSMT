#' Hard Group Features
#'
#' Divide each features into several groups with knn. This allows us obtain over-estimated and under-estimated FDP.
#' @param Org_Value a n*p matrix or a dataframe representing observed values of n items with p features; each column represents one feature.
#' @param n.group.max an integer representing the maximum number of groups for each feature.
#'
#' @return a list;
#' Grp_Value is a list where value represents the minimum value of each group, value_up represents the maximum value of each group, and ind represents the item indexes of each group;
#' Grp_Info is a n*p matrix where i-j element of Grp_Info represents which group the j-th feature of i-th item belongs to;
#' Grp_Val is a n.group.max*p matrix i-j element of where Grp_Val represents the minimum value of i-th group of j-th feature;
#' Grp_Val_up is a n.group.max*p matrix where i-j element of Grp_Val_up represents the maximum value of i-th group of j-th feature;
#' Grp_Val_summarize is the unique values of Grp_Val (in general, there will be n.group.max^p rows);
#' Grp_Val_up_summarize is the unique values of Grp_Val_up (in general, there will be n.group.max^p rows);
#' Grp_Ind_summarize is list where i-th term is an integer vector indicating the indexes of the items belongs to i-th row of Grp_Val_summarize/Grp_Val_up_summarize
#' is.exact.group is a logical value representing whether the group is exact (for example, if the value is discrete, the group would be exact)
#' @importFrom stats quantile kmeans uniroot
#' @export
#'
#' @examples
#' Org_Value <- matrix(runif(100),50,2)
#' res <- Hard_Group(Org_Value,n.group.max=3)
Hard_Group <- function(Org_Value,
                       n.group.max = 10){
  is.exact.group = T
  Grp_Value <- list()
  if(is.matrix(Org_Value)){
    Org_Value <- as.data.frame(Org_Value)
  }
  n.item <- dim(Org_Value)[1]
  n.val <- dim(Org_Value)[2]
  n.group <- rep(0,n.val)

  Grp_Info <- matrix(data=NA,nrow = n.item, ncol=n.val)
  Grp_Val_item <- matrix(data=NA,nrow = n.item, ncol=n.val)
  Grp_Val_up_item <- matrix(data=NA,nrow = n.item, ncol=n.val)
  Grp_Val <- matrix(data = NA, nrow = n.group.max, ncol = n.val)
  Grp_Val_up <- matrix(data = NA, nrow = n.group.max, ncol = n.val)
  if(length(n.group.max)==1){
    n.group.max <- rep(n.group.max,n.val)
  }

  for(j in 1:n.val){
    grp_value <- list() #To store the temporate output
    org_value <- Org_Value[,j]
    # Check whether the sequence is grouped
    order.index <- order(org_value)
    sort.value <- org_value[order.index]
    jump.ind <- which(diff(sort.value)>0)

    # The jump.ind indicates the exct group of sequence
    #If the groups of the sequence are too many, we need to assign them a group
    if(length(jump.ind)>=n.group.max[j]){
      # Basic1: Split to let the difference be large
      #jump.ind <- order(diff(sort.value),decreasing=TRUE)[1:(n.group.max[j]-1)]
      #
      # Basic2: Split to make them distributed unifomly
      #jump.ind <- order(diff(sort.value),
      #                  decreasing=TRUE)[floor(seq(1,n.item,length.out =n.group.max[j]))[1:(n.group.max[j]-1)]]
      # Basic3: Use kmeans to split
      is.exact.group = F # We use approximation

      kmeans.res <- NULL
      n.group.max[j] <- n.group.max[j]+1
      while(is.null(kmeans.res) & n.group.max[j]>0){
        n.group.max[j] <- n.group.max[j]-1
        kmeans.center <- unique(quantile(org_value,seq(0,1,length.out=n.group.max[j])))
        #print(kmeans.center)
        kmeans.res <- tryCatch(kmeans(org_value, centers = kmeans.center),
                               error = function(e) {NULL})
      }
      len.group <- length(kmeans.center)
      #kmeans.center <- unique(quantile(org_value,seq(0,1,length.out=n.group.max[j])))
      #len.group <- length(kmeans.center)
      #kmeans.res <- kmeans(org_value,
      #                     center = kmeans.center)
      #center = seq(sort.value[1],sort.value[n.item],length.out=n.group.max[j]))

      grp_ind <- list()
      grp_value$value <- numeric(len.group)
      for(k in 1:len.group){
        grp_ind[[k]] <- which(kmeans.res$cluster==k)
        grp_value$value[k] <- min(org_value[grp_ind[[k]]])
        grp_value$value_up[k] <- max(org_value[grp_ind[[k]]])
      }
      grp_value$ind <- grp_ind

      n.group[j] <- len.group
      Grp_Info[,j] <- kmeans.res$cluster
      Grp_Val[1:len.group,j] <- grp_value$value
      Grp_Val_up[1:len.group,j] <- grp_value$value_up
      Grp_Val_item[,j] <- grp_value$value[kmeans.res$cluster]
      Grp_Val_up_item[,j] <- grp_value$value_up[kmeans.res$cluster]
      # We can also let the index to be spreaded out as more as possible
      # There is a trade-off, but hasn't been detailedly considered yet.
    }else{
      # The case the weight can be grouped
      grp_value$value <- sort.value[c(1,jump.ind+1)]
      jump.ind.extend <- c(0,jump.ind,n.item)
      grp_ind <- list()
      for(k in 1:(length(jump.ind.extend)-1)){
        grp_ind[[k]] <- order.index[(jump.ind.extend[k]+1):jump.ind.extend[k+1]]
        Grp_Info[grp_ind[[k]],j] <- k
      }

      n.group[j] <- length(jump.ind)+1
      grp_value$ind <- grp_ind
      Grp_Val[1:(length(jump.ind)+1),j] <- grp_value$value
      Grp_Val_up[1:(length(jump.ind)+1),j] <- grp_value$value
      Grp_Val_item[,j] <- grp_value$value[Grp_Info[,j]]
      Grp_Val_up_item[,j] <- grp_value$value[Grp_Info[,j]]
    }

    Grp_Value[[j]] <- grp_value
  }

  #Grp_Value$Grp_Info <- Grp_Info
  #Grp_Value$Grp_Val <- Grp_Val

  #Generate the total group
  Grp_Val_summarize <- unique(Grp_Val_item)
  Grp_Val_up_summarize <- unique(Grp_Val_up_item)
  Grp_Ind_summarize <- list()
  for(i in 1:nrow(Grp_Val_summarize)){
    grp_val_single <- Grp_Val_summarize[i,]
    Grp_Ind_summarize[[i]] <- which(apply(Grp_Val_item,1,function(x){all(x==grp_val_single)}))
  }

  return(list(Grp_Value=Grp_Value,
              Grp_Info = Grp_Info,
              Grp_Val = Grp_Val,
              Grp_Val_up = Grp_Val_up,
#              Grp_Val_item = Grp_Val_item,
              Grp_Val_summarize = Grp_Val_summarize,
              Grp_Val_up_summarize = Grp_Val_up_summarize,
              Grp_Ind_summarize = Grp_Ind_summarize,
              is.exact.group = is.exact.group
  ))
}
