## Optimizes the deviance using multiple optim starts
## 
## May 8th, 2012

GP_fit <- function(X,Y,control=c(200*d,80*d,2*d),nug_thres=20,
                   trace=FALSE,maxit = 100,
                   corr = list(type="exponential",power=1.95), 
                   optim_start = NULL){
  if (is.matrix(X) == FALSE){
    X = as.matrix(X)
  }
  
  ## Checking to make sure that the X design is scaled between 0 and 1
  if ((min(X) <0) | (max(X)>1) | ((max(X) - min(X)) <=0.5)){
    cat("The optimization routine assumes that inputs are in [0,1].\n")
  }
  
  n = nrow(X)
  d = ncol(X)
  ## Checking the dimensions between the different inputs
  if (n != length(Y)){
    stop("The dimensions of X and Y do not match. \n")
  }
  if (n != length(Y)){
    stop("The dimension of the training points and simulator values are different. \n")
  }
  if (nug_thres < 10 | nug_thres > 25) {
    warning("nug_thres is outside of the normal range of [10, 25].")
  }
  ## Need to make sure that the control values are relatively higher / lower and need 
  ## to print a warning if they are not
  if (length(control) != 3){
    stop("control is defined incorrectly. Wrong number of arguments. \n")
  }
  if (control[1] < control[2]){ 
    stop("control is defined incorrectly. Need control[1] >= control[2]. \n")
  }
  if (control[1] < control[3]){
    stop("control is defined incorrectly. Need control[1] >= control[3]. \n")
  }
  if (control[2] < control[3]){
    stop("control is defined incorrectly. Need control[2] >= control[3]. \n")
  }
  
  ## Checking to see if the vigen starting values for the optim runs are in the correct 
  ## format
  if (!is.null(optim_start)){
    if(!is.matrix(optim_start)){
      if(length(optim_start)/d != floor(length(optim_start)/d)){
        stop("The dimension of optim_start does not match the dimension
					of the problem \n")
      }
      optim_start = matrix(optim_start, byrow = TRUE, ncol = d)
    } else if(ncol(optim_start)!= d){
      stop("The dimension of optim_start does not match the dimension
				of the problem \n")
    }
  }
  
  param_search = control[1L];
  param_percent = control[2L];
  param_clust = control[3L];
  ###########################################################
  ## Using a grid of control[1] points as default, need to find the 
  ## control[2] corresponding values with the lowest deviance
  ## starting out based on the defined range for beta
  beta_range = if (corr$type == "exponential") {
    corr$power - log10(d) + c(-4, log10(5))
  } else if (corr$type == "matern") {
    2 - log10(d) + c(0, log10(5))
  }
  param_init_200d = 
    lhs::maximinLHS(param_search, d)*diff(beta_range)+beta_range[1L];
  param_lower = rep(-10, d);
  param_upper = -param_lower

  ##-------------------------------------------------------##
  ## Need to evaluate the deviance for the control[1] points
  deviance1 = 
    do.call(rbind, lapply(seq_len(param_search), function(i) {
      parami = param_init_200d[i, ]
      c(GP_deviance(parami, X, Y, nug_thres, corr = corr), parami)
    }))
  ## Need to order the initial values based on their deviance
  deviance2 = deviance1[order(deviance1[ , 1L]), ];
  
  ##-------------------------------------------------------##
  ## Taking the control[2] smallest deviance values
  deviance3 = deviance2[seq_len(param_percent), ];
  
  ##-------------------------------------------------------##
  ## Going to cluster these control[2] "best" observations into control[3]
  ## groups using K-means, over 5 iterations and the smallest WSS as 
  ## the criterion
  points_percen = deviance3[ , 2L:(d+1L)];
  ## Taking the best of 5 different runs of K-means
  k = lapply(integer(5L), function(...) kmeans(points_percen, param_clust))
  kmax = k[[which.max(sapply(k, function(kk) sum(kk$withinss)))]]

  ##---------------------------------------------------------##
  ## Need to set the control[3] cluster centers as starting values
  
  #param_init = k[[1]]$centers;
  
  
  ##-------------------------------------------------------------------##
  ## Need to set the control[3] best points in cluster as starting values
  
  param_init = matrix(nrow = param_clust, ncol = d)
  for(i in seq_len(param_clust)){
    ID = which(kmax$cluster == i)
    fID = which.min(deviance3[ID, 1L])
    param_init[i, ] = deviance3[ID[fID], 2L:(d+1L)]
  }
  
  
  
  
  
  
  #############################################################
  ## For the diagonal search: use 3 points along the diagonal, 
  ## 1 near each of the ends and one near the middle of the range, 
  ## this is only necessary above 1 dimension
  if (d >= 2) {
    param_wrap = c(0.20, .5, .8)*diff(beta_range) + beta_range[1L];
    
    ## Need to run optim() on the wrapped function the 3 times in order to find
    ## the lowest starting ponit of the 3 values
    dev = do.call(rbind, lapply(param_wrap, function(pp) {
      with(optim(pp, dev_wrapper, X=X, Y=Y, nug_thres=nug_thres,
                 corr=corr, method="L-BFGS-B",
                 lower=param_lower, upper=param_upper, 
                 control = c(maxit = maxit)),
           c(par, value, pp))
    }));
    ## Take the best of the 3 based on the likelohood value
    dev = dev[order(dev[ , 2L]), ];
    
    ##---------------------------------------------------------##
    ## Combining the 2*d centers from the clusters with the single point from 
    ## the diagonal search
    param_init = rbind(rep(dev[1L], d), param_init);
  }
  
  #############################################################
  ## We now have the control[3]+1 initializing points for the optim search, as well
  ## as any points that are specified by the user
  param_init = rbind(param_init, optim_start)
  dev_val = do.call(rbind, lapply(seq_len(nrow(param_init)), function(ii) {
    with(optim(param_init[i, ], GP_deviance, X=X, Y=Y,
               nug_thres=nug_thres, corr=corr,
               method="L-BFGS-B", lower=param_lower, upper=param_upper, 
               control = c(maxit = maxit)),
         c(par, value))
  }))
  #############################################################
  ## Making a print statement of what the progress of the optimizer
  ## only if trace == TRUE
  if (trace == TRUE){
    optim_result = cbind(param_init,dev_val)
    col_name = NULL
    if(d==1){
      row_name = NULL
    } else {
      row_name = c("Diagonal Search")
    }
    for (i in 1:d){
      col_name = cbind(col_name, paste("Beta", as.character(i), "Start"))
    }
    for (i in 1:d){
      col_name = cbind(col_name, paste("Beta", as.character(i), "Final"))
    }
    col_name = cbind(col_name, "Deviance Value")
    for (i in 1:(param_clust)){
      row_name = cbind(row_name, paste("Start", as.character(i)))
    }
    colnames(optim_result) = col_name;
    rownames(optim_result) = row_name;
    print(optim_result)
  }
  dev_val = dev_val[order(dev_val[,(d+1)]),];
  beta = (dev_val[1,1:d]);
  
  dim(beta) = c(d,1);
  R = corr_matrix(X,beta,corr);
  temp = eigen(R,symmetric = TRUE, only.values = TRUE);
  eig_val = temp$values;
  condnum = kappa(R,triangular = TRUE,exact=TRUE);
  max_eigval = eig_val[1];
  delta = max(c(0,abs(max_eigval)*(condnum-exp(nug_thres))/(condnum*(exp(nug_thres)-1))));
  
  One = rep(1,n);
  dim(Y) = c(n,1);
  
  LO = diag(n);
  Sig = R + delta*LO;
  
  L = chol(Sig);
  
  Sig_invOne = solve(L,solve(t(L),One));	
  Sig_invY = solve(L,solve(t(L),Y));	
  
  mu_hat = solve(t(One)%*%Sig_invOne,t(One)%*%Sig_invY);
  Sig_invb = solve(L,solve(t(L),(Y-One%*%mu_hat)));
  sig2 = t(Y-One%*%mu_hat)%*%Sig_invb/n;
  
  GP = list(X = X, Y = Y, sig2 = as.vector(sig2),
            beta = beta, delta = delta, 
            nugget_threshold_parameter = nug_thres,
            correlation_param = corr)
  
  class(GP) = "GP"
  return(GP)
}