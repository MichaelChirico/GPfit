## Calculates the deviance values 
##
## May 8th, 2012

GP_deviance <- 
  function(beta, X, Y, nug_thres=20,
           corr = list(type="exponential", power=1.95)){
  if (is.matrix(X) == FALSE){
    X = as.matrix(X)
  }
  n = nrow(X);
  d = ncol(X);
  if (n != length(Y)){
    stop("The dimensions of X and Y do not match. \n")
  }
  ## Checking the dimensions of the parameters
  if (d != length(beta)){
    stop("The dimensions of beta and X do not match \n")
  }
  
  if (nug_thres < 10 | nug_thres > 25){
    warning("nug_thres outside of the normal range of [10, 25]")
  }
  ex_nug = exp(nug_thres)
  
  One = rep(1L, n);
  dim(Y) = c(n, 1L);
  
  dim(beta) = c(1L, d);
  
  R = corr_matrix(X, beta, corr);
  condnum = kappa(R, triangular = TRUE, exact=TRUE);
  
  if (ex_nug > condnum) {
    delta = 0 
  } else {
    delta = 
      abs(eigen(R, symmetric = TRUE, only.values = TRUE)$values[1L]) *
      (condnum-ex_nug)/condnum/(ex_nug-1)
  }
  
  #inverse of variance + delta*I_n
  Sig = R + diag(rep(delta, n))
  Sig_inv  = solve(Sig)
  #difference between observed Y and mu_hat
  mu_hat_diff = Y - sum(Sig_inv) / sum(Sig_inv %*% Y)

  part1 = 2*sum(log(abs(eigen(chol(Sig), only.values = TRUE)$values)))
  part2 = n * log(drop(t(mu_hat_diff) %*% Sig_inv %*% mu_hat_diff))
  devval = part1 + part2
  
  if (!is.finite(devval))
  {
    stop('Infinite values of the Deviance Function, 
		unable to find optimum parameters \n')
  }
  return(devval)
}
