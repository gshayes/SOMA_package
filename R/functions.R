#========================================================================================================================
#FUNCTIONS

#' Fixed Effects SO MA of overlapping studies
#'
#' @param ES Effect sizes from the meta analyses
#' @param SE Standard errors corresponding to the effect sizes
#' @param cor_mat Covariance matrix of the effect sizes
#'
#' @returns SOMA estimated effect size, standard error, Q statistic, degrees of freedom, and p-value
#' @export
fixed <- function(ES, SE, cor_mat) {
  
  # number of MAs
  k = length(ES)
  
  # vector of 1s
  O = matrix(nrow = k, ncol = 1, data = 1)
  O_t = t(O)
  
  # initialize Sigma matrix
  Sig = matrix(nrow = k, ncol = k)
  
  # fill in Sigma values with variances and covariances
  for (i in 1:k) {
    for (j in 1:k) {
      if(i==j) { Sig[i,j] = SE[i]^2} else { 
        Sig[i,j] = cor_mat[i,j] * SE[i] * SE[j] }
    }
  }
  
  # solve for inverse Sigma
  Sig_inv = solve(Sig)
  
  # compute mean effect size
  mu_hat = (O_t %*% Sig_inv %*% ES) / (O_t %*% Sig_inv %*% O)
  
  # compute variance of the estimator
  V_mu_hat = 1 / (O_t %*% Sig_inv %*% O)
  
  # compute standard error of the estimator
  SE_mu_hat = sqrt(V_mu_hat)
  
  # compute Q
  Q = t(ES) %*% Sig_inv %*% ES - (O_t %*% Sig_inv %*% ES)^2 / (O_t %*% Sig_inv %*% O)
  
  # compute degrees of freedom for chi-square distribution
  df = k-1
  
  # compute p-value for heterogeneity test
  pvalue = 1 - pchisq(q=Q, df=df)
  
  # print out results
  ANS <- matrix(nrow=1,ncol=5)
  colnames(ANS) <- c("Mean ES", "SE","Q","df", "p" )
  
  ANS[1] <- mu_hat
  ANS[2] <- SE_mu_hat
  ANS[3] <- Q
  ANS[4] <- df
  ANS[5] <- pvalue
  return(ANS)
}

#' Random Effects SO MA of overlapping studies
#'
#' @param ES Effect sizes from the meta analyses
#' @param SE Standard errors corresponding to the effect sizes
#' @param cor_mat Covariance matrix of the effect sizes
#' @param iter Number of iterations
#' @param type ML or REML
#'
#' @returns SOMA estimated effect size, standard error, Q statistic, degrees of freedom, p-value, and estimated tau-squared.
#' @export
random <- function(ES, SE, cor_mat, iter, type) {
  
  # initial fixed-effect results from TDotFM
  fix_results  = fixed(ES, SE, cor_mat)
  mu  = fix_results[1]
  Q   = fix_results[3]
  df  = fix_results[4]
  p   = fix_results[5]
  
  # number of meta-analyses
  k = length(ES)
  
  # vector of 1s
  O = matrix(nrow = k, ncol = 1, data = 1)
  O_t = t(O)
  
  # starting value for tau^2
  tau2 = 0
  
  # create covariance matrix
  Sig = matrix(0, nrow = k, ncol = k)
  
  # fill in Sigma values with variances and covariances
  for (i in 1:k) {
    for (j in 1:k) {
      if(i==j) { Sig[i,j] = SE[i]^2 + tau2} else { 
        Sig[i,j] = cor_mat[i,j] * SE[i] * SE[j] }
    }
  }
  
  # store iteration estimates
  iter_store = matrix(nrow = iter, ncol = 3)
  
  # begin iteration
  for (t in 1:iter) {
    
    Sig_inv = solve(Sig)
    V_mu = 1 / (O_t %*% Sig_inv %*% O)
    W = O_t %*% Sig_inv
    
    # REML correction term
    REML_correction = 0
    if (type == "REML") {
      REML_correction = 1 / (W %*% O)
    }
    
    # update tau^2 estimator
    tau2 = sum(W^2 * ((ES - mu)^2 - SE^2)) / sum(W^2) + REML_correction
    
    # update SS and Sig
    if (tau2 > 0) {
      for (i in 1:k) {
        for (j in 1:k) {
          if(i==j) { Sig[i,j] = SE[i]^2 + tau2} else { 
            Sig[i,j] = cor_mat[i,j] * SE[i] * SE[j] }
        }
      }
    }
    
    # updated diagonal SEs for fixed call
    new_SE = numeric(k)
    for (i in 1:k) {
      new_SE[i] = sqrt(Sig[i, i])
    }
    
    # recompute effects
    TD_update = fixed(ES, new_SE, cor_mat)
    mu  = TD_update[1]
    SE_mu = TD_update[2]
    
    # store progress
    tau = ifelse(tau2 >= 0, sqrt(tau2), 0)
    iter_store[t, 1] = mu
    iter_store[t, 2] = SE_mu
    iter_store[t, 3] = tau
  }
  # output results
  ANS = matrix(nrow = 1, ncol = 6)
  colnames(ANS) = c("Mean ES", "SE", "Q", "df", "p", "tau")
  
  ANS[1] = mu
  ANS[2] = SE_mu
  ANS[3] = Q
  ANS[4] = df
  ANS[5] = p
  ANS[6] = tau
  return(ANS)
}