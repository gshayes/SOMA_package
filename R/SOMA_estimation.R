# Second Order Meta Analysis
library(metafor)

# Fixed Effects SO MA of overlapping studies
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
  
  ANS[1] <-mu_hat
  ANS[2] <- SE_mu_hat
  ANS[3] <- Q
  ANS[4] <- df
  ANS[5] <- pvalue
  return(ANS)
}

# Random Effects SO MA of overlapping studies
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

# This creates an example of 2 independent meta-analyses
T1 <- c(.2,.5)
T2 <- c(.3,.9)
TT <- c(T1,T2)
S1 <- c(.1,.2)
S2 <- c(.12,.15)
ST <- c(S1,S2)

# The common weights, WC, is 0 since the MAs are independent
WC <- 0
WT1 <- sum(S1^-2) 
WT2 <- sum((S2^-2)) 
rho <- 1/sqrt((1 + WT1/WC)*(1 + WT2/WC))

# This generates the correlation matrix between meta-analyses
PP <- matrix(nrow=2,ncol=2,data=rho)
for (i in 1:2){PP[i,i]<-1}

# These are the mean effect sizes (TP) and standard errors (SP) for the two
# primary meta-analyses
D1 <- sum(S1^-2*T1)/WT1
SD1 <- sqrt(1/WT1)
D2 <- sum(S2^-2*T2)/WT2
SD2 <- sqrt(1/WT2)
TP <-c(D1,D2)
SP <-c(SD1,SD2)

# This generates the comprehensive (COMP) meta-analysis 
DT <- sum(ST^-2*TT)/sum(ST^-2)
SE <- sqrt(1/sum(ST^-2))
COMP <- c(DT,SE)

# compare results from fixed function above and comprehensive MA
fixed(ES=TP, SE=SP, cor_mat=PP)
print(COMP,digits=8)

# compare results from random function above to metafor
rma(yi=TP, vi=SP^2, method = "REML")
random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="REML")

rma(yi=TP, vi=SP^2, method = "ML")
random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="MLE")

# This creates the example of 2 overlapping meta-analyses
# T1 and T2 with TT are the 2 primary and the comprehensive effect sizes
# and S1, S2, and ST are the 2 primary and the comprehensive SEs
T1 <- c(.2,.5)
T2 <- c(T1[2],.9)
TT <- c(T1,T2[2])
S1 <- c(.1,.2)
S2 <- c(S1[2],.15)
ST <- c(S1,S2[2])

# WC is the shared weight. WT1 and WT2 are the total weights, 
# WT1U and WT2U are the unique (non-overlapping) weights,
# and rho is the correlation between the two meta-analyses
WC <- S1[2]^-2 
WT1 <- sum(S1^-2) 
WT2 <- sum((S2^-2)) 
WT1U <- WT1 - WC
WT2U <- WT2 - WC
rho <- 1/sqrt((1 + WT1U/WC)*(1 + WT2U/WC))

# This generates the correlation matrix between meta-analyses
PP <- matrix(nrow=2,ncol=2,data=rho)
for (i in 1:2){PP[i,i]<-1}

# These are the mean effect sizes (TP) and standard errors (SP) for the two
# primary meta-analyses
D1 <- sum(S1^-2*T1)/WT1
SD1 <- sqrt(1/WT1)
D2 <- sum(S2^-2*T2)/WT2
SD2 <- sqrt(1/WT2)
TP <-c(D1,D2)
SP <-c(SD1,SD2)

# This generates the two primary (PRIM) and comprehensive (COMP) meta-analyses 
DT <- sum(ST^-2*TT)/sum(ST^-2)
SE <- sqrt(1/sum(ST^-2))
PRIM <- c(TP,SP)
COMP <- c(DT,SE)

fixed(ES=TP, SE=SP, cor_mat=PP)
print(COMP)

random(ES=TP, SE=SP, cor_mat=PP, iter=20, type="REML")

random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="MLE")

# This creates the example of 2 overlapping meta-analyses
# T1 and T2 with TT are the 2 primary and the comprehensive effect sizes
# and S1, S2, and ST are the 2 primary and the comprehensive SEs
T1 <- c(.2, .5, .4)
T2 <- c(T1[2], .9, .6)
TT <- c(T1, T2[2], T2[3])
S1 <- c(.1, .2, .12)
S2 <- c(S1[2], .15, 0.9)
ST <- c(S1, S2[2], S2[3])

# WC is the shared weight. WT1 and WT2 are the total weights, 
# WT1U and WT2U are the unique (non-overlapping) weights,
# and rho is the correlation between the two meta-analyses
WC <- S1[2]^-2 
WT1 <- sum(S1^-2) 
WT2 <- sum((S2^-2)) 
WT1U <- WT1 - WC
WT2U <- WT2 - WC
rho <- 1/sqrt((1 + WT1U/WC)*(1 + WT2U/WC))

# This generates the correlation matrix between meta-analyses
PP <- matrix(nrow=2,ncol=2,data=rho)
for (i in 1:2){PP[i,i]<-1}

# These are the mean effect sizes (TP) and standard errors (SP) for the two
# primary meta-analyses
D1 <- sum(S1^-2*T1)/WT1
SD1 <- sqrt(1/WT1)
D2 <- sum(S2^-2*T2)/WT2
SD2 <- sqrt(1/WT2)
TP <-c(D1,D2)
SP <-c(SD1,SD2)

# This generates the two primary (PRIM) and comprehensive (COMP) meta-analyses 
DT <- sum(ST^-2*TT)/sum(ST^-2)
SE <- sqrt(1/sum(ST^-2))
PRIM <- c(TP,SP)
COMP <- c(DT,SE)

fixed(ES=TP, SE=SP, cor_mat=PP)
print(COMP)

random(ES=TP, SE=SP, cor_mat=PP, iter=20, type="REML")

random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="MLE")
