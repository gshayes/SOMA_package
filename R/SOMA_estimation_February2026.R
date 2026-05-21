rm(list = ls())

#========================================================================================================================
#FUNCTIONS
# Second Order Meta Analysis
library(metafor)
library(dplyr)

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

random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="REML")

random(ES=TP, SE=SP, cor_mat=PP, iter=10, type="MLE")
#========================================================================================================================
#APPLICATION: 5 Meta-Analyses for Augmented Reality Interventions in Education
#CONSOLIDATED primary ES and SE
# MA1, Chang et al. (2022). Ten years of augmented reality in education: A meta-analysis of (quasi-) experimental studies to investigate the impact. Computers & Education, 191, 104641.
T1 <- c(0.635, 2.971, 1.038, 2.819, -0.273, 0.223, 0.118, 1.819, -0.266, 0.381, 0.704, 0.478, 0.9108, 0.992, 2.012, 1.214, 1.1333, 0.593, 1.056, 0.469, -0.284, 0.868, 0.257, 0.34, 1.831, 0.195, 0.105, -0.276, 0.703, 0.3085, 0.531, 0.2381, -0.467, 1.428, 0.789, 0.286, 0.311, 0.851, 0.697, -0.342, 1.0369, 0.673, 1.083, 0.121, 0.731, -0.02, 0.45, 0.105, 0.185, 0.736, 0.781, 0.157, -0.906, 0.717, 0.079, 0.836, 0.396, 0.707, 0.433, 0.371, 1.755, 0.666, 0.717, 0.573, 0.915, 0.3146, 2.142, 1.184, -0.322, 0.0045, 1.105, 1.326, 3.303, 1.642, 0.808, 0.986, -0.694, 0.565, 0.704, 0.25, -0.332, 0.897, 2.91, 0.793, 0.262, 1.09, -0.188, 0.736, 0.271, 2.4655, 1.829,0.183, 0.525, 0.737, 0.239, 0.183, 0.371, 0.847, -0.346)
S1 <- c(0.2441, 0.289, 0.242, 0.4656, 0.2628, 0.1635, 0.1661, 0.282, 0.3974, 0.1747, 0.2801, 0.4801, 0.1565, 0.1903, 0.363, 0.2304, 0.1469, 0.2352, 0.3319, 0.2696, 0.2944, 0.2079, 0.2745, 0.2508, 0.4273, 0.199, 0.2005, 0.1901, 0.2528, 0.108, 0.2661, 0.2098, 0.2365, 0.2798, 0.261, 0.3179, 0.2505, 0.2663, 0.4423, 0.348, 0.2342, 0.1462, 0.2352, 0.318, 0.2288, 0.4196, 0.3071, 0.377, 0.2117, 0.2444, 0.2985, 0.1643, 0.3684, 0.2115, 0.0844, 0.2579, 0.4145, 0.2018, 0.1997, 0.0964, 0.1936, 0.2431, 0.2995, 0.2413, 0.2878, 0.099, 0.2513, 0.2768, 0.1679, 0.2291, 0.2166, 0.2821, 0.45, 0.2934, 0.1459, 0.2681, 0.2625, 0.2689, 0.2719, 0.2798, 0.2324, 0.2556, 0.439, 0.3028, 0.4747, 0.1416, 0.1875, 0.3676, 0.3594,0.2774, 0.2645, 0.2403, 0.3013, 0.2051, 0.2599, 0.2515, 0.212, 0.2855, 0.3472)

# MA2 , Wang et al. (2024). Impacts of augmented reality-supported STEM education on students' achievement: A meta-analysis of selected SSCI publications from 2010 to 2023. Education and Information Technologies, 29(15), 20547-20585.
T2 <- c(1.038, 1.705, 0.554, 1.819, 0.783, 0.297, 0.469, 0.941, -0.005, 0.531, 0.299, -0.289, 1.083, 0.121, 1.339, 0.416, 0.157, 0.079, 0.836, 0.125, 0.563, 1.766, 0.721, 0.732, -0.0015, 0.409, 1.105, -0.332, 0.281, 0.341, 0.239, 0.183, 0.6152)
S2 <- c(0.242, 0.285, 0.284, 0.282, 0.208, 0.19, 0.2696, 0.109, 0.221, 0.2661, 0.301, 0.312, 0.2352, 0.318, 0.258, 0.307, 0.177, 0.0844, 0.2579, 0.502, 0.261, 0.243, 0.271, 0.321, 0.2369, 0.189, 0.2166, 0.2324, 0.188, 0.231, 0.2599, 0.2515, 0.1685)

# MA3 , Faria & Miranda (2024). Augmented reality in natural sciences and biology teaching: systematic literature review and meta-analysis. Emerging Science Journal, 8(4), 1666-1685.
T3 <- c(1.105, 0.469, 0.195, 0.703, -0.139, -0.131, 1.511, 0.909, 1.755, 1.0686, 5.846, -0.015, 0.128, 0.185, 4.391)
S3 <- c(0.275, 0.2696, 0.199, 0.2528, 0.299, 0.317, 0.321, 0.201, 0.1936, 0.1793, 0.746, 0.233, 0.234, 0.243, 0.503)

# MA4 , Kalemkuş & Kalemkuş (2023). Effect of the use of augmented reality applications on academic achievement of student in science education: meta analysis review. Interactive Learning Environments, 31(9), 6017-6034.
T4 <- c(0.331, 0.61, 0.654, 0.391, -0.284, 0.703, -0.139, -0.147, -0.289, 1.083, 1.755, -0.077, 0.856, 3.681, 1.105, 0.649)
S4 <- c(0.299, 0.301, 0.262, 0.268, 0.295, 0.2528, 0.299, 0.279, 0.312, 0.2352, 0.1936, 0.249, 0.267, 0.449, 0.2166, 0.301)
#========================================================
# Total weights
WT1 <- sum(S1^-2) 
WT2 <- sum((S2^-2)) 
WT3 <- sum(S3^-2) 
WT4 <- sum((S4^-2)) 
#========================================================
#Overlapping weights
WC12 <- sum(S1[c(81, 55, 44, 96, 95, 20, 31, 56, 3, 43, 71, 8)]^-2)
WC13 <- sum(S1[c(26, 29, 61  )]^-2)
WC14 <- sum(S1[c(43,29,71 , 61 )]^-2)

WC23 <- sum(S2[c(7 )]^-2)
WC24 <- sum(S2[c(13,27 )]^-2)

WC34 <- sum(S3[c(5, 4,9 )]^-2)
#=======================================================
#Unique weights
WT1U <- WT1 - WC12 - WC13 - WC14 
WT2U <- WT2 - WC12 - WC23 - WC24
WT3U <- WT3 - WC13 - WC23 - WC34
WT4U <- WT4 - WC14 - WC24 - WC34
#=======================================================
#CORRELATION MATRIX

rho_fun <- function(WTi, WTj, WCij) {
  if (WCij <= 0) return(0)
  WTiU <- WTi - WCij
  WTjU <- WTj - WCij
  1 / sqrt((1 + WTiU/WCij) * (1 + WTjU/WCij))
}


rho12 <- rho_fun(WT1, WT2, WC12)
rho13 <- rho_fun(WT1, WT3, WC13)
rho14 <- rho_fun(WT1, WT4, WC14)

rho23 <- rho_fun(WT2, WT3, WC23)
rho24 <- rho_fun(WT2, WT4, WC24)

rho34 <- rho_fun(WT3, WT4, WC34)


c(rho12=rho12, rho13=rho13, rho14=rho14,
  rho23=rho23, rho24=rho24, 
  rho34=rho34
)


PP <- matrix(0, 4, 4)
diag(PP) <- 1

PP[1,2] <- rho12; PP[2,1] <- rho12
PP[1,3] <- rho13; PP[3,1] <- rho13
PP[1,4] <- rho14; PP[4,1] <- rho14

PP[2,3] <- rho23; PP[3,2] <- rho23
PP[2,4] <- rho24; PP[4,2] <- rho24

PP[3,4] <- rho34; PP[4,3] <- rho34

PP
#===================================================================================
#SYNTHESIS AT THE META-ANALYSIS LEVEL
#WITH FIXED EFFECTS -> As if the MA reported results from a FE MA
#mean effect sizes (TP) and standard errors (SP) for the primary meta-analyses
D1 <- sum(S1^-2*T1)/WT1
D2 <- sum(S2^-2*T2)/WT2
D3 <- sum(S3^-2*T3)/WT3
D4 <- sum(S4^-2*T4)/WT4

SD1 <- sqrt(1/WT1)
SD2 <- sqrt(1/WT2)
SD3 <- sqrt(1/WT3)
SD4 <- sqrt(1/WT4)

TP_FE <-c(D1,D2, D3, D4)
SP_FE <-c(SD1,SD2, SD3, SD4)

#TP_FE and SP_FE vectors can now be used for SOFF and SOFR
#===================================================================================
#WITH RANDOM EFFECTS -> As if the MA reported results from a RE MA

# Helper: within-meta-analysis RANDOM (REML) summary using random()
ma_re_summary <- function(T, S, iter = 50) {
  PP_I <- diag(length(T))  # independence within the meta-analysis
  res  <- random(ES = T, SE = S, cor_mat = PP_I, iter = iter, type = "REML")
  c(mean = as.numeric(res[1, "Mean ES"]),
    se   = as.numeric(res[1, "SE"]),
    tau  = as.numeric(res[1, "tau"]))
}

# Compute RE (REML) summaries for each primary meta-analysis
re1 <- ma_re_summary(T1, S1)
re2 <- ma_re_summary(T2, S2)
re3 <- ma_re_summary(T3, S3)
re4 <- ma_re_summary(T4, S4)

#  vectors to feed into second-order meta-analysis
TP_RE <- c(re1["mean"], re2["mean"], re3["mean"], re4["mean"])
SP_RE <- c(re1["se"],   re2["se"],   re3["se"],   re4["se"])

# (Optional) within-MA taus
tau_within <- c(re1["tau"], re2["tau"], re3["tau"], re4["tau"])

# Print
print("Random-effects (REML) MA-level inputs:")
print(round(cbind(TP_RE = TP_RE, SP_RE = SP_RE, tau_within = tau_within), 6))

#TP_RE and SP_RE vectors can now be used for SORF and SORR

#-----------------------------------------------------
#ALT 1
# Helper: within-meta-analysis RANDOM-effects (REML) summary using rma.uni()
ma_re_summary_rma <- function(T, S) {
  fit <- rma.uni(yi = T, vi = S^2, method = "REML")  # RE model, REML
  
  c(
    mean = as.numeric(fit$b),        # pooled mean
    se   = as.numeric(fit$se),       # SE of pooled mean
    tau  = sqrt(as.numeric(fit$tau2))# tau = sqrt(tau^2)
  )
}
# Compute RE (REML) summaries for each primary meta-analysis
re1_a <- ma_re_summary_rma(T1, S1)
re2_a <- ma_re_summary_rma(T2, S2)
re3_a <- ma_re_summary_rma(T3, S3)
re4_a <- ma_re_summary_rma(T4, S4)

re1_a; re2_a; re3_a; re4_a

#ALT 2
re1_a <- rma(yi=T1, vi=S1^2, method = "REML")
re2_a <- rma(yi=T2, vi=S2^2, method = "REML")
re3_a <- rma(yi=T3, vi=S3^2, method = "REML")
re4_a <- rma(yi=T4, vi=S4^2, method = "REML")

re1_a; re2_a; re3_a; re4_a

#  vectors to feed into second-order meta-analysis
TP_RE_a <- c(
  coef(re1_a),
  coef(re2_a),
  coef(re3_a),
  coef(re4_a)
)

SP_RE_a <- c(
  re1_a$se,
  re2_a$se,
  re3_a$se,
  re4_a$se
)
TP_RE_a <- as.numeric(TP_RE_a)  
SP_RE_a <- as.numeric(SP_RE_a)
TP_RE_a
SP_RE_a

#===================================================================================

#(2) SECOND ORDER META-ANALYSIS
#FE INPUT
#Fixed and random effects SO MA taking dependence/overlapping primary effect sizes into account
fixed(ES=TP_FE, SE=SP_FE, cor_mat=PP)

random(ES=TP_FE, SE=SP_FE, cor_mat=PP, iter=50, type="REML")
#---------------------------------------------------------------------------------------------------
#RE INPUT
#Fixed and random effects SO MA taking dependence/overlapping primary effect sizes into account
fixed(ES=TP_RE, SE=SP_RE, cor_mat=PP)

random(ES=TP_RE, SE=SP_RE, cor_mat=PP, iter=50, type="REML")
#---------------------------------------------------------------
#ALT
random(ES=TP_RE_a, SE=SP_RE_a, cor_mat=PP, iter=50, type="REML")

#===================================================================================
#Second-order MA ignoring dependence/assuming independence
#FE inputs
m <- length(TP_FE)
# Naive (assume independence)
PP_I <- diag(1, m)
fixed(ES=TP_FE, SE=SP_FE, cor_mat=PP_I)
random(ES=TP_FE, SE=SP_FE, cor_mat=PP_I, iter=50, type="REML")
#===========================================================================#
#RE inputs
m <- length(TP_RE)
# Naive (assume independence)
PP_I <- diag(1, m)

fixed(ES=TP_RE, SE=SP_RE, cor_mat=PP_I)
random(ES=TP_RE, SE=SP_RE, cor_mat=PP_I, iter=50, type="REML")

#===========================================================================#
#### SENSITIVITY ANALYSIS ###
#RHO = 0.25

# Fixed between-MA correlation
rho_fixed <- 0.25
m <- 4

PP_rho25 <- matrix(rho_fixed, nrow = m, ncol = m)
diag(PP_rho25) <- 1

colnames(PP_rho25) <- rownames(PP_rho25) <- c("MA1","MA2","MA3","MA4")

PP_rho25
#---------------------------------------------------------------------------------------------------
#FE SOMA, FE inputs
SOF_rho25 <- fixed(
  ES = TP_FE,
  SE = SP_FE,
  cor_mat = PP_rho25
)

SOF_rho25
#---------------------------------------------------------------------------------------------------
#RE SOMA, FE inputs
SOR_rho25 <- random(
  ES = TP_FE,
  SE = SP_FE,
  cor_mat = PP_rho25,
  iter = 50,
  type = "REML"
)

SOR_rho25
#---------------------------------------------------------------------------------------------------
#FE SOMA, RE inputs
SOF_RE_rho25 <- fixed(
  ES = TP_RE,
  SE = SP_RE,
  cor_mat = PP_rho25
)

SOF_RE_rho25

#---------------------------------------------------------------------------------------------------
#RE SOMA, RE inputs
SOR_RE_rho25 <- random(
  ES = TP_RE,
  SE = SP_RE,
  cor_mat = PP_rho25,
  iter = 50,
  type = "REML"
)

SOR_RE_rho25

#---------------------------------------------------------------------------------------------------
#===========================================================================#
#COMPREHENSIVE MA
TT_noMA5 <- c(
  -0.9060, -0.6940, -0.4670, -0.3460, -0.3420, -0.3320, -0.3220, -0.2890,
  -0.2840, -0.2840, -0.2760, -0.2730, -0.2660, -0.1880, -0.1470, -0.1390,
  -0.1310, -0.0770, -0.0200, -0.0150, -0.0050, -0.0015,  0.0045,  0.0790,
  0.1050,  0.1050,  0.1180,  0.1210,  0.1250,  0.1280,  0.1570,  0.1570,
  0.1830,  0.1830,  0.1850,  0.1850,  0.1950,  0.2230,  0.2381,  0.2390,
  0.2500,  0.2570,  0.2620,  0.2710,  0.2810,  0.2860,  0.2970,  0.2990,
  0.3085,  0.3110,  0.3146,  0.3310,  0.3400,  0.3410,  0.3710,  0.3710,
  0.3810,  0.3910,  0.3960,  0.4090,  0.4160,  0.4330,  0.4500,  0.4690,
  0.4780,  0.5310,  0.5540,  0.5630,  0.5650,  0.5730,  0.5930,  0.6100,
  0.6152,  0.6350,  0.6490,  0.6540,  0.6660,  0.6730,  0.6970,  0.7030,
  0.7040,  0.7040,  0.7070,  0.7170,  0.7170,  0.7210,  0.7310,  0.7320,
  0.7360,  0.7360,  0.7370,  0.7810,  0.7830,  0.7890,  0.7930,  0.8080,
  0.8360,  0.8470,  0.8510,  0.8560,  0.8680,  0.8970,  0.9090,  0.9150,
  0.9410,  0.9860,  0.9920,  1.0369,  1.0380,  1.0560,  1.0686,  1.0830,
  1.0900,  1.1050,  1.1050,  1.1333,  1.1840,  1.2140,  1.3260,  1.3390,
  1.4280,  1.5110,  1.6420,  1.7050,  1.7550,  1.7660,  1.8190,  1.8290,
  1.8310,  2.0120,  2.1420,  2.8190,  2.9100,  2.9710,  3.3030,  3.6810,
  4.3910,  5.8460
)

ST_noMA5 <- c(
  0.3684, 0.2625, 0.2365, 0.3472, 0.3480, 0.2324, 0.1679, 0.3120,
  0.2944, 0.2950, 0.1901, 0.2628, 0.3974, 0.1875, 0.2790, 0.2990,
  0.3170, 0.2490, 0.4196, 0.2330, 0.2210, 0.2369, 0.2291, 0.0844,
  0.2005, 0.3770, 0.1661, 0.3180, 0.5020, 0.2340, 0.1643, 0.1770,
  0.2403, 0.2515, 0.2117, 0.2430, 0.1990, 0.1635, 0.2098, 0.2599,
  0.2798, 0.2745, 0.4747, 0.3594, 0.1880, 0.3179, 0.1900, 0.3010,
  0.1080, 0.2505, 0.0990, 0.2990, 0.2508, 0.2310, 0.0964, 0.2120,
  0.1747, 0.2680, 0.4145, 0.1890, 0.3070, 0.1997, 0.3071, 0.2696,
  0.4801, 0.2661, 0.2840, 0.2610, 0.2689, 0.2413, 0.2352, 0.3010,
  0.1685, 0.2441, 0.3010, 0.2620, 0.2431, 0.1462, 0.4423, 0.2528,
  0.2801, 0.2719, 0.2018, 0.2115, 0.2995, 0.2710, 0.2288, 0.3210,
  0.2444, 0.3676, 0.2051, 0.2985, 0.2080, 0.2610, 0.3028, 0.1459,
  0.2579, 0.2855, 0.2663, 0.2670, 0.2079, 0.2556, 0.2010, 0.2878,
  0.1090, 0.2681, 0.1903, 0.2342, 0.2420, 0.3319, 0.1793, 0.2352,
  0.1416, 0.2166, 0.2750, 0.1469, 0.2768, 0.2304, 0.2821, 0.2580,
  0.2798, 0.3210, 0.2934, 0.2850, 0.1936, 0.2430, 0.2820, 0.2645,
  0.4273, 0.3630, 0.2513, 0.4656, 0.4390, 0.2890, 0.4500, 0.4490,
  0.5030, 0.7460
)
TT <- TT_noMA5
ST <- ST_noMA5
#Primary studies that appear in different MAs but with different primary Effect Sizes:
#Chang and Hwang (2018) 
# ES      SE
# 1.1333	0.1469  MA1
#	0.297	  0.19    MA2
#Chien et al. (2019)
# ES      SE
# 0.2381	0.2098  MA1
#-0.139	0.299	    MA3
#-0.139	0.299     MA4
#Giasiranis and Sofos (2017)
# ES    SE
# 0.45	  0.3071  MA1
# 0.416	  0.307   MA2
#Wang (2017)b
# ES        SE
# 0.736 	0.3676  MA1
# 0.281	  0.188   MA2

## Helper: inverse-variance pooled ES/SE for a set of (es,se)
iv_pool <- function(es, se) {
  w <- 1 / (se^2)
  es_bar <- sum(w * es) / sum(w)
  se_bar <- sqrt(1 / sum(w))
  c(es = es_bar, se = se_bar)
}
## ---- Flagged studies (as ES/SE pairs) ----
# 1) Chang & Hwang (2018)
chg_hwg_es <- c(1.1333, 0.297)
chg_hwg_se <- c(0.1469, 0.19)

# 2) Chien et al. (2019)
chien_es <- c(0.2381, -0.139, -0.139)
chien_se <- c(0.2098, 0.299, 0.299)

# 3) Giasiranis & Sofos (2017)
gias_es <- c(0.45, 0.416)
gias_se <- c(0.3071, 0.307)

# 4) Wang (2017b)
wang_es <- c(0.736, 0.281)
wang_se <- c(0.3676, 0.188)
## Bundle for looping
flagged <- list(
  ChangHwang2018 = list(es = chg_hwg_es, se = chg_hwg_se),
  Chien2019      = list(es = chien_es,   se = chien_se),
  GiasSofos2017  = list(es = gias_es,    se = gias_se),
  Wang2017b      = list(es = wang_es,    se = wang_se)
)

## Helper: find indices in TT/ST matching (es,se) pairs
find_pair_indices <- function(TT_noMA5, ST_noMA5, es_vec, se_vec, tol = 1e-10) {
  idx_all <- integer(0)
  for (k in seq_along(es_vec)) {
    idx <- which(abs(TT - es_vec[k]) < tol & abs(ST - se_vec[k]) < tol)
    if (length(idx) == 0) {
      warning(sprintf("No match found for pair (%.10g, %.10g)", es_vec[k], se_vec[k]))
    }
    idx_all <- c(idx_all, idx)
  }
  idx_all
}
## ============================================================
## 1) COLLAPSED SET: replace each flagged set with IV-pooled ES/SE
## ============================================================
TT_avg <- TT_noMA5
ST_avg <- ST_noMA5

for (nm in names(flagged)) {
  es_vec <- flagged[[nm]]$es
  se_vec <- flagged[[nm]]$se
  
  idx <- find_pair_indices(TT_avg, ST_avg, es_vec, se_vec)
  
  # If duplicates exist (same pair appears multiple times), idx may be longer than length(es_vec).
  # That's fine: remove *all* matched rows (most conservative).
  idx <- sort(unique(idx))
  
  pooled <- iv_pool(es_vec, se_vec)
  
  # Remove matched entries
  TT_avg <- TT_avg[-idx]
  ST_avg <- ST_avg[-idx]
  
  # Append pooled entry
  TT_avg <- c(TT_avg, pooled["es"])
  ST_avg <- c(ST_avg, pooled["se"])
}

# Optional: keep sorted
ord_avg <- order(TT_avg)
TT_avg <- TT_avg[ord_avg]
ST_avg <- ST_avg[ord_avg]

#--------------------------------------------------------------------------------------
# 2) Identity correlation matrix (comprehensive unique ES assumed independent)
PP_I_comp <- diag(length(TT_avg))

# 3) Comprehensive FIXED (CF)
#--- AVERAGE EFFECT SIZE PER PRIMARY STUDY ---
COMP_FIXED_avg <- fixed(ES = TT_avg, SE = ST_avg, cor_mat = PP_I_comp)

# 4) Comprehensive RANDOM (REML) (CR)
COMP_RANDOM_REML_avg <- random(ES = TT_avg, SE = ST_avg,
                           cor_mat = PP_I_comp, iter = 50, type = "REML")

# 5) Print results
cat("\nUPDATED COMPREHENSIVE FIXED (no MA5):\n")
print(COMP_FIXED_avg)

cat("\nUPDATED COMPREHENSIVE RANDOM (REML) (no MA5):\n")
print(COMP_RANDOM_REML_avg)

#ALT
rma(yi=TT_avg, vi=ST_avg^2, method = "FE")
rma(yi=TT_avg, vi=ST_avg^2, method = "REML")


#========================================================
#TABLE 4
# Bundle MAs
Ts <- list(MA1 = T1, MA2 = T2, MA3 = T3, MA4 = T4)
Ss <- list(MA1 = S1, MA2 = S2, MA3 = S3, MA4 = S4)
MA_names <- names(Ts)
k <- length(Ts)
# ------------------------------------------------------------------------------
# 1) Build per-MA data frames with keys so overlap is definable via (T,S)
# ------------------------------------------------------------------------------
df_list <- lapply(MA_names, function(ma) {
  data.frame(
    ma  = ma,
    idx = seq_along(Ts[[ma]]),
    T   = Ts[[ma]],
    S   = Ss[[ma]],
    key = paste0(Ts[[ma]], " | ", Ss[[ma]]),
    stringsAsFactors = FALSE
  )
})
pairs_df <- do.call(rbind, df_list)

# Helper: overlap keys between two MAs (exact match of (T,S))
overlap_keys <- function(ma_i, ma_j) {
  intersect(pairs_df$key[pairs_df$ma == ma_i],
            pairs_df$key[pairs_df$ma == ma_j])
}

# ------------------------------------------------------------------------------
# 2) FIXED-effects weights and correlations
# ------------------------------------------------------------------------------
pairs_df$w_FE <- 1 / (pairs_df$S^2)

WT_FE <- tapply(pairs_df$w_FE, pairs_df$ma, sum)

WC_FE <- matrix(0, k, k, dimnames = list(MA_names, MA_names))
for (i in 1:k) {
  for (j in 1:k) {
    if (i == j) next
    ki <- overlap_keys(MA_names[i], MA_names[j])
    if (length(ki) == 0) next
    # sum of weights of overlapping studies as they appear in MA i
    WC_FE[i, j] <- sum(pairs_df$w_FE[pairs_df$ma == MA_names[i] &
                                       pairs_df$key %in% ki])
  }
}

RHO_FE <- matrix(0, k, k, dimnames = list(MA_names, MA_names))
diag(RHO_FE) <- 1
for (i in 1:k) {
  for (j in 1:k) {
    if (i == j) next
    RHO_FE[i, j] <- WC_FE[i, j] / sqrt(WT_FE[i] * WT_FE[j])
  }
}

# ------------------------------------------------------------------------------
# 3) RANDOM-effects weights (REML within each MA) and correlations
# ------------------------------------------------------------------------------
tau2_hat <- sapply(MA_names, function(ma) {
  fit <- rma.uni(yi = Ts[[ma]], sei = Ss[[ma]], method = "REML")
  as.numeric(fit$tau2)
})

pairs_df$tau2 <- tau2_hat[pairs_df$ma]
pairs_df$w_RE <- 1 / (pairs_df$S^2 + pairs_df$tau2)

WT_RE <- tapply(pairs_df$w_RE, pairs_df$ma, sum)

WC_RE <- matrix(0, k, k, dimnames = list(MA_names, MA_names))
for (i in 1:k) {
  for (j in 1:k) {
    if (i == j) next
    ki <- overlap_keys(MA_names[i], MA_names[j])
    if (length(ki) == 0) next
    WC_RE[i, j] <- sum(pairs_df$w_RE[pairs_df$ma == MA_names[i] &
                                       pairs_df$key %in% ki])
  }
}

RHO_RE <- matrix(0, k, k, dimnames = list(MA_names, MA_names))
diag(RHO_RE) <- 1
for (i in 1:k) {
  for (j in 1:k) {
    if (i == j) next
    RHO_RE[i, j] <- WC_RE[i, j] / sqrt(WT_RE[i] * WT_RE[j])
  }
}

# ------------------------------------------------------------------------------
# 4) Combine into Table 4 layout:
#    - lower triangle = fixed
#    - upper triangle = random
# ------------------------------------------------------------------------------
Table4 <- matrix("", k, k, dimnames = list(MA_names, MA_names))

for (i in 1:k) {
  for (j in 1:k) {
    if (i == j) {
      Table4[i, j] <- "1.000"
    } else if (i > j) {
      # below diagonal: fixed
      Table4[i, j] <- sprintf("%.3f", RHO_FE[i, j])
    } else {
      # above diagonal: random
      Table4[i, j] <- sprintf("%.3f", RHO_RE[i, j])
    }
  }
}

cat("\nTable 4 (lower=Fixed, upper=Random REML):\n")
print(Table4, quote = FALSE)





#================================================================================





















