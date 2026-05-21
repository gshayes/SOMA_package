library(metafor)
library(dplyr)
# load(functions.R)

# EXAMPLE 1 ===================================================================
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

# EXAMPLE 2 ===================================================================

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