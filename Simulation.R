#----------------------------------------------------------------------------------------------------------------------------
#   Simulation program for multilevel ordinal outcome with missing data
#
# This program runs the simultation study presented in 
# "Multiple Imputation Methods for Missing Multilevel Ordinal Outcomes"
# Authors: Mei Dong, Aya Mitani
#
# Code Written by: Mei Dong
# need to specifiy a few parameters: N (sample size), truenu (scale for ICS), truetau (scale for ICC), rate (missing rate), pattern (missing pattern: MAR, MNAR, MCAR), 
# The R version we used is 4.1.2

library(geepack)
library(bridgedist)
library(MASS)
library(Matrix)
library(tidyverse)
library(reshape2)
library(gee)
library(multgee)
library(jomo) ## version of jomo package is 2.7.2
library(mitools)

#-----------------------------
# simulation parameters
#-----------------------------

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
print(args)


# ## args is now a list of character vectors
# ## First check to see if arguments are passed.
# ## Then cycle through each element of the list and evaluate the expressions.
#
if(length(args)==0){
  print("No arguments supplied.")
  ## supply default values
  N=50
  #ICS
  truenu=0.4
  #teeth correlation level
  truetau=0.3
  #missing rate
  rate=0.2
  #missing pattern
  pattern="MAR"
  nteeth=28
  
  rep=1
}else{
  for (i in (1:length(args))) eval(parse(text=args[[i]]))
}

SS <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(SS)) SS <- 1

print(paste("N=",N))
print(paste("SS=",SS))
print(paste("rep=",rep))

set.seed(rep)

## load function for generating missing data
source("generate_missing.R")



##define parameters
#betas
truebetac_y1=c(-0.4, 0.8, 1.6)
truebetac_y2=c(-1, 0.5, 1)
truebetac_y3=c(1, 2, 3)
truebetac_y4=c(-2, -0.5, 1)

truebeta1=-0.2
truebeta2=-0.5

### define more parameters
truebetac_y_mat <- rbind(truebetac_y1,truebetac_y2, truebetac_y3,truebetac_y4)
numcat <- length(truebetac_y1) + 1
catvec <- c(1:numcat)

truebeta <- c(truebeta1, truebeta2)


nteeth1 <- nteeth - 1
tau <- sin( pi * truetau / 2)
p <- length(truebeta)
phi <- 0.5

xi <- rnorm(N, 0, 1)
zi <- rbinom(N,1, 0.3)

##sample w from multivariate normal with mean 0, covariance thesigma
themu <- rep(0, nteeth)

##correlation matrix for teeth
thesigma <- toeplitz(c(1, rep(tau, nteeth1)))

maxnstayiter <- 1000

nstay <- 0

ally <- vector("list", N)
nteethpp <- baserisk <- rep(NA, N)
teethpp <- vector("list", N)

#######################################     Begin simulation     #########################################

## remove data that can't converge for full data
est_wgee_full_int1 <- est_wgee_full_int5 <- 100
se_wgee_full_int1 <- se_wgee_full_int5 <- 100
count <- 0
while(abs(est_wgee_full_int1)>10 | se_wgee_full_int1==0 | se_wgee_full_int1>5 | abs(est_wgee_full_int5)>10 | se_wgee_full_int5==0 | se_wgee_full_int5>5){
  count <- count+1
  y1_table <- c(1,1,1,1)
  while(any(y1_table<10)){
    for(i in 1:N){
      
      nstayiter <- maxnstayiter
      
      while(nstayiter == maxnstayiter){
        #step1
        Z <- mvrnorm(n = 1, mu = themu, Sigma = thesigma)
        #step2
        u <- pnorm(Z)
        #step3
        b <- ( 1 / phi ) * log( sin( phi * pi * u) / sin( phi * pi * ( 1 - u ) ) )
        #step4
        lambda <- exp(truenu * mean(b)) / ( 1 + exp(truenu * mean(b)) )
        baserisk[i] <- lambda
        
        #step5, cluster size
        nstayiter <- 0
        
        while(nstay < 2 & nstayiter < maxnstayiter){
          nstay <- rbinom(n = 1, size = nteeth, prob = lambda)
          #nstay <- rbinom(n = 1, size = nteeth, prob = 0.25)
          nstayiter <- nstayiter + 1
        }
      }
      
      nteethpp[i] <- nstay
      ##sample which teeth to stay
      tstay <- sort(sample(c(1:nteeth), nstay))
      ##teeth per participant
      teethpp[[i]] <- tstay
      
      bstay <- b[tstay]
      
      yvec <- rep(NA, nstay)
      ymat <- matrix(NA, nrow=nstay, ncol = 4)
      
      #step6
      truexb <- truebeta[1] * xi[i] + truebeta[2] * zi[i] 
      
      for (m in 1:4){
        for (t in 1:nstay ){
          
          est1 <- exp( bstay[t] + ( truebetac_y_mat[m, 1] + truexb ) / phi ) / (1 + exp( bstay[t] + ( truebetac_y_mat[m, 1] + truexb ) / phi ))
          est2 <- exp( bstay[t] + ( truebetac_y_mat[m, 2] + truexb ) / phi ) / (1 + exp( bstay[t] + ( truebetac_y_mat[m, 2] + truexb ) / phi ))
          est3 <- exp( bstay[t] + ( truebetac_y_mat[m, 3] + truexb ) / phi ) / (1 + exp( bstay[t] + ( truebetac_y_mat[m, 3] + truexb ) / phi ))
          
          prob1 <- est1
          prob2 <- est2-est1
          prob3 <- est3-est2
          prob4 <- 1-est3
          
          theprob <- c(prob1, prob2, prob3, prob4)
          ydummy <- rmultinom( n = 1, size = 1, prob = theprob )
          
          yvec[[t]] <- catvec %*% ydummy
          
        }
        ymat[,m] <- yvec
      }
      ally[[i]] <-ymat
      nstay <- 0
      
    }
    
    y_mat <- do.call(rbind, ally)
    colnames(y_mat) <- c("Y1", "Y2", "Y3", "Y4")
    y1_table <- table(y_mat[,1])
  }
  
  
  x_mat <- cbind(rep(xi, times = nteethpp), rep(zi, times= nteethpp))
  colnames(x_mat) <- c("X", "Z")
  
  subject <- rep(1:N, times = nteethpp)
  
  toothlist <- vector("list", N)
  for(i in 1:N) toothlist[[i]] <- rep(1:nteethpp[i])
  tooth <- unlist(toothlist)
  
  blnteethlist <- vector("list", N)
  for(i in 1:N) blnteethlist[[i]] <- rep(nteethpp[i], each = nteethpp[i])
  blnteeth <- unlist(blnteethlist)
  
  simdat <- as.data.frame(cbind(subject, tooth, x_mat, y_mat, blnteeth))
  simdat <- simdat %>% mutate(Y1=Y1, Y2=Y2, Y3=Y3, Y4=Y4)
  
  cluster <- simdat$subject
  unit <- simdat$tooth
  Y <- simdat$Y1
  X <- cbind(rep(1, length(Y)), simdat$X, simdat$Z)
  
  #----------------------------------------------------
  # check simulated data
  #----------------------------------------------------
  ##extract the Y1 out from ally
  ally_Y1 <- rapply(ally, classes = 'matrix', how = 'list', f = function(x) x[, 1, drop = FALSE])
  ymean <- unlist(lapply(ally_Y1, mean), use.names = F)
  # 
  # plot(baserisk, ymean)
  #plot(nteethpp, ymean)
  cor.ics = cor(nteethpp,ymean, method="pearson")
  
  # hist(nteethpp)
  # 
  # meany <- mean(Y)
  # meanx <- mean(x)
  # mediannteethpp <- median(nteethpp)
  # meannteethpp <- mean(nteethpp)
  
  #----------------------------------------------------
  ###  Analysis of Full data 
  #----------------------------------------------------
  
  
  ### cluster weighted GEE
  wgee_full <-  ordgee(ordered(Y1)~X+factor(Z), data=simdat, id=subject, mean.link = "logit", weights = blnteeth^(-1), corstr = "independence", rev=T)
  summ_wgee_full <- summary(wgee_full)
  est_wgee_full <- summ_wgee_full$mean$estimate
  se_wgee_full <- summ_wgee_full$mean$san.se
  
  est_wgee_full_int1 <- est_wgee_full[1]
  se_wgee_full_int1 <- se_wgee_full[1]
  est_wgee_full_int5 <- est_wgee_full[5]
  se_wgee_full_int5 <- se_wgee_full[5]
}


#----------------------------------------------------
# Generate missing patterns (MAR, MNAR, MCAR)
#----------------------------------------------------

alpha_missing <- rep(0.5, 6)

simdat_miss <- generate_NA(data=simdat, missing = c(rate, 0.3, 0.3, 0.1), alpha=alpha_missing, pattern=pattern, start = -10,  by=0.01, diff=0.01)

simdat_misy <- simdat_miss$data_missing 
## Currently only numeric missing codes are allowed in blimp
simdat_misy_blimp <- simdat_misy
simdat_misy_blimp[is.na(simdat_misy_blimp)] <- 99999

dir.create(file.path(getwd(), "output"))
write.table(simdat_misy_blimp, file=paste0("output/missing_rate", rate,  "_rep", rep, ".txt"), row.names = F, sep=",", col.names = F)
#----------------------------------------------------
# Complete case analysis
#----------------------------------------------------
##get missing rate
true_missing <- length(which(is.na(simdat_misy$Y1)))/length(simdat_misy$Y1)

### cluster weighted GEE
wgee_cca <-  ordgee(ordered(Y1)~X+factor(Z), data=simdat_misy, id=subject, mean.link = "logit", weights = blnteeth^(-1), corstr = "independence", rev=T)
summ_wgee_cca <- summary(wgee_cca)
est_wgee_cca <- summ_wgee_cca$mean$estimate
se_wgee_cca <- summ_wgee_cca$mean$san.se


#----------------------------------------------------
# Imputation with joint model
# joint model modeling ordinal outcome by treating it as categorical data
#----------------------------------------------------
Y_jomo = simdat_misy %>% mutate(Y1=as.factor(Y1), Y2=as.factor(Y2), Y3=as.factor(Y3), Y4=as.factor(Y4))  %>% select(Y1, Y2, Y3, Y4)
#imputation without cluster size
X_jomo=data.frame(rep(1,dim(simdat_misy)[1]),simdat_misy[,c("X", "Z")])

colnames(X_jomo)[1] <- "const"

imp <- jomo(Y=Y_jomo, X=X_jomo, clus = simdat_misy$subject, nburn = 2000, nbetween = 1000, nimp = 5)

imp.list <- imputationList(split(imp, imp$Imputation)[-1])

### analysis with cluster size, weighted GEE
fit.imp.wgee <- with(data = imp.list, ordgee(ordered(Y1)~X+factor(Z), id=subject, mean.link = "logit", weights = blnteeth^(-1), corstr = "independence", rev=T))
MI_est_wgee_jomo_lv <- MIextract(fit.imp.wgee, fun = function(x) x[["beta"]])
MI_see_wgee_jomo_lv <- MIextract(fit.imp.wgee, fun = function(x) x[["vbeta"]])
# Pool results with Rubin's rules
results_wgee_jomo_lv <- MIcombine(MI_est_wgee_jomo_lv, MI_see_wgee_jomo_lv)
est_wgee_jomo_lv <- summary(results_wgee_jomo_lv)$results
se_wgee_jomo_lv <- summary(results_wgee_jomo_lv)$se


#imputation with cluster size, denoted by _wN
X_jomo_wN=data.frame(rep(1,dim(simdat_misy)[1]),simdat_misy[,c("X", "Z", "blnteeth")])

colnames(X_jomo_wN)[1] <- "const"

imp_wN <- jomo(Y=Y_jomo, X=X_jomo_wN, clus = simdat_misy$subject, nburn = 2000, nbetween = 1000, nimp = 5)

imp.list_wN <- imputationList(split(imp_wN, imp_wN$Imputation)[-1])

### analysis with cluster size, weighted GEE
fit.imp.wgee_wN <- with(data = imp.list_wN, ordgee(ordered(Y1)~X+factor(Z), id=subject, mean.link = "logit", weights = (simdat_misy$blnteeth)^(-1), corstr = "independence", rev=T))
MI_est_wgee_jomo_lv_wN <- MIextract(fit.imp.wgee_wN, fun = function(x) x[["beta"]])
MI_see_wgee_jomo_lv_wN <- MIextract(fit.imp.wgee_wN, fun = function(x) x[["vbeta"]])
# Pool results with Rubin's rules
results_wgee_jomo_lv_wN <- MIcombine(MI_est_wgee_jomo_lv_wN, MI_see_wgee_jomo_lv_wN)
est_wgee_jomo_lv_wN <- summary(results_wgee_jomo_lv_wN)$results
se_wgee_jomo_lv_wN <- summary(results_wgee_jomo_lv_wN)$se

#----------------------------------------------------
# Imputation with Fully conditional specification
# use system function to integrate blimp
# you will need to have blimp installed
# the blimp version is blimp_binary_dynlibs for linux
# the following code need to run in linux environment
#----------------------------------------------------


blimpPath <- "~/software/blimp_binary_dynlibs/blimp"

## need to write the blimp file first for the imputation model
inputPath <- "FCS.imp"

dataPath <- paste0(getwd(),"/output/missing_rate", rate,  "_rep", rep, ".txt")

outputPath <- paste0(getwd(), "/output/imputed_blimp_rate", rate,  "_rep", rep, ".txt")

##run blimp
out <- system(paste(blimpPath,inputPath,"-o",outputPath,"-s", rep,"-d",dataPath),intern = F)

imputation_fcs <- read.table(paste0(getwd(), "/output/imputed_blimp_rate", rate,  "_rep", rep, ".txt"))
colnames(imputation_fcs) <- c("Imputation", colnames(simdat))

imp.list_fcs <- imputationList(split(imputation_fcs, imputation_fcs$Imputation)[-1])


### weighted GEE
fit.imp.N_fcs <- with(data = imp.list_fcs, ordgee(ordered(Y1)~X+factor(Z), id=subject, mean.link = "logit", weights = (simdat_misy$blnteeth)^(-1), corstr = "independence", rev=T))
#summary(fit.imp.N_fcs[[1]])
MI_est_wgee_fcs <- MIextract(fit.imp.N_fcs, fun = function(x) x[["beta"]])
MI_see_wgee_fcs <- MIextract(fit.imp.N_fcs, fun = function(x) x[["vbeta"]])
# Pool results with Rubin's rules
results_wgee_fcs <- MIcombine(MI_est_wgee_fcs, MI_see_wgee_fcs)
est_wgee_fcs <- summary(results_wgee_fcs)$results
se_wgee_fcs <- summary(results_wgee_fcs)$se


###imputation with cluster size
inputPath2 <- "FCS_wN.imp"
outputPath2 <-  paste0(getwd(), "/output/imputed_blimp_wN_rate", rate,  "_rep", rep, ".txt")
##run blimp
out <- system(paste(blimpPath,inputPath2,"-o",outputPath2,"-s", 2022,"-d",dataPath),intern = F)

imputation_fcs_wN <- read.table(paste0(getwd(), "/output/imputed_blimp_wN_rate", rate,  "_rep", rep, ".txt"))
colnames(imputation_fcs_wN) <- c("Imputation", colnames(simdat))
imp.list_fcs_wN <- imputationList(split(imputation_fcs_wN, imputation_fcs_wN$Imputation)[-1])


### weighted GEE
fit.imp.N_fcs_wN <- with(data = imp.list_fcs_wN, ordgee(ordered(Y1)~X+factor(Z), id=subject, mean.link = "logit", weights = blnteeth^(-1), corstr = "independence", rev=T))
summary(fit.imp.N_fcs_wN[[1]])
MI_est_wgee_fcs_wN <- MIextract(fit.imp.N_fcs_wN, fun = function(x) x[["beta"]])
MI_see_wgee_fcs_wN <- MIextract(fit.imp.N_fcs_wN, fun = function(x) x[["vbeta"]])
# Pool results with Rubin's rules
results_wgee_fcs_wN <- MIcombine(MI_est_wgee_fcs_wN, MI_see_wgee_fcs_wN)
est_wgee_fcs_wN <- summary(results_wgee_fcs_wN)$results
se_wgee_fcs_wN <- summary(results_wgee_fcs_wN)$se


