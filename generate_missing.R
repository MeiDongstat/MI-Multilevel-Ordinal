###this is a function to generate missing values and calculate missing rate

# Code Written by: Mei Dong

generate_missing <- function(truealpha0, prob_noint){
  prob_mar <- truealpha0+prob_noint
  missprob <- exp(prob_mar) / ( 1 + exp(prob_mar) )
  missind <- rbinom(length(missprob), 1, missprob)
  return(missind)
}


## a function to generate desired missing rate for one predictor
## need to specify missing rate, missing pattern alphas, dataset, which predictor, 
generate_NA_one <- function(rate=0.1, alpha=rep(0, 6), data=simdat, predictor="Y1", start=-10,  by=0.01, diff=0.001,seed=2022){
  set.seed(seed)
  prob_noint <- as.matrix(data[, c("X", "Z", "Y1", "Y2", "Y3", "Y4")])  %*% alpha
 
  truealpha0 <- start
  missind <- generate_missing(truealpha0=truealpha0, prob_noint=prob_noint)
  missrate <- sum(missind)/length(missind)
  while(abs(missrate-rate) > diff & truealpha0 <= 20){
      truealpha0 <- truealpha0+by
      missind <- generate_missing(truealpha0=truealpha0, prob_noint=prob_noint)
      missrate <- sum(missind)/length(missind)
  }
  
  # if(truealpha0_mar>=20){
  #   by2=by/2
  #   truealpha0_mar <- start
  #   missrate <- generate_missing(truealpha0=truealpha0_mar, prob_noint=prob_noint)$missrate
  #   while(abs(missrate-rate) > diff2 & truealpha0_mar <= 20){
  #     truealpha0_mar <- truealpha0_mar+by2
  #     missrate <- generate_missing(truealpha0=truealpha0_mar, prob_noint=prob_noint)$missrate
  #     missind <- generate_missing(truealpha0=truealpha0_mar, prob_noint=prob_noint)$missind
  #   }
  # }
  print(truealpha0)
  idx_rm <- which(missind==1)
  return(idx_rm <- idx_rm)
}

generate_NA <- function(data=simdat, missing=c(0.2, 0.3, 0.3, 0.1), alpha=rep(0, 6), pattern="MCAR", start=-10,  by=0.01, diff=0.001){
  data_miss <- data  
  for(i in 1:4){
    alpha_tem <- alpha
    predictor <- paste0("Y", i)
    if(pattern=="MCAR"){
      alpha_tem <- rep(0, 6)
    } else if(pattern=="MAR"){
      alpha_tem[i+2] = 0
    } 
    idx_rm <- generate_NA_one(rate=missing[i], data=data, alpha=alpha_tem, predictor=predictor, start=start, by=by, diff=diff)
    data_miss[idx_rm, predictor] <- NA
  }
  return(list(data=data, data_missing=data_miss))
}
