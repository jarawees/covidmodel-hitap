## Codes below are adpated from Khoury & Cromer - Predicting-Efficacy-Variant-Modified-Boosters-main Github
ProbRemainUninfected <- function(logTitre, logk, C50){
  1/(1 + exp(-exp(logk)*(logTitre - C50)))
}

LogisticModel_PercentUninfected <- function(mu_titre,sig_titre,logk,C50){
  NumInteration <- max(length(mu_titre),length(sig_titre),length(logk),length(C50))
  Output <- NULL
  
  if (length(C50)==1) C50 <- rep(C50,NumInteration)
  if (length(logk)==1) logk <- rep(logk,NumInteration)
  if (length(sig_titre)==1) sig_titre <- rep(sig_titre,NumInteration)
  if (length(mu_titre)==1) mu_titre <- rep(mu_titre,NumInteration)
  
  for (i in 1:NumInteration) {
    Step <- sig_titre[i]*0.001
    IntegralVector <- seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i] <- sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}

# Covariate matrix
cov_inf <- matrix(c(0.03106460, 0.010755914, 0.010755914, 0.005727749), nrow = 2, ncol = 2) # infection/symptomatic/disease protection
cov_sev <- matrix(c(0.02999148, 0.02972476, 0.0297247, 0.09938996), nrow = 2, ncol = 2) # severe protection
cov_sev <- as.matrix(forceSymmetric(cov_sev)) # fource a symmetric matrix

# n50 is the (log10) neutralisation level that provides an individual with 50% protective efficacy of COVID-19
# k is the parameter determining the steepness of the logistic relationship
k_n50_inf <- c(1.13, log10(0.2)) # infection/symptomatic/disease protection
k_n50_sev <- c(1.12, log10(0.03)) # severe infection protection

# Generate distribution in efficacy from 10,000 bootstraps for each neutralisation level
# Estimate 95% confidence limits using the percentile method (2.5 and 97.5 percentile for a regular 95% confidence interval)

NeuToEfficacy <- function(MeanSD_NAb, n = 10000, k_n50, cov, SE_data){
  # debug
  # n = 1000
  # k_n50 <- k_n50_inf
  # cov <- cov_inf
  # MeanSD_NAb <- test
  # SE_data <- max(test$SEM)
  
  N <- n # number of bootstraps
  ModelParamtemp <- rmvnorm(N, mean = k_n50, sigma = cov)
  SDrandom <- rnorm(N, mean = PooledSD, sd = SE_PooledSD)
  
  LowerBound <- 0.025
  UpperBound <- 0.975
  
  MeanSD_NAb$Efficacy <- NA
  #MeanSD_NAb$EfficacySD <- NA
  MeanSD_NAb$EfficacyLower <- NA
  MeanSD_NAb$EfficacyUpper <- NA
  
  pb <- progress_bar$new(total = nrow(MeanSD_NAb))
  
  for (i in 1:nrow(MeanSD_NAb)) {
    MeanRandom <- rnorm(N, mean = MeanSD_NAb$MeanRatio[i], sd = SE_data)
    tempEvaluateFunction <- LogisticModel_PercentUninfected(MeanRandom, SDrandom, ModelParamtemp[,1], ModelParamtemp[,2])
    MeanSD_NAb$EfficacyLower[i] <- 100*quantile(tempEvaluateFunction, LowerBound)
    MeanSD_NAb$EfficacyUpper[i] <- 100*quantile(tempEvaluateFunction, UpperBound)
    MeanSD_NAb$Efficacy[i] <- 100*LogisticModel_PercentUninfected(MeanSD_NAb$MeanRatio[i], PooledSD, k_n50[1], k_n50[2])
    pb$tick()
  }
  return(MeanSD_NAb)
}

GetNeuOverTime <- function(starting_neut, decay_rate, time){
  if(decay_rate > 0) {decay_rate = -decay_rate}
  current_neut = starting_neut*exp(decay_rate*time)
}

GenerateBootstrap <- function(MeanSD_NAb, n = 10000, SE_data){
  N <- n # number of bootstraps
  pb <- progress_bar$new(total = nrow(MeanSD_NAb))
  
  # create blank dataframe to store bootstrap NAb and efficacy values
  col_name <- paste(MeanSD_NAb$study_id, MeanSD_NAb$group, MeanSD_NAb$hybrid, MeanSD_NAb$variant, sep = "__")
  BootNAb <- data.frame(matrix(ncol = nrow(MeanSD_NAb), nrow = N, dimnames = list(NULL, col_name)))
  BootInf <- data.frame(matrix(ncol = nrow(MeanSD_NAb), nrow = N, dimnames = list(NULL, col_name)))
  BootSev <- data.frame(matrix(ncol = nrow(MeanSD_NAb), nrow = N, dimnames = list(NULL, col_name)))
  
  for (i in 1:nrow(MeanSD_NAb)) {
    MeanRandom <- rnorm(N, mean = MeanSD_NAb$MeanRatio[i], sd = SE_data)
    BootNAb[,i] <- MeanRandom
    BootInf[,i] <- 100*LogisticModel_PercentUninfected(MeanRandom, PooledSD, k_n50_inf[1], k_n50_inf[2])
    BootSev[,i] <- 100*LogisticModel_PercentUninfected(MeanRandom, PooledSD, k_n50_sev[1], k_n50_sev[2])
    pb$tick()
  }
  
  BootNAb$label <- "NAb"
  BootInf$label <- "Infection"
  BootSev$label <- "Severe"
  
  BootCombine <- rbind(BootNAb, BootInf, BootSev)
  return(BootCombine)
}
