library(ggplot2)
library(mice)
library(sm)

setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")

set.seed(190)

#####################################################################################
## Function that generates a simulated population
#####################################################################################

gen_population_data <- function(){
  N <- 500000
  
  #Generate covariates
  smoking <- sample(c(0,1), N, replace = TRUE, prob = c(0.8, 0.2))
  gender <- sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5))
  
  #Generate pollutant given covariates: X|C
  sd_genx <- 1.15
  gamma0 <- -0.5
  gamma1 <- 1.25
  gamma2 <- 1.25
  logpoll <- gamma0 + gamma1*smoking + gamma2*gender + rnorm(N, mean = 0, sd = sd_genx)
  poll <- exp(logpoll)
  
  #Coefficients to generate Y|X,C
  expbeta1 <- 1.5
  expbeta2 <- 1.25
  expbeta3 <- 1.25
  
  expbeta0 <- exp(-2.68) #Determined by setting P(Y=1) = 0.1 and solving for Beta0
  
  #Generate Case/Control Status: Y|X,C
  y <- log(expbeta0) + log(expbeta1)*logpoll + log(expbeta2)*smoking + log(expbeta3)*gender
  
  #True Probabilities with noise
  probability <- 1/(1+exp(-y))
  case_cntrl <- rbinom(N, 1, probability)
  sum(case_cntrl)
  
  #Create Dataset and vector of cases and cntrls to sample from
  data <- as.data.frame(cbind(id=1:N,case_cntrl,logpoll,poll,smoking,gender))
  cases <- subset(data, case_cntrl == 1)
  cntrls <- subset(data, case_cntrl == 0)
  
  #Other checks
  model <- glm(case_cntrl~logpoll+smoking+gender, family="binomial")
  model2 <- lm(logpoll~smoking+gender+case_cntrl)
  
  return(list(cases, cntrls, model, model2))
}

#Generate Study population and true model for the entire population

pop <- gen_population_data()
cases <- pop[[1]]
cntrls <- pop[[2]]
truth <- pop[[4]]

###################################################################################################
## Function that checks whether the MLE estimated by CLMI coincides with the population parameters.
##
## Arguments:
##   wb: Stored simulation results when a batch indicator is included in the analysis model
##   wob: Stored simulation results when a batch indicator is not included in the analysis model
###################################################################################################

check_mle <- function(wb, wob){
  gamma0_bias <- rep(0,9)
  gamma1_bias <- rep(0,9)
  gamma2_bias <- rep(0,9)
  gamma3_bias <- rep(0,9)
  sigma2_bias <- rep(0,9)
  
  gamma0_bias_wob <- rep(0,9)
  gamma1_bias_wob <- rep(0,9)
  gamma2_bias_wob <- rep(0,9)
  gamma3_bias_wob <- rep(0,9)
  sigma2_bias_wob <- rep(0,9)
  
  gamma0_se <- rep(0,9)
  gamma1_se <- rep(0,9)
  gamma2_se <- rep(0,9)
  gamma3_se <- rep(0,9)
  sigma2_se <- rep(0,9)
  
  gamma0_se_wob <- rep(0,9)
  gamma1_se_wob <- rep(0,9)
  gamma2_se_wob <- rep(0,9)
  gamma3_se_wob <- rep(0,9)
  sigma2_se_wob <- rep(0,9)
  
  gamma0_covprob <- rep(0,9)
  gamma1_covprob <- rep(0,9)
  gamma2_covprob <- rep(0,9)
  gamma3_covprob <- rep(0,9)
  sigma2_covprob <- rep(0,9)
  
  gamma0_covprob_wob <- rep(0,9)
  gamma1_covprob_wob <- rep(0,9)
  gamma2_covprob_wob <- rep(0,9)
  gamma3_covprob_wob <- rep(0,9)
  sigma2_covprob_wob <- rep(0,9)
  
  #Add % Below LOD Variable
  wb$LOD1Per <- c(rep(15,3000),rep(30,3000),rep(60,3000))
  wb$LOD2Per <- c(rep(15,1000),rep(30,1000),rep(60,1000),
                  rep(15,1000),rep(30,1000),rep(60,1000),
                  rep(15,1000),rep(30,1000),rep(60,1000))
  
  wob$LOD1Per <- c(rep(15,3000),rep(30,3000),rep(60,3000))
  wob$LOD2Per <- c(rep(15,1000),rep(30,1000),rep(60,1000),
                  rep(15,1000),rep(30,1000),rep(60,1000),
                  rep(15,1000),rep(30,1000),rep(60,1000))
  
  lod1 <- rep(0,9)
  lod2 <- rep(0,9)
  lods <- c(15,30,60)
  cnt <- 0
  for(l1 in lods){
    for(l2 in lods){
      cnt <- cnt + 1
      lod1[cnt] <- l1
      lod2[cnt] <- l2
      sub <- subset(wb, LOD1Per == l1 & LOD2Per == l2)
      sub_wob <- subset(wob, LOD1Per == l1 & LOD2Per == l2)
      
      gamma0_bias[cnt] <- mean(sub$PMIMLE0-truth$coefficients[1])
      gamma1_bias[cnt] <- mean(sub$PMIMLE1-truth$coefficients[2])
      gamma2_bias[cnt] <- mean(sub$PMIMLE2-truth$coefficients[3])
      gamma3_bias[cnt] <- mean(sub$PMIMLE3-truth$coefficients[4])
      sigma2_bias[cnt] <- mean(sub$PMIMLE4-(summary(truth)$sigma)^2)
      
      gamma0_bias_wob[cnt] <- mean(sub_wob$PMIMLE0-truth$coefficients[1])
      gamma1_bias_wob[cnt] <- mean(sub_wob$PMIMLE1-truth$coefficients[2])
      gamma2_bias_wob[cnt] <- mean(sub_wob$PMIMLE2-truth$coefficients[3])
      gamma3_bias_wob[cnt] <- mean(sub_wob$PMIMLE3-truth$coefficients[4])
      sigma2_bias_wob[cnt] <- mean(sub_wob$PMIMLE4-(summary(truth)$sigma)^2)
      
      gamma0_se[cnt] <- mean(sub$PMIMLEStdErr0)
      gamma1_se[cnt] <- mean(sub$PMIMLEStdErr1)
      gamma2_se[cnt] <- mean(sub$PMIMLEStdErr2)
      gamma3_se[cnt] <- mean(sub$PMIMLEStdErr3)
      sigma2_se[cnt] <- mean(sub$PMIMLEStdErr4)
      
      gamma0_se_wob[cnt] <- mean(sub_wob$PMIMLEStdErr0)
      gamma1_se_wob[cnt] <- mean(sub_wob$PMIMLEStdErr1)
      gamma2_se_wob[cnt] <- mean(sub_wob$PMIMLEStdErr2)
      gamma3_se_wob[cnt] <- mean(sub_wob$PMIMLEStdErr3)
      sigma2_se_wob[cnt] <- mean(sub_wob$PMIMLEStdErr4)
      
      lower_g0 <- sub$PMIMLE0 - 1.96*sub$PMIMLEStdErr0
      upper_g0 <- sub$PMIMLE0 + 1.96*sub$PMIMLEStdErr0
      lower_g1 <- sub$PMIMLE1 - 1.96*sub$PMIMLEStdErr1
      upper_g1 <- sub$PMIMLE1 + 1.96*sub$PMIMLEStdErr1
      lower_g2 <- sub$PMIMLE2 - 1.96*sub$PMIMLEStdErr2
      upper_g2 <- sub$PMIMLE2 + 1.96*sub$PMIMLEStdErr2
      lower_g3 <- sub$PMIMLE3 - 1.96*sub$PMIMLEStdErr3
      upper_g3 <- sub$PMIMLE3 + 1.96*sub$PMIMLEStdErr3
      lower_s2 <- (sub$PMIMLE4*2*(sub$PMIMLE4/sub$PMIMLEStdErr4)^2)/qchisq(.975, df=2*(sub$PMIMLE4/sub$PMIMLEStdErr4)^2)
      upper_s2 <- (sub$PMIMLE4*2*(sub$PMIMLE4/sub$PMIMLEStdErr4)^2)/qchisq(.025, df=2*(sub$PMIMLE4/sub$PMIMLEStdErr4)^2)
      
      gamma0_covprob[cnt] <- sum(lower_g0 < truth$coefficients[1] & upper_g0 > truth$coefficients[1])/nrow(sub)
      gamma1_covprob[cnt] <- sum(lower_g1 < truth$coefficients[2] & upper_g1 > truth$coefficients[2])/nrow(sub)
      gamma2_covprob[cnt] <- sum(lower_g2 < truth$coefficients[3] & upper_g2 > truth$coefficients[3])/nrow(sub)
      gamma3_covprob[cnt] <- sum(lower_g3 < truth$coefficients[4] & upper_g3 > truth$coefficients[4])/nrow(sub)
      sigma2_covprob[cnt] <- sum(lower_s2 < (summary(truth)$sigma)^2 & upper_s2 > (summary(truth)$sigma)^2)/nrow(sub)
      
      lower_g0_wob <- sub_wob$PMIMLE0 - 1.96*sub_wob$PMIMLEStdErr0
      upper_g0_wob <- sub_wob$PMIMLE0 + 1.96*sub_wob$PMIMLEStdErr0
      lower_g1_wob <- sub_wob$PMIMLE1 - 1.96*sub_wob$PMIMLEStdErr1
      upper_g1_wob <- sub_wob$PMIMLE1 + 1.96*sub_wob$PMIMLEStdErr1
      lower_g2_wob <- sub_wob$PMIMLE2 - 1.96*sub_wob$PMIMLEStdErr2
      upper_g2_wob <- sub_wob$PMIMLE2 + 1.96*sub_wob$PMIMLEStdErr2
      lower_g3_wob <- sub_wob$PMIMLE3 - 1.96*sub_wob$PMIMLEStdErr3
      upper_g3_wob <- sub_wob$PMIMLE3 + 1.96*sub_wob$PMIMLEStdErr3
      lower_s2_wob <- (sub_wob$PMIMLE4*2*(sub_wob$PMIMLE4/sub_wob$PMIMLEStdErr4)^2)/qchisq(.975, df=2*(sub_wob$PMIMLE4/sub_wob$PMIMLEStdErr4)^2)
      upper_s2_wob <- (sub_wob$PMIMLE4*2*(sub_wob$PMIMLE4/sub_wob$PMIMLEStdErr4)^2)/qchisq(.025, df=2*(sub_wob$PMIMLE4/sub_wob$PMIMLEStdErr4)^2)
      
      gamma0_covprob_wob[cnt] <- sum(lower_g0_wob < truth$coefficients[1] & upper_g0_wob > truth$coefficients[1])/nrow(sub_wob)
      gamma1_covprob_wob[cnt] <- sum(lower_g1_wob < truth$coefficients[2] & upper_g1_wob > truth$coefficients[2])/nrow(sub_wob)
      gamma2_covprob_wob[cnt] <- sum(lower_g2_wob < truth$coefficients[3] & upper_g2_wob > truth$coefficients[3])/nrow(sub_wob)
      gamma3_covprob_wob[cnt] <- sum(lower_g3_wob < truth$coefficients[4] & upper_g3_wob > truth$coefficients[4])/nrow(sub_wob)
      sigma2_covprob_wob[cnt] <- sum(lower_s2_wob < (summary(truth)$sigma)^2 & upper_s2_wob > (summary(truth)$sigma)^2)/nrow(sub_wob)
    }
  }
  summar_mle <- as.data.frame(cbind(lod1, lod2, gamma0_bias, gamma1_bias, gamma2_bias, gamma3_bias, sigma2_bias,
                                    gamma0_se, gamma1_se, gamma2_se, gamma3_se, sigma2_se,
                                    gamma0_covprob, gamma1_covprob, gamma2_covprob, gamma3_covprob, sigma2_covprob))
  
  summar_mle_wob <- as.data.frame(cbind(lod1, lod2, gamma0_bias_wob, gamma1_bias_wob, gamma2_bias_wob, gamma3_bias_wob, sigma2_bias_wob,
                                        gamma0_se_wob, gamma1_se_wob, gamma2_se_wob, gamma3_se_wob, sigma2_se_wob,
                                        gamma0_covprob_wob, gamma1_covprob_wob, gamma2_covprob_wob, gamma3_covprob_wob, sigma2_covprob_wob))
  
  return(list(summar_mle, summar_mle_wob))
}

#####################################################################################
## Read in Raw Simulation Results and check MLE determined from CLMI
#####################################################################################

#######################################################################################
# Large Cohort - Random Batch Assignment
#######################################################################################

wb <- read.table("largecohort_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecohort_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Moderate Cohort - Random Batch Assignment
#######################################################################################

wb <- read.table("modcohort_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcohort_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Large Cohort - Outcome Dependent Batch Assignment
#######################################################################################

wb <- read.table("largecohort_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecohort_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Moderate Cohort - Outcome Dependent Batch Assignment
#######################################################################################

wb <- read.table("modcohort_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcohort_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Large Cohort - Covariate Dependent Batch Assignment
#######################################################################################

wb <- read.table("largecohort_depbatchcov_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecohort_depbatchcov_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Moderate Cohort - Covariate Dependent Batch Assignment
#######################################################################################

wb <- read.table("modcohort_depbatchcov_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcohort_depbatchcov_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Large Case-Control - Random Batch Assignment
#######################################################################################

wb <- read.table("largecc_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecc_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Moderate Case-Control - Random Batch Assignment
#######################################################################################

wb <- read.table("modcc_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcc_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Large Case-Control - Dependent Batch Assignment
#######################################################################################

wb <- read.table("largecc_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecc_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)

#######################################################################################
# Moderate Case-Control - Dependent Batch Assignment
#######################################################################################

wb <- read.table("modcc_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcc_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- check_mle(wb = wb, wob = wob)



