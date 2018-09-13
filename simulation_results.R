library(ggplot2)
library(sm)
library(xlsx)

setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")

set.seed(190)

#####################################################################################
## Regenerate Population
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
  # print("Mean of the pollutant")
  # mean(poll)
  # print("Standard Deviation of the pollutant")
  # sd(poll)
  # print("Mean of the log-adjusted pollutant")
  # mean(logpoll)
  # print("Standard Deviation of the log_adjusted pollutant")
  # sd(logpoll)
  # summary(lm(logpoll~smoking+gender))
  
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
  #summary(glm(case_cntrl~logpoll+smoking+gender, family=binomial))
  
  # casef <- factor(case_cntrl, levels = c(0,1), labels = c("Control","Case"))
  # sm.density.compare(logpoll, case_cntrl, xlab = "Log-Adjusted Concentration")
  # title(main = "Contaminant Distribution by Case/Control Status")
  # 
  # colfill <- c(2:(2+length(levels(casef))))
  # legend("topright", levels(casef), fill=colfill)
  
  #Create Dataset and vector of cases and cntrls to sample from
  data <- as.data.frame(cbind(id=1:N,case_cntrl,logpoll,poll,smoking,gender))
  cases <- subset(data, case_cntrl == 1)
  cntrls <- subset(data, case_cntrl == 0)
  
  #Other checks
  model <- glm(case_cntrl~logpoll+smoking+gender, family="binomial")
  #summary(model)
  
  model2 <- lm(logpoll~smoking+gender+case_cntrl)
  
  return(list(cases, cntrls, model, model2))
}

pop <- gen_population_data()

model <- pop[[3]]

#####################################################################################
## Function that returns simulation results
#####################################################################################

sim.results.summarize <- function(wb, wob){
  #Organize LOD
  wb$lod_name <- c(rep("15 15",1000),rep("15 30",1000),rep("15 60",1000),rep("30 15",1000),
                   rep("30 30",1000),rep("30 60",1000),rep("60 15",1000),rep("60 30",1000),
                   rep("60 60",1000))
  wob$lod_name <- c(rep("15 15",1000),rep("15 30",1000),rep("15 60",1000),rep("30 15",1000),
                    rep("30 30",1000),rep("30 60",1000),rep("60 15",1000),rep("60 30",1000),
                    rep("60 60",1000))
  
  #Store Results in vectors
  gs_bias <- rep(0,9)
  gs_bias_wob <- rep(0,9)
  cca_bias <- rep(0,9)
  cca_bias_wob <- rep(0,9)
  sq2_bias <- rep(0,9)
  sq2_bias_wob <- rep(0,9)
  mice_bias <- rep(0,9)
  mice_bias_wob <- rep(0,9)
  pmi_bias <- rep(0,9)
  pmi_bias_wob <- rep(0,9)
  
  gs_rel_bias <- rep(0,9)
  gs_rel_bias_wob <- rep(0,9)
  cca_rel_bias <- rep(0,9)
  cca_rel_bias_wob <- rep(0,9)
  sq2_rel_bias <- rep(0,9)
  sq2_rel_bias_wob <- rep(0,9)
  mice_rel_bias <- rep(0,9)
  mice_rel_bias_wob <- rep(0,9)
  pmi_rel_bias <- rep(0,9)
  pmi_rel_bias_wob <- rep(0,9)
  
  gs_emp_sd <- rep(0,9)
  gs_emp_sd_wob <- rep(0,9)
  cca_emp_sd <- rep(0,9)
  cca_emp_sd_wob <- rep(0,9)
  sq2_emp_sd <- rep(0,9)
  sq2_emp_sd_wob <- rep(0,9)
  mice_emp_sd <- rep(0,9)
  mice_emp_sd_wob <- rep(0,9)
  pmi_emp_sd <- rep(0,9)
  pmi_emp_sd_wob <- rep(0,9)
  
  gs_avg_se <- rep(0,9)
  gs_avg_se_wob <- rep(0,9)
  cca_avg_se <- rep(0,9)
  cca_avg_se_wob <- rep(0,9)
  sq2_avg_se <- rep(0,9)
  sq2_avg_se_wob <- rep(0,9)
  mice_avg_se <- rep(0,9)
  mice_avg_se_wob <- rep(0,9)
  pmi_avg_se <- rep(0,9)
  pmi_avg_se_wob <- rep(0,9)
  
  gs_mse <- rep(0,9)
  gs_mse_wob <- rep(0,9)
  cca_mse <- rep(0,9)
  cca_mse_wob <- rep(0,9)
  sq2_mse <- rep(0,9)
  sq2_mse_wob <- rep(0,9)
  mice_mse <- rep(0,9)
  mice_mse_wob <- rep(0,9)
  pmi_mse <- rep(0,9)
  pmi_mse_wob <- rep(0,9)
  
  gs_cov <- rep(0,9)
  gs_cov_wob <- rep(0,9)
  cca_cov <- rep(0,9)
  cca_cov_wob <- rep(0,9)
  sq2_cov <- rep(0,9)
  sq2_cov_wob <- rep(0,9)
  mice_cov <- rep(0,9)
  mice_cov_wob <- rep(0,9)
  pmi_cov <- rep(0,9)
  pmi_cov_wob <- rep(0,9)
  
  #Store LODs
  lods <- rep(0,9)
  lod_pairs <- unique(wb$lod_name)
  cnt <- 0
  for(l in lod_pairs){
    cnt <- cnt + 1
    lods[cnt] <- l
    sub <- subset(wb, lod_name == l)
    sub_wob <- subset(wob, lod_name == l)
    
    #Bias
    gs_bias[cnt] <- mean(sub$GSBeta1-model$coefficients[2])
    gs_bias_wob[cnt] <- mean(sub_wob$GSBeta1-model$coefficients[2])
    cca_bias[cnt] <- mean(sub$CCABeta1-model$coefficients[2])
    cca_bias_wob[cnt] <- mean(sub_wob$CCABeta1-model$coefficients[2])
    sq2_bias[cnt] <- mean(sub$SQ2Beta1-model$coefficients[2])
    sq2_bias_wob[cnt] <- mean(sub_wob$SQ2Beta1-model$coefficients[2])
    mice_bias[cnt] <- mean(sub$MICEBeta1-model$coefficients[2])
    mice_bias_wob[cnt] <- mean(sub_wob$MICEBeta1-model$coefficients[2])
    pmi_bias[cnt] <- mean(sub$PMIBeta1-model$coefficients[2])
    pmi_bias_wob[cnt] <- mean(sub_wob$PMIBeta1-model$coefficients[2])
    
    #Relative Bias
    gs_rel_bias[cnt] <- 100*mean(sub$GSBeta1-model$coefficients[2])/model$coefficients[2]
    gs_rel_bias_wob[cnt] <- 100*mean(sub_wob$GSBeta1-model$coefficients[2])/model$coefficients[2]
    cca_rel_bias[cnt] <- 100*mean(sub$CCABeta1-model$coefficients[2])/model$coefficients[2]
    cca_rel_bias_wob[cnt] <- 100*mean(sub_wob$CCABeta1-model$coefficients[2])/model$coefficients[2]
    sq2_rel_bias[cnt] <- 100*mean(sub$SQ2Beta1-model$coefficients[2])/model$coefficients[2]
    sq2_rel_bias_wob[cnt] <- 100*mean(sub_wob$SQ2Beta1-model$coefficients[2])/model$coefficients[2]
    mice_rel_bias[cnt] <- 100*mean(sub$MICEBeta1-model$coefficients[2])/model$coefficients[2]
    mice_rel_bias_wob[cnt] <- 100*mean(sub_wob$MICEBeta1-model$coefficients[2])/model$coefficients[2]
    pmi_rel_bias[cnt] <- 100*mean(sub$PMIBeta1-model$coefficients[2])/model$coefficients[2]
    pmi_rel_bias_wob[cnt] <- 100*mean(sub_wob$PMIBeta1-model$coefficients[2])/model$coefficients[2]
    
    #Empirical Standard Deviation
    gs_emp_sd[cnt] <- sqrt(var(sub$GSBeta1))
    gs_emp_sd_wob[cnt] <- sqrt(var(sub_wob$GSBeta1))
    cca_emp_sd[cnt] <- sqrt(var(sub$CCABeta1))
    cca_emp_sd_wob[cnt] <- sqrt(var(sub_wob$CCABeta1))
    sq2_emp_sd[cnt] <- sqrt(var(sub$SQ2Beta1))
    sq2_emp_sd_wob[cnt] <- sqrt(var(sub_wob$SQ2Beta1))
    mice_emp_sd[cnt] <- sqrt(var(sub$MICEBeta1))
    mice_emp_sd_wob[cnt] <- sqrt(var(sub_wob$MICEBeta1))
    pmi_emp_sd[cnt] <- sqrt(var(sub$PMIBeta1))
    pmi_emp_sd_wob[cnt] <- sqrt(var(sub_wob$PMIBeta1))
    
    #Average Standard Error
    gs_avg_se[cnt] <- mean(sqrt(sub$GSVar1))
    gs_avg_se_wob[cnt] <- mean(sqrt(sub_wob$GSVar1))
    cca_avg_se[cnt] <- mean(sqrt(sub$CCAVar1))
    cca_avg_se_wob[cnt] <- mean(sqrt(sub_wob$CCAVar1))
    sq2_avg_se[cnt] <- mean(sqrt(sub$SQ2Var1))
    sq2_avg_se_wob[cnt] <- mean(sqrt(sub_wob$SQ2Var1))
    mice_avg_se[cnt] <- mean(sqrt(sub$MICEVar1))
    mice_avg_se_wob[cnt] <- mean(sqrt(sub_wob$MICEVar1))
    pmi_avg_se[cnt] <- mean(sqrt(sub$PMIVar1))
    pmi_avg_se_wob[cnt] <- mean(sqrt(sub_wob$PMIVar1))
    
    #Mean-Squared Error
    gs_mse[cnt] <- mean((sub$GSBeta1-model$coefficients[2])^2)
    gs_mse_wob[cnt] <- mean((sub_wob$GSBeta1-model$coefficients[2])^2)
    cca_mse[cnt] <- mean((sub$CCABeta1-model$coefficients[2])^2)
    cca_mse_wob[cnt] <- mean((sub_wob$CCABeta1-model$coefficients[2])^2)
    sq2_mse[cnt] <- mean((sub$SQ2Beta1-model$coefficients[2])^2)
    sq2_mse_wob[cnt] <- mean((sub_wob$SQ2Beta1-model$coefficients[2])^2)
    mice_mse[cnt] <- mean((sub$MICEBeta1-model$coefficients[2])^2)
    mice_mse_wob[cnt] <- mean((sub_wob$MICEBeta1-model$coefficients[2])^2)
    pmi_mse[cnt] <- mean((sub$PMIBeta1-model$coefficients[2])^2)
    pmi_mse_wob[cnt] <- mean((sub_wob$PMIBeta1-model$coefficients[2])^2)
    
    #95% Coverage Probability
    lower_gs <- sub$GSBeta1 - 1.96*sqrt(sub$GSVar1)
    upper_gs <- sub$GSBeta1 + 1.96*sqrt(sub$GSVar1)
    lower_cca <- sub$CCABeta1 - 1.96*sqrt(sub$CCAVar1)
    upper_cca <- sub$CCABeta1 + 1.96*sqrt(sub$CCAVar1)
    lower_sq2 <- sub$SQ2Beta1 - 1.96*sqrt(sub$SQ2Var1)
    upper_sq2 <- sub$SQ2Beta1 + 1.96*sqrt(sub$SQ2Var1)
    lower_mice <- sub$MICEBeta1 - qt(0.975, sub$MICEVMStar1)*sqrt(sub$MICEVar1)
    upper_mice <- sub$MICEBeta1 + qt(0.975, sub$MICEVMStar1)*sqrt(sub$MICEVar1)
    lower_pmi <- sub$PMIBeta1 - qt(0.975, sub$PMIVMStar1)*sqrt(sub$PMIVar1)
    upper_pmi <- sub$PMIBeta1 + qt(0.975, sub$PMIVMStar1)*sqrt(sub$PMIVar1)
    
    lower_gs_wob <- sub_wob$GSBeta1 - 1.96*sqrt(sub_wob$GSVar1)
    upper_gs_wob <- sub_wob$GSBeta1 + 1.96*sqrt(sub_wob$GSVar1)
    lower_cca_wob <- sub_wob$CCABeta1 - 1.96*sqrt(sub_wob$CCAVar1)
    upper_cca_wob <- sub_wob$CCABeta1 + 1.96*sqrt(sub_wob$CCAVar1)
    lower_sq2_wob <- sub_wob$SQ2Beta1 - 1.96*sqrt(sub_wob$SQ2Var1)
    upper_sq2_wob <- sub_wob$SQ2Beta1 + 1.96*sqrt(sub_wob$SQ2Var1)
    lower_mice_wob <- sub_wob$MICEBeta1 - qt(0.975, sub_wob$MICEVMStar1)*sqrt(sub_wob$MICEVar1)
    upper_mice_wob <- sub_wob$MICEBeta1 + qt(0.975, sub_wob$MICEVMStar1)*sqrt(sub_wob$MICEVar1)
    lower_pmi_wob <- sub_wob$PMIBeta1 - qt(0.975, sub_wob$PMIVMStar1)*sqrt(sub_wob$PMIVar1)
    upper_pmi_wob <- sub_wob$PMIBeta1 + qt(0.975, sub_wob$PMIVMStar1)*sqrt(sub_wob$PMIVar1)
    
    gs_cov[cnt] <- sum(lower_gs < model$coefficients[2] & upper_gs > model$coefficients[2])/nrow(sub)
    gs_cov_wob[cnt] <- sum(lower_gs_wob < model$coefficients[2] & upper_gs_wob > model$coefficients[2])/nrow(sub_wob)
    cca_cov[cnt] <- sum(lower_cca < model$coefficients[2] & upper_cca > model$coefficients[2])/nrow(sub)
    cca_cov_wob[cnt] <- sum(lower_cca_wob < model$coefficients[2] & upper_cca_wob > model$coefficients[2])/nrow(sub_wob)
    sq2_cov[cnt] <- sum(lower_sq2 < model$coefficients[2] & upper_sq2 > model$coefficients[2])/nrow(sub)
    sq2_cov_wob[cnt] <- sum(lower_sq2_wob < model$coefficients[2] & upper_sq2_wob > model$coefficients[2])/nrow(sub_wob)
    mice_cov[cnt] <- sum(lower_mice < model$coefficients[2] & upper_mice > model$coefficients[2])/nrow(sub)
    mice_cov_wob[cnt] <- sum(lower_mice_wob < model$coefficients[2] & upper_mice_wob > model$coefficients[2])/nrow(sub_wob)
    pmi_cov[cnt] <- sum(lower_pmi < model$coefficients[2] & upper_pmi > model$coefficients[2])/nrow(sub)
    pmi_cov_wob[cnt] <- sum(lower_pmi_wob < model$coefficients[2] & upper_pmi_wob > model$coefficients[2])/nrow(sub_wob)
  }
  
  summar <- as.data.frame(cbind(lods, gs_bias, cca_bias, sq2_bias, mice_bias, pmi_bias,
                                gs_rel_bias, cca_rel_bias, sq2_rel_bias, mice_rel_bias, pmi_rel_bias,
                                gs_emp_sd, cca_emp_sd, sq2_emp_sd, mice_emp_sd, pmi_emp_sd,
                                gs_avg_se, cca_avg_se, sq2_avg_se, mice_avg_se, pmi_avg_se,
                                gs_mse, cca_mse, sq2_mse, mice_mse, pmi_mse,
                                gs_cov, cca_cov, sq2_cov, mice_cov, pmi_cov))
  
  summar_wob <- as.data.frame(cbind(lods, gs_bias_wob, cca_bias_wob, sq2_bias_wob, mice_bias_wob, pmi_bias_wob,
                                    gs_rel_bias_wob, cca_rel_bias_wob, sq2_rel_bias_wob, mice_rel_bias_wob, pmi_rel_bias_wob,
                                    gs_emp_sd_wob, cca_emp_sd_wob, sq2_emp_sd_wob, mice_emp_sd_wob, pmi_emp_sd_wob,
                                    gs_avg_se_wob, cca_avg_se_wob, sq2_avg_se_wob, mice_avg_se_wob, pmi_avg_se_wob,
                                    gs_mse_wob, cca_mse_wob, sq2_mse_wob, mice_mse_wob, pmi_mse_wob,
                                    gs_cov_wob, cca_cov_wob, sq2_cov_wob, mice_cov_wob, pmi_cov_wob))
  
  names(summar)[1] <- "lod_name"
  names(summar_wob)[1] <- "lod_name"
  
  return(list(summar, summar_wob))
}

#####################################################################################
## Read in Raw Simulation Results
#####################################################################################

#######################################################################################
# Large Cohort - Random Batch Assignment
#######################################################################################

wb <- read.table("largecohort_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecohort_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "largecohort_randombatch_withind.xlsx")
write.xlsx(final_res[[2]], "largecohort_randombatch_withoutind.xlsx")

#######################################################################################
# Moderate Cohort - Random Batch Assignment
#######################################################################################

wb <- read.table("modcohort_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcohort_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "modcohort_randombatch_withind.xlsx")
write.xlsx(final_res[[2]], "modcohort_randombatch_withoutind.xlsx")

#######################################################################################
# Large Cohort - Dependent Batch Assignment
#######################################################################################

wb <- read.table("largecohort_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecohort_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "largecohort_depbatch_withind.xlsx")
write.xlsx(final_res[[2]], "largecohort_depbatch_withoutind.xlsx")

#######################################################################################
# Moderate Cohort - Dependent Batch Assignment
#######################################################################################

wb <- read.table("modcohort_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcohort_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "modcohort_depbatch_withind.xlsx")
write.xlsx(final_res[[2]], "modcohort_depbatch_withoutind.xlsx")

#######################################################################################
# Large Case-Control - Random Batch Assignment
#######################################################################################

wb <- read.table("largecc_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecc_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "largecc_randombatch_withind.xlsx")
write.xlsx(final_res[[2]], "largecc_randombatch_withoutind.xlsx")

#######################################################################################
# Moderate Case-Control - Random Batch Assignment
#######################################################################################

wb <- read.table("modcc_randombatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcc_randombatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "modcc_randombatch_withind.xlsx")
write.xlsx(final_res[[2]], "modcc_randombatch_withoutind.xlsx")

#######################################################################################
# Large Case-Control - Dependent Batch Assignment
#######################################################################################

wb <- read.table("largecc_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("largecc_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "largecc_depbatch_withind.xlsx")
write.xlsx(final_res[[2]], "largecc_depbatch_withoutind.xlsx")

#######################################################################################
# Moderate Case-Control - Dependent Batch Assignment
#######################################################################################

wb <- read.table("modcc_depbatch_withind.txt", sep=",", stringsAsFactors = FALSE)
wob <- read.table("modcc_depbatch_withoutind.txt", sep=",", stringsAsFactors = FALSE)
final_res <- sim.results.summarize(wb = wb, wob = wob)
write.xlsx(final_res[[1]], "modcc_depbatch_withind.xlsx")
write.xlsx(final_res[[2]], "modcc_depbatch_withoutind.xlsx")



