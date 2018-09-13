library(ggplot2)
library(mice)
library(sm)

#Version

setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")

set.seed(190)

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

gen_lods <- function(cases, cntrls, design, large){
  if(large == TRUE & (design == "Case Control - Dependent Batch" | design == "Case Control - Random Batch")){
    ncases <- 1000
    ncntrls <- 4000
  } else if(large == FALSE & (design == "Case Control - Dependent Batch" | design == "Case Control - Random Batch")){
    ncases <- 500
    ncntrls <- 500
  } else if(large == TRUE & (design == "Cohort - Dependent Batch" | design == "Cohort - Random Batch")){
    ncohort <- 5000
  } else if(large == FALSE & (design == "Cohort - Dependent Batch" | design == "Cohort - Random Batch")){
    ncohort <- 1000
  }
  #If we have dependent batch assignment need to calculate LOD separately for batch 1 and batch 2
  if(design == "Case Control - Dependent Batch" | design == "Cohort - Dependent Batch"){
    lod_15_b1 <- rep(0,1000)
    lod_30_b1 <- rep(0,1000)
    lod_60_b1 <- rep(0,1000)
    lod_15_b2 <- rep(0,1000)
    lod_30_b2 <- rep(0,1000)
    lod_60_b2 <- rep(0,1000)
  #If we have random batch assignment can calculate LOD irrespective of batch
  } else if(design == "Case Control - Random Batch" | design == "Cohort - Random Batch"){
    lod_15 <- rep(0,1000)
    lod_30 <- rep(0,1000)
    lod_60 <- rep(0,1000)
  }
  #Create population dataset for the cohort designs
  if(design == "Cohort - Dependent Batch" | design == "Cohort - Random Batch"){
    total_pop <- rbind(cases, cntrls)
  }
  for(s in 1:1000){
    print(s)
    if(design == "Case Control - Dependent Batch"){
      cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
      cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
      dataset <- rbind(cases_samp, cntrls_samp)
      expalpha0 <- exp(-0.3) #Determined by setting P(Y=1) = 0.5 and solving for Beta0
      expalpha1 <- 3
      batch <- log(expalpha0) + log(expalpha1)*dataset$case_cntrl
      #True Probabilities with noise
      prob <- 1/(1+exp(-batch))
      batch1 <- rbinom(length(batch), 1, prob)
      dataset$batch1 <- batch1
      lod_15_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.15)
      lod_15_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.15)
      lod_30_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.3)
      lod_30_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.3)
      lod_60_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.6)
      lod_60_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.6)
    } else if(design == "Cohort - Dependent Batch"){
      dataset <- subset(total_pop, id %in% sample(total_pop$id, ncohort, replace = FALSE))
      expalpha0 <- exp(0.73) #Determined by setting P(Y=1) = 0.5 and solving for Beta0
      expalpha1 <- exp(1.5)
      expalpha2 <- exp(-2)
      batch <- log(expalpha0) + log(expalpha1)*dataset$smoking + log(expalpha2)*dataset$gender
      #True Probabilities with noise
      prob <- 1/(1+exp(-batch))
      batch1 <- rbinom(length(batch), 1, prob)
      #sum(batch1)
      dataset$batch1 <- batch1
      lod_15_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.15)
      lod_15_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.15)
      lod_30_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.3)
      lod_30_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.3)
      lod_60_b1[s] <- quantile(subset(dataset, batch1 == 1)$logpoll, probs = 0.6)
      lod_60_b2[s] <- quantile(subset(dataset, batch1 == 0)$logpoll, probs = 0.6)
    } else if(design == "Case Control - Random Batch"){
      cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
      cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
      dataset <- rbind(cases_samp, cntrls_samp)
      lod_15[s] <- quantile(dataset$logpoll, probs = 0.15)
      lod_30[s] <- quantile(dataset$logpoll, probs = 0.3)
      lod_60[s] <- quantile(dataset$logpoll, probs = 0.6)
    } else if(design == "Cohort - Random Batch"){
      lod_15[s] <- quantile(total_pop$logpoll, probs = 0.15)
      lod_30[s] <- quantile(total_pop$logpoll, probs = 0.3)
      lod_60[s] <- quantile(total_pop$logpoll, probs = 0.6)
    }
  }
  if(design == "Case Control - Dependent Batch" | design == "Cohort - Dependent Batch"){
    b1_15 <- mean(lod_15_b1)
    b2_15 <- mean(lod_15_b2)
    b1_30 <- mean(lod_30_b1)
    b2_30 <- mean(lod_30_b2)
    b1_60 <- mean(lod_60_b1)
    b2_60 <- mean(lod_60_b2)
    
    #LODs
    lod_b1 <- c(b1_15, b1_30, b1_60)
    lod_b2 <- c(b2_15, b2_30, b2_60)
    lod_mat <- matrix(rep(0,9), nrow = 9, ncol = 2)
    cnt <- 1
    for(lod1_mat in lod_b1){
      for(lod2_mat in lod_b2){
        lod_mat[cnt,] <- c(lod1_mat, lod2_mat)
        cnt <- cnt + 1
      }
    }
  } else if(design == "Case Control - Random Batch" | design == "Cohort - Random Batch"){
    #LODs
    lod <- c(mean(lod_15), mean(lod_30), mean(lod_60))
    lod_mat <- matrix(rep(0,9), nrow = 9, ncol = 2)
    cnt <- 1
    for(lod1_mat in lod){
      for(lod2_mat in lod){
        lod_mat[cnt,] <- c(lod1_mat, lod2_mat)
        cnt <- cnt + 1
      }
    }
  }
  return(lod_mat)
}

gold_standard <- function(dataset){
  
  #With batch indicator
  gs <- glm(data=dataset, case_cntrl~logpoll+smoking+gender+batch1, family="binomial")
  gs_coeff <- gs$coefficients
  gs_var <- (diag(summary(gs)$cov.unscaled)*summary(gs)$dispersion)
  
  #Without batch indicator
  gs_wob <- glm(data=dataset, case_cntrl~logpoll+smoking+gender, family="binomial")
  gs_coeff_wob <- gs_wob$coefficients
  gs_var_wob <- (diag(summary(gs_wob)$cov.unscaled)*summary(gs_wob)$dispersion)
  
  return(list(gs_coeff, gs_var, gs_coeff_wob, gs_var_wob))
}

complete_case_analysis <- function(dataset){
  
  cca_data <- na.omit(dataset)
  
  #With batch indicator
  cca <- glm(data=cca_data, case_cntrl~logpoll_cens+smoking+gender+batch1, family="binomial")
  cca_coeff <- cca$coefficients
  cca_var <- (diag(summary(cca)$cov.unscaled)*summary(cca)$dispersion)
  
  #Without batch indicator
  cca_wob <- glm(data=cca_data, case_cntrl~logpoll_cens+smoking+gender, family="binomial")
  cca_coeff_wob <- cca_wob$coefficients
  cca_var_wob <- (diag(summary(cca_wob)$cov.unscaled)*summary(cca_wob)$dispersion)
  
  return(list(cca_coeff, cca_var, cca_coeff_wob, cca_var_wob))
}

lod_sqrt2_substitution <- function(dataset, lod1, lod2){
  
  dataset$logpoll_cens[dataset$batch1 == 1 & is.na(dataset$logpoll_cens)] <- log(lod1/sqrt(2))
  dataset$logpoll_cens[dataset$batch1 == 0 & is.na(dataset$logpoll_cens)] <- log(lod2/sqrt(2))
  
  #With batch indicator
  lod2mod <- glm(data=dataset, case_cntrl~logpoll_cens+smoking+gender+batch1, family = "binomial")
  lod2mod_coeff <- lod2mod$coefficients
  lod2mod_var <- (diag(summary(lod2mod)$cov.unscaled)*summary(lod2mod)$dispersion)
  
  #Without batch indicator
  lod2mod_wob <- glm(data=dataset, case_cntrl~logpoll_cens+smoking+gender, family = "binomial")
  lod2mod_coeff_wob <- lod2mod_wob$coefficients
  lod2mod_var_wob <- (diag(summary(lod2mod_wob)$cov.unscaled)*summary(lod2mod_wob)$dispersion)
  
  return(list(lod2mod_coeff, lod2mod_var, lod2mod_coeff_wob, lod2mod_var_wob))
}

mice_imputation <- function(dataset){
  
  #Imputation Model Variables
  imputekeep <- c("case_cntrl","logpoll_cens","smoking","gender","batch1")
  micedata <- dataset[,names(dataset) %in% imputekeep]
  mi <- mice(micedata,m=5,printFlag = FALSE)
  
  #With Batch
  imputed_together <- with(mi, glm(case_cntrl~logpoll_cens+smoking+gender+batch1, family="binomial"))
  imp_tog <- pool(imputed_together)
  qbar <- imp_tog$qbar
  t <- diag(imp_tog$t)
  vmstar <- imp_tog$df
  
  #Without Batch
  imputed_together2 <- with(mi, glm(case_cntrl~logpoll_cens+smoking+gender, family="binomial"))
  imp_tog_wob <- pool(imputed_together2)
  qbar_wob <- imp_tog_wob$qbar
  t_wob <- diag(imp_tog_wob$t)
  vmstar_wob <- imp_tog_wob$df
  
  return(list(qbar, t, vmstar, qbar_wob, t_wob, vmstar_wob))
}

pmi <- function(dataset, l, l2){
  total_imp <- 5
  
  #Organize Data
  b1_obs <- subset(dataset, !is.na(logpoll_cens) & batch1 == 1)
  b1_cens <- subset(dataset, is.na(logpoll_cens) & batch1 == 1)
  b2_obs <- subset(dataset, !is.na(logpoll_cens) & batch1 == 0)
  b2_cens <- subset(dataset, is.na(logpoll_cens) & batch1 == 0)
  
  miss_b1 <- nrow(b1_cens)
  miss_b2 <- nrow(b2_cens)
  
  imp_b1 <- list(b1_cens, b1_cens, b1_cens, b1_cens, b1_cens)
  imp_b2 <- list(b2_cens, b2_cens, b2_cens, b2_cens, b2_cens)
  
  for(k1 in 1:total_imp){
    #Get Boostratpped Sample
    dataset_boot <- dataset[sample(1:nrow(dataset), nrow(dataset), replace = TRUE),]
    
    b1_obs_boot <- subset(dataset_boot, !is.na(logpoll_cens) & batch1 == 1)
    b1_cens_boot <- subset(dataset_boot, is.na(logpoll_cens) & batch1 == 1)
    b2_obs_boot <- subset(dataset_boot, !is.na(logpoll_cens) & batch1 == 0)
    b2_cens_boot <- subset(dataset_boot, is.na(logpoll_cens) & batch1 == 0)
    
    #Get MLE for Bootstrapped Sample
    dist_x_given_yc <- function(theta){
      -(sum(log(pnorm(l, mean = theta[1]+theta[2]*b1_cens_boot$smoking+theta[3]*b1_cens_boot$gender+theta[4]*b1_cens_boot$case_cntrl, sd = sqrt(theta[5])))) +
          sum(log(pnorm(l2, mean = theta[1]+theta[2]*b2_cens_boot$smoking+theta[3]*b2_cens_boot$gender+theta[4]*b2_cens_boot$case_cntrl, sd = sqrt(theta[5])))) -
          (1/(2*theta[5]))*sum((b1_obs_boot$logpoll_cens-(theta[1]+theta[2]*b1_obs_boot$smoking+theta[3]*b1_obs_boot$gender+theta[4]*b1_obs_boot$case_cntrl))^2) -
          (1/(2*theta[5]))*sum((b2_obs_boot$logpoll_cens-(theta[1]+theta[2]*b2_obs_boot$smoking+theta[3]*b2_obs_boot$gender+theta[4]*b2_obs_boot$case_cntrl))^2) -
          (nrow(b1_obs_boot)+nrow(b2_obs_boot))*(1/2)*log(2*pi*theta[5]))
    }
    
    mle <- optim(c(0,0,0,0,1), dist_x_given_yc, hessian=TRUE)
    mle_param <- mle$par
    #fisher_inf <- -solve(-mle$hessian)
    #prop_sigma <- sqrt(diag(fisher_inf))
    
    #Impute missing values in Batch 1
    for(k2 in 1:miss_b1){
      p <- pnorm(l, mean = mle_param[1]+mle_param[2]*b1_cens$smoking[k2]+mle_param[3]*b1_cens$gender[k2]+mle_param[4]*b1_cens$case_cntrl[k2], sd = sqrt(mle_param[5]))
      unif_draw <- runif(1,min=0,max=p)
      draw <- qnorm(unif_draw, mean = mle_param[1]+mle_param[2]*b1_cens$smoking[k2]+mle_param[3]*b1_cens$gender[k2]+mle_param[4]*b1_cens$case_cntrl[k2], sd = sqrt(mle_param[5]))
      imp_b1[[k1]]$logpoll_cens[k2] <- draw
    }
    
    #Impute missing values in Batch 2
    for(k2 in 1:miss_b2){
      p <- pnorm(l2, mean = mle_param[1]+mle_param[2]*b2_cens$smoking[k2]+mle_param[3]*b2_cens$gender[k2]+mle_param[4]*b2_cens$case_cntrl[k2], sd = sqrt(mle_param[5]))
      unif_draw <- runif(1,min=0,max=p)
      draw <- qnorm(unif_draw, mean = mle_param[1]+mle_param[2]*b2_cens$smoking[k2]+mle_param[3]*b2_cens$gender[k2]+mle_param[4]*b2_cens$case_cntrl[k2], sd = sqrt(mle_param[5]))
      imp_b2[[k1]]$logpoll_cens[k2] <- draw
    }
  }
  
  #Check MLE Estimates
  #Get MLE for Bootstrapped Sample
  dist_x_given_yc <- function(theta){
    -(sum(log(pnorm(l, mean = theta[1]+theta[2]*b1_cens$smoking+theta[3]*b1_cens$gender+theta[4]*b1_cens$case_cntrl, sd = sqrt(theta[5])))) +
        sum(log(pnorm(l2, mean = theta[1]+theta[2]*b2_cens$smoking+theta[3]*b2_cens$gender+theta[4]*b2_cens$case_cntrl, sd = sqrt(theta[5])))) -
        (1/(2*theta[5]))*sum((b1_obs$logpoll_cens-(theta[1]+theta[2]*b1_obs$smoking+theta[3]*b1_obs$gender+theta[4]*b1_obs$case_cntrl))^2) -
        (1/(2*theta[5]))*sum((b2_obs$logpoll_cens-(theta[1]+theta[2]*b2_obs$smoking+theta[3]*b2_obs$gender+theta[4]*b2_obs$case_cntrl))^2) -
        (nrow(b1_obs)+nrow(b2_obs))*(1/2)*log(2*pi*theta[5]))
  }
  
  mle <- optim(c(0,0,0,0,1), dist_x_given_yc, hessian=TRUE)
  mle_param <- mle$par
  fisher_inf <- -solve(-mle$hessian)
  prop_sigma <- sqrt(diag(fisher_inf))
  
  #Aggregate imputed dataset with observed datasets
  poll_imp <- imp_b1
  for(k1 in 1:length(poll_imp)){
    poll_imp[[k1]] <- rbind(imp_b1[[k1]], imp_b2[[k1]], b1_obs, b2_obs)
  }
  
  #Get estimates and standard errors for each imputed dataset
  betas_with <- matrix(rep(0,5*total_imp), nrow = total_imp, ncol = 5)
  se_with <- matrix(rep(0,5*total_imp), nrow = total_imp, ncol = 5)
  betas_without <- matrix(rep(0,4*total_imp), nrow = total_imp, ncol = 4)
  se_without <- matrix(rep(0,4*total_imp), nrow = total_imp, ncol = 4)
  for(k1 in 1:length(poll_imp)){
    withbatchsep <- glm(data=poll_imp[[k1]], case_cntrl~logpoll_cens+smoking+gender+batch1, family="binomial")
    betas_with[k1,] <- coefficients(withbatchsep)
    se_with[k1,] <- sqrt((diag(summary(withbatchsep)$cov.unscaled)*summary(withbatchsep)$dispersion))
    withoutbatchsep <- glm(data=poll_imp[[k1]], case_cntrl~logpoll_cens+smoking+gender, family="binomial")
    betas_without[k1,] <- coefficients(withoutbatchsep)
    se_without[k1,] <- sqrt((diag(summary(withoutbatchsep)$cov.unscaled)*summary(withoutbatchsep)$dispersion))
  }
  
  #Get estimates using Rubin's combination rules
  qbar_with <- apply(betas_with, 2, mean)
  qbar_without <- apply(betas_without, 2, mean)
  ubar_with <- apply(se_with^2, 2, mean)
  ubar_without <- apply(se_without^2, 2, mean)
  b_with <- apply(betas_with, 2, var)
  b_without <- apply(betas_without, 2, var)
  t_with <- ubar_with + (1+(1/total_imp))*b_with
  t_without <- ubar_without + (1+(1/total_imp))*b_without
  r_with <- ((1+(1/total_imp))*b_with)/ubar_with
  r_without <- ((1+(1/total_imp))*b_without)/ubar_without
  vm_with <- (total_imp-1)*(1+(1/r_with))^2
  vm_without <- (total_imp-1)*(1+(1/r_without))^2
  gamma_with <- ((1+(1/total_imp))*b_with)/t_with
  gamma_without <- ((1+(1/total_imp))*b_without)/t_without
  v0_with <- nrow(dataset)-5
  v0_without <- nrow(dataset)-4
  vmstar_with <- ((1/vm_with)+(1/(((1-gamma_with)*(v0_with*(v0_with+1)))/(v0_with+3))))^(-1)
  vmstar_without <- ((1/vm_without)+(1/(((1-gamma_without)*(v0_without*(v0_without+1)))/(v0_without+3))))^(-1)
  
  return(list(qbar_with, t_with, vm_with, vmstar_with,
              qbar_without, t_without, vm_without, vmstar_without,
              mle_param, prop_sigma))
}

sim.lod <- function(total_itr, lod_mat, cases, cntrls, design, large){
  #Store results
  with_store <- matrix(rep(0,9*total_itr*90), nrow = 9*total_itr, ncol = 90)
  without_store <- matrix(rep(0,9*total_itr*77), nrow = 9*total_itr, ncol = 77)
  
  #Main simulation
  i <- 0
  increm_seed <- 1736
  total_pop <- rbind(cases, cntrls)
  for(lod_combination in 1:nrow(lod_mat)){
    l <- lod_mat[lod_combination,1]
    l2 <- lod_mat[lod_combination,2]
    i <- i + 1
    print(i)
    seed_num <- 425
    
    for(j in 1:total_itr){
      print(j)
      up_seed <- seed_num + (j-1)*increm_seed
      set.seed(up_seed)
      
      #Draw sample corresponding to study design
      if(design == "Case Control - Random Batch" & large == TRUE){
        ncases <- 1000
        ncntrls <- 4000
        cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
        cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
        dataset <- rbind(cases_samp, cntrls_samp)
      } else if(design == "Case Control - Random Batch" & large == FALSE){
        ncases <- 500
        ncntrls <- 500
        cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
        cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
        dataset <- rbind(cases_samp, cntrls_samp)
      } else if(design == "Case Control - Dependent Batch" & large == TRUE){
        ncases <- 1000
        ncntrls <- 4000
        cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
        cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
        dataset <- rbind(cases_samp, cntrls_samp)
      } else if(design == "Case Control - Dependent Batch" & large == FALSE){
        ncases <- 500
        ncntrls <- 500
        cases_samp <- subset(cases, id %in% sample(cases$id, ncases, replace = FALSE))
        cntrls_samp <- subset(cntrls, id %in% sample(cntrls$id, ncntrls, replace = FALSE))
        dataset <- rbind(cases_samp, cntrls_samp)
      } else if(design == "Cohort - Random Batch" & large == TRUE){
        ncohort <- 5000
        dataset <- subset(total_pop, id %in% sample(total_pop$id, ncohort, replace = FALSE))
      } else if(design == "Cohort - Random Batch" & large == FALSE){
        ncohort <- 1000
        dataset <- subset(total_pop, id %in% sample(total_pop$id, ncohort, replace = FALSE))
      } else if(design == "Cohort - Dependent Batch" & large == TRUE){
        ncohort <- 5000
        dataset <- subset(total_pop, id %in% sample(total_pop$id, ncohort, replace = FALSE))
      } else if(design == "Cohort - Dependent Batch" & large == FALSE){
        ncohort <- 1000
        dataset <- subset(total_pop, id %in% sample(total_pop$id, ncohort, replace = FALSE))
      }
      
      #Determine Batch Assignment
      if(design == "Case Control - Dependent Batch"){
        expalpha0 <- exp(-0.3) #Determined by setting P(Y=1) = 0.5 and solving for Beta0
        expalpha1 <- 3
        batch <- log(expalpha0) + log(expalpha1)*dataset$case_cntrl
        prob <- 1/(1+exp(-batch))
        batch1 <- rbinom(length(batch), 1, prob)
        dataset$batch1 <- batch1
      } else if(design == "Cohort - Dependent Batch"){
        expalpha0 <- exp(0.73) #Determined by setting P(Y=1) = 0.5 and solving for Beta0
        expalpha1 <- exp(1.5)
        expalpha2 <- exp(-2)
        batch <- log(expalpha0) + log(expalpha1)*dataset$smoking + log(expalpha2)*dataset$gender
        #True Probabilities with noise
        prob <- 1/(1+exp(-batch))
        batch1 <- rbinom(length(batch), 1, prob)
        #sum(batch1)
        dataset$batch1 <- batch1
      } else if(design == "Case Control - Random Batch" | design == "Cohort - Random Batch"){
        batch1 <- sample(c(0,1), nrow(dataset), replace = TRUE)
        dataset$batch1 <- batch1
      }
      
      #Censor values below LOD
      dataset$logpoll_cens <- dataset$logpoll
      dataset$logpoll_cens[dataset$batch1 == 1 & dataset$logpoll_cens < l] <- NA
      dataset$logpoll_cens[dataset$batch1 == 0 & dataset$logpoll_cens < l2] <- NA
      
      #Descriptive Statistics By Batch
      #Proportion censored in each batch
      per_bel_b1 <- sum(is.na(subset(dataset, batch1 == 1)$logpoll_cens))/nrow(subset(dataset, batch1 == 1))
      per_bel_b2 <- sum(is.na(subset(dataset, batch1 == 0)$logpoll_cens))/nrow(subset(dataset, batch1 == 0))
      #Descriptives
      dataset_cases <- sum(dataset$case_cntrl)
      dataset_smokers <- sum(dataset$smoking)
      dataset_gender <- sum(dataset$gender)
      dataset_batch <- sum(dataset$batch1)
      mean_logpoll <- mean(dataset$logpoll)
      sd_logpoll <- sd(dataset$logpoll)
      #Percent of cases within each batch
      per_cases_in_b1 <- sum(dataset$case_cntrl*dataset$batch1)/sum(dataset$batch1)
      per_cases_in_b2 <- sum(dataset$case_cntrl*(1-dataset$batch1))/sum(1-dataset$batch1)
      #Estimated probability of being in batch 1 given case
      prob_b1_case <- sum(dataset$case_cntrl*dataset$batch1)/sum(dataset$case_cntrl)
      #Estimated probability of being in batch 2 given case
      prob_b2_case <- sum(dataset$case_cntrl*(1-dataset$batch1))/sum(dataset$case_cntrl)
      
      #LOD Methods
      
      #Gold Standard
      gs <- gold_standard(dataset)
      
      #Complete Case Analysis
      cca <- complete_case_analysis(dataset)
      
      #LOD/SQRT(2)
      lodsqrt2 <- lod_sqrt2_substitution(dataset, exp(l), exp(l2)) #Need to exponentiate to get LOD on raw scale
      
      #Multiple Imputation - MICE
      mice_imp <- mice_imputation(dataset)
      
      #Censored Likliehood Multiple Imputation
      pmimp <- pmi(dataset, l, l2)
      
      #Store results
      with_store[(j+(i-1)*total_itr),] <- c(l, l2, up_seed, dataset_cases, dataset_smokers, dataset_gender,
                                            dataset_batch, mean_logpoll, sd_logpoll, per_bel_b1, per_bel_b2,
                                            per_cases_in_b1, per_cases_in_b2, prob_b1_case, prob_b2_case,
                                            gs[[1]], gs[[2]], cca[[1]], cca[[2]], lodsqrt2[[1]], lodsqrt2[[2]],
                                            mice_imp[[1]], mice_imp[[2]], mice_imp[[3]],
                                            pmimp[[1]], pmimp[[2]], pmimp[[3]], pmimp[[4]],
                                            pmimp[[9]], pmimp[[10]])
      
      without_store[(j+(i-1)*total_itr),] <- c(l, l2, up_seed, dataset_cases, dataset_smokers, dataset_gender,
                                            dataset_batch, mean_logpoll, sd_logpoll, per_bel_b1, per_bel_b2,
                                            per_cases_in_b1, per_cases_in_b2, prob_b1_case, prob_b2_case,
                                            gs[[3]], gs[[4]], cca[[3]], cca[[4]], lodsqrt2[[3]], lodsqrt2[[4]],
                                            mice_imp[[4]], mice_imp[[5]], mice_imp[[6]],
                                            pmimp[[5]], pmimp[[6]], pmimp[[7]], pmimp[[8]],
                                            pmimp[[9]], pmimp[[10]])
    }
  }
  return(list(with_store, without_store))
}

export_prep <- function(wb, wob){
  wb_df <- as.data.frame(wb)
  names(wb_df) <- c("LOD1","LOD2","Seed","NumCasesTotal","NumSmkTotal","NumGenTotal","NumBatch1Total","MeanLogPoll",
                    "SDLogPoll","PerBelowB1","PerBelowB2","PerCasesB1","PerCasesB2","ProbCaseB1","ProbCaseB2",
                    "GSBeta0","GSBeta1","GSBeta2","GSBeta3","GSBeta4","GSVar0","GSVar1","GSVar2","GSVar3","GSVar4",
                    "CCABeta0","CCABeta1","CCABeta2","CCABeta3","CCABeta4","CCAVar0","CCAVar1","CCAVar2","CCAVar3","CCAVar4",
                    "SQ2Beta0","SQ2Beta1","SQ2Beta2","SQ2Beta3","SQ2Beta4","SQ2Var0","SQ2Var1","SQ2Var2","SQ2Var3","SQ2Var4",
                    "MICEBeta0","MICEBeta1","MICEBeta2","MICEBeta3","MICEBeta4","MICEVar0","MICEVar1","MICEVar2","MICEVar3","MICEVar4",
                    "MICEVMStar0","MICEVMStar1","MICEVMStar2","MICEVMStar3","MICEVMStar4",
                    "PMIBeta0","PMIBeta1","PMIBeta2","PMIBeta3","PMIBeta4","PMIVar0","PMIVar1","PMIVar2","PMIVar3","PMIVar4",
                    "PMIVM0","PMIVM1","PMIVM2","PMIVM3","PMIVM4","PMIVMStar0","PMIVMStar1","PMIVMStar2","PMIVMStar3","PMIVMStar4",
                    "PMIMLE0","PMIMLE1","PMIMLE2","PMIMLE3","PMIMLE4","PMIMLEStdErr0","PMIMLEStdErr1","PMIMLEStdErr2","PMIMLEStdErr3","PMIMLEStdErr4")
  
  wob_df <- as.data.frame(wob)
  names(wob_df) <- c("LOD1","LOD2","Seed","NumCasesTotal","NumSmkTotal","NumGenTotal","NumBatch1Total","MeanLogPoll",
                     "SDLogPoll","PerBelowB1","PerBelowB2","PerCasesB1","PerCasesB2","ProbCaseB1","ProbCaseB2",
                     "GSBeta0","GSBeta1","GSBeta2","GSBeta3","GSVar0","GSVar1","GSVar2","GSVar3",
                     "CCABeta0","CCABeta1","CCABeta2","CCABeta3","CCAVar0","CCAVar1","CCAVar2","CCAVar3",
                     "SQ2Beta0","SQ2Beta1","SQ2Beta2","SQ2Beta3","SQ2Var0","SQ2Var1","SQ2Var2","SQ2Var3",
                     "MICEBeta0","MICEBeta1","MICEBeta2","MICEBeta3","MICEVar0","MICEVar1","MICEVar2","MICEVar3",
                     "MICEVMStar0","MICEVMStar1","MICEVMStar2","MICEVMStar3",
                     "PMIBeta0","PMIBeta1","PMIBeta2","PMIBeta3","PMIVar0","PMIVar1","PMIVar2","PMIVar3",
                     "PMIVM0","PMIVM1","PMIVM2","PMIVM3","PMIVMStar0","PMIVMStar1","PMIVMStar2","PMIVMStar3",
                     "PMIMLE0","PMIMLE1","PMIMLE2","PMIMLE3","PMIMLE4","PMIMLEStdErr0","PMIMLEStdErr1","PMIMLEStdErr2","PMIMLEStdErr3","PMIMLEStdErr4")
  
  return(list(wb_df, wob_df))
}

#############################
# Main Part of Program
#############################
pop <- gen_population_data()
cases <- pop[[1]]
cntrls <- pop[[2]]

#######################################################################################
# Large Cohort - Random Batch Assignment
#######################################################################################

set.seed(58681)

design <- "Cohort - Random Batch"
large <- TRUE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
        design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "largecohort_randombatch_withind.txt", sep=",")
write.table(save_res[[2]], "largecohort_randombatch_withoutind.txt", sep=",")

#######################################################################################
# Moderate Cohort - Random Batch Assignment
#######################################################################################

set.seed(24900)

design <- "Cohort - Random Batch"
large <- FALSE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "modcohort_randombatch_withind.txt", sep=",")
write.table(save_res[[2]], "modcohort_randombatch_withoutind.txt", sep=",")

#######################################################################################
# Large Cohort - Dependent Batch Assignment
#######################################################################################

set.seed(33195)

design <- "Cohort - Dependent Batch"
large <- TRUE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "largecohort_depbatch_withind.txt", sep=",")
write.table(save_res[[2]], "largecohort_depbatch_withoutind.txt", sep=",")

#######################################################################################
# Moderate Cohort - Dependent Batch Assignment
#######################################################################################

set.seed(10183)

design <- "Cohort - Dependent Batch"
large <- FALSE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "modcohort_depbatch_withind.txt", sep=",")
write.table(save_res[[2]], "modcohort_depbatch_withoutind.txt", sep=",")

#######################################################################################
# Large Case-Control - Random Batch Assignment
#######################################################################################

set.seed(21217)

design <- "Case Control - Random Batch"
large <- TRUE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "largecc_randombatch_withind.txt", sep=",")
write.table(save_res[[2]], "largecc_randombatch_withoutind.txt", sep=",")

#######################################################################################
# Moderate Case-Control - Random Batch Assignment
#######################################################################################

set.seed(88771)

design <- "Case Control - Random Batch"
large <- FALSE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "modcc_randombatch_withind.txt", sep=",")
write.table(save_res[[2]], "modcc_randombatch_withoutind.txt", sep=",")

#######################################################################################
# Large Case-Control - Dependent Batch Assignment
#######################################################################################

set.seed(52298)

design <- "Case Control - Dependent Batch"
large <- TRUE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "largecc_depbatch_withind.txt", sep=",")
write.table(save_res[[2]], "largecc_depbatch_withoutind.txt", sep=",")

#######################################################################################
# Moderate Case-Control - Dependent Batch Assignment
#######################################################################################

set.seed(78306)

design <- "Case Control - Dependent Batch"
large <- FALSE
total_itr <- 1000
get_lods <- gen_lods(cases = cases, cntrls = cntrls, design = design, large = large)
sim_res <- sim.lod(total_itr = total_itr, lod_mat = get_lods, cases = cases, cntrls = cntrls,
                   design = design, large = large)
save_res <- export_prep(sim_res[[1]], sim_res[[2]])
write.table(save_res[[1]], "modcc_depbatch_withind.txt", sep=",")
write.table(save_res[[2]], "modcc_depbatch_withoutind.txt", sep=",")


