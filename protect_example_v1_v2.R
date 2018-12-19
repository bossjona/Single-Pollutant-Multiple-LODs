library(mice)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

#####################################################################
## Functions
#####################################################################

##############################################################################################
## Important: This is an older version of CLMI that was used for this data example. For the
## most recent version visit https://github.com/bossjona/Single-Pollutant-Multiple-LODs
## and download the clmi.R file.
##############################################################################################

clmi <- function(data, cont_name, batch_name, outcome_name, imp_vars, lod_info, total_imp = 5, seed_num, transform_imp){
  #data: Data frame with contaminant concentration, batch number, covariates used in imputation, precision variables
  #(make sure categorical variables are numeric, not factors). Can use model.matrix function to convert data frame with
  #factors to numeric design matrix and convert that matrix back into a data frame.
  #cont_name: string corresponding to name of log-transformed contaminant variable.
  #batch_name: string corresponding to name of batch indicator variable.
  #outcome_name: string corresponding to name of binary health outcome.
  #imp_vars: vector of strings corresponding to covariates associated with the contaminant (covariates used in imputation).
  #lod_info: data frame with first column corresponding to each batch (entered as a string) and the
  #second column listing the LODs (entered numerically) corresponding to each batch.
  #First column must be named "batch_info" and the second column must be named "lod".
  #total_imp: number of multiply imputed datasets (integer).
  #seed_num: seed number to ensure reproducability.
  #transform_imp: transformation applied to contaminant data to make it normally distributed, default is no transformation.
  
  set.seed(seed_num)
  
  ################################################################################
  ## Organize Data into non-detects and detects
  ################################################################################
  
  data$lod <- rep(0, nrow(data))
  for(i in 1:nrow(data)){
    data$lod[i] <- lod_info[which(lod_info$batch_info == data[[batch_name]][i]),"lod"]
  }
  
  #save_original_data <- data
  
  cont_imp_name <- paste0(cont_name,"_imputed")
  #data[[cont_imp_name]] <- transform_imp(data[[cont_name]])
  data[[cont_imp_name]] <- data[[cont_name]]
  
  above_lod <- data[which(!is.na(data[[cont_imp_name]])),]
  below_lod <- data[which(is.na(data[[cont_imp_name]])),]
  
  if(!is.null(imp_vars)){
    above_matrix <- matrix(rep(0, nrow(above_lod)*length(imp_vars)), nrow = nrow(above_lod), ncol = length(imp_vars))
    below_matrix <- matrix(rep(0, nrow(below_lod)*length(imp_vars)), nrow = nrow(below_lod), ncol = length(imp_vars))
    for(l in 1:length(imp_vars)){
      above_matrix[,l] <- above_lod[[imp_vars[l]]]
      below_matrix[,l] <- below_lod[[imp_vars[l]]]
    }
  }
  
  ################################################################################
  ## Perform Multiple Imputation
  ################################################################################
  
  imp <- rep(list(below_lod), total_imp)
  
  for(j in 1:total_imp){
    #Bootstrap dataset
    data_boot <- data[sample(1:nrow(data), nrow(data), replace = TRUE),]
    
    #Subjects with concentration above LOD
    above_lod_boot <- data_boot[which(!is.na(data_boot[[cont_imp_name]])),]
    
    #Subjects with concentration below LOD for each batch
    below_lod_boot <- data_boot[which(is.na(data_boot[[cont_imp_name]])),]
    
    #Get MLE for Bootstrapped Sample
    if(!is.null(imp_vars)){
      #Organize covariates used in imputation
      above_matrix_boot <- matrix(rep(0, nrow(above_lod_boot)*length(imp_vars)), nrow = nrow(above_lod_boot), ncol = length(imp_vars))
      below_matrix_boot <- matrix(rep(0, nrow(below_lod_boot)*length(imp_vars)), nrow = nrow(below_lod_boot), ncol = length(imp_vars))
      for(k in 1:length(imp_vars)){
        above_matrix_boot[,k] <- above_lod_boot[[imp_vars[k]]]
        below_matrix_boot[,k] <- below_lod_boot[[imp_vars[k]]]
      }
      
      #Calculate MLE
      dist_x_given_yc <- function(theta){
        -(sum(log(pnorm(transform_imp(below_lod_boot[["lod"]]), mean = theta[1]+theta[2]*below_lod_boot[[outcome_name]]+below_matrix_boot%*%theta[4:(3+length(imp_vars))], sd = sqrt(theta[3])))) -
            (1/(2*theta[3]))*sum((transform_imp(above_lod_boot[[cont_imp_name]])-(theta[1]+theta[2]*above_lod_boot[[outcome_name]]+above_matrix_boot%*%theta[4:(3+length(imp_vars))]))^2) -
            nrow(above_lod_boot)*(1/2)*log(2*pi*theta[3]))
      }
      
      mle <- optim(c(0,0,1,rep(0,length(imp_vars))), dist_x_given_yc, hessian=FALSE)
      mle_param <- mle$par
      
      #Impute missing values
      for(k1 in 1:nrow(below_lod)){
        p <- pnorm(transform_imp(below_lod[["lod"]])[k1], mean = mle_param[1]+mle_param[2]*below_lod[[outcome_name]][k1]+below_matrix[k1,]%*%mle_param[4:length(mle_param)], sd = sqrt(mle_param[3]))
        unif_draw <- runif(1,min=0,max=p)
        draw <- qnorm(unif_draw, mean = mle_param[1]+mle_param[2]*below_lod[[outcome_name]][k1]+below_matrix[k1,]%*%mle_param[4:length(mle_param)], sd = sqrt(mle_param[3]))
        imp[[j]][[cont_imp_name]][k1] <- draw
      }
    } else if(is.null(imp_vars)){
      #Calculate MLE
      dist_x_given_yc <- function(theta){
        -(sum(log(pnorm(transform_imp(below_lod_boot[["lod"]]), mean = theta[1]+theta[2]*below_lod_boot[[outcome_name]], sd = sqrt(theta[3])))) -
            (1/(2*theta[3]))*sum((transform_imp(above_lod_boot[[cont_imp_name]])-(theta[1]+theta[2]*above_lod_boot[[outcome_name]]))^2) -
            nrow(above_lod_boot)*(1/2)*log(2*pi*theta[3]))
      }
      
      mle <- optim(c(0,0,1), dist_x_given_yc, hessian=FALSE)
      mle_param <- mle$par
      
      #Impute missing values
      for(k1 in 1:nrow(below_lod)){
        p <- pnorm(transform_imp(below_lod[["lod"]])[k1], mean = mle_param[1]+mle_param[2]*below_lod[[outcome_name]][k1], sd = sqrt(mle_param[3]))
        unif_draw <- runif(1,min=0,max=p)
        draw <- qnorm(unif_draw, mean = mle_param[1]+mle_param[2]*below_lod[[outcome_name]][k1], sd = sqrt(mle_param[3]))
        imp[[j]][[cont_imp_name]][k1] <- draw
      }
    }
  }
  
  if(!is.null(imp_vars)){
    #Check MLE estimation for original dataset
    dist_x_given_yc <- function(theta){
      -(sum(log(pnorm(transform_imp(below_lod[["lod"]]), mean = theta[1]+theta[2]*below_lod[[outcome_name]]+below_matrix%*%theta[4:(3+length(imp_vars))], sd = sqrt(theta[3])))) -
          (1/(2*theta[3]))*sum((transform_imp(above_lod[[cont_imp_name]])-(theta[1]+theta[2]*above_lod[[outcome_name]]+above_matrix%*%theta[4:(3+length(imp_vars))]))^2) -
          nrow(above_lod)*(1/2)*log(2*pi*theta[3]))
    }
    
    mle <- optim(c(0,0,1,rep(0,length(imp_vars))), dist_x_given_yc, hessian=TRUE)
    mle_param <- mle$par
    fisher_inf <- -solve(-mle$hessian)
    prop_sigma <- sqrt(diag(fisher_inf))
  } else if(is.null(imp_vars)){
    #Check MLE estimation for original dataset
    dist_x_given_yc <- function(theta){
      -(sum(log(pnorm(transform_imp(below_lod[["lod"]]), mean = theta[1]+theta[2]*below_lod[[outcome_name]], sd = sqrt(theta[3])))) -
          (1/(2*theta[3]))*sum((transform_imp(above_lod[[cont_imp_name]])-(theta[1]+theta[2]*above_lod[[outcome_name]]))^2) -
          nrow(above_lod)*(1/2)*log(2*pi*theta[3]))
    }
    
    mle <- optim(c(0,0,1), dist_x_given_yc, hessian=FALSE)
    mle_param <- mle$par
    fisher_inf <- -solve(-mle$hessian)
    prop_sigma <- sqrt(diag(fisher_inf))
  }
  
  #Transform Back to Original Scale
  inverse <- function(f, lower = -1000, upper = 1000){
    function(y){
      uniroot((function(x) f(x) - y), lower = lower, upper = upper)[1]
    }
  }
  
  inverse_transform <- inverse(transform_imp, 0, max(lod_info[,"lod"]))
  
  #Aggregate imputed datasets with observed dataset
  poll_imp <- imp
  for(k1 in 1:length(poll_imp)){
    for(conv in 1:nrow(imp[[k1]])){
      imp[[k1]][[cont_imp_name]][conv] <- inverse_transform(imp[[k1]][[cont_imp_name]][conv])$root
    }
    poll_imp[[k1]] <- rbind(imp[[k1]], above_lod)
  }
  return(list("nimp" = total_imp, "imputations" = poll_imp, "mle" = mle_param, "fisher_inf" = fisher_inf,
              var_names = list("cont_name" = cont_name, "batch_name" = batch_name, "outcome_name" = outcome_name,
                               "imp_vars" = imp_vars, "lod_info" = lod_info)))
}

pool.clmi <- function(clmi_obj, transform_model, regression_model, precis_vars){
  #clmi_obj: List generated from clmi() function is used as input for this argument
  #transform_imp: transformaed contaminant used in the regression model; default is no transformation (must be a 1-1 transformation)
  #transform_model: Linear Regression or Logistic Regression
  #precision variables: vector of strings corresponding to covariates associated with the contaminant
  #Get names used in clmi() function
  poll_imp <- clmi_obj$imputations
  batch_name <- clmi_obj$var_names$batch_name
  cont_imp_name <- paste0(clmi_obj$var_names$cont_name,"_imputed")
  outcome_name <- clmi_obj$var_names$outcome_name
  imp_vars <- clmi_obj$var_names$imp_vars
  total_imp <- clmi_obj$nimp
  nobs <- nrow(poll_imp[[1]])
  
  #Transform contaminant
  for(idx in 1:total_imp){
    poll_imp[[idx]][[cont_imp_name]] <- transform_model(poll_imp[[idx]][[cont_imp_name]])
  }
  
  #Get number of coefficients to store
  data.fit <- poll_imp[[1]][, names(poll_imp[[1]]) %in% imp_vars | names(poll_imp[[1]]) %in% precis_vars |
                              names(poll_imp[[1]]) %in% cont_imp_name]
  if(is.vector(data.fit)){
    num_continuous <- 1
    num_levels <- 0
  } else if(!is.vector(data.fit)){
    #Number of continuous variables
    # num_continuous <- sum(!sapply(data.fit, is.factor))
    num_continuous <- sum(sapply(data.fit, is.numeric))
    #Number of levels in factor variables
    # if(sum(sapply(data.fit, is.factor)) == 0){
    #   num_levels <- 0
    # } else if(sum(sapply(data.fit, is.factor)) != 0){
    #   num_levels <- sum(sapply(data.fit[,sapply(data.fit, is.factor)], nlevels)) - sum(sapply(data.fit, is.factor))
    # }
  }
  #Total number of parameters in logistic regression (add one for the intercept)
  # total_param <- 1 + num_continuous + num_levels
  total_param <- 1 + num_continuous
  
  #Get estimates and standard errors for each imputed dataset
  betas <- matrix(rep(0,total_param*total_imp), nrow = total_imp, ncol = total_param)
  varcov <- rep(list(0), total_imp)
  regressions <- rep(list(0), total_imp)
  if(regression_model == "Logistic"){
    for(k1 in 1:length(poll_imp)){
      data.fit <- poll_imp[[k1]][, names(poll_imp[[k1]]) %in% imp_vars | names(poll_imp[[k1]]) %in% precis_vars |
                                   names(poll_imp[[k1]]) %in% cont_imp_name | names(poll_imp[[k1]]) %in% outcome_name]
      names(data.fit)[which(names(data.fit) == outcome_name)] <- "outcome"
      logreg <- glm(factor(outcome)~., family="binomial", data = data.fit)
      betas[k1,] <- coefficients(logreg)
      varcov[[k1]] <- summary(logreg)$cov.unscaled
      regressions[[k1]] <- logreg
    }
  } else if(regression_model == "Linear"){
    for(k1 in 1:length(poll_imp)){
      data.fit <- poll_imp[[k1]][, names(poll_imp[[k1]]) %in% imp_vars | names(poll_imp[[k1]]) %in% precis_vars |
                                   names(poll_imp[[k1]]) %in% cont_imp_name | names(poll_imp[[k1]]) %in% outcome_name]
      names(data.fit)[which(names(data.fit) == outcome_name)] <- "outcome"
      linreg <- lm(outcome~., data = data.fit)
      betas[k1,] <- coefficients(linreg)
      varcov[[k1]] <- vcov(linreg)
      regressions[[k1]] <- linreg
    }
  }
  
  #Pooled Inference for Multiply Imputed Datasets
  qbar <- apply(betas, 2, mean)
  ubar <- Reduce('+', varcov)/total_imp
  q_diff <- apply(betas, 1, function(x){x - qbar}) #Each column of this matrix is Qi - Qbar
  ui <- matrix(rep(0, ncol(betas)*ncol(betas)), nrow = ncol(betas), ncol = ncol(betas))
  for(m in 1:ncol(q_diff)){
    ui <- ui + as.matrix(q_diff[,m])%*%t(as.matrix(q_diff[,m]))
  }
  b <- ui/(total_imp-1)
  t <- ubar + (1+(1/total_imp))*b
  gamma <- ((1+(1/total_imp))*diag(b))/diag(t)
  r <- ((1+(1/total_imp))*diag(b))/diag(ubar)
  v <- (total_imp-1)*(1+(1/r))^2
  v0 <- nobs-ncol(betas)
  denom_adj_df <- ((1-gamma)*v0*(v0+1))/(v0+3)
  vstar <- ((1/v)+(1/denom_adj_df))^(-1)
  
  #Get Confidence Intervals
  LCL <- qbar - qt(0.975, vstar)*sqrt(diag(t))
  UCL <- qbar + qt(0.975, vstar)*sqrt(diag(t))
  
  #Get p-values
  p_vals <- 2*pt(abs(qbar/sqrt(diag(t))), df = vstar, lower.tail = FALSE)
  
  #Summary of each regression models on each imputed datasets
  individual_summary <- summary(regressions[[1]])$coefficients
  for(k1 in 2:length(poll_imp)){
    individual_summary <- rbind(individual_summary, summary(regressions[[k1]])$coefficients)
  }
  
  return(list("output" = data.frame("Est" = qbar, "SE" = sqrt(diag(t)), "df" = vstar, "P.Value" = p_vals, "LCL.95" = LCL, "UCL.95" = UCL), "pooled_vcov" = t,
              "reg_models" = regressions, "reg_models_summary" = individual_summary))
}

#Plot paneled graphs with common legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, c(lapply(plots, function(x)
      x + theme(legend.position="none")), nrow = nrow, ncol = ncol)),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

#Read in data
setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18\\PROTECT Example")

visit1 <- read.csv("visit1.csv", stringsAsFactors = FALSE)
visit2 <- read.csv("visit2.csv", stringsAsFactors = FALSE)

#Define covariates
visit1$currjob_cat[visit1$currjob_cat == 2] <- 0
visit2$currjob_cat[visit2$currjob_cat == 2] <- 0

visit1$childnum1 <- rep(0, nrow(visit1))
visit1$childnum1[visit1$childnum == "1"] <- 1
visit1$childnum2 <- rep(0, nrow(visit1))
visit1$childnum2[visit1$childnum == ">= 2"] <- 1

visit2$childnum1 <- rep(0, nrow(visit2))
visit2$childnum1[visit2$childnum == "1"] <- 1
visit2$childnum2 <- rep(0, nrow(visit2))
visit2$childnum2[visit2$childnum == ">= 2"] <- 1

#Convert batch to character
visit1$mehp_batch <- as.character(visit1$mehp_batch)
visit1$mcpp_batch <- as.character(visit1$mcpp_batch)
visit1$bpb_batch <- as.character(visit1$bpb_batch)
visit1$bpf_batch <- as.character(visit1$bpf_batch)
visit1$tcs_batch <- as.character(visit1$tcs_batch)
visit1$tcc_batch <- as.character(visit1$tcc_batch)

visit2$mehp_batch <- as.character(visit2$mehp_batch)
visit2$mcpp_batch <- as.character(visit2$mcpp_batch)
visit2$bpb_batch <- as.character(visit2$bpb_batch)
visit2$bpf_batch <- as.character(visit2$bpf_batch)
visit2$tcs_batch <- as.character(visit2$tcs_batch)
visit2$tcc_batch <- as.character(visit2$tcc_batch)

#Smaller Dataset for BPF and TCC
visit1_small <- subset(visit1, (!is.na(bpf) | (is.na(bpf) & bpf_lod != "")) & (!is.na(tcc) | (is.na(tcc) & tcc_lod != "")))
visit2_small <- subset(visit2, (!is.na(bpf) | (is.na(bpf) & bpf_lod != "")) & (!is.na(tcc) | (is.na(tcc) & tcc_lod != "")))

####################################################################################
## Visit 1
####################################################################################

#######################
## Summary Statistics
#######################

visit1 %>% group_by(spont_preterm, childnum) %>% tally() %>% mutate("Per" = 100*n/sum(n))
visit1 %>% group_by(spont_preterm, currjob_cat) %>% tally() %>% mutate("Per" = 100*n/sum(n))

visit1 %>% group_by(spont_preterm) %>% summarise("M" = mean(isage), "SD" = sd(isage))
visit1 %>% group_by(spont_preterm) %>% summarise("M" = mean(sg_summary), "SD" = sd(sg_summary))

###############
## MEHP
###############

#CLMI
data <- visit1
cont_name <- "mehp"
batch_name <- "mehp_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.5, 0.8))
total_imp <- 10
seed_num <- 667754
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$mehp_batch_model <- as.numeric(data$mehp_batch)
data$mehp_batch_model[data$mehp_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
     lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
          precis_vars = c("isage","currjob_cat","childnum1","childnum2","mehp_batch_model"))
res$output["mehp_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch),
    family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(mehp))"], sqrt(vcov(cca_res)["I(log(mehp))","I(log(mehp))"]), confint(cca_res)["I(log(mehp))",])

#MICE
mice_data <- visit1[, names(visit1) %in%
                      c("spont_preterm","mehp_batch","mehp","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$mehp_batch <- as.factor(mice_data$mehp_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(mehp))","est"], digits = 3), round(summary(imp)["I(log(mehp))","se"], digits = 3),
          round(summary(imp)["I(log(mehp))","lo 95"], digits = 3), round(summary(imp)["I(log(mehp))","hi 95"], digits = 3))

#SQRT 2
data$mehp[is.na(data$mehp) & data$mehp_lod == "<LOD(0.5)"] <- 0.5/sqrt(2)
data$mehp[is.na(data$mehp) & data$mehp_lod == "<LOD(0.8)"] <- 0.8/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(mehp))"], sqrt(vcov(sq2_res)["I(log(mehp))","I(log(mehp))"]), confint(sq2_res)["I(log(mehp))",])

#######################
## Plot
#######################

visit1_mehp <- ggplot(data = data, aes(x = log(mehp), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed MEHP (ng/ml)") + ggtitle("Visit 1") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.5/sqrt(2))) +
  geom_vline(xintercept = log(0.8/sqrt(2)))

###############
## MCPP
###############

#CLMI
data <- visit1
cont_name <- "mcpp"
batch_name <- "mcpp_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.2, 0.4))
total_imp <- 10
seed_num <- 943201
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$mcpp_batch_model <- as.numeric(data$mcpp_batch)
data$mcpp_batch_model[data$mcpp_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","mcpp_batch_model"))
res$output["mcpp_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(mcpp))"], sqrt(vcov(cca_res)["I(log(mcpp))","I(log(mcpp))"]), confint(cca_res)["I(log(mcpp))",])

#MICE
mice_data <- visit1[, names(visit1) %in%
                      c("spont_preterm","mcpp_batch","mcpp","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$mcpp_batch <- as.factor(mice_data$mcpp_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(mcpp))","est"], digits = 3), round(summary(imp)["I(log(mcpp))","se"], digits = 3),
  round(summary(imp)["I(log(mcpp))","lo 95"], digits = 3), round(summary(imp)["I(log(mcpp))","hi 95"], digits = 3))

#SQRT 2
data$mcpp[is.na(data$mcpp) & data$mcpp_lod == "<LOD(0.2)"] <- 0.2/sqrt(2)
data$mcpp[is.na(data$mcpp) & data$mcpp_lod == "<LOD(0.4)"] <- 0.4/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(mcpp))"], sqrt(vcov(sq2_res)["I(log(mcpp))","I(log(mcpp))"]), confint(sq2_res)["I(log(mcpp))",])

#######################
## Plot
#######################

visit1_mcpp <- ggplot(data = data, aes(x = log(mcpp), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed MCPP (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.2/sqrt(2))) +
  geom_vline(xintercept = log(0.4/sqrt(2)))

###############
## BPB
###############

#CLMI
data <- visit1
cont_name <- "bpb"
batch_name <- "bpb_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2","3","4"), "lod" = c(0.15, 0.2, 0.1, 0.4))
total_imp <- 10
seed_num <- 50032
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$bpb_batch_model2 <- as.numeric(data$bpb_batch == 2)
data$bpb_batch_model34 <- as.numeric(data$bpb_batch %in% c(3,4))

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","bpb_batch_model2","bpb_batch_model34"))
res$output["bpb_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_data <- data
cca_data$bpb_batch[cca_data$bpb_batch == 4] <- "3"
cca_res <- glm(data = cca_data, factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(bpb))"], sqrt(vcov(cca_res)["I(log(bpb))","I(log(bpb))"]), confint(cca_res)["I(log(bpb))",])

#MICE
mice_data <- visit1[, names(visit1) %in%
                      c("spont_preterm","bpb_batch","bpb","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$bpb_batch[mice_data$bpb_batch == 4] <- "3"
mice_data$bpb_batch <- as.factor(mice_data$bpb_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(bpb))","est"], digits = 3), round(summary(imp)["I(log(bpb))","se"], digits = 3),
  round(summary(imp)["I(log(bpb))","lo 95"], digits = 3), round(summary(imp)["I(log(bpb))","hi 95"], digits = 3))

#SQRT 2
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.1"] <- 0.1/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.15"] <- 0.15/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.2"] <- 0.2/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.4"] <- 0.4/sqrt(2)
data$bpb_batch[data$bpb_batch == 4] <- "3"
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(bpb))"], sqrt(vcov(sq2_res)["I(log(bpb))","I(log(bpb))"]), confint(sq2_res)["I(log(bpb))",])

#######################
## Plot
#######################

visit1_bpb <- ggplot(data = data, aes(x = log(bpb), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed BPB (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(0.15/sqrt(2))) + geom_vline(xintercept = log(0.2/sqrt(2))) +
  geom_vline(xintercept = log(0.4/sqrt(2)))

###############
## BPF
###############

#CLMI
data <- visit1_small
cont_name <- "bpf"
batch_name <- "bpf_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.2, 0.1))
total_imp <- 10
seed_num <- 84332
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$bpf_batch_model <- as.numeric(data$bpf_batch)
data$bpf_batch_model[data$bpf_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2"))
res$output["bpf_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat,
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(bpf))"], sqrt(vcov(cca_res)["I(log(bpf))","I(log(bpf))"]), confint(cca_res)["I(log(bpf))",])

#MICE
mice_data <- visit1_small[, names(visit1_small) %in%
                      c("spont_preterm","bpf","sg","isage","childnum1","childnum2","currjob_cat")]
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat, family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(bpf))","est"], digits = 3), round(summary(imp)["I(log(bpf))","se"], digits = 3),
  round(summary(imp)["I(log(bpf))","lo 95"], digits = 3), round(summary(imp)["I(log(bpf))","hi 95"], digits = 3))

#SQRT 2
data$bpf[is.na(data$bpf) & data$bpf_lod == "<0.1"] <- 0.1/sqrt(2)
data$bpf[is.na(data$bpf) & data$bpf_lod == "<0.2"] <- 0.2/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat,
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(bpf))"], sqrt(vcov(sq2_res)["I(log(bpf))","I(log(bpf))"]), confint(sq2_res)["I(log(bpf))",])

#######################
## Plot
#######################

visit1_bpf <- ggplot(data = data, aes(x = log(bpf), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed BPF (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(0.2/sqrt(2)))

###############
## TCS
###############

#CLMI
data <- visit1
cont_name <- "tcs"
batch_name <- "tcs_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(2.3, 1.7))
total_imp <- 10
seed_num <- 22560
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$tcs_batch_model <- as.numeric(data$tcs_batch)
data$tcs_batch_model[data$tcs_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","tcs_batch_model"))
res$output["tcs_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(tcs))"], sqrt(vcov(cca_res)["I(log(tcs))","I(log(tcs))"]), confint(cca_res)["I(log(tcs))",])

#MICE
mice_data <- visit1[, names(visit1) %in%
                      c("spont_preterm","tcs_batch","tcs","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$tcs_batch <- as.factor(mice_data$tcs_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(tcs))","est"], digits = 3), round(summary(imp)["I(log(tcs))","se"], digits = 3),
  round(summary(imp)["I(log(tcs))","lo 95"], digits = 3), round(summary(imp)["I(log(tcs))","hi 95"], digits = 3))

#SQRT 2
data$tcs[is.na(data$tcs) & data$tcs_lod == "<1.7"] <- 1.7/sqrt(2)
data$tcs[is.na(data$tcs) & data$tcs_lod == "<2.3"] <- 2.3/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(tcs))"], sqrt(vcov(sq2_res)["I(log(tcs))","I(log(tcs))"]), confint(sq2_res)["I(log(tcs))",])

#######################
## Plot
#######################

visit1_tcs <- ggplot(data = data, aes(x = log(tcs), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed TCS (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(1.7/sqrt(2))) +
  geom_vline(xintercept = log(2.3/sqrt(2)))

###############
## TCC
###############

#CLMI
data <- visit1_small
cont_name <- "tcc"
batch_name <- "tcc_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.1, 1.7))
total_imp <- 10
seed_num <- 11092
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$tcc_batch_model <- as.numeric(data$tcc_batch)
data$tcc_batch_model[data$tcc_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","tcc_batch_model"))
res$output["tcc_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(tcc))"], sqrt(vcov(cca_res)["I(log(tcc))","I(log(tcc))"]), confint(cca_res)["I(log(tcc))",])

#MICE
mice_data <- visit1_small[, names(visit1_small) %in%
                      c("spont_preterm","tcc_batch","tcc","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$tcc_batch <- as.factor(mice_data$tcc_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(tcc))","est"], digits = 3), round(summary(imp)["I(log(tcc))","se"], digits = 3),
  round(summary(imp)["I(log(tcc))","lo 95"], digits = 3), round(summary(imp)["I(log(tcc))","hi 95"], digits = 3))

#SQRT 2
data$tcc[is.na(data$tcc) & data$tcc_lod == "<1.7"] <- 1.7/sqrt(2)
data$tcc[is.na(data$tcc) & data$tcc_lod == "<0.1"] <- 0.1/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(tcc))"], sqrt(vcov(sq2_res)["I(log(tcc))","I(log(tcc))"]), confint(sq2_res)["I(log(tcc))",])

#######################
## Plot
#######################

visit1_tcc <- ggplot(data = data, aes(x = log(tcc), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed TCC (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(1.7/sqrt(2)))

####################################################################################
## Visit 2
####################################################################################

#######################
## Summary Statistics
#######################

visit2 %>% group_by(spont_preterm, childnum) %>% tally() %>% mutate("Per" = 100*n/sum(n))
visit2 %>% group_by(spont_preterm, currjob_cat) %>% tally() %>% mutate("Per" = 100*n/sum(n))

visit2 %>% group_by(spont_preterm) %>% summarise("M" = mean(isage), "SD" = sd(isage))
visit2 %>% group_by(spont_preterm) %>% summarise("M" = mean(sg_summary), "SD" = sd(sg_summary))

###############
## MEHP
###############

#CLMI
data <- visit2
cont_name <- "mehp"
batch_name <- "mehp_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.5, 0.8))
total_imp <- 10
seed_num <- 98743
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$mehp_batch_model <- as.numeric(data$mehp_batch)
data$mehp_batch_model[data$mehp_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","mehp_batch_model"))
res$output["mehp_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(mehp))"], sqrt(vcov(cca_res)["I(log(mehp))","I(log(mehp))"]), confint(cca_res)["I(log(mehp))",])

#MICE
mice_data <- visit2[, names(visit2) %in%
                      c("spont_preterm","mehp_batch","mehp","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$mehp_batch <- as.factor(mice_data$mehp_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(mehp))","est"], digits = 3), round(summary(imp)["I(log(mehp))","se"], digits = 3),
  round(summary(imp)["I(log(mehp))","lo 95"], digits = 3), round(summary(imp)["I(log(mehp))","hi 95"], digits = 3))

#SQRT 2
data$mehp[is.na(data$mehp) & data$mehp_lod == "<LOD(0.5)"] <- 0.5/sqrt(2)
data$mehp[is.na(data$mehp) & data$mehp_lod == "<LOD(0.8)"] <- 0.8/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(mehp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mehp_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(mehp))"], sqrt(vcov(sq2_res)["I(log(mehp))","I(log(mehp))"]), confint(sq2_res)["I(log(mehp))",])

#######################
## Plot
#######################

visit2_mehp <- ggplot(data = data, aes(x = log(mehp), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed MEHP (ng/ml)") + ggtitle("Visit 2") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.5/sqrt(2))) +
  geom_vline(xintercept = log(0.8/sqrt(2)))

###############
## MCPP
###############

#CLMI
data <- visit2
cont_name <- "mcpp"
batch_name <- "mcpp_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.2, 0.4))
total_imp <- 10
seed_num <- 889110
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$mcpp_batch_model <- as.numeric(data$mcpp_batch)
data$mcpp_batch_model[data$mcpp_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","mcpp_batch_model"))
res$output["mcpp_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(mcpp))"], sqrt(vcov(cca_res)["I(log(mcpp))","I(log(mcpp))"]), confint(cca_res)["I(log(mcpp))",])

#MICE
mice_data <- visit2[, names(visit2) %in%
                      c("spont_preterm","mcpp_batch","mcpp","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$mcpp_batch <- as.factor(mice_data$mcpp_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(mcpp))","est"], digits = 3), round(summary(imp)["I(log(mcpp))","se"], digits = 3),
  round(summary(imp)["I(log(mcpp))","lo 95"], digits = 3), round(summary(imp)["I(log(mcpp))","hi 95"], digits = 3))

#SQRT 2
data$mcpp[is.na(data$mcpp) & data$mcpp_lod == "<LOD(0.2)"] <- 0.2/sqrt(2)
data$mcpp[is.na(data$mcpp) & data$mcpp_lod == "<LOD(0.4)"] <- 0.4/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(mcpp))+sg+isage+childnum1+childnum2+currjob_cat+factor(mcpp_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(mcpp))"], sqrt(vcov(sq2_res)["I(log(mcpp))","I(log(mcpp))"]), confint(sq2_res)["I(log(mcpp))",])

#######################
## Plot
#######################

visit2_mcpp <- ggplot(data = data, aes(x = log(mcpp), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed MCPP (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.2/sqrt(2))) +
  geom_vline(xintercept = log(0.4/sqrt(2)))

###############
## BPB
###############

#CLMI
data <- visit2
cont_name <- "bpb"
batch_name <- "bpb_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2","3","4"), "lod" = c(0.15, 0.2, 0.1, 0.4))
total_imp <- 10
seed_num <- 53127
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$bpb_batch_model2 <- as.numeric(data$bpb_batch == 2)
data$bpb_batch_model34 <- as.numeric(data$bpb_batch %in% c(3,4))

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","bpb_batch_model2","bpb_batch_model34"))
res$output["bpb_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_data <- data
cca_data$bpb_batch[cca_data$bpb_batch == 4] <- "3"
cca_res <- glm(data = cca_data, factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(bpb))"], sqrt(vcov(cca_res)["I(log(bpb))","I(log(bpb))"]), confint(cca_res)["I(log(bpb))",])

#MICE
mice_data <- visit2[, names(visit2) %in%
                      c("spont_preterm","bpb_batch","bpb","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$bpb_batch[mice_data$bpb_batch == 4] <- "3"
mice_data$bpb_batch <- as.factor(mice_data$bpb_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(bpb))","est"], digits = 3), round(summary(imp)["I(log(bpb))","se"], digits = 3),
  round(summary(imp)["I(log(bpb))","lo 95"], digits = 3), round(summary(imp)["I(log(bpb))","hi 95"], digits = 3))

#SQRT 2
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.1"] <- 0.1/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.15"] <- 0.15/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.2"] <- 0.2/sqrt(2)
data$bpb[is.na(data$bpb) & data$bpb_lod == "<0.4"] <- 0.4/sqrt(2)
data$bpb_batch[data$bpb_batch == 4] <- "3"
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(bpb))+sg+isage+childnum1+childnum2+currjob_cat+factor(bpb_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(bpb))"], sqrt(vcov(sq2_res)["I(log(bpb))","I(log(bpb))"]), confint(sq2_res)["I(log(bpb))",])

#######################
## Plot
#######################

visit2_bpb <- ggplot(data = data, aes(x = log(bpb), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed BPB (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(0.15/sqrt(2))) + geom_vline(xintercept = log(0.2/sqrt(2))) +
  geom_vline(xintercept = log(0.4/sqrt(2)))

###############
## BPF
###############

#CLMI
data <- visit2_small
cont_name <- "bpf"
batch_name <- "bpf_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.2, 0.1))
total_imp <- 10
seed_num <- 63327
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$bpf_batch_model <- as.numeric(data$bpf_batch)
data$bpf_batch_model[data$bpf_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2"))
res$output["bpf_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat,
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(bpf))"], sqrt(vcov(cca_res)["I(log(bpf))","I(log(bpf))"]), confint(cca_res)["I(log(bpf))",])

#MICE
mice_data <- visit2_small[, names(visit2_small) %in%
                            c("spont_preterm","bpf","sg","isage","childnum1","childnum2","currjob_cat")]
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat, family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(bpf))","est"], digits = 3), round(summary(imp)["I(log(bpf))","se"], digits = 3),
  round(summary(imp)["I(log(bpf))","lo 95"], digits = 3), round(summary(imp)["I(log(bpf))","hi 95"], digits = 3))

#SQRT 2
data$bpf[is.na(data$bpf) & data$bpf_lod == "<0.1"] <- 0.1/sqrt(2)
data$bpf[is.na(data$bpf) & data$bpf_lod == "<0.2"] <- 0.2/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(bpf))+sg+isage+childnum1+childnum2+currjob_cat,
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(bpf))"], sqrt(vcov(sq2_res)["I(log(bpf))","I(log(bpf))"]), confint(sq2_res)["I(log(bpf))",])

#######################
## Plot
#######################

visit2_bpf <- ggplot(data = data, aes(x = log(bpf), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed BPF (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(0.2/sqrt(2)))

###############
## TCS
###############

#CLMI
data <- visit2
cont_name <- "tcs"
batch_name <- "tcs_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(2.3, 1.7))
total_imp <- 10
seed_num <- 23569
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$tcs_batch_model <- as.numeric(data$tcs_batch)
data$tcs_batch_model[data$tcs_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","tcs_batch_model"))
res$output["tcs_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(tcs))"], sqrt(vcov(cca_res)["I(log(tcs))","I(log(tcs))"]), confint(cca_res)["I(log(tcs))",])

#MICE
mice_data <- visit2[, names(visit2) %in%
                      c("spont_preterm","tcs_batch","tcs","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$tcs_batch <- as.factor(mice_data$tcs_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(tcs))","est"], digits = 3), round(summary(imp)["I(log(tcs))","se"], digits = 3),
  round(summary(imp)["I(log(tcs))","lo 95"], digits = 3), round(summary(imp)["I(log(tcs))","hi 95"], digits = 3))

#SQRT 2
data$tcs[is.na(data$tcs) & data$tcs_lod == "<1.7"] <- 1.7/sqrt(2)
data$tcs[is.na(data$tcs) & data$tcs_lod == "<2.3"] <- 2.3/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(tcs))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcs_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(tcs))"], sqrt(vcov(sq2_res)["I(log(tcs))","I(log(tcs))"]), confint(sq2_res)["I(log(tcs))",])

#######################
## Plot
#######################

visit2_tcs <- ggplot(data = data, aes(x = log(tcs), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed TCS (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(1.7/sqrt(2))) +
  geom_vline(xintercept = log(2.3/sqrt(2)))

###############
## TCC
###############

#CLMI
data <- visit2_small
cont_name <- "tcc"
batch_name <- "tcc_batch"
outcome_name <- "spont_preterm"
imp_vars <- c("sg")
lod_info <- data.frame("batch_info" = c("1", "2"), "lod" = c(0.1, 1.7))
total_imp <- 10
seed_num <- 40043
transform_imp <- function(x){log(x)}

#Add batch indicator to model
data$tcc_batch_model <- as.numeric(data$tcc_batch)
data$tcc_batch_model[data$tcc_batch_model == 2] <- 0

imp <- clmi(data = data, cont_name = cont_name, batch_name = batch_name, outcome_name = outcome_name, imp_vars = imp_vars,
            lod_info = lod_info, total_imp = total_imp, seed_num = seed_num, transform_imp = transform_imp)
res <- pool.clmi(clmi_obj = imp, transform_model = function(x){log(x)}, regression_model = "Logistic",
                 precis_vars = c("isage","currjob_cat","childnum1","childnum2","tcc_batch_model"))
res$output["tcc_imputed",c("Est","SE","LCL.95","UCL.95")]

#CCA
cca_res <- glm(data = data, factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch),
               family = binomial)
summary(cca_res)
c(cca_res$coefficients["I(log(tcc))"], sqrt(vcov(cca_res)["I(log(tcc))","I(log(tcc))"]), confint(cca_res)["I(log(tcc))",])

#MICE
mice_data <- visit2_small[, names(visit2_small) %in%
                            c("spont_preterm","tcc_batch","tcc","sg","isage","childnum1","childnum2","currjob_cat")]
mice_data$tcc_batch <- as.factor(mice_data$tcc_batch)
mi <- mice(mice_data, m = 10, printFlag = FALSE)
comb <- with(mi, glm(factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch), family="binomial"))
imp <- pool(comb)
c(round(summary(imp)["I(log(tcc))","est"], digits = 3), round(summary(imp)["I(log(tcc))","se"], digits = 3),
  round(summary(imp)["I(log(tcc))","lo 95"], digits = 3), round(summary(imp)["I(log(tcc))","hi 95"], digits = 3))

#SQRT 2
data$tcc[is.na(data$tcc) & data$tcc_lod == "<1.7"] <- 1.7/sqrt(2)
data$tcc[is.na(data$tcc) & data$tcc_lod == "<0.1"] <- 0.1/sqrt(2)
sq2_res <- glm(data = data, factor(spont_preterm)~I(log(tcc))+sg+isage+childnum1+childnum2+currjob_cat+factor(tcc_batch),
               family = binomial)
summary(sq2_res)
c(sq2_res$coefficients["I(log(tcc))"], sqrt(vcov(sq2_res)["I(log(tcc))","I(log(tcc))"]), confint(sq2_res)["I(log(tcc))",])

#######################
## Plot
#######################

visit2_tcc <- ggplot(data = data, aes(x = log(tcc), color = factor(spont_preterm))) + geom_line(stat = "density", lwd = 1, aes(linetype = factor(spont_preterm))) +
  ylab("Density") + xlab("Log-Transformed TCC (ng/ml)") +
  scale_color_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("orange","darkblue")) +
  scale_linetype_manual(labels=c("Full Term Birth","Spontaneous Preterm Birth"), values = c("solid","dashed")) +
  labs(color = "", linetype = "") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = log(0.1/sqrt(2))) +
  geom_vline(xintercept = log(1.7/sqrt(2)))

#####################################################################
## Final Plot
#####################################################################

grid_arrange_shared_legend(visit1_mehp, visit2_mehp, visit1_mcpp, visit2_mcpp, visit1_bpb, visit2_bpb,
                           visit1_bpf, visit2_bpf, visit1_tcs, visit2_tcs, visit1_tcc, visit2_tcc, nrow = 6, ncol = 2)

png(file = "./cont_dist.png", width = 8, height = 11, res = 400, units = 'in')

grid_arrange_shared_legend(visit1_mehp, visit2_mehp, visit1_mcpp, visit2_mcpp, visit1_bpb, visit2_bpb,
                           visit1_bpf, visit2_bpf, visit1_tcs, visit2_tcs, visit1_tcc, visit2_tcc, nrow = 6, ncol = 2)

dev.off()

