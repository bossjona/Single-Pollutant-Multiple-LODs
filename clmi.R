clmi <- function(data, cont_name, batch_name, outcome_name, imp_vars, lod_info, total_imp = 5, seed_num, transform_imp = function(x){x}){
  #data: Data frame with contaminant concentration, batch number, covariates used in imputation, precision variables
  #(make sure categorical variables are numeric, not factors). Can use model.matrix function to convert data frame with
  #factors to numeric design matrix and convert that matrix back into a data frame.
  #cont_name: string corresponding to name of untransformed contaminant variable.
  #batch_name: string corresponding to name of batch indicator variable.
  #outcome_name: string corresponding to name of continuous or binary health outcome.
  #imp_vars: vector of strings corresponding to covariates associated with the contaminant (covariates used in imputation).
  #lod_info: data frame with first column corresponding to each batch (entered as a string) and the
  #second column listing the LODs (entered numerically) corresponding to each batch.
  #First column must be named "batch_info" and the second column must be named "lod".
  #total_imp: number of multiply imputed datasets (positive integer).
  #seed_num: seed number to ensure reproducability.
  #transform_imp: transformation applied to contaminant data to make it normally distributed, default is no transformation.
  #Resulting model will be fit on the transformed contaminant.
  
  set.seed(seed_num)
  
  ################################################################################
  ## Organize Data into non-detects and detects
  ################################################################################
  
  data$lod <- rep(0, nrow(data))
  for(i in 1:nrow(data)){
    data$lod[i] <- lod_info[which(lod_info$batch_info == data[[batch_name]][i]),"lod"]
  }
  
  cont_imp_name <- paste0(cont_name,"_transform","_imputed")
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
    
    mle <- optim(c(0,0,1), dist_x_given_yc, hessian=TRUE)
    mle_param <- mle$par
    fisher_inf <- -solve(-mle$hessian)
    prop_sigma <- sqrt(diag(fisher_inf))
  }
  
  #Aggregate imputed datasets with observed dataset
  poll_imp <- imp
  above_lod[[cont_imp_name]] <- transform_imp(above_lod[[cont_imp_name]])
  for(k1 in 1:length(poll_imp)){
    poll_imp[[k1]] <- rbind(imp[[k1]], above_lod)
  }
  return(list("nimp" = total_imp, "imputations" = poll_imp, "transform" = transform_imp,
              "mle" = mle_param, "fisher_inf" = fisher_inf,
              var_names = list("cont_name" = cont_name, "batch_name" = batch_name, "outcome_name" = outcome_name,
                               "imp_vars" = imp_vars, "lod_info" = lod_info, "cont_imp_name" = cont_imp_name)))
}

pool.clmi <- function(clmi_obj, regression_model, precis_vars){
  #clmi_obj: List generated from clmi() function is used as input for this argument
  #regression_model: Linear Regression or Logistic Regression
  #precis_vars: vector of strings corresponding to adjustment covariates
  
  #Get names used in clmi() function
  poll_imp <- clmi_obj$imputations
  batch_name <- clmi_obj$var_names$batch_name
  cont_imp_name <- clmi_obj$var_names$cont_imp_name
  outcome_name <- clmi_obj$var_names$outcome_name
  imp_vars <- clmi_obj$var_names$imp_vars
  total_imp <- clmi_obj$nimp
  nobs <- nrow(poll_imp[[1]])
  
  #Get number of coefficients to store
  data.fit <- poll_imp[[1]][, names(poll_imp[[1]]) %in% imp_vars | names(poll_imp[[1]]) %in% precis_vars |
                              names(poll_imp[[1]]) %in% cont_imp_name]
  if(is.vector(data.fit)){
    num_continuous <- 1
  } else if(!is.vector(data.fit)){
    num_continuous <- sum(sapply(data.fit, is.numeric))
  }
  #Total number of parameters in logistic regression (add one for the intercept)
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

