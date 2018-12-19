setwd("Enter your working directory")

source("clmi.R")

test <- read.table("toy_dataset.txt", stringsAsFactors = FALSE)

save_imp <- clmi(data = test, cont_name = "poll", batch_name = "batch1", outcome_name = "case_cntrl",
                 imp_vars = c("smoking","gender"),
                 lod_info = data.frame(batch_info = c("1","0"), lod = c(0.8,0.65)),
                 total_imp = 5, seed_num = 12345, transform_imp = function(x){log(x)})

results <- pool.clmi(clmi_obj = save_imp, regression_model = "Logistic", precis_vars = NULL)

results$output