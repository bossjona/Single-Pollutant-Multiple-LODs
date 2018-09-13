library(ggplot2)
library(grid)
library(gridExtra)
library(xlsx)
library(reshape2)
library(tidyr)
library(ggpubr)

setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")

#Function that merges ggplot panels with the same legend
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

gen.sim.plot <- function(design, large){
  #Read in files corresponding to each study design
  if(design == "Cohort" & large == TRUE){
    randbatch_withind <- read.xlsx("largecohort_randombatch_withind.xlsx", sheetIndex = 1)
    randbatch_withoutind <- read.xlsx("largecohort_randombatch_withoutind.xlsx", sheetIndex = 1)
    depbatch_withind <- read.xlsx("largecohort_depbatch_withind.xlsx", sheetIndex = 1)
    depbatch_withoutind <- read.xlsx("largecohort_depbatch_withoutind.xlsx", sheetIndex = 1)
  } else if(design == "Cohort" & large == FALSE){
    randbatch_withind <- read.xlsx("modcohort_randombatch_withind.xlsx", sheetIndex = 1)
    randbatch_withoutind <- read.xlsx("modcohort_randombatch_withoutind.xlsx", sheetIndex = 1)
    depbatch_withind <- read.xlsx("modcohort_depbatch_withind.xlsx", sheetIndex = 1)
    depbatch_withoutind <- read.xlsx("modcohort_depbatch_withoutind.xlsx", sheetIndex = 1)
  } else if(design == "Case Control" & large == TRUE){
    randbatch_withind <- read.xlsx("largecc_randombatch_withind.xlsx", sheetIndex = 1)
    randbatch_withoutind <- read.xlsx("largecc_randombatch_withoutind.xlsx", sheetIndex = 1)
    depbatch_withind <- read.xlsx("largecc_depbatch_withind.xlsx", sheetIndex = 1)
    depbatch_withoutind <- read.xlsx("largecc_depbatch_withoutind.xlsx", sheetIndex = 1)
  } else if(design == "Case Control" & large == FALSE){
    randbatch_withind <- read.xlsx("modcc_randombatch_withind.xlsx", sheetIndex = 1)
    randbatch_withoutind <- read.xlsx("modcc_randombatch_withoutind.xlsx", sheetIndex = 1)
    depbatch_withind <- read.xlsx("modcc_depbatch_withind.xlsx", sheetIndex = 1)
    depbatch_withoutind <- read.xlsx("modcc_depbatch_withoutind.xlsx", sheetIndex = 1)
  }
  
  #######################################################################################
  # Bias
  #######################################################################################
  
  #Extract Bias
  bias_randbatch_withind_tmp <- randbatch_withind[, names(randbatch_withind) %in%
                                                c("lod_name","gs_rel_bias","cca_rel_bias","sq2_rel_bias","mice_rel_bias","pmi_rel_bias")]
  bias_randbatch_withoutind_tmp <- randbatch_withoutind[, names(randbatch_withoutind) %in%
                                                          c("lod_name","gs_rel_bias_wob","cca_rel_bias_wob","sq2_rel_bias_wob","mice_rel_bias_wob","pmi_rel_bias_wob")]
  bias_depbatch_withind_tmp <- depbatch_withind[, names(depbatch_withind) %in%
                                                  c("lod_name","gs_rel_bias","cca_rel_bias","sq2_rel_bias","mice_rel_bias","pmi_rel_bias")]
  bias_depbatch_withoutind_tmp <- depbatch_withoutind[, names(depbatch_withoutind) %in%
                                                        c("lod_name","gs_rel_bias_wob","cca_rel_bias_wob","sq2_rel_bias_wob","mice_rel_bias_wob","pmi_rel_bias_wob")]
  #Stack Bias
  bias_randbatch_withind <- gather(bias_randbatch_withind_tmp[, !(names(bias_randbatch_withind_tmp) %in% c("gs_rel_bias"))],
                                   method, bias, cca_rel_bias:pmi_rel_bias, factor_key = TRUE)
  bias_randbatch_withoutind <- gather(bias_randbatch_withoutind_tmp[, !(names(bias_randbatch_withoutind_tmp) %in% c("gs_rel_bias_wob"))],
                                      method, bias, cca_rel_bias_wob:pmi_rel_bias_wob, factor_key = TRUE)
  bias_depbatch_withind <- gather(bias_depbatch_withind_tmp[, !(names(bias_depbatch_withind_tmp) %in% c("gs_rel_bias"))],
                                   method, bias, cca_rel_bias:pmi_rel_bias, factor_key = TRUE)
  bias_depbatch_withoutind <- gather(bias_depbatch_withoutind_tmp[, !(names(bias_depbatch_withoutind_tmp) %in% c("gs_rel_bias_wob"))],
                                     method, bias, cca_rel_bias_wob:pmi_rel_bias_wob, factor_key = TRUE)
  
  #Covert to numeric
  bias_randbatch_withind$bias <- as.numeric(bias_randbatch_withind$bias)
  bias_randbatch_withoutind$bias <- as.numeric(bias_randbatch_withoutind$bias)
  bias_depbatch_withind$bias <- as.numeric(bias_depbatch_withind$bias)
  bias_depbatch_withoutind$bias <- as.numeric(bias_depbatch_withoutind$bias)
  
  if(design == "Cohort" & large == TRUE){
    bias_bounds <- c(-15, 20, 5)
    mse_bounds <- c(0, 0.01, 0.005)
  } else if(design == "Cohort" & large == FALSE){
    bias_bounds <- c(-20, 20, 5)
    mse_bounds <- c(0, 0.05, 0.01)
  } else if(design == "Case Control" & large == TRUE){
    bias_bounds <- c(-60, 50, 10)
    mse_bounds <- c(0, 0.065, 0.01)
  } else if(design == "Case Control" & large == FALSE){
    bias_bounds <- c(-75, 60, 15)
    mse_bounds <- c(0, 0.11, 0.02)
  }
  
  #Plot Bias
  b1 <- ggplot(data = bias_randbatch_withoutind, aes(x = lod_name, y = bias, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("A)          Without Batch Indicator - Independent") + 
    ylab("Relative Bias (%)") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks=round(seq(from = bias_bounds[1], to = bias_bounds[2], by = bias_bounds[3]), digits = 2), limits=c(bias_bounds[1],bias_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  b2 <- ggplot(data = bias_randbatch_withind, aes(x = lod_name, y = bias, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("B)             With Batch Indicator - Independent") + 
    ylab("Relative Bias (%)") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks=round(seq(from = bias_bounds[1], to = bias_bounds[2], by = bias_bounds[3]), digits = 2), limits=c(bias_bounds[1],bias_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  b3 <- ggplot(data = bias_depbatch_withoutind, aes(x = lod_name, y = bias, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("C)           Without Batch Indicator - Dependent") + 
    ylab("Relative Bias (%)") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks=round(seq(from = bias_bounds[1], to = bias_bounds[2], by = bias_bounds[3]), digits = 2), limits=c(bias_bounds[1],bias_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  b4 <- ggplot(data = bias_depbatch_withind, aes(x = lod_name, y = bias, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("D)              With Batch Indicator - Dependent") + 
    ylab("Relative Bias (%)") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks=round(seq(from = bias_bounds[1], to = bias_bounds[2], by = bias_bounds[3]), digits = 2), limits=c(bias_bounds[1],bias_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18\\Final Plots")
  
  if(design == "Cohort" & large == TRUE){
    png(file = "./grid_large_cohort_rel_bias.png", width = 9, height = 9, res = 400, units = 'in')
    grid_arrange_shared_legend(b1,b2,b3,b4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Cohort" & large == FALSE){
    pdf(file = "./grid_mod_cohort_rel_bias.pdf", width = 9, height = 9, onefile=FALSE)
    grid_arrange_shared_legend(b1,b2,b3,b4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Case Control" & large == TRUE){
    png(file = "./grid_large_cc_rel_bias.png", width = 9, height = 9, res = 400, units = 'in')
    grid_arrange_shared_legend(b1,b2,b3,b4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Case Control" & large == FALSE){
    pdf(file = "./grid_mod_cc_rel_bias.pdf", width = 9, height = 9, onefile=FALSE)
    grid_arrange_shared_legend(b1,b2,b3,b4, ncol = 2, nrow = 2)
    dev.off()
  }
  
  setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")
  
  #######################################################################################
  # Mean-Squared Error
  #######################################################################################
  
  #Extract Mean-Squared Error
  mse_randbatch_withind_tmp <- randbatch_withind[, names(randbatch_withind) %in%
                                                    c("lod_name","gs_mse","cca_mse","sq2_mse","mice_mse","pmi_mse")]
  mse_randbatch_withoutind_tmp <- randbatch_withoutind[, names(randbatch_withoutind) %in%
                                                          c("lod_name","gs_mse_wob","cca_mse_wob","sq2_mse_wob","mice_mse_wob","pmi_mse_wob")]
  mse_depbatch_withind_tmp <- depbatch_withind[, names(depbatch_withind) %in%
                                                  c("lod_name","gs_mse","cca_mse","sq2_mse","mice_mse","pmi_mse")]
  mse_depbatch_withoutind_tmp <- depbatch_withoutind[, names(depbatch_withoutind) %in%
                                                        c("lod_name","gs_mse_wob","cca_mse_wob","sq2_mse_wob","mice_mse_wob","pmi_mse_wob")]
  
  #Stack MSE
  mse_randbatch_withind <- gather(mse_randbatch_withind_tmp[, !(names(mse_randbatch_withind_tmp) %in% c("gs_mse"))],
                                   method, mse, cca_mse:pmi_mse, factor_key = TRUE)
  mse_randbatch_withoutind <- gather(mse_randbatch_withoutind_tmp[, !(names(mse_randbatch_withoutind_tmp) %in% c("gs_mse_wob"))],
                                      method, mse, cca_mse_wob:pmi_mse_wob, factor_key = TRUE)
  mse_depbatch_withind <- gather(mse_depbatch_withind_tmp[, !(names(mse_depbatch_withind_tmp) %in% c("gs_mse"))],
                                  method, mse, cca_mse:pmi_mse, factor_key = TRUE)
  mse_depbatch_withoutind <- gather(mse_depbatch_withoutind_tmp[, !(names(mse_depbatch_withoutind_tmp) %in% c("gs_mse_wob"))],
                                     method, mse, cca_mse_wob:pmi_mse_wob, factor_key = TRUE)
  
  #Convert to Numeric
  mse_randbatch_withind$mse <- as.numeric(mse_randbatch_withind$mse)
  mse_randbatch_withoutind$mse <- as.numeric(mse_randbatch_withoutind$mse)
  mse_depbatch_withind$mse <- as.numeric(mse_depbatch_withind$mse)
  mse_depbatch_withoutind$mse <- as.numeric(mse_depbatch_withoutind$mse)
  
  #Get Gold Standard MSE for Comparison
  comp1 <- as.numeric(as.character(randbatch_withoutind$gs_mse))[1]
  comp2 <- as.numeric(as.character(randbatch_withind$gs_mse))[1]
  comp3 <- as.numeric(as.character(depbatch_withoutind$gs_mse))[1]
  comp4 <- as.numeric(as.character(depbatch_withind$gs_mse))[1]
  
  #Plot MSE
  m1 <- ggplot(data = mse_randbatch_withoutind, aes(x = lod_name, y = mse, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("A)          Without Batch Indicator - Independent") + 
    ylab("MSE") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = comp1) +
    scale_y_continuous(breaks=round(seq(from = mse_bounds[1], to = mse_bounds[2], by = mse_bounds[3]), digits = 4), limits=c(mse_bounds[1],mse_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  m2 <- ggplot(data = mse_randbatch_withind, aes(x = lod_name, y = mse, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("B)             With Batch Indicator - Independent") + 
    ylab("MSE") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = comp2) +
    scale_y_continuous(breaks=round(seq(from = mse_bounds[1], to = mse_bounds[2], by = mse_bounds[3]), digits = 4), limits=c(mse_bounds[1],mse_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  m3 <- ggplot(data = mse_depbatch_withoutind, aes(x = lod_name, y = mse, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("C)           Without Batch Indicator - Dependent") + 
    ylab("MSE") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = comp3) +
    scale_y_continuous(breaks=round(seq(from = mse_bounds[1], to = mse_bounds[2], by = mse_bounds[3]), digits = 4), limits=c(mse_bounds[1],mse_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  m4 <- ggplot(data = mse_depbatch_withind, aes(x = lod_name, y = mse, color = factor(method), shape = factor(method))) +
    geom_point(size = 3) +
    ggtitle("D)             With Batch Indicator - Dependent") + 
    ylab("MSE") + xlab("LOD Pairs") +
    #scale_shape_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c(15,16,17,18)) +
    #scale_color_manual(name = "Method:", labels = c("CCA", expression(LOD/sqrt(2)), "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    scale_shape_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c(15,16,17,18)) +
    scale_color_manual(name = "Method:", labels = c("CCA", "SQ2", "MICE", "CLMI"), values = c("blue","forestgreen","red","goldenrod")) +
    geom_hline(yintercept = comp4) +
    scale_y_continuous(breaks=round(seq(from = mse_bounds[1], to = mse_bounds[2], by = mse_bounds[3]), digits = 4), limits=c(mse_bounds[1],mse_bounds[2])) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle=45, vjust=0.5), plot.title = element_text(hjust = -0.3, vjust = 2, size = 12))
  
  #######################################################################################
  # Generate Bias and MSE Plot
  #######################################################################################
  
  setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18\\Final Plots")
  
  if(design == "Cohort" & large == TRUE){
    png(file = "./grid_large_cohort_mse.png", width = 9, height = 9, res = 400, units = 'in')
    grid_arrange_shared_legend(m1,m2,m3,m4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Cohort" & large == FALSE){
    pdf(file = "./grid_mod_cohort_mse.pdf", width = 9, height = 9, onefile=FALSE)
    grid_arrange_shared_legend(m1,m2,m3,m4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Case Control" & large == TRUE){
    png(file = "./grid_large_cc_mse.png", width = 9, height = 9, res = 400, units = 'in')
    grid_arrange_shared_legend(m1,m2,m3,m4, ncol = 2, nrow = 2)
    dev.off()
  } else if(design == "Case Control" & large == FALSE){
    pdf(file = "./grid_mod_cc_mse.pdf", width = 9, height = 9, onefile=FALSE)
    grid_arrange_shared_legend(m1,m2,m3,m4, ncol = 2, nrow = 2)
    dev.off()
  }
  
  setwd("C:\\Users\\John Boss\\Desktop\\Apples and Oranges Backup\\Jonathan\\Documents\\GSRA Presentations\\Single Pollutant Multiple LOD\\Final Results - Updated 6-1-18")
  
}

#######################################################################################
# Large Cohort
#######################################################################################

design <- "Cohort"
large <- TRUE
gen.sim.plot(design = design, large = large)

#######################################################################################
# Moderate Cohort
#######################################################################################

design <- "Cohort"
large <- FALSE
gen.sim.plot(design = design, large = large)

#######################################################################################
# Large Case-Control
#######################################################################################

design <- "Case Control"
large <- TRUE
gen.sim.plot(design = design, large = large)

#######################################################################################
# Small Case-Control
#######################################################################################

design <- "Case Control"
large <- FALSE
gen.sim.plot(design = design, large = large)

