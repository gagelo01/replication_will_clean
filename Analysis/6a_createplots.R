#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(ckbplotr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
res_class <- fread( "Data/Modified/res_class.txt")
harm_class <- fread("Data/Modified/harm_class.txt")
###plotting Figure 1 and supplementary figure 1

###plotting Figure 1 and supplementary figure 1
exp_name = list(c("BMI_UKB", "WC_UKB", "whradjbmi_ukbonly"), c("GIANT_2015_BMI", "GIANT_2015_WC","WHRadjBMI"))
out_name = list(c("NAFLD","Fat_Liver"), c("NAFLD", "Fat_Liver"))
file_name <- c("Figure1", "SupplementaryFigure1")

for(i in 1:length(exp_name)) {
  doA <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[[i]][1]]
  # doA <- doA[!(method == "MR Egger")]
  doA[,exposure := gsub("BMI_UKB", "Body mass index", exposure)]
  doA[,exposure := gsub("WC_UKB", "Waist circumference", exposure)]
  doA[,exposure := gsub("GIANT_2015_BMI", "Body mass index", exposure)]
  doA[,exposure := gsub("GIANT_2015_WC", "Waist circumference", exposure)]
  doA[, exposure := gsub("whradjbmi_ukbonly", "WHRadjBMI", exposure)]
  doA[,colour := "black"]
  doB <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[[i]][2]]
  
  resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                         estimate = round(doA$b, digits =2),
                         lci =  round(doA$lci, digits = 2),
                         uci =  round(doA$uci, digits = 2),
                         n = doA$nsnp,
                         P_value = formatC(doA$pval, format = "e", digits = 1),
                         colour = doA$colour )
  
  resultsB <- data.frame(variable = as.character(1:nrow(doB)),
                         estimate = round(doB$b, digits =2),
                         lci =  round(doB$lci, digits = 2),
                         uci =  round(doB$uci, digits = 2),
                         n = doB$nsnp,
                         P_value = formatC(doB$pval, format = "e", digits = 1),
                         colour = doA$colour )
  
  mylabels <- data.frame(heading1 = doA$exposure,
                         heading2 = doA$method,
                         heading3 = as.character(NA),
                         variable = as.character(1:nrow(doA)))
  
  
  k <- make_forest_plot(panels = list(resultsA, resultsB),
                        col.key = "variable",
                        row.labels = mylabels,
                        exponentiate = FALSE,
                        pointsize = 2, 
                        rows = unique(mylabels$heading1),
                        col.stderr = NULL,
                        col.lci = "lci",
                        col.uci = "uci",
                        col.left         = c("n"),
                        col.left.heading = c("n SNPs"),
                        col.right = "P_value",
                        col.right.heading = c("Effect (95% CI)", "P-value"),
                        xlab = c("Effect on NAFLD log(OR)", "Effect on Liver Fat (SD)"),
                        blankrows = c(0,0,0,0),
                        col.right.hjust = 1,
                        nullval = 0,
                        colour           = "colour",
                        panel.headings = NULL)
  k  
  ggsave(paste0("Results/", file_name[i], ".png"),
         width=936/72,height=529/72, units="in", scale=1,  dpi = 1000,
         device = "png")
}
#Figure 2 clean
mvmr_object <- list(c("WC_UKB-and-BMI_UKB-on-NAFLD", "WC_UKB-and-BMI_UKB-on-Fat_Liver"), c("GIANT_2015_WC-and-GIANT_2015_BMI-on-NAFLD", "GIANT_2015_WC-and-GIANT_2015_BMI-on-Fat_Liver"))
file_name <- c("Figure2", "SupplementaryFigure2")
for(i in 1:length(file_name)) {
  mvmr_results <- lapply(as.list(mvmr_object[[i]]), function(x) resmvmr[[x]]) %>% rbindlist(.)
  
  # k  <- gsub("-and-|-on-", ",", mvmr_object[[i]])
  # k <- strsplit(k, split = ",")   %>% unlist
  # uni <- res_univariate[exposure %in% k[1:2] & outcome == k[3], ]
  uni <- res_univariate[exposure %in% mvmr_results$exposure & outcome %in% mvmr_results$outcome, ]
  mvmr_results <-  mvmr_results[clump_exposure=="none", ]
  mvmr_results <- rbindlist(list(uni, mvmr_results), fill = TRUE)
  
  unimeth<-"Inverse variance weighted" 
  multimeth<- c("Multivariable IVW", "Multivariable Median",
                "Multivariable Lasso", "Multivariable Egger")  
  
  data <- mvmr_results[method %in% c(unimeth, multimeth),]
  data[, Category_other := ifelse(method %in% unimeth, "Univariable", "Multivariable")]
  data[, Category_other := factor(Category_other, levels = c("Univariable", "Multivariable"))]
  data[,name := gsub("BMI_UKB|GIANT_2015_BMI", "BMI", exposure) %>% gsub("WC_UKB|GIANT_2015_WC", "Waist circumference", .)]
  data[Category_other == "Multivariable", name := name %>% ifelse(. == "BMI", "BMI adjusted for Waist circumference", .) %>% ifelse(. == "Waist circumference", "Waist circumference adjusted for BMI", .)]
  data[, outcome := factor(outcome, levels = c("NAFLD", "Fat_Liver"))]
  data<-data[order(Category_other,exposure,outcome)]
  ggforestplot::forestplot(
    df = data,
    name = name,
    se = se,
    estimate = b,
    pvalue = pval,
    psignif = 0.05,
    xlab = "Effect size (SD) per 1-SD\nincrease in adiposity",
    ci = 0.95,
    logodds = FALSE,
    colour = method
  ) + facet_grid( ~ outcome ) +
    theme(legend.position = "right") +
    theme(text = element_text(size = 10)) +
    theme(legend.position="right") #+
  # ggforce::facet_col(
  #   # facets = ~outcome_category,
  #   facets = ~Category_other,
  #   scales = "free_y",
  #   space = "free"
  # ) + facet_grid( ~ outcome ) 
  
  ggsave(paste0("Results/", file_name[i],"colour", ".png"),
         width=483/72,height=250/72, units="in", scale=1, dpi = 1000,
         device = "png")
  #### figure 2
  
  doA <- data[outcome=="NAFLD",]
  doB<-  data[outcome=="Fat_Liver",]
  resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                         estimate = round(doA$b, digits =2),
                         lci =  round(doA$lci, digits = 2),
                         uci =  round(doA$uci, digits = 2),
                         n = doA$nsnp,
                         P_value = formatC(doA$pval, format = "e", digits = 1))
  
  resultsB <- data.frame(variable = as.character(1:nrow(doB)),
                         estimate = round(doB$b, digits =2),
                         lci =  round(doB$lci, digits = 2),
                         uci =  round(doB$uci, digits = 2),
                         n = doB$nsnp,
                         P_value = formatC(doB$pval, format = "e", digits = 1))
  
  mylabels <- data.frame(heading1 = doA$name,
                         heading2 = doA$method,
                         heading3 = as.character(NA),
                         variable = as.character(1:nrow(doA)))
  
  ckbplotr::make_forest_plot(panels = list(resultsA,resultsB),
                             col.key = "variable",
                             row.labels = mylabels,
                             exponentiate = FALSE,
                             pointsize = 2, 
                             rows = unique(mylabels$heading1),
                             col.stderr = NULL,
                             col.lci = "lci",
                             col.uci = "uci",
                             col.left         = c("n"),
                             col.left.heading = c("n SNPs"),
                             col.right = "P_value",
                             col.right.heading = c("Effect (95% CI)", "P-value"),
                             xlab = c("Effect on NAFLD log(OR)", "Effect on Liver fat (SD)"),
                             blankrows = c(0,0,0,0),
                             col.right.hjust = 1,
                             panel.headings = NULL,
                             nullval = 0,
                             scalepoints = FALSE)
  
  ggsave(paste0("Results/", file_name[i],"", ".png"),
         width=840/72,height=346/72, units="in", scale=1,  dpi = 1000,
         device = "png")
}

#Figure 3
harm_class <- harm_class[class != "WHRonly-"] #there is no value at doing the analysis on WHRonly-
harm_class[,class := NULL]
res_class <- separate(res_class, col = "exposure", into = c("class", "exposure"), extra = "merge", sep = "_", remove =  TRUE) %>% as.data.table(.)
harm_class <- separate(harm_class, col = "exposure", into = c("class", "exposure"), extra = "merge", sep = "_", remove =  TRUE) %>% as.data.table(.)
res_class[, id.exposure := exposure]
harm_class[, id.exposure := exposure]

mr_results <-  res_class[exposure == "BMI_UKB" & outcome %in% c("NAFLD", "Fat_Liver"),]
dat <- harm_class[exposure == "BMI_UKB" & outcome %in% c("NAFLD", "Fat_Liver"),]

mr_results[,align := class]
dat[,align := class]
dat[, Locus := ""]
dat[, align := factor(align, levels = c("BMI+WHR+","BMI+WHR-","BMIonly+"))]
dat[, outcome := factor(outcome, levels = c("NAFLD","Fat_Liver"))]
dat <- dat[order(align,outcome),]
mr_results[, align := factor(align, levels = levels(dat$align))]
mr_results[, outcome := factor(outcome, levels = levels(dat$outcome))]
mr_results <- mr_results[order(align,outcome),]

source("/mnt/sda/gagelo01/Projects/small_MR_exploration/FI_BMI/Analysis/my_mr_scatter_plot.R")
k <- my_mr_scatter_plot( dat = dat, mr_results = mr_results, equation_facet_grid = "outcome ~  align", legend.position = "top")
k + xlab("SNP effect on BMI") +ylab("SNP effect on liver traits")

ggsave(paste0("Results/", "Figure3", ".png"),
       width=530/72,height=418/72, units="in", scale=1, dpi = 1000,
       device = "png")


#Figure 4
exp_name <- list(c("WC_UKB", "Fat_Liver"), c("WC_UKB", "Fat_Liver"), c("GIANT_2015_WC", "Fat_Liver"), c("GIANT_2015_WC", "Fat_Liver"))
out_name <-  c("Mahajan_Type2diabetes", "van_der_Harst_CAD", "Mahajan_Type2diabetes", "van_der_Harst_CAD") 
study<-c("ukb", "giant")
file_name <- c("Figure4", "Figure4_CAD", "SupplementaryFigure3", "SupplementaryFigure3_CAD")

for(i in 1:length(exp_name)) {
  
  doA<-res_univariate[outcome == out_name[i] & exposure %in% exp_name[[i]]]
  
  
  format_toforest <- function(do, multivariate) {
    do[,heading := "Univariable"]
    multivariate <- multivariate[method != "Multivariable Egger Intercept"]
    multivariate[,heading := "Multivariable"]
    multivariate <- multivariate[clump_exposure == "Fat_Liver", ][order(exposure)]
    MVMR<- rbindlist(list(do, multivariate), use.names = TRUE, fill = TRUE)
    MVMR[, name := exposure %>% gsub("WC_UKB|GIANT_2015_WC", "Waist circumference", . ) %>% gsub("Fat_Liver", "Liver fat", .)]
    return(MVMR)
  }
  
  MVMRA <- format_toforest(do = doA, multivariate = resmvmr[[paste0(exp_name[[i]][1],"-and-Fat_Liver-on-",out_name[i])]])
  MVMRA[heading == "Multivariable", name := name %>% ifelse(. == "Liver fat", "Liver fat adjusted for Waist circumference", .) %>% ifelse(. == "Waist circumference", "Waist circumference adjusted for Liver fat", .)]
  
  resultsA <- data.frame(variable = LETTERS[1:nrow(MVMRA)],
                         estimate = round(MVMRA$b, digits =2),
                         stderr = round(MVMRA$se, digits = 2),
                         lci =  round(MVMRA$lci, digits = 2),
                         uci =  round(MVMRA$uci, digits = 2),
                         P_value = formatC(MVMRA$pval, format = "e", digits = 1),
                         nsnp = MVMRA$nsnp,
                         colour = "black")
  
  mylabels <- data.frame(heading = as.character(MVMRA$name),
                         subheading = MVMRA$method,
                         label = as.character(NA),
                         # label = as.character(MVMRA$exposure),
                         variable = LETTERS[1:nrow(MVMRA)])
  
  make_forest_plot(panels = list(resultsA),
                   col.key = "variable",
                   row.labels = mylabels,
                   row.labels.levels = c("heading", "subheading", "label"),
                   rows             = unique(mylabels$heading),
                   panel.headings = NULL,
                   exponentiate = TRUE,
                   pointsize = 2,
                   col.stderr = NULL,
                   col.lci = "lci",
                   col.uci = "uci",
                   col.right = "P_value",
                   col.right.heading = c("OR (95% CI)", "P-value"),
                   col.left         = c("nsnp"),
                   col.left.heading = c("n SNPs"),
                   xlab = paste0("Effect of 1-SD increase in Liver fat/waist circumference on ",ifelse(out_name[i] == "van_der_Harst_CAD", "CAD", "T2D") ," (OR)"),
                   colour           = "colour",
                   blankrows = c(0,0,0,0))
  
  
  ggsave(paste0("Results/", file_name[i], ".png"),
         width=798/72,height=583/72, units="in", scale=1, dpi = 1000,
         device = "png")
}

#########Figures for presentation univariable WC + liver fat ~ CAD

exp_name = list(c("WC_UKB", "Fat_Liver"))
out_name = list("van_der_Harst_CAD")
file_name <- c("Presentationliver_WC_CAD")

for(i in 1:length(exp_name)) {
  doA <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[[i]][1]]
  # doA <- doA[!(method == "MR Egger")]
  doA[,exposure := gsub("BMI_UKB", "Body mass index", exposure)]
  doA[,exposure := gsub("WC_UKB", "Waist circumference", exposure)]
  doA[,exposure := gsub("GIANT_2015_BMI", "Body mass index", exposure)]
  doA[,exposure := gsub("GIANT_2015_WC", "Waist circumference", exposure)]
  doA[,exposure := gsub("Fat_Liver", "Liver Fat", exposure)]
  doA[,colour := "black"]
  
  resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                         estimate = round(doA$b, digits =2),
                         lci =  round(doA$lci, digits = 2),
                         uci =  round(doA$uci, digits = 2),
                         n = doA$nsnp,
                         P_value = formatC(doA$pval, format = "e", digits = 1),
                         colour = doA$colour )
  
  
  mylabels <- data.frame(heading1 = doA$exposure,
                         heading2 = doA$method,
                         heading3 = as.character(NA),
                         variable = as.character(1:nrow(doA)))
  
  
  k <- make_forest_plot(panels = list(resultsA),
                        col.key = "variable",
                        row.labels = mylabels,
                        exponentiate = FALSE,
                        pointsize = 2, 
                        rows = unique(mylabels$heading1),
                        col.stderr = NULL,
                        col.lci = "lci",
                        col.uci = "uci",
                        col.left         = c("n"),
                        col.left.heading = c("n SNPs"),
                        col.right = "P_value",
                        col.right.heading = c("Effect (95% CI)", "P-value"),
                        xlab = c("Effect on CAD log(OR)"),
                        blankrows = c(0,0,0,0),
                        col.right.hjust = 1,
                        nullval = 0,
                        colour           = "colour",
                        panel.headings = NULL)
  k  
  ggsave(paste0("Results/", file_name[i], ".png"),
         width=936/72,height=529/72, units="in", scale=1, dpi = 1000,
         device = "png")
}

message("This script finished without errors")

