library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(ckbplotr)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
res_class <- fread( "Data/Modified/res_class.txt")
harm_class <- fread("Data/Modified/harm_class.txt")
###plotting Figure 1 and supplementary figure 1
exp_name = list(c("BMI_UKB", "WC_UKB", "WHRadjBMI"), c("GIANT_2015_BMI", "GIANT_2015_WC" ))
out_name = c("NAFLD", "NAFLD")
file_name <- c("Figure1", "Supplementary_figure1")
for(i in 1:length(exp_name)) {
  doA <- res_univariate[exposure %in% exp_name[[i]] & outcome %in% out_name[i]]
  # doA <- doA[!(method == "MR Egger")]
  doA[,exposure := gsub("BMI_UKB", "Body mass index", exposure)]
  doA[,exposure := gsub("WC_UKB", "Waist circumference", exposure)]
  doA[,exposure := gsub("GIANT_2015_BMI", "Body mass index", exposure)]
  doA[,exposure := gsub("GIANT_2015_WC", "Waist circumference", exposure)]
  doA[,colour := "black"]
  # doB <- doB[!(method == "MR Egger")]
  
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
  
  
 make_forest_plot(panels = list(resultsA),
                   col.key = "variable",
                   row.labels = mylabels,
                   exponentiate = TRUE,
                   pointsize = 2, 
                   rows = unique(mylabels$heading1),
                   col.stderr = NULL,
                   col.lci = "lci",
                   col.uci = "uci",
                   col.left         = c("n"),
                   col.left.heading = c("n SNP"),
                   col.right = "P_value",
                   col.right.heading = c("Effect (95% CI)", "P-value"),
                   xlab = c("Effect of 1 SD increase in anthropometric measures on NAFLD (OR)"),
                   blankrows = c(0,0,0,0),
                   col.right.hjust = 1,
                   nullval = 0,
                   colour           = "colour",
                   panel.headings = NULL)

  ggsave(paste0("Results/", file_name[i], ".png"),
         width=1000/72,height=472/72, units="in", scale=1,
         device = "png")
}

  #Figure 2 clean
mvmr_object <- c("WC_UKB-and-BMI_UKB-on-NAFLD", "GIANT_2015_WC-and-GIANT_2015_BMI-on-NAFLD")
file_name <- c("Figure2", "Supplementary_figure2")
for(i in 1:length(file_name)) {
 mvmr_results <- resmvmr[[mvmr_object[i]]]
  
k  <- gsub("-and-|-on-", ",", mvmr_object[i])
k <- strsplit(k, split = ",")   %>% unlist
uni <- res_univariate[exposure %in% k[1:2] & outcome == k[3], ]
mvmr_results <-  mvmr_results[clump_exposure=="none", ]
mvmr_results <- rbindlist(list(uni, mvmr_results), fill = TRUE)

unimeth<-"Inverse variance weighted" 
  multimeth<- c("Multivariable IVW", "Multivariable Median",
                "Multivariable Lasso", "Multivariable Egger")  

  data <- mvmr_results[method %in% c(unimeth, multimeth),]
  data[, Category_other := ifelse(method %in% unimeth, "Univariable", "Multivariable")]
  data[, Category_other := factor(Category_other, levels = c("Univariable", "Multivariable"))]
  data[,name := gsub("BMI_UKB|GIANT_2015_BMI", "BMI", exposure) %>% gsub("WC_UKB|GIANT_2015_WC", "WC", .)]
  
 ggforestplot::forestplot(
    df = data,
    name = name,
    se = se,
    estimate = b,
    pvalue = pval,
    psignif = 0.05,
    xlab = "Effect size (SD) per 1-SD\nincrease in adiposity",
    ci = 0.95,
    logodds = TRUE,
    colour = method,
    xlim = data[method != "MR Egger", round(c(min(exp(lci)),max(exp(uci))), digits = 1)]
  ) + theme(legend.position = "right") +
    theme(text = element_text(size = 10)) +
    scale_x_continuous(breaks = seq(from = round(min(exp(data$lci)), digits = 1), to = round(max(exp(data$uci)), digits = 1), by = 0.2)) +
    theme(legend.position="right") +
    ggforce::facet_col(
      # facets = ~outcome_category,
      facets = ~Category_other,
      scales = "free_y",
      space = "free"
    )
  
  ggsave(paste0("Results/", file_name[i],"colour", ".png"),
         width=483/72,height=250/72, units="in", scale=1,
         device = "png")
  #### figure 2
  
  data[Category_other == "Multivariable", name := name %>% ifelse(. == "BMI", "BMI adjusted for WC", .) %>% ifelse(. == "WC", "WC adjusted for BMI", .)]
  doA <- data
  resultsA <- data.frame(variable = as.character(1:nrow(doA)),
                         estimate = round(doA$b, digits =2),
                         lci =  round(doA$lci, digits = 2),
                         uci =  round(doA$uci, digits = 2),
                         n = doA$nsnp,
                         P_value = formatC(doA$pval, format = "e", digits = 1))
  
  
  mylabels <- data.frame(heading1 = doA$name,
                         heading2 = doA$method,
                         heading3 = as.character(NA),
                         variable = as.character(1:nrow(doA)))
  
  ckbplotr::make_forest_plot(panels = list(resultsA),
                   col.key = "variable",
                   row.labels = mylabels,
                   exponentiate = TRUE,
                   pointsize = 2, 
                   rows = unique(mylabels$heading1),
                   col.stderr = NULL,
                   col.lci = "lci",
                   col.uci = "uci",
                   col.left         = c("n"),
                   col.left.heading = c("n SNP"),
                   col.right = "P_value",
                   col.right.heading = c("Effect (95% CI)", "P-value"),
                   xlab = c("Effect of 1 SD increase in WC/BMI on NAFLD (OR)"),
                   blankrows = c(0,0,0,0),
                   col.right.hjust = 1,
                   panel.headings = NULL,
                   scalepoints = FALSE)
 
 ggsave(paste0("Results/", file_name[i],"Clean", ".png"),
        width=483/72,height=250/72, units="in", scale=1,
        device = "png")
}

#Figure 3
harm_class <- harm_class[class != "WHRonly-"] #there is no value at doing the analysis on WHRonly-
res_class <- separate(res_class, col = "exposure", into = c("class", "exposure"), sep = "_", remove =  TRUE) %>% as.data.table(.)
harm_class[,exposure := exposure %>% gsub(".*_", "", .) ]
res_class[, id.exposure := exposure]
harm_class[, id.exposure := exposure]

mr_results <-  res_class[exposure == "ukb-b-19953" & outcome == "NAFLD",]
dat <- harm_class[exposure == "ukb-b-19953" & outcome == "NAFLD",]

mr_results[,align := class]
dat[,align := class]
dat[, Locus := ""]
source("/mnt/sde/gagelo01/Projects/small_MR_exploration/FI_BMI/Analysis/my_mr_scatter_plot.R")
debug(my_mr_scatter_plot)
k <- my_mr_scatter_plot( dat = dat, mr_results = mr_results)
k + xlab("SNP effect on BMI") +ylab("SNP effect on NAFLD")

ggsave(paste0("Results/", "Figure3_possibly", ".png"),
       width=657/72,height=385/72, units="in", scale=1,
       device = "png")


#Figure 4
exp_name <- list(c("WC_UKB", "Fat_Liver"), c("GIANT_2015_WC", "Fat_Liver"))
study<-c("ukb", "giant")
file_name <- c("Figure3", "Supplementary_figure3")

for(i in 1:length(exp_name)) {
  
  doA<-res_univariate[outcome == "Mahajan_Type2diabetes" & exposure %in% exp_name[[i]]]
  
  
  format_toforest <- function(do, multivariate) {
    do[,heading := "Univariable"]
    multivariate <- multivariate[method != "Multivariable Egger Intercept"]
    multivariate[,heading := "Multivariable"]
    multivariate <- multivariate[clump_exposure == "Fat_Liver", ][order(exposure)]
    MVMR<- rbindlist(list(do, multivariate), use.names = TRUE, fill = TRUE)
    MVMR[, name := exposure %>% gsub("WC_UKB|GIANT_2015_WC", "WC", . ) %>% gsub("Fat_Liver", "Liver fat", .)]
    return(MVMR)
  }
  
  MVMRA <- format_toforest(do = doA, multivariate = resmvmr[[paste0(exp_name[[i]][1],"-and-Fat_Liver-on-Mahajan_Type2diabetes")]])
  MVMRA[heading == "Multivariable", name := name %>% ifelse(. == "Liver fat", "Liver fat adjusted for WC", .) %>% ifelse(. == "WC", "WC adjusted for Liver fat", .)]
  
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
                   col.left.heading = c("n SNP"),
                   xlab = c("Effect of 1-SD increase in Liver fat/WC on NAFLD (OR)"),
                   colour           = "colour",
                   blankrows = c(0,0,0,0),
                   xticks = round(c(1, 2), digits = 0))
  
  
  ggsave(paste0("Results/", file_name[i], "Complete.png"),
         width=798/72,height=583/72, units="in", scale=1,
         device = "png")
}
