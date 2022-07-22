#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(MendelianRandomization)

# objects <- c("res_univariate", "res_bmiwcukb_nafld_eloi", "res_wcnafldukb_t2d_eloi", "res_wcnafldukb_cad_eloi",
#              "res_bmiwcgiant_nafld_eloi", "res_wcnafldgiant_t2d_eloi", "res_wcnafldgiant_cad_eloi")
# 
# for(i in 1:length(objects)) {
#   assign(objects[i] , fread(paste0("Results/", objects[i], ".txt"))) 
# }
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
FandQ <- fread( "Data/Modified/FandQ_univariate.txt")
res_class <- fread("Data/Modified/res_class.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
harm_class <- fread("Data/Modified/harm_class.txt")

#Abstract
return_format_data<-function(data) {
return(data[, paste0(round(exp(b), digits = 2), " 95% CI=", round(exp(lci), digits = 2), "-",  round(exp(uci), digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}

return_format_data(resmvmr$`WC_UKB-and-BMI_UKB-on-NAFLD`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "none",])
return_format_data(resmvmr$`WC_UKB-and-BMI_UKB-on-NAFLD`[method == "Multivariable IVW" & exposure == "BMI_UKB" & clump_exposure == "none",])
return_format_data(resmvmr$`WC_UKB-and-Fat_Liver-on-Mahajan_Type2diabetes`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "Fat_Liver",])
return_format_data(resmvmr$`WC_UKB-and-Fat_Liver-on-van_der_Harst_CAD`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "Fat_Liver",])

#Results
res_univariate[exposure == "WC_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "WC_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[exposure == "WC_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",])
res_univariate[exposure == "BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[exposure == "BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",])

#para 2
res_univariate[exposure == "whradjbmi_ukbonly" & outcome == "NAFLD" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "whradjbmi_ukbonly" & outcome == "NAFLD" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[exposure == "whradjbmi_ukbonly" & outcome == "NAFLD" & method == "Inverse variance weighted",])


#std of BMI is 4.77168 in BMI points; std of wc is 13.4238 in cm
# res_univariate[exposure == "WC_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",
#                c(exp(c(b, lci, uci)/13.4238), pval)]
# 
# res_univariate[exposure == "BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted",
#                c(exp(c(b, lci, uci)/4.77168), pval)]

sum(inst_all_sign_clump[exposure == "BMI_UKB", unique(SNP)] %in% inst_all_sign_clump[exposure == "WC_UKB", unique(SNP)])
return_format_data(resmvmr$`WC_UKB-and-BMI_UKB-on-NAFLD`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "none",])
return_format_data(resmvmr$`WC_UKB-and-BMI_UKB-on-NAFLD`[method == "Multivariable IVW" & exposure == "BMI_UKB" & clump_exposure == "none",])

resmvmr$`WC_UKB-and-BMI_UKB-on-NAFLD`[clump_exposure == "none",.(exposure, F_stastistics)]

#para 3
res_class[exposure == "BMI+WHR+_BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted", ] %>% 
  return_format_data()
res_class[exposure == "BMI+WHR-_BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted", ] %>% 
  return_format_data()
res_class[exposure == "BMIonly+_BMI_UKB" & outcome == "NAFLD" & method == "Inverse variance weighted", ] %>% 
  return_format_data()


res_univariate[exposure == "WC_UKB" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "WC_UKB" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "WC_UKB" & outcome == "Mahajan_Type2diabetes",])
return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "WC_UKB" & outcome == "van_der_Harst_CAD",])

res_univariate[exposure == "NAFLD" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "NAFLD" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "NAFLD" & outcome == "Mahajan_Type2diabetes",])

res_univariate[exposure == "Fat_Liver" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted",nsnp]
FandQ[exposure == "Fat_Liver" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted", ]
return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "Fat_Liver" & outcome == "Mahajan_Type2diabetes",])
return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "Fat_Liver" & outcome == "van_der_Harst_CAD",])


return_format_data(resmvmr$`WC_UKB-and-Fat_Liver-on-Mahajan_Type2diabetes`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "Fat_Liver",])
return_format_data(resmvmr$`WC_UKB-and-Fat_Liver-on-Mahajan_Type2diabetes`[method == "Multivariable IVW" & exposure == "Fat_Liver" & clump_exposure == "Fat_Liver",])

# res_univariate[exposure == "WC_UKB" & outcome == "van_der_Harst_CAD" & method == "Inverse variance weighted",nsnp]
# FandQ[exposure == "WC_UKB" & outcome == "van_der_Harst_CAD" & method == "Inverse variance weighted", ]
# return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "WC_UKB" & outcome == "van_der_Harst_CAD",])
# res_univariate[exposure == "WC_UKB" & outcome == "van_der_Harst_CAD" & method == "Inverse variance weighted",nsnp]
# FandQ[exposure == "NALFD" & outcome == "van_der_Harst_CAD" & method == "Inverse variance weighted", ]
# return_format_data(res_univariate[method == "Inverse variance weighted" & exposure == "NAFLD" & outcome == "van_der_Harst_CAD",])
# return_format_data(resmvmr$`WC_UKB-and-NAFLD-on-van_der_Harst_CAD`[method == "Multivariable IVM" & exposure == "WC_UKB" & clump_exposure == "NAFLD",])
# return_format_data(resmvmr$`WC_UKB-and-NAFLD-on-van_der_Harst_CAD`[method == "Multivariable IVM" & exposure == "NAFLD" & clump_exposure == "NAFLD",])

#mediation
thetad <- resmvmr$`WC_UKB-and-Fat_Liver-on-Mahajan_Type2diabetes`[method == "Multivariable IVW" & exposure == "WC_UKB" & clump_exposure == "Fat_Liver",]$b
thetat <- res_univariate[exposure == "WC_UKB" & outcome == "Mahajan_Type2diabetes" & method == "Inverse variance weighted", b]
1 - (thetad/thetat)

thetad/resmvmr$`WC_UKB-and-Fat_Liver-on-Mahajan_Type2diabetes`[method == "Multivariable IVW" & exposure == "Fat_Liver" & clump_exposure == "Fat_Liver",]$b

##########Discussion##############
res_univariate[exposure == "NAFLD" & outcome ==  "Mahajan_Type2diabetes",][method != "MR Egger",][pval < 0.05, .N]

###########Method###################
to_exclude = c("APOE", "ABO", "HLA-A")
window = 2e+06
gencode <- fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
list <- vector(mode = "list", length = length(to_exclude))
for (i in 1:length(to_exclude)) {
  bon <- gencode[gene_name == to_exclude[i], ]
  list[[i]] <- data.table(chr = bon[1, ]$chr, start = min(bon$start) - 
                            window/2, end = max(bon$end) + window/2, gene_name = bon[1, 
                            ]$gene_name)

  }
list

message("This script finished without errors")
