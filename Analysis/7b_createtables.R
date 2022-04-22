#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library("xlsx")
library(writexl)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
inst_all_sign_clump  <- inst_all_sign_clump[!(exposure %in% c("Volume_VAT", "Volume_ASAT")), ]
resmvmr <- readRDS( "Data/Modified/res_mvmr.rds")
res_univariate <- fread("Data/Modified/res_univariate.txt")
res_univariate <- res_univariate[!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
FandQ <- fread( "Data/Modified/FandQ_univariate.txt")
FandQ  <- FandQ [!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
egger_intercept <- fread("Data/Modified/egger_intercept.txt")
egger_intercept <- egger_intercept[!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
harm_univariate <- harm_univariate[!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
res_class <- fread("Data/Modified/res_class.txt")

#supptable 1 description of cohorts
ao <- fread("/mnt/sdf/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ID_mrbase_exp <- c("ieu-a-61", "ieu-a-835", "ukb-b-19953", "ukb-b-9405","ieu-a-7","ieu-a-79" ) 
ao[id %in% ID_mrbase_exp,]

df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
ID_server_exp <- c("dis-2-1", "dis-9-1", "dis-8-1", "dis-6-1", "trait-14-8")
df_index[id %in% ID_server_exp,]

dataset <- rbindlist(list(ao[id %in% ID_mrbase_exp,], df_index[id %in% ID_server_exp,]), fill = TRUE)
dataset <- dataset[,.(trait,group_name, year, author,  consortium, sex, population, unit, nsnp, sample_size,ncase, ncontrol,pmid)]
dataset[, trait := trait %>% gsub("Waist-to-hip ratio", "Waist-to-hip ratio adjusted for BMI", . )]
dattrait <- data.frame(trait = dataset$trait,
url = c("https://gwas.mrcieu.ac.uk/files/ukb-b-19953/ukb-b-19953.vcf.gz",
  "https://gwas.mrcieu.ac.uk/files/ukb-b-9405/ukb-b-9405.vcf.gz",
  "https://gwas.mrcieu.ac.uk/files/ieu-a-7/ieu-a-7.vcf.gz",
  "https://gwas.mrcieu.ac.uk/files/ieu-a-835/ieu-a-835.vcf.gz",
  "https://gwas.mrcieu.ac.uk/files/ieu-a-61/ieu-a-61.vcf.gz",
  "https://gwas.mrcieu.ac.uk/files/ieu-a-79/ieu-a-79.vcf.gz",
  "https://www.ebi.ac.uk/gwas/publications/34841290",
  "http://diagram-consortium.org/downloads.html",
  "https://data.mendeley.com/datasets/gbbsrpx6bs/1",
  "http://diagram-consortium.org/downloads.html",
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016676/"))
  
dataset <- merge(dataset, dattrait, "trait")  
#supp table 
for(i in 1:length(resmvmr)) {
  resmvmr[[i]][, exposure_outcome_analyzed := names(resmvmr)[i]]
}

dt_resmvmr <- rbindlist(resmvmr, fill = TRUE)
dt_resmvmr[,F_stastistics := NULL]
inst_all_sign_clump <- inst_all_sign_clump[,.(exposure, SNP, chr.exposure, pos.exposure, other_allele.exposure,
                       effect_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure)]

########Supplementary table titles and description
dt_title <- data.table(title = paste0("Supplementary Table ", 1:7),
                       caption = c("Description of the datasets used.",
                                   "Harmonised data sets.",
                                   "Univariable Mendelian Randomization results",
                                   "Univariable Mendelian Randomization Egger's intercept",
                                   "Multivariable Mendelian randomization results.",
                                   "Class specific Mendelian randomization results",
                                   "Instrument strength for univariable MR"))


writexl::write_xlsx(x = list("Tables captions and titles" = dt_title,
                            "Supplementary Table 1" = dataset,
                            "Supplementary Table 2" = harm_univariate,
                            "Supplementary Table 3" = res_univariate, 
                            "Supplementary Table 4" = egger_intercept, 
                            "Supplementary Table 5" = dt_resmvmr, 
                            "Supplementary Table 6" = res_class,
                            "Supplementary Table 7" = FandQ),
                    path = "Results/supplementary_tables_clean.xlsx")


