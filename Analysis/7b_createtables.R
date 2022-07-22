#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library("xlsx")
library(writexl)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
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
NAFLDsensi <-fread("Data/Modified/NAFLD_exposure_p5e-6.txt")
NAFLDsensi[,c("id.exposure", "id.outcome") := NULL]
egger_intercept5e6 <- fread("Data/Modified/NAFLD_egger_intercept_p5e-6.txt")

#supptable 1 description of cohorts
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ID_mrbase_exp <- c("ieu-a-61", "ieu-a-835", "ukb-b-19953", "ukb-b-9405","ieu-a-7","ieu-a-79" ) 

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ID_server_exp <- c("dis-2-1", "dis-9-1", "dis-8-1", "dis-6-1", "trait-14-8", "trait-10-2")

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
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016676/",
  "https://zenodo.org/record/3341464#.Ytgqs3bMKUk"),
download_date = c(rep("2021-16-06", 6), "2021-11-03", rep("2021-16-06", 3), "2022-04-15", "2022-07-20"))
  
dataset <- merge(dataset, dattrait, "trait")  
#
harm_univariate[,c("id.outcome", "pval_origin.outcome", "action", "pval_origin.exposure", "id.exposure", "exposure_outcome", "mr_keep") := NULL]
#
egger_intercept[,c("id.exposure","id.outcome") := NULL] 
#supp table 
for(i in 1:length(resmvmr)) {
  resmvmr[[i]][, exposure_outcome_analyzed := names(resmvmr)[i]]
}


dt_resmvmr <- rbindlist(resmvmr, fill = TRUE)
inst_all_sign_clump <- inst_all_sign_clump[,.(exposure, SNP, chr.exposure, pos.exposure, other_allele.exposure,
                       effect_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure)]

convert <- data.table(id = c("ieu-a-61", "ieu-a-835", "UKB-b-19953", "UKB-b-9405","ieu-a-7", "ieu-a-79"),
                      trait = c("GIANT_2015_WC", "GIANT_2015_BMI", "BMI_UKB", "WC_UKB", "Nikpay_CAD", "WHRadjBMI")) 
egger_intercept5e6 <- merge(egger_intercept5e6, convert, by.x= "outcome", by.y = "id", all.x = TRUE)
egger_intercept5e6[!is.na(trait), outcome := trait]
egger_intercept5e6[, trait := NULL]
########Supplementary table titles and description
dt_title <- data.table(title = paste0("Supplementary Table ", 1:9),
                       caption = c("Description of the datasets used.",
                                   "Harmonised data sets.",
                                   "Univariable Mendelian Randomization results",
                                   "Univariable Mendelian Randomization Egger's intercept",
                                   "Multivariable Mendelian randomization results.",
                                   "Group specific Mendelian randomization results",
                                   "NAFLD effect on T2D and CAD using pvalue threshold of 5e-6 and LD clump of R2<0.001",
                                   "NAFLD egger intercept on T2D and CAD using pvalue threshold of 5e-6 and LD clump of R2<0.001",
                                   "Instrument strength and heterogeneity statistics for univariable MR"),
                       column_legend = c("trait = UNique trait identifier; group_name = the name of the study group;  year = year data was published
; author = author of the data; consortium = conortium for the sample;  
sex = sex of the sample (in this study always Males and Females); 
population = Ancestry of the sample; 
unit = The unit either standard deviation (SD) or log Odds ratio (log(orr)); nsnp = number of SNPs in the summary statistics;
sample_size = The maximum sample size of the meta-analysis; ncase = number of cases; ncontrol = number of controls; pmid = Pubmed ID",
                                         "SNP = the rsid;  effect_allele.exposure = effect allele of the exposure; other_allele.exposure = non effect allele of the exposure
effect_allele.outcome = effect allele of the outcome; other_allele.outcome = non effect allele of the outcome; 
beta.exposure =  SNP effect on exposure; beta.outcome = SNP effect on outcome; eaf.exposure = effect allele frequency in the exposure;           
eaf.outcome = effect allele frequency in the outcome; remove = should you remove the SNP; palindromic = is the SNP palindromic;
ambiguous = is the SNP A/T or C/G; outcome = Unique identifier for the outcome; chr.outcome = chromosome of the SNP integer from 1 to 22;
pos.outcome = position on the chromosome (GRCH37); se.outcome = standard error of the SNP outcome association; pval.outcome = p-value of the SNP outcome association;
samplesize.outcome = Maximum sample size of the outcome GWAS; ncase.outcome = maximum nmber of case of the outcome GWAS;       
ncontrol.outcome = Maximum nmber of controls of the outcome GWAS; mr_keep.outcome = Should you remove this genetic instruments because the allele frequency does not match;
exposure = Unique name for exposure;  chr.exposure = SNP Chromosome; pos.exposure =  SNP position on the chromosome (GRCH37);
se.exposure = standard error of the SNP exposure association; pval.exposure = p-value of the SNP exposure association;        
samplesize.exposure = Maximum sample size of the exposure GWAS; ncase.exposure = maximum nmber of case of the exposure GWAS;       
ncontrol.exposure = Maximum nmber of controls of the exposure GWAS; mr_keep.exposure = Should you remove this genetic instruments because the allele frequency does not match;",
                                         "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome;",
"nsnp = The number of SNPs genetic instruments used to compute the estimate; b = the estimate scaled per SD or log(orr); 
se =  standard error of the estimate; pval = the p-value of the estimate; lci = the 95% confidence intervall lower bound of the effect;
uci = the 95% confidence intervall upper bound of the effect; type_of_test = categorisation of the method based on Slob and Burgess pmid = 32249995",
"outcome = unique name for the outcome; exposure = unique name for the exposure; egger_intercept = the egger intecept estimate;
se = the egger intecept standard error; pval = the egger intecept p-value", 
"exposure = unique name for the exposure; outcome = unique name for the outcome; b = the direct effect (adjusted for the other exposures);
se = the standard error of the direct effect; lci = The 95% confidence intervall lower bound of the direct effect; 
uci = The 95% confidence intervall upper bound of the direct effect; pval = The p-value of the direct effect;
cochranQ = The cochran's Q statitics; cochranQpval =The cochran's Q statitics  p-value; nsnp = the number of SNPs (genetic instruments);
method = The name of the multivarialble MR method used;   F_stastistics = The conditionnal F-statistics; 
clump_exposure = The clumping was based on the highest pvalue of which exposure. None = both exposures;         
exposure_outcome_analyzed = The name of the exposures and the outcomes analysed in the form paste0(exposure1,'-and-',exposure2,'-on-',outcome)",
"Same as Supplementary Table 3",
"Same as Supplementary Table 4",
"outcome = Unique name for the outcome; exposure = Unique name for the exposure; method = the method to compute the cochran's Q test;
Q = cochran's Q value; Q_df = The cochran's Q statistics number of degree of freedom; Q_pval = The Cochran's Q  p-value;
fstat = the F-statistics of the instrumental variables; rsq = the variance explained of the exposure by the instrumental variables."))


writexl::write_xlsx(x = list("Tables captions and titles" = dt_title,
                            "Supplementary Table 1" = dataset,
                            "Supplementary Table 2" = harm_univariate,
                            "Supplementary Table 3" = res_univariate, 
                            "Supplementary Table 4" = egger_intercept, 
                            "Supplementary Table 5" = dt_resmvmr, 
                            "Supplementary Table 6" = res_class,
                            "Supplementary Table 7" = NAFLDsensi,
                            "Supplementary Table 8" = egger_intercept5e6,
                            "Supplementary Table 9" = FandQ),
                    path = "Results/supplementary_tables_clean.xlsx")


message("This script finished without errors")
