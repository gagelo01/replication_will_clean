#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref <-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

ao <- fread("/mnt/sdf/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao[id %in% list.files("/mnt/sdf/gagelo01/Vcffile/MRBase_vcf/")]
df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
df_index[category %in% c("Trait", "Disease"), ][group_name != "Dietary habits",][trait != "Chronotype_Morningness", ]


ID_mrbase_exp <- c("ieu-a-61", "ieu-a-835", "ukb-b-19953", "ukb-b-9405","ieu-a-7", "ieu-a-79" ) 
exp_mrbase <- paste0("/mnt/sdf/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_exp, "/", ID_mrbase_exp, ".vcf.gz")
ID_server_exp <- c("dis-2-1", "dis-9-1", "dis-8-1", "dis-6-1", "trait-14-6", "trait-14-7", "trait-14-8")
exp_server <- paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")


options(future.globals.maxSize= 5e9)
plan(multisession, workers = 9, gc = TRUE) #I should try using multicore

inst <- future_map(as.list(c(exp_mrbase,exp_server)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  GagnonMR::get_inst(x)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


out <- future_map(as.list(c(exp_mrbase,exp_server)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr <- gwasvcf::query_gwas(vcf = x, rsid = inst[,unique(SNP)], proxies = "yes", bfile = ldref, tag_r2 = 0.8)
  outr <- outr %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


out$outcome %>% unique
convert <- data.table(id = c("ieu-a-61", "ieu-a-835", "UKB-b-19953", "UKB-b-9405","ieu-a-7", "ieu-a-79"),
                     trait = c("GIANT_2015_WC", "GIANT_2015_BMI", "BMI_UKB", "WC_UKB", "Nikpay_CAD", "WHRadjBMI")) 

inst <- merge(inst, convert, by.x = "exposure", by.y = "id", all.x = TRUE)
inst[,exposure := ifelse(is.na(trait), exposure, trait)]
inst[,trait := NULL]
inst$exposure %>% unique
out <- merge(out, convert, by.x = "outcome", by.y = "id", all.x = TRUE)
out[,outcome := ifelse(is.na(trait), outcome, trait)]
out[,trait := NULL]
out$outcome %>% unique
exp <- TwoSampleMR::convert_outcome_to_exposure(out)

#save
fwrite(exp, "Data/Modified/all_inst_mvmr")
fwrite(out, "Data/Modified/all_outcome_mvm.txt" )
fwrite(inst, "Data/Modified/inst_all_sign_clump.txt")

message("This script finished without errors")
