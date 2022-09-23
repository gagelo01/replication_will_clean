#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref <-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/")]
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index[category %in% c("Trait", "Disease"), ][group_name != "Dietary habits",][trait != "Chronotype_Morningness", ]


ID_mrbase_exp <- c("ieu-a-7") 
exp_mrbase <- paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_exp, "/", ID_mrbase_exp, ".vcf.gz")
ID_server_exp <- c("dis-9-1", "dis-8-1", "dis-6-1")
exp_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")


inst <- GagnonMR::get_inst("/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-2-1/dis-2-1.vcf.gz", pval = 5e-6, r2 = 0.001)

options(future.globals.maxSize= 5e9)
plan(multisession, workers = 4, gc = TRUE) #I should try using multicore

out <- future_map(as.list(c(exp_mrbase,exp_server)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr <- gwasvcf::query_gwas(vcf = x, rsid = inst[,unique(SNP)], proxies = "yes", bfile = ldref, tag_r2 = 0.8)
  outr <- outr %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

out[outcome %in% "ieu-a-7", outcome := "Nikpay_CAD"]
harm <- TwoSampleMR::harmonise_data(inst, out, action = 1)
setDT(harm)
harm[,exposure_outcome := paste0(exposure, "_", outcome)]
harm_split <- split(harm, harm$exposure_outcome)
res <- future_map(harm_split, function(x) {GagnonMR::all_mr_methods(x)}, .options = furrr_options(seed = TRUE)) %>%
  rbindlist(., fill = TRUE)
egger_intercept <- TwoSampleMR::mr_pleiotropy_test(harm)
setDT(egger_intercept)
egger_intercept[, c("id.exposure", "id.outcome") := NULL]
res[, c("id.exposure", "id.outcome") := NULL]
fwrite(res, "Data/Modified/NAFLD_exposure_p5e-6.txt")
fwrite(egger_intercept, "Data/Modified/NAFLD_egger_intercept_p5e-6.txt")
message("This script finished without errors")
