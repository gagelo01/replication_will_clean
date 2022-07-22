#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref <-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

instclass <- readxl::read_excel("Data/Raw/41467_2018_4124_MOESM4_ESM.xlsx", skip = 3) %>% as.data.table(.)
arguments <- instclass[!is.na(rsid),.(rsid, Class)]
ID_mrbase_exp <- c("ukb-b-19953", "ieu-a-835") #BMI ukb, BMI GIANT

arguments <- tidyr::crossing(arguments, data.frame(id = ID_mrbase_exp)) %>% as.data.table(.)
arguments[, class_id := paste0(Class, "_", id)]
ldmat <- ieugwasr::ld_matrix_local(unique(arguments$rsid), plink_bin = genetics.binaRies::get_plink_binary(), 
                                   bfile = ldref, with_alleles = FALSE)
arguments <- arguments[rsid %in% rownames(ldmat), ]
arguments_split <- split(arguments, arguments$class_id)


inst_class <- map(arguments_split, function(x) {
tsmr <- gwasvcf::query_gwas(paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", x[1,]$id, "/", x[1,]$id, ".vcf.gz"), rsid =  x$rsid, proxies = "yes", tag_r2 = 0.8, bfile = ldref) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(.) %>% 
  as.data.table(.)
tsmr[,exposure := x[1,]$id]
tsmr[,id.exposure := exposure]
tsmr[, class := x[1,]$Class]
return(tsmr)}) %>% rbindlist(., fill = TRUE)

inst_class[,  exposure := exposure %>% gsub("ukb-b-19953", "BMI_UKB", . ) %>% gsub("ieu-a-835", "GIANT_2015_BMI", . )]
inst_class[,class_exposure := paste0(class, "_", exposure)]

ID_server_exp <- c("dis-2-1", "trait-14-8")
exp_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")

out <- future_map(as.list(c(exp_server)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr <- gwasvcf::query_gwas(vcf = x, rsid = inst_class[,unique(SNP)], proxies = "yes", bfile = ldref, tag_r2 = 0.8)
  outr <- outr %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table(.)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

out[,id.outcome := outcome]
inst_class[,exposure := class_exposure]
inst_class[,id.exposure := exposure]

harm <- TwoSampleMR::harmonise_data(exposure_dat = inst_class, outcome_dat = out, action = 1)
setDT(harm)
harm[, exposure_outcome := paste0(exposure, "_", outcome)]
fwrite(harm, "Data/Modified/harm_class.txt")

message("This script finished without errors")