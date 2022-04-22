#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")

code <- fread("/home/couchr02/Mendel_Commun/Christian/GWAS/IRM_Liu/GCST_code.tsv")
code[,trait_organ := paste0(trait, "_", organ)]
traduction = fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]

newrow <- data.table(id = paste0("trait-14-", 1:nrow(code)), trait = code[,trait_organ], group_name = "public",
                     year = 2021,author = "Liu Yi",consortium = "UKBiobank",sex = "Males and Females",population = "European",
                     initial_build = "HG19/GRCh37", unit = "SD", nsnp = NA, sample_size = 38000,
                     category = "Trait", pmid = 34128465, ncase = NA,
                     sd = 1, note = NA, ncontrol = NA)

df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
df_index <- rbind(df_index, newrow)

code <- merge( code,newrow[,.(trait, id)],  by.x = "trait_organ", by.y = "trait")
code_list<-split(code, 1:nrow(code))
df_index_copy <- df_index
traduction_copy <- traduction

options(future.globals.maxSize= 1e10)
plan(multicore, workers = 6)

formatvcf_wrapper <- function(codevec, df_index, traduction) {
data <- fread(paste0("/home/couchr02/Mendel_Commun/Christian/GWAS/IRM_Liu/",codevec$GCST_code,"_buildGRCh37.tsv.gz"))
GagnonMR::formattovcf_createindex(all_out = data,
                                  snp_col = "variant_id",
                                  outcome_name = codevec$trait_organ,
                                  beta_col = "beta",
                                  se_col = "standard_error",
                                  pval_col = "p_value",
                                  eaf_col = "effect_allele_frequency",
                                  effect_allele_col = "effect_allele",
                                  other_allele_col =  "other_allele",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = 38000,
                                  chr_col = "chromosome",
                                  pos_col = "base_pair_location",
                                  units = "SD",
                                  traduction = traduction,
                                  out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = "public",
                                  year = 2021,
                                  author = "Liu Yi",
                                  consortium = "UKBiobank",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 34128465,
                                  note = NA,
                                  should_create_id = FALSE,
                                  ID = codevec$id )
}


future_map(code_list, function(x) {formatvcf_wrapper( codevec = x, df_index = df_index_copy, traduction = traduction_copy)},
           .options = furrr_options(seed = TRUE))

df_index[pmid == 34128465 & grepl("Volume", trait) & trait != "Volume_Pancreas", sample_size := 32860]
df_index[pmid == 34128465 & trait == "Volume_Pancreas", sample_size := 31758]
df_index[pmid == 34128465 & trait %in% c("Fat_Pancreas", "Iron_Pancreas"), sample_size := 25617]
df_index[pmid == 34128465 & trait %in% c("Fat_Liver", "Iron_Liver"), sample_size := 32858]

fwrite(df_index, "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")

message("This script finished without error")



