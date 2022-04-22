#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")

t2d <- fread("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean/Data/Raw/METAANALYSIS_DIAGRAM_SE1.txt")
t2d <- separate(t2d, "Chr:Position", sep = ":", into = c("Chr", "Position"))
t2d[, Chr := as.integer(Chr)]
t2d[, Position := as.integer(Position)]

traduction = fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]

debugonce(GagnonMR::formattovcf_createindex)
GagnonMR::formattovcf_createindex(all_out = t2d,
                                  snp_col = NULL,
                                  outcome_name = "Scott_Type2Diabetes",
                                  beta_col = "Effect",
                                  se_col = "StdErr",
                                  pval_col = "P-value",
                                  eaf_col = NULL,
                                  effect_allele_col = "Allele1",
                                  other_allele_col =  "Allele2",
                                  ncase_col = 26676,
                                  ncontrol_col = 132532,
                                  samplesize_col = "TotalSampleSize",
                                  chr_col = "Chr",
                                  pos_col = "Position",
                                  units = "log odds",
                                  traduction = traduction,
                                  out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf",
                                  df_index = fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt"),
                                  group_name = "public",
                                  year = 2017,
                                  author = "Scott RA",
                                  consortium = "DIAGRAMM",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Disease",
                                  pmid = 28566273,
                                  note = NA,
                                  should_create_id = TRUE)


message("This script finished without error")
