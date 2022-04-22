#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref <-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

harm <- fread("Data/Modified/harm_class.txt")
harm <- harm[class != "WHRonly-"] #there is no value at doing the analysis on WHRonly-
all_mr_methods_safely <- safely(GagnonMR::all_mr_methods)


options(future.globals.maxSize= 5e9)
plan(multisession, workers = 6, gc = TRUE) #I should try using multicore

res <- future_map(split(harm, harm$exposure_outcome),function(x) {all_mr_methods(x)}, .options = furrr_options(seed = TRUE)) %>%
  rbindlist(.,fill = TRUE)

fwrite(res, "Data/Modified/res_class.txt")
message("This script finished without errors")