#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(MendelianRandomization)
library(writexl)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/replication_will_clean")
# all_inst_proxies <- fread( "Data/Modified/all_inst_proxies.txt")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr")
all_inst_mvmr[, id.exposure := exposure]
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvm.txt" )
all_outcome_mvmr [, id.outcome := outcome]
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
inst_all_sign_clump[, id.exposure := exposure]

#####Liu
# liu <- c("rs738409", "rs58542926")
# LDlinkR::LDproxy_batch(snp = liu, token = Sys.getenv("LDLINK_TOKEN"), 
#                        append = TRUE)
# liu_LD <- fread("combined_query_snp_list.txt")
# rs_prox <- liu_LD[RS_Number %in% inst_all_sign_clump[exposure == "NAFLD",]$SNP, ]$RS_Number
# inst_NAFLD <- inst_all_sign_clump[exposure == "NAFLD",][pval.exposure < 5e-8,]
# index_to_remove <- inst_all_sign_clump[exposure == "NAFLD"]$SNP[!(inst_all_sign_clump[exposure == "NAFLD"]$SNP %in% inst_NAFLD[pval.exposure < 5e-8,]$SNP)]
# inst_all_sign_clump <- inst_all_sign_clump[!(SNP %in% index_to_remove)]

###ivw-mr on everything
# exposure_name_df <- data.frame(exposure_name = inst_all_sign_clump$exposure %>% unique)
# outcome_name_df <- data.frame(outcome_name = all_outcome_mvmr$outcome %>% unique)
# arguments_df <- tidyr::crossing(exposure_name_df, outcome_name_df)
# harm_wrapper<-function(exposure_name, outcome_name) {
#   dat <- harmonise_data(inst_all_sign_clump[exposure == exposure_name], all_outcome_mvmr[outcome == outcome_name], action = 2)
# return(dat)
# }
# list_harm <- pmap(arguments_df, harm_wrapper)
# steiger_filtering_safely <- safely(steiger_filtering)
# list_harm_s <- lapply(list_harm, function(x) steiger_filtering_safely)
# index <- sapply(list_harm_s, function(x) is.null(x$error))
# list_harm[!index]
# harm_ivw <- rbindlist(list_harm, fill = TRUE)
# 
# 
# res_ivw <- mr(dat = harm_ivw, method_list = "mr_ivw")
# setDT(res_ivw)
# res_ivw <- res_ivw[outcome != exposure]
# res_ivw[exposure == "NAFLD",][pval < 0.05]$outcome
# res_ivw[exposure == "Mahajan_Type2diabetes",][pval < 0.05]$outcome
# res_ivw[exposure == "scott_t2d",][pval < 0.05]$outcome
###univariate

# inst_all_sign_clump <- inst_all_sign_clump[exposure %in% c("BMI_UKB", "WC_UKB", "NAFLD", "Mahajan_Type2diabetes",  "GIANT_2015_WC",  "GIANT_2015_BMI") ]
# all_outcome_mvmr <- all_outcome_mvmr[outcome %in% c("Mahajan_Type2diabetes", "BMI_UKB", "WC_UKB", "NAFLD", "van_der_Harst_CAD", "GIANT_2015_WC",  "GIANT_2015_BMI"), ]

inst_all_sign_clump_split <- split(inst_all_sign_clump, inst_all_sign_clump$exposure)
harm_univariate <- lapply(inst_all_sign_clump_split, function(x) {harmonise_data(x, all_outcome_mvmr, action = 1)}) %>% rbindlist(.,fill = TRUE)


setDT(harm_univariate)
harm_univariate <- harm_univariate[!(exposure == outcome), ]
harm_univariate <- harm_univariate[!(exposure %in% c("BMI_UKB", "WC_UKB", "GIANT_2015_WC", "GIANT_2015_BMI") & outcome %in% c("BMI_UKB", "WC_UKB", "GIANT_2015_WC", "GIANT_2015_BMI"))]
harm_univariate[, exposure_outcome := paste0(exposure,"_", outcome)]
harm_univariate <- harm_univariate[!(exposure == outcome),]
harm_univariate <-  harm_univariate[exposure != "Volume_ASAT",]
harm_univariate <- harm_univariate[!(exposure_outcome %in% c("NAFLD_GIANT_2015_BMI", "NAFLD_GIANT_2015_WC", "NAFLD_WHRadjBMI")), ]
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
egger_intercept <- mr_pleiotropy_test(harm_univariate[!(exposure == "NAFLD" & outcome %in% c("GIANT_2015_BMI", "GIANT_2015_WC")),])
list_harm_univariate <- split(harm_univariate, harm_univariate$exposure_outcome)


options(future.globals.maxSize= 5e9)
plan(multicore, workers = 9, gc = TRUE) #I should try using multicore

list_res_univariate <- future_map(list_harm_univariate, function(x) {GagnonMR::all_mr_methods(x)}, .options = furrr_options(seed = TRUE))
res_univariate <- rbindlist(list_res_univariate, fill = TRUE)
res_univariate[, c("id.exposure", "id.outcome") := NULL]

FandQ <- lapply(list_harm_univariate, function(x) {
  res <- TwoSampleMR::mr_heterogeneity(x)
  x <- TwoSampleMR::add_rsq(x)
  res$fstat<-GagnonMR::fstat_fromdat(list(x))
  res$rsq <- sum(x$rsq.exposure)
  return(res)
}) %>% rbindlist(.,fill = TRUE)

FandQ[,c("id.exposure", "id.outcome") := NULL]
fwrite(FandQ, "Data/Modified/FandQ_univariate.txt")
fwrite(res_univariate, "Data/Modified/res_univariate.txt")

#### multivariate with TwoSample MR
#"res_bmiwc_nafld"
performmvmr <- function(exposure_vec, outcome_vec, pval_threshold = 1, 
                        clump_exp_arg = "none" ) { #clump_exp_arg either none, first, or second
  exposure_dat <- inst_all_sign_clump[exposure %in% exposure_vec,]
  d1 <- all_inst_mvmr[(exposure %in% exposure_vec) & (SNP %in% unique(exposure_dat$SNP)),] 
  
  if(clump_exp_arg == "none") {clump_exp<-NULL} else if(clump_exp_arg == "first"){clump_exp<-exposure_vec[1]} else if(clump_exp_arg == "second"){clump_exp<-exposure_vec[2]}
  inst_mvmr <- prepare_for_mvmr(exposure_dat = exposure_dat, d1 =d1, pval_threshold = pval_threshold, clump_exp = clump_exp)
  
  exposure_outcome_harmonized <- mv_harmonise_data(exposure_dat = inst_mvmr,
                                                   outcome_dat = all_outcome_mvmr[outcome == outcome_vec,],
                                                   harmonise_strictness = 1)
  mvmr_results <- GagnonMR::mv_multiple_MendelianRandomization(exposure_outcome_harmonized = exposure_outcome_harmonized)
  mvmr_results[,clump_exposure := clump_exp %>% ifelse(is.null(.), "none", .)]
  return(mvmr_results)
}

args <- list(list(exposure_vec = c("WC_UKB",  "BMI_UKB"), outcome_vec = "NAFLD"),
             list(exposure_vec = c("WC_UKB",  "BMI_UKB"), outcome_vec = "Fat_Liver"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Nikpay_CAD"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Nikpay_CAD"),
             list(exposure_vec = c("GIANT_2015_WC",  "GIANT_2015_BMI"), outcome_vec = "NAFLD"),
             list(exposure_vec = c("GIANT_2015_WC",  "GIANT_2015_BMI"), outcome_vec = "Fat_Liver"),
             list(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "Nikpay_CAD"),
             list(exposure_vec = c("GIANT_2015_WC",  "Fat_Liver"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("GIANT_2015_WC",  "Fat_Liver"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("GIANT_2015_WC",  "Fat_Liver"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("GIANT_2015_WC",  "Fat_Liver"), outcome_vec = "Nikpay_CAD"))


resmvmr <- map(args, function(x) {
  resnone<-performmvmr(exposure_vec = x$exposure_vec, outcome_vec = x$outcome_vec, clump_exp_arg = "none")
  resfirst <- performmvmr(exposure_vec = x$exposure_vec, outcome_vec = x$outcome_vec, clump_exp_arg = "first")
  ressecond <-performmvmr(exposure_vec = x$exposure_vec, outcome_vec = x$outcome_vec, clump_exp_arg = "second")
  return(rbindlist(list(resnone, resfirst, ressecond), fill = TRUE))
})

names(resmvmr) <- sapply(args, function(x) paste0(paste(x[[1]],collapse = "-and-"), "-on-", x[[2]]))

saveRDS(resmvmr, "Data/Modified/res_mvmr.rds")
setDT(egger_intercept)
egger_intercept <- egger_intercept[!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
fwrite(egger_intercept, "Data/Modified/egger_intercept.txt")

message("This script finished without errors")

