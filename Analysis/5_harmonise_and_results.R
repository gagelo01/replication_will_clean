#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(MendelianRandomization)
library(writexl)
library(furrr)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
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
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Nikpay_CAD"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Mahajan_Type2diabetes"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "van_der_Harst_CAD"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Scott_Type2Diabetes"),
             list(exposure_vec = c("WC_UKB",  "Fat_Liver"), outcome_vec = "Nikpay_CAD"),
             list(exposure_vec = c("GIANT_2015_WC",  "GIANT_2015_BMI"), outcome_vec = "NAFLD"),
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
egger_intercept <- egger_intercept[!(exposure %in% c("Volume_VAT", "Volume_ASAT") |outcome %in% c("Volume_VAT", "Volume_ASAT")), ]
fwrite(egger_intercept, "Data/Modified/egger_intercept.txt")

message("This script finished without errors")

# performmvmr <- function(exposure_vec, outcome_vec, pval_threshold = 1, clump_exp = NULL) {
# exposure_dat <- inst_all_sign_clump[exposure %in% exposure_vec,]
# d1 <- all_inst_mvmr[(exposure %in% exposure_vec) & (SNP %in% unique(exposure_dat$SNP)),] 
# 
# inst_mvmr <- prepare_for_mvmr(exposure_dat = exposure_dat, d1 =d1, pval_threshold = pval_threshold, clump_exp = clump_exp)
# 
# exposure_outcome_harmonized <- mv_harmonise_data(exposure_dat = inst_mvmr,
#                                                  outcome_dat = all_outcome_mvmr[outcome == outcome_vec,],
#                                                  harmonise_strictness = 2)
# mvmr_results <- mv_multiple(exposure_outcome_harmonized)$result
# setDT(mvmr_results)
# mvmr_results[,method := "IVW"]
# mvmr_results[,CILower := b-1.96*se]
# mvmr_results[,CIUpper := b+1.96*se]
# mvmr_results<-mvmr_results[,.(method, exposure, outcome, b,se,CILower,CIUpper,pval, nsnp)]
# mvmr_results %>% setnames(., colnames(.), c("method", "Exposure", "Outcome",  "Estimate", "StdError", "CILower", "CIUpper",  "Pvalue", "nsnp"))
# return(mvmr_results)
# 
# # mvharm<- harmonise_data(exposure_dat = inst_mvmr,
# #                         outcome_dat = all_outcome_mvmr[outcome == outcome_vec,],
# #                         action = 2)
# # 
# # mv_input <- GagnonMR::mr_mvinput_wrapper(split(mvharm, mvharm$exposure))
# # 
# # res <- rbind(
# #   cbind(method = "IVW", GagnonMR::mr_mvivw_wrapper(mv_input)),
# #   cbind(method = "Egger", GagnonMR::mr_mvegger_wrapper(mv_input)))
# # 
# # return(res)
# }
# 
# res_bmiwcukb_nafld_eloi <- performmvmr(exposure_vec = c("WC_UKB",  "BMI_UKB"), outcome_vec = "NAFLD")
# res_wcnafldukb_t2d_eloi <- performmvmr(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "Mahajan_Type2diabetes", pval_threshold = 1, clump_exp = "NAFLD")
# res_wcnafldukb_cad_eloi <- performmvmr(exposure_vec = c("WC_UKB",  "NAFLD"), outcome_vec = "van_der_Harst_CAD", pval_threshold = 1, clump_exp = "NAFLD")
# res_bmiwcgiant_nafld_eloi <- performmvmr(exposure_vec = c("GIANT_2015_WC",  "GIANT_2015_BMI"), outcome_vec = "NAFLD")
# res_wcnafldgiant_t2d_eloi <- performmvmr(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "Mahajan_Type2diabetes", pval_threshold = 1, clump_exp = "NAFLD")
# res_wcnafldgiant_cad_eloi <- performmvmr(exposure_vec = c("GIANT_2015_WC",  "NAFLD"), outcome_vec = "van_der_Harst_CAD", pval_threshold = 1, clump_exp = "NAFLD")
# 
# 
# 
# objects <- c("res_univariate", "res_bmiwcukb_nafld_eloi", "res_wcnafldukb_t2d_eloi", "res_wcnafldukb_cad_eloi",
#              "res_bmiwcgiant_nafld_eloi", "res_wcnafldgiant_t2d_eloi", "res_wcnafldgiant_cad_eloi")
# 
# for(i in 1:length(objects)) {
#   fwrite(get(objects[i]), paste0("Results/", objects[i], ".txt")) 
# }
# 





# ########Three sample MR
# #UKB -> GIANT -> NAFLD
# snp_to_get <- inst_all_sign_clump[grepl("UKB", exposure), ]$SNP %>% unique
# outcome <- TwoSampleMR::extract_outcome_data(snps = snp_to_get, outcomes = c("ieu-a-2", "ieu-a-61"))
# exposure2 <- TwoSampleMR::convert_outcome_to_exposure(outcome_dat = outcome)
# outcome <- all_outcome_mvmr[outcome == "NAFLD" & SNP %in% exposure2$SNP,][order(SNP)]
# list_all_inst <- split(exposure2, exposure2$exposure)
# all_inst_mvmr <- GagnonMR::align_list_exposure_mvmr(list_all_inst, action = 2)
# 
# mvharm<- harmonise_data(exposure_dat = all_inst_mvmr, 
#                         outcome_dat = outcome, 
#                         action = 2)
# 
# mv_input <- GagnonMR::mr_mvinput_wrapper(split(mvharm, mvharm$exposure))
# 
# res <- rbind(
#   cbind(method = "IVW", GagnonMR::mr_mvivw_wrapper(mv_input)), 
#   cbind(method = "Egger", GagnonMR::mr_mvegger_wrapper(mv_input)))
# 
# #GIANT -> UKB -> NAFLD
# snp_to_get <- inst_all_sign_clump[grepl("GIANT", exposure), ]$SNP %>% unique
# exposure2 <- all_inst_mvmr[grepl("UKB", exposure) & SNP %in% snp_to_get,]
# exposure2 <- GagnonMR::prepare_for_mvmr(exposure_dat = exposure2)
# outcome <- all_outcome_mvmr[outcome == "NAFLD" & SNP %in% exposure2$SNP,][order(SNP)]
# list_all_inst <- split(exposure2, exposure2$exposure)
# all_inst_mvmr <- GagnonMR::align_list_exposure_mvmr(list_all_inst, action = 2)
# 
# mvharm<- harmonise_data(exposure_dat = all_inst_mvmr, 
#                         outcome_dat = outcome, 
#                         action = 2)
# 
# mv_input <- GagnonMR::mr_mvinput_wrapper(split(mvharm, mvharm$exposure))
# 
# res <- rbind(
#   cbind(method = "IVW", GagnonMR::mr_mvivw_wrapper(mv_input)), 
#   cbind(method = "Egger", GagnonMR::mr_mvegger_wrapper(mv_input)))
# 
# 
# 
# 
# # 
# # unique(inst_all_sign_clump$exposure)
# # instrument_bmiwc <- unique(inst_all_sign_clump[exposure %in% c("GIANT_2015_WC",  "GIANT_2015_BMI"),]$SNP)
# # 
# # mvharm<- mv_harmonise_data(exposure_dat = all_inst_mvmr[SNP %in% instrument_bmiwc & exposure %in% c("GIANT_2015_WC",  "GIANT_2015_BMI")], 
# #                            outcome_dat = all_outcome_mvmr[outcome == "NAFLD",], 
# #                            harmonise_strictness = 2)
# # res_bmiwc_nafld <- mv_multiple(mvharm)
# # res_bmiwc_nafld
# # 
# # #"res_wcnafld_cad"
# # unique(inst_all_sign_clump$exposure)
# # unique(all_outcome_mvmr$outcome)
# # instrument_wcnafld <- unique(inst_all_sign_clump[exposure %in% c("GIANT_2015_WC",  "NAFLD"),]$SNP)
# # 
# # mvharm<- mv_harmonise_data(exposure_dat = all_inst_mvmr[SNP %in% instrument_wcnafld & exposure %in% c("GIANT_2015_WC",  "NAFLD")], 
# #                            outcome_dat = all_outcome_mvmr[outcome == "van_der_Harst_CAD",], 
# #                            harmonise_strictness = 2)
# # res_wcnafld_cad <- mv_multiple(mvharm)
# # res_wcnafld_cad  #only one snp because select pvalue < 5*10^-8
# # 
# # #"res_wcnafld_t2d"
# # mvharm<- mv_harmonise_data(exposure_dat = all_inst_mvmr[SNP %in% instrument_wcnafld & exposure %in% c("GIANT_2015_WC",  "NAFLD")], 
# #                            outcome_dat = all_outcome_mvmr[outcome == "Mahajan_Type2diabetes",], 
# #                            harmonise_strictness = 2)
# # res_wcnafld_t2d <- mv_multiple(mvharm)
# # res_wcnafld_t2d  #only one snp because select pvalue < 5*10^-8
# # 
# # 
# # 
# # #####MVMR my way
# # 
# # #"res_bmiwc_nafld"
# # list_all_list_mvmr  <- split(all_inst_mvmr, all_inst_mvmr$exposure) 
# # list_inst_bmiwc <- lapply(list_all_list_mvmr[ c("GIANT_2015_WC","GIANT_2015_BMI")], function(x) x[SNP %in% instrument_bmiwc, ])
# # 
# # 
# # list_harm_bmiwc_nafld <- lapply(list_inst_bmiwc, function(x) {
# #   harmonise_data(x, all_outcome_mvmr[outcome == "NAFLD",], 2)
# # })
# # lapply(list_harm_bmiwc_nafld, setDT)
# # res_bmiwc_nafld_eloi <- GagnonMR::mr_mvinput_wrapper(list_harm_bmiwc_nafld) %>% GagnonMR::mr_mvivw_wrapper(.)
# # 
# # #"res_nalfdwc_cad"
# # list_inst_wcnafld <- lapply(list_all_list_mvmr[ c("GIANT_2015_WC","NAFLD")], function(x) x[SNP %in% instrument_wcnafld, ])
# # list_harm_wcnafld_cad <- lapply(list_inst_wcnafld, function(x) {
# #   harmonise_data(x, all_outcome_mvmr[outcome == "van_der_Harst_CAD",], 2)
# # })
# # 
# # res_wcnafld_cad_eloi <- GagnonMR::mr_mvinput_wrapper(list_harm_wcnafld_cad) %>% GagnonMR::mr_mvivw_wrapper(.)
# # 
# # #"res_wcnafld_t2d"
# # list_harm_wcnafld_t2d <- lapply(list_inst_wcnafld, function(x) {
# #   harmonise_data(x, all_outcome_mvmr[outcome == "Mahajan_Type2diabetes",], 2)
# # })
# # 
# # res_wcnafld_t2d_eloi <- GagnonMR::mr_mvinput_wrapper(list_harm_wcnafld_t2d) %>% GagnonMR::mr_mvivw_wrapper(.)
# # 
# # 
# # ###Compare method éloi with TwoSampleMR
# # res_bmiwc_nafld
# # res_bmiwc_nafld_eloi #more or less the same
# # 
# # res_wcnafld_cad
# # res_wcnafld_cad_eloi #more or less the same (so really it is the nsnp column that is not ok)
# # 
# # res_wcnafld_t2d
# # res_wcnafld_t2d_eloi
# 
# 
# ##write everything
