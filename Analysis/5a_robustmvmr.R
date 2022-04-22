#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(MendelianRandomization)
library(writexl)
library(robustMVMR)

setwd("/mnt/sde/gagelo01/Projects/small_MR_exploration/replication_will_clean")
all_inst_proxies <- fread( "Data/Modified/all_inst_proxies.txt")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr")
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvm.txt" )
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")


robustMVMR_wrapper <- function(inpoot, pval_threshold = 1e-05){
  fit <- robustMVMR::robustMVMR(betaGX = as.matrix(inpoot[, colnames(inpoot)[grepl("beta.exposure",colnames(inpoot))]]),
                                sebetaGX = as.matrix(inpoot[, colnames(inpoot)[grepl("se.exposure",colnames(inpoot))]]),
                                pvalbetaGX = as.matrix(inpoot[, colnames(inpoot)[grepl("pval.exposure",colnames(inpoot))]]),
                                betaGY = inpoot[, colnames(inpoot)[grepl("beta.outcome",colnames(inpoot))]],
                                sebetaGY = inpoot[, colnames(inpoot)[grepl("se.outcome",colnames(inpoot))]], 
                                pvalbetaGY = inpoot[, colnames(inpoot)[grepl("pval.outcome",colnames(inpoot))]],
                                pval_threshold = pval_threshold, plot = FALSE)
  return(fit)

}

getinput_robustmvmr <- function (list_harm3)  {
  lapply(list_harm3, setDT)
  rsid_include <- Reduce(intersect, lapply(list_harm3, function(x) x$SNP))
  list_harm3 <- lapply(list_harm3, function(x) x[SNP %in% rsid_include, 
  ])
  stopifnot(all(list_harm3[[1]]$SNP == list_harm3[[2]]$SNP & 
                  list_harm3[[1]]$SNP == list_harm3[[length(list_harm3)]]$SNP))
  bon <- map(list_harm3, `[`, "beta.exposure")
  bx <- matrix(unlist(bon), ncol = length(list_harm3), nrow = nrow(bon[[1]]), 
               byrow = FALSE)
  bon <- map(list_harm3, `[`, "se.exposure")
  bxse <- matrix(unlist(bon), ncol = length(list_harm3), nrow = nrow(bon[[1]]), 
                 byrow = FALSE)
  bon <- map(list_harm3, `[`, "pval.exposure")
  bxpval <- matrix(unlist(bon), ncol = length(list_harm3), nrow = nrow(bon[[1]]), 
                   byrow = FALSE)
  
  stopifnot(all(lapply(list_harm3, function(x) x$beta.outcome == 
                         list_harm3[[1]]$beta.outcome) %>% unlist(.)))
  by <- list_harm3[[1]]$beta.outcome
  byse <- list_harm3[[1]]$se.outcome
  bypval <- list_harm3[[1]]$pval.outcome
  stopifnot(all(list_harm3[[1]]$effect_allele == list_harm3[[2]]$effect_allele & 
                  list_harm3[[1]]$effect_allele == list_harm3[[length(list_harm3)]]$effect_allele) & 
              all(list_harm3[[1]]$other_allele == list_harm3[[2]]$other_allele & 
                    list_harm3[[1]]$other_allele == list_harm3[[length(list_harm3)]]$other_allele) & 
              all(list_harm3[[1]]$eaf.exposure == list_harm3[[2]]$eaf.exposure & 
                    list_harm3[[1]]$eaf.exposuree == list_harm3[[length(list_harm3)]]$eaf.exposure))
  
  colnames(bx) <-   paste0(sapply(list_harm3, function(x) x[1, ]$exposure), ".beta.", "exposure")
  colnames(bxse) <-   paste0(sapply(list_harm3, function(x) x[1, ]$exposure), ".se.", "exposure")
  colnames(bxpval) <-   paste0(sapply(list_harm3, function(x) x[1, ]$exposure), ".pval.", "exposure")
  
  df_outcome <-data.frame(by, byse, bypval)
  colnames(df_outcome) <-  paste0(sapply(list_harm3, function(x) x[1, ]$outcome), c(".beta.", ".se.", ".pval."),  "outcome")
  
  
  inpoot<- cbind(bx, bxse, bxpval, df_outcome)
  
  return(inpoot)
}

exposure_vec <- c("GIANT_2015_WC",  "GIANT_2015_BMI")
outcome_vec <- "NAFLD"

performrobustMVMR <- function(exposure_vec, outcome_vec) {
  instrument <- unique(inst_all_sign_clump[exposure %in% exposure_vec,]$SNP)
  exposure_dat <- all_inst_mvmr[SNP %in% instrument & exposure %in% exposure_vec]
  
  exposure_dat <- GagnonMR::prepare_for_mvmr(exposure_dat = exposure_dat)
  
  list_inst_mvmr <- split(exposure_dat, exposure_dat$exposure)
  inst_mvmr <- GagnonMR::align_list_exposure_mvmr(list_inst_mvmr = list_inst_mvmr)
  
  if(all(exposure_vec %in% c("GIANT_2015_WC",  "Mahajan_Type2diabetes")) & outcome_vec == "NAFLD") {
    inst_mvmr<- inst_mvmr[SNP != "rs62007683"]
  } 
  
  mvharm<- harmonise_data(exposure_dat = inst_mvmr, 
                          outcome_dat = all_outcome_mvmr[outcome == outcome_vec,], 
                          action = 2)
  
  mv_input <- getinput_robustmvmr(split(mvharm, mvharm$exposure))
  res <- robustMVMR_wrapper(mv_input)
  return(res)
}


fit <- performrobustMVMR(exposure_vec = c("GIANT_2015_WC",  "GIANT_2015_BMI"), outcome_vec = "NAFLD")


## -- Main results of the robust MVMR
fit$mvMRResult_heter_robust
## -- The modified Q-statistic for testing instrument validity
fit$Q_pleiotropy_test
## -- The pair-wise conditional F-statistic matrix
fit$Conditional_F_statistic_matrix
## -- The correlation matrix of the exposures
fit$rho_Exposures


