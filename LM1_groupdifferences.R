#LM1
#running linear regression to assess sensitivity of each method to capture group differences (FTD variants vs healthy controls)
#model run separately for each method
#Measure ~ 1 + Diagnostic Group + Age + Sex + TIV

library(lmerTest)
library(lme4)
library(multcomp)
library(dplyr)
library(tidyr)
library(lubridate)

#read data
data_lme <- read.csv("/data/dadmah/metame/NIFD_aim1/data/LME_data_251020.csv", stringsAsFactors=TRUE)

#get regional columns (CerebrA regions)
measures <- tail(colnames(data_lme), 90)

#DBM has signal for cortical grey matter, subcortical grey matter, white matter, ventricles
#CT only for cortical grey matter
#VBM only for cortical and subcortical grey matter
#so we need to differentiate measures depending on which methods have them

measures_no_ct <- measures[
  sapply(measures, function(v) all(is.na(data_lme[data_lme$method == "CT", v])))
]
measures_no_vbm <- measures[
  sapply(measures, function(v) all(is.na(data_lme[data_lme$method == "VBMdir", v])))
]
measures_all <- measures[!measures %in% c(measures_no_ct)]
measures_with_vbm <- measures[!measures %in% c(measures_no_vbm,measures_all)]

#exclude rows that are fully NA, 'other' diagnoses, and (only for subset analysis) only keep visits that passed QC for all methods
data_lme <- data_lme %>% 
  filter(rowSums(!is.na(data_lme[, measures])) > 0)  %>%
  filter(DX!='OTHER') %>%
  group_by(LONI_ID, visit_mri) %>% #only for subset analysis
  filter(n() >= 5) %>%  #only for subset analysis              
  ungroup() #only for subset analysis

#set healthy controls as reference
data_lme$DX <- factor(data_lme$DX, levels = c("CON", "BV", "SV", "PNFA"))
data_lme$GENDER <- as.factor(data_lme$GENDER)

results <- data.frame(
  measure = character(),
  EffectType = character(),
  pValue = numeric(),
  tStat = numeric(),
  est = numeric(),
  p_adj = numeric(),
  FDR = numeric(),
  method = character(),
  stringsAsFactors = FALSE
)

#run lm for measures that are available in cortical thickness

for (i in unique(data_lme$method)) {
  
  data_tmp <- data_lme %>%
    filter(method == i)
  
  tmp1 <- data.frame()
  
  for (variable in measures_all) {
    
    mylme <- lm(
      as.formula(paste(variable, "~ DX + agez + GENDER + TIV_adjz")),
      data = data_tmp
    )
    
    coefs <- summary(mylme)$coefficients
    
    BVp   <- coefs["DXBV", "Pr(>|t|)"]
    BVt   <- coefs["DXBV", "t value"]
    BVest <- coefs["DXBV", "Estimate"]
    
    SVp   <- coefs["DXSV", "Pr(>|t|)"]
    SVt   <- coefs["DXSV", "t value"]
    SVest <- coefs["DXSV", "Estimate"]
    
    PNFAp   <- coefs["DXPNFA", "Pr(>|t|)"]
    PNFAt   <- coefs["DXPNFA", "t value"]
    PNFAest <- coefs["DXPNFA", "Estimate"]
    
    
    tmp2 <- data.frame(
      measure    = rep(variable, 3),
      EffectType = c("BV", "SV", "PNFA"),
      pValue     = c(BVp, SVp, PNFAp),
      tStat      = c(BVt, SVt, PNFAt),
      est        = c(BVest, SVest, PNFAest),
      FDR        = NA_real_,
      stringsAsFactors = FALSE
    )
    tmp1 <- bind_rows(tmp1,tmp2)
    tmp1$method <- i
  }
  results <- bind_rows(results,tmp1)
}

#run lm for measures that are not available in cortical thickness

for (i in unique(data_lme$method[data_lme$method!='CT'])) {
  
  data_tmp <- data_lme %>%
    filter(method == i)
  
  tmp1 <- data.frame()
  
  for (variable in measures_with_vbm) {
    
    mylme <- lm(
      as.formula(paste(variable, "~ DX + agez + GENDER + TIV_adjz")),
      data = data_tmp
    )
    
    coefs <- summary(mylme)$coefficients
    
    BVp   <- coefs["DXBV", "Pr(>|t|)"]
    BVt   <- coefs["DXBV", "t value"]
    BVest <- coefs["DXBV", "Estimate"]
    
    SVp   <- coefs["DXSV", "Pr(>|t|)"]
    SVt   <- coefs["DXSV", "t value"]
    SVest <- coefs["DXSV", "Estimate"]
    
    PNFAp   <- coefs["DXPNFA", "Pr(>|t|)"]
    PNFAt   <- coefs["DXPNFA", "t value"]
    PNFAest <- coefs["DXPNFA", "Estimate"]
    
    
    tmp2 <- data.frame(
      measure    = rep(variable, 3),
      EffectType = c("BV", "SV", "PNFA"),
      pValue     = c(BVp, SVp, PNFAp),
      tStat      = c(BVt, SVt, PNFAt),
      est        = c(BVest, SVest, PNFAest),
      FDR        = NA_real_,
      stringsAsFactors = FALSE
    )
    tmp1 <- bind_rows(tmp1,tmp2)
    tmp1$method <- i
  }
  results <- bind_rows(results,tmp1)
}

#add fdr correction
results <- results %>%
  group_by(EffectType) %>%
  mutate(p_adj = p.adjust(pValue, method = "fdr")) %>%
  ungroup() %>%
  mutate(FDR = ifelse(p_adj <= .05, 1, 0))

write.csv(results, 'LM1_results_sub.csv')
