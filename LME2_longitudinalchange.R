#LME2
#running linear mixed effects model to assess sensitivity of each method to capture longitudinal changes in FTD variants vs healthy controls
#model run separately for each method
#Measure ~ 1 + Diagnostic Group * Time from Baseline + Diagnostic Group * Age at Baseline + Sex + TIV + (1|Participant_ID)

library(lmerTest)
library(lme4)
library(multcomp)
library(dplyr)
library(tidyr)
library(lubridate)

#read data
data_lme <- read.csv(data_path, stringsAsFactors=TRUE)

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
#exclude participants with only baseline visit
#exclude visits over 3.9years after baseline (90th percentile)

data_lme <- data_lme %>% 
  filter(rowSums(!is.na(data_lme[, measures])) > 0)  %>%
  filter(DX!='OTHER') %>%
  filter(timefrombaseline < 1425)  %>%
  group_by(LONI_ID) %>%
  filter(n_distinct(visit_mri) > 1) %>%
  ungroup() %>%
  group_by(LONI_ID, visit_mri) %>% #only for subset analysis
  filter(n() >= 6) %>%  #only for subset analysis               
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

#run lme for measures that are available in cortical thickness

for (i in unique(data_lme$method)) {
  
  data_tmp <- data_lme %>%
    filter(method == i)
  
  tmp1 <- data.frame()
  
  for (variable in measures_all) {
    
    mylme <- lmer(
      as.formula(paste(variable, "~ DX*timefrombaselinez + DX*age_blz + GENDER + TIV_adjz + (1|LONI_ID)")),
      data = data_tmp
    )
    
    coefs <- summary(mylme)$coefficients
    
    BVp   <- coefs["DXBV:timefrombaselinez", "Pr(>|t|)"]
    BVt   <- coefs["DXBV:timefrombaselinez", "t value"]
    BVest <- coefs["DXBV:timefrombaselinez", "Estimate"]
    
    SVp   <- coefs["DXSV:timefrombaselinez", "Pr(>|t|)"]
    SVt   <- coefs["DXSV:timefrombaselinez", "t value"]
    SVest <- coefs["DXSV:timefrombaselinez", "Estimate"]
    
    PNFAp   <- coefs["DXPNFA:timefrombaselinez", "Pr(>|t|)"]
    PNFAt   <- coefs["DXPNFA:timefrombaselinez", "t value"]
    PNFAest <- coefs["DXPNFA:timefrombaselinez", "Estimate"]
    
    
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

#run lme for measures that are not available in cortical thickness

for (i in unique(data_lme$method[data_lme$method!='CT'])) {
  
  data_tmp <- data_lme %>%
    filter(method == i)
  
  tmp1 <- data.frame()
  
  for (variable in measures_with_vbm) {
    
    mylme <- lmer(
      as.formula(paste(variable, "~ DX*timefrombaselinez + DX*age_blz + GENDER + TIV_adjz + (1|LONI_ID)")),
      data = data_tmp
    )
    
    coefs <- summary(mylme)$coefficients
    
    BVp   <- coefs["DXBV:timefrombaselinez", "Pr(>|t|)"]
    BVt   <- coefs["DXBV:timefrombaselinez", "t value"]
    BVest <- coefs["DXBV:timefrombaselinez", "Estimate"]
    
    SVp   <- coefs["DXSV:timefrombaselinez", "Pr(>|t|)"]
    SVt   <- coefs["DXSV:timefrombaselinez", "t value"]
    SVest <- coefs["DXSV:timefrombaselinez", "Estimate"]
    
    PNFAp   <- coefs["DXPNFA:timefrombaselinez", "Pr(>|t|)"]
    PNFAt   <- coefs["DXPNFA:timefrombaselinez", "t value"]
    PNFAest <- coefs["DXPNFA:timefrombaselinez", "Estimate"]
    
    
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

write.csv(results, 'LME2_results.csv')
