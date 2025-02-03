#load packages
library(tidyverse)
library(haven)
library(lubridate)
library(ggthemes)
library(survival)
library(survminer)
library(Boruta)
library(timeROC)
library(rms)
library(broom)

###########################################
#HOVON AML studies
###########################################

#final dataset
complete_data <- dev_dat

##########################################
# Breems et al score
##########################################

# predictors in final model: ~age + rfi + cytogenetics(inv(16), t(16;16), t(8;21)) + allo and auto SCT
dat_breems <- complete_data %>% select(pid, trialnr, inv16, t8_21, sct_au, sct_allo, time, event, efs_rel, efsi_rel, rfs2, rfsi2, ttcr2, ttcri2) %>% mutate(
  cyto = case_when(
    inv16 == 1 ~ 1,
    t8_21 == 1 ~ 1,
    inv16 == 0 & t8_21 ==0 ~ 0
  ),
  sct_au_al = ifelse(sct_au == 1 | sct_allo == 1, 1, 0),
  age_rel = rel_data$age_rel,
  rfi = rel_data$rfi
)

#risk groups
dat_breems <- dat_breems %>% mutate(
  rfi_scr = case_when(
    between(rfi, 0, 6) ~ 5,
    between(rfi, 7, 18) ~ 3,
    between(rfi, 18, Inf) ~ 0,
  ),
  cyt_scr = case_when(
    inv16 == 1 ~ 0,
    t8_21 == 1 ~ 3,
    inv16 == 0 & t8_21 == 0 ~ 5,
  ),
  age_scr = case_when(
    between(age_rel, 0, 35) ~ 0,
    between(age_rel, 35, 45) ~ 1,
    between(age_rel, 45, Inf) ~ 2,
  ),
  sct_scr = case_when(
    sct_au_al == 1 ~ 2,
    sct_au_al == 0 ~ 0
  ),
  tot_scr = rfi_scr + cyt_scr + age_scr + sct_scr,
  group  = factor(ifelse(tot_scr <= 6, "Favorable", 
                         ifelse(tot_scr <= 9, "Intermediate", "Poor")), 
                  levels = c('Favorable', 'Intermediate', 'Poor'))
)



#plot risk groups
km_breems <- survfit(Surv(time, event %in% 1) ~ group, dat_breems)

ggsurvplot(km_breems,
           palette = c("#009E73", "#F0E442", "darkred"),
           risk.table = T,
           xlim = c(0,12),
           censor = F,
           break.time.by = 3) 

# 1-yr survival
summary(survfit(Surv(time, event) ~ group, dat_breems), times = 12)

#ROC
timeROC(dat_breems$time, dat_breems$event, 
        marker = as.numeric(as.factor(dat_breems$group)), 
        cause=1,
        weighting="marginal", 
        times=11.9, 
        iid=TRUE)

##############################################
# GOELAMS score
##############################################

# predictors in final model: ~ rfi + FLT3-ITD + cyto
dat_GOELAMS <- rel_data %>% 
  select(pid, trialnr, flt3itd, t8_21, t15_17, inv16, del5q, m5, m7, ck,
         starts_with("abn3"), abn03_q27, time, event, efs_rel, efsi_rel, rfs2, rfsi2, ttcr2, ttcri2) %>%
  mutate(
    rfi = complete_data$rfi,
    cyto = case_when(
      t8_21 == 1 | t15_17 == 1 | inv16 == 1 ~ "Favorable/intermediate",
      m5 == 1 | m7 == 1 | del5q == 1 | abn3q26 == 1 | ck2017 == 1 |
        abn3_invq21q26 == 1 | abn3_q25 == 1 | abn3_q21q26 == 1 | 
        abn03_q27 == 1 ~ "Adverse", 
      t8_21 == 0 & t15_17 == 0 & inv16 == 0 &  
        m5 == 0 & m7 == 0 & del5q == 0 & abn3q26 == 0 & 
        abn3_invq21q26 == 0 & abn3_q25 == 0 & abn3_q21q26 == 0 & abn03_q27 == 0 ~ "Favorable/intermediate"
    ),
    cyto = factor(cyto, levels = c("Favorable/intermediate", "Adverse"))
  ) %>%
  select(pid, trialnr, rfi, flt3itd, cyto, time, event, efs_rel, efsi_rel, rfs2, rfsi2, ttcr2, ttcri2)

#no adverse markers, so favorable/intermediate
dat_GOELAMS$cyto[is.na(dat_GOELAMS$cyto)] <- "Favorable/intermediate"

#risk groups
dat_GOELAMS <- dat_GOELAMS %>% mutate(
  rfi_scr = ifelse(rfi <= 12, 1, 0),
  flt_scr = ifelse(flt3itd == 2, 1, 0),
  cyto_scr = ifelse(cyto == "Adverse", 1, 0),
  tot_scr = rfi_scr + flt_scr + cyto_scr,
  group  = factor(ifelse(tot_scr == 0, "Favorable", 
                         ifelse(tot_scr == 1, "Intermediate", "Poor")), 
                  levels = c('Favorable', 'Intermediate', 'Poor'))
)

#plot risk groups
km_GOELAMS <- survfit(Surv(time, event) ~ group, dat_GOELAMS)

ggsurvplot(km_GOELAMS,
           palette = c("#009E73", "#F0E442", "darkred"),
           risk.table = T,
           xlim = c(0,12),
           censor = TRUE,
           break.time.by = 3) 

# 2-yr survival
summary(survfit(Surv(time, event) ~ group, dat_GOELAMS), times = 12)

#timeROC
timeROC(dat_GOELAMS$time, dat_GOELAMS$event, 
        marker = as.numeric(dat_GOELAMS$group), 
        cause=1,
        weighting="marginal", 
        times=11.9, 
        iid=TRUE)

######################################################
# Make a Random Survival Forest
######################################################

names_dat_com <- names(boruta$finalDecision[boruta$finalDecision == "Confirmed"])

rf_data <- complete_data %>% 
  select(all_of(names_dat_com))

#plot VIMP on z-score
set.seed(321)
z <- Surv(complete_data$time + 0.01, complete_data$event)
boruta <- Boruta(z ~., data = rf_data, doTrace = 2, maxRuns = 1000, mcAdj = T)

##################################################################
# Perform a Cox proportional hazards regression
##################################################################

names <- names(boruta$finalDecision[boruta$finalDecision == "Confirmed"])

#select confirmed predictors and preprocess data
data_cox <- complete_data %>%
  select(all_of(names)) %>%
  mutate(time = complete_data$time,
         event = complete_data$event)

#model formula
survival_formula <- paste0("Surv(time + 0.01, event)", " ~ ", 
                           paste0(names, collapse = " + "))
survival_formula <- as.formula(survival_formula)

#final model
(v <- rms::fastbw(cph(survival_formula, data_cox, x = T, y = T), rule = "p", sls = 0.1))

#fit final model in survival package
fit_dev <- coxph(Surv(time, event) ~ age_cat + rfi_cat + wbc_cat + sct_allo + ch_11 + ck_mk + flt3itd + tp53 + inv16_t8_21, data = complete_data, x = TRUE, y = TRUE)

#use tidy function from broom to get nice overview with 95% CI and P value
broom::tidy(fit_dev, conf.int = T, exponentiate = T)
