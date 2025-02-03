set.seed(123) #seed for reproducibility
library(survival)
library(survminer)
library(tidyverse)
library(caret)
library(timeROC)

#load data
test_data 

#select patients in cr and relapse. Remove certain IDs with other treatments
test_data <- test_data %>%
  filter(BestResponse %in% c("CR-Complete remission", "CRi-Complete remission with Incomplete Recovery")) %>% # Best response CR/CRi
  drop_na(rel_date) # drop NA relapse date
  

test_data <- test_data %>%
  mutate(
    sct_allo = ifelse(AlloSCTin1stCR == "Yes" & DateAlloSCTin1stCR <= rfs_date, 1, 0),
    rfi = test_data2$rfi,
    age_rel = Age + RFSurvival,
    os_rel = rfs_date %--% DateLastSeen %/% ddays(1)/30.4368499 + 0.1, 
    osi_rel = Death %in% "Died",
    flt3itd = ifelse(FLT3ITD == "Yes", 1, 0)
  ) %>% 
  filter(!ID %in% c(2268, 2274, 2281, 2292, 2295, 2300, 2301, 2309, 2318, 2320, 2325, 2328, 2341, 2350, 2360)) #IDs to filter

#set up new data
test <- test_data %>% 
  mutate(
    tp53 = factor(ifelse(TP53 == "Yes", 1,0)),
    inv16_t8_21 = ifelse(inv_16 == "Abnormality present" | t8_21 == "Abnormality present", 0, 1),
    ck_mk = ifelse(ck == "Abnormality present" | mk == "Abnormality present", 1, 0),
    ch_11 = ifelse(abn11q23 == "Abnormality present" | tv_11 == "Abnormality present", 1, 0),
    rfi_cat = factor(cut(rfi, breaks = c(-Inf, 12, Inf)), labels = c(1, 0)),
    age_cat = factor(cut(age_rel, breaks = c(-Inf, 60, Inf), labels = c(0,1))),
    wbc_cat = factor(cut(WBC, breaks = c(-Inf, 10, Inf), labels = c(0,1))),
    rfi_cat = factor(rfi_cat, levels = c("0", "1")),
    wbc = WBC,
    rfi = as.numeric(rfi)
  ) %>% select(sct_allo, ch_11, flt3itd, ck_mk, tp53, rfi_cat, inv16_t8_21, age_cat, wbc_cat, os_rel, osi_rel) %>%
  mutate(across(1:8, as.factor))

summary(test)

#predict
test <- test %>% 
  mutate(
    pred = predict(cox_points, test),
    risk  =  ifelse(pred <= 6, 1, ifelse(pred <= 7, 2, 3))
  )

#km plot
ggsurvplot(
  survfit(Surv(os_rel, osi_rel) ~ risk, test),
  risk.table = TRUE,
  break.time.by = 6,
  xlim = c(0,48),
  xlab = "Time (months)",
  censor = T
)

#survival at 1 year
summary(survfit(Surv(os_rel, osi_rel) ~ risk, test), times = 42)

# Calculate c-index at 1 year
(c_r <- timeROC(T=test$os_rel, 
        delta=test$osi_rel %in% 1,
        marker=as.numeric(test$risk), 
        cause=1,
        weighting="marginal", 
        times=12, 
        iid=TRUE))
