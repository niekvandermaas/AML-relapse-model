library(broom)
library(timeROC)
library(riskRegression)
library(tidyverse)

#define points
(coef <- as.data.frame(fit_dev$coefficients) %>%
  mutate(
    points = round(fit_dev$coefficients*4)
  )
)

#survival version
cox_points <- coxph(Surv(time, event) ~ age_cat + rfi_cat + wbc_cat + sct_allo + ch_11 + ck_mk + flt3itd + tp53 + inv16_t8_21,
                    complete_data_2, x = T, y = T)

#replace coef with points
cox_points$coefficients <- coef$points

#make prediction on dev data
data_pm <- complete_data_2 %>% 
  mutate(pred = predict(cox_points, complete_data_2))

#make clusters
data_pm$risk <- factor(ifelse(data_pm$pred <= 6, "Favorable", 
                              ifelse(data_pm$pred <=7, "Intermediate","Poor")), 
                       levels = c("Favorable", "Intermediate",  "Poor"))

#histogram of prognostic index
hist_pi <- data_pm %>%
  ggplot(aes(pred, fill = risk)) +
  geom_bar() +
  theme_fivethirtyeight() +
  theme(legend.position = "right",
        legend.direction='vertical',
        legend.title = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.title = element_text()) +
  scale_x_continuous(breaks = 0:18) + 
  xlab("Prognostic index score") +
  ylab("Number of patients")

set_palette(hist_pi, palette = c("#20854EFF", "#EFC000FF", "#CD534CFF"))


#km plot
ggsurvplot(
  survfit(Surv(complete_data_2$time, complete_data_2$event) ~ risk, data = data_pm),
  risk.table = TRUE,
  break.time.by = 6,
  xlim = c(0,12),
  xlab = "Time (months)",
  censor = F
)

#survival at 1 year
summary(survfit(Surv(complete_data_2$time, complete_data_2$event) ~ pred, data = data_pm), times = 12)

#calculate c-index
timeROC(T=complete_data_2$time, 
        delta=complete_data_2$event %in% 1,
        marker=as.numeric(data_pm$risk), 
        cause=1,
        weighting="marginal", 
        times=11.9, 
        iid=TRUE)
