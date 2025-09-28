########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

# set studydrug variable that will be used to define reference group
studydrug_var = paste0("studydrug", main)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

# pARR = S0(t)^HR - S0(t)

cohort <- cohort %>% 
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100,
         )

temp <- cohort

reference_group = levels(temp[[studydrug_var]])[1]
############################1 ASSESSING CALIBRATION OF RISK SCORE IN PRESERVED EGFR################################################################

## remove double overlapping entries for DPP4 and SU that overlap (take one only) and select observations with preserved eGFR only
cohort <- temp %>% filter(.imp != 0 ) %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup() %>% filter(preegfr >= 60)

# check number of subjects
table(cohort[[studydrug_var]])

# make variable for risk deciles
cohort$risk_decile <- ntile(cohort$ckdpc_40egfr_score, n.quantiles)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort %>%
  group_by(risk_decile, !!sym(studydrug_var)) %>%
  summarise(mean_40egfr_pred=mean(ckdpc_40egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_40egfr_pred=mean(ckdpc_40egfr_score)/100)

## Find actual observed probabilities by risk score category and studydrug

EST.ref <- SE.ref <-
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_ref <- tibble() %>% mutate(
  observed_ref=NA,
  lower_ci_ref=NA,
  upper_ci_ref=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:n.quantiles) {
  for (i in 1:n.imp) {
    
    observed_ref_40egfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
                                     data=cohort[cohort$.imp == i & 
                                                   cohort$risk_decile == k &
                                                   cohort[[studydrug_var]]==reference_group,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))
    
    EST.ref[k,i] <- observed_ref_40egfr$estimate
    SE.ref[k,i] <- observed_ref_40egfr$std.error
    
    
    # observed_all_40egfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile, 
    #                               data=cohort[cohort$.imp == i & 
    #                                             cohort$risk_decile == k,]) %>%
    #   tidy() %>%
    #   # group_by(strata) %>%
    #   filter(time==max(time))
    # 
    # EST.all[k,i] <- observed_all_40egfr$estimate
    # SE.all[k,i] <- observed_all_40egfr$std.error
  }
  
  est.ref <- pool.rubin.KM(EST.ref[k,], SE.ref[k,], n.imp)
  observed_ref[k,] <- observed_ref[k,] %>% 
    mutate(
      observed_ref=est.ref[1],
      lower_ci_ref=est.ref[2],
      upper_ci_ref=est.ref[3],
      strata=k
    )

  
  # est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  # observed_all[k,] <- observed_all[k,] %>% 
  #   mutate(
  #     observed=est.all[1],
  #     lower_ci=est.all[2],
  #     upper_ci=est.all[3],
  #     strata=k
  #   )
  
}


events_ref <- cohort %>%
  filter(!!sym(studydrug_var)==reference_group & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(events=round(n()/n.imp, 0))


obs_v_pred <- cbind((predicted %>% filter(!!sym(studydrug_var)==reference_group)), observed_ref) %>%
  mutate(observed=observed_ref,
         lower_ci=lower_ci_ref,
         upper_ci=upper_ci_ref)

events_table <- data.frame(t(events_ref)) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_40egfr_pred=NA, risk_decile=0)

max_y_value <- ((obs_v_pred$upper_ci*1000) %>% max(na.rm = T) %>% ceiling())/10

## FINAL PLOT
p_uncal_bydeciles_presegfr_ref <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_40egfr_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_ref*100,ymin=lower_ci_ref*100, color=!!sym(studydrug_var)),width=0.1,linewidth=1) +
  geom_point(aes(y = observed_ref*100, group=!!sym(studydrug_var), color=!!sym(studydrug_var)), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Predicted 3-year risk of kidney disease progression (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100), breaks = seq(0, max_y_value, 2))+
  scale_y_continuous(limits=c(-1,100), breaks = seq(0, max_y_value, 2)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.35)),
        axis.title=element_text(size=rel(1.35)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  #ggtitle("Predicted risk of kidney disease progression") +
  coord_cartesian(xlim = c(0,max_y_value), ylim = c(0,max_y_value)) +
  ggtitle("CKD-PC risk score (≥60mL/min/1.73m2)")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_uncalibrated_risk_score_calibration_presegfr.tiff"), width=6, height=5.5, units = "in", res=800) 
print(p_uncal_bydeciles_presegfr_ref)
dev.off()

## C-stat
cohort <- cohort %>%
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100)

raw_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=cohort[cohort[[studydrug_var]] == reference_group,], method="breslow")
cstat_est <- summary(raw_mod)$concordance[1]
cstat_est_ll <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_est_ul <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])

## AUC
ROC_raw <- roc(cohort, ckd_egfr40_censvar, ckdpc_40egfr_score)
auc(ROC_raw)
ci.auc(ROC_raw)


## brier score for raw risk score
brier_raw <- rep(NA, n.imp)
brier_raw_se <- rep(NA, n.imp)


for (i in 1:n.imp) {
  print(paste("Imputation ", i))

  temp1 <- cohort %>% filter(.imp == i) %>% select(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar, ckdpc_40egfr_survival)
  temp1 <- temp1 %>% # error if times = 3, therefore adding extra row with time beyond t=3
    rbind(    # adds rows below your dataset
      temp1 %>%
        slice(1) %>% # this selects the first patients in dataset
        mutate(ckd_egfr40_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  raw_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival,
                   data=temp1, x=T)


  score_raw <-
    Score(object = list(raw_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ 1, # null model
          data = temp1,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T)

  brier_raw[i] <- score_raw$Brier$score$Brier[2]
  brier_raw_se[i] <- score_raw$Brier$score$se[2]

  rm(temp1)
}

brier_raw_se_pooled <- sqrt(mean(brier_raw_se^2) + (1+1/n.imp)*var(brier_raw))

## calibration slope:

# create empty vectors to store values of each calculation in multiple imputations:
bh_new <- rep(NA, n.imp)
se_bh_new <- rep(NA, n.imp)
cal_slope <- rep(NA, n.imp)
slope_optimism_presegfr <- rep(NA, n.imp)
var_slope <- rep(NA, n.imp)

options(datadist=NULL)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  # fit model with linear predictor as only variable
  recal_mod2 <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ (ckdpc_40egfr_lin_predictor), 
                    data = cohort %>% filter(.imp == i & !!sym(studydrug_var) == reference_group), x = TRUE, y = TRUE, surv = TRUE)
  
  # Obtain baseline survival estimates for model
  x <- summary(survfit(recal_mod2),time=3)
  bh_new[i] <- x$surv
  se_bh_new[i] <- x$std.err
  
  # bootstrap internal validation
  boot <- validate(recal_mod2)
  slope_optimism_presegfr[i] <- boot[3,5]
  
  # store calibration slope
  cal_slope[i] <- recal_mod2$coefficients * boot[3,5]
  var_slope[i] <- recal_mod2$var
}

# pool baseline hazard results
bh_recal <- mean(bh_new)
bh_recal_se <- sqrt(mean(se_bh_new)^2 + (1+1/n.imp)*var(bh_new))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_recal, ", 95% CI ", bh_recal-1.96*bh_recal_se, "-", bh_recal+1.96*bh_recal_se))

# pool calibration slope results
coef_recal <- mean(cal_slope)
coef_recal_se <- sqrt(mean(var_slope) + (1+1/n.imp)*var(cal_slope))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_recal, ", 95% CI ", coef_recal-1.96*coef_recal_se, "-", coef_recal+1.96*coef_recal_se))

# c statistic:
paste0("C statistic in subjects with preserved eGFR: ", round(cstat_est, 4), ", 95% CI ", round(cstat_est_ll, 4), "-", round(cstat_est_ul,4))

#pool and print brier score
print(paste0("Brier score for uncalibrated risk score in subjects with preserved eGFR: ", mean(brier_raw), ", 95% CI ", mean(brier_raw)-1.96*brier_raw_se_pooled, "-", mean(brier_raw)+1.96*brier_raw_se_pooled))

model_metrics_presegfr <- data.frame(
  model = "Preserved eGFR",
  baseline_hazard = round(bh_recal, 4),
  baseline_hazard_ci = paste0(round(bh_recal-1.96*bh_recal_se, 4), "-", round(bh_recal+1.96*bh_recal_se, 4)),
  calibration_slope = round(coef_recal, 4),
  calibration_slope_ci = paste0(round(coef_recal-1.96*coef_recal_se, 4), "-", round(coef_recal+1.96*coef_recal_se, 4)),
  c_statistic = round(cstat_est, 4),
  c_statistic_ci = paste0(round(cstat_est_ll, 4), "-", round(cstat_est_ul,4)),
  brier_score = mean(brier_raw),
  brier_score_ci = paste0(mean(brier_raw)-1.96*brier_raw_se_pooled, "-", mean(brier_raw)+1.96*brier_raw_se_pooled)
)


############################2 ASSESSING CALIBRATION OF RISK SCORE IN REDUCED EGFR################################################################

## remove double overlapping entries for DPP4 and SU that overlap (take one only) and select observations with reduced eGFR only
cohort <- temp %>% filter(.imp != 0 ) %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>%
  arrange(dstartdate) %>%
  filter(!duplicated(!!sym(studydrug_var))) %>%
  ungroup() %>% filter(preegfr < 60)

# check number of subjects
table(cohort[[studydrug_var]])

# make variable for risk deciles
cohort$risk_decile <- ntile(cohort$ckdpc_40egfr_score, n.quantiles/2)

## Get mean predicted probabilities by risk decile and studydrug
predicted <- cohort %>%
  group_by(risk_decile, !!sym(studydrug_var)) %>%
  summarise(mean_40egfr_pred=mean(ckdpc_40egfr_score)/100)

# get mean predicted probabilities by risk decile (not by studydrug)
predicted_all <- cohort %>%
  group_by(risk_decile) %>%
  summarise(mean_40egfr_pred=mean(ckdpc_40egfr_score)/100)

## Find actual observed probabilities by risk score category and studydrug

EST.ref <- SE.ref <-
  EST.all <- SE.all <-
  matrix(data = NA, nrow = n.quantiles, ncol = n.imp)

observed_ref <- tibble() %>% mutate(
  observed_ref=NA,
  lower_ci_ref=NA,
  upper_ci_ref=NA,
  strata=NA
)

observed_all <- tibble() %>% mutate(
  observed=NA,
  lower_ci=NA,
  upper_ci=NA,
  strata=NA
)

for (k in 1:(n.quantiles/2)) {
  for (i in 1:n.imp) {

    observed_ref_40egfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile,
                                     data=cohort[cohort$.imp == i &
                                                   cohort$risk_decile == k &
                                                   cohort[[studydrug_var]]==reference_group,]) %>%
      tidy() %>%
      # group_by(strata) %>%
      filter(time==max(time))

    EST.ref[k,i] <- observed_ref_40egfr$estimate
    SE.ref[k,i] <- observed_ref_40egfr$std.error


    # observed_all_40egfr <- survfit(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ risk_decile,
    #                               data=cohort[cohort$.imp == i &
    #                                             cohort$risk_decile == k,]) %>%
    #   tidy() %>%
    #   # group_by(strata) %>%
    #   filter(time==max(time))
    #
    # EST.all[k,i] <- observed_all_40egfr$estimate
    # SE.all[k,i] <- observed_all_40egfr$std.error
  }

  est.ref <- pool.rubin.KM(EST.ref[k,], SE.ref[k,], n.imp)
  observed_ref[k,] <- observed_ref[k,] %>%
    mutate(
      observed_ref=est.ref[1],
      lower_ci_ref=est.ref[2],
      upper_ci_ref=est.ref[3],
      strata=k
    )

  # est.all <- pool.rubin.KM(EST.all[k,], SE.all[k,], n.imp)
  # observed_all[k,] <- observed_all[k,] %>%
  #   mutate(
  #     observed=est.all[1],
  #     lower_ci=est.all[2],
  #     upper_ci=est.all[3],
  #     strata=k
  #   )

}


events_ref <- cohort %>%
  filter(!!sym(studydrug_var)==reference_group & ckd_egfr40_censvar==1) %>%
  group_by(risk_decile) %>%
  summarise(events=round(n()/n.imp, 0))

obs_v_pred <- cbind((predicted %>% filter(!!sym(studydrug_var)==reference_group)), observed_ref)

events_table <- data.frame(t(events_ref)) %>%
  rownames_to_column() %>%
  filter(rowname!="risk_decile")

dodge <- position_dodge(width=0.3)

empty_tick <- obs_v_pred %>%
  filter(risk_decile==1) %>%
  mutate(observed=NA, lower_ci=NA, upper_ci=NA, mean_40egfr_pred=NA, risk_decile=0)

max_y_value <- ((obs_v_pred$upper_ci_ref*1000) %>% max(na.rm = T) %>% ceiling())/10

## FINAL PLOT
p_uncal_bydeciles_redegfr_ref <- ggplot(data=bind_rows(empty_tick,obs_v_pred), aes(x=mean_40egfr_pred*100)) +
  geom_errorbar(aes(ymax=upper_ci_ref*100,ymin=lower_ci_ref*100, color=!!sym(studydrug_var)),width=0.6,size=1) +
  geom_point(aes(y = observed_ref*100, group=!!sym(studydrug_var), color=!!sym(studydrug_var)), shape=18, size=3) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_bw() +
  xlab("Predicted 3-year risk of kidney disease progression (%)") + ylab("Observed risk (%)")+
  scale_x_continuous(limits=c(0,100))+
  scale_y_continuous(limits=c(-10,100)) +
  scale_colour_manual(values = cols) +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
  theme(axis.text=element_text(size=rel(1.35)),
        axis.title=element_text(size=rel(1.35)),
        plot.title=element_text(hjust = 0.5),
        plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
        legend.position = "none") +
  #ggtitle("Predicted risk of kidney disease progression") +
  coord_cartesian(xlim = c(0,max_y_value), ylim = c(0,max_y_value)) +
  ggtitle("CKD-PC risk score (<60mL/min/1.73m2)")


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_uncalibrated_risk_score_calibration_redegfr.tiff"), width=6, height=5.5, units = "in", res=800)
print(p_uncal_bydeciles_redegfr_ref)
dev.off()

## C-stat

raw_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival, data=cohort[cohort[[studydrug_var]] == reference_group,], method="breslow")
cstat_est_redegfr <- summary(raw_mod)$concordance[1]
cstat_est_ll_redegfr <- summary(raw_mod)$concordance[1]-(1.96*summary(raw_mod)$concordance[2])
cstat_est_ul_redegfr <- summary(raw_mod)$concordance[1]+(1.96*summary(raw_mod)$concordance[2])

## AUC
ROC_raw <- roc(cohort, ckd_egfr40_censvar, ckdpc_40egfr_score)
auc(ROC_raw)
ci.auc(ROC_raw)


## brier score for raw risk score
brier_raw_redegfr <- rep(NA, n.imp)
brier_raw_se_redegfr <- rep(NA, n.imp)


for (i in 1:n.imp) {
  print(paste("Imputation ", i))

  temp1 <- cohort %>% filter(.imp == i) %>% select(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar, ckdpc_40egfr_survival)
  temp1 <- temp1 %>% # error if times = 3, therefore adding extra row with time beyond t=3
    rbind(    # adds rows below your dataset
      temp1 %>%
        slice(1) %>% # this selects the first patients in dataset
        mutate(ckd_egfr40_censtime_yrs = 3.5)   # changes censored time to 3.5
    )
  raw_mod <- coxph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ ckdpc_40egfr_survival,
                   data=temp1, x=T)


  score_raw <-
    Score(object = list(raw_mod), # need to pass cox model to Score() as a list in order for it to be processed
          formula = Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ 1, # null model
          data = temp1,
          summary = "ibs", # statistic of interest is integrated brier score
          times = 3,
          splitMethod = "bootcv",  # use bootstrapping for confidence intervals
          B = n.bootstrap,
          verbose = T)

  brier_raw_redegfr[i] <- score_raw$Brier$score$Brier[2]
  brier_raw_se_redegfr[i] <- score_raw$Brier$score$se[2]

  rm(temp1)
}

brier_raw_se_pooled_redegfr <- sqrt(mean(brier_raw_se_redegfr^2) + (1+1/n.imp)*var(brier_raw_redegfr))

### calibration slope

# create empty vectors to store values of each calculation in multiple imputations:
bh_new_redegfr <- rep(NA, n.imp)
se_bh_new_redegfr <- rep(NA, n.imp)
cal_slope_redegfr <- rep(NA, n.imp)
slope_optimism_redegfr <- rep(NA, n.imp)
var_slope_redegfr <- rep(NA, n.imp)

options(datadist=NULL)

for (i in 1:n.imp) {
  print(paste0("Calculations in imputation ", i))
  # fit model with linear predictor as only variable
  recal_mod2 <- cph(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ (ckdpc_40egfr_lin_predictor),
                    data = cohort %>% filter(.imp == i & !!sym(studydrug_var) == reference_group), x = TRUE, y = TRUE, surv = TRUE)

  # Obtain baseline survival estimates for model
  x <- summary(survfit(recal_mod2),time=3)
  bh_new_redegfr[i] <- x$surv
  se_bh_new_redegfr[i] <- x$std.err

  # bootstrap internal validation
  boot <- validate(recal_mod2)
  slope_optimism_redegfr[i] <- boot[3,5]

  # store calibration slope
  cal_slope_redegfr[i] <- recal_mod2$coefficients * boot[3,5]
  var_slope_redegfr[i] <- recal_mod2$var
}

# pool baseline hazard results
bh_recal_redegfr <- mean(bh_new_redegfr)
bh_recal_se_redegfr <- sqrt(mean(se_bh_new_redegfr)^2 + (1+1/n.imp)*var(bh_new_redegfr))
# print baseline hazard with 95% CI
print(paste0("Baseline hazard ", bh_recal_redegfr, ", 95% CI ", bh_recal_redegfr-1.96*bh_recal_se_redegfr, "-", bh_recal_redegfr+1.96*bh_recal_se_redegfr))

# pool calibration slope results
coef_recal_redegfr <- mean(cal_slope_redegfr)
coef_recal_se_redegfr <- sqrt(mean(var_slope_redegfr) + (1+1/n.imp)*var(cal_slope_redegfr))
# print calibration slope with 95% CI
print(paste0("Calibration slope ", coef_recal_redegfr, ", 95% CI ", coef_recal_redegfr-1.96*coef_recal_se_redegfr, "-", coef_recal_redegfr+1.96*coef_recal_se_redegfr))


# summary statistics for discimination and overall performance:
paste0("C statistic in subjects with reduced eGFR: ", round(cstat_est_redegfr, 4), ", 95% CI ", round(cstat_est_ll_redegfr, 4), "-", round(cstat_est_ul_redegfr,4))
print(paste0("Brier score for uncalibrated risk score in subjects with reduced eGFR: ", mean(brier_raw_redegfr), ", 95% CI ", mean(brier_raw_redegfr)-1.96*brier_raw_se_pooled_redegfr, "-", mean(brier_raw_redegfr)+1.96*brier_raw_se_pooled_redegfr))


model_metrics_redegfr <- data.frame(
  model = "Reduced eGFR",
  baseline_hazard = round(bh_recal_redegfr, 4),
  baseline_hazard_ci = paste0(round(bh_recal_redegfr-1.96*bh_recal_se_redegfr, 4), "-", round(bh_recal_redegfr+1.96*bh_recal_se_redegfr, 4)),
  calibration_slope = round(coef_recal_redegfr, 4),
  calibration_slope_ci = paste0(round(coef_recal_redegfr-1.96*coef_recal_se_redegfr, 4), "-", round(coef_recal_redegfr+1.96*coef_recal_se_redegfr, 4)),
  c_statistic = round(cstat_est_redegfr, 4),
  c_statistic_ci = paste0(round(cstat_est_ll_redegfr, 4), "-", round(cstat_est_ul_redegfr,4)),
  brier_score = mean(brier_raw_redegfr),
  brier_score_ci = paste0(mean(brier_raw_redegfr)-1.96*brier_raw_se_pooled_redegfr, "-", mean(brier_raw_redegfr)+1.96*brier_raw_se_pooled_redegfr)
)

############################3 STORE MODEL METRICS################################################################

model_metrics <- rbind(model_metrics_presegfr, model_metrics_redegfr)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
write.csv2(model_metrics, file=paste0(today, "_ckdpc_egfr40_score_performance.csv"))
