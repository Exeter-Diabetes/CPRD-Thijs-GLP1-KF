########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_with_observed_surv.Rda"))

studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)

############################1 COMPARE pARR TO EGFR/ALBUMINURIA CRITERIA################################################################
 
# calculate predicted SGLT2 benefit (absolute risk reduction = ARR):
# pARR = S0(t)^HR - S0(t)
trial_hr_kf_GLP1 <- 0.81

cohort <- cohort %>% 
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100,
         # added benefit of GLP1 on top of SGLT2
         `ckdpc_40egfr_survival_SGLT2 + GLP1`=ckdpc_40egfr_survival^trial_hr_kf_GLP1,
         `ckdpc_40egfr_SGLT2 + GLP1_benefit`=`ckdpc_40egfr_survival_SGLT2 + GLP1` - ckdpc_40egfr_survival,
  )

# The FLOW trial (N Engl J Med. 2024 Jul 11;391(2):109-121.) randomised participants to GLP1 and found a significant reduction kidney disease progression
# inclusion criteria were eGFR 20-50 and uACR 11.3-565mg/mmol OR eGFR 50-75 and uACR 33.9-565mg/mmol

# In this example we will compare a treatment strategy based on the above (without upper limits for uACR) to a strategy based on pARR targeting a similar proportion of the population

cohort <- cohort %>% 
  mutate(
    flow_criteria = ifelse(
      preegfr < 50 & uacr >= 11.3, T, ifelse(
        preegfr >= 50 & uacr >= 33.9, T, F
      )
    )
  )

# the main question will be which people will benefit from adding GLP1 to SGLT2, which is why we will evaluate the pARR cutoff for SGLT2 + GLP1
n = "SGLT2 + GLP1"
benefit_variable = paste0("ckdpc_40egfr_", n, "_benefit")

# get pARR threshold targeting same proportion of population as FLOW criteria 
cutoff <- cohort[[benefit_variable]] %>% quantile(
  1 -
    (nrow(cohort[cohort$flow_criteria == T,]) / 
       nrow(cohort))
) %>% as.numeric()
  
cutoff_equivalent_40egfr_score <- cohort  %>% .$ckdpc_40egfr_score %>% quantile(
  1 -
    (nrow(cohort[cohort$flow_criteria == T,]) / 
       nrow(cohort))
) %>% as.numeric()

cohort <- cohort %>% mutate(
  treat_model = ifelse(`ckdpc_40egfr_SGLT2 + GLP1_benefit` > cutoff, T, F),        # pARR strategy (matching FLOW criteria treatment proportion)
  recommendation_group = ifelse(flow_criteria == F,
                                ifelse(treat_model == F, "FN_MN", "FN_MY"),
                                ifelse(treat_model == F, "FY_MN", "FY_MY"))
)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv
vars <- c(vars, studydrug_var)
factors <- c(factors, studydrug_var)
table <- CreateTableOne(vars = vars, strata = "recommendation_group", data = cohort,
                         factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)

write.csv2(tabforprint, file = paste0(today, "_baseline_table_by_parr.csv"))

############################2 NUMBERS TREATED AND EVENTS AVOIDED################################################################


## estimate projected impacted of different treatment strategies over 3 years 

#Treat none
#estimated number of events with no one treated
events.notx <- round(nrow(cohort)*mean(cohort$ckdpc_40egfr_score/100)) 

#estimated number if events if we were to treat all these people with SGLT2 + GLP1
cohort <- cohort %>% mutate(`ckdpc_40egfr_score.SGLT2 + GLP1.all` = 
                                     100*(1-cohort$`ckdpc_40egfr_survival_SGLT2 + GLP1`)
)

`events.SGLT2 + GLP1.all` <- round(
  nrow(cohort) * mean(cohort$`ckdpc_40egfr_score.SGLT2 + GLP1.all`/100)
)

#estimated number of events if we were to treat SGLT2 as per guidelines and added GLP1 only to those per FLOW criteria
cohort <- cohort %>% mutate(`ckdpc_40egfr_score.SGLT2 + GLP1.flow` = 
                              ifelse(flow_criteria == F, ckdpc_40egfr_score,
                                     100*(1-cohort$`ckdpc_40egfr_survival_SGLT2 + GLP1`))
)

`events.SGLT2 + GLP1.flow` <- round(
  nrow(cohort) * mean(cohort$`ckdpc_40egfr_score.SGLT2 + GLP1.flow`/100)
)

#estimated number of events if we were to treat SGLT2 as per guidelines and added GLP1 based on pARR
cohort <- cohort %>% mutate(`ckdpc_40egfr_score.SGLT2 + GLP1.model` = 
                              ifelse(treat_model == F, ckdpc_40egfr_score,
                                     100*(1-cohort$`ckdpc_40egfr_survival_SGLT2 + GLP1`))
)

`events.SGLT2 + GLP1.model` <- round(
  nrow(cohort) * mean(cohort$`ckdpc_40egfr_score.SGLT2 + GLP1.model`/100)
)

# collate in dataframe
predicted_events <- data.frame(
  strategy = "treat none",
  n_treated = 0, 
  perc_treated = 0, 
  events = round(events.notx/n.imp), 
  events_avoided = round((events.notx - events.notx)/n.imp), 
  events_avoided_perc = round((events.notx - events.notx) / events.notx * 100, 1), 
  nnt = Inf
)

predicted_events <- predicted_events %>% rbind(
  data.frame(
    strategy = "treat all",
    n_treated = round(nrow(cohort)/n.imp),
    perc_treated = round(nrow(cohort)/nrow(cohort) * 100, 1),
    events = round(`events.SGLT2 + GLP1.all`/n.imp),
    events_avoided = round((`events.SGLT2 + GLP1.all` - events.notx)/n.imp),
    events_avoided_perc = round((`events.SGLT2 + GLP1.all` - events.notx) / events.notx * 100, 1),
    nnt = round(1/
                  (mean(cohort$ckdpc_40egfr_score/100) - 
                     mean(cohort$`ckdpc_40egfr_score.SGLT2 + GLP1.all`/100)))
    
  )
)

predicted_events <- predicted_events %>% rbind(
  data.frame(
    strategy = "FLOW criteria",
    n_treated = round(nrow(cohort %>% filter(flow_criteria == T))/n.imp),
    perc_treated = round(nrow(cohort %>% filter(flow_criteria == T))/nrow(cohort) * 100, 1),
    events = round(`events.SGLT2 + GLP1.flow`/n.imp),
    events_avoided = round((`events.SGLT2 + GLP1.flow` - events.notx)/n.imp),
    events_avoided_perc = round((`events.SGLT2 + GLP1.flow` - events.notx) / events.notx * 100, 1),
    nnt = round(1/
                  (mean(cohort[cohort$flow_criteria == T,]$ckdpc_40egfr_score/100) - 
                     mean(cohort[cohort$flow_criteria == T,]$`ckdpc_40egfr_score.SGLT2 + GLP1.flow`/100)))
    
  )
)

predicted_events <- predicted_events %>% rbind(
  data.frame(
    strategy = "Risk score",
    n_treated = round(nrow(cohort %>% filter(treat_model == T))/n.imp),
    perc_treated = round(nrow(cohort %>% filter(treat_model == T))/nrow(cohort) * 100, 1),
    events = round(`events.SGLT2 + GLP1.model`/n.imp),
    events_avoided = round((`events.SGLT2 + GLP1.model` - events.notx)/n.imp),
    events_avoided_perc = round((`events.SGLT2 + GLP1.model` - events.notx) / events.notx * 100, 1),
    nnt = round(1/
                  (mean(cohort[cohort$treat_model == T,]$ckdpc_40egfr_score/100) - 
                     mean(cohort[cohort$treat_model == T,]$`ckdpc_40egfr_score.SGLT2 + GLP1.model`/100)))
    
  )
)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
write.csv2(predicted_events, file = paste0(today, "_nnt_table_by_strategy.csv"))

############################3 OBSERVED BENEFIT DURING FOLLOW UP TO 5-YEARS################################################################

## compare outcomes in 5-year extended observational analyses of treatment strategy based on FLOW criteria vs pARR

# formula for use in weighted cox models:
ps.formula <- formula(paste0("studydrug", main, " ~ ", paste(covariates, collapse=" + ")))


## in order for the weights to be valid by subgroup, re-calculate weights within each subgruop

cohort$overlap_within_recommendation_group <- NA

for (p in levels(as.factor(cohort$recommendation_group))) {
  print(paste0("Overlap weights within group for ", p))
  
  overlap <- SumStat(ps.formula=ps.formula, data = as.data.frame(cohort[cohort$recommendation_group == p,]))
  
  cohort[cohort$recommendation_group == p,]$overlap_within_recommendation_group <- overlap$ps.weights$overlap
  
}


options(datadist=NULL)

# Initialize a list to store pooled results
pooled_results <- list()

# Loop through imputations
for (i in 1:n.imp) {
  cohort_imp <- cohort %>% filter(.imp == i)
  
  for (p in levels(as.factor(cohort_imp$recommendation_group))) {
    print(p)
    
    for (n in levels(as.factor(cohort_imp[[studydrug_var]]))) {
      df_name <- paste0("df_", n)
      
      fit <- survfit(Surv(ckd_egfr40_5y_censtime_yrs, ckd_egfr40_5y_censvar) ~ 1, 
                     data = cohort_imp %>% filter(recommendation_group == p & !!sym(studydrug_var) == n), 
                     weights = cohort_imp %>% filter(recommendation_group == p & !!sym(studydrug_var) == n) %>% .$overlap_within_recommendation_group)
      
      df <- data.frame(time = fit$time, surv = fit$surv, std.err = fit$std.err)
      df <- rbind(data.frame(time = 0, surv = 1, std.err = 0), df) # Add row for t=0
      
      assign(df_name, df)
    }
    
    df_diff_name <- paste0("df_diff_", p)
    
    df_diff <- merge(`df_SGLT2 + DPP4/SU`, `df_SGLT2 + GLP1`, by = "time", suffixes = c("_SGLT2 + DPP4/SU", "_SGLT2 + GLP1"))
    df_diff <- df_diff %>% mutate(
      difference = ifelse(time == 0, 0, `surv_SGLT2 + DPP4/SU` - `surv_SGLT2 + GLP1`),
      se_difference = sqrt(`std.err_SGLT2 + DPP4/SU`^2 + `std.err_SGLT2 + GLP1`^2),
      lower_ci = ifelse(time == 0, 0, difference - 1.96 * se_difference),
      upper_ci = ifelse(time == 0, 0, difference + 1.96 * se_difference),
      recommendation_group = p
    )
    
    df_diff$.imp <- i  # Add imputation identifier
    pooled_results[[paste0(p, "_", i)]] <- df_diff
  }
}

# Combine all results across imputations
pooled_data <- bind_rows(pooled_results)

# Pool the results by averaging across imputations
pooled_summary <- pooled_data %>%
  group_by(recommendation_group, time) %>%
  summarise(
    pooled_difference = mean(difference, na.rm = TRUE),
    pooled_se = sqrt(mean(se_difference^2, na.rm = TRUE)),
    pooled_lower_ci = pooled_difference - 1.96 * pooled_se,
    pooled_upper_ci = pooled_difference + 1.96 * pooled_se
  ) %>%
  mutate(
    recommendation_group_label = case_when(
      recommendation_group == "FN_MN" ~ paste0(recommendation_group,
        "(n=", 
        if (nrow(cohort %>% filter(recommendation_group == "FN_MN")) %% 10 == 5) {round(nrow(cohort %>% filter(recommendation_group == "FN_MN")) / n.imp) - 1} 
        else {round(nrow(cohort %>% filter(recommendation_group == "FN_MN")) / n.imp)}, 
        ")"
      ),
      recommendation_group == "FY_MN" ~ paste0(recommendation_group,
        "(n=", 
        if ( nrow(cohort %>% filter(recommendation_group == "FY_MN")) %% 10 == 5) {round( nrow(cohort %>% filter(recommendation_group == "FY_MN")) / n.imp) - 1} 
        else {round( nrow(cohort %>% filter(recommendation_group == "FY_MN")) / n.imp)}, 
        ")"
      ),
      recommendation_group == "FY_MY" ~ paste0(recommendation_group,
        "(n=", 
        round( nrow(cohort %>% filter(recommendation_group == "FY_MY")) / n.imp), 
        ")"
      ),
      recommendation_group == "FN_MY" ~ paste0(recommendation_group,
        "(n=", 
        round( nrow(cohort %>% filter(recommendation_group == "FN_MY")) / n.imp), 
        ")"
      )
    )
  )

# Split data by recommendation_group
pooled_summary_split <- pooled_summary %>% group_split(recommendation_group)

# Apply loess smoothing for each subgroup
smoothed_results <- lapply(pooled_summary_split, function(data) {
  loess_diff <- loess(pooled_difference ~ time, data = data, span = 0.1)
  loess_lower <- loess(pooled_lower_ci ~ time, data = data, span = 0.1)
  loess_upper <- loess(pooled_upper_ci ~ time, data = data, span = 0.1)
  
  data %>%
    mutate(
      smoothed_diff = predict(loess_diff, newdata = data),
      smoothed_lower_ci = predict(loess_lower, newdata = data),
      smoothed_upper_ci = predict(loess_upper, newdata = data)
    )
})

# Combine the smoothed results back into one dataframe
df_combined <- bind_rows(smoothed_results)

## cumulative absolute risk reduction plot by recommendation group:
ci_plot <- df_combined %>% ggplot(aes(x = time, y = smoothed_diff*100, color = recommendation_group_label)) +
  geom_line(size = 1) +
  # geom_ribbon(aes(ymin = smoothed_lower_ci*100, ymax = smoothed_upper_ci*100, fill = recommendation_group_label), alpha = 0.2, color = NA) +
  labs(#title = "Kidney protection benefit by SGLT2i treatment recommendation",
    x = "Time (years)",
    y = "Observed ARR in kidney disease progression") +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#D55E00", "#E69F00", "grey40", "grey15")) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "grey40", "grey15")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(legend.position = c(0.44, 0.875), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5), # General axis line style
        
        # Remove top and right axes lines
        axis.line.x.top = element_blank(),    # No line on the top
        axis.line.y.right = element_blank(),
        panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0,3.5), xlim = c(0, 4.8))


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_arr_by_treatment_recommendation.tiff"), width=7.5, height=6, units = "in", res=800)
ci_plot
dev.off()


############################4 DECISION CURVE ANALYSIS################################################################

# define absolute risk of outcome
dca_data <- cohort %>% mutate(ckdpc_40egfr_score_risk = ckdpc_40egfr_score/100)

dca_results <- dca(Surv(ckd_egfr40_censtime_yrs, ckd_egfr40_censvar) ~ flow_criteria + ckdpc_40egfr_score_risk,
                thresholds = seq(0, 0.20, by = 0.001),
                time = 3,
                data=dca_data)

# get true positive / negative rate for example
threshold_data <- dca_results$dca %>%
  filter(threshold == round(cutoff_equivalent_40egfr_score/100, 1))

# difference in true positive rate for pARR strategy and albuminuria threshold, per 100,000
(threshold_data$tp_rate[4] - threshold_data$tp_rate[3])*100000

# difference in true negative rate for pARR strategy and albuminuria threshold, per 100,000
((1-threshold_data$fp_rate[4]) - (1-threshold_data$fp_rate[3]))*100000


p_dca <- as_tibble(dca_results) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label, linetype = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", 
              span = 0.2, size = 1.25) +
  coord_cartesian(ylim = c(-0.00105984276971715, 0.0105984276971715
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Risk tolerance\n(of 3-year absolute risk of kidney disease progression)", y = "Net utility", color = "") +
  theme_minimal() +
  scale_color_manual(
    values = c(
      "black", "grey40", "#0072B2", "#E69F00"
    ),
    labels = c("Treat all with SGLT2 + GLP1", 
               "Treat none with GLP1 (with SGLT2 only)", 
               "Treat with SGLT2 + GLP1 according to FLOW criteria", 
               "Treat with SGLT2 + GLP1 according to risk score\n(threshold varying by risk tolerance)")) + 
  scale_linetype_manual(
    values = c(
      "solid", 
      "solid", 
      "solid", 
      "solid"),
    labels = c("Treat all with SGLT2 + GLP1", 
               "Treat none with GLP1 (with SGLT2 only)", 
               "Treat with SGLT2 + GLP1 according to FLOW criteria", 
               "Treat with SGLT2 + GLP1 according to risk score\n(threshold varying by risk tolerance)")) + 
  guides(color = guide_legend("Treatment strategy"), 
         linetype = guide_legend("Treatment strategy")) +
  theme(legend.position = c(0.65, 0.8),
        panel.border = element_blank(),
        
        # Add custom axis lines
        axis.line = element_line(color = "black", size = 0.5), # General axis line style
        
        # Remove top and right axes lines
        axis.line.x.top = element_blank(),    # No line on the top
        axis.line.y.right = element_blank(),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0.0012,0.12), ylim = c(0,0.015))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_decision_curve_analysis.tiff"), width=6, height=5, units = "in", res=800) 
print(p_dca)
dev.off()
