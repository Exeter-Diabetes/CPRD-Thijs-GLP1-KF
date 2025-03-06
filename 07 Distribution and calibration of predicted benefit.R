# in this script we will compare predicted with observed absolute risk reductions and create main figures
# we will also compare treatment strategies based on Predicted 3-year absolute risk reductions with the current albuminuria treshold

########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_GLP1_imputed_data_with_observed_surv.Rda"))

m = 2
studydrug_var = paste0("studydrug", m)
weights_overlap = paste0("overlap", m)

############################1 DISTRIBUTION/CALIBRATION of pARR################################################################

# calculate predicted SGLT2 benefit (absolute risk reduction = ARR):
# pARR = S0(t)^HR - S0(t)
trial_hr_kf_SGLT2 <- 0.62
trial_hr_kf_GLP1 <- 0.81

cohort <- cohort %>% 
  mutate(ckdpc_50egfr_survival=(100-ckdpc_50egfr_score)/100,
         ckdpc_50egfr_survival_SGLT2=ckdpc_50egfr_survival^trial_hr_kf_SGLT2,
         ckdpc_50egfr_SGLT2_benefit=ckdpc_50egfr_survival_SGLT2 - ckdpc_50egfr_survival,
         ckdpc_50egfr_survival_GLP1=ckdpc_50egfr_survival^trial_hr_kf_GLP1,
         ckdpc_50egfr_GLP1_benefit=ckdpc_50egfr_survival_GLP1 - ckdpc_50egfr_survival,
         # added benefit of GLP1 on top of SGLT2
         `ckdpc_50egfr_survival_GLP1 + SGLT2`=ckdpc_50egfr_survival_SGLT2^trial_hr_kf_GLP1,
         `ckdpc_50egfr_GLP1 + SGLT2_benefit`=`ckdpc_50egfr_survival_GLP1 + SGLT2` - ckdpc_50egfr_survival_SGLT2,
         #create variabels for egfr and albuminuria categories
         egfr_cat = ifelse(preegfr < 45, "20-45", ifelse(preegfr < 60, "45-60", "≥60")),
         egfr_cat = factor(egfr_cat),
         albuminuria_cat = ifelse(uacr >30, "≥30", ifelse(uacr > 3, "3-30", "<3")),
         albuminuria_cat = factor(albuminuria_cat),
         ckd_cat = paste(egfr_cat, albuminuria_cat, sep = "_"),
         ckd_cat = factor(ckd_cat))

## FIGURE: histogram of predicted benefit

histogram_list <- list()

proportional_heights <- table(cohort$ckd_cat) / sum(table(cohort$ckd_cat)) 
# make separate plots for different ckd categories

max_x_value = 0

for (c in levels(cohort$ckd_cat)) {

  data <- cohort %>% 
    filter(ckd_cat == c) %>% 
    select(contains("SGLT2_benefit")) %>%
    pivot_longer(cols = everything(), 
                 names_to = "studydrug", 
                 values_to = "benefit") %>%
    mutate(studydrug = case_when(
      grepl("SGLT2", studydrug, ignore.case = TRUE) ~ "SGLT2",
      grepl("GLP1", studydrug, ignore.case = TRUE) ~ "GLP1 + SGLT2",
      TRUE ~ NA_character_)) %>%
    select(benefit, studydrug)

benefit_histogram <- ggplot(data %>%
                                    mutate(predicted_benefit_percent = benefit * 100), 
                                  aes(x = predicted_benefit_percent)) +
  geom_histogram(data = filter(data, studydrug == "SGLT2"),
                 aes(y = ..count.. / n.imp, fill = studydrug),
                 binwidth = 0.02, color = "black") +  
  geom_histogram(data = filter(data, studydrug == "GLP1/SGLT2"),
                 aes(y = ..count.. / n.imp, fill = studydrug),
                 binwidth = 0.02, color = "black") +  
  scale_fill_manual(values = cols, labels = names(cols)) +
  # geom_vline(xintercept = cutoff1*100, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 2.5, y = proportional_heights[c]*900, label = paste0("eGFR ", egfr_value, "mL/min/1.73m2 and albuminuria ", albuminuria_value, "mg/mmol (", round(100*nrow(cohort %>% filter(ckd_cat == c))/nrow(cohort), 1), "%, n=",nrow(cohort %>% filter(ckd_cat == c))/n.imp,")"), vjust = 0, hjust = 1, angle = 0, size = 4, color = "black") +
  labs(x = "", y = "Frequency") +
  theme_bw() +
  theme(
    # Remove panel border (default borders around the plot area)
    panel.border = element_blank(),
    
    # Add custom axis lines
    axis.line = element_line(color = "black", size = 0.5), # General axis line style
    
    # Remove top and right axes lines
    axis.line.x.top = element_blank(),    # No line on the top
    axis.line.y.right = element_blank(),   # No line on the right
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.position = c(0.725, 0.25),
        legend.title = element_blank())

max_x_value <- pmax(max_x_value, max(data_benefit_percent))

histogram_list[[c]] <- benefit_histogram

}

for (c in levels(cohort$ckd_cat)) {
  # ensure all histograms have same x-axis range and proportional heights to their size
  histogram_list[[c]] <- histogram_list[[c]] + 
    coord_cartesian(xlim=c(0,max_x_value), ylim=c(0,proportional_heights[c]*1000), expand = F)

}

# combine all plots in one
combined_histogram <- wrap_plots(histogram_list, ncol = 1) + 
  plot_layout(heights = proportional_heights)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_predicted_benefit_histogram.tiff"), width=6, height=12, units = "in", res=800) 
combined_histogram
dev.off()



##FIGURE: calibration plot of predicted vs observed absolute risk reductions
benefit_calibration_plot_list <- list()

for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  benefit_variable = paste0("ckdpc_50egfr_", n, "_benefit")
  
  obs_v_pred_for_plot <- cohort %>%
    # group predicted benefit by decile
    mutate(benefit_decile = ntile(!!sym(benefit_variable), n.quantiles)) %>%
    group_by(benefit_decile) %>%
    summarise(median_predicted_benefit=median(!!sym(benefit_variable), na.rm=T),
              mean_predicted_benefit=mean(!!sym(benefit_variable), na.rm=T),
              mean_benefit=mean(survdiff_ckd_egfr50),
              se_benefit=mean(se_survdiff_ckd_egfr50),
              median_benefit=median(survdiff_ckd_egfr50),
              lq_benefit=quantile(survdiff_ckd_egfr50, prob=c(.25)),
              uq_benefit=quantile(survdiff_ckd_egfr50, prob=c(.75)),
              upper_ci=mean_benefit + 1.96*se_benefit,
              lower_ci=mean_benefit - 1.96*se_benefit)
  
  empty_tick <- data.frame(matrix(NA, nrow = 1, ncol = length(obs_v_pred_for_plot)))
  names(empty_tick) <- names(obs_v_pred_for_plot)
  empty_tick <- empty_tick %>%
    mutate(benefit_decile=0)z
  
  p_benefit_calibration <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=median_predicted_benefit*100)) +
    geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100, color= cols[n]),width=0.1,size=1) +
    geom_point(aes(y = median_benefit*100, color=cols[n]), shape=18, size=3) +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    theme_bw() +
    xlab("Model-predicted absolute risk reduction (%)") + ylab("Counterfactual absolute risk reduction (%)\nestimated from observed data")+
    scale_colour_manual(values = cols[n]) +
    theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
          plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
    theme(axis.text=element_text(size=rel(1.5)),
          axis.title=element_text(size=rel(1.5)),
          plot.title=element_text(hjust = 0.5),
          plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
          legend.position = "none") +
    coord_cartesian(xlim = c(0,3.5), ylim = c(-.1,3.56))
  
  benefit_calibration_plot_list[[n]] <- p_benefit_calibration
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_predicted_benefit_", n, "_calibration.tiff"), width=6, height=5.5, units = "in", res=800) 
  print(benefit_calibration_plot_list[[n]])
  dev.off()
  
}





for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  contrast_drug = ifelse(n == "GLP1 + SGLT2", "SGLT2", "DPP4 + SU")
  survdiff_var = paste0("survdiff_ckd_egfr50_", n, "_", contrast_drug)
  benefit_var = paste0("ckdpc_50egfr_", n, "_benefit")
  
  # calculate calibration slope and Brier score
  slope <- slope_se <- brier <- brier_se <- rep(NA, n.imp)
  
  for (i in 1:n.imp) {
    print(paste("Imputation ", i))
    
    data <- cohort %>% filter(.imp == i)
    
    # calculate calibration slope
    x <- lm(!!sym(survdiff_var) ~ !!sym(benefit_var), data = data)
    slope[i] <- coef(x)[2]
    slope_se[i] <- vcov(x)[2,2]
    rm(x)
    
    bootstrap_brier <- rep(NA, n.bootstrap)
    
    for (b in 1:n.bootstrap) {
      # Resample the combined data with replacement
      bootstrap_sample <- data %>% sample_frac(size = 1, replace = TRUE)
      
      # Compute squared errors
      bootstrap_sample <- bootstrap_sample %>%
        mutate(
          squared_error = (!!sym(benefit_var) - !!sym(survdiff_var))^2,
          weighted_error = !!sym(weights_overlap) * squared_error
        )
      
      # Compute weighted Brier score for this bootstrap sample
      bootstrap_brier[b] <- sum(bootstrap_sample$weighted_error, na.rm = TRUE) / sum(bootstrap_sample[[weights_overlap]], na.rm = TRUE)
      rm(bootstrap_sample)
    }
    
    # Calculate mean Brier score
    brier[i] <- mean(bootstrap_brier)
    
    # Calculate standard error (SE)
    brier_se[i] <- sd(bootstrap_brier)
    rm(data)
  }
  
  slope_mean <- mean(slope)
  B <- var(slope)            # Between-imputation variance
  W <- mean(slope_se^2)      # Within-imputation variance
  TV <- W + (1 + 1 / n.imp) * B  # Total variance
  # Confidence interval for pooled slope
  slope_lc <- slope_mean - 1.96 * sqrt(TV)
  slope_uc <- slope_mean + 1.96 * sqrt(TV)
  print(paste0("Calibration slope for pARR (", n, "): ", round(slope_mean, 3), ", 95% CI ", round(slope_lc, 3), "-", round(slope_uc,3)))
  
  #pool and print brier score
  brier_se_pooled <- sqrt(mean(brier_se^2) + (1+1/n.imp)*var(brier))
  print(paste0("Brier score for pARR  (", n, ")", mean(brier), ", 95% CI ", mean(brier)-1.96*brier_se_pooled, "-", mean(brier)+1.96*brier_se_pooled))
  
}