# in this script we will compare predicted with counterfactual absolute risk reductions (estimated based on observed data) and create main figures

########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_with_observed_surv.Rda"))

studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)
drug_levels <- levels(cohort[[studydrug_var]])
############################1 CALCULATE PREDICTED BENEFIT################################################################

# calculate predicted benefit (absolute risk reduction = ARR):
# pARR = S0(t)^HR - S0(t)
 
# GLP1-RA trial meta-analysis HR: HR 0·81, 95% CI 0·72–0·92 (Lancet Diabetes Endocrinol. 2025 Jan;13(1):15-28.)
trial_hr_kf_GLP1 <- 0.81

cohort <- cohort %>% 
  mutate(ckdpc_40egfr_survival=(100-ckdpc_40egfr_score)/100,
         ckdpc_40egfr_survival_GLP1=ckdpc_40egfr_survival^trial_hr_kf_GLP1,
         `ckdpc_40egfr_SGLT2 + GLP1_benefit`=ckdpc_40egfr_survival_GLP1 - ckdpc_40egfr_survival,
  )


# print overall benefit per drug and by CKD category
for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  predicted_benefit_var = paste0("ckdpc_40egfr_", n, "_benefit")
  
  print(paste0("Overall median pARR for ", n, ": ", sprintf("%.2f", median(cohort[[predicted_benefit_var]]*100)), "% (IQR ", sprintf("%.2f", quantile(cohort[[predicted_benefit_var]]*100, 0.25)), "-", sprintf("%.2f", quantile(cohort[[predicted_benefit_var]]*100, 0.75)), ")"))
  
}

############################2 DISTRIBUTION OF PREDICTED BENEFIT################################################################


# make separate plots for different ckd categories



risk_histogram <- ggplot(cohort, 
                         aes(x = ckdpc_40egfr_score)) +
  geom_histogram(aes(y = ..count.. / n.imp, fill = as.logical(oha)),#predicted_benefit_percent > cutoff1*100),
                 binwidth = 0.25, color = "black") +  
  scale_fill_manual(values = c("grey")) + 

  labs(x = "Risk score", y = "Frequency") +
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
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold"),
        axis.title = element_text(size = rel(0.9))) + theme(plot.margin = margin()) +
  guides(fill = "none") +
  theme(legend.position = c(0.725, 0.25),
        legend.title = element_blank()) + 
  coord_cartesian(xlim=c(0,32), ylim=c(0,3250), expand = F)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_risk_histogram_grey.tiff"), width=6, height=3, units = "in", res=600) 
print(risk_histogram)
dev.off()

benefit_histogram <- ggplot(cohort, 
                         aes(x = `ckdpc_40egfr_SGLT2 + GLP1_benefit`*100)) +
  geom_histogram(aes(y = ..count.. / n.imp, fill = as.logical(oha)),#predicted_benefit_percent > cutoff1*100),
                 binwidth = 0.06640625, color = "black") +  
  scale_fill_manual(values = c("grey")) + 
  
  labs(x = "Predicted benefit from adding GLP1 treatment", y = "Frequency") +
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
  theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
        plot.title = element_text(size = rel(1.5), face = "bold"),
        axis.title = element_text(size = rel(0.9))) + theme(plot.margin = margin()) +
  guides(fill = "none") +
  theme(legend.position = c(0.725, 0.25),
        legend.title = element_blank()) + 
  coord_cartesian(xlim=c(0,8.5), ylim=c(0,5050), expand = F)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_benefit_histogram_grey.tiff"), width=6, height=3, units = "in", res=600) 
print(benefit_histogram)
dev.off()
############################3 CALIBRATION PLOTS OF PREDICTED BENEFIT################################################################


for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  benefit_calibration_plot_list <- list()
  
  predicted_benefit_var = paste0("ckdpc_40egfr_", n, "_benefit")
  
  reference = levels(cohort[[studydrug_var]])[1]
  
  counterfactual_benefit_var = paste0("survdiff_", n, "_", reference, "_ckd_egfr40")
  counterfactual_benefit_var_se = paste0("se_survdiff_", n, "_", reference, "_ckd_egfr40")
  
  max_x_value = 0
  max_y_value = 0

  
    obs_v_pred_for_plot <- cohort %>% 

      mutate(benefit_decile = ntile(!!sym(predicted_benefit_var), n.quantiles)) %>%
      group_by(benefit_decile) %>%
      summarise(median_predicted_benefit=median(!!sym(predicted_benefit_var), na.rm=T),
                mean_predicted_benefit=mean(!!sym(predicted_benefit_var), na.rm=T),
                mean_benefit=mean(!!sym(counterfactual_benefit_var)),
                se_benefit=mean(!!sym(counterfactual_benefit_var_se)),
                median_benefit=median(!!sym(counterfactual_benefit_var)),
                lq_benefit=quantile(!!sym(counterfactual_benefit_var), prob=c(.25)),
                uq_benefit=quantile(!!sym(counterfactual_benefit_var), prob=c(.75)),
                upper_ci=mean_benefit + 1.96*se_benefit,
                lower_ci=mean_benefit - 1.96*se_benefit,
                group = n)
    
    empty_tick <- data.frame(matrix(NA, nrow = 1, ncol = length(obs_v_pred_for_plot)))
    names(empty_tick) <- names(obs_v_pred_for_plot)
    empty_tick <- empty_tick %>%
      mutate(benefit_decile=0)
    
    max_x_value <- pmax(max_x_value, max(obs_v_pred_for_plot$median_predicted_benefit*100)) * 1.1
    max_y_value <- pmax(max_y_value, max(obs_v_pred_for_plot$uq_benefit*100))    
    
    p_benefit_calibration <- ggplot(data=bind_rows(empty_tick,obs_v_pred_for_plot), aes(x=median_predicted_benefit*100)) +
      geom_errorbar(aes(ymax=uq_benefit*100,ymin=lq_benefit*100, color= group),width=0.015*max_x_value,size=1) +
      geom_point(aes(y = median_benefit*100, color=group), shape=18, size=3) +
      geom_abline(intercept = 0, slope = 1, lty = 2) +
      theme_bw() +
      xlab("Model-predicted absolute risk reduction (%)") + ylab("Counterfactual absolute risk reduction (%)\nestimated from observed data")+
      scale_colour_manual(values = cols[n]) +
      theme(panel.border=element_blank(), panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            axis.line.x=element_line(colour = "black"), axis.line.y=element_line(colour="black"),
            plot.title = element_text(size = rel(1.25), face = "bold")) + theme(plot.margin = margin()) +
      theme(axis.text=element_text(size=rel(1.25)),
            axis.title=element_text(size=rel(1.25)),
            plot.title=element_text(hjust = 0.5),
            plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
            legend.position = "none")  +

      coord_cartesian(xlim = c(0,max_x_value), ylim = c(-.1,max_y_value))


  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_predicted_benefit_calibration_", n, ".tiff"), width=5.5, height=6, units = "in", res=800) 
  print(p_benefit_calibration)
  dev.off()
  
}


############################4 CALIBRATION SLOPE OF PREDICTED BENEFIT################################################################



for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  print(paste0("Calibration slope for ", n, " benefit"))
  
  contrast_drug = levels(cohort[[studydrug_var]])[1]
  counterfactual_benefit_var = paste0("`survdiff_", n, "_", contrast_drug, "_ckd_egfr40`")
  predicted_benefit_var = paste0("`ckdpc_40egfr_", n, "_benefit`")
  
  # calculate calibration slope and Brier score
  slope <- slope_se <- brier <- brier_se <- rep(NA, n.imp)
  
  for (i in 1:n.imp) {
    print(paste("Imputation ", i))
    
    data <- cohort %>% filter(.imp == i)
    
    # calculate calibration slope
    x <- lm(as.formula(paste0(counterfactual_benefit_var, " ~ ", predicted_benefit_var)), data = data)
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
          squared_error = (eval(str2lang(predicted_benefit_var)) - eval(str2lang(counterfactual_benefit_var)))^2,
          weighted_error = eval(str2lang(weights_overlap)) * squared_error
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


model_metrics <- data.frame(
  model = "Benefit model",
  calibration_slope = round(slope_mean, 4),
  calibration_slope_ci = paste0(round(slope_lc, 4), "-", round(slope_uc, 4)),
  brier_score = mean(brier),
  brier_score_ci = paste0(mean(brier)-1.96*brier_se_pooled, "-", mean(brier)+1.96*brier_se_pooled)
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
write.csv2(model_metrics, file=paste0(today, "_benefit_model_performance.csv"))
