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
         # create extra variable to differentiate between reduced/preserved eGFR for plots below
           egfr_cat2 = ifelse(preegfr < 60, "<60", "≥60"),
           egfr_cat2 = factor(egfr_cat2),
         albuminuria_cat = ifelse(uacr >30, "≥30", ifelse(uacr > 3, "3-30", "<3")),
         albuminuria_cat = factor(albuminuria_cat, levels = c("<3", "3-30", "≥30")),
         ckd_cat = paste(egfr_cat2, albuminuria_cat, sep = "_"),
         ckd_cat = factor(ckd_cat, levels = c(
           "≥60_<3",
           "≥60_3-30",
           "≥60_≥30",
           "<60_<3",
           "<60_3-30",
           "<60_≥30"
         )))


# print overall benefit per drug and by CKD category
for (n in drug_levels[-1]) {
  
  predicted_benefit_variable = paste0("ckdpc_50egfr_", n, "_benefit")
  
  print(paste0("Overall median pARR for ", n, ": ", sprintf("%.2f", median(cohort[[predicted_benefit_variable]]*100)), "% (IQR ", sprintf("%.2f", quantile(cohort[[predicted_benefit_variable]]*100, 0.25)), "-", sprintf("%.2f", quantile(cohort[[predicted_benefit_variable]]*100, 0.75)), ")"))
  
  for (c in levels(cohort$ckd_cat)) {
    
    print(paste0("Median pARR for ", n, " in CKD category ", c, ": ", sprintf("%.2f", median(cohort[cohort$ckd_cat == c,][[predicted_benefit_variable]]*100)), "% (IQR ", sprintf("%.2f", quantile(cohort[cohort$ckd_cat == c,][[predicted_benefit_variable]]*100, 0.25)), "-", sprintf("%.2f", quantile(cohort[cohort$ckd_cat == c,][[predicted_benefit_variable]]*100, 0.75)), ")"))
  }
}

############################2 DISTRIBUTION OF PREDICTED BENEFIT################################################################

histogram_list <- list()

# make separate plots for different ckd categories


for (e in levels(cohort$egfr_cat2)) {
  
  max_x_value = 0
  
  proportional_heights <- table(cohort[cohort$egfr_cat2==e,]$albuminuria_cat) / sum(table(cohort[cohort$egfr_cat2==e,]$albuminuria_cat)) 
  
  
  for (a in levels(cohort$albuminuria_cat)) {
    
    data <- cohort %>% 
      filter(egfr_cat2 == e) %>%       
      filter(albuminuria_cat == a) %>% 
      select(contains("SGLT2_benefit")) %>%
      pivot_longer(cols = everything(), 
                   names_to = "studydrug", 
                   values_to = "benefit") %>%
      mutate(studydrug = case_when(
        grepl("GLP1", studydrug, ignore.case = TRUE) ~ "GLP1 + SGLT2",      
        grepl("SGLT2", studydrug, ignore.case = TRUE) ~ "SGLT2",
        TRUE ~ NA_character_),
        studydrug = factor(studydrug)) %>%
      select(benefit, studydrug) %>%
      mutate(predicted_benefit_percent = benefit * 100)
    
    max_x_value <- pmax(max_x_value, max(data$predicted_benefit_percent))
    
    
    benefit_histogram <- ggplot(data=data, aes(x = predicted_benefit_percent)) 
    
    
    histogram_list[[a]] <- benefit_histogram
    
  }
  
  for (a in levels(cohort$albuminuria_cat)) {
    
    data <- cohort %>% 
      filter(egfr_cat2 == e) %>%       
      filter(albuminuria_cat == a) %>% 
      select(contains("SGLT2_benefit")) %>%
      pivot_longer(cols = everything(), 
                   names_to = "studydrug", 
                   values_to = "benefit") %>%
      mutate(studydrug = case_when(
        grepl("GLP1", studydrug, ignore.case = TRUE) ~ "GLP1 + SGLT2",      
        grepl("SGLT2", studydrug, ignore.case = TRUE) ~ "SGLT2",
        TRUE ~ NA_character_),
        studydrug = factor(studydrug)) %>%
      select(benefit, studydrug) %>%
      mutate(predicted_benefit_percent = benefit * 100)
    
    # ensure all histograms have same x-axis range and proportional heights to their size
    histogram_list[[a]] <- histogram_list[[a]] +
      geom_histogram(data = filter(data, studydrug == "SGLT2"),
                     aes(y = ..count.. / n.imp, fill = studydrug),
                     binwidth = max_x_value/50, color = "black") +  
      geom_histogram(data = filter(data, studydrug == "GLP1 + SGLT2"),
                     aes(y = ..count.. / n.imp, fill = studydrug),
                     binwidth = max_x_value/50, color = "black", alpha = 0.5) +  
      scale_fill_manual(values = cols[rev(levels(data$studydrug))], labels = names(cols[rev(levels(data$studydrug))])) +
      # geom_vline(xintercept = cutoff1*100, linetype = "dashed", color = "black", size = 1) +
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
      guides(fill = if (a == levels(cohort$albuminuria_cat)[1]) {guide_legend(reverse = TRUE)} else {"none"}) +
      theme(legend.position = c(0.725, 0.25),
            legend.title = element_blank()) + 
      annotate("text", x = max_x_value-1, y = proportional_heights[a]*nrow(cohort[cohort$egfr_cat2 == e,])*.75*.5, label = paste0("eGFR ", e, "mL/min/1.73m2 and albuminuria ", a, "mg/mmol\n(", round(100*nrow(cohort %>% filter(albuminuria_cat == a & egfr_cat2 == e))/nrow(cohort %>% filter(egfr_cat2 == e)), 1), "%, n=",nrow(cohort %>% filter(albuminuria_cat == a & egfr_cat2 == e))/n.imp,")"), vjust = 0, hjust = 1, angle = 0, size = 4, color = "black") +
      coord_cartesian(xlim=c(0,max_x_value), ylim=c(0,proportional_heights[a]*nrow(cohort[cohort$egfr_cat2 == e,])*.5), expand = F)
    
  }
  
  # combine all plots in one
  combined_histogram <- wrap_plots(histogram_list, ncol = 1) + 
    plot_layout(heights = proportional_heights)
  
  
  egfr_value <- gsub("<60", "reduced", e)
  egfr_value <- gsub("≥60", "preserved", egfr_value)

  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_predicted_benefit_histogram_", egfr_value, "_egfr.tiff"), width=9, height=15, units = "in", res=800) 
  print(combined_histogram)
  dev.off()
}


############################3 CALIBRATION OF PREDICTED BENEFIT################################################################


for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  benefit_calibration_plot_list <- list()
  
  predicted_benefit_variable = paste0("ckdpc_50egfr_", n, "_benefit")
  
  if (n == levels(cohort[[studydrug_var]])[4]) {
    reference = "SGLT2"
  } else {
    reference = "DPP4 + SU"
  }
  counterfactual_benefit_variable = paste0("survdiff_", n, "_", reference, "_ckd_egfr50")
  counterfactual_benefit_variable_se = paste0("se_survdiff_", n, "_", reference, "_ckd_egfr50")
  
  max_x_value = 0
  max_y_value = 0
  
  for (p in 1:length(levels(cohort$egfr_cat2))) {
    
    q = levels(cohort$egfr_cat2)[p]
    
    obs_v_pred_for_plot <- cohort %>% 
      filter(egfr_cat2 == q) %>%
      # group predicted benefit by decile
      mutate(benefit_decile = ntile(!!sym(predicted_benefit_variable), n.quantiles)) %>%
      group_by(benefit_decile) %>%
      summarise(median_predicted_benefit=median(!!sym(predicted_benefit_variable), na.rm=T),
                mean_predicted_benefit=mean(!!sym(predicted_benefit_variable), na.rm=T),
                mean_benefit=mean(!!sym(counterfactual_benefit_variable)),
                se_benefit=mean(!!sym(counterfactual_benefit_variable_se)),
                median_benefit=median(!!sym(counterfactual_benefit_variable)),
                lq_benefit=quantile(!!sym(counterfactual_benefit_variable), prob=c(.25)),
                uq_benefit=quantile(!!sym(counterfactual_benefit_variable), prob=c(.75)),
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
            plot.title = element_text(size = rel(1.5), face = "bold")) + theme(plot.margin = margin()) +
      theme(axis.text=element_text(size=rel(1.5)),
            axis.title=element_text(size=rel(1.5)),
            plot.title=element_text(hjust = 0.5),
            plot.subtitle=element_text(hjust = 0.5,size=rel(1.2)),
            legend.position = "none") +
      # coord_cartesian(xlim = c(0,3.5), ylim = c(-.1,3.56)) +
      ggtitle(if (q == "<60") {"Reduced eGFR (<60mL/min/1.73m2)"} else {"Preserved eGFR (≥60mL/min/1.73m2)"})
    

    
    benefit_calibration_plot_list[[p]] <- p_benefit_calibration
    
  }
  
  for (p in 1:length(levels(cohort$egfr_cat2))) {
    
    benefit_calibration_plot_list[[p]] <- benefit_calibration_plot_list[[p]] +
      coord_cartesian(xlim = c(0,max_x_value), ylim = c(-.1,max_y_value))
    
    }
  
  # wrap calibration plot for reduced and preserved eGFR together side by side
  combined_calibration_plot <- wrap_plots(benefit_calibration_plot_list, ncol = nlevels(cohort$egfr_cat2)) 
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_predicted_benefit_calibration_", n, ".tiff"), width=18, height=5.5, units = "in", res=800) 
  print(combined_calibration_plot)
  dev.off()
  
}


############################4 CALIBRATION SLOPE OF PREDICTED BENEFIT################################################################



for (n in levels(cohort[[studydrug_var]])[-1]) {
  
  print(paste0("Calibration slope for ", n, " benefit"))
  
  contrast_drug = ifelse(n == "GLP1 + SGLT2", "SGLT2", "DPP4 + SU")
  survdiff_var = paste0("`survdiff_", n, "_", contrast_drug, "_ckd_egfr50`")
  benefit_var = paste0("`ckdpc_50egfr_", n, "_benefit`")
  
  # calculate calibration slope and Brier score
  slope <- slope_se <- brier <- brier_se <- rep(NA, n.imp)
  
  for (i in 1:n.imp) {
    print(paste("Imputation ", i))
    
    data <- cohort %>% filter(.imp == i)
    
    # calculate calibration slope
    x <- lm(as.formula(paste0(survdiff_var, " ~ ", benefit_var)), data = data)
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
          squared_error = (eval(str2lang(benefit_var)) - eval(str2lang(survdiff_var)))^2,
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
