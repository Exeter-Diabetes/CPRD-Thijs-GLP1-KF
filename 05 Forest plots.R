########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_hrs.Rda"))
load(paste0(today, "_factor_hrs.Rda"))

# set default studydrug variable
studydrug_var = paste0("studydrug", main)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

drug_reference = levels(cohort[[studydrug_var]])[1]
drug_of_interest = levels(cohort[[studydrug_var]])[2]

############################1 FOREST PLOT FOR MAIN OUTCOMES################################################################

# GLP1-RA trial meta-analysis HR: 0.79, 95% CI 0.66-0.95 (Circulation. 2024 Nov 26;150(22):1781-1790.)

# other GLP1-RA trial meta-analysis HR: HR 0·81, 95% CI 0·72–0·92 (Lancet Diabetes Endocrinol. 2025 Jan;13(1):15-28.)

hrs <- hrs %>% filter(!grepl("interaction", variable))


n.studydrug.vars <- hrs %>% .$variable %>% as.factor() %>% nlevels()

# main outcomes
for (m in 1:n.studydrug.vars) {
  
  temp <- hrs
  
  if (m == 1) {
    temp <- temp %>%
      filter(!str_detect(contrast, "GLP1"))
  }
  
  # Collect datasets across outcomes
  all_outcomes <- lapply(outcomes_per_drugclass, function(k) {
    temp <- temp %>%
      filter(analysis=="ow",
             outcome == k,
             variable == paste0("studydrug", m)) %>%
      mutate(HR = ifelse(is.na(string), 1, HR),
             string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)),
             analysis_label = case_when(analysis == "ow" ~ ""),
             outcome = k)
    
    temp[c("outcome","contrast","variable","analysis")] <- 
      lapply(temp[c("outcome","contrast","variable","analysis")], factor)
    
    temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))
    
    # Add a heading row per outcome
    heading <- tibble(
      outcome = k,
      analysis_label = "",
      drug = k,
      HR = NA, LB = NA, UB = NA,
      string = ifelse(k == outcomes_per_drugclass[1], 
                      "Hazard ratio (95% CI)", 
                      ""), 
      nN = ifelse(k == outcomes_per_drugclass[1], 
                  "Events/subjects", 
                  "") 
    )
    
    bind_rows(heading, temp)
  })
  
  plot_df <- bind_rows(all_outcomes)
  
  # Define pretty labels
  outcome_labels <- c(
    "ckd_egfr40" = "≥40% eGFR decline/ESKD",
    "ckd_egfr50" = "≥50% eGFR decline/ESKD",
    "mace"       = "MACE (incl. CV death)",
    "hf"         = "Hospitalisation for HF",
    "death"      = "All-cause mortality"
  )
  
  xmin = 0.15
  xmax = 2.5
  # Apply mapping before building plot_df
  plot_df <- plot_df %>%
    mutate(outcome = recode(outcome, !!!outcome_labels),
           # also update headings if you used outcome as heading text
           drug = ifelse(drug %in% outcomes_per_drugclass, outcome, drug))
  
  # Set y order (reverse for top-down order)
  plot_df$y_order <- rev(seq_len(nrow(plot_df)))
  
  # Forest plot with tabulated extras
  forest_plot <- ggplot(plot_df, aes(y = y_order)) +
    # Forest CI + point
    geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
    geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    # Outcome headings on the left
    geom_text(aes(x = 0.15, label = drug,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 0) +
    
    # Events/subjects in middle-left
    geom_text(aes(x = 0.475, label = nN,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 1) +
    
    # HR text on right-hand side
    geom_text(aes(x = 1.7, label = string,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 0) +
    
    scale_x_continuous(trans = "log10",
                       breaks = c(0.5, 1, 1.5),
                       limits = c(xmin, xmax)) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = "Hazard Ratio (95% CI)", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(
            hjust = (log10(1.0) - log10(xmin)) / (log10(xmax) - log10(xmin))
          )
    )
  
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_HR_", m, "_main_outcomes.tiff"), width=10, height=length(outcomes_per_drugclass)*1.2, units = "in", res=800)
  print(forest_plot)
  dev.off()
}

############################2 FOREST PLOT FOR SAFETY AND CONTROL OUTCOMES################################################################

for (m in 1:n.studydrug.vars) {
  
  temp <- hrs
  
  if (m == 1) {
    temp <- temp %>%
      filter(!str_detect(contrast, "GLP1"))
  }
  
  # Collect datasets across outcomes
  all_outcomes <- lapply(safety_outcomes, function(k) {
    temp <- temp %>%
      filter(analysis=="ow",
             outcome == k,
             variable == paste0("studydrug", m)) %>%
      mutate(HR = ifelse(is.na(string), 1, HR),
             string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)),
             analysis_label = case_when(analysis == "ow" ~ ""),
             outcome = k)
    
    temp[c("outcome","contrast","variable","analysis")] <- 
      lapply(temp[c("outcome","contrast","variable","analysis")], factor)
    
    temp$drug <- as.factor(sub(" vs.*", "", temp$contrast))
    
    # Add a heading row per outcome
    heading <- tibble(
      outcome = k,
      analysis_label = "",
      drug = k,
      HR = NA, LB = NA, UB = NA,
      string = ifelse(k == safety_outcomes[1], 
                      "Hazard ratio (95% CI)", 
                      ""), 
      nN = ifelse(k == safety_outcomes[1], 
                  "Events/subjects", 
                  "") 
    )
    
    bind_rows(heading, temp)
  })
  
  plot_df <- bind_rows(all_outcomes)
  
  # Define pretty labels
  outcome_labels <- c(
    "acutepancreatitis"      = "Acute pancreatitis",
    "retinopathy"            = "Diabetic retinopathy",
    "lowerlimbfracture"      = "Lower limb fracture",
    "upperlimbfracture"      = "Upper limb fracture"
  )
  
  xmin = 0.15
  xmax = 2.5
  # Apply mapping before building plot_df
  plot_df <- plot_df %>%
    mutate(outcome = recode(outcome, !!!outcome_labels),
           # also update headings if you used outcome as heading text
           drug = ifelse(drug %in% safety_outcomes, outcome, drug))
  
  # Set y order (reverse for top-down order)
  plot_df$y_order <- rev(seq_len(nrow(plot_df)))
  
  # Forest plot with tabulated extras
  forest_plot <- ggplot(plot_df, aes(y = y_order)) +
    # Forest CI + point
    geom_linerange(aes(xmin = LB, xmax = UB), na.rm = TRUE) +
    geom_point(aes(x = HR), shape = 15, size = 3, na.rm = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    # Outcome headings on the left
    geom_text(aes(x = 0.15, label = drug,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 0) +
    
    # Events/subjects in middle-left
    geom_text(aes(x = 0.475, label = nN,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 1) +
    
    # HR text on right-hand side
    geom_text(aes(x = 1.7, label = string,
                  fontface = ifelse(HR %>% is.na(), "bold", "plain")),
              hjust = 0) +
    
    scale_x_continuous(trans = "log10",
                       breaks = c(0.5, 1, 1.5),
                       limits = c(xmin, xmax)) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = "Hazard Ratio (95% CI)", y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(
            hjust = (log10(1.0) - log10(xmin)) / (log10(xmax) - log10(xmin))
          )
    )
  
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_HR_", m, "_safety_outcomes.tiff"), width=10, height=length(safety_outcomes)*1.2, units = "in", res=800)
  print(forest_plot)
  dev.off()
}


############################3 FOREST PLOTs FOR HR BY COMORBIDITIES################################################################



factors <- c("malesex", "white_ethnicity", "predrug_cvd", "predrug_heartfailure")

factor_hrs$condition_label <- case_when(
  factor_hrs$factor == "malesex" & factor_hrs$condition == TRUE ~ "Sex: male",
  factor_hrs$factor == "malesex" & factor_hrs$condition == FALSE ~ "Sex: female",
  
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == TRUE ~ "Ethnicity: white",
  factor_hrs$factor == "white_ethnicity" & factor_hrs$condition == FALSE ~ "Ethnicity: non-white",
  
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == TRUE ~ "ASCVD: present",
  factor_hrs$factor == "predrug_cvd" & factor_hrs$condition == FALSE ~ "ASCVD: absent",
  
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == TRUE ~ "Heart failure: present",
  factor_hrs$factor == "predrug_heartfailure" & factor_hrs$condition == FALSE ~ "Heart failure: absent",
  
  TRUE ~ NA_character_
) 

factor_hrs <- factor_hrs %>% mutate(
  condition_label = factor(condition_label, levels = c(
    "ASCVD: absent", "ASCVD: present", "Heart failure: absent", "Heart failure: present", 
    "Sex: male", "Sex: female", "Ethnicity: white", "Ethnicity: non-white"
  ))
)

factor_hrs_all <- factor_hrs


for (k in outcomes_per_drugclass) {
  
  # get counts per category added: (so totals added but by condition instead):
  factor_hrs <- factor_hrs_all %>%
    filter(outcome == k & contrast %in% c(paste0(drug_reference, " vs ", drug_reference), paste0(drug_of_interest, " vs ", drug_reference)))
  
  
  test <- factor_hrs %>% group_by(factor, condition) %>% summarise(events_number = sum(as.numeric(events_number)), count = sum(as.numeric(count))) %>%
    
    mutate(
      `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
    )
  
  factor_hrs <- factor_hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference)) %>% select(-nN, -count, -events_number) %>% left_join(
    test, by = c("factor", "condition")
  )
  
  HR <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$HR)
  LB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$LB)
  UB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$UB)
  
  
  
  temp_hrs <- factor_hrs %>%
    filter(factor %in% factors) %>%
    filter(analysis=="ow") %>%
    filter(outcome == k) %>%
    filter(variable == paste0("studydrug", main)) %>%
    filter(contrast == paste0(drug_of_interest, " vs ", drug_reference)) %>%
    mutate(HR = ifelse(is.na(string), 1, HR),
           string = ifelse(is.na(string), "1.00 (ref.)", as.character(string)),
    )
  
  
  
  temp_hrs[c("outcome", "contrast", "variable", "analysis", "condition_label")] <- lapply(temp_hrs[c("outcome", "contrast", "variable", "analysis", "condition_label")], factor)
  
  temp_hrs$drug <- as.factor(sub(" vs.*", "", temp_hrs$contrast))
  # 
  # set factor reference and order of levels
  temp_hrs <- temp_hrs %>% mutate(
    drug = factor(drug, levels = c(drug_reference))
  )
  
  # Prepare labels
  labels_plot <- temp_hrs
  
  
  labels <- data.frame(matrix("", nrow = 1, ncol = length(labels_plot)))
  names(labels) <- names(labels_plot)
  labels <- labels %>% mutate(string = "Hazard Ratio (95% CI)",
                              nN = "Events/subjects",
                              drug = "Treatment arm"
  )
  
  for (n in levels(temp_hrs$drug)) {
    
    j = paste0(drug_of_interest, " vs ", drug_reference)
    
    labels_temp <- labels %>% mutate(
      condition_label = paste0(" Overall"),
      contrast = j,
      HR = !!HR,
      LB = !!LB,
      UB = !!UB,
      string = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", LB), ", ", sprintf("%.2f", UB), ")"),
    )
    
    labels_temp <- labels_temp %>%
      mutate(
        condition_label = fct_relevel(condition_label, paste0(" Overall"))
      )
    
    
    labels_plot <- rbind(labels_temp, labels_plot)
  }
  
  # have to coerce HR and CI to class numeric as they sometimes default to character
  
  class(labels_plot$HR) <- class(labels_plot$LB) <- class(labels_plot$UB) <- "numeric"
  
  
  
  plot_expression <- ""
  
  for (n in rev(levels(temp_hrs$drug))) {
    
    j = paste0(drug_of_interest, " vs ", drug_reference)
    
    labels_temp <- labels_plot %>% filter(contrast == j)
    
    p_counts <- labels_temp %>%
      mutate(., condition_label = fct_rev(fct_relevel(condition_label, paste0(" Overall")))) %>%
      ggplot(aes(y = (fct_relevel(condition_label, paste0(" Overall"), after = Inf)))) +
      geom_text(aes(x = 1, label = nN), hjust = 1,
                colour = ifelse(labels_temp$nN == labels$nN,
                                "white", "black"),
                fontface = ifelse(labels_temp$nN == labels$nN,
                                  "bold", "plain")) +
      annotate("text", x = 1, hjust = 1,
               y = length(levels(droplevels(labels_temp$condition_label))) + 1,
               label = ifelse(j==levels(as.factor(temp_hrs$contrast))[1], "Events/subjects", ""),
               fontface = "bold") +
      theme_void() +
      coord_cartesian(xlim = c(-2, 3), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(levels(droplevels(labels_temp$condition_label))) + 1, length(levels(labels_temp$condition_label)))))
    
    p_hr <- labels_temp %>%
      mutate(., condition_label = fct_rev(fct_relevel(condition_label, paste0(" Overall")))) %>%
      
      ggplot(aes(y = condition_label)) +
      
      scale_x_continuous(trans = "log10", breaks = c(0.15, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)) +
      coord_cartesian(ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(levels(droplevels(labels_temp$condition_label))) + 1, length(levels(labels_temp$condition_label)))),
                      xlim=c(0.15, 3.0)) +
      theme_classic() +
      geom_point(aes(x=HR), shape=15, size=3) +
      geom_linerange(aes(xmin=LB, xmax=UB)) +
      geom_vline(xintercept = 1, linetype="dashed") +
      annotate("text", x = .5,
               y = length(levels(droplevels(labels_temp$condition_label))) + 1,
               label = ifelse(n==levels(temp_hrs$drug)[1], paste0("Favours ", drug_of_interest), "")) +
      annotate("text", x = 2.0,
               y = length(levels(droplevels(labels_temp$condition_label))) + 1,
               label = ifelse(n==levels(temp_hrs$drug)[1], paste0("Favours ", drug_reference, collapse = ""), "")) +
      labs(x=ifelse(n==levels(temp_hrs$drug)[nlevels(as.factor(temp_hrs$drug))], "Hazard ratio", ""), y="") +
      theme(axis.line.y = element_blank(),
            axis.ticks.y= element_blank(),
            axis.text.y= element_blank(),
            axis.title.y= element_blank(),
            axis.line.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
            axis.text.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
            axis.ticks.x = if (n!=levels(temp_hrs$drug)[nlevels(temp_hrs$drug)]) {element_blank()},
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    p_left <- labels_temp %>%
      mutate(., condition_label = fct_rev(fct_relevel(condition_label, paste0(" Overall")))) %>%
      
      ggplot(aes(y = (fct_relevel(condition_label, paste0(" Overall"), after = Inf)))) +
      
      geom_text(
        aes(x = 1, label = condition_label),
        hjust = 0,
        fontface = ifelse(labels_temp$
                            condition_label %in% paste0(" Overall"), "bold", "plain"),
        colour = "black"
        #colour = ifelse(n==levels(temp_hrs$drug)[1] & !labels_temp$condition_label %in% paste0(" ", labels_temp$drug), "white", "black")
      ) +
      theme_void() +
      coord_cartesian(xlim = c(0, 4), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(levels(droplevels(labels_temp$condition_label))) + 1, length(levels(labels_temp$condition_label)))))
    
    p_right <- labels_temp %>%
      mutate(., condition_label = fct_rev(fct_relevel(condition_label, paste0(" Overall")))) %>%
      
      ggplot(aes(y = condition_label)) +
      
      geom_text(
        aes(x = 0, label = string),
        hjust = 0,
        fontface = ifelse(labels_temp$string == "Hazard Ratio (95% CI)", "bold", "plain"),
        colour = ifelse(labels_temp$string == "Hazard Ratio (95% CI)",
                        "white", "black")) +
      annotate("text", x = 0, hjust = 0,
               y = length(levels(droplevels(labels_temp$condition_label))) + 1,
               label = ifelse(n==levels(temp_hrs$drug)[1], "Hazard Ratio (95% CI)", ""),
               fontface = "bold") +
      theme_void() +
      coord_cartesian(xlim = c(0, 4), ylim=c(1,ifelse(n==levels(temp_hrs$drug)[1], length(levels(droplevels(labels_temp$condition_label))) + 1, length(levels(labels_temp$condition_label)))))
    
    assign(paste0("p_counts_", n), p_counts)
    assign(paste0("p_hr_", n), p_hr)
    assign(paste0("p_left_", n), p_left)
    assign(paste0("p_right_", n), p_right)
    
    plot_expression <- paste0(plot_expression, "`p_counts_", n, "` + `p_left_", n, "` + `p_hr_", n, "` + `p_right_", n, "` + ")
  }
  
  n.plots <- nlevels(as.factor(temp_hrs$contrast))
  
  # layout for plots below
  
  i <- 1
  
  height_all_plots <- 24
  height_first_plot <- 24
  
  layout <- paste("area(t = ",(i-1)*height_all_plots, ", l = 3, b = ",(i-1)*height_all_plots+height_first_plot,", r = 6),
                    area(t = ",(i-1)*height_all_plots, ", l = 0, b = ",(i-1)*height_all_plots+height_first_plot,", r = 3),
                    area(t = ",(i-1)*height_all_plots, ", l = 6, b = ",(i-1)*height_all_plots+height_first_plot,", r = 16),
                    area(t = ",(i-1)*height_all_plots, ", l = 17, b = ", (i-1)*height_all_plots+height_first_plot,", r = 19)")
  
  
  layout <- paste0("c(", layout, ")")
  
  layout <- eval(str2lang(layout))
  
  # Final plot arrangement
  
  final_plot_expression <- paste0(plot_expression, "plot_layout(design = layout)")
  
  plot_for_saving <- eval(str2lang(final_plot_expression))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_HR_for_", k, "_as_factors.tiff"), width=14, height=6, units = "in", res=800)
  print(plot_for_saving)
  dev.off()
}


############################4 SPLINE PLOTS################################################################

cohort <- cohort %>% 
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup()

ddist <- cohort %>% datadist()
options(datadist = "ddist") 


numerical_covariates <- c("dstartdate_age", "prebmi", "preegfr", "uacr", "prehba1c", "ckdpc_40egfr_score")

var_labels <- c(
  dstartdate_age = "Age (years)",
  prehba1c = "HbA1c (mmol/mol)",
  prebmi = "BMI (kg/m²)",
  preegfr = "eGFR (mL/min/1.73m²)",
  uacr = "uACR (mg/mmol)",
  presbp = "Systolic blood pressure (mmHg)",
  ckdpc_40egfr_score = "CKD-PC risk score"
)



for (k in outcomes_per_drugclass) {
  print(paste0("Spline plots for outcome ", k))
  
  HR <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$HR)
  LB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$LB)
  UB <- (hrs %>% filter(contrast == paste0(drug_of_interest, " vs ", drug_reference) & outcome == k & analysis == "ow") %>% .$UB)
  
  
  for (q in numerical_covariates) {
    
    r <- var_labels[[q]]
    
    # optimal number of knots is 3 for all variables
    
    # print(paste0("Calculating optimal number of knots for variable ", r))
    # # Define the range of knots to test
    # k_range <- 3:5
    # 
    # # Initialize empty vectors to store results
    # bic_values <- numeric(length(k_range))
    # 
    # # Loop over each value of k
    # for (i in seq_along(k_range)) {
    #   k <- k_range[i]
    #   
    #   # Fit the model with k knots
    #   model <- cph(
    #     as.formula(paste0(
    #       "Surv(ckd_egfr50_censtime_yrs, ckd_egfr50_censvar) ~ studydrug2*rcs(", q, ",", k, ") + ",
    #       paste(setdiff(covariates, q), collapse=" + ") 
    #     )),
    #     data = cohort %>% filter(.imp == n.imp), x = TRUE, y = TRUE
    #   )
    #   
    #   # Store the BIC values
    #   bic_values[i] <- BIC(model)
    # }
    # 
    # # Find the optimal k based on minimum BIC
    # optimal_k_bic <- k_range[which.min(bic_values)]
    # 
    # print(paste0("Optimal number of knots for variable ", r, ": ", optimal_k_bic))
    
    optimal_k_bic = 3
    
    print(paste0("Spline plot for variable ", r))
    
    contrast_dataframe <- data.frame()
    
    for (i in 1:n.imp) {
      
      print(paste0("Imputation number ", i))
      
      # fit model with optimal number of knots
      final_model <- cph(
        as.formula(paste0(
          "Surv(", k, "_censtime_yrs, ", k, "_censvar) ~ studydrug2*rcs(", q, ",", optimal_k_bic, ") + ",
          paste(setdiff(covariates, q), collapse=" + "))
        ),
        data = cohort %>% filter(.imp == i), x = TRUE, y = TRUE
      )
      
      
      anova(final_model)
      p_value_non_linear <- anova(final_model)[2,3] # p value for non-linear interaction term
      
      print(paste0("p-value for non-linear interaction term with ", r, ": ", p_value_non_linear))
      
      q_vals <- seq(
        quantile(cohort[[q]], 0.025, na.rm = TRUE),
        quantile(cohort[[q]], 0.975, na.rm = TRUE),
        by = 0.05
      )
      
      contrast_spline <- contrast(
        final_model, 
        setNames(list("SGLT2 + GLP1", q_vals), c("studydrug2", q)),
        setNames(list("SGLT2 + DPP4/SU",  q_vals), c("studydrug2", q))
      )
      # extract beta and SE
      contrast_dataframe <- rbind(contrast_dataframe, as.data.frame(contrast_spline[c(q,'Contrast','SE')]) %>% mutate(.imp = i))
      
    }
    
    # pool results
    contrast_spline_df <- contrast_dataframe %>% 
      group_by(!!sym(q)) %>% 
      summarise(mean.coef = mean(Contrast),
                W = mean(SE^2),
                B = var(Contrast),
                T.var = W + (1+1/n.imp)*B,
                se.coef = sqrt(T.var),
                Upper = mean.coef + se.coef*1.96,
                Lower = mean.coef - se.coef*1.96) %>%
      select(-c(W, B, T.var, se.coef)) %>%
      rename(Contrast = mean.coef)
    
    
    # plot
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/output/")
    
    # define scale
    x_scale <- if (q == "uacr") {
      scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30))
    } else {
      scale_x_continuous(breaks = function(x) sort(unique(c(0, pretty(x)))))
    }
    
    # set binwidth for histogram
    range_q <- range(if (q == "uacr") {log10(contrast_spline_df[[q]])} else (contrast_spline_df[[q]]), na.rm = TRUE)
    binwidth <- (range_q[2] - range_q[1]) / 30
    
    p_spline <- ggplot(data=contrast_spline_df, aes(x=.data[[q]], y=exp(Contrast))) +
      geom_line(data=contrast_spline_df,aes(x=.data[[q]], y=exp(Contrast)), size=1) +
      xlab(r) +
      ylab("Hazard ratio") +
      x_scale +
      scale_y_log10(breaks = c(0.25, 0.50, 0.75, 1.0, 1.50, 2.0)) +
      geom_ribbon(data=contrast_spline_df, aes(x=.data[[q]], ymin=exp(Lower), ymax=exp(Upper)), alpha=0.2) +
      geom_hline(yintercept = 1, linetype = "dashed")  +
      geom_hline(aes(yintercept = HR, linetype = "hr", size="hr"), color="#56B4E9")  +
      geom_hline(aes(yintercept = LB, linetype = "hr_95", size="hr_95"), color="#56B4E9")  +
      geom_hline(aes(yintercept = UB, linetype = "hr_95", size="hr_95"), color="#56B4E9")  +
      annotate("text", x = mean(range(contrast_spline_df[[q]])), y = 0.35, 
               label = "Favours SGLT2 + GLP1", 
               size = 5, hjust = 0.5, parse = F) +
      annotate("text", x = mean(range(contrast_spline_df[[q]])), y = 1.5, 
               label = "Favours SGLT2 + DPP4/SU", 
               size = 5, hjust = 0.5, parse = F) +
      theme_bw() +
      theme(text = element_text(size = 18),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="top",
            legend.title = element_text(size=14, face = "italic"),
            legend.text = element_text(face="italic"),
            # Add custom axis lines
            axis.line = element_line(color = "black", size = 0.5), # General axis line style
            
            # Remove top and right axes lines
            axis.line.x.top = element_blank(),    # No line on the top
            axis.line.y.right = element_blank(),) +
      scale_linetype_manual(values = c(hr = "twodash", hr_95 = "twodash"), labels = c(hr = sprintf("%.2f", HR), hr_95 = paste0("95% CI ", sprintf("%.2f", LB), "-", sprintf("%.2f", UB))), name="Overall hazard ratio") +
      scale_size_manual(values = c(hr = 1, hr_95 = 0.5), labels = c(hr = sprintf("%.2f", HR), hr_95 = paste0("95% CI ", sprintf("%.2f", LB), "-", sprintf("%.2f", UB))), name="Overall hazard ratio") +
      coord_cartesian(xlim = c(min(contrast_spline_df[[q]]), max(contrast_spline_df[[q]])), ylim = c(0.25, 2.0), expand = F)
    
    
    histogram <- ggplot(cohort %>% filter(.imp == n.imp), aes(x = .data[[q]])) +
      geom_histogram(binwidth = binwidth, fill = "grey70", color = "black") +
      xlab(NULL) +
      ylab(NULL) +
      theme_bw() +
      theme(
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0, 5, 5, 5)
      ) +
      scale_y_continuous(breaks = function(x) {
        top <- (ceiling(max(x)/ 1000) * 1000)
        mid <- ceiling(top / 2)
        c(0, top)
      }) +
      x_scale +
      coord_cartesian(xlim = c(min(contrast_spline_df[[q]]), max(contrast_spline_df[[q]])), expand = FALSE)
    
    spline_and_histogram <- p_spline / histogram + plot_layout(heights = c(6, 1)) 
    
    tiff(paste0(today, "_HR_for_", k, "_by_", q, ".tiff"), width=10, height=4, units = "in", res=600) 
    print(spline_and_histogram)
    dev.off()
    
  }
}

