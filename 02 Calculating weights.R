########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data.Rda"))

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end

# weight variables for each studydrug variable
n.studydrug.vars <- sum(grepl("studydrug", colnames(temp)))

# for computational speed, keep minimal dataset only for each studydrug variable
# first rename full dataset
temp_all <- temp

for (m in 1:n.studydrug.vars) {
  
  studydrug_var = paste0("studydrug", m)
  
  temp <- temp_all %>% 
    filter(.imp != 0) %>% 
    select("patid", ".imp", contains("studydrug"), all_of(covariates), "dstartdate", "predrug_pancreatitis", "predrug_retinopathy", "with_hes", "ncurrtx2", "weight_pct_change")
  
  # create empty variables for weights
  temp[[paste0("IPTW", m, collapse = "")]] <-  # IPTW
    temp[[paste0("overlap", m, collapse = "")]] <- # overlap weights
    temp[[paste0("ps", m, collapse = "")]] <- # raw propensity scores
    # weights for subsets without prior retinopathy / pancreatitis
    temp[[paste0("IPTW", m, "_retinopathy", collapse = "")]] <- temp[[paste0("overlap", m, "_retinopathy", collapse = "")]] <- NA
  
  # remove duplicate drug episodes
  temp <- temp %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(temp, file=paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", m, ".Rda"))
  rm(temp)
}  
rm(temp_all)
gc()

# calculate weights
for (m in 1:n.studydrug.vars) {
  
  studydrug_var = paste0("studydrug", m)
  
  print(paste0("Processing data for variable studydrug", m, collapse = ""))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", m, ".Rda"))
  
  gc()
  
  ps.formula <- formula(paste0("studydrug", m, " ~ ", paste(covariates, collapse=" + ")))
  
  # force temp to be data.frame() for SumStat function
  temp <- temp %>% as.data.frame(temp)
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp[temp$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 2nd and 98th percentile
    w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       ifelse(
                                         w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         w.overlap$ps.weights$IPW
                                       ))
    
    
    ## Add weights to data frame
    weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
    temp[temp$.imp == i,][[paste0("IPTW", m, collapse = "")]]  <- weights$IPW
    temp[temp$.imp == i,][[paste0("overlap", m, collapse = "")]] <- weights$overlap

    # save propensity score additionally for diagnostics
    temp[temp$.imp == i,][[paste0("ps", m, collapse = "")]] <- w.overlap$propensity[,2]
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  ## check distribution of propensity score before and after overlap weighting 
  studydrug_var <- paste0("studydrug", m)
  weights_overlap <- paste0("overlap", m)
  ps_var <- paste0("ps", m)
  
  # prep data for plot
  weight_plot_data <- temp %>% 
    group_by(.imp, !!sym(studydrug_var)) %>%
    mutate(weight_before = 1 / n(), # overlap weights are proportional to total number in group; standardise unweighted observations similarly to get proportion
           weight_after  = !!sym(weights_overlap)) %>%
    ungroup() %>%
    select(weight_before, weight_after, !!sym(ps_var), !!sym(studydrug_var)) %>%
    pivot_longer(cols = starts_with("weight_"),
                 names_to = "stage",
                 values_to = "weight") %>%
    mutate(stage = recode(stage,
                          weight_after  = "After weighting",
                          weight_before = "Before weighting"),
           stage = as.factor(stage),
           stage = relevel(stage, ref = "Before weighting"))
  
  # histogram with propensity score by treatment arm before and after weighting
  weight_plot <- ggplot(weight_plot_data, aes(x = !!sym(ps_var), fill = !!sym(studydrug_var))) +
    geom_histogram(data = weight_plot_data %>% filter(stage == "Before weighting"),
                   aes(weight = weight), alpha = 0.5, colour = "grey20", bins = 30, position = "identity") +
    geom_histogram(data = weight_plot_data %>% filter(stage == "After weighting"),
                   aes(weight = weight), alpha = 0.5, colour = "grey20", bins = 30, position = "identity") +
    scale_fill_manual(values = cols) +
    facet_wrap(stage ~ ., nrow = 1, scales = "free_x") +
    labs(x = "Propensity score", y = "Percent")+
    theme_bw() +
    theme(legend.title=element_blank(),
          legend.text = element_text(size=10),
          strip.text = element_text(size = 12),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  
  #save plot
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  tiff(paste0(today, "_propensity_score_distribution_plot_", m, ".tiff"), width=8, height=5, units = "in", res=800)
  print(weight_plot)
  dev.off()
  
  # calculate separate weights for subset without prior retinopathy / pancreatitis for outcomes pancreatitis and diabetic retinopathy
  
  ## retinopathy
  temp_no_retinopathy <- temp %>% filter(predrug_retinopathy == F)
  temp_retinopathy <- temp %>% filter(predrug_retinopathy == T)
  
  temp_no_retinopathy <- temp_no_retinopathy %>% as.data.frame()
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for those without retinopathy in imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp_no_retinopathy[temp_no_retinopathy$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 2nd and 98th percentile
    w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       ifelse(
                                         w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         w.overlap$ps.weights$IPW
                                       ))
    
    
    ## Add weights to data frame
    weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
    temp_no_retinopathy[temp_no_retinopathy$.imp == i,][[paste0("IPTW", m, "_retinopathy", collapse = "")]]  <- weights$IPW
    temp_no_retinopathy[temp_no_retinopathy$.imp == i,][[paste0("overlap", m, "_retinopathy", collapse = "")]] <- weights$overlap
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  temp <- rbind(temp_retinopathy, temp_no_retinopathy)
  rm(temp_retinopathy)
  rm(temp_no_retinopathy)
  
  ## pancreatitis
  
  temp_no_pancreatitis <- temp %>% filter(predrug_pancreatitis == F)
  temp_pancreatitis <- temp %>% filter(predrug_pancreatitis == T)
  
  temp_no_pancreatitis <- temp_no_pancreatitis %>% as.data.frame()
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for those without pancreatitis in imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp_no_pancreatitis[temp_no_pancreatitis$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 2nd and 98th percentile
    w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                       ifelse(
                                         w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                         w.overlap$ps.weights$IPW
                                       ))
    
    
    ## Add weights to data frame
    weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
    temp_no_pancreatitis[temp_no_pancreatitis$.imp == i,][[paste0("IPTW", m, "_acutepancreatitis", collapse = "")]]  <- weights$IPW
    temp_no_pancreatitis[temp_no_pancreatitis$.imp == i,][[paste0("overlap", m, "_acutepancreatitis", collapse = "")]] <- weights$overlap
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  temp <- rbind(temp_pancreatitis, temp_no_pancreatitis)
  rm(temp_pancreatitis)
  rm(temp_no_pancreatitis)
  

  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(temp, file=paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_withweights_", m, ".Rda"))
  rm(temp)
  gc()
  
}

############################2 SAVE WEIGHTS################################################################

# add weights to imputed dataset with all variables (load this one first)
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data.Rda"))

temp_all <- temp %>% filter(.imp != 0) 
rm(temp)
gc()
for (m in 1:n.studydrug.vars) {
  
  
  studydrug_var = paste0("studydrug", m)
  print(studydrug_var)
  
  # load dataset with weights for merging
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_withweights_", m, ".Rda"))
  
  cohort <- temp %>% select(.imp, patid, !!sym(studydrug_var), 
                                                dstartdate, contains("overlap"), contains("IPTW")) %>% 
                     left_join(temp_all, by = c(".imp", "patid", "dstartdate", studydrug_var))
  
  cohort <- cohort %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    ungroup()
  
  rm(temp)
  
  # save dataset with weights so this can be used in subsequent scripts
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  
}
rm(temp_all)


############################3 WEIGHTED BASELINE TABLE################################################################

# factor and nonnormal variables are already specified as factor_vars and nonnormal
# specify which ones are continuous and normally distributed
continuous_vars <- setdiff(vars, factor_vars)
normal_vars <- setdiff(continuous_vars, nonnormal)

for (m in 1:n.studydrug.vars) {
  
  print(paste0("Making weighted table for variable studydrug", m))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  studydrug_var <- paste0("studydrug", m)
  weights_overlap <- paste0("overlap", m)
  treat_var <- sym(studydrug_var)
  weight_var <- sym(weights_overlap)
  
  # write function to get weighted table
  weighted_summary <- function(data, varname, weight = NULL, group, type = "normal") {
    
    if (is.null(weight)) {
      data <- data %>% mutate(wt = 1)
    } else {
      data <- data %>% mutate(wt = .data[[weight]])
    }
    
    var_sym <- sym(varname)
    
    if (type == "normal") {
      data %>%
        group_by({{group}}) %>%
        summarise(
          mean = weighted.mean(!!var_sym, wt, na.rm = TRUE),
          sd   = sqrt(Hmisc::wtd.var(!!var_sym, weights = wt, na.rm = TRUE)),
          .groups = "drop"
        ) %>%
        mutate(stat = sprintf("%.2f (%.2f)", mean, sd)) %>%
        select({{group}}, stat)
    } else {
      data %>%
        group_by({{group}}) %>%
        summarise(
          med = Hmisc::wtd.quantile(!!var_sym, weights = wt, probs = 0.5, na.rm = TRUE),
          q1  = Hmisc::wtd.quantile(!!var_sym, weights = wt, probs = 0.25, na.rm = TRUE),
          q3  = Hmisc::wtd.quantile(!!var_sym, weights = wt, probs = 0.75, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(stat = sprintf("%.2f [%.2f–%.2f]", med, q1, q3)) %>%
        select({{group}}, stat)
    }
  }
  
  
  weighted_table_cat <- function(data, varname, weight = NULL, group) {
    
    
    if (is.null(weight)) {
      data <- data %>% mutate(wt = 1)
    } else {
      data <- data %>% mutate(wt = .data[[weight]])
    }
    
    var_sym <- sym(varname)
    
    data %>%
      group_by({{group}}, !!var_sym) %>%
      summarise(n = sum(wt, na.rm = TRUE), .groups = "drop") %>%
      group_by({{group}}) %>%
      mutate(percent = 100 * n / sum(n)) %>%
      ungroup() %>%
      mutate(stat = sprintf("%f (%.1f%%)", n, percent)) %>%
      select({{group}}, !!var_sym, stat)
  }
  
  
  # numerical variables
  tbl_num_norm_w <- lapply(normal_vars, function(v)
    weighted_summary(cohort, var = v, weight = weight_var, group = !!treat_var, type = "normal") %>%
      mutate(variable = v, type = "Weighted")
  ) %>% bind_rows()
  
  
  tbl_num_skew_w <- lapply(nonnormal, function(v)
    weighted_summary(cohort, v, weight_var, !!treat_var, type = "nonnormal") %>%
      mutate(variable = v, type = "Weighted")
  ) %>% bind_rows()
  
  tbl_num <- bind_rows(tbl_num_norm_w,
                       tbl_num_skew_w)
  
  # categorical
  
  tbl_cat_w <- lapply(factor_vars, function(v)
    weighted_table_cat(cohort, v, weight_var, !!treat_var) %>%
      mutate(variable = v, type = "Weighted")
  ) %>% bind_rows()
  
  # smd
  
  fml <- reformulate(vars, response = studydrug_var)
  
  bal <- bal.tab(
    x = fml,                 
    data = cohort,
    weights = cohort[[weight_var]], 
    method = "weighting",
    estimand = "ATO",        # "ATO" for overlap, "ATE" for IPTW
    s.d.denom = "pooled",
    quick = TRUE
  )
  
  # Extract balance summary
  bal_df <- bal$Balance %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable")
  
  # Collapse factor levels by variable name (everything before ":")
  col_un <- grep("Diff\\.Un", names(bal_df), value = TRUE)
  col_adj <- grep("Diff\\.Adj", names(bal_df), value = TRUE)
  
  bal_summary <- bal_df %>%
    mutate(variable = sub(":.*", "", variable)) %>%
    group_by(variable) %>%
    summarise(
      SMD_weighted   = mean(.data[[col_adj]], na.rm = TRUE)
    ) %>%
    arrange(desc(abs(SMD_weighted)))
  
  smds <- bal$Balance %>%
    tibble::rownames_to_column("variable") %>%
    select(variable,
           SMD_weighted   = .data[[col_adj]])
  
  # join together
  table_num_fmt <- tbl_num %>%
    pivot_wider(names_from = !!treat_var, values_from = stat) %>%
    select(variable, type, everything())
  
  table_cat_fmt <- tbl_cat_w %>%
    pivot_wider(names_from = !!treat_var, values_from = stat) %>%
    select(variable, type, everything())
  
  table1 <- bind_rows(table_num_fmt, table_cat_fmt) %>%
    mutate(across(all_of(factor_vars), as.character)) %>%
    mutate(factor_level = coalesce(!!!syms(factor_vars))) %>%
    select(-all_of(factor_vars)) 
  
  table1 <- table1 %>%
    # join SMDs
    mutate(
      variable2 = ifelse(
        is.na(factor_level) | factor_level %in% c(0, 1, "TRUE", "FALSE", T, F),
        variable,
        paste0(variable, "_", factor_level)
      )
    ) %>%
    left_join(smds %>% rename(variable2 = variable), by = "variable2") %>%
    arrange(type, variable) %>%
    select(-type, -variable2) %>%
    # keep max SMD per variable
    group_by(variable) %>%
    mutate(SMD_weighted = max(SMD_weighted, na.rm = TRUE)) %>%
    ungroup() %>%
    # order vars in correct order of table for display
    mutate(variable = factor(variable, levels = vars)) %>%
    arrange(variable)
  
  
  for (q in levels(cohort[[studydrug_var]])) {
    
    table1 <- table1 %>% 
      # multiply weighted proportions with treatment group size and adjust for n.imp
      mutate(
        extracted_value = as.numeric(str_extract(!!sym(q), "^[0-9\\.]+")),
        N   = round(extracted_value / n.imp * nrow(cohort[cohort[[studydrug_var]] == q,])/n.imp, 2),
        pct = ifelse(is.na(factor_level), NA, extracted_value / n.imp * 100)
      ) %>%
      mutate(N = ifelse(is.na(factor_level), !!sym(q), as.numeric(N))) %>%
      select(-extracted_value, -!!sym(q)) %>%
      rename(
        !!q := N,
        !!paste0(q, "_pct") := pct
      )
  }
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  write.table(as.data.frame(table1), file = paste0(today, "_weighted_table_", m, ".csv"), sep = ";", dec = ",", row.names = F)
}
