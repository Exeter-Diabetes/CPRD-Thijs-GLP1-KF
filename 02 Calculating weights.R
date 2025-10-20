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
    select("patid", ".imp", contains("studydrug"), all_of(covariates), "dstartdate", "predrug_pancreatitis", "predrug_retinopathy", "with_hes", "ncurrtx2")
  
  # create empty variables for weights
  temp[[paste0("IPTW", m, collapse = "")]] <- temp[[paste0("overlap", m, collapse = "")]] <- 
    temp[[paste0("IPTW", m, "_retinopathy", collapse = "")]] <- temp[[paste0("overlap", m, "_retinopathy", collapse = "")]] <- NA
  
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
    
    # truncate IPTW at 5th and 95th percentile
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
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  # calculate separate weights for subset without retinopathy for outcomes pancreatitis and diabetic retinopathy
  
  # ---- retinopathy
  temp_no_retinopathy <- temp %>% filter(predrug_retinopathy == F)
  temp_retinopathy <- temp %>% filter(predrug_retinopathy == T)
  
  temp_no_retinopathy <- temp_no_retinopathy %>% as.data.frame()
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for those without retinopathy in imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp_no_retinopathy[temp_no_retinopathy$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 5th and 95th percentile
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
  
  # ----- pancreatitis
  
  temp_no_pancreatitis <- temp %>% filter(predrug_pancreatitis == F)
  temp_pancreatitis <- temp %>% filter(predrug_pancreatitis == T)
  
  temp_no_pancreatitis <- temp_no_pancreatitis %>% as.data.frame()
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for those without pancreatitis in imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp_no_pancreatitis[temp_no_pancreatitis$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 5th and 95th percentile
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
continuous_vars <- setdiff(vars, factors)
normal_vars <- setdiff(continuous_vars, nonnormal)

for (m in 1:n.studydrug.vars) {
  
  print(paste0("Making weighted table for variable studydrug", m))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  studydrug_var <- paste0("studydrug", m)
  weights_overlap <- paste0("overlap", m)
  treat_var <- sym(studydrug_var)
  weight_var <- sym(weights_overlap)
  
  if (m == 1) {
    cohort <- cohort %>% filter(!!sym(studydrug_var) != "SGLT2i + GLP1-RA" ) %>% mutate(!!sym(studydrug_var) := droplevels(!!sym(studydrug_var)))
  }
  
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
  
  tbl_cat_w <- lapply(factors, function(v)
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
    mutate(across(all_of(factors), as.character)) %>%
    mutate(factor_level = coalesce(!!!syms(factors))) %>%
    select(-all_of(factors)) 
  
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
