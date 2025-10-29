## SENSITIVITY ANALYSIS INCLUDING INDIVIDUALS WITH MISSING BASELINE UACR (WITH IMPUTED VALUES)


########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/sens/")
source("00 Setup_sens.R")


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data.Rda"))

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end


# for computational speed, keep minimal dataset only for each studydrug variable
# first rename full dataset
temp_all <- temp

for (m in 2) {
  
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
for (m in 2) {
  
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
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  

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
for (m in 2) {
  
  
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
continuous_vars <- setdiff(vars, factor_vars)
normal_vars <- setdiff(continuous_vars, nonnormal)

for (m in 2) {
  
  print(paste0("Making weighted table for variable studydrug", m))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  studydrug_var <- paste0("studydrug", m)
  weights_overlap <- paste0("overlap", m)
  treat_var <- sym(studydrug_var)
  weight_var <- sym(weights_overlap)
  
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

############################4 HAZARD RATIOS OVERALL################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs <- data.frame()


# main dataset is large - for speed of computation we will only load in dataset we need each time
gc()


# calculate hazard ratios

# for every studydrug variable:
for (m in 2) {
  
  print(paste0("Loading data for variable studydrug", m, collapse = ""))
  
  # define studydrug variable and weights variables to be used
  studydrug_var = paste0("studydrug", m)
  weights_overlap = paste0("overlap", m)
  weights_iptw = paste0("IPTW", m)
  
  # load data
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  gc()
  
  # define which drugs are evaluated with current studydrug variable
  drug_levels <- levels(cohort[[studydrug_var]])
  
  
  # remove double overlapping entries (take one only)
  cohort <- cohort %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  cohort_all <- cohort
  
  # calculate hazard ratios per outcome
  for (k in outcomes) {
    
    cohort <- cohort_all
    
    if (k == "retinopathy") {
      # for outcome retinopathy, only retain those without a history of this
      cohort <- cohort %>% filter(predrug_retinopathy == F)
      weights_overlap = paste0(weights_overlap, "_", k)
      weights_iptw = paste0(weights_iptw, "_",k) 
      
    }
    
    if (k == "acutepancreatitis") {
      # for outcome pancreatitis, only retain those without a history of this
      cohort <- cohort %>% filter(predrug_pancreatitis == F)
      weights_overlap = paste0(weights_overlap, "_", k)
      weights_iptw = paste0(weights_iptw, "_",k) 
      
    }
    
    
    analysis_approaches <-  c("unadj", "adj", "ow", "iptw")
    
    
    # for other studydrug variables use regular censoring variables
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    
    
    print(paste0("Calculating event numbers per drug level for outcome ", k))
    
    # calculate number of subjects in each group
    count <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var)) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", m, "}_count"),
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var)) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", m, "}_followup"),
                  values_from=time)
    
    # summarise number of events per group
    events <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var)) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(!!sym(studydrug_var), events) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", m, "}_events"),
                  values_from=events)
    
    
    # write formulas for adjusted and unadjusted analyses
    f2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, collapse = ""))
    
    f_adjusted2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, " + ", paste(covariates, collapse=" + "), collapse = ""))
    
    # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
    
    for (n in drug_levels[-1]) {
      for (c in c("COEFS", "SE")) {
        for (d in analysis_approaches) {
          assign(paste0(c, ".", n, ".", d), rep(NA, n.imp))
          
        }
      }
    }
    
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      #overlap weighted analyses
      fit.ow <- coxph(f_adjusted2, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])

      #store coefficients and standard errors from this model
      for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coefficients to extract
        
        drug_name <- drug_levels[n+1]
        
        # for every analysis approach, extract coefficients
        for (d in analysis_approaches) {
          
          # get model for this analysis approach
          model <- get(paste0("fit.", d))
          
          # write commands as strings that will dynamically extract coefficient variables for this analysis approach + drug level:
          
          # coef_vector[i] <- model$coefficients[n]
          coef_statement <- paste0("`COEFS.", drug_name, ".", d, "`[", i, "] <- model$coefficients[", n, "]", collapse = "")
          # se_vector[i] <- sqrt(model$var[n,n])
          se_statement <- paste0("`SE.", drug_name, ".", d, "`[", i, "] <- sqrt(model$var[", n, ",", n, "])", collapse = "")
          
          # execute commands
          eval(str2lang(coef_statement))
          eval(str2lang(se_statement))
          
          rm(model)
          
        }
      }
      
    }
    
    
    ## loop to pool and store results
    for (n in 1:length(drug_levels)) {
      
      # get drug name
      drug_name <- drug_levels[n]
      
      
      for (d in analysis_approaches) {
        if (n == 1) {
          # if drug level is reference category then HR will be NA
          pooled_hr <- c(NA, NA, NA)
          pooled_hr_string <- NA
          
        } else {
          
          # define names for objects containing pooled hr + 95% ci (as vector and as string)
          pooled_hr_name <- paste0(d, "_", drug_name, "_hr")
          pooled_hr_string_name <- paste0(d, "_", drug_name, "_string")
          
          # define names for coefficient vectors
          coef_name <- paste0("COEFS.", drug_name, ".", d)
          se_name <- paste0("SE.", drug_name, ".", d)
          
          coef_vector <- get(coef_name)
          se_vector <- get(se_name)
          
          # assign appropriate name to pooled hr + 95% ci (as vector)
          assign(pooled_hr_name, pool.rubin.HR(coef_vector, se_vector, n.imp))
          # get dynamic handle to object
          pooled_hr <- get(pooled_hr_name)
          
          # assign appropriate name to pooled hr + 95% ci (as string)
          assign(pooled_hr_string_name, paste0(sprintf("%.2f", round(pooled_hr[1], 2)), " (", sprintf("%.2f", round(pooled_hr[2], 2)), ", ", sprintf("%.2f", round(pooled_hr[3], 2)), ")"))
          # get dynamic handle to object
          pooled_hr_string <- get(pooled_hr_string_name)
          
          
        }
        
        # create dataframe containing events, follow-up etc
        outcome_hr <- data.frame(outcome = k, 
                                 count = as.numeric(count[n]),
                                 followup = as.numeric(followup[n]), 
                                 events = as.character(events[n]),
                                 contrast = paste0(drug_name, " vs ", drug_levels[1], collapse = ""),
                                 variable = paste0("studydrug", m, collapse = ""),
                                 analysis = d,
                                 HR = pooled_hr[1],
                                 LB = pooled_hr[2],
                                 UB = pooled_hr[3],
                                 string = pooled_hr_string)
        
        # combine results by each analysis approach within each drug type (within studydrug variable)
        hrs <- rbind(hrs, outcome_hr)
        
        
      }
      
    }
    
  }    
  
  
}

rm(cohort_all)

# ensure hazard ratio and 95% ci are stored as numeric variables
class(hrs$HR) <- class(hrs$LB) <- class(hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
hrs <- hrs %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

# store hazard ratios
hrs_incl_missing_uacr <- hrs
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs_incl_missing_uacr, file=paste0(other_day, "_hrs_incl_missing_uacr.Rda"))
