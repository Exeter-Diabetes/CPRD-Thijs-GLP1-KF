########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

############################1 CALCULATE WEIGHTS AND HAZARD RATIOS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in themain, then stack them again at the end
gc()

# calculate weights
for (n in 1:4) {
  

  studydrug_var = paste0("studydrug", main)
  
  print(paste0("Processing data for variable studydrug", main, collapse = ""))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", main, ".Rda"))
  
    if (n == 1) {
    print(paste0("Remove multiple episodes per subject"))
      #subjects can theoretically be in the non-GLP1 group, and later on have a 2nd episode in the GLP1 group
      temp <- temp %>%
        group_by(.imp, patid) %>%
        filter(
          # Keep "SGLT2 + DPP4/SU" if patient also has "SGLT2 + GLP1"
          !(studydrug2 == "SGLT2 + GLP1" & "SGLT2 + DPP4/SU" %in% studydrug2)
        ) %>%
        ungroup()
      q = "single_episodes"
    }
  
  if (n == 2) {
    print(paste0("Remove subjects on insulin"))
    temp <- temp %>% filter(INS == F)
    q = "no_insulin"
  }
  
  if (n == 3) {
    print(paste0("Remove subjects on more than 4 current diabetes treatments"))
    temp <- temp %>% filter(
      ncurrtx2 < 5
    ) %>% select(-ncurrtx2)
    q = "ncurrtx_upto_4"
  }
  
  if (n == 4) {
    print(paste0("All previous exclusion rules combined"))
    temp <- temp %>% filter(
      ncurrtx2 < 5 & INS == F
    ) %>% select(-ncurrtx2) %>%
      group_by(.imp, patid) %>%
      filter(
        !(studydrug2 == "SGLT2 + GLP1" & "SGLT2 + DPP4/SU" %in% studydrug2)
      ) %>%
      ungroup()
    q = "all_3_exclusions"
  }
  
  gc()
  
  ps.formula <- formula(paste0("studydrug", main, " ~ ", paste(covariates, collapse=" + ")))
  
  # force temp to be data.frame() for SumStat function
  temp <- temp %>% as.data.frame(temp)
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp[temp$.imp == i,],
                         weight = c("IPW", "overlap"))
    
    # truncate IPTW at 5th and 95th percentile
    w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.05, 0.95))[1],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.05, 0.95))[1],
                                       ifelse(
                                         w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.05, 0.95))[2],
                                         quantile(w.overlap$ps.weights$IPW, probs = c(0.05, 0.95))[2],
                                         w.overlap$ps.weights$IPW
                                       ))
    
    # truncate overlap weights at 1st and 99th percentile
    w.overlap$ps.weights$overlap <- ifelse(w.overlap$ps.weights$overlap < quantile(w.overlap$ps.weights$overlap, probs = c(0.01, 0.99))[1],
                                           quantile(w.overlap$ps.weights$overlap, probs = c(0.01, 0.99))[1],
                                           ifelse(
                                             w.overlap$ps.weights$overlap > quantile(w.overlap$ps.weights$overlap, probs = c(0.01, 0.99))[2],
                                             quantile(w.overlap$ps.weights$overlap, probs = c(0.01, 0.99))[2],
                                             w.overlap$ps.weights$overlap
                                           ))
    
    
    ## Add weights to data frame
    weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
    temp[temp$.imp == i,][[paste0("IPTW", main, collapse = "")]]  <- weights$IPW
    temp[temp$.imp == i,][[paste0("overlap", main, collapse = "")]] <- weights$overlap
    
    rm(w.overlap)
    rm(weights)
    gc()
  }
  
  temp_weights <- temp
  
  #add to weights to full dataset
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data.Rda"))
  
  temp <- temp %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  if (n == 1) {
    print(paste0("Remove multiple episodes per subject"))
    temp <- temp %>%
      group_by(.imp, patid) %>%
      filter(
        # Keep "SGLT2 + DPP4/SU" if patient also has "SGLT2 + GLP1"
        !(studydrug2 == "SGLT2 + GLP1" & "SGLT2 + DPP4/SU" %in% studydrug2)
      ) %>%
      ungroup()
    
  }
  
  if (n == 2) {
    print(paste0("Remove subjects on insulin"))
    temp <- temp %>% filter(INS == F)
  }
  
  if (n == 3) {
    print(paste0("Remove subjects on more than 4 current diabetes treatments"))
    temp <- temp %>% filter(
      ncurrtx2 < 5
    ) %>% select(-ncurrtx2)
  }
  
  if (n == 4) {
    print(paste0("All previous exclusion rules combined"))
    temp <- temp %>% filter(
      ncurrtx2 < 5 & INS == F
    ) %>% select(-ncurrtx2) %>%
      group_by(.imp, patid) %>%
      filter(
        !(studydrug2 == "SGLT2 + GLP1" & "SGLT2 + DPP4/SU" %in% studydrug2)
      ) %>%
      ungroup()
  }
  
  temp_all <- temp %>% filter(.imp != 0) 
  
  temp_weights <- temp_weights %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    ungroup()
  
  # ensure temp_all is in same row order as dataset with studydrug variable it will be merged with
  temp2 <- temp_all %>%
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  cohort <- temp2 %>% cbind(temp_weights %>% select(contains("overlap"), contains("IPTW")))
  
  rm(temp)
  rm(temp2)
  rm(temp_all)
  rm(temp_weights)
  
  
  # create empty data frame to which we can append the hazard ratios once calculated
  hrs <- data.frame()
  
  
  # main dataset is large - for speed of computation we will only load in dataset we need each time
  gc()
  
  
  # calculate hazard ratios
  
  # for every studydrug variable:

    print(paste0("Loading data for variable studydrug", main, collapse = ""))
    
    # define studydrug variable and weights variables to be used
    studydrug_var = paste0("studydrug", main)
    weights_overlap = paste0("overlap", main)
    weights_iptw = paste0("IPTW", main)

    gc()
    
    # define which drugs are evaluated with current studydrug variable
    drug_levels <- levels(cohort[[studydrug_var]])
    
    
    # remove double overlapping entries (take one only)
    cohort <- cohort %>% 
      group_by(.imp, patid, !!sym(studydrug_var)) %>% 
      arrange(dstartdate) %>% 
      filter(!duplicated(!!sym(studydrug_var))) %>% 
      ungroup()
    
    
    # calculate hazard ratios per outcome
    k = "ckd_egfr40"      
      
      # for other studydrug variables use regular censoring variables
      censvar_var=paste0(k, "_censvar")
      censtime_var=paste0(k, "_censtime_yrs")
      
      
      
      print(paste0("Calculating event numbers per drug level for outcome ", k))
      
      # calculate number of subjects in each group
      count <- cohort[cohort$.imp != 0,] %>%
        group_by(!!sym(studydrug_var)) %>%
        summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
        pivot_wider(names_from=!!sym(studydrug_var),
                    names_glue=paste0("{studydrug", main, "}_count"),
                    values_from=count)
      
      # calculate median follow up time (years) per group
      followup <- cohort[cohort$.imp != 0,] %>%
        group_by(!!sym(studydrug_var)) %>%
        summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
        pivot_wider(names_from=!!sym(studydrug_var),
                    names_glue=paste0("{studydrug", main, "}_followup"),
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
                    names_glue=paste0("{studydrug", main, "}_events"),
                    values_from=events)
      
      
      # write formulas for adjusted and unadjusted analyses
      f2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, collapse = ""))
      
      f_adjusted2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, " + ", paste(covariates, collapse=" + "), collapse = ""))
      
      # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
      
      for (n in drug_levels[-1]) {
        for (c in c("COEFS", "SE")) {
          for (d in c("unadj", "adj", "ow", "iptw")) {
            assign(paste0(c, ".", n, ".", d), rep(NA, n.imp))
            
          }
        }
      }
      
      
      for (i in 1:n.imp) {
        print(paste0("Analyses in imputed dataset number ", i))
        
        #unadjusted analyses first
        fit.unadj <- coxph(f2, cohort[cohort$.imp == i,])
        #adjusted analyses
        fit.adj <- coxph(f_adjusted2, cohort[cohort$.imp == i,])      
        #overlap weighted analyses
        fit.ow <- coxph(f_adjusted2, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
        #inverse probability of treatment weighted analyses
        fit.iptw <- coxph(f_adjusted2, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_iptw]])      
        
        #store coefficients and standard errors from this model
        for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coefficients to extract
          
          drug_name <- drug_levels[n+1]
          
          # for every analysis approach, extract coefficients
          for (d in c("unadj", "adj", "ow", "iptw")) {
            
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
        
        
        for (d in c("unadj", "adj", "ow", "iptw")) {
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
                                   variable = paste0("studydrug", main, collapse = ""),
                                   analysis = d,
                                   HR = pooled_hr[1],
                                   LB = pooled_hr[2],
                                   UB = pooled_hr[3],
                                   string = pooled_hr_string)
          
          # combine results by each analysis approach within each drug type (within studydrug variable)
          hrs <- rbind(hrs, outcome_hr)
          
          
        }
        
      }
      
  
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
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(hrs, file=paste0(today, "_hrs_sens_", q, ".Rda"))
  
  rm(hrs)
  rm(cohort)
}
