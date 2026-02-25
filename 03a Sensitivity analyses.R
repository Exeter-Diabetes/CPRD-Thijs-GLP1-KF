########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

## sensitivity analyses - excluding several subgroups
## sensitivity analysis including those with missing uacr is in separate codes in folder /sens/

# primary outcome only
outcomes_sensitivity = outcomes[1]

# primary analysis approach only
analysis_approaches = "ow"
############################1 CALCULATE WEIGHTS AND HAZARD RATIOS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in themain, then stack them again at the end
gc()

# calculate weights
for (r in 1:5) {
  

  studydrug_var = paste0("studydrug", main)
  
  print(paste0("Processing data for variable studydrug", main, collapse = ""))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", main, ".Rda"))
  
    if (r == 1) {
    print(paste0("Remove multiple episodes per subject"))
      #subjects can  be in the comparator group, and later on have a 2nd episode in the GLP1 group
      temp <- temp %>%
        group_by(.imp, patid) %>%
        filter(
          # only keep "SGLT2i + DPP4i/SU" if patient also has "SGLT2i + GLP1-RA"
          !(studydrug2 == "SGLT2i + GLP1-RA" & "SGLT2i + DPP4i/SU" %in% studydrug2)
        ) %>%
        ungroup()
      q = "single_episodes"
    }
  
  if (r == 2) {
    print(paste0("Remove subjects on insulin"))
    temp <- temp %>% filter(INS == F)
    q = "no_insulin"
  }
  
  if (r == 3) {
    print(paste0("Remove subjects on more than 4 current diabetes treatments"))
    temp <- temp %>% filter(
      ncurrtx2 < 5
    ) %>% select(-ncurrtx2)
    q = "ncurrtx_upto_4"
  }
  
  if (r == 4) {
    print(paste0("All previous exclusion rules combined"))
    temp <- temp %>% filter(
      ncurrtx2 < 5 & INS == F
    ) %>% select(-ncurrtx2) %>%
      group_by(.imp, patid) %>%
      filter(
        !(studydrug2 == "SGLT2i + GLP1-RA" & "SGLT2i + DPP4i/SU" %in% studydrug2)
      ) %>%
      ungroup()
    q = "all_3_exclusions"
  }
  
  if (r == 5) {
    print(paste0("Remove subjects without linkage to secondary care data"))
    temp <- temp %>% filter(with_hes == 1)
    q = "hes_data_only"
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
                         weight = c("overlap"))
    

    ## Add weights to data frame
    weights <- w.overlap$ps.weights # note that these do not contain an index variable but are in the same order as our data frame
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
  

  if (r == 1) {
    print(paste0("Remove multiple episodes per subject"))
    temp <- temp %>%
      group_by(.imp, patid) %>%
      filter(
        # only keep "SGLT2i + DPP4i/SU" if patient also has "SGLT2i + GLP1-RA"
        !(studydrug2 == "SGLT2i + GLP1-RA" & "SGLT2i + DPP4i/SU" %in% studydrug2)
      ) %>%
      ungroup()
    
  }
  
  if (r == 2) {
    print(paste0("Remove subjects on insulin"))
    temp <- temp %>% filter(INS == F)
  }
  
  if (r == 3) {
    print(paste0("Remove subjects on more than 4 current diabetes treatments"))
    temp <- temp %>% filter(
      ncurrtx2 < 5
    ) %>% select(-ncurrtx2)
  }
  
  if (r == 4) {
    print(paste0("All previous exclusion rules combined"))
    temp <- temp %>% filter(
      ncurrtx2 < 5 & INS == F
    ) %>% select(-ncurrtx2) %>%
      group_by(.imp, patid) %>%
      filter(
        !(studydrug2 == "SGLT2i + GLP1-RA" & "SGLT2i + DPP4i/SU" %in% studydrug2)
      ) %>%
      ungroup()
  
  }
  
  if (r == 5) {
    print(paste0("Remove subjects without linkage to secondary care data"))
    temp <- temp %>% filter(with_hes == 1)
    
    outcomes_sensitivity <- outcomes_per_drugclass
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
    
    
    # for safety remove double overlapping entries (should be redundant code)
    cohort <- cohort %>% 
      group_by(.imp, patid, !!sym(studydrug_var)) %>% 
      arrange(dstartdate) %>% 
      filter(!duplicated(!!sym(studydrug_var))) %>% 
      ungroup()
    
    
    # calculate hazard ratios per outcome
    for (k in outcomes_sensitivity) {      

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
          for (d in c(analysis_approaches)) {
            assign(paste0(c, ".", n, ".", d), rep(NA, n.imp))
            
          }
        }
      }
      
      
      for (i in 1:n.imp) {
        print(paste0("Analyses in imputed dataset number ", i))
        
        #overlap weighted model only
        fit.ow <- coxph(f_adjusted2, cluster = patid, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
 
        #store coefficients and standard errors from this model
        for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coefficients to extract
          
          drug_name <- drug_levels[n+1]
          
          # for every analysis approach, extract coefficients
          for (d in c(analysis_approaches)) {
            
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
        
        
        for (d in c(analysis_approaches)) {
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


############################2 POST-HOC ANALYSIS: ABSOLUTE RISK DIFFERENCE################################################################


k  <- "ckd_egfr40"
censvar_var  <- paste0(k, "_censvar")
censtime_var <- paste0(k, "_censtime_yrs")

studydrug_var   <- paste0("studydrug", main)
weights_overlap <- paste0("overlap", main)

outcome_variables <- paste0("(", paste("ckd_egfr40", collapse = "|"), ")(?!.*(5y|sens|pp)).*(_censtime_yrs|_censvar)$")


cohort <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         # hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year))  %>%
  #select relevant variables only
  select(patid, .imp, !!sym(studydrug_var), !!sym(weights_overlap),
         matches(outcome_variables, perl = TRUE), 
         all_of(covariates)) 


print(paste0("Starting G-computation of absolute risk difference for outcome: ", k))

# Ensure studydrug variable is factor variable
cohort[[studydrug_var]] <- as.factor(cohort[[studydrug_var]])
drug_levels <- levels(cohort[[studydrug_var]])

# Formula for double-robust model
f_adjusted <- as.formula(
  paste0("Surv(", censtime_var, ", ", censvar_var, ") ~ ", 
         studydrug_var, " + ", paste(covariates, collapse = " + "))
)

# Create empty vectors to store risks and their SEs per imputation
# G-computation requires storing the risk and its variance from each imputed set
assign(paste0("arr_imp.", k), rep(NA, n.imp))
assign(paste0("arr_se.", k),  rep(NA, n.imp))

for (i in 1:n.imp) {
  print(paste0("ARD calculation for imputation ", i, " of ", n.imp))
  
  # Fit the doubly robust outcome model
  fit <- coxph(f_adjusted, 
               data = cohort %>% filter(.imp == i), 
               weights = cohort[cohort$.imp == i,][[weights_overlap]])
  
  # Counterfactual 1: Everyone Treated (GLP1-RA)
  sf1 <- survfit(fit, 
                 newdata = cohort %>% filter(.imp == i) %>% mutate(!!sym(studydrug_var) := drug_levels[2]), 
                 weights = cohort[cohort$.imp == i,][[weights_overlap]])
  
  s1_val <- weighted.mean(summary(sf1, times = 3)$surv, 
                          cohort[cohort$.imp == i,][[weights_overlap]], 
                          na.rm = TRUE)
  
  # Counterfactual 0: Everyone Untreated (DPP4i/SU)
  sf0 <- survfit(fit, 
                 newdata = cohort %>% filter(.imp == i) %>% mutate(!!sym(studydrug_var) := drug_levels[1]), 
                 weights = cohort[cohort$.imp == i,][[weights_overlap]])
  
  s0_val <- weighted.mean(summary(sf0, times = 3)$surv, 
                          cohort[cohort$.imp == i,][[weights_overlap]], 
                          na.rm = TRUE)
  
  # Calculate risk difference for this imputation
  # (1-s1) - (1-s0)
  arr_val <- (1 - s1_val) - (1 - s0_val)
  
  # Extract SE for the difference (standard error of the survival estimate)
  # Using the delta method logic or simple combined SE
  
  se_val <- sqrt(
    (sum(cohort[cohort$.imp == i,][[weights_overlap]]^2 * summary(sf0, times = 3)$std.err^2, na.rm = TRUE) / sum(cohort[cohort$.imp == i,][[weights_overlap]], na.rm = TRUE)^2) + 
      (sum(cohort[cohort$.imp == i,][[weights_overlap]]^2 * summary(sf1, times = 3)$std.err^2, na.rm = TRUE) / sum(cohort[cohort$.imp == i,][[weights_overlap]], na.rm = TRUE)^2)
  )  
  
  # assign names and store
  arr_vec <- get(paste0("arr_imp.", k))
  arr_vec[i] <- arr_val
  assign(paste0("arr_imp.", k), arr_vec)
  
  se_vec <- get(paste0("arr_se.", k))
  se_vec[i] <- se_val
  assign(paste0("arr_se.", k), se_vec)
  
  rm(sf1, sf0, fit, arr_vec, se_vec)
  if(i %% 10 == 0) gc() 
}

# save results
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(arr_imp.ckd_egfr40, file=paste0(today, "_arr.Rda"))
save(arr_se.ckd_egfr40, file=paste0(today, "_arr_se.Rda"))


# pool results using Rubin's rules
COEFS <- get(paste0("arr_imp.", k))
SE    <- get(paste0("arr_se.", k))

mean.coef <- mean(COEFS)
W <- mean(SE^2)
B <- var(COEFS)
T.var <- W + (1+1/n.imp)*B
se.coef <- sqrt(T.var)

# calculate confidence intervals
LB.CI <- mean.coef - (se.coef*1.96)
UB.CI <- mean.coef + (se.coef*1.96)

# summary of results
ard_results <- data.frame(
  outcome = k,
  ARD     = mean.coef,
  LB      = LB.CI,
  UB      = UB.CI,
  NNT     = 1 / abs(mean.coef)
)

print(ard_results)