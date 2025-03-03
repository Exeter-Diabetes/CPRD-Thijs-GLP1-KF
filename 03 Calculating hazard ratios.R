########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights.Rda"))

############################1 HAZARD RATIOS################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs <- data.frame()

# number of studydrug variables
n.studydrug.vars <- sum(grepl("studydrug", colnames(cohort)))

# main dataset is large - for speed of computation we will only load in dataset we need each time
rm(cohort)
gc()


# calculate hazard ratios

# for every studydrug variable:
for (m in 1:n.studydrug.vars) {
  
  print(paste0("Loading data for variable studydrug", m, collapse = ""))
  
  # define studydrug variable and weights variables to be used
  studydrug_var = paste0("studydrug", m)
  weights_overlap = paste0("overlap", m)
  weights_iptw = paste0("IPTW", m)
  
  # reload data
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_imputed_data_withweights.Rda"))
  gc()
  
  cohort <- cohort %>% mutate(across(starts_with("studydrug"), as.factor))
  
  # define which drugs are evaluated with current studydrug variable
  drug_levels <- levels(cohort[[studydrug_var]])
  
  
  # remove double overlapping entries (take one only)
  cohort <- cohort %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  
  # calculate hazard ratios per outcome
  for (k in outcomes_per_drugclass) {
    
    print(paste0("Calculating event numbers per drug level for outcome ", k))
    
    if (studydrug_var == "studydrug1") {
      
      # studydrug1 is meant to check validity of taking DPP4/SU as one group as well as the difference between different GLP1 types
      # therefore we will use censoring variables _pp_ which will censor observations if switching between those
      censvar_var=paste0(k, "_pp_censvar")
      censtime_var=paste0(k, "_pp_censtime_yrs")
      
    } else {
      
      # for other studydrug variables use regular censoring variables
      censvar_var=paste0(k, "_censvar")
      censtime_var=paste0(k, "_censtime_yrs")
      
    }
    
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
    
    for (b in drug_levels[-1]) {
      for (c in c("COEFS", "SE")) {
        for (d in c("unadj", "adj", "ow", "iptw")) {
          assign(paste0(c, ".", b, ".", d), rep(NA, n.imp))
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
      for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coeffficients to extract
        
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
          # if drug level is SU (reference category) then HR will be NA
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

# ensure hazard ratio and 95% ci are stored as numeric variables
class(hrs$HR) <- class(hrs$LB) <- class(hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
hrs <- hrs %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
     mutate(
       `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
       `nN` = paste0(`events_number`, "/", `count`),
     )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs, file=paste0(today, "_hrs.Rda"))
