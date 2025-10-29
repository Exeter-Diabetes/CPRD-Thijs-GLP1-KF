########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

# set default studydrug variable
studydrug_var = paste0("studydrug", main)
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

# number of studydrug variables
n.studydrug.vars <- sum(grepl("studydrug", colnames(cohort)))
rm(cohort)
############################1 HAZARD RATIOS OVERALL################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs <- data.frame()


# main dataset is large - for speed of computation we will only load in dataset we need each time
gc()


# calculate hazard ratios

# for every studydrug variable:
for (m in 1:n.studydrug.vars) {
    
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
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs, file=paste0(today, "_hrs.Rda"))
############################2 HAZARD RATIOS FOR 5 YEAR FOLLOW-UP################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs_5y <- data.frame()


# main dataset is large - for speed of computation we will only load in dataset we need each time
gc()


# calculate hazard ratios

m = main
  
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
  for (k in outcomes[1]) {
    
    cohort <- cohort_all
  
    
    if (k == "retinopathy") {
      # for outcome retinopathy, only retain those without a history of this
      cohort <- cohort %>% filter(predrug_retinopathy == F)
      weights_overlap = paste0(weights_overlap, "_", k)
      weights_iptw = paste0(weights_iptw, "_", k) 
      
    }
    
    if (k == "acutepancreatitis") {
      # for outcome pancreatitis, only retain those without a history of this
      cohort <- cohort %>% filter(predrug_pancreatitis == F)
      weights_overlap = paste0(weights_overlap, "_", k)
      weights_iptw = paste0(weights_iptw, "_",k) 
      
    }

      analysis_approaches <-  c("ow")

    
    # for other studydrug variables use regular censoring variables
    censvar_var=paste0(k, "_5y_censvar")
    censtime_var=paste0(k, "_5y_censtime_yrs")
    
    
    
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
        for (d in c("ow")) {
          
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
      
      
      for (d in c("ow")) {
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
        hrs_5y <- rbind(hrs_5y, outcome_hr)
        
        
      }
      
    }
    
  }    
  
  

# ensure hazard ratio and 95% ci are stored as numeric variables
class(hrs_5y$HR) <- class(hrs_5y$LB) <- class(hrs_5y$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
hrs_5y <- hrs_5y %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs_5y, file=paste0(today, "_hrs_5y.Rda"))
############################3 HAZARD RATIOS FOR PER PROTOCOL ANALYSES################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs_pp <- data.frame()


# main dataset is large - for speed of computation we will only load in dataset we need each time
gc()


# calculate hazard ratios

m = main

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
for (k in outcomes[1]) {
  
  cohort <- cohort_all
  
  if (k == "retinopathy") {
    # for outcome retinopathy, only retain those without a history of this
    cohort <- cohort %>% filter(predrug_retinopathy == F)
    weights_overlap = paste0(weights_overlap, "_", k)
    weights_iptw = paste0(weights_iptw, "_", k) 
    
  }
  
  if (k == "acutepancreatitis") {
    # for outcome pancreatitis, only retain those without a history of this
    cohort <- cohort %>% filter(predrug_pancreatitis == F)
    weights_overlap = paste0(weights_overlap, "_", k)
    weights_iptw = paste0(weights_iptw, "_",k) 
    
  }
  
  analysis_approaches <-  c("ow")
  
  
  # for other studydrug variables use regular censoring variables
  censvar_var=paste0(k, "_pp_censvar")
  censtime_var=paste0(k, "_pp_censtime_yrs")
  
  
  
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
      for (d in c("ow")) {
        
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
    
    
    for (d in c("ow")) {
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
      hrs_pp <- rbind(hrs_pp, outcome_hr)
      
      
    }
    
  }
  
}    



# ensure hazard ratio and 95% ci are stored as numeric variables
class(hrs_pp$HR) <- class(hrs_pp$LB) <- class(hrs_pp$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
hrs_pp <- hrs_pp %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs_pp, file=paste0(today, "_hrs_pp.Rda"))
############################4 COMPETING RISK REGRESSION (FINE-GRAY)################################################################

# load data
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))
gc()

hrs_fg_ow <- data.frame()

studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)
weights_iptw = paste0("IPTW", main)

# analyse main outcome
k = "ckd_egfr40"
censvar_var=paste0(k, "_censvar")
censtime_var=paste0(k, "_censtime_yrs")

# create new censoring variables
cohort <- cohort %>% mutate(
  status = ifelse(!!sym(censvar_var) == 1, 1, ifelse(death_censvar == 1 & death_censtime_yrs <= !!sym(censtime_var), 2, 0)),
  status = factor(status, levels = c("0", "1", "2")),
  censtime_yrs = ifelse(status == 2, death_censtime_yrs, !!sym(censtime_var))
)

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


drug_levels <- levels(cohort[[studydrug_var]])

fg_data <- finegray(Surv(censtime_yrs, status) ~ ., data = cohort)

fg_data <- fg_data %>% mutate(fg_ow = fgwt * !!sym(weights_overlap))

f_fg_adjusted <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ studydrug", main, " + ", paste0(covariates, collapse = " + ")))

for (n in drug_levels[-1]) {
  for (c in c("COEFS", "SE")) {
    for (d in c("fg_ow")) {
      assign(paste0(c, ".", n, ".", d), rep(NA, n.imp))
      
    }
  }
}


for (i in 1:n.imp) {
  print(paste0("Analyses in imputed dataset number ", i))
  
  fit.fg_ow <- coxph(f_fg_adjusted, data = fg_data %>% filter(.imp == i), weights = fg_ow)
  
  
  #store coefficients and standard errors from this model
  for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coefficients to extract
    
    drug_name <- drug_levels[n+1]
    
    # for every analysis approach, extract coefficients
    for (d in c("fg_ow")) {
      
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
  
  
  for (d in c("fg_ow")) {
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
    
    
    
  }
  
  hrs_fg_ow <- rbind(hrs_fg_ow, outcome_hr)
  
}

# ensure hazard ratio and 95% ci are stored as numeric variables
class(hrs_fg_ow$HR) <- class(hrs_fg_ow$LB) <- class(hrs_fg_ow$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
hrs_fg_ow <- hrs_fg_ow %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

# save result
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs_fg_ow, file=paste0(today, "_hrs_fg_ow.Rda"))


############################5 HAZARD RATIOS BY PRESENCE OF COMORBIDITIES################################################################

# consider following factors that may potentially alter effect of GLP1-RA on kidney protection:
factors <- c("malesex", "white_ethnicity", "low_egfr", "albuminuria", "predrug_cvd", "predrug_heartfailure")

# calculate hazard ratios

# define studydrug variable and weights variables to be used
studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)
weights_iptw = paste0("IPTW", main)

# reload data
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))
gc()

cohort <- cohort %>% mutate(across(starts_with("studydrug"), as.factor),
                            white_ethnicity = ifelse(ethnicity_4cat == "White", T, F),
                            low_egfr = ifelse(preegfr <60, T, F),
                            albuminuria = as.logical(albuminuria),
                            malesex = as.logical(malesex),
                            predrug_cvd = as.logical(predrug_cvd), 
                            predrug_heartfailure = as.logical(predrug_heartfailure),
)



# define which drugs are evaluated with current studydrug variable
drug_levels <- levels(cohort[[studydrug_var]])


# remove double overlapping entries (take one only)
cohort <- cohort %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>%
  arrange(dstartdate) %>%
  filter(!duplicated(!!sym(studydrug_var))) %>%
  ungroup()



# create empty data frame to which we can append the hazard ratios once calculated
factor_hrs <- data.frame()

for (k in outcomes[1]) {
  
  print(paste0("Calculating event numbers per drug level for outcome ", k))
  
  # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
  for (n in drug_levels) {
    for (c in c("COEFS", "SE")) {
      for (l in factors) {
        for (q in levels(as.factor(cohort[[l]]))) {
          assign(paste0(c, ".", n, ".", l, q), rep(NA, n.imp))
        }
      }
    }
  }
  
  # calculate hazard ratios per outcome
  for (l in factors) {
    
    d = "ow"
    print(paste0("Analyses by ", l))
    
    
    # for other studydrug variables use regular censoring variables
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    
    # calculate number of subjects in each group
    count <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_count"),
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_followup"),
                  values_from=time)
    
    # summarise number of events per group
    events <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(!!sym(studydrug_var), !!sym(l), events) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_events"),
                  values_from=events)
    
    # to avoid collinearity, remove variables ethnicity_4cat when relevant
    covariates_all = covariates
    
    if (grepl("ethnicity", l)) {
      covariates <- covariates %>% setdiff(covariates[grepl("ethnicity", covariates)])
    }
    
    # write formula for analyses with interaction
    f_adjusted_interaction <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, "*", l, " + ", paste(covariates, collapse=" + "), collapse = ""))
    
    # empty vector for wald statistic
    wald = rep(NA, n.imp)
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      # analyses with interaction
      model <- coxph(f_adjusted_interaction, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
      
      # save wald statistic for interaction term
      wald[i] = car::Anova(model, type = 3, test = "Wald")[paste0("studydrug", main, ":", l),"Chisq"]
      
      
      for (n in 2:length(drug_levels)) {   # skip reference drug (n=1)
        drug_name <- drug_levels[n]
        coef_name <- paste0("studydrug", main, drug_name)
        
        for (q in levels(as.factor(cohort[[l]]))) {
          ref_level <- levels(as.factor(cohort[[l]]))[1]  # reference of the factor
          coef_drug <- model$coefficients[coef_name]      # main effect for drug vs ref
          var_drug  <- vcov(model)[coef_name, coef_name]
          
          if (q == ref_level) {
            # HR for drug at reference level of factor
            beta <- coef_drug
            se   <- sqrt(var_drug)
            
          } else {
            # interaction term name, e.g. "drugX:factorY"
            inter_term <- paste0(coef_name, ":", l, q)
            if (!(inter_term %in% names(model$coefficients))) {
              # fallback if R encodes interaction in other order (factor:drug)
              inter_term <- paste0(l, q, ":", coef_name)
            }
            
            coef_int <- model$coefficients[inter_term]
            var_int  <- vcov(model)[inter_term, inter_term]
            covar    <- vcov(model)[coef_name, inter_term]
            
            beta <- coef_drug + coef_int
            se   <- sqrt(var_drug + var_int + 2*covar)
          }
          
          # store results
          assign(paste0("COEFS.", drug_name, ".", l, q), 
                 replace(get(paste0("COEFS.", drug_name, ".", l, q)), i, beta))
          
          assign(paste0("SE.", drug_name, ".", l, q), 
                 replace(get(paste0("SE.", drug_name, ".", l, q)), i, se))
        }
      }
      
      
    }
    
    
    ## loop to pool and store results
    for (n in 1:length(drug_levels)) {
      
      # get drug name
      drug_name <- drug_levels[n]
      
      
       for (q in levels(as.factor(cohort[[l]])))  {
        
        if (n == 1 & q == levels(as.factor(cohort[[l]]))[1]) {
          # if drug level is reference category then HR will be NA
          pooled_hr <- c(1, NA, NA)
          pooled_hr_string <- "1.00 (ref.)"
          
        } else {
          
          # define names for objects containing pooled hr + 95% ci (as vector and as string)
          pooled_hr_name <- paste0(drug_name, "_", l, q, "_hr")
          pooled_hr_string_name <- paste0(drug_name, "_", l, q, "_string")
          
          # define names for coefficient vectors
          coef_name <- paste0("COEFS.", drug_name, ".", l, q)
          se_name <- paste0("SE.", drug_name, ".", l, q)
          
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
                                 factor = l,
                                 condition = q,
                                 count = as.numeric(count %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_count"))),
                                 followup = as.numeric(followup %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_followup"))),
                                 events = as.character(events %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_events"))),
                                 contrast = paste0(drug_name, " vs ", drug_levels[1], collapse = ""),
                                 variable = paste0("studydrug", main, collapse = ""),
                                 analysis = d,
                                 HR = pooled_hr[1],
                                 LB = pooled_hr[2],
                                 UB = pooled_hr[3],
                                 string = pooled_hr_string)
        
        # combine results by each analysis approach within each drug type (within studydrug variable)
        factor_hrs <- rbind(factor_hrs, outcome_hr)
      }
      
    }
    
    
    # Calculate the p-value from the Wald statistic (using normal approximation)
    
    m <- length(wald)
    Qbar <- mean(wald)
    B <- var(wald)
    Tstat <- Qbar + (1 + 1/m) * B
    df_interaction <- nlevels(as.factor(cohort[[l]])) - 1
    Fstat <- Tstat / df_interaction
    r <- (1 + 1/m) * B / Qbar
    nu <- (m - 1) * (1 + 1 / r)^2
    p_value_interaction <- 1 - pf(Fstat, df1 = df_interaction, df2 = nu)
    
    print(paste0("P value for ", l, " interaction: ", p_value_interaction))
    assign(paste0("p_value_interaction_", l), p_value_interaction)
    
    covariates <- covariates_all
    rm(covariates_all)
    
  }
  
}

# ensure hazard ratio and 95% ci are stored as numeric variables
class(factor_hrs$HR) <- class(factor_hrs$LB) <- class(factor_hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
factor_hrs <- factor_hrs %>%
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

factor_hrs$p_value_interaction <- NA  

for (l in unique(factor_hrs$factor)) {
  pval <- get(paste0("p_value_interaction_", l))
  factor_hrs$p_value_interaction[factor_hrs$factor == l] <- pval
}


# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(factor_hrs, file=paste0(today, "_factor_hrs.Rda"))


############################6 HAZARD RATIOS BY WEIGHT LOSS TERTILES################################################################

# consider following tertile that may potentially alter effect of GLP1-RA on kidney protection:
tertile <- c("weight_loss_tertiles")

# calculate hazard ratios

# define studydrug variable and weights variables to be used
studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)
weights_iptw = paste0("IPTW", main)

# reload data
print(paste0("Processing data for variable studydrug", main, collapse = ""))

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", main, ".Rda"))
gc()


# keep those with weight loss data at 6 months only
temp <- temp %>% filter(!is.na(weight_pct_change))

# re-estimate propensity score weights

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
  w.overlap$ps.weights$IPW <- ifelse(w.overlap$ps.weights$IPW < quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                     quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[1],
                                     ifelse(
                                       w.overlap$ps.weights$IPW > quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       quantile(w.overlap$ps.weights$IPW, probs = c(0.02, 0.98))[2],
                                       w.overlap$ps.weights$IPW
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


# keep those with weight loss data at 6 months only
temp <- temp %>% filter(!is.na(weight_pct_change))

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

# get tertiles of weight loss

tertile_values <- quantile(cohort$weight_pct_change, c(0:3/3))

cohort <- cohort %>% mutate(across(starts_with("studydrug"), as.factor),
)

cohort$weight_loss_tertiles <- with(cohort, cut(weight_pct_change, 
                                                tertile_values, 
                                                include.lowest = T, 
                                                labels = c("Low", "Medium", "High")))

cohort$weight_loss_tertiles = as.factor(cohort$weight_loss_tertiles)


# define which drugs are evaluated with current studydrug variable
drug_levels <- levels(cohort[[studydrug_var]])


# remove double overlapping entries (take one only)
cohort <- cohort %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>%
  arrange(dstartdate) %>%
  filter(!duplicated(!!sym(studydrug_var))) %>%
  ungroup()



# create empty data frame to which we can append the hazard ratios once calculated
weight_loss_tertile_hrs <- data.frame()

for (k in outcomes[1]) {
  
  print(paste0("Calculating event numbers per drug level for outcome ", k))
  
  # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
  for (n in drug_levels) {
    for (c in c("COEFS", "SE")) {
      for (l in tertile) {
        for (q in levels(as.factor(cohort[[l]]))) {
          assign(paste0(c, ".", n, ".", l, q), rep(NA, n.imp))
        }
      }
    }
  }
  
  # calculate hazard ratios per outcome
  for (l in tertile) {
    
    d = "ow"
    print(paste0("Analyses by ", l))
    
    
    # for other studydrug variables use regular censoring variables
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    
    # calculate number of subjects in each group
    count <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_count"),
                  values_from=count)
    
    # calculate median follow up time (years) per group
    followup <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_followup"),
                  values_from=time)
    
    # summarise number of events per group
    events <- cohort[cohort$.imp != 0,] %>%
      group_by(!!sym(studydrug_var), !!sym(l)) %>%
      summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
                drug_count=round(n()/n.imp, 0)) %>%
      mutate(events_perc=round(event_count*100/drug_count, 1),
             events=paste0(event_count, " (", events_perc, "%)")) %>%
      select(!!sym(studydrug_var), !!sym(l), events) %>%
      pivot_wider(names_from=!!sym(studydrug_var),
                  names_glue=paste0("{studydrug", main, "}_events"),
                  values_from=events)
    
    # to avoid collinearity, remove variables ethnicity_4cat and dstartdate_age when relevant
    covariates_all = covariates
    
    if (grepl("ethnicity", l)) {
      covariates <- covariates %>% setdiff(covariates[grepl("ethnicity", covariates)])
    }
    
    if (grepl("age", l)) {
      covariates <- covariates %>% setdiff(covariates[grepl("age", covariates)])
    }
    
    # write formula for analyses with interaction
    f_adjusted_interaction <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, "*", l, " + ", paste(covariates, collapse=" + "), collapse = ""))
    
    # empty vector for wald statistic
    wald = rep(NA, n.imp)
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      # analyses with interaction
      model <- coxph(f_adjusted_interaction, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
      
      # save wald statistic for interaction term
      wald[i] = car::Anova(model, type = 3, test = "Wald")[paste0("studydrug", main, ":", l),"Chisq"]
      
      
      for (n in 2:length(drug_levels)) {   # skip reference drug (n=1)
        drug_name <- drug_levels[n]
        coef_name <- paste0("studydrug", main, drug_name)
        
        for (q in levels(as.factor(cohort[[l]]))) {
          ref_level <- levels(as.factor(cohort[[l]]))[1]  # reference of the factor
          coef_drug <- model$coefficients[coef_name]      # main effect for drug vs ref
          var_drug  <- vcov(model)[coef_name, coef_name]
          
          if (q == ref_level) {
            # HR for drug at reference level of factor
            beta <- coef_drug
            se   <- sqrt(var_drug)
            
          } else {
            # interaction term name, e.g. "drugX:factorY"
            inter_term <- paste0(coef_name, ":", l, q)
            if (!(inter_term %in% names(model$coefficients))) {
              # fallback if R encodes interaction in other order (factor:drug)
              inter_term <- paste0(l, q, ":", coef_name)
            }
            
            coef_int <- model$coefficients[inter_term]
            var_int  <- vcov(model)[inter_term, inter_term]
            covar    <- vcov(model)[coef_name, inter_term]
            
            beta <- coef_drug + coef_int
            se   <- sqrt(var_drug + var_int + 2*covar)
          }
          
          # store results
          assign(paste0("COEFS.", drug_name, ".", l, q), 
                 replace(get(paste0("COEFS.", drug_name, ".", l, q)), i, beta))
          
          assign(paste0("SE.", drug_name, ".", l, q), 
                 replace(get(paste0("SE.", drug_name, ".", l, q)), i, se))
        }
      }
      
      
    }
    
    
    ## loop to pool and store results
    for (n in 1:length(drug_levels)) {
      
      # get drug name
      drug_name <- drug_levels[n]
      
      
      for (q in levels(as.factor(cohort[[l]])))  {
        
        if (n == 1 & q == levels(as.factor(cohort[[l]]))[1]) {
          # if drug level is reference category then HR will be NA
          pooled_hr <- c(1, NA, NA)
          pooled_hr_string <- "1.00 (ref.)"
          
        } else {
          
          # define names for objects containing pooled hr + 95% ci (as vector and as string)
          pooled_hr_name <- paste0(drug_name, "_", l, q, "_hr")
          pooled_hr_string_name <- paste0(drug_name, "_", l, q, "_string")
          
          # define names for coefficient vectors
          coef_name <- paste0("COEFS.", drug_name, ".", l, q)
          se_name <- paste0("SE.", drug_name, ".", l, q)
          
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
                                 tertile = l,
                                 condition = q,
                                 count = as.numeric(count %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_count"))),
                                 followup = as.numeric(followup %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_followup"))),
                                 events = as.character(events %>% filter(!!sym(l)==q) %>% select(paste0(drug_name, "_events"))),
                                 contrast = paste0(drug_name, " vs ", drug_levels[1], collapse = ""),
                                 variable = paste0("studydrug", main, collapse = ""),
                                 analysis = d,
                                 HR = pooled_hr[1],
                                 LB = pooled_hr[2],
                                 UB = pooled_hr[3],
                                 string = pooled_hr_string)
        
        # combine results by each analysis approach within each drug type (within studydrug variable)
        weight_loss_tertile_hrs <- rbind(weight_loss_tertile_hrs, outcome_hr)
      }
      
    }
    
    
    # Calculate the p-value from the Wald statistic (using normal approximation)

    m <- length(wald)
    Qbar <- mean(wald)
    B <- var(wald)
    Tstat <- Qbar + (1 + 1/m) * B
    df_interaction <- nlevels(as.factor(cohort[[l]])) - 1
    Fstat <- Tstat / df_interaction
    r <- (1 + 1/m) * B / Qbar
    nu <- (m - 1) * (1 + 1 / r)^2
    p_value_interaction <- 1 - pf(Fstat, df1 = df_interaction, df2 = nu)
    
    print(paste0("P value for ", l, " interaction: ", p_value_interaction))
    assign(paste0("p_value_interaction_", l), p_value_interaction)
    
    covariates <- covariates_all
    rm(covariates_all)
    
  }
  
}

# ensure hazard ratio and 95% ci are stored as numeric variables
class(weight_loss_tertile_hrs$HR) <- class(weight_loss_tertile_hrs$LB) <- class(weight_loss_tertile_hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
weight_loss_tertile_hrs <- weight_loss_tertile_hrs %>%
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
  )

weight_loss_tertile_hrs$p_value_interaction <- NA  

for (l in unique(weight_loss_tertile_hrs$tertile)) {
  pval <- get(paste0("p_value_interaction_", l))
  weight_loss_tertile_hrs$p_value_interaction[weight_loss_tertile_hrs$tertile == l] <- pval
}


# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(weight_loss_tertile_hrs, file=paste0(today, "_weight_loss_tertile_hrs.Rda"))
