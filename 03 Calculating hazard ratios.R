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

main = 2
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
  
  
  # calculate hazard ratios per outcome
  for (k in outcomes) {
    

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
       `nN` = paste0("  ", format(as.numeric(`events_number`), big.mark = ",", scientific = F), " / ", format(`count`, big.mark = ",", scientific = F)),
     )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(hrs, file=paste0(today, "_hrs.Rda"))


############################2 HAZARD RATIOS BY PRESENCE OF COMORBIDITIES################################################################

# consider following factors that may potentially alter effect of GLP1 on kidney protection:
factors <- c("malesex", "white_ethnicity", "predrug_cvd", "predrug_heartfailure")

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
                            age_above_60 = ifelse(dstartdate_age >=60, T, F),
                            white_ethnicity = ifelse(ethnicity_5cat == "White", T, F),
                            malesex = as.logical(malesex),
                            predrug_cvd = as.logical(predrug_cvd), 
                            predrug_heartfailure = as.logical(predrug_heartfailure),
                            albuminuria = as.logical(albuminuria))



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

for (k in outcomes) {
  
  print(paste0("Calculating event numbers per drug level for outcome ", k))
  
  # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
  for (n in drug_levels) {
    for (c in c("COEFS", "SE")) {
      for (l in factors) {
        for (q in c("TRUE", "FALSE")) {
          assign(paste0(c, ".", n, ".", l, q), rep(NA, n.imp))
        }
      }
    }
  }
  
  # for interaction term p value
  for (c in c("COEFS", "SE")) {
    for (l in factors) {
      assign(paste0(c, ".interaction.term.only.", l), rep(NA, n.imp))
      
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
    
    
    # write formulas for analyses with and without interaction
    f_adjusted <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, " + ", paste(covariates, collapse=" + "), collapse = ""))
    
    f_adjusted_interaction <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", main, "*", l, " + ", paste(covariates, collapse=" + "), collapse = ""))
    
    for (i in 1:n.imp) {
      print(paste0("Analyses in imputed dataset number ", i))
      
      # analyses without interaction
      fit.no_interaction <- coxph(f_adjusted, cohort[cohort$.imp == i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
      # analyses with interaction
      fit.interaction <- coxph(f_adjusted_interaction, cohort[cohort$.imp ==i,], weights = cohort[cohort$.imp ==i,][[weights_overlap]])
      
      
      #store coefficients and standard errors from this model
      for (n in 1:length(drug_levels[-1])) { # 1st is reference category so no coefficients to extract
        
        drug_name <- drug_levels[n+1]
        
        for (q in c("FALSE", "TRUE")) {
          
          p = ifelse(q == "TRUE", 2, 1)
          
          model <- fit.interaction
          
          if (q == "FALSE") {
            # coef_vector[i] <- model$coefficients[n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", l, q, "`[", i, "] <- model$coefficients[", n, "]", collapse = "")
            # se_vector[i] <- sqrt(model$var[n,n])
            se_statement <- paste0("`SE.", drug_name, ".", l, q, "`[", i, "] <- sqrt(model$var[", n, ",", n, "])", collapse = "")
          } else {
            # coef = coef1 + coef2
            # coef_vector[i] <- model$coefficients[n] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(p-1)+n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", l, q, "`[", i, "] <- model$coefficients[", n, "] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(", p, "-1) + ", n, "]", collapse = "")
            
            # se = sqrt(var1 + var2 + cov1,2)
            # se_vector[i] <- sqrt(abs(model$var[n,n]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(p-1)+n,nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(p-1)+n]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(p-1)+n]))
            se_statement <- paste0("`SE.", drug_name, ".", l, q, "`[", i, "] <- sqrt(abs(model$var[", n, ",", n, "]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(", p, "-1) + ", n, ",nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(", p, "-1) + ", n, "]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(as.factor(cohort[[l]]))[-1]))/(", p, "-1) + ", n, "]))", collapse = "")
          }
          
          # for interaction term itself only
          coef_statement2 <- paste0("`COEFS.interaction.term.only.", l, "`[", i, "] <- model$coefficients[length(model$coefficients)]", collapse = "")
          se_statement2 <- paste0("`SE.interaction.term.only.", l, "`[", i, "] <- sqrt(abs(model$var[length(model$coefficients), length(model$coefficients)]))", collapse = "")
          
          # execute commands
          eval(str2lang(coef_statement))
          eval(str2lang(se_statement))
          
          eval(str2lang(coef_statement2))
          eval(str2lang(se_statement2))
          
          
          
          
        }
        
      }
      
      
    }
    
    
    ## loop to pool and store results
    for (n in 1:length(drug_levels)) {
      
      # get drug name
      drug_name <- drug_levels[n]
      
      
      for (q in c("FALSE", "TRUE")) {
        
        if (n == 1 & q == "FALSE") {
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
          
          
          # test significance interaction
          coef_name2 <- paste0("COEFS.interaction.term.only.", l)
          se_name2 <- paste0("SE.interaction.term.only.", l)
          
          coef_vector <- get(coef_name2)
          se_vector <- get(se_name2)
          # Calculate the Wald statistic
          wald_stat <- mean(coef_vector) / sqrt(mean(se_vector^2) + var(coef_vector) / n.imp)
          
          # Calculate the p-value from the Wald statistic (using normal approximation)
          p_value_interaction <- 2 * (1 - pnorm(abs(wald_stat)))
          print(paste0("P value for ", l, " interaction (", levels(factor(cohort[[studydrug_var]]))[length(levels(factor(cohort[[studydrug_var]])))], " with ", levels(factor(cohort[[studydrug_var]]))[1], "): ", p_value_interaction))
          
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

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(factor_hrs, file=paste0(today, "_factor_hrs.Rda"))
