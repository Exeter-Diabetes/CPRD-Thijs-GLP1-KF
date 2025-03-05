########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights.Rda"))

############################1 HAZARD RATIOS OVERALL################################################################

# create empty data frame to which we can append the hazard ratios once calculated
hrs <- data.frame()

# add extra studydrug variable with SGLT2 as reference
cohort$studydrug4 <- relevel(cohort$studydrug2, ref = "SGLT2")

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
    
    if (m == 1 | m == 3) {
      
      # studydrug1 is meant to check validity of taking DPP4/SU as one group and studydrug 3 the difference between different GLP1 types
      # therefore we will use censoring variables _sens1/_sens3 which will censor observations if switching between those
      censvar_var=paste0(k, "_sens", m, "_censvar")
      censtime_var=paste0(k, "_sens", m, "_censtime_yrs")
      
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


############################2 HAZARD RATIOS BY EGFR CATEGORY################################################################

# create empty data frame to which we can append the hazard ratios once calculated
egfr_hrs <- data.frame()

# main dataset is large - for speed of computation we will only load in dataset we need each time
rm(cohort)
gc()


# calculate hazard ratios

# use variable studydrug2
m = 2

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

# create variable for eGFR category
cohort <- cohort %>% mutate(
  egfr_cat = ifelse(preegfr < 45, "20-45", ifelse(preegfr < 60, "45-60", "≥60")),
  egfr_cat = factor(egfr_cat)
)

p_value_interaction_egfr <- rep(NA, n.imp)

# calculate hazard ratios per outcome
for (k in outcomes_per_drugclass) {
  
  print(paste0("Calculating event numbers per drug level for outcome ", k))
  
  
  
  # for other studydrug variables use regular censoring variables
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), egfr_cat) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_count"),
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), egfr_cat) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_followup"),
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), egfr_cat) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(!!sym(studydrug_var), egfr_cat, events) %>%
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_events"),
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, "*egfr_cat", collapse = ""))
  
  f_adjusted2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, "*egfr_cat", " + ", paste(covariates, collapse=" + "), collapse = ""))
  
  # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
  
  for (b in drug_levels[-1]) {
    for (c in c("COEFS", "SE")) {
      for (d in c("adj", "ow", "iptw")) {
        for (q in levels(as.factor(cohort$egfr_cat))) {
          assign(paste0(c, ".", b, ".", d, ".", q), rep(NA, n.imp))
        }
      }
    }
  }
  
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
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
      for (d in c("adj", "ow", "iptw")) {
        
        # get model for this analysis approach
        model <- get(paste0("fit.", d))
        
        for (p in 1:length(levels(cohort$egfr_cat))) {
          
          q = levels(as.factor(cohort$egfr_cat))[p]
          
          if (q == levels(as.factor(cohort$egfr_cat))[1]) {
            # coef_vector[i] <- model$coefficients[n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", d, ".", q, "`[", i, "] <- model$coefficients[", n, "]", collapse = "")
            # se_vector[i] <- sqrt(model$var[n,n])
            se_statement <- paste0("`SE.", drug_name, ".", d, ".", q, "`[", i, "] <- sqrt(model$var[", n, ",", n, "])", collapse = "")
          } else {
            # coef = coef1 + coef2
            # coef_vector[i] <- model$coefficients[n] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(p-1)+n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", d, ".", q, "`[", i, "] <- model$coefficients[", n, "] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(", p, "-1) + ", n, "]", collapse = "")
            
            # se = sqrt(var1 + var2 + cov1,2)
            # se_vector[i] <- sqrt(abs(model$var[n,n]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(p-1)+n,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(p-1)+n]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(p-1)+n]))
            se_statement <- paste0("`SE.", drug_name, ".", d, ".", q, "`[", i, "] <- sqrt(abs(model$var[", n, ",", n, "]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(", p, "-1) + ", n, ",nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(", p, "-1) + ", n, "]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$egfr_cat)[-1]))/(", p, "-1) + ", n, "]))", collapse = "")
          }          
          # execute commands
          eval(str2lang(coef_statement))
          eval(str2lang(se_statement))
          
          
          
          
        }
        
        rm(model)
      }
      
    }
    
    # save interaction significance
    if (k == "ckd_egfr50") {
      f_adjusted3 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, " + ", paste(covariates, collapse=" + "), collapse = ""))
      fit.ow.no_interaction <- coxph(f_adjusted3, cohort[cohort$.imp == i,], weights = overlap2)
      chi <- 2 * fit.ow$loglik[2] - 2 * fit.ow.no_interaction$loglik[2]
      p_value_interaction_egfr[i] <- 1 - pchisq(chi, df = 2)
    }
    
  }
  
  
  ## loop to pool and store results
  for (n in 1:length(drug_levels)) {
    
    # get drug name
    drug_name <- drug_levels[n]
    
    
    for (d in c("adj", "ow", "iptw")) {
      
      for (p in 1:length(levels(cohort$egfr_cat))) {
        q = levels(as.factor(cohort$egfr_cat))[p]
        if (n == 1) {
          # if drug level is reference category then HR will be NA
          pooled_hr <- c(1, NA, NA)
          pooled_hr_string <- "1.00 (ref.)"
          
        } else {
          
          # define names for objects containing pooled hr + 95% ci (as vector and as string)
          pooled_hr_name <- paste0(d, "_", drug_name, "_", q, "_hr")
          pooled_hr_string_name <- paste0(d, "_", drug_name, "_", q, "_string")
          
          # define names for coefficient vectors
          coef_name <- paste0("COEFS.", drug_name, ".", d, ".", q)
          se_name <- paste0("SE.", drug_name, ".", d, ".", q)
          
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
                                 egfr_cat = q,
                                 count = as.numeric(count %>% filter(egfr_cat==q) %>% select(paste0(drug_name, "_count"))),
                                 followup = as.numeric(followup %>% filter(egfr_cat==q) %>% select(paste0(drug_name, "_followup"))), 
                                 events = as.character(events %>% filter(egfr_cat==q) %>% select(paste0(drug_name, "_events"))),
                                 contrast = paste0(drug_name, " vs ", drug_levels[1], collapse = ""),
                                 variable = paste0("studydrug", m, collapse = ""),
                                 analysis = d,
                                 HR = pooled_hr[1],
                                 LB = pooled_hr[2],
                                 UB = pooled_hr[3],
                                 string = pooled_hr_string)
        
        # combine results by each analysis approach within each drug type (within studydrug variable)
        egfr_hrs <- rbind(egfr_hrs, outcome_hr)
      }
      
    }
    
  }
  
}    


# ensure hazard ratio and 95% ci are stored as numeric variables
class(egfr_hrs$HR) <- class(egfr_hrs$LB) <- class(egfr_hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
egfr_hrs <- egfr_hrs %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0(`events_number`, "/", `count`),
  )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(egfr_hrs, file=paste0(today, "_egfr_hrs.Rda"))

############################3 HAZARD RATIOS BY ALBUMINURIA CATEGORY################################################################

# create empty data frame to which we can append the hazard ratios once calculated
albuminuria_hrs <- data.frame()

# main dataset is large - for speed of computation we will only load in dataset we need each time
rm(cohort)
gc()


## calculate hazard ratios

# use variable studydrug2
m = 2

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

# create variable for albuminuria category
cohort <- cohort %>% mutate(
  albuminuria_cat = ifelse(uacr >30, "≥30", ifelse(uacr > 3, "3-30", "<3")),
  albuminuria_cat = factor(albuminuria_cat)
)

p_value_interaction_albuminuria <- rep(NA, n.imp)

# calculate hazard ratios per outcome
for (k in outcomes_per_drugclass) {
  
  print(paste0("Calculating event numbers per drug level for outcome ", k))
  
  
  
  # for other studydrug variables use regular censoring variables
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")
  
  
  # calculate number of subjects in each group
  count <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), albuminuria_cat) %>%
    summarise(count=round(n()/n.imp, 0)) %>% # the total number of subjects in the stacked imputed datasets has to be divided by the number of imputed datasets
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_count"),
                values_from=count)
  
  # calculate median follow up time (years) per group
  followup <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), albuminuria_cat) %>%
    summarise(time=round(median(!!sym(censtime_var)), 2)) %>%
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_followup"),
                values_from=time)
  
  # summarise number of events per group
  events <- cohort[cohort$.imp != 0,] %>%
    group_by(!!sym(studydrug_var), albuminuria_cat) %>%
    summarise(event_count=round(sum(!!sym(censvar_var))/n.imp, 0),
              drug_count=round(n()/n.imp, 0)) %>%
    mutate(events_perc=round(event_count*100/drug_count, 1),
           events=paste0(event_count, " (", events_perc, "%)")) %>%
    select(!!sym(studydrug_var), albuminuria_cat, events) %>%
    pivot_wider(names_from=!!sym(studydrug_var),
                names_glue=paste0("{studydrug", m, "}_events"),
                values_from=events)
  
  
  # write formulas for adjusted and unadjusted analyses
  f2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, "*albuminuria_cat", collapse = ""))
  
  f_adjusted2 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, "*albuminuria_cat", " + ", paste(covariates, collapse=" + "), collapse = ""))
  
  # create empty vectors to store the coefficients and standard errors of the hazard ratios from every imputed dataset
  
  for (b in drug_levels[-1]) {
    for (c in c("COEFS", "SE")) {
      for (d in c("adj", "ow", "iptw")) {
        for (q in levels(as.factor(cohort$albuminuria_cat))) {
          assign(paste0(c, ".", b, ".", d, ".", q), rep(NA, n.imp))
        }
      }
    }
  }
  
  
  for (i in 1:n.imp) {
    print(paste0("Analyses in imputed dataset number ", i))
    
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
      for (d in c("adj", "ow", "iptw")) {
        
        # get model for this analysis approach
        model <- get(paste0("fit.", d))
        
        for (p in 1:length(levels(cohort$albuminuria_cat))) {
          
          q = levels(as.factor(cohort$albuminuria_cat))[p]
          
          if (q == levels(as.factor(cohort$albuminuria_cat))[1]) {
            # coef_vector[i] <- model$coefficients[n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", d, ".", q, "`[", i, "] <- model$coefficients[", n, "]", collapse = "")
            # se_vector[i] <- sqrt(model$var[n,n])
            se_statement <- paste0("`SE.", drug_name, ".", d, ".", q, "`[", i, "] <- sqrt(model$var[", n, ",", n, "])", collapse = "")
          } else {
            # coef = coef1 + coef2
            # coef_vector[i] <- model$coefficients[n] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(p-1)+n]
            coef_statement <- paste0("`COEFS.", drug_name, ".", d, ".", q, "`[", i, "] <- model$coefficients[", n, "] + model$coefficients[length(model$coefficients)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(", p, "-1) + ", n, "]", collapse = "")
            
            # se = sqrt(var1 + var2 + cov1,2)
            # se_vector[i] <- sqrt(abs(model$var[n,n]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(p-1)+n,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(p-1)+n]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(p-1)+n]))
            se_statement <- paste0("`SE.", drug_name, ".", d, ".", q, "`[", i, "] <- sqrt(abs(model$var[", n, ",", n, "]) + abs(model$var[nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(", p, "-1) + ", n, ",nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(", p, "-1) + ", n, "]) + abs(2 * vcov(model)[1,nrow(model$var)-(length(drug_levels[-1])*length(levels(cohort$albuminuria_cat)[-1]))/(", p, "-1) + ", n, "]))", collapse = "")
          }          
          # execute commands
          eval(str2lang(coef_statement))
          eval(str2lang(se_statement))
          
          
          
          
        }
        
        rm(model)
      }
      
    }
    
    # save interaction significance
    if (k == "ckd_albuminuria50") {
      f_adjusted3 <- as.formula(paste0("Surv(", censtime_var, ", ", censvar_var, ") ~  studydrug", m, " + ", paste(covariates, collapse=" + "), collapse = ""))
      fit.ow.no_interaction <- coxph(f_adjusted3, cohort[cohort$.imp == i,], weights = overlap2)
      chi <- 2 * fit.ow$loglik[2] - 2 * fit.ow.no_interaction$loglik[2]
      p_value_interaction_albuminuria[i] <- 1 - pchisq(chi, df = 2)
    }
    
  }
  
  
  ## loop to pool and store results
  for (n in 1:length(drug_levels)) {
    
    # get drug name
    drug_name <- drug_levels[n]
    
    
    for (d in c("adj", "ow", "iptw")) {
      
      for (p in 1:length(levels(cohort$albuminuria_cat))) {
        q = levels(as.factor(cohort$albuminuria_cat))[p]
        if (n == 1) {
          # if drug level is reference category then HR will be NA
          pooled_hr <- c(1, NA, NA)
          pooled_hr_string <- "1.00 (ref.)"
          
        } else {
          
          # define names for objects containing pooled hr + 95% ci (as vector and as string)
          pooled_hr_name <- paste0(d, "_", drug_name, "_", q, "_hr")
          pooled_hr_string_name <- paste0(d, "_", drug_name, "_", q, "_string")
          
          # define names for coefficient vectors
          coef_name <- paste0("COEFS.", drug_name, ".", d, ".", q)
          se_name <- paste0("SE.", drug_name, ".", d, ".", q)
          
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
                                 albuminuria_cat = q,
                                 count = as.numeric(count %>% filter(albuminuria_cat==q) %>% select(paste0(drug_name, "_count"))),
                                 followup = as.numeric(followup %>% filter(albuminuria_cat==q) %>% select(paste0(drug_name, "_followup"))), 
                                 events = as.character(events %>% filter(albuminuria_cat==q) %>% select(paste0(drug_name, "_events"))),
                                 contrast = paste0(drug_name, " vs ", drug_levels[1], collapse = ""),
                                 variable = paste0("studydrug", m, collapse = ""),
                                 analysis = d,
                                 HR = pooled_hr[1],
                                 LB = pooled_hr[2],
                                 UB = pooled_hr[3],
                                 string = pooled_hr_string)
        
        # combine results by each analysis approach within each drug type (within studydrug variable)
        albuminuria_hrs <- rbind(albuminuria_hrs, outcome_hr)
      }
      
    }
    
  }
  
}    


# ensure hazard ratio and 95% ci are stored as numeric variables
class(albuminuria_hrs$HR) <- class(albuminuria_hrs$LB) <- class(albuminuria_hrs$UB) <- "numeric"

# create separate variables for events per number of drug initiations (nN)
albuminuria_hrs <- albuminuria_hrs %>% 
  separate(`events`, into = c("events_number", "events_percentage"), sep = " \\(", remove = FALSE) %>%
  mutate(
    `events_percentage` = str_replace(`events_percentage`, "\\)", ""),
    `nN` = paste0(`events_number`, "/", `count`),
  )

# store hazard ratios
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(albuminuria_hrs, file=paste0(today, "_albuminuria_hrs.Rda"))
