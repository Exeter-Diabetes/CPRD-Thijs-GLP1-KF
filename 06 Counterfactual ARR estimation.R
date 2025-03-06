########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")

########################1 PREPARE DATASET####################################################################

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_recalibrated.Rda"))

# select studydrug variable
m = 2
studydrug_var = paste0("studydrug", m)

centre_and_reference <- function(df, covariates) {
  df %>%
    mutate(across(all_of(covariates), ~ if(is.numeric(.)) . - mean(., na.rm = TRUE) else . ),  # Center numeric variables
           across(all_of(covariates), ~ if(is.logical(.)) . == names(sort(table(.), decreasing = TRUE))[1] else . ) ,  # Set most frequent level as reference for logical
           across(all_of(covariates), ~ if(is.factor(.)) relevel(., ref = names(sort(table(.), decreasing = TRUE))[1]) else . ))  # Set most frequent level as reference for factor
}

# create regex pattern of censoring variables to select
outcome_variables <- paste0("(", paste(outcomes_per_drugclass, collapse = "|"), ")(?!.*(5y|sens)).*(_censtime_yrs|_censvar)$")


cohort <- cohort %>%
  mutate(across(contains("predrug_"), as.logical),
         hosp_admission_prev_year=as.logical(hosp_admission_prev_year),
         INS=as.logical(INS),
         MFN=as.logical(MFN),
         malesex=as.factor(malesex),
         initiation_year=as.factor(initiation_year))  %>%
  #select relevant observations only (non-duplicated within studydrug variable)
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup() %>%
  #select relevant variables only
  select(patid, .imp, !!sym(studydrug_var), overlap2,
         matches(outcome_variables, perl = TRUE), 
         all_of(covariates)) %>% 
  centre_and_reference(covariates)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(cohort, file=paste0(today, "_data_centred_predictors.Rda"))

############################2 FIT WEIGHTED COX MODEL AND ESTIMATE COUNTERFACTUAL OBSERVED SURVIVAL################################################################

# clear R memory to ensure memory limit not exceeded
rm(list = setdiff(ls(), c("n.imp", "covariates", "today", "m")))

# if evaluating other risk scores, other relevant outcomes may be added
outcomes_per_drugclass <- "ckd_egfr50"

for (k in outcomes_per_drugclass) {
  
  for (i in 1:n.imp) {
    # load minimal dataset
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
    load(paste0(today, "_data_centred_predictors.Rda"))
    
    cohort <- cohort %>% filter(.imp == i)
    
    studydrug_var = paste0("studydrug", m)
    weights_overlap = paste0("overlap", m)
    drug_levels <- levels(cohort[[studydrug_var]])
    
    # clear R memory
    gc()
    
    censvar_var=paste0(k, "_censvar")
    censtime_var=paste0(k, "_censtime_yrs")
    
    # fit overlap-weighted cox model
    f_adjusted <- as.formula(paste("Surv(", censtime_var, ", ", censvar_var, ") ~  ", studydrug_var, " + ", paste(covariates, collapse=" + "))) 
    model <- cph(f_adjusted, data=cohort, x=TRUE, y=TRUE, surv=TRUE, weights=weights_overlap)
    
    for (n in 1:length(drug_levels)) { 
      
    drug_name <- drug_levels[n]

    observed_data_name <- paste0("observed_", drug_name)
    estimate_name <- paste0("estimate_", drug_name)
    se_name <- paste0("se_", drug_name)
    
    # create counterfactual dataset where everyone is treated with drug_name
    obs_data <- cohort %>%
      # create variable for factual study arm for reference later
      mutate(studydrug_original=!!sym(studydrug_var),
             !!sym(studydrug_var) := drug_name,
             rowno=row_number())
    
    print(paste0("Survival estimates for", drug_name, " in imputation ", i, "  (outcome ", k, ")"))
    
    # get counterfactual ("observed" weighted) survival
    observed <- survfit(model, newdata=as.data.frame(obs_data)) %>%
      tidy() %>%
      filter(time == 3) %>%
      pivot_longer(cols=-c(time, n.risk, n.event, n.censor), names_to = c(".value", "group"), names_pattern = "(.*)\\.(.*)") %>%
      select(group, estimate, std.error) %>%
      mutate(group=as.numeric(group)) %>%
      inner_join(obs_data, by=c("group"="rowno")) %>% 
      mutate(
        !!sym(estimate_name) := estimate,
        !!sym(se_name) := std.error) %>% 
      select(-c(estimate, std.error))
    
    # rename dataset
    assign(observed_data_name, observed)
    observed_data <- get(observed_data_name)
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
    # save results
    save(observed_data, file=paste0(today, "_adjusted_surv_",k,"_SGLT2i_imp.", i, ".Rda"))
    
    rm(obs_data)
    rm(observed)
    rm(observed_data)
    }
    # clear environment
    rm(list = setdiff(ls(), c("n.imp", "covariates", "k", "today", "m", "drug_levels", "outcomes_per_drugclass")))
  }
  

  # after separately estimating counterfactual survival probabilities in each imputation, we will now combine these in the main dataset

  for (n in 1:length(drug_levels)) { 
    
    drug_name <- drug_levels[n]
    
    observed_data_name <- paste0("observed_", drug_name)
    
    # create empty dataframes to append each imputed dataset to
    temp <- data.frame()
    
    # join estimates from each imputation in one dataframe
    for (i in 1:n.imp) {
      setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
      load(paste0(today, "_adjusted_surv_",k,"_", drug_name, "_imp.", i, ".Rda"))
      
      observed_data <- get(observed_data_name)
      
      temp <- temp %>% rbind(observed_data)
      
      rm(observed_data)
      
    }
    
    # rename dataset
    assign(observed_data_name, temp)
    observed_data <- get(observed_data_name)
    rm(temp)

    if (n == 1) {
      studydrug_var = paste0("studydrug", m)
      
      benefits <- observed_data
      benefits[[studydrug_var]] <- benefits$studydrug_original
      benefits <- benefits %>% select(-studydrug_original)

    } else {

      # combine counterfactual observed survival probabilities
      benefits <- observed_data %>% select(group, .imp, contains("estimate_"), contains("se_")) %>%
        inner_join(benefits, by = c("group", ".imp")) 
      
    }
  
  }
  
  
  # create variables for survival differences between each drug combination
  combinations <- combn(drug_levels, 2, simplify = FALSE)
  
  for (comb in combinations) {
    # obtain drug names in combination
    drug_name_1 <- comb[1]
    drug_name_2 <- comb[2]
    
    # define variable names of survival estimate and standard error per drug
    estimate_col_1 <- paste0("estimate_", drug_name_1)
    estimate_col_2 <- paste0("estimate_", drug_name_2)
    se_col_1 <- paste0("se_", drug_name_1)
    se_col_2 <- paste0("se_", drug_name_2)
    
    # add new variables for difference between two drug levels
    benefits <- benefits %>%
      mutate(
        !!paste0("survdiff_", drug_name_1, "_", drug_name_2) := 
          !!sym(estimate_col_1) - !!sym(estimate_col_2),
        
        !!paste0("se_survdiff_", drug_name_1, "_", drug_name_2) := 
          sqrt(!!sym(se_col_1)^2 + !!sym(se_col_2)^2)
      )
  }
  
  #rename variables by outcome name
  benefits <- benefits %>% rename_with(
    ~ paste0(.x, "_", k),
    contains("survdiff")
  )
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(benefits, file=paste0(today, "_adjusted_surv_",k,".Rda"))
  rm(benefits)
}

rm(list = setdiff(ls(), c("n.imp", "k", "today", "m", "outcomes_per_drugclass")))

### add counterfactual observed absolute risk reductions to main dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_recalibrated.Rda"))

studydrug_var = paste0("studydrug", m)

for (k in outcomes_per_drugclass) {
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_adjusted_surv_", k, ".Rda"))
  cohort <- cohort %>% left_join(benefits %>%
                                   select(.imp, patid, !!sym(studydrug_var), contains("survdiff")), 
                                 by=c(".imp", "patid", studydrug_var))
  rm(benefits)
}

cohort[[studydrug_var]] <- as.factor(cohort[[studydrug_var]])

# save dataset

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_with_observed_surv.Rda"))


