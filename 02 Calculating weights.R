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
  
  temp <- temp_all %>% filter(.imp != 0) %>% select("patid", ".imp", contains("studydrug"), all_of(covariates), "dstartdate")
  
  # create empty variables for weights
  temp[[paste0("IPTW", m, collapse = "")]] <- temp[[paste0("overlap", m, collapse = "")]] <- NA
  
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
  
  if (m != 2) {
    # only analyse dual treatment group in main analysis variable
    temp <- temp %>% filter(!!sym(studydrug_var) != "GLP1 + SGLT2") %>% 
      mutate(
        !!sym(studydrug_var) := droplevels(!!sym(studydrug_var))
      )
  }
  
  # force temp to be data.frame() for SumStat function
  temp <- temp %>% as.data.frame(temp)
  
  for (i in 1:n.imp) {
    
    print(paste0("Calculating weights for imputed dataset number ", i))
    
    # calculate overlap weights
    w.overlap <- SumStat(ps.formula=ps.formula,
                         data = temp[temp$.imp == i,],
                         weight = c("IPW", "overlap"))

    # truncate IPTW at 2nd and 98% percentile
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
  
  temp <- temp %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    ungroup()
  
  # ensure temp_all is in same row order as dataset with studydrug variable it will be merged with
  temp2 <- temp_all %>%
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    {if (m != 2) filter(., !!sym(studydrug_var) != "GLP1 + SGLT2") else .} %>%
    ungroup()
  
  cohort <- temp2 %>% cbind(temp %>% select(contains("overlap"), contains("IPTW")))
  
  rm(temp)
  rm(temp2)
  
  if (m == 1) {
    cohort <- cohort %>% mutate(
      !!sym(studydrug_var) := factor(!!sym(studydrug_var), levels = c("SU", "DPP4", "SGLT2", "GLP1"))
    )
  }
  
  if (m == 2) {
    
    cohort <- cohort %>% mutate(
      !!sym(studydrug_var) := factor(!!sym(studydrug_var), levels = c("DPP4 + SU", "GLP1", "SGLT2", "GLP1 + SGLT2"))
    )
    
    new_studydrug_var <- paste0("studydrug", n.studydrug.vars+1)
    new_weights_overlap = paste0("overlap", n.studydrug.vars+1)
    new_weights_iptw = paste0("IPTW", n.studydrug.vars+1)
    cohort[[new_studydrug_var]] <- relevel(cohort[[paste0("studydrug", m)]], ref = "SGLT2")
    cohort[[new_weights_overlap]] <- cohort[[paste0("overlap", m)]]
    cohort[[new_weights_iptw]] <- cohort[[paste0("IPTW", m)]]
    
    setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
    save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", n.studydrug.vars+1, ".Rda"))
  }
  
  if (m==3) {  
    cohort <- cohort %>% mutate(
      !!sym(studydrug_var) := factor(!!sym(studydrug_var), levels = c("DPP4 + SU", "SGLT2", "Oral semaglutide", "Subcutaneous semaglutide", "Other GLP1"))
    )
  }
  
  # save dataset with weights so this can be used in subsequent scripts
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  
}
rm(temp_all)
