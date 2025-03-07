########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data.Rda"))

############################1 CALCULATE WEIGHTS################################################################

# 1 calculate weights

# as the data have been imputed, take each imputed dataset, calculate weights in them, then stack them again at the end

# for computational speed, keep minimal dataset only
temp <- temp %>% filter(.imp != 0) %>% select("patid", ".imp", contains("studydrug"), all_of(covariates), "dstartdate")

# weight variables for each studydrug variable
n.studydrug.vars <- sum(grepl("studydrug", colnames(temp)))

for (m in 1:n.studydrug.vars) {
  
  # create empty variables for weights
  temp[[paste0("IPTW", m, collapse = "")]] <- temp[[paste0("overlap", m, collapse = "")]] <- NA
  
}  
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(temp, file=paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_0.Rda"))
rm(temp)
gc()

# calculate weights
for (m in 1:n.studydrug.vars) {
  
  studydrug_var = paste0("studydrug", m)
  
  print(paste0("Processing data for variable studydrug", m, collapse = ""))
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", m-1, ".Rda"))
  
  # remove observations that are duplicate for this studydrug level - these should not be included in weight calculations
  # (they will be kept NA and added at the end)
  temp_non_used_observations <- temp %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  temp <- temp %>% 
    group_by(.imp, patid, !!sym(studydrug_var)) %>% 
    arrange(dstartdate) %>% 
    filter(!duplicated(!!sym(studydrug_var))) %>% 
    ungroup()
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(temp_non_used_observations, file=paste0(today, "_t2d_glp1_redundant_observations_", m, ".Rda"))
  rm(temp_non_used_observations)
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
  load(paste0(today, "_t2d_glp1_redundant_observations_", m, ".Rda"))
  temp <- temp %>% rbind(temp_non_used_observations)
  
  save(temp, file=paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", m, ".Rda"))
  
  rm(temp)
  rm(temp_non_used_observations)
  gc()
}

# rename dataset with weights for merging
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_minimal_dataset_for_weight_calculation_", n.studydrug.vars, ".Rda"))
temp2 <- temp
rm(temp)

# add weights to imputed dataset with all variables (load this one first)
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data.Rda"))

# unimputed dataset
temp3 <- temp %>% filter(.imp == 0)
# imputed datasets but with all variables
temp <- temp %>% filter(.imp != 0)

for (m in 1:n.studydrug.vars) {
  
  # add empty variables for weights to unimputed dataset so that they can be combined later
  temp3[[paste0("IPTW", m, collapse = "")]] <- temp3[[paste0("overlap", m, collapse = "")]] <- NA
  
}  

# ensure rows in temp and temp2 (imputed datasets without and with weights) are in the same order
temp <- temp %>% 
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  ungroup()

temp2 <- temp2 %>% 
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  ungroup()

# add IPTW and overlap weights to full imputed dataset
temp <- temp %>% cbind(temp2 %>% select(contains("overlap"), contains("IPTW")))
  
# and combine with unimputed dataset
cohort <- rbind(temp3, temp)

rm(temp)
rm(temp2)
rm(temp3)

# save dataset with weights so this can be used in subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_withweights.Rda"))
#load(paste0(today, "_t2d_glp1_imputed_data_withweights.Rda"))