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
  
  temp <- temp_all %>% 
    filter(.imp != 0) %>% 
    select("patid", ".imp", contains("studydrug"), all_of(covariates), "dstartdate", "ncurrtx2")
  
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
    ungroup()
  
  cohort <- temp2 %>% cbind(temp %>% select(contains("overlap"), contains("IPTW")))
  
  rm(temp)
  rm(temp2)
  
  
  # save dataset with weights so this can be used in subsequent scripts
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  save(cohort, file=paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", m, ".Rda"))
  
  
}
rm(temp_all)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
load(paste0(today, "_t2d_glp1_imputed_data_withweights_studydrug", main, ".Rda"))

# get weighted baseline table

studydrug_var = paste0("studydrug", main)
weights_overlap = paste0("overlap", main)

# save weighted cohort
weighted_cohort <- cohort %>% filter(.imp != 0) %>%
  mutate(malesex=as.logical(malesex),
         statin=as.logical(statin),
         INS = as.logical(INS),
         ACE_or_ARB = as.logical(ACE_or_ARB),
         predrug_cvd = as.logical(predrug_cvd),
         predrug_heartfailure = as.logical(predrug_heartfailure),
         predrug_af = as.logical(predrug_af),
         hosp_admission_prev_year = as.logical(hosp_admission_prev_year)) %>%
  select(!!sym(studydrug_var), !!sym(weights_overlap), all_of(vars)) %>%
  mutate(count=round(!!sym(weights_overlap)*10000000),0) %>%
  uncount(count) 

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(weighted_cohort, file=paste0(today, "_weighted_imputed_data_for_table.Rda"))

rm(weighted_cohort)
rm(cohort)

# write function to extract summary statistics

cont <- function(x, var_name) {
  if (var_name %in% c("uacr", "dstartdate_dm_dur_all")) {
    # For specific variables, calculate Median (IQR)
    with(stats.apply.rounding(stats.default(x)), 
         c("Median (IQR)" = sprintf("%s (%s-%s)", 
                                    round_pad(as.numeric(MEDIAN), 1), 
                                    round_pad(as.numeric(Q1), 1), 
                                    round_pad(as.numeric(Q3), 1))))
  } else {
    # For other continuous variables, calculate Mean (SD)
    with(stats.apply.rounding(stats.default(x)), 
         c("Mean (SD)" = sprintf("%s (%s)", 
                                 round_pad(as.numeric(MEAN), 1), 
                                 round_pad(as.numeric(SD), 1))))
  }
}
missing <- function(x, ...) {
  with(stats.apply.rounding(stats.default(x)), c("Missing"=sprintf("%s", prettyNum(NMISS, big.mark=","))))
}


rndr <- function(x, name, ...) {
  if (is.logical(x)) {
    y <- render.default(x, name, ...)
    y[2]
  } else if (is.numeric(x)) {
    cont(x, name)  # pass both x and variable name to your cont() function
  } else {
    render.default(x, name, ...)
  }
}


strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s</span>", 
          label, prettyNum(n, big.mark=","))
}

cat <- function(x, ...) {
  vals <- stats.default(x)  # get raw stats without rounding
  c("", sapply(vals, function(y) {
    # assume y$PCT is numeric; convert if necessary
    sprintf("%.4f%%", as.numeric(y$PCT)) # 4 decimals
  }))
}

gc()


# Create chunks of 5 covariates
covariates_chunks <- split(covariates, ceiling(seq_along(covariates) / 4))

i = 1

for (covariates_chunk in covariates_chunks) {
  
  columns_to_load <- c(studydrug_var, covariates_chunk)
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
  load(paste0(today, "_weighted_imputed_data_for_table.Rda"))
  
  # Read just the needed columns
  data_chunk <- weighted_cohort %>% select(all_of(columns_to_load))
  
  rm(weighted_cohort)
  gc()
  
  # Build formula for table1
  formula <- as.formula(paste("~ ", paste(covariates_chunk, collapse = " + "), "| ", studydrug_var))
  
  # Run the table1 function on the selected chunk
  table <- table1(formula, data=data_chunk, overall=F, render=rndr, render.categorical=cat, render.continuous=cont, render.strat=strat)
  
  # calculate SMDs
  smd_table <- CreateTableOne(vars = covariates_chunk,
                              strata = studydrug_var,
                              data = data_chunk,
                              test = FALSE)
  
  # Extract SMDs
  smd <- ExtractSmd(smd_table)
  
  # save the output
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  write.table(as.data.frame(table), file = paste0(today, "_weighted_table_part_", i, ".csv"), sep = ";", dec = ",", row.names = F)
  write.table(as.data.frame(smd), file = paste0(today, "_weighted_smd_part_", i, ".csv"), sep = ";", dec = ",", row.names = F)
  
  # Clear the data from memory and run garbage collection to free memory
  rm(data_chunk)
  rm(table)
  i <- i + 1
  gc()
  
}
