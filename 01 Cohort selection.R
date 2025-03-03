########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")
########################1 COHORT SELECTION####################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition_kf function for details)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
load(paste0(today, "_t2d_1stinstance_a.Rda"))
load(paste0(today, "_t2d_1stinstance_b.Rda"))

t2d_1stinstance <- rbind(t2d_1stinstance_a, t2d_1stinstance_b)
rm(t2d_1stinstance_a)
rm(t2d_1stinstance_b)

load(paste0(today, "_t2d_all_drug_periods.Rda"))

# add variable for age and diabetes duration
t2d_1stinstance <- t2d_1stinstance %>% mutate(
  dstartdate_age=as.numeric(difftime(dstartdate, dob, units = "days")/365.25),
  dstartdate_dm_dur_all=as.numeric(difftime(dstartdate, dm_diag_date_all, units = "days")/365.25),
  malesex=ifelse(gender==1, T, F),
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/scripts/CPRD-Thijs-GLP1-KF/Functions/")
source("cohort_definition_kf.R")
cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

# "Number of subjects on GIPGLP1: 413 (not included due to small numbers)"
# "Number of subjects of starting an GLP1/SGLT2 or comparator drugs DPP4/SU between 2014-2024: 524510"
# "Number of drug episodes of starting an GLP1/SGLT2 or comparator drugs DPP4/SU between 2014-2024: 1251436"
# "Number of drug episodes excluded with established CVD: 307965"
# "Number of drug episodes excluded with established HF: 35635"
# "Number of drug episodes excluded with unknown CKD status: 142544"
# "Number of drug episodes excluded with established eGFR <20 mL/min/1.73m2 or ESKD: 2591"
# "Number of drug episodes removed (e.g. subsequent episode of starting DPP4/SU after episode of SGLT2 or GLP1): 34195"
# "Number of subjects included: 297069"
# "Number of drug episodes included: 728528"

table(cohort$studydrug1)




## B Make variables for survival analysis of all endpoints (see survival_variables_kf function for details)

source("survival_variables_kf.R")

cohort <- add_surv_vars(cohort, main_only=FALSE) # add per-protocol survival variables as well

rm(list=setdiff(ls(), c("cohort", "today", "vars", "factors", "nonnormal")))

# create acr variable for ckdpc risk scores that uses further source of acr if acr not available
cohort <- cohort %>% 
  mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         uacr=ifelse(uacr<0.6, 0.6, uacr),
         #and create variable to code whether someone is on oral hyperglycaemic agents
         oha=ifelse(ncurrtx > 1 | ncurrtx == 1 & INS == 0, 1L, 0L),
         statin=!is.na(predrug_latest_statins),
         ACEi=!is.na(predrug_latest_ace_inhibitors),
         ARB=!is.na(predrug_latest_arb),
         BB=!is.na(predrug_latest_beta_blockers),
         finerenone=!is.na(predrug_latest_finerenone),
         CCB=!is.na(predrug_latest_ca_channel_blockers),
         # ThZD=!is.na(predrug_latest_thiazide_diuretics),
         # loopD=!is.na(predrug_latest_loop_diuretics),
         # MRA=!is.na(predrug_latest_ksparing_diuretics),
         # steroids=!is.na(predrug_latest_oralsteroids),
         # immunosuppr=!is.na(predrug_latest_immunosuppressants),
         # osteoporosis=!is.na(predrug_latest_osteoporosis),
         # genital_infection=as.logical(predrug_medspecific_gi)
  )


cohort <- cohort %>%
  #default hba1c variable is from previous 6 months to index date
  #take hba1c within window of 2 years prior and 7 days post, similar as other biomarkers
  mutate(prehba1c = prehba1c2yrs) %>%
  select(patid, malesex, ethnicity_5cat, ethnicity_qrisk2, 
         # imd2015_10,  tds_2011,
         regstartdate, 
         # gp_end_date, death_date, 
         gp_end_date,
         drug_class, contains("studydrug"), dstartdate, dstopdate_class, drugline_all, drug_substance, ncurrtx, DPP4, GLP1, 
         MFN, SGLT2, SU, TZD, INS, dstartdate_age, dstartdate_dm_dur_all, preweight, height, prehba1c, prebmi, 
         prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, predbp, preegfr, preckdstage, 
         preacr, uacr, 
         # qrisk2_smoking_cat, 
         contains("cens"), last_sglt2_stop, last_glp1_stop, oha,
         ##add variables necessary to calculate qrisk2/qhdf and ckdpc scores
         predrug_fh_premature_cvd, predrug_af, predrug_rheumatoidarthritis,
         predrug_angina, predrug_myocardialinfarction, predrug_stroke, predrug_revasc,
         predrug_heartfailure, predrug_hypertension, 
         # predrug_acutepancreatitis,
         predrug_earliest_ace_inhibitors, predrug_earliest_arb,
         predrug_earliest_beta_blockers, predrug_latest_ace_inhibitors, 
         predrug_latest_arb, predrug_latest_beta_blockers, 
         predrug_earliest_ca_channel_blockers, predrug_latest_ca_channel_blockers, 
         #predrug_earliest_thiazide_diuretics, predrug_latest_thiazide_diuretics,
         # predrug_dka, predrug_falls, predrug_dementia, 
         # hosp_admission_prev_year,
         statin, ACEi, ARB, BB, finerenone, CCB, 
         # ThZD, loopD, MRA, steroids, immunosuppr, 
         # osteoporosis, genital_infection,
         ckd_egfr50_outcome_type, preacr_confirmed, preacr_previous, preacr_previous_date, preacr_next, preacr_next_date
  )


# # smoking status still missing
# cohort$qrisk2_smoking_cat <- as.factor(cohort$qrisk2_smoking_cat)

# set SU as reference group
cohort$studydrug1 <- relevel(as.factor(cohort$studydrug1), ref = "SU")

# create variable for year of treatment initiation
cohort$initiation_year <- substring(as.character(cohort$dstartdate), 1, 4)

# ethnicity cannot be calculated in the imputation model due to it being a constant variable
# for the sake of imputation, we will class missing as a separate category "missing" (5-cat ethnicity: 5; QRISK2: 10)
cohort <- cohort %>%
  mutate(ethnicity_qrisk2=ifelse(is.na(ethnicity_qrisk2), "10", ethnicity_qrisk2),
         ethnicity_5cat=ifelse(is.na(ethnicity_5cat), "5", ethnicity_5cat),
         ethnicity_5cat=factor(ethnicity_5cat,
                               levels = c(0, 1, 2, 3, 4, 5),
                               labels = c("White", "South Asian", "Black", "Other", "Mixed", "Not stated/Unknown"))) %>% 
  relocate(ethnicity_5cat, .after = last_col())

# # imd_2015 is not a continuous variable - we will categorise this as quantiles
# cohort <- cohort %>% mutate(
#   imd2015_10 = ifelse(imd2015_10 %in% c(1,2), "1/2",
#                       ifelse(imd2015_10 %in% c(3,4), "3/4",
#                              ifelse(imd2015_10 %in% c(5,6), "5/6",
#                                     ifelse(imd2015_10 %in% c(7,8), "7/8",
#                                            ifelse(imd2015_10 %in% c(9,10), "9/10", NA))))),
#   imd2015_10 = factor(imd2015_10)
# )

# variable preacr_confirmed indicates whether a person had their presence of albuminuria (3mg/mmol) confirmed on 2 readings
# this shows as NA if no second reading available to confirm - replace with NA
cohort <- cohort %>% mutate(preacr_confirmed = ifelse(is.na(preacr_confirmed), F, preacr_confirmed),
                            preacr_confirmed = ifelse(uacr<3, F, preacr_confirmed))


# save dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(cohort, file=paste0(today, "_t2d_glp1_cohort_data.Rda"))
#load(paste0(today, "_t2d_glp1_cohort_data.Rda"))
########################2 MULTIPLE IMPUTATION####################################################################

#dry run
ini <- mice(cohort, seed = 123, maxit = 0)

# remove observations if ncurrtx is missing (needs to be fixed in combo_start_stop)
cohort <- cohort %>% filter(!is.na(ncurrtx))

# there are a couple of variables that we do not need to impute, we can tell mice not to impute these
meth <- ini$meth

meth[c(
  # "death_date", 
  "preacr", "last_sglt2_stop", "last_glp1_stop", "preckdstage", 
  "dstopdate_class",
  "predrug_earliest_ace_inhibitors", 
       "predrug_earliest_arb",
       "predrug_earliest_beta_blockers", "predrug_earliest_ca_channel_blockers",
       "predrug_latest_ace_inhibitors", 
       "predrug_latest_arb",
       "predrug_latest_beta_blockers", "predrug_latest_ca_channel_blockers",
       "ethnicity_qrisk2", 
       # "predrug_earliest_thiazide_diuretics", "predrug_latest_thiazide_diuretics",
       "ckd_egfr50_outcome_type", "preacr_confirmed", 
       "preacr_previous", "preacr_previous_date", "preacr_next", "preacr_next_date")] <- ""

# # smoking status and deprivation missing at present
# meth[c("qrisk2_smoking_cat", "imd2015_10")] <- "polyreg"

meth[c("preweight", "height")] <- "pmm"

meth["prebmi"] <- "~ I( preweight / (height/100)^2)"

# use quickpred function to build predictor matrix
# we can specify which variables to definitely include (inlist) and which ones to leave out (outlist)

inlist <- c("malesex", "studydrug1",  "dstartdate_age",                # main sociodemographic factors
            # "imd2015_10",                                           # treatment variable
            "dstartdate_dm_dur_all", "prebmi", "pretotalcholesterol", # laboratory and vital sign measurements
            "presbp", "preegfr", "uacr", 
            # "qrisk2_smoking_cat", 
            "ckd_egfr50_censvar"
            #"death_censvar"     # outcome variables
)


#list variables that are 100% complete and are not interesting for the imputation model
complete_vars <- names(ini$nmis[ini$nmis == 0])
#inspect complete_vars by printing it > print(complete_vars) then choose variables that we want to omit
outlist1 <- c("patid", "gp_end_date", "drug_class", "drugline_all", "ncurrtx", 
              "DPP4", "GLP1", "SGLT2", "SU", "INS", 
              names(cohort)[grep("cens", names(cohort))],
              "oha", "predrug_angina", "predrug_myocardialinfarction", "predrug_stroke", 
              "predrug_revasc", "predrug_heartfailure", "initiation_year", "ethnicity_5cat",
              "last_glp1_stop", "last_sglt2_stop", "dstopdate_class",
              "ckd_egfr50_outcome_type", "preacr_confirmed", "preacr_previous", 
              "preacr_previous_date", "preacr_next", "preacr_next_date", 
              names(cohort)[apply(cohort, 2, function(x) any(is.na(x))) & meth == ""]) # any variables that have missing in them that are not being imputed

#list variables with outflux <0.5 
#outflux is an indicator of the potential usefulness for imputing other variables - 
#outflux depends on the proportion of missing data of the variable: 
#outflux of a completely observed variable is equal to 1, 
#whereas outflux of a completely missing variable is equal to 0
#for two variables having the same proportion of missing data, 
#the variable with higher outflux is better connected to the missing data, 
#and thus potentially more useful for imputing other variables)
fx <- flux(cohort)
outlist2 <- row.names(fx)[fx$outflux < 0.5]

#identify problematic variables from initial run (constant/collinear variables)
outlist3 <- as.character(ini$loggedEvents[, "out"])

#combine above variables in one list
outlist <- unique(c(outlist1, outlist2, outlist3))

pred <- quickpred(cohort, include = inlist, exclude = outlist)


# limit imputations to plausible range
post <- ini$post
post["preegfr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(60, 120))"
post["prebmi"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(20, 40))"
post["prehba1c"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(42, 97))"
post["uacr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0.6, 56.5))"
post["presbp"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(80, 180))"

n.imp <- 10

imp <- mice(data = cohort, 
            meth = meth, 
            pred = pred, 
            post = post, 
            m=n.imp, 
            seed = 123)

# check imputed vs original values (disabled as time-consuming)
#densityplot(x = imp, data = ~ imd2015_10 + dstartdate_dm_dur_all + preweight + height + prehba1c + prebmi + 
#              prehdl + preldl + pretriglyceride + pretotalcholesterol + prealt + presbp + predbp + preegfr + 
#              uacr + qrisk2_smoking_cat + tds_2011)

#extract imputations so we can add variables
temp <- complete(imp, action = "long", include = T)

#the post-processing of prebmi imputations does not work as the method is passive imputation
#I will set all imputed prebmi values that are <20 at 20 which is what the post-processing procedure would otherwise do
z <- cohort[is.na(cohort$prebmi),]$patid
temp <- temp %>% mutate(prebmi = ifelse(
  patid %in% z & prebmi < 20, 20, prebmi
))
rm(z)


# add ckd pc risk score

temp <- temp %>%
  
  mutate(preckdstage=ifelse(is.na(preckdstage), ifelse(preegfr<15, "stage_5", 
                                                       ifelse(preegfr<30, "stage_4", 
                                                              ifelse(preegfr<45, "stage_3b", 
                                                                     ifelse(preegfr<60, "stage_3a", 
                                                                            ifelse(preegfr<90, "stage_2", 
                                                                                   "stage_1"))))), preckdstage),
         
         black_ethnicity=ifelse(!is.na(ethnicity_qrisk2) & (ethnicity_qrisk2 == 6 | ethnicity_qrisk2 == 7), 
                                1L, 
                                ifelse(is.na(ethnicity_qrisk2), NA, 0L)),
         
         cvd=predrug_myocardialinfarction==1 | predrug_revasc==1 | predrug_heartfailure==1 | predrug_stroke==1,
         
         # ever_smoker=ifelse(!qrisk2_smoking_cat == 0, 1L, 0L),
         # current_smoker=ifelse(qrisk2_smoking_cat==2 | qrisk2_smoking_cat == 3 | qrisk2_smoking_cat == 4, 1L, 0L),
         # ex_smoker=ifelse(qrisk2_smoking_cat==1, 1L, 0L),
         current_smoker = 0L,
         ex_smoker=0L,
         
         latest_bp_med=pmax(
           ifelse(is.na(predrug_latest_ace_inhibitors),as.Date("1900-01-01"),predrug_latest_ace_inhibitors),
           ifelse(is.na(predrug_latest_arb),as.Date("1900-01-01"),predrug_latest_arb),
           ifelse(is.na(predrug_latest_beta_blockers),as.Date("1900-01-01"),predrug_latest_beta_blockers),
           ifelse(is.na(predrug_latest_ca_channel_blockers),as.Date("1900-01-01"),predrug_latest_ca_channel_blockers),
  #         ifelse(is.na(predrug_latest_thiazide_diuretics),as.Date("1900-01-01"),predrug_latest_thiazide_diuretics),
           na.rm=TRUE
         ) %>% as.Date(),
         
         bp_meds_ckdpc=ifelse(latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=183, 1L, 0L),
         
         hypertension=ifelse((!is.na(presbp) & presbp>=140) | (!is.na(predbp) & predbp>=90) | bp_meds_ckdpc==1, 1L,0L),
         
         chd=predrug_myocardialinfarction==1 | predrug_revasc==1,
         
         )

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/Functions/")
source("calculate_ckdpc_50egfr_risk.R")

temp <- temp %>% 
  
  mutate(sex=ifelse(malesex == T, "male", ifelse(malesex==F, "female", NA))) %>%
  # calculate_ckdpc_50egfr_risk(age=dstartdate_age, 
  #                             sex=sex, 
  #                             egfr=preegfr, 
  #                             acr=uacr, 
  #                             sbp=presbp, 
  #                             bp_meds=bp_meds_ckdpc, 
  #                             hf=predrug_heartfailure, 
  #                             chd=chd, 
  #                             af=predrug_af, 
  #                             current_smoker=current_smoker, 
  #                             ex_smoker=ex_smoker, 
  #                             bmi=prebmi, 
  #                             hba1c=prehba1c, 
  #                             oha=oha, 
  #                             insulin=INS) 
  calculate_ckdpc_50egfr_risk(age=dstartdate_age, 
                              sex=sex, 
                              egfr=preegfr, 
                              acr=uacr, 
                              sbp=presbp, 
                              bp_meds=bp_meds_ckdpc, 
                              hf=predrug_heartfailure, 
                              chd=chd, 
                              af=predrug_af, 
                              current_smoker=current_smoker, 
                              ex_smoker=ex_smoker, 
                              bmi=prebmi, 
                              hba1c=prehba1c, 
                              oha=oha, 
                              insulin=INS) 


temp <- temp %>%
  
  mutate(across(starts_with("ckdpc_50egfr"),
                ~ifelse(dstartdate_age>=20 & dstartdate_age<=80 &
                          prebmi>=20, .x, NA)))

q <- temp %>% filter(.imp !=0 & (dstartdate_age<20 | dstartdate_age>80 | prebmi < 20))

# those left out:
q1 <- q %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q1))

q2 <- q %>% nrow()
print(paste0("Number of drug episodes excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q2/n.imp))


# retain those with available risk scores only
temp <- temp %>% filter(!is.na(ckdpc_50egfr_score))


temp <- temp %>% mutate(
  obesity = ifelse(prebmi < 30, F, T),
#  smoking_hx = ifelse(qrisk2_smoking_cat == 0, F, T),
#  smoking_status = ifelse(qrisk2_smoking_cat == 0, "never", ifelse(qrisk2_smoking_cat == 1, "ex", "current")),
  albuminuria_unconfirmed = ifelse(uacr < 3, F, T),
  albuminuria = preacr_confirmed,        # 
  ACEi_or_ARB = ifelse(temp$ACEi + temp$ARB > 0, T, F),
  ncurrtx = ifelse(ncurrtx==1, "1.", ifelse(ncurrtx==2, "2.", "3+"))
)


q <- temp %>% filter(!.imp == 0) %>% nrow()
p <- temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
  group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup() %>% nrow()
print(paste0("Number of duplicate drug episodes removed ", (q-p)/n.imp))
print(paste0("Number of drug episodes in study population ", p/n.imp))
rm(p)
q <- temp %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects in study population ", q))

# save imputed dataset so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(temp, file=paste0(today, "_t2d_glp1_imputed_data.Rda"))


# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)

table <- CreateTableOne(vars = vars, strata = "studydrug2", data = temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
                          group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup(),  ## these "duplicate" episodes will be removed after we have done the drug-specific analyses
                        factorVars = factors, test = F)

tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
## Save to a CSV file
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv

write.csv2(tabforprint, file = paste0(today, "_baseline_table.csv"))

# table by different type of semaglutide
table2 <- CreateTableOne(vars = vars, strata = "studydrug1", data = temp %>% filter(!.imp == 0) %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
                          group_by(.imp, patid) %>% filter(!duplicated(studydrug1)) %>% ungroup(),  ## these "duplicate" episodes will be removed after we have done the drug-specific analyses
                        factorVars = factors, test = F)

tabforprint2 <- print(table2, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = T)
## Save to a CSV file
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
#my computer is set to continental settings, therefore I am using write.csv2 instead of write.csv

write.csv2(tabforprint2, file = paste0(today, "_baseline_table2.csv"))



# events rates (sum of events divided by sum of person-years) by studydrug
outcomes <- c("ckd_egfr40", "ckd_egfr50", "ckd_egfr50_5y"#, "macroalb", "dka", "side_effect", "death", "amputation"
              )

for (k in outcomes) {
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")  
  
  for (m in levels(as.factor(temp$studydrug1))) {
    events <- temp %>% filter(.imp !=0 & studydrug1 == m) %>% select(censvar_var) %>% sum()
    pyears <- temp %>% filter(.imp !=0 & studydrug1 == m) %>% select(censtime_var) %>% sum()
    print(paste0(m, " event rate for ", k, ": ", round(events/pyears*1000,1), " per 1000 patient-years"))
    rm(events)
    rm(pyears)
  }
  rm(censvar_var)
  rm(censtime_var)
}

sensitivity_outcomes <- "ckd_egfr50_pp"

for (k in sensitivity_outcomes) {
  censvar_var=paste0(k, "_censvar")
  censtime_var=paste0(k, "_censtime_yrs")  
  
  for (m in levels(as.factor(temp$studydrug1))) {
    events <- temp %>% filter(.imp !=0 & studydrug1 == m) %>% select(censvar_var) %>% sum()
    pyears <- temp %>% filter(.imp !=0 & studydrug1 == m) %>% select(censtime_var) %>% sum()
    print(paste0(m, " event rate for ", k, ": ", round(events/pyears*1000,1), " per 1000 patient-years"))
    rm(events)
    rm(pyears)
  }
  rm(censvar_var)
  rm(censtime_var)
}
