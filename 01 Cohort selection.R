########################0 SETUP####################################################################
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-GLP1-KF/")
source("00 Setup.R")
########################1 COHORT SELECTION####################################################################

# 1 Cohort selection and variable setup

## A Cohort selection (see cohort_definition_kf function for details)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
load(paste0("2025-10-01_t2d_1stinstance_a.Rda"))
load(paste0("2025-10-01_t2d_1stinstance_b.Rda"))

t2d_1stinstance <- rbind(t2d_1stinstance_a, t2d_1stinstance_b)
rm(t2d_1stinstance_a)
rm(t2d_1stinstance_b)

load(paste0("2025-10-01_t2d_all_drug_periods.Rda"))

# add variable for age and diabetes duration
t2d_1stinstance <- t2d_1stinstance %>% mutate(
  dstartdate_age=as.numeric(difftime(dstartdate, dob, units = "days")/365.25),
  dstartdate_dm_dur_all=as.numeric(difftime(dstartdate, dm_diag_date_all, units = "days")/365.25),
  malesex=ifelse(gender==1, T, F),
)

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/scripts/CPRD-Thijs-GLP1-KF/Functions/")
source("cohort_definition_kf.R")
cohort <- define_cohort(t2d_1stinstance, t2d_all_drug_periods)

# [1] "Number of subjects on GIPGLP1: 259 (not included due to small numbers)" 
# [1] "Number of subjects already on SGLT2 starting an GLP1 or comparator drugs DPP4/SU between 2013-2023: 39890" 
# [1] "Number of drug episodes of already on SGLT2 starting an GLP1 or comparator drugs DPP4/SU between 2014-2023: 49876" 
# [1] "Number of drug episodes excluded with unknown CKD status: 10716" 
# [1] "Number of drug episodes excluded with established eGFR <20 mL/min/1.73m2 or ESKD: 92" 
# [1] "Number of drug episodes removed (e.g. subsequent episode of starting DPP4/SU after already taking the other): 427" 
# [1] "Number of subjects included: 31651" 
# [1] "Number of drug episodes included: 38643"


table(cohort$studydrug1)

# SGLT2 + SU SGLT2 + DPP4 SGLT2 + GLP1 
# 51390        84820       200390

rm(t2d_1stinstance)
rm(t2d_all_drug_periods)
gc()

## B Make variables for survival analysis of all endpoints (see survival_variables_kf function for details)

source("survival_variables_kf.R")

cohort <- add_surv_vars(cohort, main_only=FALSE) # add per-protocol survival variables as well

rm(list=setdiff(ls(), c("cohort", "today", "vars", "factors", "nonnormal", "main", "outcomes")))

# create acr variable for ckdpc risk scores that uses further source of acr if acr not available
cohort <- cohort %>% 
  mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         uacr=ifelse(uacr<0.6, 0.6, uacr),
         #and create variable to code whether someone is on oral hyperglycaemic agents
         oha=ifelse(ncurrtx > 1 | ncurrtx == 1 & INS == 0, 1L, 0L),
         statin=!is.na(predrug_latest_statins),
         ACE=!is.na(predrug_latest_ace_inhibitors),
         ARB=!is.na(predrug_latest_arb),
         BB=!is.na(predrug_latest_beta_blockers),
         finerenone=!is.na(predrug_latest_finerenone),
         CCB=!is.na(predrug_latest_ca_channel_blockers),
         ThZD=!is.na(predrug_latest_thiazide_diuretics),
         loopD=!is.na(predrug_latest_loop_diuretics),
         MRA=!is.na(predrug_latest_ksparing_diuretics),
         predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1, 1, 0),
  )


cohort <- cohort %>%
  #default hba1c variable is from previous 6 months to index date
  #take hba1c within window of 2 years prior and 7 days post, similar as other biomarkers
  mutate(prehba1c = prehba1c2yrs) %>%
  select(patid, malesex, ethnicity_5cat, ethnicity_qrisk2, 
         imd_decile, tds_2011,
         regstartdate, 
         gp_end_date, death_date, 
         gp_end_date,
         drug_class, contains("studydrug"), dstartdate, dstopdate_class, drugline_all, drug_substance, ncurrtx, ncurrtx2,
         DPP4, GLP1, MFN, SGLT2, SU, TZD, INS, 
         dstartdate_age, dstartdate_dm_dur_all, preweight, height, prehba1c, prebmi, 
         prehdl, preldl, pretriglyceride, pretotalcholesterol, prealt, presbp, predbp, preegfr, preckdstage, 
         preacr, uacr, 
         qrisk2_smoking_cat, 
         contains("cens"), last_sglt2_stop, last_glp1_stop, oha, timeprevcombo_class,
         ##add variables necessary to calculate qrisk2/qhdf and ckdpc scores
         predrug_af,
         predrug_angina, predrug_myocardialinfarction, predrug_stroke, predrug_revasc,
         predrug_heartfailure, predrug_cvd,
         predrug_hypertension, 
         predrug_acutepancreatitis,
         predrug_earliest_ace_inhibitors, predrug_earliest_arb,
         predrug_earliest_beta_blockers, predrug_latest_ace_inhibitors, 
         predrug_latest_arb, predrug_latest_beta_blockers, 
         predrug_earliest_ca_channel_blockers, predrug_latest_ca_channel_blockers, 
         predrug_earliest_thiazide_diuretics, predrug_latest_thiazide_diuretics,
         predrug_rheumatoidarthritis, predrug_fh_premature_cvd, 
         hosp_admission_prev_year, predrug_efi_score,
         statin, ACE, ARB, BB, finerenone, CCB, 
         ThZD, loopD, MRA, 
         ckd_egfr40_outcome_type, preacr_confirmed, preacr_previous, preacr_previous_date, preacr_next, preacr_next_date
  )

# set SU as reference group
cohort$studydrug1 <- relevel(as.factor(cohort$studydrug1), ref = "SGLT2 + SU")
cohort$studydrug2 <- relevel(as.factor(cohort$studydrug2), ref = "SGLT2 + DPP4/SU")
cohort$studydrug3 <- relevel(as.factor(cohort$studydrug3), ref = "SGLT2 + DPP4/SU")

# create variable for year of treatment initiation
cohort$initiation_year <- substring(as.character(cohort$dstartdate), 1, 4)
cohort$initiation_year <- as.numeric(cohort$initiation_year)

# Set start_year to the minimum year in the data
start_year <- min(cohort$initiation_year)

# Create a new variable that groups the years into two-year periods
cohort$initiation_year <- paste0(
  floor((cohort$initiation_year - start_year) / 2) * 2 + start_year, "/", 
  floor((cohort$initiation_year - start_year) / 2) * 2 + start_year + 1
)

# ethnicity cannot be calculated in the imputation model due to it being a constant variable
# for the sake of imputation, we will class missing as a separate category "missing" (5-cat ethnicity: 5; QRISK2: 10)
cohort <- cohort %>%
  mutate(ethnicity_qrisk2=ifelse(is.na(ethnicity_qrisk2), "10", as.character(ethnicity_qrisk2)),
         ethnicity_4cat=ifelse(is.na(ethnicity_5cat) | ethnicity_5cat == "4", "3", as.character(ethnicity_5cat)),
         ethnicity_4cat=factor(ethnicity_4cat,
                               levels = c(0, 1, 2, 3),
                               labels = c("White", "South Asian", "Black", "Other or unknown"))) %>% 
  relocate(ethnicity_4cat, .after = last_col())

# imd_decile is not a continuous variable - we will categorise this as quintiles
cohort <- cohort %>% mutate(
  imd_decile = ifelse(imd_decile %in% c(1,2), "1/2",
                      ifelse(imd_decile %in% c(3,4), "3/4",
                             ifelse(imd_decile %in% c(5,6), "5/6",
                                    ifelse(imd_decile %in% c(7,8), "7/8",
                                           ifelse(imd_decile %in% c(9,10), "9/10", NA))))),
  imd_decile = factor(imd_decile)
)

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
  "death_date", 
  "preacr", "last_sglt2_stop", "last_glp1_stop", "preckdstage", 
  "dstopdate_class",
  "predrug_earliest_ace_inhibitors", 
       "predrug_earliest_arb",
       "predrug_earliest_beta_blockers", "predrug_earliest_ca_channel_blockers",
       "predrug_latest_ace_inhibitors", 
       "predrug_latest_arb",
       "predrug_latest_beta_blockers", "predrug_latest_ca_channel_blockers",
       "ethnicity_qrisk2", 
       "predrug_earliest_thiazide_diuretics", "predrug_latest_thiazide_diuretics",
       "ckd_egfr40_outcome_type", "preacr_confirmed", 
       "preacr_previous", "preacr_previous_date", "preacr_next", "preacr_next_date")] <- ""

# # smoking status and deprivation missing at present
# meth[c("qrisk2_smoking_cat", "imd_decile")] <- "polyreg"

meth[c("preweight", "height")] <- "pmm"

meth["prebmi"] <- "~ I( preweight / (height/100)^2)"

# use quickpred function to build predictor matrix
# we can specify which variables to definitely include (inlist) and which ones to leave out (outlist)

inlist <- c("malesex",  "dstartdate_age",  "imd_decile",  "tds_2011",            # main sociodemographic factors
            paste0("studydrug", main),                                           # treatment variable
            "dstartdate_dm_dur_all", "prebmi", "pretotalcholesterol",            # laboratory and vital sign measurements
            "presbp", "preegfr", "uacr", 
            "qrisk2_smoking_cat", 
            "ckd_egfr40_censvar",
            "death_censvar"     # outcome variables
)


#list variables that are 100% complete and are not interesting for the imputation model
complete_vars <- names(ini$nmis[ini$nmis == 0])
#inspect complete_vars by printing it > print(complete_vars) then choose variables that we want to omit
outlist1 <- c("patid", "gp_end_date", "drug_class", "drugline_all", "ncurrtx", 
              "DPP4", "GLP1", "SGLT2", "SU", "INS", 
              names(cohort)[grep("cens", names(cohort))],
              "oha", "predrug_angina", "predrug_myocardialinfarction", "predrug_stroke", 
              "predrug_revasc", "predrug_heartfailure", "initiation_year", "ethnicity_4cat",
              "last_glp1_stop", "last_sglt2_stop", "dstopdate_class",
              "ckd_egfr40_outcome_type", "preacr_confirmed", "preacr_previous", 
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
#densityplot(x = imp, data = ~ imd_decile + dstartdate_dm_dur_all + preweight + height + prehba1c + prebmi + 
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

temp <- temp %>% filter(!.imp == 0) %>% 
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
         
         ever_smoker=ifelse(!qrisk2_smoking_cat == 0, 1L, 0L),
         current_smoker=ifelse(qrisk2_smoking_cat==2 | qrisk2_smoking_cat == 3 | qrisk2_smoking_cat == 4, 1L, 0L),
         ex_smoker=ifelse(qrisk2_smoking_cat==1, 1L, 0L),

         latest_bp_med=pmax(
           ifelse(is.na(predrug_latest_ace_inhibitors),as.Date("1900-01-01"),predrug_latest_ace_inhibitors),
           ifelse(is.na(predrug_latest_arb),as.Date("1900-01-01"),predrug_latest_arb),
           ifelse(is.na(predrug_latest_beta_blockers),as.Date("1900-01-01"),predrug_latest_beta_blockers),
           ifelse(is.na(predrug_latest_ca_channel_blockers),as.Date("1900-01-01"),predrug_latest_ca_channel_blockers),
           ifelse(is.na(predrug_latest_thiazide_diuretics),as.Date("1900-01-01"),predrug_latest_thiazide_diuretics),
           na.rm=TRUE
         ) %>% as.Date(),
         
         bp_meds_ckdpc=ifelse(latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=183, 1L, 0L),
         
         hypertension=ifelse((!is.na(presbp) & presbp>=140) | (!is.na(predbp) & predbp>=90) | bp_meds_ckdpc==1, 1L,0L),
         
         chd=predrug_myocardialinfarction==1 | predrug_revasc==1,
         
         )




setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2023/scripts/CPRD-Thijs-SGLT2-KF-scripts/Functions/")

source("calculate_ckdpc_40egfr_risk.R")

temp <- temp %>%
  mutate(sex=ifelse(malesex == T, "male", ifelse(malesex==F, "female", NA))) %>%
  
  calculate_ckdpc_40egfr_risk(age=dstartdate_age,
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


source("calculate_ckdpc_50egfr_risk.R")

temp <- temp %>% 
  
  mutate(sex=ifelse(malesex == T, "male", ifelse(malesex==F, "female", NA))) %>%

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


# temp <- temp %>%
#   
#   mutate(across(starts_with("ckdpc_50egfr"),
#                 ~ifelse(dstartdate_age>=20 & dstartdate_age<=80 &
#                           prebmi>=20, .x, NA)))
# 

q <- temp %>% filter(.imp !=0 & (dstartdate_age<20 | dstartdate_age>80 | prebmi < 20))

# those left out:
q1 <- q %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q1))

q2 <- q %>% nrow()
print(paste0("Number of drug episodes excluded with missing ckdpc risk scores due to age/BMI/HbA1c/uACR/SBP out of range: ", q2/n.imp))


# retain those with available CKD risk scores only
temp <- temp %>% filter(!is.na(ckdpc_50egfr_score))


# add qrisk2 score

temp <- temp %>% mutate(ethnicity_qrisk2 = ifelse(ethnicity_qrisk2 == "10", "9", ethnicity_qrisk2),
                        precholhdl=pretotalcholesterol/prehdl,
                        ckd45=ifelse(is.na(preckdstage), ifelse(preegfr < 30, T, F), 
                                     ifelse(preckdstage == "stage_4" | preckdstage == "stage_5", T, F)),
                        cvd=predrug_myocardialinfarction==1 | predrug_angina==1 | predrug_stroke==1,
                        sex=ifelse(malesex==1, "male", "female"),
                        dm_duration_cat=ifelse(dstartdate_dm_dur_all<=1, 0L,
                                               ifelse(dstartdate_dm_dur_all<4, 1L,
                                                      ifelse(dstartdate_dm_dur_all<7, 2L,
                                                             ifelse(dstartdate_dm_dur_all<11, 3L, 4L)))),
                        
                        earliest_bp_med=pmin(
                          if_else(is.na(predrug_earliest_ace_inhibitors),as.Date("2050-01-01"),predrug_earliest_ace_inhibitors),
                          if_else(is.na(predrug_earliest_arb),as.Date("2050-01-01"),predrug_earliest_arb),
                          if_else(is.na(predrug_earliest_beta_blockers),as.Date("2050-01-01"),predrug_earliest_beta_blockers),
                          if_else(is.na(predrug_earliest_ca_channel_blockers),as.Date("2050-01-01"),predrug_earliest_ca_channel_blockers),
                          if_else(is.na(predrug_earliest_thiazide_diuretics),as.Date("2050-01-01"),predrug_earliest_thiazide_diuretics),
                          na.rm=TRUE
                        ),
                        latest_bp_med=pmax(
                          if_else(is.na(predrug_latest_ace_inhibitors),as.Date("1900-01-01"),predrug_latest_ace_inhibitors),
                          if_else(is.na(predrug_latest_arb),as.Date("1900-01-01"),predrug_latest_arb),
                          if_else(is.na(predrug_latest_beta_blockers),as.Date("1900-01-01"),predrug_latest_beta_blockers),
                          if_else(is.na(predrug_latest_ca_channel_blockers),as.Date("1900-01-01"),predrug_latest_ca_channel_blockers),
                          if_else(is.na(predrug_latest_thiazide_diuretics),as.Date("1900-01-01"),predrug_latest_thiazide_diuretics),
                          na.rm=TRUE
                        ),
                        bp_meds_qrisk2=ifelse(earliest_bp_med!=as.Date("2050-01-01") & latest_bp_med!=as.Date("1900-01-01") & difftime(dstartdate, latest_bp_med, units="days")<=28 & earliest_bp_med!=latest_bp_med, 1L, 0L),
                        
                        type1=0L,
                        type2=1L,
                        surv_5yr=5L)

# calculate qrisk score

temp <- temp %>%
  
  mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  
  calculate_qrisk2(sex=sex2, age=dstartdate_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=predrug_fh_premature_cvd, renal=ckd45, af=predrug_af, rheumatoid_arth=predrug_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds_qrisk2, town=tds_2011, surv=surv_5yr) %>%
  
  rename(qrisk2_score_5yr=qrisk2_score) 



temp <- temp %>% mutate(
  obesity = ifelse(prebmi < 30, F, T),
  smoking_hx = ifelse(qrisk2_smoking_cat == 0, F, T),
  smoking_status = ifelse(qrisk2_smoking_cat == 0, "never", ifelse(qrisk2_smoking_cat == 1, "ex", "current")),
  albuminuria_unconfirmed = ifelse(uacr < 3, F, T),
  albuminuria = preacr_confirmed,        # 
  ACE_or_ARB = ifelse(temp$ACE + temp$ARB > 0, T, F),
  ncurrtx2 = ncurrtx,
  ncurrtx = ifelse(ncurrtx==1, "1.", ifelse(ncurrtx==2, "2.", ifelse(ncurrtx == 3, "3.", "4+"))),
  ncurrtx = relevel(as.factor(ncurrtx), ref = "3."),
  predrug_efi_cat = case_when(
    predrug_efi_score < 0.12 ~ "fit",
    predrug_efi_score >= 0.12 & predrug_efi_score < 0.24 ~ "mild",
    predrug_efi_score >= 0.24 & predrug_efi_score < 0.36 ~ "moderate",
    predrug_efi_score >= 0.36 ~ "severe"
  )
)


q <- temp %>% nrow()
p <- temp %>%  ## dataset at present contains separate drug episodes if a subject started a DPP4i and later a sulfonylurea
  group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup() %>% nrow()
print(paste0("Number of duplicate drug episodes removed ", (q-p)/n.imp))
print(paste0("Number of drug episodes in study population ", p/n.imp))
rm(p)
q <- temp %>% .$patid %>% unique() %>% length()
print(paste0("Number of subjects in study population ", q))

studydrug_var = paste0("studydrug", main)
# save imputed dataset so this can be used in the subsequent scripts
temp <- temp %>%
       mutate(across(starts_with("studydrug"), as.factor),
              egfr_cat = ifelse(preegfr < 45, "20-45", ifelse(preegfr < 60, "45-60", "≥60")),
              egfr_cat = factor(egfr_cat),
              albuminuria_cat = ifelse(uacr >30, "≥30", ifelse(uacr > 3, "3-30", "<3")),
              albuminuria_cat = factor(albuminuria_cat)) %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup()


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(temp, file=paste0(today, "_t2d_glp1_imputed_data.Rda"))


# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)
n.studydrug.vars <- sum(grepl("studydrug", colnames(temp)))

for (m in 1:n.studydrug.vars) {
  
  studydrug_var = paste0("studydrug", m)
  
  table <- CreateTableOne(vars = vars, strata = studydrug_var, data = temp  %>% 
                            group_by(.imp, patid) %>% filter(!duplicated(!!sym(studydrug_var))) %>% ungroup(),  
                          factorVars = factors, test = F)
  
  tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, printToggle = T)
  
  ## Save to a CSV file
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  write.csv2(tabforprint, file = paste0(today, "_baseline_table_studydrug", m, ".csv"))
  
}


# events rates (sum of events divided by sum of person-years) by studydrug

for (m in 1:n.studydrug.vars) {  
  
  studydrug_var = paste0("studydrug", m)
  
  print(paste0("Event rates for ", studydrug_var))
  
    
    for (k in outcomes) {
      
      if (m == 1) {
        k = paste0(k, "_sens1")
      }
      
      censvar_var=paste0(k, "_censvar")
      censtime_var=paste0(k, "_censtime_yrs")  
      
      for (p in levels(temp[[studydrug_var]])) {
        events <- temp %>% filter(.imp !=0 & !!sym(studydrug_var) == p) %>% select(censvar_var) %>% sum()
        pyears <- temp %>% filter(.imp !=0 & !!sym(studydrug_var) == p) %>% select(censtime_var) %>% sum()
        print(paste0(p, " event rate for ", k, ": ", round(events/pyears*1000,1), " per 1000 patient-years"))
        rm(events)
        rm(pyears)
      }
      rm(censvar_var)
      rm(censtime_var)
    }
    
}
