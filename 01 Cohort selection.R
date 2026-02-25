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

# [1] "Number of subjects on GIP-GLP1-RA: 259 (not included due to small numbers)"
# [1] "Number of subjects already on SGLT2i initiating a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: 39888"
# [1] "Number of drug episodes of already on SGLT2i initiating a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: 49874"
# [1] "Number of drug episodes excluded with unknown eGFR or uACR: 10715"
# [1] "Number of drug episodes excluded with established eGFR <20 mL/min/1.73m2 or ESKD: 92"
# [1] "Number of drug episodes removed (e.g. subsequent episode of starting DPP4i/SU after already taking the other): 427"
# [1] "Number of subjects included: 31650"
# [1] "Number of drug episodes included: 38642"
# this number includes "duplicate" drug episodes - these will be removed later on

rm(t2d_1stinstance)
rm(t2d_all_drug_periods)
gc()

## B Make variables for survival analysis of all endpoints (see survival_variables_kf function for details)

source("survival_variables_kf.R")

cohort <- add_surv_vars(cohort, main_only=FALSE) # add per-protocol survival variables as well

rm(list=setdiff(ls(), c("cohort", "today", "vars", "factor_vars", "nonnormal", "main", "outcomes", "n.imp")))

## set reference group
# list of studydrug vars
studydrug_vars <- grep("^studydrug", names(cohort), value = TRUE)

# Loop through studydrug vars and relevel
for (i in seq_along(studydrug_vars)) {
  ref_level <- if (i == 1) "SGLT2i + SU" else "SGLT2i + DPP4i/SU"
  cohort[[studydrug_vars[i]]] <- relevel(as.factor(cohort[[studydrug_vars[i]]]), ref = ref_level)
}

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


## prepare dataset for imputation: create other relevant variables and deselect redundant variables 
cohort <- cohort %>% 
  mutate(
         # coalesce acr and acr_from_separate into single variable
         uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
         # lowest acr value possible is 0.6
         uacr=ifelse(uacr<0.6, 0.6, uacr),
         # create additional variables that indicate non-glucose-lowering drug prescriptions
         oha=ifelse(ncurrtx > 1 | ncurrtx == 1 & INS == 0, 1L, 0L),
         statin=!is.na(predrug_latest_statins),
         ACE=!is.na(predrug_latest_ace_inhibitors),
         ARB=!is.na(predrug_latest_arb),
         BB=!is.na(predrug_latest_beta_blockers),
         finerenone=!is.na(predrug_latest_finerenone), # all FALSE as finerenone was not in use in the UK until after March 2023
         CCB=!is.na(predrug_latest_ca_channel_blockers),
         ThZD=!is.na(predrug_latest_thiazide_diuretics),
         loopD=!is.na(predrug_latest_loop_diuretics),
         MRA=!is.na(predrug_latest_ksparing_diuretics),
         predrug_cvd=ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1 | predrug_tia==1, 1, 0),
         # default hba1c variable is from previous 6 months to index date
         # take hba1c within window of 2 years prior and 7 days post, similar as other biomarkers
         prehba1c = prehba1c2yrs,
         # create variables reflecting number of eGFR measurements in 12 months prior to and post drug initiation
         preegfr_count_12m = coalesce(preegfr_count_12m, 0),
         egfr_count_12m = coalesce(egfr_count_12m, 0),
         initiation_year = case_when(
           initiation_year %in% c("2013/2014", "2015/2016") ~ "2013/2016",
           TRUE ~ initiation_year),
         # ethnicity cannot be calculated in the imputation model due to it being a constant variable
         # furthermore, the number of individuals in the categories "other" and "mixed" are too small - merge into one miscellaneous category
         ethnicity_qrisk2=ifelse(is.na(ethnicity_qrisk2), "10", as.character(ethnicity_qrisk2)),
         ethnicity_4cat=ifelse(is.na(ethnicity_5cat) | ethnicity_5cat == "4", "3", as.character(ethnicity_5cat)),
         ethnicity_4cat=factor(ethnicity_4cat,
                               levels = c(0, 1, 2, 3),
                               labels = c("White", "South Asian", "Black", "Other or unknown")),
         # imd_decile is not a continuous variable - we will categorise this as quintiles
         imd_decile = ifelse(imd_decile %in% c(1,2), "1/2",
                      ifelse(imd_decile %in% c(3,4), "3/4",
                             ifelse(imd_decile %in% c(5,6), "5/6",
                                    ifelse(imd_decile %in% c(7,8), "7/8",
                                           ifelse(imd_decile %in% c(9,10), "9/10", NA))))),
         imd_decile = factor(imd_decile),
         qrisk2_smoking_cat = factor(qrisk2_smoking_cat),
         ) %>% 
         relocate(ethnicity_4cat, .after = last_col()
         ) %>%
  # select relevant variables only
  select(patid, malesex, ethnicity_4cat, ethnicity_qrisk2, 
         imd_decile, tds_2011,
         regstartdate, 
         gp_end_date, death_date, 
         gp_end_date,
         drug_class, contains("studydrug"), dstartdate, dstopdate_class, drugline_all, drug_substance, ncurrtx,
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
         predrug_retinopathy, predrug_acutepancreatitis, predrug_chronicpancreatitis,
         hosp_admission_prev_year, predrug_efi_score,
         statin, ACE, ARB, BB, finerenone, CCB, 
         ThZD, loopD, MRA, 
         ckd_egfr40_outcome_type, with_hes, preacr_confirmed,
         weightresp6m, weightresp12m, preegfr_count_12m, egfr_count_12m,
         initiation_year,
         ) 


# save dataset
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(cohort, file=paste0(today, "_t2d_glp1_cohort_data.Rda"))
########################2 MULTIPLE IMPUTATION####################################################################
#load(paste0(today, "_t2d_glp1_cohort_data.Rda"))

#dry run
ini <- mice(cohort, seed = 123, maxit = 0)

# remove observations if ncurrtx is missing (needs to be fixed in combo_start_stop)
cohort <- cohort %>% filter(!is.na(ncurrtx))

# specify imputation method for specific variables in ini$meth
# default method is predictive mean matching ("pmm")
meth <- ini$meth

# variables not to be imputed
meth[c(
  "death_date", 
  "preacr", "last_sglt2_stop", "last_glp1_stop", "preckdstage", 
  "dstopdate_class", "timeprevcombo_class", "prehdl", "preldl", "pretriglyceride", "prealt", "predbp",
  "predrug_earliest_ace_inhibitors", 
       "predrug_earliest_arb",
       "predrug_earliest_beta_blockers", "predrug_earliest_ca_channel_blockers",
       "predrug_latest_ace_inhibitors", 
       "predrug_latest_arb",
       "predrug_latest_beta_blockers", "predrug_latest_ca_channel_blockers",
       "ethnicity_qrisk2", 
       "predrug_earliest_thiazide_diuretics", "predrug_latest_thiazide_diuretics",
       "ckd_egfr40_outcome_type", "preacr_confirmed",
       "weightresp6m", "weightresp12m")] <- ""

# set imputation method for smoking status and deprivation as polynomial regression
meth[c("qrisk2_smoking_cat", "imd_decile")] <- "polyreg"

# ensure weight and height are imputed using "pmm"
meth[c("preweight", "height")] <- "pmm"
# BMI should not be imputed directly but rather imputed based on weight and height values
meth["prebmi"] <- "~ I( preweight / (height/100)^2)"

## build prediction matrix - which variables will predict imputed values of which missing variables
# we will use the quickpred function to build the predictor matrix
# we can specify which variables to definitely include (inlist) and which ones to leave out (outlist)

inlist <- c("malesex",  "dstartdate_age",  "imd_decile",  "tds_2011",            # main sociodemographic factors
            paste0("studydrug", main),                                           # treatment variable
            "dstartdate_dm_dur_all", "prebmi", "pretotalcholesterol",            # laboratory and vital sign measurements
            "presbp", "preegfr", "uacr", 
            "qrisk2_smoking_cat", 
            "ckd_egfr40_censvar", "with_hes",
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
              "ckd_egfr40_outcome_type", "preacr_confirmed", 
              "weightresp6m", "weightresp12m", "preegfr_count_12m", "egfr_count_12m",
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


## make final prediction matrix
pred <- quickpred(cohort, include = inlist, exclude = outlist)


# limit imputations to plausible range
post <- ini$post
post["preegfr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(60, 120))"
post["prebmi"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(20, 40))"
post["prehba1c"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(42, 97))"
post["uacr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0.6, 56.5))"
post["presbp"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(80, 180))"


# save % missing and imputation methods per variable imputed:
vars_to_impute <- rownames(pred)[apply(pred, 1, function(x) any(x != 0))]
missing_summary <- data.frame(
  variable = vars_to_impute,
  n_missing = sapply(cohort[vars_to_impute], function(x) sum(is.na(x))),
  n_total = nrow(cohort)
) %>%
  mutate(
    pct_missing = round(100 * n_missing / n_total, 1)
  ) %>%
  arrange(desc(pct_missing))
missing_summary$method <- meth[missing_summary$variable]

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
write.csv2(missing_summary, paste0(today, "_missing_data_summary.csv"), row.names = F)


## impute data
imp <- mice(data = cohort, 
            meth = meth, 
            pred = pred, 
            post = post, 
            m=n.imp, 
            seed = 123)

#check distribution imputed vs original values

density_plot <- densityplot(x = imp, data = ~ imd_decile + prebmi + presbp + pretotalcholesterol +
                              prehba1c + dstartdate_dm_dur_all + qrisk2_smoking_cat + hosp_admission_prev_year)


setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
tiff(paste0(today, "_mice_density_plot.tiff"), width=10, height=6, units = "in", res=800)
print(density_plot)
dev.off()

#extract imputations so we can add variables
temp <- complete(imp, action = "long", include = T)

#the post-processing of prebmi imputations (squeeze) did not compute as expected
#manually force this squeeze
z <- cohort[is.na(cohort$prebmi),]$patid
temp <- temp %>% mutate(prebmi = ifelse(
  patid %in% z & prebmi < 20, 20, prebmi
))
rm(z)


## add CKD-PC risk score
# add variables for risk score calculation
temp <- temp %>% 
  # remove un-imputed dataset
  filter(!.imp == 0) %>% 
  mutate(black_ethnicity=ifelse(!is.na(ethnicity_qrisk2) & (ethnicity_qrisk2 == 6 | ethnicity_qrisk2 == 7), 
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



# tidy variables
temp <- temp %>% mutate(
  # simplify smoking status for covariates
  smoking_status = ifelse(qrisk2_smoking_cat == 0, "never", ifelse(qrisk2_smoking_cat == 1, "ex", "current")),
  # coalesce ACE or ARB use into variable reflecting RAASi use
  ACE_or_ARB = ifelse(temp$ACE + temp$ARB > 0, T, F),
  # categorise number of current glucose-lowering treatments
  ncurrtx2 = ncurrtx,
  ncurrtx = ifelse(ncurrtx==1, "1.", ifelse(ncurrtx==2, "2.", ifelse(ncurrtx == 3, "3.", "4+"))),
  ncurrtx = relevel(as.factor(ncurrtx), ref = "3."),
  # categorise eFI for table
  predrug_efi_cat = case_when(
    predrug_efi_score < 0.12 ~ "fit",
    predrug_efi_score >= 0.12 & predrug_efi_score < 0.24 ~ "mild",
    predrug_efi_score >= 0.24 & predrug_efi_score < 0.36 ~ "moderate",
    predrug_efi_score >= 0.36 ~ "severe"
  ),
  # variable to represent prior history of pancreatitis
  predrug_pancreatitis = ifelse(predrug_acutepancreatitis == T | predrug_chronicpancreatitis == T, T, F),
  # weight change variables
  weight_pct_change_6m = weightresp6m / preweight * 100,
  weight_pct_change_12m = weightresp12m / preweight * 100,
  weight_pct_change = ifelse(is.na(weight_pct_change_12m), weight_pct_change_6m, weight_pct_change_12m),
  across(starts_with("studydrug"), as.factor),
  # categories of eGFR and albuminuria
  egfr_cat = ifelse(preegfr < 60, "<60", "≥60"),
  egfr_cat = factor(egfr_cat),
  albuminuria_cat = ifelse(uacr >=30, "≥30", ifelse(uacr > 3, "3-30", "<3")),
  albuminuria_cat = factor(albuminuria_cat),
  albuminuria_cat2 = ifelse(uacr >=3, "≥3", "<3"),
  albuminuria_cat2 = factor(albuminuria_cat2),
  # categorise eGFR counts for table
  preegfr_count_12m_cat = ifelse(preegfr_count_12m >= 4, "≥4", as.character(preegfr_count_12m)),
  egfr_count_12m_cat = ifelse(egfr_count_12m >= 4, "≥4", as.character(egfr_count_12m)))


q <- temp %>% nrow()
p <- temp %>%  ## in case data contains separate drug episodes if a subject started a DPP4i and later an SU
  group_by(.imp, patid) %>% filter(!duplicated(studydrug2)) %>% ungroup() %>% nrow()
print(paste0("Number of duplicate drug episodes removed ", (q-p)/n.imp))

print(paste0("Number of drug episodes in study population ", p/n.imp))

rm(p)
q <- temp %>% .$patid %>% unique() %>% length()
# print(paste0("Number of subjects in study population ", q))

# remove duplicate episodes
studydrug_var = paste0("studydrug", main)
temp <- temp  %>%
  group_by(.imp, patid, !!sym(studydrug_var)) %>% 
  arrange(dstartdate) %>% 
  filter(!duplicated(!!sym(studydrug_var))) %>% 
  ungroup()

# save imputed dataset so this can be used in the subsequent scripts
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Processed data/")
save(temp, file=paste0(today, "_t2d_glp1_imputed_data.Rda"))


# create table one: this will be an average of the imputed datasets (n to be divided by n.imp)
n.studydrug.vars <- sum(grepl("studydrug", colnames(temp)))

for (m in 1:n.studydrug.vars) {
  
  studydrug_var = paste0("studydrug", m)
  
  table <- CreateTableOne(vars = vars, strata = studydrug_var, data = temp  %>% 
                            group_by(.imp, patid) %>% filter(!duplicated(!!sym(studydrug_var))) %>% ungroup(),  
                          factorVars = factor_vars, test = F)
  
  tabforprint <- print(table, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, printToggle = T)
  
  # save
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
  write.csv2(tabforprint, file = paste0(today, "_baseline_table_studydrug", m, ".csv"))
  
}

## get summary of weight change at 6 / 12 months for studydrug2 variable
temp %>% filter(studydrug2 == "SGLT2i + GLP1-RA") %>% .$weight_pct_change %>% summary()
temp %>% filter(studydrug2 != "SGLT2i + GLP1-RA") %>% .$weight_pct_change %>% summary()


# events rates (sum of events divided by sum of person-years) by studydrug

# prep empty data frame
event_rates <- data.frame(
  studydrug_var = character(),
  drug_level = character(),
  outcome = character(),
  n_at_risk = numeric(),
  events = numeric(),
  pyears = numeric(),
  event_rate_1000yrs = numeric(),
  stringsAsFactors = F
)

temp_all <- temp

for (m in 1:n.studydrug.vars) {  
  
  studydrug_var <- paste0("studydrug", m)
  print(paste0("Event rates for ", studydrug_var))
  
  
  for (k in outcomes) {
    
    print(k)
    
    temp <- temp_all 
    
    if (k == "retinopathy") {
      temp <- temp_all %>% filter(predrug_retinopathy == F)
    }
    
    if (k == "acutepancreatitis") {
      temp <- temp_all %>% filter(predrug_pancreatitis == F)
    }
    
    if (m == 1) {
      # use different censoring variables for studydrug1 (censoring if starting DPP4i or SU if not in that respective arm)
      k <- paste0(k, "_sens1")
    }
 
    censvar_var  <- paste0(k, "_censvar")
    censtime_var <- paste0(k, "_censtime_yrs")  
    
    # iterate over each treatment group
    for (n in levels(temp[[studydrug_var]])) {
      
      temp_sub <- temp %>% filter(.imp != 0 & !!sym(studydrug_var) == n)
      
      n_at_risk <- nrow(temp_sub)
      events  <- sum(temp_sub[[censvar_var]], na.rm = T)
      pyears  <- sum(temp_sub[[censtime_var]], na.rm = T)
      
      event_1000yrs <- round(events / pyears * 1000, 1)
      
      # store result
      event_rates <- rbind(
        event_rates,
        data.frame(
          studydrug_var = studydrug_var,
          drug_level = n,
          outcome = k,
          n_at_risk = n_at_risk/n.imp,
          events = events/n.imp,
          pyears = pyears/n.imp,
          event_rate_1000yrs = event_1000yrs,
          stringsAsFactors = F
        )
      )
    }
  }
}

# save event rates 
setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Output/")
write.csv2(event_rates, paste0(today, "_event_rates_table.csv"), row.names = FALSE)

# give split of event types for ckd_egfr40
temp %>% filter(.imp == n.imp) %>% 
  filter(ckd_egfr40_censvar == 1) %>% 
  .$ckd_egfr40_outcome_type %>% 
  as.factor() %>%
  summary()

## No of confirmed vs 40% at end
temp %>% filter(.imp == n.imp) %>% 
  filter(ckd_egfr40_censvar == 1 & ckd_egfr40_outcome_type == "40% eGFR decline") %>% 
  .$post_egfr_40_decline_date_next_egfr %>% 
  is.na() %>% 
  summary()
