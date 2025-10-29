
# Extract dataset of all first instance drug periods (i.e. the first time patient has taken this particular drug class) for ALL DIABETES/T2Ds ONLY WITH HES LINKAGE
## Exclude drug periods starting within 90 days of registration
## Set drugline to missing where diagnosed before registration

## Do not exclude where first line
## Do not exclude where patient is on insulin at drug initiation
## Do not exclude where only 1 prescription (dstartdate=dstopdate)

## Set hosp_admission_prev_year to 0/1 rather than NA/1

# Also extract all drug start and stop dates so that you can see if people later initiate SGLT2is/GLP1s etc.


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

analysis = cprd$analysis("mm")

############################################################################################

# Today's date for table names

# today <- format(Sys.Date(), "%Y%m%d")
today <- "20251019"

############################################################################################

# Get handles to pre-existing data tables

## Cohort and patient characteristics including Townsend scores and death causes
analysis = cprd$analysis("all")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")
townsend_score <- townsend_score %>% analysis$cached("patid_townsend_score")
death_causes <- death_causes %>% analysis$cached("death_causes")

## Drug info
analysis = cprd$analysis("mm")
drug_start_stop <- drug_start_stop %>% analysis$cached("drug_start_stop")
combo_start_stop <- combo_start_stop %>% analysis$cached("combo_start_stop")

## Biomarkers inc. CKD
#baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
analysis = cprd$analysis("thijs_glp1")
response_biomarkers <- response_biomarkers %>% analysis$cached("response_biomarkers") #includes baseline biomarker values for first instance drug periods so no need to use baseline_biomakers table
analysis = cprd$analysis("mm")
ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages")

# number of eGFR counts in 12 months post baseline
egfr_counts_12m <- egfr_counts_12m %>% analysis$cached("response_biomarkers_egfr_count_12m") 
preegfr_counts_12m <- preegfr_counts_12m %>% analysis$cached("response_biomarkers_preegfr_count_12m") 

## Comorbidities and eFI
analysis = cprd$analysis("thijs_glp1")
comorbidities <- comorbidities %>% analysis$cached("comorbidities")

analysis = cprd$analysis("mm")
efi <- efi %>% analysis$cached("efi")

## Non-diabetes meds
non_diabetes_meds <- non_diabetes_meds %>% analysis$cached("non_diabetes_meds")

## Smoking status at drug start
smoking <- smoking %>% analysis$cached("smoking")

#Alcohol at drug start
alcohol <- alcohol %>% analysis$cached("alcohol")

## Discontinuation
discontinuation <- discontinuation %>% analysis$cached("discontinuation")

## Glycaemic failure
glycaemic_failure <- glycaemic_failure %>% analysis$cached("glycaemic_failure")


############################################################################################

# Make first instance drug period dataset

## Define all diabetes cohort (1 line per patient)
## Add in Townsend Deprivation Scores
all_diabetes <- diabetes_cohort %>%
  left_join((townsend_score %>% select(patid, tds_2011)), by="patid") %>%
  relocate(tds_2011, .after=imd_decile)


## Get info for first instance drug periods for cohort (1 line per patid-drug_substance period)
### Make new drugline variable which is missing where diagnosed before registration or within 90 days following

all_diabetes_drug_periods <- all_diabetes %>%
  inner_join(drug_start_stop, by="patid") %>%
  inner_join(combo_start_stop, by=c("patid", c("dstartdate"="dcstartdate"))) %>%
  mutate(drugline=ifelse(dm_diag_date_all<regstartdate | is.na(dm_diag_date), NA, drugline_all)) %>%
  relocate(drugline, .after=drugline_all) %>%
  analysis$cached(paste0(today, "_all_1stinstance_interim_1"), indexes=c("patid", "dstartdate", "drug_class", "drug_substance"))

all_diabetes_drug_periods %>% distinct(patid) %>% count()
# 1,752,846


### Keep first instance only
all_diabetes_1stinstance <- all_diabetes_drug_periods %>%
  filter(drug_instance==1)

all_diabetes_1stinstance %>% distinct(patid) %>% count()
# 1,752,846 as above


### Exclude drug periods starting within 90 days of registration
all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  filter(datediff(dstartdate, regstartdate)>90)

all_diabetes_1stinstance %>% count()
# 3,004,907

all_diabetes_1stinstance %>% distinct(patid) %>% count()
# 1,346,921


### Exclude drug periods starting after or on same day as death
all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  filter(is.na(death_date) | death_date>dstartdate)

all_diabetes_1stinstance %>% count()
# 3,004,813

all_diabetes_1stinstance %>% distinct(patid) %>% count()
# 1,346,887


## Merge in biomarkers, comorbidities, eFI, non-diabetes meds, smoking status, alcohol
### Could merge on druginstance too, but quicker not to
### Remove some variables to avoid duplicates
### Now in two stages to speed it up

analysis = cprd$analysis("thijs_glp1")


all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  inner_join((response_biomarkers %>% select(-c(drug_class, drug_instance, timeprevcombo_class))), by=c("patid", "dstartdate", "drug_substance")) %>%
  inner_join((ckd_stages %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  inner_join((comorbidities %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  inner_join((efi %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  analysis$cached(paste0(today, "_all_1stinstance_interim_2"), indexes=c("patid", "dstartdate", "drug_substance"))

all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  inner_join((non_diabetes_meds %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  inner_join((smoking %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  inner_join((alcohol %>% select(-c(drug_class, drug_instance))), by=c("patid", "dstartdate", "drug_substance")) %>%
  left_join(death_causes, by="patid") %>%
  mutate(hosp_admission_prev_year=ifelse(is.na(hosp_admission_prev_year) & with_hes==1, 0L,
                                         ifelse(hosp_admission_prev_year==1, 1L, NA)),
         hosp_admission_prev_year_count=ifelse(is.na(hosp_admission_prev_year_count) & with_hes==1, 0L, hosp_admission_prev_year_count)) %>%
  analysis$cached(paste0(today, "_all_1stinstance_interim_3"), indexes=c("patid", "dstartdate", "drug_substance"))


## Merge in glycaemic failure and discontinuation - both are currently by drug class (not substance) only
### Glycaemic failure doesn't include all drug periods - only those with HbA1cs
### Make new variables: age at drug start, diabetes duration at drug start

all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  left_join((glycaemic_failure %>% select(-c(dstopdate, timetochange, timetoaddrem, nextdrugchange, nextdcdate, prehba1c, prehba1cdate, threshold_7.5, threshold_8.5, threshold_baseline, threshold_baseline_0.5))), by=c("patid", "dstartdate", "drug_class")) %>%
  left_join((discontinuation %>% select(-c(drugline_all, dstopdate_class, drug_instance, timeondrug, nextremdrug, timetolastpx))), by=c("patid", "dstartdate"="dstartdate_class", "drug_class")) %>%
  mutate(dstartdate_age=datediff(dstartdate, dob)/365.25,
         dstartdate_dm_dur_all=datediff(dstartdate, dm_diag_date_all)/365.25,
         dstartdate_dm_dur=datediff(dstartdate, dm_diag_date)/365.25) %>%
  analysis$cached(paste0(today, "_all_1stinstance_interim_4"), indexes=c("patid", "dstartdate", "drug_substance"))

all_diabetes_1stinstance <- all_diabetes_1stinstance %>%
  left_join((egfr_counts_12m), by = c("patid", "dstartdate", "drug_substance")) %>%
  left_join((preegfr_counts_12m), by = c("patid", "dstartdate", "drug_substance")) %>%
  analysis$cached(paste0(today, "_all_1stinstance_interim_5"), indexes=c("patid", "dstartdate", "drug_substance"))

# Check counts

all_diabetes_1stinstance %>% count()
# 3,004,813

all_diabetes_1stinstance %>% distinct(patid) %>% count()
# 1,346,887

############################################################################################

## Filter just type 2s
t2d_1stinstance <- all_diabetes_1stinstance %>% filter(diabetes_type=="type 2") %>%
  analysis$cached(paste0(today, "_t2d_1stinstance"), indexes=c("patid", "dstartdate", "drug_class", "drug_substance"))

### Check unique patid count
t2d_1stinstance %>% distinct(patid) %>% count()
#1,269,977


############################################################################################

# Make dataset of all drug starts so that can see whether people later initiate SGLT2i/GLP1-RA etc.
## Not cleaned to remove those close to reg start / after death

## Just T2s
t2d_all_drug_periods <- all_diabetes %>%
  filter(diabetes_type=="type 2") %>%
  select(patid) %>%
  inner_join(drug_start_stop, by="patid") %>%
  analysis$cached(paste0(today, "_t2d_all_drug_periods"))

