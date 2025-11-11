
# Inclusion/exclusion criteria:
## a) Subjects with T2D
## b) 1st instance of initiating glucose-lowering drug
## c) Exclude if initiating drug within 91 days of registration

## d) Aged 18+
## e) Already receiving SGLT2i initiating GLP1-RA or initiating DPP4i/SU
## f) Initiated between 31/03/2013 and end of HES data (31/03/2023)
## g) No missing eGFR/uACR
## h) No eGFR < 20 mL/min/1.73m2 or ESKD before index date
## j) Remove further episodes of initiating DPP4i/SU if already taking SGLT2i or GLP1-RA in previous episode


# Use "t2d_1stinstance" cohort_dataset which already has a)-d) applied
# all_drug_periods_dataset is used to define later start dates of meds for censoring


define_cohort <- function(cohort_dataset, all_drug_periods_dataset) {
  
  # d)-f) Keep those aged >=18 and within study period
  cohort <- cohort_dataset %>%
    filter(dstartdate_age>=18)
  
  q <- cohort %>% filter(drug_class=="GIPGLP1") %>% nrow()
  print(paste0("Number of subjects on GIPGLP1: ", q, " (not included due to small numbers)"))
  
  cohort <- cohort %>% filter(
    (  drug_class=="GLP1" | 
       drug_class=="DPP4" | 
       drug_class=="SU") &
      dstartdate>=as.Date("2013-03-31") &
      dstartdate<=as.Date("2023-03-31")
  ) %>%
    mutate(studydrug1=ifelse(drug_class=="GLP1", "GLP1-RA", ifelse(drug_class=="DPP4", "DPP4i", "SU")))
  
  # Remove new drug episodes of DPP4i/SU if already initiated GLP1-RA
  cohort <- cohort %>%
    filter(
      !((drug_class=="DPP4" | drug_class=="SU") & (GLP1==1 | GLP1 == 1 & SGLT2==1)))
  
  
  
  ### only keep drug episodes if individuals were already receiving SGLT2i before initiating DPP4i or SU or GLP1-RA
  cohort <- cohort %>% filter(
    ((drug_class=="DPP4" | drug_class=="SU") & SGLT2 == 1) | 
      (drug_class == "GLP1" & SGLT2 == 1)
  )
  
  
  # create variables that denote whether this is first/only/last episode
  cohort <- cohort %>%
    group_by(patid) %>%
    mutate(earliest_start = min(dstartdate, na.rm=TRUE),
           latest_start   = max(dstartdate, na.rm=TRUE),
           episode_order  = ifelse(earliest_start == latest_start, "only", 
                                   ifelse(earliest_start == dstartdate, "first", 
                                          ifelse(latest_start == dstartdate, "last", "other")))) 
  
  # remove drug episodes that are only a change in drug substance, with the exception of initiating a GLP1-RA
  cohort <- cohort %>% 
    filter(
      (drug_class_start == 0 & grepl("tide", drug_substance, ignore.case = T)) | drug_class_start == 1
    ) %>% mutate(
      # create variables to define study drug: 
      ## GLP1-RA vs DPP4i vs SU (sensitivity analysis)
      studydrug1 = ifelse(studydrug1 == "GLP1-RA",
                          "SGLT2i + GLP1-RA", ifelse(studydrug1 == "DPP4i", "SGLT2i + DPP4i", "SGLT2i + SU")),
      ## GLP1-RA vs DPP4i/SU (primary analysis)
      studydrug2 = ifelse(studydrug1 != "SGLT2i + GLP1-RA", "SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA"),
      ## sc semaglutide/dulaglutide vs other GLP1-RA vs DPP4i/SU (sensitivity analysis)
      studydrug3 = ifelse(studydrug2 == "SGLT2i + GLP1-RA",
                          ifelse(grepl("semaglutide", drug_substance, ignore.case = T) & !grepl("oral", drug_substance, ignore.case = T) |    # sc semaglutide
                                   grepl("dulaglutide", drug_substance, ignore.case = T) |                                                    # dulaglutide
                                   grepl("efpeglenatide", drug_substance, ignore.case = T),                                                   # efpeglenatide - not used in the UK
                                 "GLP1-RA with direct kidney outcome evidence", "Other GLP1-RA"), studydrug2)

    )
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects already receiving SGLT2i initiating a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes of already receiving SGLT2i initiating a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: ", q))
  
  
  
  # g) Remove if missing eGFR/uACR
  
  # create uacr variable that calculates acr manually if acr value not available
  cohort <- cohort %>% 
    mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
           uacr=ifelse(uacr<0.6, 0.6, uacr))
  
  q <- cohort %>% filter(is.na(preckdstage) | is.na(preegfr) | is.na(uacr)) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with unknown eGFR/uACR: ", q))
  
  cohort <- cohort %>%
    filter(
      !(is.na(preckdstage) & is.na(preegfr) | is.na(uacr))
    )
  
  # h) Remove if ESKD before index date or eGFR <20
  
  q <- cohort %>% filter(preckdstage=="stage_5" | predrug_ckd5_code == 1 | preegfr < 20) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established eGFR <20 mL/min per 1.73m2 or ESKD: ", q))
  
  cohort <- cohort %>%
    filter(!(preegfr < 20 | preckdstage=="stage_5" | predrug_ckd5_code == 1) )
  
  
  # j) Remove further episodes of initiating DPP4i/SU if already initiated GLP1-RA in previous episode (these episodes would overlap)
  #    or episodes of taking DPP4i following episode of SU and vice versa
  q <- cohort %>% filter(
    episode_order %in% c("last", "other") & (
        (drug_class == "DPP4" & SU == 1) |
        (drug_class == "SU" & DPP4 == 1))
  ) %>% nrow()
  
  print(paste0("Number of drug episodes removed (e.g. subsequent episode of initiating DPP4i/SU after already taking the other): ", q))
  
  cohort <- cohort %>%
    filter(!(
      episode_order %in% c("last", "other") & (
          (drug_class == "DPP4" & SU == 1) |
          (drug_class == "SU" & DPP4 == 1))
    )
    )
  
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects included: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes included: ", q))
  
  rm(q)
  
  ### get post-initiation dates of DPP4i, SU, and GLP1-RA for censoring variables
  
  later_dpp4 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="DPP4") %>%
                  select(patid, next_dpp4=dstartdate)), by="patid") %>%
    filter(next_dpp4>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_dpp4_start=min(next_dpp4, na.rm=TRUE)) %>%
    ungroup()
  
  later_su <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="SU") %>%
                  select(patid, next_su=dstartdate)), by="patid") %>%
    filter(next_su>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_su_start=min(next_su, na.rm=TRUE)) %>%
    ungroup()
  
  later_glp1 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="GLP1") %>%
                  select(patid, next_glp1=dstartdate)), by="patid") %>%
    filter(next_glp1>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_glp1_start=min(next_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  last_dpp4_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="DPP4") %>%
                  select(patid, last_dpp4=dstopdate_class)), by="patid") %>%
    filter(last_dpp4<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_dpp4_stop=min(last_dpp4, na.rm=TRUE)) %>%
    ungroup()
  
  last_su_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="SU") %>%
                  select(patid, last_su=dstopdate_class)), by="patid") %>%
    filter(last_su<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_su_stop=min(last_su, na.rm=TRUE)) %>%
    ungroup()
  
  last_glp1_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="GLP1") %>%
                  select(patid, last_glp1=dstopdate_class)), by="patid") %>%
    filter(last_glp1<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_glp1_stop=min(last_glp1, na.rm=TRUE)) %>%
    ungroup()
  
  
  cohort <- cohort %>%
    left_join(later_dpp4, by=c("patid", "dstartdate")) %>%
    left_join(later_su, by=c("patid", "dstartdate")) %>%
    left_join(later_glp1, by=c("patid", "dstartdate")) %>%
    left_join(last_dpp4_stop, by=c("patid", "dstartdate")) %>%
    left_join(last_su_stop, by=c("patid", "dstartdate")) %>%
    left_join(last_glp1_stop, by=c("patid", "dstartdate")) 
  
  
  # Tidy up gender, ncurrtx, drugline and ethnicity variables
  ## Also code up death cause variables
  
  cohort <- cohort %>%
    
    mutate(malesex=ifelse(gender==1, 1, 0),
           
           ncurrtx=DPP4+GLP1+MFN+SU+SGLT2+TZD+INS,          #ignore Acarbose and Glinide
           
           drugline_all=as.factor(ifelse(drugline_all>=5, 5, drugline_all)),
           
           drug_substance=ifelse(grepl("&", drug_substance), NA, drug_substance),
           
           ethnicity_5cat_decoded=case_when(ethnicity_5cat==0 ~"White",
                                            ethnicity_5cat==1 ~"South Asian",
                                            ethnicity_5cat==2 ~"Black",
                                            ethnicity_5cat==3 ~"Other",
                                            ethnicity_5cat==4 ~"Mixed") ,
           
           cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           kf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(kf_death_any_cause) & kf_death_any_cause==1, death_date, as.Date(NA)),
           kf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(kf_death_primary_cause) & kf_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA))
    )
  
  return(cohort)
  
}