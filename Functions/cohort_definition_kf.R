
# Inclusion/exclusion criteria:
## a) Subjects with T2D
## b) With HES linkage
## c) 1st instance
## d) Exclude if start drug within 91 days of registration

## e) Aged 18+
## f) SGLT2 + GLP1 or SGLT2 + DPP4i/SU
## g) Initiated between 01/01/2013 and end of HES data (31/03/2023)
## j) No missing eGFR/uACR
## k) No advanced CKD (egfr < 20 mL/min/1.73m2 or CKD stage 5) before index date
## l) Remove further episodes of starting DPP4i/SU if already taking SGLT2i or GLP1 in previous episode (these episodes would overlap)


# Use "t2d_1stinstance" cohort_dataset which already has a)-d) applied
# all_drug_periods_dataset is used to define later start dates of meds for censoring


define_cohort <- function(cohort_dataset, all_drug_periods_dataset) {
  
  # e)-g) Keep those aged >=18 and within study period
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
  
  # Remove new drug episodes of DPP4i/SU if already on GLP1-RA or GLP1-RA + SGLT2i
  cohort <- cohort %>%
    filter(
      !((drug_class=="DPP4" | drug_class=="SU") & (GLP1==1 | GLP1 == 1 & SGLT2==1)))
  
  
  
  ### only keep drug episodes if individuals were already on SGLT2i before starting DPP4i or SU or GLP1
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
  
  # remove drug episodes that are only a change in drug substance, with the exception of starting a GLP1-RA
  cohort <- cohort %>% 
    filter(
      (drug_class_start == 0 & grepl("tide", drug_substance, ignore.case = T)) | drug_class_start == 1
    ) %>% mutate(
      # create other variables to define study drug: 
      # create combination group for dual therapy with SGLT2i and GLP1-RA (studydrug1)
      ## combine DPP4i/SU as one group (studydrug2) 
      ## create distinct level for oral vs sc semaglutide semaglutide vs other GLP1-RAs (studydrug3)
      studydrug1 = ifelse(studydrug1 == "GLP1-RA",
                          "SGLT2i + GLP1-RA", ifelse(studydrug1 == "DPP4i", "SGLT2i + DPP4i", "SGLT2i + SU")),
      studydrug2 = ifelse(studydrug1 != "SGLT2i + GLP1-RA", "SGLT2i + DPP4i/SU", "SGLT2i + GLP1-RA"),
      studydrug3 = ifelse(studydrug2 == "SGLT2i + GLP1-RA",
                          ifelse(grepl("semaglutide", drug_substance, ignore.case = T) |#& !grepl("oral", drug_substance, ignore.case = T) |    # sc semaglutide
                                   grepl("dulaglutide", drug_substance, ignore.case = T) |                                                    # dulaglutide
                                   grepl("efpeglenatide", drug_substance, ignore.case = T),                                                   # efpeglenatide
                                 "GLP1-RA with direct kidney outcome evidence", "Other GLP1-RA"), studydrug2),
      studydrug4 = ifelse(studydrug2 == "SGLT2i + GLP1-RA",
                          ifelse(grepl("semaglutide", drug_substance, ignore.case = T), 
                                 "Semaglutide", 
                                 ifelse(grepl("dulaglutide", drug_substance, ignore.case = T),
                                        "Dulaglutide", "Other GLP1-RA")), studydrug2)

    )
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects already on SGLT2i starting a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes of already on SGLT2i starting a GLP1-RA or comparator drugs DPP4i/SU between 2013-2023: ", q))
  
  
  
  # j) Remove if missing CKD status
  
  # create acr variable for ckdpc risk scores that uses further source of acr if acr not available
  cohort <- cohort %>% 
    mutate(uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
           uacr=ifelse(uacr<0.6, 0.6, uacr))
  
  q <- cohort %>% filter(is.na(preckdstage) | is.na(preegfr) | is.na(uacr)) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with unknown CKD status: ", q))
  
  cohort <- cohort %>%
    filter(
      !(is.na(preckdstage) & is.na(preegfr) | is.na(uacr))
    )
  
  # k) Remove if ESKD before index date or eGFR <60
  
  q <- cohort %>% filter(preckdstage=="stage_5" | predrug_ckd5_code == 1 | preegfr < 20) %>% nrow()
  
  print(paste0("Number of drug episodes excluded with established eGFR <20 mL/min/1.73m2 or ESKD: ", q))
  
  cohort <- cohort %>%
    filter(!(preegfr < 20 | preckdstage=="stage_5" | predrug_ckd5_code == 1) )
  
  
  # m) Remove further episodes of starting DPP4i/SU if already taking SGLT2i or GLP1-RA in previous episode (these episodes would overlap)
  #    or episodes of taking DPP4i following episode of SU and vice versa
  q <- cohort %>% filter(
    episode_order %in% c("last", "other") & (
      # ((drug_class == "DPP4" | drug_class == "SU") & (SGLT2 == 1 | GLP1-RA == 1)) |
        (drug_class == "DPP4" & SU == 1) |
        (drug_class == "SU" & DPP4 == 1))
  ) %>% nrow()
  
  print(paste0("Number of drug episodes removed (e.g. subsequent episode of starting DPP4i/SU after already taking the other): ", q))
  
  cohort <- cohort %>%
    filter(!(
      episode_order %in% c("last", "other") & (
        # ((drug_class == "DPP4" | drug_class == "SU") & (SGLT2 == 1 | GLP1 == 1)) |
          (drug_class == "DPP4" & SU == 1) |
          (drug_class == "SU" & DPP4 == 1))
    )
    )
  
  
  
  q <- cohort %>% .$patid %>% unique() %>% length()
  print(paste0("Number of subjects included: ", q))
  
  q <- cohort %>% nrow()
  print(paste0("Number of drug episodes included: ", q))
  
  rm(q)
  
  ## Use all SGLT2i, GLP1-RA, DPP4i, and SU starts to code up later censoring
  
  #
  ### Also get latest GLP1-RA and SGLT2i stop dates before drug start for DPP4i/SU arms 
  ### we need this for sensitivity analysis where we exclude people who tried SGLT2i in the year before starting a DPP4i/SU.
  
  later_sglt2 <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="SGLT2") %>%
                  select(patid, next_sglt2=dstartdate)), by="patid") %>%
    filter(next_sglt2>dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE)) %>%
    ungroup()
  
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
  
  last_sglt2_stop <- cohort %>%
    select(patid, dstartdate) %>%
    inner_join((all_drug_periods_dataset %>%
                  filter(drug_class=="SGLT2") %>%
                  select(patid, last_sglt2=dstopdate_class)), by="patid") %>%
    filter(last_sglt2<dstartdate) %>%
    group_by(patid, dstartdate) %>%
    summarise(last_sglt2_stop=min(last_sglt2, na.rm=TRUE)) %>%
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
    left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
    left_join(later_dpp4, by=c("patid", "dstartdate")) %>%
    left_join(later_su, by=c("patid", "dstartdate")) %>%
    left_join(later_glp1, by=c("patid", "dstartdate")) %>%
    left_join(last_sglt2_stop, by=c("patid", "dstartdate")) %>%
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