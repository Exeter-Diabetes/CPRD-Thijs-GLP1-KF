
# Produce survival variables for all endpoints (including for sensitivity analysis)
## All censored at 5 years post drug start (3 years for 'ckd_egfr40') / end of GP records / death / starting a different diabetes med which affects CV risk (GLP1 + SGLT2), and also drug stop date + 6 months for per-protocol analysis

# Main analysis:
## 'ckd_egfr50': decline in eGFR of >=50% from baseline or onset of CKD stage 5 OR death from renal causes

# Sensitivity analysis:
## '{outcome}_sens': all of main analysis but with different groupings of drugs
## intention to treat: censoring if starting an SGLT2 inhibitor or GLP1 agonist (if in DPP4 + SU arm).
## sensitivity: censoring if starting any other treatment arm in analyses to compare DPP4i with SU arm or different GLP1s.


add_surv_vars <- function(cohort_dataset, main_only=FALSE) {
  
  # Add survival variables for outcomes for main analysis
  # main_outcomes <- c("ckd_egfr40", "ckd_egfr50", "death", "macroalb", "dka", "amputation", "side_effect")
  # keep ckd_egfr50 as only outcome for now
  main_outcomes <- c("ckd_egfr40", "ckd_egfr50")
  
  cohort <- cohort_dataset %>%
    
    mutate(cens_itt_5_yrs=pmin(dstartdate+(365.25*5),
                         
                         gp_end_date,
                         # death_date,
                         if_else(studydrug2!="GLP1" & studydrug2!="GLP1 + SGLT2", next_glp1_start, as.Date("2050-01-01")),
                         if_else(studydrug2!="SGLT2" & studydrug2!="GLP1 + SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                         na.rm=TRUE),

           
           cens_itt_3_yrs=pmin(dstartdate+(365.25*3),
                               
                               gp_end_date,
                               # death_date,
                               if_else(studydrug2!="GLP1" & studydrug2!="GLP1 + SGLT2", next_glp1_start, as.Date("2050-01-01")),
                               if_else(studydrug2!="SGLT2" & studydrug2!="GLP1 + SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                               na.rm=TRUE),
           
           cens_sens1_3_yrs=pmin(dstartdate+(365.25*3),
                                
                                gp_end_date,
                                # death_date,
                                if_else(studydrug1!="GLP1"  & studydrug1!="GLP1 + SGLT2", next_glp1_start, as.Date("2050-01-01")),
                                if_else(studydrug1!="SGLT2" & studydrug1!="GLP1 + SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                                # censor subjects if starting Su or DPP4 only if in those groups
                                if_else(studydrug1=="DPP4", next_su_start, as.Date("2050-01-01")),
                                if_else(studydrug1=="SU", next_dpp4_start, as.Date("2050-01-01")),
                                #    dstopdate_class+183,
                                na.rm=TRUE),
           
           cens_sens3_3_yrs=pmin(dstartdate+(365.25*3),
                                
                                gp_end_date,
                                # death_date,           
                                if_else(studydrug3!="Other GLP1", next_other_glp1_start, as.Date("2050-01-01")),
                                if_else(studydrug3!="Subcutaneous semaglutide", next_semaglutide_subcut_start, as.Date("2050-01-01")),
                                if_else(studydrug3!="Oral semaglutide", next_semaglutide_oral_start, as.Date("2050-01-01")),
                                if_else(studydrug3!="SGLT2", next_sglt2_start, as.Date("2050-01-01")),
                                #    dstopdate_class+183,
                                na.rm=TRUE),
           

           
           
           
           ckd_egfr40_outcome=pmin(egfr_40_decline_date,
                                   postckdstage5date,
                            #       kf_death_date_any_cause,
                                   na.rm=TRUE),
           
           ckd_egfr50_outcome=pmin(egfr_50_decline_date,
                                   postckdstage5date,
                            #       kf_death_date_any_cause,
                                   na.rm=TRUE),
           
           ckd_egfr50_outcome_type=ifelse(
             is.na(ckd_egfr50_outcome), 
             NA,
             ifelse(ckd_egfr50_outcome == egfr_50_decline_date, 
                    "50% eGFR decline",
                    ifelse(ckd_egfr50_outcome == postckdstage5date, 
                           "ESKD", 
                           "Kidney-related death")
             )
           ),
           
           # macroalb_outcome=macroalb_date,
           # 
           # dka_outcome=postdrug_first_dka,
           # 
           # amputation_outcome=postdrug_first_amputation,
           # 
           # side_effect_outcome=pmin(postdrug_first_medspecific_gi,
           #                          na.rm=TRUE),
           # 
           # death_outcome=death_date
           )
  
  
  for (i in main_outcomes) {
    
    # outcome variable at any time
    outcome_var=paste0(i, "_outcome")
    
    # outcome/censor and censoring time variables (within 3 years)
    censdate_var=paste0(i, "_censdate")
    censvar_var=paste0(i, "_censvar")
    censtime_var=paste0(i, "_censtime_yrs")
    
    # outcome/censor and censoring time variables (within 5 years)
    censdate_var_5y=paste0(i, "_5y_censdate")
    censvar_var_5y=paste0(i, "_5y_censvar")
    censtime_var_5y=paste0(i, "_5y_censtime_yrs")
    
    
    
    cohort <- cohort %>%
        mutate(
          {{censdate_var}}:=pmin(!!sym(outcome_var), cens_itt_3_yrs, na.rm=TRUE),
          {{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
          {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25,
          
          {{censdate_var_5y}}:=pmin(!!sym(outcome_var), cens_itt_5_yrs, na.rm=TRUE),
          {{censvar_var_5y}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var_5y)==!!sym(outcome_var), 1, 0),
          {{censtime_var_5y}}:=as.numeric(difftime(!!sym(censdate_var_5y), dstartdate, unit="days"))/365.25,
          
           )

    
  }
  
  if (main_only==TRUE) {
    message(paste("survival variables for", paste(main_outcomes, collapse=", "), "added"))
  }
  
  
  # Add survival variables for outcomes for sensitivity analyses
  
  else {
    
    # sensitivity analyses
    sensitivity_outcomes <- c("ckd_egfr50_sens1", "ckd_egfr50_sens3")
    
    
    for (i in sensitivity_outcomes) {
      
      sens = substr(i, nchar(i)-4, nchar(i))
      censoring_var_im=paste0("cens_", sens, "_3_yrs")
      censdate_var=paste0(i, "_censdate")
      censvar_var=paste0(i, "_censvar")
      censtime_var=paste0(i, "_censtime_yrs")
      
      
      outcome_var=paste0(substr(i, 1,  nchar(i)-6), "_outcome")
      

        cohort <- cohort %>%
          mutate({{censdate_var}}:=pmin(!!sym(outcome_var), !!sym(censoring_var_im), na.rm=TRUE))
      
      
      
      cohort <- cohort %>%
        mutate({{censvar_var}}:=ifelse(!is.na(!!sym(outcome_var)) & !!sym(censdate_var)==!!sym(outcome_var), 1, 0),
               {{censtime_var}}:=as.numeric(difftime(!!sym(censdate_var), dstartdate, unit="days"))/365.25)
      
    }
    
    if (main_only==FALSE) {
      message(paste("survival variables for", paste(main_outcomes, collapse=", "), ",", paste(unlist(sensitivity_outcomes), collapse=", "), "added"))
    }
    
  }
  
  return(cohort) 
  
}