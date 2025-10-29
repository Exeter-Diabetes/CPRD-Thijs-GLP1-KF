
# Produce survival variables for all endpoints (including for sensitivity analysis)
## All censored at 3 years post drug initiation (5 years for sensitivity analysis) / end of GP records / death / starting a GLP1-RA if not in the treatment group, and also drug stop date + 6 months for per-protocol analysis

# Main analysis:
## 'ckd_egfr40': decline in eGFR of >=40% from baseline OR onset of CKD stage 5 OR death from kidney-related causes

# Sensitivity analysis:
## '{outcome}_sens': all of main analysis but with different groupings of drugs
## intention to treat: censoring if starting an SGLT2i or GLP1-RA (if in DPP4i/SU arm).


add_surv_vars <- function(cohort_dataset, main_only=FALSE) {
  
  # Add survival variables for outcomes for main analysis
  
  
  main_outcomes <- c("ckd_egfr40", "death", "mace", "hf", "lowerlimbfracture", "retinopathy", "acutepancreatitis")
  
  
  cohort <- cohort_dataset %>%
    
    mutate(gp_end_date=pmax(gp_end_date, as.Date("2023-03-31"), na.rm=TRUE),
           
           
           cens_itt_3_yrs=pmin(dstartdate+(365.25*3),
                               hes_end_date,
                               gp_end_date,
                               death_date,
                               if_else(studydrug2!="SGLT2 + GLP1", next_glp1_start, as.Date("2050-01-01")),
                               na.rm=TRUE),
           
           
           cens_itt_5_yrs=pmin(dstartdate+(365.25*5),
                               hes_end_date,
                               gp_end_date,
                               death_date,
                               if_else(studydrug2!="SGLT2 + GLP1", next_glp1_start, as.Date("2050-01-01")),
                               na.rm=TRUE),
           
           cens_sens1_3_yrs=pmin(dstartdate+(365.25*3),
                                 hes_end_date,
                                 gp_end_date,
                                 death_date,
                                 if_else(studydrug1!="SGLT2 + GLP1", next_glp1_start, as.Date("2050-01-01")),
                                 # censor subjects if starting Su or DPP4 only if in those groups
                                 if_else(studydrug1=="DPP4", next_su_start, as.Date("2050-01-01")),
                                 if_else(studydrug1=="SU", next_dpp4_start, as.Date("2050-01-01")),
                                 na.rm=TRUE),
           
           cens_pp_3_yrs=pmin(dstartdate+(365.25*3),
                                 hes_end_date,
                                 gp_end_date,
                                 death_date,
                                 if_else(studydrug2!="SGLT2 + GLP1", next_glp1_start, as.Date("2050-01-01")),
                                 dstopdate_class+183,
                                 na.rm=TRUE),
           
           
           
           # if eGFR measurement showing a 40% decline is a final measurement, use this
           egfr40_decline_date = ifelse(
             is.na(cohort$post_egfr_40_decline_date_next_egfr) | is.null(cohort$post_egfr_40_decline_date_next_egfr),
             cohort$egfr_40_decline_date,
             cohort$egfr_40_decline_date_confirmed
           ),
           
           ckd_egfr40_outcome=pmin(egfr40_decline_date,
                                   postckdstage5date,
                                   kf_death_date_any_cause,
                                   na.rm=TRUE),
           

           ckd_egfr40_outcome_type=ifelse(
             is.na(ckd_egfr40_outcome), 
             NA,
             ifelse(ckd_egfr40_outcome == egfr_40_decline_date, 
                    "40% eGFR decline",
                    ifelse(ckd_egfr40_outcome == postckdstage5date, 
                           "ESKD", 
                           "Kidney-related death")
             )
           ),
           
           mace_outcome=pmin(postdrug_first_myocardialinfarction,
                             postdrug_first_stroke,
                             cv_death_date_any_cause,
                             na.rm=TRUE),
           
           hf_outcome=pmin(postdrug_first_heartfailure,
                           hf_death_date_any_cause,
                           na.rm=TRUE),
           
           acutepancreatitis_outcome=pmin(postdrug_first_acutepancreatitis,
                                          na.rm=TRUE),
           
           retinopathy_outcome=pmin(postdrug_first_retinopathy,
                                    na.rm=TRUE),
           
           lowerlimbfracture_outcome=pmin(postdrug_first_lowerlimbfracture,
                                          na.rm=TRUE),
           
           death_outcome=death_date
           
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
    
    for (q in c("sens1", "pp")) {
      # sensitivity analyses
      sensitivity_outcomes <- paste0(main_outcomes, "_", q)
      
      
      for (i in sensitivity_outcomes) {
        
        censoring_var_im=paste0("cens_", q, "_3_yrs")
        censdate_var=paste0(i, "_censdate")
        censvar_var=paste0(i, "_censvar")
        censtime_var=paste0(i, "_censtime_yrs")
        
        
        outcome_var=paste0(substr(i, 1,  nchar(i)-(nchar(q) + 1)), "_outcome")
        
        
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
  }
  
  return(cohort) 
  
}