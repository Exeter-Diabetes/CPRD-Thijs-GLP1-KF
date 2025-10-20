############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

analysis = cprd$analysis("thijs_glp1")


############################################################################################

# Today's date for table names

today <- as.character(Sys.Date(), format="%Y%m%d")


############################################################################################

############################################################################################

## Filter just type 2s
t2d_1stinstance <- t2d_1stinstance %>% 
  analysis$cached("20251019_t2d_1stinstance", indexes=c("patid", "dstartdate", "drug_class", "drug_substance"))

### Check unique patid count
# t2d_1stinstance %>% distinct(patid) %>% count()
#1,025,666

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")

t2d_1stinstance_a <- collect(t2d_1stinstance %>% filter(patid<2000000000000) %>% mutate(patid=as.character(patid)))

is.integer64 <- function(x){
  class(x)=="integer64"
}

t2d_1stinstance_a <- t2d_1stinstance_a %>%
  mutate_if(is.integer64, as.integer)

save(t2d_1stinstance_a, file=paste0(today, "_t2d_1stinstance_a.Rda"))

rm(t2d_1stinstance_a)


t2d_1stinstance_b <- collect(t2d_1stinstance %>% filter(patid>=2000000000000) %>% mutate(patid=as.character(patid)))

t2d_1stinstance_b <- t2d_1stinstance_b %>%
  mutate_if(is.integer64, as.integer)

save(t2d_1stinstance_b, file=paste0(today, "_t2d_1stinstance_b.Rda"))

rm(t2d_1stinstance_b)



############################################################################################

# Make dataset of all drug starts so that can see whether people later initiate SGLT2i/GLP1 etc.
## Add in discontinuation variables

## Just T2s
t2d_all_drug_periods <- t2d_all_drug_periods %>%
  analysis$cached("20251019_t2d_all_drug_periods")


t2d_all_drug_periods <- collect(t2d_all_drug_periods %>% mutate(patid=as.character(patid)))

save(t2d_all_drug_periods, file=paste0(today, "_t2d_all_drug_periods.Rda"))