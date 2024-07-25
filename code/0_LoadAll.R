## Load required packages
require(pacman)
require(tidyverse)
require(httr) 
require(jsonlite)
require(countrycode) 
require(data.table) 
require(socialmixr) 
require(imputeTS)
require(lubridate)
require(mgcv) 
require(DEoptim) 
require(magrittr)
require(progress)
require(readxl)
require(Rcpp)
require(here) 
require(testthat)

##### load covidm #####
data_path <- "data/"
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))

# A. Population structure
fread(paste0(data_path, "pop_str_2021.csv")) %>%
  gather(key = "sex", value = "pop", both) %>%
  mutate(pop = parse_number(pop)) |> 
  mutate(age_group = age %/% 5 + 1,
         # age = paste0(age_group * 5, "-", age_group * 5 + 4),
         # age = replace(age, age_group>=18, "90+"),
         # age = factor(age, levels = limits_to_agegroups(seq(0, 90, by = 5))),
         age_group = if_else(age_group > 16, 16, age_group)) |> 
  group_by(age_group) |> 
  summarise(pop_age = sum(pop)) -> pop_TH

cm_populations |> 
  filter(name == "Thailand") |> 
  separate(age, into = c("age_LL", "age_UL")) |> 
  mutate(age_group = as.numeric(age_LL) %/%5 + 1,
         age_group = if_else(age_group > 16, 16, age_group)) |> 
  group_by(age_group) |> 
  summarise(tot_age = (f + m) * 1000) |> 
  ungroup() |> mutate(tot = sum(tot_age)) -> popTH_cm

# B. Vaccine uptake
source("code/0_4_Vaccinations.R")

# B. load custom functions
source("code/0_0_util.R")



## Load required data
# A. Covid-19 deaths
# B. Covid-19 cases/hospitalisations
# home isolation policy were proclaimed on 4 Jan 2022, thus, new cases afterwards were mixed between hospitalised and home isolated.

# res_round1to2 <- GET("https://covid19.ddc.moph.go.th/api/Cases/round-1to2-all")
# df_round1to2 <- fromJSON(rawToChar(res_round1to2$content)) 
# write.csv(df_round1to2, file = "data/epi_round1to2.csv", row.names = F)
epi_round1to <- fread(paste0(data_path, "epi_round1to2.csv")) # snapshot from 2020-01-12 to 2021-03-31

# download the data file from DDC MOPH TH if it doesn't exist in your directory
# if(!file.exists(paste0("data/epi_update.csv"))){
#   res_update <- GET("https://covid19.ddc.moph.go.th/api/Cases/timeline-cases-all")
#   df_update <- fromJSON(rawToChar(res_update$content))[,1:9] 
#   write.csv(df_update, file = "data/epi_update.csv", row.names = F)
# }
 
# update the data file from DDC MOPH TH if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/epi_update.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   res_update <- GET("https://covid19.ddc.moph.go.th/api/Cases/timeline-cases-all")
#   df_update <- fromJSON(rawToChar(res_update$content))[,1:9] 
#   write.csv(df_update, file = "data/epi_update.csv", row.names = F)
# }

epi_update <- fread(paste0(data_path, "epi_update.csv")) # time-series from 2021-04-01 onward

epi <- bind_rows(epi_round1to, epi_update) %>%
  dplyr::select(txn_date, new_case_excludeabroad, new_death) %>% 
  distinct()
  
rm(epi_round1to, epi_update)

# I. Stringency index
# source("code/0_2_StringencyIndex.R")
# H. Google mobility data
# source("code/0_3_Mobility.R")
contact_schedule <- read_rds(paste0(data_path,"c_schedule.rds"))

#### K. Epi parameters ####
source("code/0_5_EpiParams.R")

#### L. Burden processes #### 
country_list <- read_rds(paste0(data_path, "country_list_thailand.rds"))
HSR_cleaned <- read_rds(paste0(data_path, "HSR_cleaned_thailand.rds"))
source("code/0_6_HealthCareSystem.R")

#### Vaccine Market ####
source("code/0_7_vaccine_market.R")

# birth rate
cbr <- read_rds(paste0(data_path, "cbr.rds"))

# death rate
mu_weighted_16 <- read_rds(paste0(data_path, "mu_weighted_16.rds"))

# fitted table (test)
fitted_table_baseline <- read_rds(paste0(data_path, "fitted_table_baseline_THA.rds"))

# states
compartment_pop <- c("S", "Sv_l", "Sv_m", "Sv_h",
                     "E", "Ev_l", "Ev_m", "Ev_h",
                     "Ip", "Ip_l", "Ip_m", "Ip_h",
                     "Ia", "Ia_l", "Ia_m", "Ia_h",
                     "Is", "Is_l", "Is_m", "Is_h",
                     "R", "Rv_l", "Rv_m", "Rv_h")

compartment_process <- c("cases", "cases_reported",
                         "subclinical",
                         "foi", "foiv_l", "foiv_m", "foiv_h")

compartment_process_voc <- c("severe", "critical", "death")

# other index

voc_phases <- read_rds(paste0(data_path, "voc_phases_thailand.rds"))
voc_phases_imputation_index <-  read_rds(paste0(data_path, "voc_phases_imputation_index_thailand.rds"))
voc_features_test <- read_rds(paste0(data_path, "voc_features_test_thailand.rds"))
load(paste0(data_path, "severe_strain_thailand.rdata"))