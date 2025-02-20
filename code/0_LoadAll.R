## Load required packages
if(!require(pacman)) install.packages("pacman")
p_load(tidyverse, httr, jsonlite, countrycode, data.table, socialmixr, imputeTS,
       lubridate, mgcv, DEoptim, magrittr, progress, readxl, Rcpp, here, testthat)

##### load covidm #####
data_path <- "/Users/yangliu/Dropbox/Github_Data/HITAP_CovidM/"
# "C:/Users/eideyliu/Dropbox/Github_Data/HITAP_CovidM/"
#"D:/GitHub/covidmodel-hitap/data/"
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

pop_proj <- read_rds(paste0(data_path, "pop_101.rds")) %>% dplyr::filter(country_code == "THA")
# pop_proj %>% dplyr::filter(year == 2030, age_from >= 75) %>% pull(value) %>% sum


# B. Vaccine uptake
source("code/0_4_Vaccinations.R")

# load HSR Cleaned
HSR <- read_rds(paste0(data_path, "HealthSystemRates_by_country.rds")) %>%
  group_by(country_code) %>% group_split()
HSR_labels <- (HSR %>% map(pull, country_code) %>% map(unique) %>% unlist)
HSR %>%
  map(dplyr::select,
      age_group,
      ihr_by_age_group,
      ifr,
      picu_by_age_group) %>%
  map(rename,
      ihr = ihr_by_age_group,
      picu = picu_by_age_group) %>%
  setNames(HSR_labels) %>%
  map(mutate,
      P.critical = ihr*picu,
      P.severe = ihr*(1-picu),
      P.death = ifr,
      P.hosp = ihr) %>%
  setNames(HSR_labels) -> HSR_cleaned

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

# C. PCR positivity rate
# download the data file from COLAB-2 if it doesn't exist in your directory
# if(!file.exists(paste0("data/thailand_covid-19_testing_data.csv"))){
#   download.file("https://data.go.th/dataset/9f6d900f-f648-451f-8df4-89c676fce1c4/resource/0092046c-db85-4608-b519-ce8af099315e/download/",
#                 paste0("data/thailand_covid-19_testing_data.csv"))
# }
# 
# # update the data file from COLAB-2 if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/thailand_covid-19_testing_data.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   download.file("https://data.go.th/dataset/9f6d900f-f648-451f-8df4-89c676fce1c4/resource/0092046c-db85-4608-b519-ce8af099315e/download/",
#                 paste0("data/thailand_covid-19_testing_data.csv"))
# }

pcr_rate <- fread(paste0(data_path, "thailand_covid-19_testing_data.csv"))[,1:3] %>%
  filter(!is.na(positive)) %>%
  mutate_at(vars(Date), ~lubridate::dmy(.)) %>%
  rename("total_test" = "Total Testing") %>%
  mutate(pos_rate = positive/total_test)

# D. Seroprevalence
# snapshot over the past 2 years, surveillance at Nov 2021
# sero <- fread("data/serosurveillance65.csv")

# E. Genotype frequencies
# download the data file from COLAB-2 if it doesn't exist in your directory
# if(!file.exists(paste0("data/sars-cov-2-variants-dmsc.csv"))){
#   download.file("https://data.go.th/dataset/5b1fb1cf-7ddf-4194-89be-c2658fdcd7a8/resource/152ed762-3c69-465e-a5ae-e592540559d8/download/",
#                 paste0("data/sars-cov-2-variants-dmsc.csv"))
# }
# 
# # update the data file from COLAB-2 if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/sars-cov-2-variants-dmsc.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   download.file("https://data.go.th/dataset/5b1fb1cf-7ddf-4194-89be-c2658fdcd7a8/resource/152ed762-3c69-465e-a5ae-e592540559d8/download/",
#                 paste0("data/sars-cov-2-variants-dmsc.csv"))
# }

# geno_freq <- fread(paste0(data_path, "sars-cov-2-variants-dmsc.csv"))[,1:6] %>%
#   mutate_at(vars(Date_Start, Date_End), ~as.Date(., format = "%d/%m/%Y")) %>%
#   rename("Alpha" = "B.1.1.7 (Alpha)",
#          "Delta" = "B.1617.2 (Delta)",
#          "Beta" = "B.1.351 (Beta)",
#          "Omicron" = "B.1.1.529 (Omicron") %>%
#   filter(!is.na(Alpha))

# G. Contact matrices
# load("data/contact_all.rdata")
# contact_all <- contact_all["THA"]

# I. Stringency index
source("code/0_2_StringencyIndex.R")
# H. Google mobility data
# source("code/0_3_Mobility.R")
contact_schedule <- read_rds(paste0(data_path, "c_schedule.rds")) 

#### K. Epi parameters ####
source("code/0_5_EpiParams.R")

#### L. Burden processes #### 
country_list <- read_rds(paste0(data_path, "country_list_thailand.rds"))
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
voc_features_test <- read_rds(paste0(data_path, "voc_features_test_thailand.rds")) %>% 
  mutate(change_u = if_else(voc_name == "omicron", 1, change_u),
         change_severity = if_else(voc_name == "omicron", 0.4, change_severity))

load(paste0(data_path, "severe_strain_thailand.rdata"))

#### move these things here to facilitate parallel ####
date_start = "2020-03-01"
fit_results_dir <- paste0("fit/")

out_all <- paste0(fit_results_dir, list.files(fit_results_dir, pattern = ".rds")) %>% 
  map(read_rds) %>% 
  map(., "optim") %>% 
  map(., "bestmem") %>% 
  bind_rows() %>% 
  mutate(fit_end_threshold = list.files(fit_results_dir, pattern = ".rds") %>% 
           gsub("fit_", "", .) %>% 
           gsub(".rds", "", .) %>% 
           as.numeric(),
         seed_raw = ymd(date_start) - 30 + par2 - ymd("2020-01-01"),
         seed_20200101 = as.numeric(seed_raw),
         country = "Thailand",
         continent = "Asia",
         country_code = "THA") %>% 
  rename(R0_assumed_2 = par1)

# Create panels for baseline (no vaccination), WHO scenario, annual scenarios
panel_WHO <- expand.grid(cov_2024 = c(seq(0.1, 0.8, 0.1)), 
                         start_age_annual = 60,
                         start_age_6m = 75) %>%
  mutate(scenario = "WHO")

panel_additional <- expand.grid(cov_2024 = c(seq(0.2, 0.8, 0.1)), 
                                start_age_annual = seq(0,75,by=5),
                                start_age_6m = 80) %>% # i.e. only annual vaccination
  mutate(scenario = paste(as.character(start_age_annual),"y+"))

panel_baseline <- data.frame(cov_2024 = 0,
                             start_age_annual = 80,
                             start_age_6m = 80,
                             scenario = "base_case")

panel_final <- bind_rows(panel_baseline,panel_WHO,panel_additional) %>%
  arrange(scenario, cov_2024)

grid_table <- CJ(fit_table_index = 1:nrow(out_all),
                 panel_final_index = 1:nrow(panel_final)) %>% 
  rownames_to_column(var = "grid_table_index")

label_age <- data.frame(group = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                                  "30-34", "35-39", "40-44", "45-49", "50-54", "55"))
