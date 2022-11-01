## Load required packages
if(!require(pacman)) install.packages("pacman")
library(pacman)
p_load(tidyverse, httr, jsonlite, countrycode, data.table, socialmixr,
       lubridate, mgcv)

##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- T
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))
source("code/0_1_util_functions.R")

## Load required data
# A. Covid-19 deaths
# B. Covid-19 cases/hospitalisations
# home isolation policy were proclaimed on 4 Jan 2022, thus, new cases afterwards were mixed between hospitalised and home isolated.

# res_round1to2 <- GET("https://covid19.ddc.moph.go.th/api/Cases/round-1to2-all")
# df_round1to2 <- fromJSON(rawToChar(res_round1to2$content)) 
# write.csv(df_round1to2, file = "data/epi_round1to2.csv", row.names = F)
epi_round1to <- fread("data/epi_round1to2.csv") # snapshot from 2020-01-12 to 2021-03-31

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

epi_update <- fread("data/epi_update.csv") # time-series from 2021-04-01 onward

epi <- bind_rows(epi_round1to, epi_update) %>%
  select(txn_date, new_case_excludeabroad, new_death)

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

pcr_rate <- fread("data/thailand_covid-19_testing_data.csv")[,1:3] %>%
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

geno_freq <- fread("data/sars-cov-2-variants-dmsc.csv")[,1:6] %>%
  mutate_at(vars(Date_Start, Date_End), ~as.Date(., format = "%d/%m/%Y")) %>%
  rename("Alpha" = "B.1.1.7 (Alpha)",
         "Delta" = "B.1617.2 (Delta)",
         "Beta" = "B.1.351 (Beta)",
         "Omicron" = "B.1.1.529 (Omicron") %>%
  filter(!is.na(Alpha))

# F. Vaccine uptake
source("code/0_4_Vaccinations.R")

# G. Contact matrices
# load("data/contact_all.rdata")
# contact_all <- contact_all["THA"]

# I. Stringency index
source("code/0_2_StringencyIndex.R")
# H. Google mobility data
source("code/0_3_Mobility.R")

# J. Population structure
# popTH <- fread("data/pop_str_2021.csv")

# popTH <- popTH %>%
#   gather(key = "sex", value = "pop", male, female, both) %>%
#   mutate_at(vars(pop), ~as.numeric(gsub("[^\\d]+", "", ., perl=TRUE))) %>%
#   mutate(age_group = age %/% 5,
#          age = paste0(age_group * 5, "-", age_group * 5 + 4),
#          age = replace(age, age_group>=18, "90+"),
#          age = factor(age, levels = limits_to_agegroups(seq(0, 90, by = 5))))

#### K. Epi parameters ####
source("code/0_5_EpiParams.R")

#### L. Burden processes #### 
source("code/0_6_HealthCareSystem.R")
