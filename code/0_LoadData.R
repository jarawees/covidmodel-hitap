## Load required packages
if(!require(pacman)) install.packages("pacman")
library(pacman)
p_load(tidyverse, httr, jsonlite, countrycode, gsheet)

## Load required data
# A. Covid-19 deaths
res_round1to2 <- GET("https://covid19.ddc.moph.go.th/api/Cases/round-1to2-all")
df_round1to2 <- fromJSON(rawToChar(res_round1to2$content)) # snapshot from 2020-01-12 to 2021-03-31

res_uptodate <- GET("https://covid19.ddc.moph.go.th/api/Cases/timeline-cases-all")
df_uptodate <- fromJSON(rawToChar(res_uptodate$content)) # time-series from 2021-04-01 onward

# B. Covid-19 hospitalisation

# C. PCR positivity rate
pcr_rate <- gsheet2tbl('https://docs.google.com/spreadsheets/d/13j9SZ01F4ATJT9PEPaVEKhAhouuU6Wo_/edit#gid=1001738166')[,1:3]

# D. Seroprevalence
# snapshot over the past 2 years, surveillance at Nov 2021
sero <- read_csv("data/serosurveillance65.csv")

# E. Genotype frequencies
geno_freq <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1edvjCkHvbaymMaXCJe1ezhfwwTEQSyNF/edit#gid=674137022')

# F. Vaccine uptake
file <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv"
owid_vac <- read_csv(url(file)) %>% 
  mutate(wb = countrycode(location, "country.name", "wb")) %>%
  filter(location == "Thailand")

# Google mobility data
# download the data file from Google if it doesn't exist in your directory
if(!file.exists(paste0("data/gm.csv"))){
  download.file("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
                paste0("data/gm.csv"))
}

# update the data file from google if the time difference is greater than a week
if(as.numeric(abs(as.Date(file.info(paste0("data/gm.csv"))$mtime) -
                  as.Date(Sys.time()))) > 7){
  download.file("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
                paste0("data/gm.csv"))
}

gm <- fread("data/gm.csv")
