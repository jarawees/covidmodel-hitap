options(scipen=999)

## Load DDC MOPH data
# download the data file from DDC MOPH TH if it doesn't exist in your directory
if(!file.exists(paste0("data/epi_round1to2.csv"))){
  res_round1to2 <- GET("https://covid19.ddc.moph.go.th/api/Cases/round-1to2-all")
  df_round1to2 <- fromJSON(rawToChar(res_round1to2$content))[,1:10]
  write.csv(df_round1to2, file = "data/epi_round1to2.csv", row.names = F)
}

if(!file.exists(paste0("data/epi_round3.csv"))){
  res_round3 <- GET("https://covid19.ddc.moph.go.th/api/Cases/report-round-3-y21-line-lists")
  df_round3 <- fromJSON(rawToChar(res_round3$content))[["data"]][,1:10]
  write.csv(df_round3, file = "data/epi_round3.csv", row.names = F)
}

if(!file.exists(paste0("data/epi_update.csv"))){
  res_update <- GET("https://covid19.ddc.moph.go.th/api/Cases/timeline-cases-all")
  df_update <- fromJSON(rawToChar(res_update$content))[,1:10]
  write.csv(df_update, file = "data/epi_update.csv", row.names = F)
}

# update the data file from DDC MOPH TH if the time difference is greater than a week
if(as.numeric(abs(as.Date(file.info(paste0("data/epi_update.csv"))$mtime) -
                  as.Date(Sys.time()))) > 7){
  res_update <- GET("https://covid19.ddc.moph.go.th/api/Cases/timeline-cases-all")
  df_update <- fromJSON(rawToChar(res_update$content))[,1:10]
  write.csv(df_update, file = "data/epi_update.csv", row.names = F)
}

epi_round1to2 <- fread("data/epi_round1to2.csv") # weekly data from 2020-01-12 to 2021-03-31
epi_round3 <- fread("data/epi_round3.csv") # weekly data from 2021-04-01 to 2021-12-31
epi_update <- fread("data/epi_update.csv") # weekly data from 2022-01-01 onward

epi_ddc <- bind_rows(epi_round1to2, epi_round3, epi_update)
rm(epi_round1to2, epi_round3, epi_update)


## Load model results
source("code/0_LoadAll_hitap.R")

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31", "2025-01-01")

panel_baseline <- data.frame(cov_2024 = 0,
                             start_age_annual = 80,
                             start_age_6m = 80,
                             scenario = "base_case")

setting_valid <- parameterise_setting(
  start_age_annual = panel_baseline$start_age_annual,
  start_age_6m = panel_baseline$start_age_6m,
  cov_2024 = panel_baseline$cov_2024)

res_valid <- cm_simulate(setting_valid)$dynamics %>% 
  filter(grepl("case|sever|critical|death", compartment)) %>% 
  filter(!grepl("_p|reported", compartment)) %>% 
  mutate(date = t + ymd("2021-02-01"),
         year = year(date))  %>% 
  pivot_wider(names_from = compartment, values_from = value) %>% 
  mutate(severe_all = case_when(date <= date_switch[1] ~ severe_i,
                                date > date_switch[1] & date <= date_switch[2] ~ severe_voc1_i,
                                date > date_switch[2] & date <= date_switch[3] ~ severe_voc2_i,
                                date > date_switch[3] & date <= date_switch[4] ~ severe_voc3_i,
                                date > date_switch[4] ~ severe_voc4_i),
         critical_all = case_when(date <= date_switch[1] ~ critical_i,
                                  date > date_switch[1] & date <= date_switch[2] ~ critical_voc1_i,
                                  date > date_switch[2] & date <= date_switch[3] ~ critical_voc2_i,
                                  date > date_switch[3] & date <= date_switch[4] ~ critical_voc3_i,
                                  date > date_switch[4] ~ critical_voc4_i),
         death_all = case_when(date <= date_switch[1] ~ death_o,
                               date > date_switch[1] & date <= date_switch[2] ~ death_voc1_o,
                               date > date_switch[2] & date <= date_switch[3] ~ death_voc2_o,
                               date > date_switch[3] & date <= date_switch[4] ~ death_voc3_o,
                               date > date_switch[4] ~ death_voc4_o)) %>% 
  dplyr::select(date, year, group, cases, ends_with("all"))

validation_plot <- function(y = 2022){
df_res_plot <- res_valid %>%
  mutate(weeknum = epiweek(date)) %>% 
  filter(year == y) %>% 
  group_by(weeknum, year) %>% 
  summarise(cases = sum(cases),
            hospital_noICU = sum(severe_all),
            ICU = sum(critical_all),
            death = sum(death_all))

df_ddc_plot <- epi_ddc %>% 
  select(year, weeknum, new_case, new_death) %>% 
  filter(year == y) %>% 
  pivot_longer(cols = new_case:new_death) %>% 
  mutate(name = factor(name, levels = c("new_case", "new_death"), labels = c("Actual cases", "Actual deaths")))

df_res_plot %>% 
  select(year, weeknum, cases, death) %>% 
  pivot_longer(cols = cases:death) %>%
  mutate(name = factor(name, levels = c("cases", "death"), labels = c("Predicted cases", "Predicted deaths"))) %>% 
  ggplot(., aes(x = weeknum, y = value, col = name)) +
  geom_point() +
  labs(x = "Week", y = "Cummulative weekly outcomes", col = "Health outcomes", title = y) +
  theme_light() +
  geom_point(data = df_ddc_plot, aes(x = weeknum, y = value, col = name), alpha = 0.5)
}

validation_plot(2021)
