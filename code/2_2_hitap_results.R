source("code/0_LoadAll_hitap.R")

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")


# PANEL for WHO scenario
panel_WHO <- expand.grid(cov_2024 = c(0.00001, seq(0.1, 0.8, 0.1)), 
                         start_age_annual = 55,
                         start_age_6m = 75)


# PANEL for annual scenarios
panel_additional <- expand.grid(cov_2024 = c(0.00001, seq(0.1, 0.8, 0.1)), 
                                start_age_annual = c(40,45,50,55,60,65,70,75),
                                start_age_6m = 80) # i.e. only annual vaccination

# PANEL merged for all scenarios
panel_final <- bind_rows(panel_WHO,panel_additional)


# Parameters for each scenario
setting_list <- list()
for(i in 1:nrow(panel_final)) {
  setting_list[[i]]<- parameterise_setting(
    start_age_annual = panel_final$start_age_annual[i],
    start_age_6m = panel_final$start_age_6m[i],
    cov_2024 = panel_final$cov_2024[i])
}


### Populate the model with the parameters in setting_list ###
res_all <- list()
for(i in 1:length(setting_list)){
  cm_simulate(setting_list[[i]])$dynamics %>% 
    filter(grepl("case|sever|critical|death", compartment)) %>% 
    filter(!grepl("_p|reported", compartment)) %>% 
    mutate(date = t + ymd("2021-02-01"),
           year = year(date))  %>% 
    pivot_wider(names_from = compartment, values_from = value) %>% 
    mutate(severe_all = case_when(date <= date_switch[1] ~ severe_i,
                                  date > date_switch[1] & date <= date_switch[2] ~ severe_voc1_i,
                                  date > date_switch[2] & date <= date_switch[3] ~ severe_voc2_i,
                                  date > date_switch[3] ~ severe_voc3_i),
           critical_all = case_when(date <= date_switch[1] ~ critical_i,
                                    date > date_switch[1] & date <= date_switch[2] ~ critical_voc1_i,
                                    date > date_switch[2] & date <= date_switch[3] ~ critical_voc2_i,
                                    date > date_switch[3] ~ critical_voc3_i),
           death_all = case_when(date <= date_switch[1] ~ death_o,
                                 date > date_switch[1] & date <= date_switch[2] ~ death_voc1_o,
                                 date > date_switch[2] & date <= date_switch[3] ~ death_voc2_o,
                                 date > date_switch[3] ~ death_voc3_o)) %>% 
    dplyr::select(date, year, group, cases, ends_with("all")) -> res_all[[i]]
}


### OUTPUT TABLE ###
# IMPORTANT: Estimates are from 2023 onwards (i.e. not Oct-Dec 2022)

panel_cua <- panel_final
pop_2024 <- 71886000 # estimate from UN data for 2024 https://data.un.org/Data.aspx?q=thailand&d=PopDiv&f=variableID%3A12%3BcrID%3A764
                  # ! update population from CovidM if possible

# Calculate total cases, hospitalisations, deaths for each scenario
temp_res <- data.frame(cases = numeric(0),
                       hospital_noICU = numeric(0),
                       ICU = numeric(0),
                       death = numeric(0))

for (row in 1:nrow(panel_cua)){
  
  # calculate total cases, hospitalisations, deaths for each scenario
  res_all[[row]] %>%
    filter(date >= "2024-01-01") %>%
    filter(date < "2025-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all)) -> new_row_temp
  temp_res <- rbind(temp_res,new_row_temp)
  
}

panel_cua <- cbind(panel_cua,temp_res) %>%
  mutate(p_dead = death/pop_2024,
         p_ICU = ICU/pop_2024,
         p_hosp = hospital_noICU/pop_2024,
         p_case = (cases - hospital_noICU - ICU - death)/pop_2024)
  
  
  
  
