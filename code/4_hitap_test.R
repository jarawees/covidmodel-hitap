source("code/0_LoadAll.R")

out <- read_rds("data/out.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")

#### Scenarios for 5y+ up to 60y+ ####
ages <- c("60y+","55y+","50y+","45y+","40y+","35y+","30y+","25y+","20y+","15y+","10y+","5y+")

# add a list of NAs if the age group is given a booster vaccine, 1 if the group is vaccinated
age_list <- list()
for(j in 1:length(ages)){
  na_j <- 13 - j
  ones_j <- 3 + j
  age_list[[j]] <- c(rep(NA,na_j),rep(1,ones_j))
}

panel <- expand.grid(f = c(1:2), boosting_level = c(0.00001, seq(0.1, 0.9, 0.1)), prioritisation = age_list)
setting_list <- list()
for(i in 1:nrow(panel)) {
      setting_list[[i]]<- parameterise_setting(
        f = panel$f[i],
        prioritisation_followup = panel$prioritisation[[i]],
        boosting_level = panel$boosting_level[i]
      )
    }


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
