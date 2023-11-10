source("code/0_LoadAll_hitap.R")
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
