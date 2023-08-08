source("code/0_LoadAll.R")

out <- read_rds("data/out.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")

#### Scenarios for 5y+ up to 60y+ ####
ages <- c("60y+","55y+","50y+","45y+","40y+","35y+","30y+","25y+","20y+","15y+","10y+","05y+")

# add a list of NAs if the age group is given a booster vaccine, 1 if the group is vaccinated
age_list <- list()
for(j in 1:length(ages)){
  na_j <- 13 - j
  ones_j <- 3 + j
  age_list[[j]] <- c(rep(NA,na_j),rep(1,ones_j))
}

### Set up a panel of the different scenarios ###
# f = 1 annual, f = 2 biannual, boosting_level sets coverage from 10% to 90%, prioritisation is age group
panel <- expand.grid(f = c(1:2), boosting_level = c(0.00001, seq(0.1, 0.9, 0.1)), prioritisation = age_list)


### Parameters for each of the scenarios being tested ###
setting_list <- list()
for(i in 1:nrow(panel)) {
      setting_list[[i]]<- parameterise_setting(
        f = panel$f[i],
        prioritisation_followup = panel$prioritisation[[i]],
        boosting_level = panel$boosting_level[i]
      )
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


### Create a list of matrices for the heatmap ###
# IMPORTANT: Estimates are from 2023 onwards (i.e. not Oct-Dec 2022)

panel_heatmap <- rownames_to_column(panel) # adds a column with row number for reference


# Calculate total cases, hospitalisations, deaths for each scenario
temp_res <- data.frame(cases = numeric(0),
                       hospitalisation = numeric(0),
                       death = numeric(0))

for (row in 1:nrow(panel_heatmap)){
  
  # calculate total cases, hospitalisations, deaths for each scenario
  res_all[[row]] %>%
    filter(date >= "2023-01-01") %>%
    summarise(cases = sum(cases),
              hospitalisation = sum(severe_all) + sum(critical_all),
              death = sum(death_all)) -> new_row_temp
  temp_res <- rbind(temp_res,new_row_temp)
  
}

panel_heatmap <- cbind(panel_heatmap,temp_res)


# Prepare the dataframe for heatmap (HOSPITALISATION ONLY)

hosp_heatmap <- panel_heatmap %>%
  select(boosting_level, f, prioritisation, hospitalisation) %>%
  pivot_wider(names_from = prioritisation, values_from = hospitalisation) %>%
  rename_with(~ages, 3:14) %>%
  slice(3:20) %>%
  pivot_longer(cols = all_of(ages),
               names_to = "target_pop",
               values_to = "hospn") 
  

# Heatmap graph
p_hosp <- ggplot(hosp_heatmap, aes(x = target_pop, y = boosting_level, fill = hospn)) +
  geom_tile() + scale_fill_distiller(palette = "BuGn") +
  labs(title = "Hospitalisations (2023-2030)",
       x = "Target population", 
       y = "Coverage (as a percentage of the eligible population)") +
  facet_grid(rows = hosp_heatmap$f, labeller = as_labeller(c(
    "1"= "Annual", "2" = "Biannual")))
  

# Prepare the dataframe for heatmap (DEATH ONLY)

dead_heatmap <- panel_heatmap %>%
  select(boosting_level, f, prioritisation, death) %>%
  pivot_wider(names_from = prioritisation, values_from = death) %>%
  rename_with(~ages, 3:14) %>%
  slice(3:20) %>%
  pivot_longer(cols = all_of(ages),
               names_to = "target_pop",
               values_to = "dead") 


# Heatmap graph
p_dead <- ggplot(dead_heatmap, aes(x = target_pop, y = boosting_level, fill = dead)) +
  geom_tile() + scale_fill_distiller(palette = "YlOrRd") +
  labs(title = "Deaths (2023-2030)",
       x = "Target population", 
       y = "Coverage (as a percentage of the eligible population)") +
  facet_grid(rows = hosp_heatmap$f, labeller = as_labeller(c(
    "1"= "Annual", "2" = "Biannual")))

