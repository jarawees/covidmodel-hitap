source("code/0_LoadAll_hitap.R")

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")


# PANEL for baseline (no vaccination), WHO scenario, annual scenarios
panel_WHO <- expand.grid(cov_2024 = c(seq(0.1, 0.8, 0.1)), 
                         start_age_annual = 55,
                         start_age_6m = 75) %>%
  mutate(scenario = "WHO")

panel_additional <- expand.grid(cov_2024 = c(seq(0.1, 0.8, 0.1)), 
                                start_age_annual = c(40,45,50,55,60,65,70,75),
                                start_age_6m = 80) %>% # i.e. only annual vaccination
  mutate(scenario = paste(as.character(start_age_annual),"y+"))

panel_baseline <- data.frame(cov_2024 = 0,
                             start_age_annual = 80,
                             start_age_6m = 80,
                             scenario = "base_case")

# PANEL for all scenarios
panel_final <- bind_rows(panel_baseline,panel_WHO,panel_additional)


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
write_rds(res_all, "data/20231211_res_all.rds")

### Extract future population from CovidM ###
pop_future <- cm_simulate(setting_list[[1]])$dynamics %>% 
  filter(!grepl("case|sever|critical|death", compartment)) %>% 
  filter(!grepl("_p|reported", compartment)) %>% 
  filter(!grepl("foi", compartment)) %>%
  filter(!grepl("subclinical", compartment)) %>% 
  mutate(date = t + ymd("2021-02-01"),
         year = year(date)) %>% 
  group_by(t, group) %>% 
  mutate(total_pop = sum(value)) %>% 
  distinct(t, total_pop, .keep_all = T) %>% 
  select(-c(compartment, population, value))


### OUTPUT TABLE 2024 to 2030 ###

panel_cua <- panel_final

# Calculate total cases, hospitalisations, deaths for each scenario
temp_res <- data.frame(cases = numeric(0),
                       hospital_noICU = numeric(0),
                       ICU = numeric(0),
                       death = numeric(0),
                       year = (0))

for (row in 1:nrow(panel_cua)){
  
  # calculate total cases, hospitalisations, deaths for each scenario
  res_all[[row]] %>%
    filter(date >= "2024-01-01") %>%
    filter(date < "2025-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2024) -> new_row_temp24
  
  res_all[[row]] %>%
    filter(date >= "2025-01-01") %>%
    filter(date < "2026-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2025) -> new_row_temp25
  
  res_all[[row]] %>%
    filter(date >= "2024-01-01") %>%
    filter(date < "2025-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2024) -> new_row_temp24
  
  res_all[[row]] %>%
    filter(date >= "2024-01-01") %>%
    filter(date < "2025-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2024) -> new_row_temp24
  
  res_all[[row]] %>%
    filter(date >= "2024-01-01") %>%
    filter(date < "2025-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2024) -> new_row_temp24
  
  temp_res <- rbind(temp_res,new_row_temp)
  
}

 
panel_cua <- cbind(panel_cua,temp_res)


### OUTPUT TABLE 2024-2030 ###




### CHART 2024 ###


# Prepare the dataframe (HOSPITALISATION ONLY)

hosp_chart <- panel_cua %>%
  select(cov_2024, scenario, hospital_noICU) %>%
  filter(cov_2024 > 0) %>%
  pivot_wider(names_from = cov_2024, values_from = hospital_noICU)%>%
  pivot_longer(cols = c(2:9),
               names_to = "coverage",
               values_to = "hospn") %>%
  mutate(coverage = as.numeric(coverage))
hosp_base <- panel_cua[1,"hospital_noICU"]

# Scatter plot
p_hosp <- ggplot(hosp_chart, aes(x=scenario, y=hospn, col=coverage)) +
  geom_point(size=3) +
  scale_shape_manual(values=4) +
  scale_color_gradient2(midpoint=0.5, low="indianred4", mid="khaki1", high = "mediumseagreen", 
                        name = "Coverage") +
  labs(x = "Scenario", 
       y = "Number of hospitalisations (2024)") + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), # set scale to millions
                     limits = c(9800000,10200000),
                     breaks = seq(8000000,11000000,100000)) +
  theme_minimal_hgrid() +
  theme(text = element_text(size = 12),
        axis.line.x.bottom=element_line(color="white")) +
  geom_hline(yintercept=hosp_base, linetype="dashed", color = "grey")
p_hosp


# Prepare the dataframe (DEATH ONLY)

death_chart <- panel_cua %>%
  select(cov_2024, scenario, death) %>%
  filter(cov_2024 > 0) %>%
  pivot_wider(names_from = cov_2024, values_from = death)%>%
  pivot_longer(cols = c(2:9),
               names_to = "coverage",
               values_to = "dead") %>%
  mutate(coverage = as.numeric(coverage))
dead_base <- panel_cua[1,"death"]

# Scatter plot
p_dead <- ggplot(death_chart, aes(x=scenario, y=dead, col=coverage)) +
  geom_point(size=3) +
  scale_shape_manual(values=4) +
  scale_color_gradient2(midpoint=0.5, low="indianred4", mid="khaki1", high = "mediumseagreen", 
                        name = "Coverage") +
  labs(x = "Scenario", 
       y = "Number of deaths (2024)") + 
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), # set scale to millions
                     limits = c(4620000,4670000),
                     breaks = seq(4620000,4670000,10000)) +
  theme_minimal_hgrid() +
  theme(text = element_text(size = 12),
        axis.line.x.bottom=element_line(color="white")) +
  geom_hline(yintercept=dead_base, linetype="dashed", color = "grey")
p_dead
