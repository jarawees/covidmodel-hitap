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
panel_final <- bind_rows(panel_baseline,panel_WHO,panel_additional) %>%
  arrange(scenario, cov_2024)


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

years <- c(2024,2025,2026,2027,2028,2029,2030)
panel_cua <- panel_final[rep(1:nrow(panel_final),times = length(years)),]

# Calculate total cases, hospitalisations, deaths for each scenario
temp_res <- data.frame(cases = numeric(0),
                       hospital_noICU = numeric(0),
                       ICU = numeric(0),
                       death = numeric(0),
                       year = numeric(0))

for (row in 1:nrow(panel_final)){
  
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
    filter(date >= "2026-01-01") %>%
    filter(date < "2027-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2026) -> new_row_temp26
  
  res_all[[row]] %>%
    filter(date >= "2027-01-01") %>%
    filter(date < "2028-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2027) -> new_row_temp27
  
  res_all[[row]] %>%
    filter(date >= "2028-01-01") %>%
    filter(date < "2029-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2028) -> new_row_temp28
  
  res_all[[row]] %>%
    filter(date >= "2029-01-01") %>%
    filter(date < "2030-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2029) -> new_row_temp29
  
  res_all[[row]] %>%
    filter(date >= "2030-01-01") %>%
    filter(date < "2031-01-01") %>%
    summarise(cases = sum(cases),
              hospital_noICU = sum(severe_all),
              ICU = sum(critical_all),
              death = sum(death_all),
              year = 2030) -> new_row_temp30
  
  temp_res <- rbind(temp_res,new_row_temp24,new_row_temp25,new_row_temp26,
                    new_row_temp27,new_row_temp28,new_row_temp29,new_row_temp30)
  
}

temp_res <- temp_res %>% arrange(year)

panel_cua <- cbind(panel_cua,temp_res) %>%
  select(-c(start_age_annual,start_age_6m)) %>%
  arrange(year,scenario,cov_2024) %>%
  rename(coverage = cov_2024)


### HOSPITALISATION CHARTS 2024 ###
p_load(scales)

#1 Scatter chart for 2024

# Prepare the dataframe
hosp_chart <- panel_cua %>%
  filter(year == 2024, coverage > 0) %>%
  select(coverage, scenario, hospital_noICU) %>%
  mutate(coverage = as.numeric(coverage)) %>%
  rename(hospn = hospital_noICU)

hosp_base <- with(panel_cua, hospital_noICU[which (scenario == "base_case" )])[1] # base case (no vaccination) for the year 2024

# Scatter plot
p_hosp <- ggplot(hosp_chart, aes(x=scenario, y=hospn, col=coverage)) +
  geom_point(size=3) +
  scale_shape_manual(values=4) +
  scale_color_gradient2(midpoint=0.5, low="indianred4", mid="khaki1", high = "mediumseagreen", 
                        name = "Coverage") +
  labs(x = "Scenario", 
       y = "Number of hospitalisations (2024)") + 
  
  theme_minimal_hgrid() +
  theme(text = element_text(size = 12),
        axis.line.x.bottom=element_line(color="white")) +
  geom_hline(yintercept=hosp_base, linetype="dashed", color = "grey")
p_hosp

#2 Bar chart for 2024-2030

# scenarios to compare
comparison_bar <- c("WHO","50 y+")
hosp_bar <- panel_cua %>%
  filter(scenario %in% comparison_bar, coverage == 0.5) %>%
  select(scenario, year, hospital_noICU)

# Bar chart
bar_hosp <- ggplot(hosp_bar, aes(x=year,y=hospital_noICU,fill=scenario)) +
  geom_bar(stat="identity",position=position_dodge())  +
  labs(x = "Year", 
       y = "Number of hospitalisations") + 
  scale_y_continuous(labels = label_comma())
bar_hosp
