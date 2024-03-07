source("code/0_LoadAll_hitap.R")

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31", "2024-01-01")

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
  setting_list[[i]] <- gen_country_basics(country = "Thailand",
                                          R0_assumed = out$optim$bestmem[1],
                                          date_start = "2021-02-01",
                                          date_end = "2030-12-31",
                                          contact = contact_schedule %>% mutate(work = if_else(date >= "2023-01-01",1,work), other = if_else(date >= "2023-01-01",1,other)),
                                          processes = gen_burden_processes(VE = efficacy_all),
                                          period_wn  = 3*365, # duration, waning of natural immunity
                                          period_wv_m2l = 1*365, # duration, waning from medium to low levels vaccine induced 
                                          period_wv_h2m = 1*365, # duration, waning from high to medium levels vaccine induced 
                                          prob_v_p_2l = 1,
                                          prob_v_p_2m = 0,
                                          prob_v_b_l2m = 0.5,
                                          deterministic = TRUE,
                                          scenario_primary = scenario3_primary,
                                          scenario_booster = scenario3_booster %>% mutate(prob_v_b_l2m = if_else(date >= "2024-01-01", 0, prob_v_b_l2m)),
                                          seed = out$optim$bestmem[2]) %>% 
    update_u_y(para = .,
               date_switch = c("2021-01-15", "2021-07-05", "2021-12-31", "2024-01-01"),
               rc_u = c(1, 1.5, 1, 1), # relative changes in u
               rc_y = c(1, 1, 1, 1), # relative changes in y
               # rc_ve = c(1, 0.79, 0.79, 0.58), # update rc_ve using meta-analysis from Dec2023
               rc_ve = c(1, 0.79, 0.79, 1), # base case
               efficacy_baseline = efficacy_all
    ) %>%
    emerge_VOC_burden(para = .,
                      # rc_severity = c(1, 1.5, 0.7, 0.7), # relative change in ihr and ifr
                      rc_severity = c(1, 1.5, 0.7, 1), # base case
                      efficacy_baseline = efficacy_all) %>%
    vaccinate_primary(para = .,
                      vac_data = owid_vac,
                      values = primary_allocation_plan) %>%
    vaccinate_additional(para = .,
                         vac_data = owid_vac,
                         booster_plan = booster_allocation_plan,
                         start_age_annual = panel_final$start_age_annual[i],
                         start_age_6m = panel_final$start_age_6m[i],
                         cov_2024 = panel_final$cov_2024[i],
                         month_annual = c(5:6),
                         month_6m = c(11:12))
    # parameterise_setting(
    # start_age_annual = panel_final$start_age_annual[i],
    # start_age_6m = panel_final$start_age_6m[i],
    # cov_2024 = panel_final$cov_2024[i])
}

### Populate the model with the parameters in setting_list ###
res_all <- list()

for(i in 1:length(setting_list)){
# for(i in 1:2){
  cm_simulate(setting_list[[i]])$dynamics %>% 
    mutate(date = t + ymd("2021-02-01"),
           year = year(date)) %>% 
    pivot_wider(names_from = compartment, values_from = value) %>% 
    mutate(
      severe_all = case_when(date <= date_switch[1] ~ severe_i,
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
                                 date > date_switch[4] ~ death_voc4_o),
           heathly_all = S + Sv_l + Sv_m + Sv_h + E + Ev_l + Ev_m + Ev_h + R + Rv_l + Rv_m + Rv_h,
           population_all = S + Sv_l + Sv_m + Sv_h + E + Ev_l + Ev_m + Ev_h + Ip + Ip_l + Ip_m + Ip_h + Is + Is_l + Is_m + Is_h + Ia + Ia_l + Ia_m + Ia_h + R + Rv_l + Rv_m + Rv_h) %>% 
    dplyr::select(date, year, group, cases, ends_with("all")) -> res_all[[i]]
}

# write_rds(res_all, "data/20240129_res_all.rds")

## scenario 40 y+ coverage 80%
# res_all[[8]] %>% 
#   pivot_longer(cols = c("cases", "severe_all", "critical_all", "death_all")) %>% 
#   group_by(year, name) %>% 
#   summarise(value = sum(value)) %>% 
#   pivot_wider(names_from = name,
#               values_from = value)
# 
# tmp <- cm_simulate(setting_list[[i]])$dynamics %>% 
#   mutate(date = t + ymd("2021-02-01"),
#          year = year(date))
# 
# compartments_status <- c("S", "Sv_l", "Sv_m", "Sv_h",
#                          "E", "Ev_l", "Ev_m", "Ev_h",
#                          "Ip", "Ip_l", "Ip_m", "Ip_h",
#                          "Is", "Is_l", "Is_m", "Is_h",
#                          "Ia", "Ia_l", "Ia_m", "Ia_h",
#                          "R", "Rv_l", "Rv_m", "Rv_h")
# 
# tmp %>% 
#   dplyr::filter(compartment %in% compartments_status) %>% 
#   dplyr::filter(group == "60-64") %>% 
#   mutate(compartment_general = substr(compartment, 1, 1)) %>% 
#   ggplot(., aes(x = date, y = value, color = compartment, fill = compartment)) +
#   geom_bar(position = "stack", stat = "identity") +
#   facet_wrap(~compartment, scales = "free")


### Extract future population from CovidM ###
# res_all[[1]] %>% 
#   ggplot(., aes(x = date, y = severe_all, group = group)) +
#   geom_line() +
#   facet_wrap(~group, )
# 
# pop_future <- cm_simulate(setting_list[[1]])$dynamics %>% 
#   filter(!grepl("case|sever|critical|death", compartment)) %>% 
#   filter(!grepl("_p|reported", compartment)) %>% 
#   filter(!grepl("foi", compartment)) %>%
#   filter(!grepl("subclinical", compartment)) %>% 
#   mutate(date = t + ymd("2021-02-01"),
#          year = year(date)) %>% 
#   group_by(t, group) %>% 
#   mutate(total_pop = sum(value)) %>% 
#   distinct(t, total_pop, .keep_all = T) %>% 
#   select(-c(compartment, population, value))
# 
# pop_future %>% 
#   ggplot(., aes(x = date, y = total_pop, group = group)) +
#   geom_line() +
#   facet_wrap(~group)


### OUTPUT TABLE 2024 to 2030 ###
# res_all <- readRDS("data/20240129_res_all.rds")

years <- c(2023:2030)
# panel_cua <- panel_final[rep(1:nrow(panel_final), times = length(years)),]

# years <- sort(rep(years, times = nrow(panel_final)))
# 
# panel_cua <- panel_cua %>% 
#   mutate(year = years)

# Calculate total cases, hospitalisations, deaths for each scenario
temp_res <- data.frame(
  year = numeric(0),
  healthy_prob = numeric(0),
  cases_prob = numeric(0),
  hospital_noICU_prob = numeric(0),
  ICU_prob = numeric(0),
  death_prob = numeric(0),
  total_pop = numeric(0),
  healthy = numeric(0),
  cases = numeric(0),
  hospital_noICU = numeric(0),
  ICU = numeric(0),
  death = numeric(0)
)

for(row in 1:nrow(panel_final)){
  for(y in years){
    res_all[[row]] %>%
      filter(year == y) %>% 
      group_by(date) %>% 
      summarise(healthy = sum(heathly_all),
                cases = sum(cases),
                hospital_noICU = sum(severe_all),
                ICU = sum(critical_all),
                death = sum(death_all),
                pop = sum(population_all)) %>% 
      mutate(healthy_prob = healthy/pop,
             cases_prob = cases/pop,
             hospital_noICU_prob = hospital_noICU/pop,
             ICU_prob = ICU/pop,
             death_prob = death/pop,
             year = y) %>% 
      ungroup() %>% 
      group_by(year) %>% 
      summarise(total_pop = mean(pop),
                healthy = mean(healthy),
                cases = sum(cases),
                hospital_noICU = sum(hospital_noICU),
                ICU = sum(ICU),
                death = sum(death),
                healthy_prob = mean(healthy_prob),
                cases_prob = mean(cases_prob),
                hospital_noICU_prob = mean(hospital_noICU_prob),
                ICU_prob = mean(ICU_prob),
                death_prob = mean(death_prob)) %>% 
      ungroup() %>% 
      mutate(cov_2024 = panel_final$cov_2024[row],
             start_age_annual = panel_final$start_age_annual[row],
             start_age_6m = panel_final$start_age_6m[row],
             scenario = panel_final$scenario[row]) -> temp_res_year
    temp_res <- bind_rows(temp_res, temp_res_year)
  }
}

result_cua <- temp_res %>% 
  arrange(year,scenario,cov_2024) %>% 
  select(year, cov_2024, start_age_annual, start_age_6m, scenario, everything()) %>% 
  filter(year != 2023 | (year == 2023 & cov_2024 == 0.1 & scenario == "40 y+"))
  
write_csv(result_cua, "data/20230307_Panel for CUA_base case.csv")


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
