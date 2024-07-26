source("code/0_LoadAll_hitap.R")

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")

# Create panels for baseline (no vaccination), WHO scenario, annual scenarios
panel_WHO <- expand.grid(cov_2024 = c(seq(0.1, 0.8, 0.1)), 
                         start_age_annual = 60,
                         start_age_6m = 75) %>%
  mutate(scenario = "WHO")

panel_additional <- expand.grid(cov_2024 = c(seq(0.2, 0.8, 0.1)), 
                                start_age_annual = seq(0,75,by=5),
                                start_age_6m = 80) %>% # i.e. only annual vaccination
  mutate(scenario = paste(as.character(start_age_annual),"y+"))

panel_baseline <- data.frame(cov_2024 = 0,
                             start_age_annual = 80,
                             start_age_6m = 80,
                             scenario = "base_case")

panel_final <- bind_rows(panel_baseline,panel_WHO,panel_additional) %>%
  arrange(scenario, cov_2024)

# Create set of parameter settings based on all scenarios
setting_list <- list()

for(i in 1:nrow(panel_final)){
  setting_list[[i]] <- gen_country_basics( date_start = "2020-01-01",
                            date_end = "2030-12-31",
                            processes_set = burden_processes_all,
                            prob_v_p_2l = 0.33,
                            prob_v_p_2m = 0.33,
                            prob_v_b_l2m = 0,
                            # seed = fitted_table_baseline$seed_20200101,
                            # R0_assumed = fitted_table_baseline$R0_assumed_2,
                            # period_wn = fitted_table_baseline$wn,
                            fitted_table_tmp = fitted_table_baseline,
                            period_wv_h2m = 1*365, 
                            period_wv_m2l = 1*365, 
                            deterministic = TRUE) %>% 
  update_u_y(para = .,
             country_tmp = "Thailand",
             country_code_tmp = "THA",
             detection_threshold = 0.3,
             voc_features_inuse = voc_features_test,
             future_severe = F,
             future_severe_lvl = "mean", # sensitivity analysis 4
             efficacy_baseline = efficacy_all# vaccine efficacy
  ) %>% 
  emerge_voc_burden(para = .,
                    detection_threshold = 0.3,
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
}

# # debug
# res <- cm_simulate(para)
# 
# res$dynamics %>% 
#   dplyr::filter(compartment %in% compartment_pop) %>% 
#   group_by(t, population, compartment) %>% 
#   summarise(value = sum(value)) %>% 
#   ggplot(., aes(x = t, y = value, color = compartment)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~compartment)

# Populate the model with the parameters in setting_list and store # of daily severe, critical, death, healthy, population
res_all <- list()

for(i in 1:length(setting_list)){
  cm_simulate(setting_list[[i]])$dynamics %>% 
    mutate(date = t + ymd("2021-02-01"),
           year = year(date)) %>% 
    pivot_wider(names_from = compartment, values_from = value) %>% 
    mutate(
      severe_all = case_when(date <= date_switch[1] ~ severe_i,
                                  date > date_switch[1] & date <= date_switch[2] ~ severe_voc_alpha_i,
                                  date > date_switch[2] & date <= date_switch[3] ~ severe_voc_delta_i,
                                  date > date_switch[3] ~ severe_voc_omicron_i),
           critical_all = case_when(date <= date_switch[1] ~ critical_i,
                                    date > date_switch[1] & date <= date_switch[2] ~ critical_voc_alpha_i,
                                    date > date_switch[2] & date <= date_switch[3] ~ critical_voc_delta_i,
                                    date > date_switch[3] ~ critical_voc_omicron_i),
           death_all = case_when(date <= date_switch[1] ~ death_o,
                                 date > date_switch[1] & date <= date_switch[2] ~ death_voc_alpha_o,
                                 date > date_switch[2] & date <= date_switch[3] ~ death_voc_delta_o,
                                 date > date_switch[3] ~ death_voc_omicron_o),
      heathly_all = S + Sv_l + Sv_m + Sv_h + E + Ev_l + Ev_m + Ev_h + R + Rv_l + Rv_m + Rv_h,
      population_all = S + Sv_l + Sv_m + Sv_h + E + Ev_l + Ev_m + Ev_h + Ip + Ip_l + Ip_m + Ip_h + Is + Is_l + Is_m + Is_h + Ia + Ia_l + Ia_m + Ia_h + R + Rv_l + Rv_m + Rv_h) %>% 
    dplyr::select(date, year, group, cases, ends_with("all")) -> res_all[[i]]
}

# Calculate total cases, hospitalisations, deaths for each scenario
years <- c(2023:2030)

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
      group_by(date, group) %>% 
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
      group_by(year, group) %>% 
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
  filter(year != 2023 | (year == 2023 & cov_2024 == 0.2 & scenario == "40 y+"))
