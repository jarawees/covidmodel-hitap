source("code/0_LoadAll.R")

date_start = "2021-02-15"
fit_results_dir <- paste0("fit/")

out_all <- paste0(fit_results_dir, list.files(fit_results_dir, pattern = "fit_0")) %>% 
  map(read_rds) %>% 
  map(., "optim") %>% 
  map(., "bestmem") %>% 
  bind_rows() %>% 
  mutate(fit_end_threshold = list.files(fit_results_dir, pattern = "fit_0") %>% 
           gsub("fit_", "", .) %>% 
           gsub(".rds", "", .) %>% 
           as.numeric(),
         seed_raw = ymd(date_start) - 30 + par2 - ymd("2020-01-01"),
         seed_20200101 = as.numeric(seed_raw),
         country = "Thailand",
         continent = "Asia",
         country_code = "THA") %>% 
  rename(R0_assumed_2 = par1,
         wn = par4)

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
                                           fitted_table_tmp = out_all[1,], # this needs to vary over the possible fitted results
                                           period_wv_h2m = 1*365, 
                                           period_wv_m2l = 1*365, 
                                           deterministic = TRUE) %>% 
    update_u_y(para = .,
               voc_features_inuse = voc_features_test %>% mutate(change_u = 1),
               future_severe = F,
               future_severe_lvl = "mean",
               efficacy_baseline = efficacy_all 
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

# Populate the model with the parameters in setting_list and store # of daily severe, critical, death, healthy, population
res_all <- list()

for(i in 1:length(setting_list)){
  cm_simulate(setting_list[[i]])$dynamics %>% 
    aggregate_results(dynamics_tmp = .) -> res_all[[i]]
}



res_all[[i]] %>% 
  dplyr::select(-prop, -cohort_all) %>% 
  pivot_wider(names_from = compartment,
              values_from = incidence) %>% 
  dplyr::filter(critical < death, group_index < 12)

res_all %>% 
  bind_rows(.id = "scenario_id") %>% 
  left_join(panel_final %>% 
              rownames_to_column(var = "scenario_id"),
            by = "scenario_id") -> output

output %>% 
  group_by(scenario_id, cov_2024, start_age_annual, start_age_6m, scenario, year, compartment) %>% 
  summarise(incidence = sum(incidence),
            cohort_all = sum(cohort_all)) %>% 
  dplyr::filter(compartment == "death") %>% 
  ggplot(., aes(x = year, y = incidence/cohort_all, colour = cov_2024, group = cov_2024)) +
  geom_line() +
  facet_wrap(~scenario, scales = "free")
