# PANEL for baseline (no vaccination), WHO scenario, annual scenarios
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


para <- gen_country_basics( date_start = "2020-01-01",
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


res <- cm_simulate(para)

res$dynamics %>% 
  dplyr::filter(compartment %in% compartment_pop) %>% 
  group_by(t, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  ggplot(., aes(x = t, y = value, color = compartment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~compartment)
