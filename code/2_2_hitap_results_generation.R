source("code/0_LoadAll.R")

# Create set of parameter settings based on all scenarios
gen_results_table <- function(i){
  
    index1 <- grid_table[i, ]$fit_table_index
    index2 <- grid_table[i, ]$panel_final_index
    
    tmp <- gen_country_basics(date_start = "2020-01-01",
                              date_end = "2030-12-31",
                              processes_set = burden_processes_all,
                              prob_v_p_2l = 0.33,
                              prob_v_p_2m = 0.33,
                              prob_v_b_l2m = 0,
                              fitted_table_tmp = out_all[index1,], # this needs to vary over the possible fitted results
                              period_wv_h2m = 1*365, 
                              period_wv_m2l = 1*365, 
                              deterministic = TRUE) %>% 
      update_u_y(para = .,
                 voc_features_inuse = voc_features_test %>% mutate(change_u = 1), # we have been told that the original version of voc_features_test overestimate changes in transmissibility 
                 future_severe = F,
                 future_severe_lvl = "mean",
                 efficacy_baseline = efficacy_all 
      ) %>% 
      emerge_voc_burden(para = .,
                        detection_threshold = 0.3,
                        split_E = T,
                        voc_features_inuse = voc_features_test %>% mutate(change_u = 1), #
                        efficacy_baseline = efficacy_all) %>% 
      vaccinate_primary(para = .,
                        vac_data = owid_vac,
                        values = primary_allocation_plan) %>%
      vaccinate_additional(para = .,
                           vac_data = owid_vac,
                           booster_plan = booster_allocation_plan,
                           start_age_annual = panel_final$start_age_annual[index2],
                           start_age_6m = panel_final$start_age_6m[index2],
                           cov_2024 = panel_final$cov_2024[index2],
                           month_annual = c(5:6),
                           month_6m = c(11:12))
    
    cm_simulate(tmp)$dynamics -> tmp_dynamics
    aggregate_results(dynamics_tmp = tmp_dynamics, by = "year") -> res
    
    tmp_fn <- paste0("~/Dropbox/Github_Data/HITAP_CovidM/results/", i, ".csv")
    write_csv(res, file = tmp_fn)
      
    # return(res)

}

cl <- makeCluster(5)
clusterEvalQ(cl, source("code/0_LoadAll.R"))
parLapply(cl, 201:1573, gen_results_table) 

res_all[[i]] %>% 
  dplyr::filter(compartment == "death") %>% 
  ggplo

tmp %>% 
  dplyr::filter(compartment == "death_o") %>% 
  group_by(t) %>% 
  summarise(value = sum(value)) %>% 
  mutate(date = ymd("2020-01-01") + t) %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_line()

res_all %>% 
  bind_rows(.id = "scenario_id") %>% 
  left_join(panel_final %>% 
              rownames_to_column(var = "scenario_id"),
            by = "scenario_id") -> output

output %>% 
  dplyr::filter(compartment == "cases")  %>% 
  group_by(scenario_id, year, compartment, scenario) %>% 
  summarise(incidence = sum(incidence),
            cov_2024 = unique(cov_2024)) %>% 
  ggplot(., aes(x = year, 
                y = incidence, 
                color = cov_2024,
                group = interaction(compartment, scenario, scenario_id))) +
  geom_line() +
  facet_wrap(compartment~scenario, scales = "free")

tmp %>%
  dplyr::filter(compartment %in% compartment_pop) %>%
  mutate(date = ymd("2020-01-01") + t,
         year = year(date),
         compartment_broad = substr(compartment, 1, 1)) %>%
  dplyr::filter(date <= ymd("2026-12-31"), date >= ymd("2021-01-01")) %>%
  dplyr::filter(group == "20-24") %>% 
  ggplot(., aes(x = date, y = value, colour = compartment, fill = compartment, group = compartment)) +
  # geom_line() + geom_point() +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~compartment, scales = "free")
  facet_grid(rows = vars(compartment), scales = "free") #+
  # geom_vline(xintercept = ymd("2023-01-01"), linetype = 2)+
  # geom_vline(xintercept = ymd("2024-12-31"), linetype = 2)

# res_all[[4]] %>% 
#   group_by(year, compartment) %>% 
#   summarise(incidence = sum(incidence)) %>% 
#   pivot_wider(names_from = compartment,
#               values_from = incidence)

  
  # check your roll-out strategy
  # setting_list[[i]]$schedule$primary_course$values %>%
  #   unlist %>%
  #   matrix(., ncol = 16, byrow = T) %>%
  #   data.table %>%
  #   mutate(t =  setting_list[[i]]$schedule$primary_course$times,
  #          date = lubridate::ymd("2020-01-01") + t) %>%
  #   pivot_longer(cols = starts_with("V"),
  #                names_to = "age_group",
  #                values_to = "vaccinated") %>%
  #   mutate(age_group = factor(age_group, levels = paste0("V",1:16))) %>%
  #   ggplot(., aes(x = date, y = vaccinated, group = age_group, color = age_group)) +
  #   geom_point() +
  #   facet_wrap(~age_group, scales = "free") -> p_primary
  # 
  # setting_list[[i]]$schedule$booster$values %>%
  #   unlist %>%
  #   matrix(., ncol = 16, byrow = T) %>%
  #   data.table %>%
  #   mutate(t =  setting_list[[i]]$schedule$booster$times,
  #          date = lubridate::ymd("2020-01-01") + t) %>%
  #   pivot_longer(cols = starts_with("V"),
  #                names_to = "age_group",
  #                values_to = "vaccinated") %>%
  #   mutate(age_group = factor(age_group, levels = paste0("V",1:16))) %>%
  #   ggplot(., aes(x = date, y = vaccinated, group = age_group, color = age_group)) +
  #   geom_line() +
  #   facet_wrap(~age_group, scales = "free") -> p_booster

# i = 1
# tmp <- cm_simulate(setting_list[[i]])$dynamics

# 
# tmp %>% 
#   dplyr::filter(compartment %in% compartment_pop) %>% 
#   mutate(date = ymd("2020-01-01") + t,
#          year = year(date),
#          compartment_broad = substr(compartment, 1, 1)) %>% 
#   dplyr::filter(date <= ymd("2022-10-04"),
#                 date >= ymd("2022-10-01")) %>% 
#   dplyr::select(-t, -year, -compartment_broad) %>% 
#   group_by(compartment, date) %>% summarise(value = sum(value)) %>% 
#   pivot_wider(names_from = date, values_from = value) %>% 
#   dplyr::filter(compartment %in% c("Sv_l", "Sv_m", "Sv_h", "Rv_l", "Rv_m", "Rv_h")) 
# 
# setting_list[[i]]$schedule$primary_course$values %>% 
#   map(t) %>%   
#   map(data.table) %>% 
#   rbindlist() %>% 
#   mutate(t = setting_list[[i]]$schedule$primary_course$times) %>% 
#   # rownames_to_column(var = "t") %>% 
#   melt(., id.vars = "t") %>% 
#   mutate(group = parse_number(as.character(variable)),
#          date = ymd("2020-01-01") + as.numeric(t)) %>% 
#   dplyr::filter(date <= ymd("2022-12-04"),
#                 date >= ymd("2022-05-30")) %>% 
#   ggplot(., aes(x = date, y = value, group = group, color = group)) +
#   geom_point() +
#   facet_wrap(~group) +
#   geom_vline(xintercept = ymd("2023-01-01"), linetype = 2)+
#   geom_vline(xintercept = ymd("2024-12-31"), linetype = 2)
#   
# setting_list[[i]]$schedule$booster$values %>% 
#   map(t) %>%   
#   map(data.table) %>% 
#   rbindlist() %>% 
#   mutate(t = setting_list[[i]]$schedule$booster$times) %>% 
#   # rownames_to_column(var = "t") %>% 
#   melt(., id.vars = "t") %>% 
#   mutate(group = parse_number(as.character(variable)),
#          date = ymd("2020-01-01") + as.numeric(t)) %>% 
#   dplyr::filter(date <= ymd("2022-10-05"),
#                 date >= ymd("2022-09-28")) %>% 
#   dplyr::select(-t) %>% 
#   pivot_wider(names_from = date, values_from = value) %>% View()
#   
# ggplot(., aes(x = t, y = value, group = group, color = group)) +
#   geom_line() +
#   facet_wrap(~group) +
#   geom_vline(xintercept = ymd("2023-01-01"), linetype = 2)+
#   geom_vline(xintercept = ymd("2024-12-31"), linetype = 2)
# 
# setting_list[[i]]$schedule$booster
# 
# res_all[[i]] %>% 
#   dplyr::select(-prop, -cohort_all) %>% 
#   pivot_wider(names_from = compartment,
#               values_from = incidence) %>% 
#   dplyr::filter(critical < death, group_index < 12)
# output %>% 
#   group_by(scenario_id, cov_2024, start_age_annual, start_age_6m, scenario, year, compartment) %>% 
#   summarise(incidence = sum(incidence),
#             cohort_all = sum(cohort_all)) %>% 
#   dplyr::filter(compartment == "death") %>% 
#   ggplot(., aes(x = year, y = incidence/cohort_all, colour = cov_2024, group = cov_2024)) +
#   geom_line() +
#   facet_wrap(~scenario, scales = "free")
