# process_parameters <- data.frame(source = unlist(lapply(para$processes, "[[", "source"))) 
# process_parameters[["type"]] <- unlist(lapply(para$processes, "[[", "type"))
# process_parameters[["names"]] <- unlist(lapply(para$processes, "[[", "names") %>% 
#                                           map(paste, collapse = ", "))
# process_parameters[["report"]] <- unlist(lapply(para$processes, "[[", "report") %>% 
#                                           map(paste, collapse = ", "))
# 
# process_parameters
# 
# para$processes[[1]]$prob
# para$processes[[5]]$prob

para <- gen_country_basics(country = "Thailand",
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
                           prob_v_b_l2m = 0,
                           deterministic = TRUE,
                           scenario_primary = scenario3_primary,
                           scenario_booster = scenario3_booster %>% mutate(prob_v_b_l2m = if_else(date >= "2024-01-01", 0, prob_v_b_l2m)),
                           seed = out$optim$bestmem[2]) %>% 
  update_u_y(para = .,
             date_switch = c("2021-01-15", "2021-07-05", "2021-12-31", "2025-01-01"),
             rc_u = c(1, 1.5, 1, 1), # relative changes in u
             rc_y = c(1, 1, 1, 1), # relative changes in y
             rc_ve = c(1, 0.9, 0.9, 1/0.81), # update rc_ve using meta-analysis from Dec2023
             efficacy_baseline = efficacy_all
  ) %>%
  emerge_VOC_burden(para = .,
                    rc_severity = c(1, 1.5, 0.7, 0.7), # relative change in ihr and ifr
                    efficacy_baseline = efficacy_all) %>%
  vaccinate_primary(para = .,
                    vac_data = owid_vac,
                    values = primary_allocation_plan) %>%
  vaccinate_additional(para = .,
                       vac_data = owid_vac,
                       booster_plan = booster_allocation_plan,
                       start_age_annual = 40,
                       start_age_6m = 80,
                       cov_2024 = 0.9,
                       month_annual = c(5:6),
                       month_6m = c(11:12))

res <- cm_simulate(para)$dynamics %>% 
  mutate(date = t + ymd("2021-02-01"),
         year = year(date))

res_by_compartment <- res %>% 
  dplyr::filter(compartment %in% compartments_status) %>% 
  mutate(compartment_general = substr(compartment, 1, 1))

res_by_compartment %>% 
  dplyr::filter(compartment == "Ev_l",
                date <= "2023-01-01") %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_point() +
  facet_wrap(~group)

res_by_outcome <- res %>% 
  filter(grepl("case|sever|critical|death", compartment)) %>% 
  filter(!grepl("_p|reported", compartment)) %>% 
  pivot_wider(names_from = compartment, values_from = value) %>% 
  mutate(severe_all = case_when(date <= date_switch[1] ~ severe_i,
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
                               date > date_switch[4] ~ death_voc4_o)) %>% 
  dplyr::select(date, year, group, cases, ends_with("all"))

res_by_outcome %>% 
  pivot_longer(cols = c("cases", "severe_all", "critical_all", "death_all")) %>% 
  group_by(year, name) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = name,
              values_from = value)

res_by_outcome %>% 
  dplyr::filter(group == "75+") %>% 
  ggplot(., aes(x = date)) +
  geom_line(aes(y = severe_all), color = "yellow")+
  geom_line(aes(y = critical_all), color = "orange")+
  geom_line(aes(y = death_all), color = "red")

res_by_compartment %>% 
  dplyr::filter(group == "75+") %>% 
  ggplot(., aes(x = date, y = value, color = compartment, fill = compartment)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~compartment, scales = "free")

res_by_outcome %>% 
  pivot_longer(cols = c("cases", "severe_all", "critical_all", "death_all")) %>% 
  dplyr::filter(group == "75+") %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~name, scales = "free")



# res_by_compartment %>% 
#   group_by(date) %>% 
#   summarise(value = sum(value)) %>% 
#   ggplot(., aes(x = date, y = value)) + 
#   geom_line()

# para$schedule$booster$values %>% bind_cols() -> tmp_b
# para$schedule$primary_course$values %>% bind_cols() -> tmp_p
# 
# # paste0("f = ", panel[i,1], 
# #        "; boosting level = ", panel[i,2],
# #        "; prioritisation = ", panel[i,3]) -> tmp_title
# 
# tmp_b %>% 
#   t %>% 
#   data.table %>% 
#   bind_cols(date_grid1) %>% 
#   right_join(date_grid2, by = c("date", "t_within")) %>% 
#   arrange(date) %>% 
#   mutate_at(vars(starts_with("V")),
#             na_locf) %>% 
#   pivot_longer(cols = starts_with("V"),
#                names_to = "age_group") %>% 
#   mutate(age_group = parse_number(age_group)) %>% 
#   left_join(pop_TH,
#             by = "age_group") %>% 
#   mutate(daily_cov = value/pop_age,
#          age_group_broad = case_when(age_group == 1 ~ "infants/ toddlers",
#                                      age_group %in% 2 ~ "children",
#                                      age_group %in% c(3,4) ~ "adolescents",
#                                      age_group %in% c(5:12) ~ "adults",
#                                      age_group %in% c(13:16) ~ "older adults"),
#          age_group_broad = factor(age_group_broad,
#                                   levels = c("infants/ toddlers",
#                                              "children",
#                                              "adolescents",
#                                              "adults",
#                                              "older adults"))) %>% 
#   ggplot(., aes(x = date, y = daily_cov, group = age_group, color = age_group_broad)) +
#   geom_line() +
#   geom_vline(xintercept = seq(ymd("2021-01-01"),
#                               ymd("2031-01-01"),
#                               by = "year"),
#              linetype = 2) + 
#   facet_wrap(~age_group)
