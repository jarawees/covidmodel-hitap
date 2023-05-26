source("code/0_LoadAll.R")

# params <- cm_parameters_SEI3R("Thailand")
# res <- cm_simulate(params)

out <- read_rds("data/out_20230328.rds")
date_switch <- c("2021-01-15", "2021-07-05", "2021-12-31")

#### Simple example ####
panel <- expand.grid(f = c(1:2), boosting_level = c(0.00001, seq(0.1, 0.9, 0.1)), prioritisation = c("OA only",
                                                                                                     "OA then A",
                                                                                                     "OA and A"))
setting_list <- list()
for(i in 1:nrow(panel)) {
  if (panel$prioritisation[i] == "OA then A") {
    setting_list[[i]] <- parameterise_setting(
      f = panel$f[i],
      prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
      boosting_level = panel$boosting_level[i]
    )
  }
  if (panel$prioritisation[i] == "OA and A") {
    setting_list[[i]] <- parameterise_setting(
      f = panel$f[i],
      prioritisation_followup = c(NA,rep(1,15)),
      boosting_level = panel$boosting_level[i]
    )
  }
  if (panel$prioritisation[i] == "OA only") {
    setting_list[[i]] <- parameterise_setting(
      f = panel$f[i],
      prioritisation_followup = c(rep(NA,12),rep(1,4)),
      boosting_level = panel$boosting_level[i]
    )
  }
}

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

write_rds(res_all, "data/res_all_low_transmissibility.rds")
  # vaccinate_booster_annual(para = .,
  #                          vac_data = owid_vac,
  #                          uptake_by_existing = 0.9,
  #                          prioritisation_initial = c(rep(NA, 4), rep(1,12)),
  #                          prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
  #                          campaign_month = c(10:12,1:2))

# check_vaccination_program(type = "booster",
#                           para = params) -> p_booster
# 
# check_vaccination_program(type = "primary_course",
#                           para = params) -> p_primary


# check schedule objects generated


# res <- cm_simulate(params)
# 
# res$dynamics %>%
#   filter(compartment == "death_o") %>% 
#   mutate(date = ymd(ymd("2020-10-01") + out$optim$bestmem[2]) + t,
#          year = year(date)) %>% 
#   group_by(year) %>% 
#   summarise(deaths = sum(value))
# 
# res$dynamics %>% 
#   filter(compartment == "death_o") %>% 
#   ggplot(., aes(x = t, y = value)) +
#   geom_line() +
#   facet_wrap(~group)

# res$dynamics |> 
#   filter(t == 458) |> 
#   filter(grepl("v_l",compartment))

# res$dynamics |> filter(compartment == "Sv_l") |> filter(value > 0) |> View()
# all_labels <- unique(res$dynamics$compartment)
# compartment_labels <- all_labels[1:24]
# outcome_labels <-  all_labels[25:51]

# testing
# (1) population level unchanged

# res$dynamics |> 
#   filter(compartment %in% compartment_labels) |> 
#   group_by(t) |> 
#   summarise(tot = sum(value)) |> 
#   mutate(diff = c(0,diff(tot))) |> 
#   filter(diff > 0.1) |> 
#   nrow() -> tmp
# 
# if(tmp != 0) print("population unbalanced!")

# (2) people moving to vaccinated stages
# res$dynamics |>
#   filter(compartment %in% c("Rv_l", "Rv_m", "R", "Rv_h", "S", "Sv_l", "Sv_m", "Sv_h")) |>
#   # filter(compartment %in% compartment_labels) |>
#   group_by(t, compartment) |> summarise(value = sum(value)) |>
#   mutate(date = ymd("2020-01-01") + t) |>
#   ggplot(aes(x = date, y = value, group = compartment, color = compartment, fill = compartment)) +
#   # geom_bar(position = "stack", stat = "identity")
#   geom_point() +
#   facet_wrap(~compartment)
# 
# res$dynamics |>
#   filter(compartment %in% c("Sv_l", "Rvl")) |>
#   group_by(t) |> summarise(value = sum(value)) |> 
#   ggplot(aes(x = t, y = value)) +geom_point()
# 
# res$dynamics |>
#   filter(compartment %in% c("Sv_l")) |> 
#   ggplot(aes(x = t, y = value)) +
#   geom_point() +
#   facet_wrap(~group)

# res$dynamics |> 
#   filter( grepl("Sv|Rv|R|S",compartment)) |> 
#   pivot_wider(names_from = compartment,
#               values_from = value) |> 
#   filter(group == "60-64") |> View()
# 
# 
# # everything fall on day 501
# res$dynamics |>
#   filter(compartment %in% c("Rv_l", "Rv_m", "R", "Rv_h", "S", "Sv_l", "Sv_m", "Sv_h")) |>
#   # group_by(t, compartment) |> summarise(value = sum(value))  |> 
#   pivot_wider(names_from = compartment,
#               values_from = value) |> 
#   mutate(date = ymd("2020-01-01") + t) -> tmp
# 
# res$dynamics |> 
#   filter(compartment %in% compartment_labels) |> 
#   filter(t %in% c(501, 500)) |> 
#   pivot_wider(names_from = t, values_from = value) |> 
#   mutate(diff = `501` - `500`) |> 
#   filter(group == "20-24") |> View()
# 
#   group_by(t, compartment, group) |> summarise(value = sum(value)) |> 
#   pivot_wider(names_from = compartment, 
#               values_from = value) |> 
#   mutate(pde = S + E + Ip + Ia + R,
#          p_S = S/pde,
#          p_R = R/pde) |> 
#   filter(t == 501) |> View()
# 
# 
# 
# res$dynamics$
# 
# res$dynamics  |> 
#   filter(compartment %in% c("S", "Sv_l","Sv_m", "Sv_h","R","Rv_l","Rv_m","Rv_h")) |> 
#   # group_by(t, compartment) |> 
#   # summarise(value = sum(value)) |>
#   pivot_wider(names_from = compartment,
#               values_from = value) |> 
#   mutate(date = ymd("2020-01-01") + t) -> tmp
# 
# tmp |> 
#   ungroup() |> 
#   mutate(all_s_vaxxed = Sv_l + Sv_m + Sv_h,
#          all_r_vaxxed = Rv_l + Rv_m + Rv_h,
#          tot = all_s_vaxxed + all_r_vaxxed) |> 
#   filter(t >= 500) |> View()
#   ggplot(aes(x = date, y = tot_diff)) +
#   geom_point() +
#   geom_hline(aes(yintercept = doses))
# 
# 
# # |> 
#   #filter(compartment %in% compartment_labels) |> 
#   group_by(t) |> 
#   summarise(tot = sum(value)) |> 
#   mutate(diff = tot - 69799978) |> 
#   filter(t > 500)
# 

# # 
# res$dynamics |>
#   # filter(compartment == "Rv_m") |>
#   # filter(grepl("_m", compartment)) |>
#   filter(compartment %in% compartment_labels) |>
#   group_by(t, compartment) |> summarise(value = sum(value)) |>
#   mutate(date = ymd("2020-01-01") + t)  |>
#   ggplot(aes(x = date, y = value)) +
#   geom_point() +
#   facet_wrap(~compartment, scales = "free")
#   # filter(value > 0, t == 502) |>
#   # pull(value) |> sum()
# 
# res$dynamics |> 
#   filter(grepl("_m", compartment)) |> 
#   filter(compartment %in% compartment_labels) |> 
#   pivot_wider(names_from = compartment, values_from = value) |> 
#   filter(t > 500)
