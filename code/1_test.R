# if(!require(pacman)) install.packages("pacman")
# library(pacman)
# p_load(tidyverse)

# require(tidyverse)

source("code/0_LoadAll.R")

compartment_labels <- c("S", "Sv_l", "Sv_m",
                        "E", "Ev_l", "Ev_m",
                        "Ip", "Ip_l", "Ip_m",
                        "Is", "Is_l", "Is_m",
                        "Ia", "Ia_l", "Ia_m",
                        "R", "Rv_l", "Rv_m")
outcome_labels <- c("cases", "severe_p", "critical_p", "severe_i", "critical_i", "death_o")

#### Simple example ####
params <- gen_country_basics(country = "Thailand",
                             R0_assumed = 2.7,
                             date_start = "2020-01-01",
                             date_end = "2022-12-31",
                             contact = contact_schedule,
                             period_wn = 3*365,
                             period_wv_ml = 1*365,
                             processes = burden_processes_az,
                             prob_v_p_l = 0.999999) %>%
  update_u_y(para = .,
                  date_switch = c("2021-01-15", "2021-04-15", "2021-12-15"),
                  rc_u = c(1, 1.5, 0.5), # relative changes in u
                  rc_y = c(1, 0.5, 0.5), # relative changes in y
                  rc_ve = c(1, 0.9, 0.7),
                  efficacy_baseline = ve_az,
                  efficacy_weights = efficacy_weights_test
                  ) %>%
  emerge_VOC_burden(para = .,
    rc_severity = c(1, 1.5,1.5), # relative change in ihr and ifr
    efficacy_baseline = ve_az) %>%
  vaccinate_primary(para = .) # %>%
  # vaccinate_booster(para = .,
  #                   program_start = "2021-12-31",
  #                   program_end = "2022-06-30")
  
res <- cm_simulate(params) 

# plot(x = params$schedule$yv_l_scaler$times,
#      y = params$schedule$yv_l_scaler$values[[1]],
#      type = "l")

# params$pop[[1]]$ur <- rep(0, 16)

# print(unique(res$dynamics$compartment))
# res$dynamics |>
#   # filter(!compartment %in% c("cases", "cases_reported", "foi","foiv_l", "foiv_m", "subclinical")) |>
#   filter(compartment %in% compartment_labels) |>
#   group_by(t, compartment) |> summarise(value = sum(value)) |>
#   # group_by(t) |>
#   # summarise(value = sum(value)) |>
#   # ggplot(aes(x = t, y = value)) + geom_line()
#   ggplot(aes(x = t, y = value, group = compartment, color = compartment, fill = compartment)) +
#   # geom_line() +
#   geom_bar(position = "stack", stat = "identity")
# # # 
# res$dynamics |>
#   # filter(compartment %in% compartment_labels) |>
#   filter(compartment %in% c("Rv_l", "Rv_m", "R", "S", "Sv_l", "Sv_m")) |>
#   # filter(!compartment %in% c("cases", "cases_reported", "foi","foiv_l", "foiv_m", "subclinical")) |>
#   # filter(compartment %in% outcome_labels[c(1, 4:6)])|>
#   # filter(grepl("_m", compartment)) |> 
#   group_by(t, compartment) |> summarise(value = sum(value)) |>
#   mutate(date = ymd("2020-01-01") + t) |>
#   ggplot(aes(x = date, y = value, group = compartment, color = compartment, fill = compartment)) +
#   # geom_bar(position = "stack", stat = "identity")
#   geom_point() +
#   facet_wrap(~compartment, scales = "free", ncol = 3)
# 
# 
res$dynamics |>
  # filter(compartment == "Rv_m") |>
  # filter(grepl("_m", compartment)) |>
  filter(compartment %in% compartment_labels) |>
  group_by(t, compartment) |> summarise(value = sum(value)) |>
  mutate(date = ymd("2020-01-01") + t)  |>
  ggplot(aes(x = date, y = value)) +
  geom_point() +
  facet_wrap(~compartment, scales = "free")
  # filter(value > 0, t == 502) |>
  # pull(value) |> sum()

res$dynamics |> 
  filter(grepl("_m", compartment)) |> 
  filter(compartment %in% compartment_labels) |> 
  pivot_wider(names_from = compartment, values_from = value) |> 
  filter(t > 500)
