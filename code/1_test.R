source("code/0_LoadAll.R")

# params <- cm_parameters_SEI3R("Thailand")
# res <- cm_simulate(params)

#### Simple example ####
params <- gen_country_basics(country = "Thailand",
                             R0_assumed = 2.7,
                             date_start = "2020-01-01",
                             date_end = "2030-12-31",
                             contact = contact_schedule,
                             processes = gen_burden_processes(VE = efficacy_all),
                             period_wn  = 3*365, # duration, waning of natural immunity
                             period_wv_m2l = 1*365, # duration, waning from medium to low levels vaccine induced 
                             period_wv_h2m = 1*365, # duration, waning from medium to low levels vaccine induced 
                             prob_v_p_2l = 1,
                             prob_v_p_2m = 0,
                             prob_v_b_l2m = 0.5,
                             deterministic = TRUE) %>% 
  update_u_y(para = .,
             date_switch = c("2021-01-15", "2021-04-15", "2021-12-15"),
             rc_u = c(1, 1.5, 0.5), # relative changes in u
             rc_y = c(1, 0.5, 0.5), # relative changes in y
             rc_ve = c(1, 0.9, 0.7), # relative evasiveness 
             efficacy_baseline = efficacy_all,
             efficacy_weights = efficacy_weights_test
  ) %>%
  emerge_VOC_burden(para = .,
    rc_severity = c(1, 1.5,1.5), # relative change in ihr and ifr
    efficacy_baseline = efficacy_all) %>%
  vaccinate_primary(para = .,
                    vac_data = owid_vac,
                    values = primary_allocation_plan) %>%
  vaccinate_booster_annual(para = .,
                           vac_data = owid_vac,
                           uptake_by_existing = 0.9,
                           prioritisation_initial = c(rep(NA, 4), rep(1,12)),
                           prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
                           campaign_month = c(10:12,1:2))

check_vaccination_program(type = "booster",
                          para = params) -> p_booster

check_vaccination_program(type = "primary_course",
                          para = params) -> p_primary


# check schedule objects generated
params$schedule$primary_course$values |> 
  map(data.frame) |> map(t) |> map(data.frame) |> bind_rows() |> 
  mutate(rn = 1:n()) |> set_rownames(NULL) |> 
  pivot_longer(starts_with("X")) |> 
  mutate(name = factor(name, levels = paste0("X",1:16))) |> 
  ggplot(aes(x = rn, y = value, group = name)) +
  geom_point() +
  facet_wrap(~name)

params$schedule$booster$times

params$schedule$booster$values |> 
  map(data.frame) |> map(t) |> map(data.frame) |> bind_rows() |> 
  mutate(rn = 1:n()) |> set_rownames(NULL) |> 
  pivot_longer(starts_with("X")) |> 
  mutate(name = factor(name, levels = paste0("X",1:16))) |> 
  ggplot(aes(x = rn, y = value, group = name)) +
  geom_point() +
  facet_wrap(~name)

res <- cm_simulate(params)
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
res$dynamics |>
  filter(compartment %in% c("Rv_l", "Rv_m", "R", "Rv_h", "S", "Sv_l", "Sv_m", "Sv_h")) |>
  # filter(compartment %in% compartment_labels) |>
  group_by(t, compartment) |> summarise(value = sum(value)) |>
  mutate(date = ymd("2020-01-01") + t) |>
  ggplot(aes(x = date, y = value, group = compartment, color = compartment, fill = compartment)) +
  # geom_bar(position = "stack", stat = "identity")
  geom_point() +
  facet_wrap(~compartment)

res$dynamics |>
  filter(compartment %in% c("Sv_l", "Rvl")) |>
  group_by(t) |> summarise(value = sum(value)) |> 
  ggplot(aes(x = t, y = value)) +geom_point()

res$dynamics |>
  filter(compartment %in% c("Sv_l")) |> 
  ggplot(aes(x = t, y = value)) +
  geom_point() +
  facet_wrap(~group)

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
