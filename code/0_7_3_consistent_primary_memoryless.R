

vaccine_daily %>% 
  filter(dose_index == "second") %>% 
  dplyr::select(-dose_index) %>% 
  left_join(vp_levels,
            by = "vac_type") %>% 
  mutate(Vl_doses = Vl*value,
         Vm_doses = Vm*value,
         Vh_doses = Vh*value) %>% 
  group_by(date) %>% 
  summarise(Vl_doses = sum(Vl_doses),
            Vm_doses = sum(Vm_doses),
            Vh_doses = sum(Vh_doses)) %>% 
  mutate(all_doses = Vl_doses + Vm_doses + Vh_doses,
         prob_v_p_2l = na_locf(Vl_doses/all_doses),
         prob_v_p_2m = na_locf(Vm_doses/all_doses)) -> scenario3_primary

# scenario3_primary %>% 
#   dplyr::select(-ends_with("doses")) %>%
#   mutate(prob_v_p_2h = 1 - prob_v_p_2l - prob_v_p_2m) %>% 
#   pivot_longer(starts_with("prob")) %>% 
#   ggplot(., aes(x = date, y = value, group = name, color = name)) +
#   geom_line()

vaccine_daily %>% 
  filter(dose_index == "boost") %>% 
  dplyr::select(-dose_index) %>% 
  rename(booster_type = vac_type) %>% 
  left_join(bp_levels,
            by = "booster_type") %>% 
  .[complete.cases(.),] %>% 
  group_by(date) %>% 
  summarise(weighted_Vl2m = sum(Vl2m*value),
            weighted_Vl2h = sum(Vl2h*value),
            unweighted = sum(value)) %>% 
  mutate(prob_v_b_l2m = weighted_Vl2m/unweighted,
         prob_v_b_l2m = na_locf(prob_v_b_l2m)) -> scenario3_booster


### FOR 2024 ONWARDS - additional booster doses ###
# NOTE: a lot of this is hard-coded, need to update

# assume 50% Pfizer and 50% Moderna for future additional booster doses
phase3_prob <- as.numeric(bp_levels[2,"Vl2m"])*0.5 + as.numeric(bp_levels[3,"Vl2m"])*0.5
date_df3 <- as.Date("2022-10-29")

require(lubridate)
df_phase3 <- data.frame(
  date = seq(date_df3, date_end, by = "day")) |>
  mutate(prob_v_b_l2m = phase3_prob)


scenario3_booster <- scenario3_booster |> 
  select(date,prob_v_b_l2m) |> 
  bind_rows(df_phase3)
  

  
  