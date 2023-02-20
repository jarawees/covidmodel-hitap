owid_vac |> 
  dplyr::select(people_vaccinated, 
                people_fully_vaccinated,
                total_boosters,
                date) |> 
  mutate(dose1_daily = c(0, diff(people_vaccinated)),
         dose2_daily = c(0, diff(people_fully_vaccinated)),
         boost_daily = c(0,diff(total_boosters)),
         date = ymd(date)) |> 
  left_join(seg_primary |> 
              mutate(vac_type = paste0(vac_type,"_primary")) |> 
              pivot_wider(names_from = vac_type,
                          values_from = value) |> 
              dplyr::select(-vac_type2),
            by = "date") |> 
  mutate_at(vars(ends_with("_primary")), imputeTS::na_locf, na_remaining = "rev") |> 
  mutate(sv_dose1 = dose1_daily*sv_primary,
         az_dose1 = dose1_daily*az_primary,
         sp_dose1 = dose1_daily*sp_primary,
         pf_dose1 = dose1_daily*pf_primary,
         md_dose1 = dose1_daily*moderna_primary,
         sv_dose2 = dose2_daily*sv_primary,
         az_dose2 = dose2_daily*az_primary,
         sp_dose2 = dose2_daily*sp_primary,
         pf_dose2 = dose2_daily*pf_primary,
         md_dose2 = dose2_daily*moderna_primary,
         sv_vac_dose1 = NA,
         sv_vac_dose1_parked = NA,
         sv_vac_full = NA) |> 
  data.frame() -> vaccine_allocation

# vaccine_allocation |> 
#   dplyr::select(date, ends_with("primary")) |> 
#   pivot_longer(ends_with("primary")) |> 
#   ggplot(aes(x = date, y = value, color = name, fill = name)) +
#   geom_bar(position = "stack",
#            stat = "identity")

for(t in 1:(nrow(vaccine_allocation))){
  if(t == 1){
    vaccine_allocation$sv_vac_dose1_parked[t] <- vaccine_allocation$sv_dose1[t]
    vaccine_allocation$sv_vac_dose1[t] <- 0
    vaccine_allocation$sv_vac_full[t] <- 0
  } 
  
  if(t > 1 & t <= 21){
    vaccine_allocation$sv_vac_dose1_parked[t] <- vaccine_allocation$sv_vac_dose1_parked[t-1] + vaccine_allocation$sv_dose1[t]
    vaccine_allocation$sv_vac_dose1[t] <- 0
    vaccine_allocation$sv_vac_full[t] <-  0
  } 
  
  if(t > 21){
    vaccine_allocation$sv_vac_dose1_parked[t] <- vaccine_allocation$sv_vac_dose1_parked[t-1] + vaccine_allocation$sv_dose1[t] - vaccine_allocation$sv_dose1[t-21]
    vaccine_allocation$sv_vac_dose1[t] <- vaccine_allocation$sv_vac_dose1[t-1] + vaccine_allocation$sv_dose1[t-21] - vaccine_allocation$sv_dose2[t]
    vaccine_allocation$sv_vac_full[t] <-  vaccine_allocation$sv_vac_full[t-1] + vaccine_allocation$sv_dose2[t]
  }
}

vaccine_allocation  %>% 
  ggplot(., aes(x = date)) +
  geom_line(aes(y = sv_vac_dose1_parked))+
  geom_line(aes(y = sv_vac_dose1), color = "red") +
  geom_line(aes(y = sv_vac_full), color = "blue")
  