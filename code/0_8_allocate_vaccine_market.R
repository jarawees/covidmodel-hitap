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
         sv_full = as.numeric(NA),
         az_full = as.numeric(NA),
         sp_full = as.numeric(NA),
         pf_full = as.numeric(NA),
         md_full = as.numeric(NA),
         sv_remain = as.numeric(NA),
         az_remain = as.numeric(NA),
         sp_remain = as.numeric(NA),
         pf_remain = as.numeric(NA),
         md_remain = as.numeric(NA)) |> 
  data.frame() -> vaccine_allocation

# vaccine_allocation |> 
#   dplyr::select(date, ends_with("primary")) |> 
#   pivot_longer(ends_with("primary")) |> 
#   ggplot(aes(x = date, y = value, color = name, fill = name)) +
#   geom_bar(position = "stack",
#            stat = "identity")

for(t in 1:(nrow(vaccine_allocation)-28)){
  val_dose2 <- vaccine_allocation |> dplyr::select(ends_with("_dose2")) |> slice(t+28) 
  val_dose1 <- vaccine_allocation |> dplyr::select(ends_with("_dose1")) |> slice(t)

  
  for(i in 1:length(val_dose2)){
    tar_full <- gsub("dose2", "full",colnames(val_dose2)[i])
    tar_remain <- gsub("dose2", "remain",colnames(val_dose2)[i])
    
    if(val_dose1[,i] <= val_dose2[,i]){
      vaccine_allocation[,tar_full][t+28] <- val_dose1[,i]
      vaccine_allocation[,tar_remain][t+28] <- val_dose2[,i] - val_dose1[,i]
      vaccine_allocation[,tar_remain][t] <- 0
    }
    
    if(val_dose1[,i] > val_dose2[,i]){
      vaccine_allocation[,tar_full][t+28] <- val_dose2[,i]
      vaccine_allocation[,tar_remain][t] <- val_dose1[,i] - val_dose2[,i]
      vaccine_allocation[,tar_remain][t+28] <- 0
    }
  }
}
  
  
vaccine_allocation |> slice(30,58)

vaccine_allocation |> 
  ggplot() +
  geom_point(aes(x = date, y = sv_full)) +
  geom_point(aes(x = date, y = sv_remain), color = "red")

