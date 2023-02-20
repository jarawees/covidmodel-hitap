doc[[2]] %>% 
  dplyr::select(date, ends_with("second")) %>% 
  .[,1:7] %>% 
  dplyr::filter(!is.na(total_second)) %>% 
  mutate(date = lubridate::ymd(date),
         check_second = sv_second + az_second + sp_second + pf_second + moderna_second) %>% 
  arrange(date) %>% 
  filter(!is.na(check_second), check_second != 0) %>% 
  mutate(sv_second_prop = sv_second/check_second,
         az_second_prop = az_second/check_second,
         sp_second_prop = sp_second/check_second,
         pf_second_prop = pf_second/check_second,
         moderna_second_prop = moderna_second/check_second) %>% 
  add_row(date = ymd("2021-02-28")) %>%
  replace(., is.na(.), 0) %>% 
  mutate(sv_second_prop = if_else(date == "2021-02-28",
                                  1,
                                  sv_second_prop),
         observed = T) %>%
  right_join(data.frame(date = seq(ymd("2021-02-28"),
                                   ymd("2022-10-28"),
                                   "day")),
             by = "date") %>% 
  arrange(date) %>% 
  dplyr::select(-sv_second, -az_second, -sp_second, -pf_second, -moderna_second, -check_second, -total_second) %>% 
  mutate(observed = if_else(is.na(observed), F, observed)) %>% 
  mutate_at(vars(ends_with("prop")), function(x) imputeTS::na_interpolation(x, option = "spline")) %>% 
  mutate(sv_second_prop = if_else(sv_second_prop > 1, 1, sv_second_prop),
         az_second_prop = if_else(az_second_prop > 1, 1, az_second_prop),
         sp_second_prop = if_else(sp_second_prop > 1, 1, sp_second_prop),
         pf_second_prop = if_else(pf_second_prop > 1, 1, pf_second_prop),
         moderna_second_prop = if_else(moderna_second_prop > 1, 1, moderna_second_prop),
         sv_second_prop = if_else(sv_second_prop < 0, 0, sv_second_prop),
         az_second_prop = if_else(az_second_prop < 0, 0, az_second_prop),
         sp_second_prop = if_else(sp_second_prop < 0, 0, sp_second_prop),
         pf_second_prop = if_else(pf_second_prop < 0, 0, pf_second_prop),
         moderna_second_prop = if_else(moderna_second_prop < 0, 0, moderna_second_prop)) %>% 
  left_join(doc[[2]] %>% 
              dplyr::select(date, total_second) %>% 
              mutate(date = ymd(date)),
            by = "date") %>% 
  mutate(total_second = imputeTS::na_interpolation(total_second, option = "stine"),
         total_second = if_else(total_second < 0, 0, total_second),
         second_daily = c(0, diff(total_second)),
         sv_full = sv_second_prop*second_daily,
         az_full = az_second_prop*second_daily,
         sp_full = sp_second_prop*second_daily,
         pf_full = pf_second_prop*second_daily,
         md_full = moderna_second_prop*second_daily) %>% 
  dplyr::select(-ends_with("prop"), -total_second, -second_daily, -observed) %>% 
  setNames(c("date", "2_sv_sv", "2_az_az","2_sp_sp","2_pf_pf","2_md_md")) %>% 
  mutate_at(vars(starts_with("2_")), cumsum) -> vaccinate_primary_bybrand

combo %>% 
  filter(applicability == "Y",
         dose1 == dose2) %>% 
  mutate(label = paste0("3_", dose1,"_",dose2,"_",dose3)) %>% 
  pull(label) -> new_cols
vaccinate_primary_bybrand[,new_cols] <- as.numeric(NA)
vaccinate_primary_bybrand[1,6:17] <- 0

for(d in 2:nrow(vaccinate_primary_bybrand)){
  date_tmp <- vaccinate_primary_bybrand$date[d]
  dose3_tmp <- list()
  for(v in c("az", "md", "pf")){
    vac_type_tmp <- v
    labels_target <- data.frame(
      new = combo %>%
        filter(applicability == "Y",
               dose1 == dose2,
               dose3 == v) %>%
        mutate(label = paste0("3_", dose1, "_", dose2, "_", dose3)) %>%
        pull(label)
    ) %>%
      distinct() %>%
      mutate(
        previous = substr(new, 1, 7),
        previous = gsub("3_", "2_", previous),
        previous_r = paste0(previous, "_r")
      )
    
    vaccinate_primary_bybrand[vaccinate_primary_bybrand$date == date_tmp-1, 
                              labels_target$previous] %>% mutate(tot = sum(.)) -> dose3_tmp[[vac_type_tmp]]
    
    for (m in 1:nrow(labels_target)){
      dose3_tmp[[vac_type_tmp]][, labels_target$previous_r[m]] <-
        if_else(dose3_tmp[[vac_type_tmp]]$tot == 0,
                0,
                (dose3_tmp[[vac_type_tmp]] %>% pull(labels_target$previous[m]) %>% unlist()) /
                  dose3_tmp[[vac_type_tmp]]$tot) 
    }
    
    for (m in 1:nrow(labels_target)){
      dose3_tmp[[vac_type_tmp]][, labels_target$new[m]] <- ((vaccine_daily %>% 
                                                               filter(date %in% date_tmp, dose_index == "boost", vac_type == vac_type_tmp) %>% 
                                                               pull(value)) * dose3_tmp[[vac_type_tmp]]) %>% 
        pull(labels_target$previous_r[m]) %>% unlist()
    }
    vaccinate_primary_bybrand[vaccinate_primary_bybrand$date == date_tmp, labels_target$new] <-   vaccinate_primary_bybrand[vaccinate_primary_bybrand$date == date_tmp - 1, labels_target$new] + dose3_tmp[[vac_type_tmp]][,labels_target$new]
    vaccinate_primary_bybrand[vaccinate_primary_bybrand$date == date_tmp, labels_target$previous] <- vaccinate_primary_bybrand[vaccinate_primary_bybrand$date == date_tmp, labels_target$previous] -  dose3_tmp[[vac_type_tmp]][,labels_target$new] 
  }
}  

vaccinate_primary_bybrand %>% 
  tail(1) %>% 
  pivot_longer(colnames(.)[2:17]) %>% 
  pull(value) %>% sum()

# doc[[2]] %>% 
#   dplyr::select(date, ends_with("boost")) %>% 
#   dplyr::filter(total_boost >= 0) %>% 
#   replace(., is.na(.), 0) %>% 
#   mutate(date = ymd(date),
#          check_boost = sv_boost + az_boost + sp_boost + pf_boost + moderna_boost,
#          sv_boost_prop = sv_boost/check_boost,
#          az_boost_prop = az_boost/check_boost,
#          sp_boost_prop = sp_boost/check_boost,
#          pf_boost_prop = pf_boost/check_boost,
#          md_boost_prop = moderna_boost/check_boost) %>% 
#   right_join(data.frame(date = seq(min(.$date),
#                                    max(.$date),
#                                    "day")),
#              by = "date") %>% 
#   arrange(date) %>% 
#   dplyr::select(-sv_boost, 
#                 -sp_boost, 
#                 -az_boost, 
#                 -pf_boost,
#                 -check_boost,
#                 -moderna_boost) %>% 
#   mutate_at(vars(ends_with("prop")), function(x) imputeTS::na_interpolation(x, option = "spline")) %>% 
#   mutate_at(vars(ends_with("prop")), function(x) replace(x, x < 0, 0)) %>% 
#   mutate_at(vars(ends_with("prop")), function(x) replace(x, x > 1, 1)) %>% 
#   mutate(total_boost = imputeTS::na_interpolation(total_boost, option = "stine"),
#          boost_daily = c(0, diff(total_boost)),
#          sv_boost = sv_boost_prop*boost_daily,
#          az_boost = az_boost_prop*boost_daily,
#          sp_boost = sp_boost_prop*boost_daily,
#          pf_boost = pf_boost_prop*boost_daily,
#          md_boost = md_boost_prop*boost_daily) %>% 
#   dplyr::select(-ends_with("prop")) -> vaccinate_boost_bybrand

# vaccinate_primary_bybrand[,c(1,5:9)] -> res_identical
# vaccinate_primary_bybrand %>% 
#   dplyr::select(-total_second, -second_daily) %>% 
#   pivot_longer(ends_with("full")) %>% 
#   pull(value) %>% sum()
# 
# vaccinate_primary_full_bybrand %>%
#   dplyr::select(-total_second, -second_daily) %>% 
#   #ggplot(., aes(x = date, y = second_daily)) +geom_point()
#   pivot_longer(ends_with("full")) %>% 
#   ggplot(., aes(x = date, y = value, color = name, fill = name)) +
#   geom_point() +
#   labs(y = "Market Share") +
#   facet_wrap(~name)