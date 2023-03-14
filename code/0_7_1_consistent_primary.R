
dose2 %>% 
  dplyr::select(date, ends_with("daily")) %>% 
  setNames(c("date", "2_sv_sv", "2_az_az","2_sp_sp","2_pf_pf","2_md_md")) -> dose2_tmp

combo %>% 
  filter(applicability == "Y",
         dose1 == dose2) %>% 
  mutate(label = paste0("3_", dose1,"_",dose2,"_",dose3)) %>% 
  pull(label) -> new_cols
tmp <- data.frame(date = dose2_tmp$date)
tmp[,c(colnames(dose2_tmp)[2:ncol(dose2_tmp)], new_cols)] <- as.numeric(NA)
tmp[1,2:17] <- 0

for(d in 2:nrow(tmp)){
# for(d in 2:156){
  date_tmp <- tmp$date[d]
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
    
    tmp[tmp$date == date_tmp-1, labels_target$previous] %>% mutate(tot = sum(.)) -> dose3_tmp[[vac_type_tmp]]
    
    for (m in 1:nrow(labels_target)){
      dose3_tmp[[vac_type_tmp]][, labels_target$previous_r[m]] <-
        if_else(dose3_tmp[[vac_type_tmp]]$tot == 0,
                0,
                (dose3_tmp[[vac_type_tmp]] %>% pull(labels_target$previous[m]) %>% unlist()) /
                  dose3_tmp[[vac_type_tmp]]$tot) 
    }
    
    for (m in 1:nrow(labels_target)){
      dose3_tmp[[vac_type_tmp]][, labels_target$new[m]] <- ((vaccine_daily %>% 
                                                               filter(date == date_tmp, dose_index == "boost", vac_type == vac_type_tmp) %>% 
                                                               pull(value)) * dose3_tmp[[vac_type_tmp]]) %>% 
        pull(labels_target$previous_r[m]) %>% unlist()
    }
    
    tmp[tmp$date == date_tmp, labels_target$new] <-
      tmp[tmp$date == date_tmp - 1, labels_target$new] + 
      dose3_tmp[[vac_type_tmp]][, labels_target$new]
    
    tmp[tmp$date == date_tmp, labels_target$previous] <-
      tmp[tmp$date == date_tmp - 1, labels_target$previous] +
      dose2_tmp[dose2_tmp$date == date_tmp, labels_target$previous] -  
      dose3_tmp[[vac_type_tmp]][, labels_target$new] 
  }
}  

# write_rds(tmp, paste0(data_path,"vaccine_market_identical_results.rds"))

tmp %>%
  pivot_longer(colnames(.)[2:17]) %>%
  ggplot(., aes(x = date, y = value, group = name)) +
  geom_line() +
  facet_wrap(~name, scales = "free")

# tmp %>%
#   mutate(tot = rowSums(.[,2:17])) %>% 
#   left_join(vaccine_daily %>% 
#               filter(dose_index == "second") %>% 
#               group_by(date) %>% 
#               summarise(dose2 = sum(value)) %>% 
#               mutate(dose2_cs = cumsum(dose2)),
#             by = "date") %>% 
#   left_join(vaccine_daily %>% 
#               filter(dose_index == "boost") %>% 
#               group_by(date) %>% 
#               summarise(dose3 = sum(value)) %>% 
#               mutate(dose3_cs = cumsum(dose3)),
#             by = "date") %>% 
#   mutate(doses_sum = dose3_cs + dose2_cs) %>% 
#   filter(round(tot) != round(doses_sum)) %>% 
#   mutate(doses_sum - tot)


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