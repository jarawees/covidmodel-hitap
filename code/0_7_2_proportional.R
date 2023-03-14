
vaccine_daily[,1] %>% distinct() -> tmp
tmp[,vaccine_compartments] <- as.numeric(NA)
tmp[1,2:ncol(tmp)] <- 0

# for(d in 2:nrow(tmp)){
  for(d in 2:30){
  date_tmp <- tmp$date[d]
  tmp[tmp$date == date_tmp, "1_sv"] <- tmp[tmp$date == date_tmp - 1, "1_sv"] + vaccine_daily %>% dplyr::filter(date %in% date_tmp, dose_index == "first", vac_type == "sv") %>% pull(value)
  tmp[tmp$date == date_tmp, "1_az"] <- tmp[tmp$date == date_tmp - 1, "1_az"] + vaccine_daily %>% dplyr::filter(date %in% date_tmp, dose_index == "first", vac_type == "az") %>% pull(value)
  tmp[tmp$date == date_tmp, "1_sp"] <- tmp[tmp$date == date_tmp - 1, "1_sp"] + vaccine_daily %>% dplyr::filter(date %in% date_tmp, dose_index == "first", vac_type == "sp") %>% pull(value)
  tmp[tmp$date == date_tmp, "1_pf"] <- tmp[tmp$date == date_tmp - 1, "1_pf"] + vaccine_daily %>% dplyr::filter(date %in% date_tmp, dose_index == "first", vac_type == "pf") %>% pull(value)
  tmp[tmp$date == date_tmp, "1_md"] <- tmp[tmp$date == date_tmp - 1, "1_md"] + vaccine_daily %>% dplyr::filter(date %in% date_tmp, dose_index == "first", vac_type == "md") %>% pull(value)
  
  # allocate second doses
  dose2_tmp <- list()
  for(v in vac_type_list){
    vac_type_tmp <- v
    labels_target <- data.frame(
      new = combo %>% 
        dplyr::select(-combo_index) %>% 
        filter(applicability == "Y", dose2 == vac_type_tmp) %>% 
        unique() %>% 
        mutate(labels = paste0("2_",dose1,"_",dose2)) %>% 
        pull(labels)
    ) %>% 
      distinct() %>% 
      mutate(previous = substr(new,1,4),
             previous = gsub("2_","1_",previous),
             previous_r = paste0(previous, "_r"))
    
    tmp[tmp$date == date_tmp-1, labels_target$previous] %>% mutate(tot = sum(.)) -> dose2_tmp[[vac_type_tmp]]
    
    for (m in 1:nrow(labels_target)){
      dose2_tmp[[vac_type_tmp]][, labels_target$previous_r[m]] <-
        if_else(dose2_tmp[[vac_type_tmp]]$tot == 0,
                0,
                (dose2_tmp[[vac_type_tmp]] %>% pull(labels_target$previous[m]) %>% unlist()) /
                  dose2_tmp[[vac_type_tmp]]$tot) 
    }
    
    for (m in 1:nrow(labels_target)){
      dose2_tmp[[vac_type_tmp]][, labels_target$new[m]] <- ((vaccine_daily %>% 
                                                               filter(date %in% date_tmp, dose_index == "second", vac_type == vac_type_tmp) %>% 
                                                               pull(value)) * dose2_tmp[[vac_type_tmp]]) %>% 
        pull(labels_target$previous_r[m]) %>% unlist()
    }
    
    tmp[tmp$date == date_tmp, labels_target$new] <-   tmp[tmp$date == date_tmp - 1, labels_target$new] + dose2_tmp[[vac_type_tmp]][,labels_target$new]
    tmp[tmp$date == date_tmp, labels_target$previous] <- tmp[tmp$date == date_tmp, labels_target$previous] -  dose2_tmp[[vac_type_tmp]][,labels_target$new]
  }

  # allocate third doses
   dose3_tmp <- list()
   for(v in vac_type_list){
     vac_type_tmp <- v
     if(v != "sp"){
       labels_target <- data.frame(
         new = combo %>% 
           dplyr::select(-combo_index) %>% 
           filter(applicability == "Y", dose3 == vac_type_tmp) %>% 
           unique() %>% 
           mutate(labels = paste0("3_",dose1,"_",dose2,"_",vac_type_tmp)) %>% 
           pull(labels)
       ) %>% 
         distinct() %>% 
         mutate(previous = substr(new,1,7),
                previous = gsub("3_","2_",previous),
                previous_r = paste0(previous, "_r"))
       
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
                                                                  filter(date %in% date_tmp, dose_index == "boost", vac_type == vac_type_tmp) %>% 
                                                                  pull(value)) * dose3_tmp[[vac_type_tmp]]) %>% 
           pull(labels_target$previous_r[m]) %>% unlist()
       }
       tmp[tmp$date == date_tmp, labels_target$new] <-   tmp[tmp$date == date_tmp - 1, labels_target$new] + dose3_tmp[[vac_type_tmp]][,labels_target$new]
       tmp[tmp$date == date_tmp, labels_target$previous] <- tmp[tmp$date == date_tmp, labels_target$previous] -  dose3_tmp[[vac_type_tmp]][,labels_target$new] 
     }
   }
}

tmp %>% 
  dplyr::select(date, "2_sv_sv", "2_az_az", "2_sp_sp","2_pf_pf","2_md_md") -> res_proportional

# write_rds(tmp, paste0(data_path,"vaccine_market_proportional_results.rds"))
read_rds(paste0(data_path,"vaccine_market_proportional_results.rds"))%>% 
  dplyr::select(date, "2_sv_sv", "2_az_az", "2_sp_sp","2_pf_pf","2_md_md") %>% 
  ggplot(., aes(x = date, y = `2_sv_sv`)) +
  geom_point()

# tmp %>% 
#   dplyr::select(date, starts_with(c("1_", "2_","3_"))) %>% 
#   pivot_longer(starts_with(c("1_", "2_","3_"))) %>% 
#   ggplot(., aes(x= date, y = value, group = name, color = name)) +
#   geom_bar(stat = "identity",
#            position = "stack")
# 
# 
# tmp %>% 
#   tail(1) %>% 
#   pivot_longer(cols = colnames(.)[2:40]) %>% pull(value) %>% sum(., na.rm = T)
