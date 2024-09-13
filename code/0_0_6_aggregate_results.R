aggregate_results <- function(dynamics_tmp = NULL,
                              by = NULL,
                              country_tmp = "Thailand",
                              country_code_tmp = "THA",
                              detection_threshold = 0.3){
  # debug
  # country_tmp = "Thailand"
  # country_code_tmp = "THA"
  # detection_threshold = 0.3
  # dynamics_tmp <- res$dynamics
  
  # extract voc information VOC
  # healthy_state_compartments <- c("S", "Sv_l", "Sv_m", "Sv_h", "E", "Ev_l", "Ev_m", 
  #                          "Ev_h", "R", "Rv_l", "Rv_m", "Rv_h")
  # all_state_compartments <- c(healthy_state_compartments, 
  #                             "Ip", "Ip_l", "Ip_m", "Ip_h", 
  #                             "Is", "Is_l", "Is_m", "Is_h",
  #                             "Ia", "Ia_l", "Ia_m", "Ia_h")
  
  voc_phases_imputation_index  %>% 
    dplyr::filter(country_code == country_code_tmp) %>% 
    pull(voc_phases_source) -> source_country_tmp
  country_list %>% 
    dplyr::filter(country == source_country_tmp) %>% 
    pull(country_code) -> source_country_code_tmp
  voc_phases %>% 
    dplyr::filter(country_code == source_country_code_tmp,
                  threshold == detection_threshold) -> voc_phases_tmp
  voc_phases_tmp %>% arrange(week_min) %>% pull(week_min) -> date_switch
  voc_phases_tmp %>% arrange(week_min) %>% pull(voc_name) -> voc_phases_names
  
  dynamics_tmp[, date := ymd("2020-01-01") + t] %>% 
    .[, voc_phases := cut(date,
                          lubridate::ymd(c("2020-01-01", as.character(date_switch), "2101-01-01")),
                          labels = c("wildtype", voc_phases_names))]
  
  dynamics_tmp <- dynamics_tmp[!grepl("_p", compartment)]
  dynamics_tmp %>% 
    .[, bin_death := grepl("death", compartment)] %>% 
    .[, bin_critical := grepl("critical", compartment)] %>% 
    .[, bin_severe := grepl("severe", compartment)] %>% 
    .[, bin_all_states := compartment %in% compartment_pop]
  
  unique(dynamics_tmp[,"compartment"]) %>% 
    .[, c("seg1", "seg2", "voc_phase_simulated", "seg4") := tstrsplit(compartment, "_")] %>% 
    .[, voc_phase_simulated := if_else(is.na(voc_phase_simulated) & seg1 %in% c("severe", "critical", "death"), "wildtype", voc_phase_simulated)] %>% 
    .[, .(compartment, voc_phase_simulated)] %>% 
    .[!is.na(voc_phase_simulated)] -> index_voc_phase
  
  index_voc_phase[dynamics_tmp, on = "compartment"] %>% 
    .[(voc_phases == voc_phase_simulated | is.na(voc_phase_simulated))] -> dynamics_tmp
  
  dynamics_tmp[, compartment := if_else(bin_death == TRUE, "death", compartment)]
  dynamics_tmp[, compartment := if_else(bin_severe == TRUE, "severe", compartment)]
  dynamics_tmp[, compartment := if_else(bin_critical == TRUE, "critical", compartment)]
  dynamics_tmp[, cohort_all := sum(value*bin_all_states), by = .(date, group)]
  
  dynamics_tmp <- dynamics_tmp[compartment %in% c("cases", "death", "severe", "critical")]
  dynamics_tmp[, year := lubridate::year(date)]
  
  if(by == "year"){
    dynamics_tmp <- dynamics_tmp[, keyby = .(year, compartment, group), .(incidence = sum(value), cohort_all = mean(cohort_all))]
  }
  
  if(by == "day"){
    dynamics_tmp <- dynamics_tmp[, keyby = .(date, compartment, group), .(incidence = sum(value), cohort_all = mean(cohort_all))]
  }
  

  dynamics_tmp[, prop := incidence/cohort_all]
  data.table(group=unique(dynamics_tmp[,group]),
             group_index = 0:15) -> group_index
  output <- group_index[dynamics_tmp, on = "group"]

  return(output)
}


