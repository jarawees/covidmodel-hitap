emerge_voc_burden <- function(para = NULL,
                              country_tmp = "Thailand",
                              country_code_tmp = "THA",
                              detection_threshold = 0.3,
                              future_severe = F,
                              future_severe_lvl = "mean", # sensitivity analysis 4
                              efficacy_baseline = NULL){
  # debug 
  # country_tmp = country_list$country[i]
  # country_code_tmp = country_list$country_code[i]
  # detection_threshold = 0.3
  # efficacy_baseline = efficacy_all
  # future_severe = T
  
  if(!exists("voc_features_test")) stop("voc_features_test is not loaded!")
  if(!exists("HSR_cleaned")) stop("HSR_cleaned is not loaded!")
  
  # put toghether voc stages
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
  
  voc_features_tmp <- list()
  for(j in 1:length(voc_phases_names)){
    voc_features_tmp[[j]] <- voc_features_test %>% dplyr::filter(voc_name == voc_phases_names[j])
  }
  voc_features_tmp %<>% bind_rows()
  
  if(future_severe == T){
    if(future_severe_lvl == "mean"){
      voc_features_tmp  %<>% 
        bind_rows(severe_strain_definition)
    }
    
    if(future_severe_lvl == "low"){
      voc_features_tmp  %<>% 
        bind_rows(severe_strain_definition_low)
    }
    
    if(future_severe_lvl == "high"){
      voc_features_tmp  %<>% 
        bind_rows(severe_strain_definition_high)
    }
    
    date_switch <- c(date_switch, "2024-01-01")
    voc_phases_names <- c(voc_phases_names, severe_strain_definition$voc_name)
  }
  
  #
  change_severity <- voc_features_tmp %>% pull(change_severity)
  
  n_voc = length(change_severity)
  expect_equal(length(change_severity), length(voc_phases_names))
  expect_equal(length(change_severity), length(date_switch))
  compartments_E <- c("newE", "newEv_l", "newEv_m", "newEv_h")
  P.death <- HSR_cleaned[[country_code_tmp]]$P.death
  P.severe <- HSR_cleaned[[country_code_tmp]]$P.severe
  P.critical <- HSR_cleaned[[country_code_tmp]]$P.critical
  P.hosp <- HSR_cleaned[[country_code_tmp]]$P.hosp
  
  
  # generate death processes
  generate_death_processes <- function(source_compartment = NULL,
                                       voc_index = NULL){
    multiplier1 <- change_severity[voc_index]
    if(source_compartment == "newE") multiplier2 <- 1
    if(source_compartment == "newEv_l") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[1]
    if(source_compartment == "newEv_m") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[2]
    if(source_compartment == "newEv_h") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[3]
    tmp_var <- paste0("death_voc_", voc_phases_names[voc_index])
    
    tmp_process <- cm_multinom_process(source_compartment,
                                       data.frame(P.death*multiplier1*multiplier2) |> setNames(tmp_var),
                                       delays = data.frame(delay_2death) |> setNames(tmp_var),
                                       report = "o")
    
    return(tmp_process)
  }
  
  output_table <- data.table::CJ(source_compartment = compartments_E,
                                 voc_index = 1:n_voc)
  
  death_processes <- lapply(1:nrow(output_table), 
                            function(x) {generate_death_processes(source_compartment = output_table$source_compartment[x],
                                                                  voc_index = output_table$voc_index[x])}
  )
  
  # generate severe and critical, intermediate processes
  generate_intermediate_processes <- function(source_compartment = NULL,
                                              voc_index = NULL){
    multiplier1 <- change_severity[voc_index]
    
    if(source_compartment == "newE") {
      multiplier2_severe <- multiplier2_critical <- 1
    }
    
    if(source_compartment == "newEv_l") {
      multiplier2_severe <- 1 - efficacy_baseline$v_severe_condition[1]
      multiplier2_critical <- 1 - efficacy_baseline$v_critical_condition[1]
    }
    
    if(source_compartment == "newEv_m") {
      multiplier2_severe <- 1 - efficacy_baseline$v_severe_condition[2]
      multiplier2_critical <- 1 - efficacy_baseline$v_critical_condition[2]
    }
    
    if(source_compartment == "newEv_h") {
      multiplier2_severe <- 1 - efficacy_baseline$v_severe_condition[3]
      multiplier2_critical <- 1 - efficacy_baseline$v_critical_condition[3]
    }
    
    tmp_var_to_severe <- paste0("to_severe_voc_", voc_phases_names[voc_index])
    tmp_var_to_critical <- paste0("to_critical_voc_", voc_phases_names[voc_index])
    tmp_var_severe <- paste0("severe_voc_", voc_phases_names[voc_index])
    tmp_var_critical <- paste0("critical_voc_", voc_phases_names[voc_index])
    
    tmp_process <-
      cm_multinom_process(source_compartment,
                          data.frame(var_to_severe = P.severe*multiplier1*multiplier2_severe,
                                     var_to_critical = P.critical*multiplier1*multiplier2_critical) |> 
                            setNames(c(tmp_var_to_severe, tmp_var_to_critical)),
                          delays = data.frame(var_to_severe = delay_2severe,
                                              var_to_critical = delay_2severe) |> 
                            setNames(c(tmp_var_to_severe, tmp_var_to_critical)))
    
    
    return(tmp_process)
  }
  
  intermediate_processes <- lapply(1:nrow(output_table), 
                                   function(x) {generate_intermediate_processes(source_compartment = output_table$source_compartment[x],
                                                                                voc_index = output_table$voc_index[x])}
  ) 
  
  intermediate_processes_agg <- list()
  
  for(j in 1:n_voc){
    tmp_var_to_severe <- paste0("to_severe_voc_", voc_phases_names[j])
    tmp_var_to_critical <- paste0("to_critical_voc_", voc_phases_names[j])
    tmp_var_severe <- paste0("severe_voc_", voc_phases_names[j])
    tmp_var_critical <- paste0("critical_voc_", voc_phases_names[j])
    
    list(  
      cm_multinom_process(tmp_var_to_severe, 
                          data.frame(rep(1,16)) |> setNames(tmp_var_severe),
                          delays = data.frame(delay_2hosp)|> 
                            setNames(tmp_var_severe),   report = "ip"),
      
      cm_multinom_process(tmp_var_to_critical,
                          data.frame(rep(1,16))|> setNames(tmp_var_critical),
                          delays = data.frame(delay_2hosp_critical)|> 
                            setNames(tmp_var_critical),   report = "ip")
    ) -> intermediate_processes_agg[[j]]
  }
  
  intermediate_processes_agg |> purrr::flatten() -> intermediate_processes_agg
  to_attach <- c(death_processes, 
                 intermediate_processes,
                 intermediate_processes_agg)
  
  para$processes <- c(para$processes, to_attach)
  return(para)
}

