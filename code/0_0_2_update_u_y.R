update_u_y <- function(para = NULL,
                       country_tmp = "Thailand",
                       country_code_tmp = "THA",
                       detection_threshold = 0.3,
                       voc_features_inuse = voc_features_test,
                       future_severe = F,
                       future_severe_lvl = "mean", # sensitivity analysis 4
                       efficacy_baseline = NULL# vaccine efficacy
){
  # debug
  # country_tmp = country_list$country[country_index]
  # country_code_tmp = country_list$country_code[country_index]
  # detection_threshold = 0.3
  # efficacy_baseline = efficacy_all
  # future_severe = F
  # voc_features_inuse = voc_features_test

  if(!exists("country_list")) stop("country_list is not loaded!")
  if(!exists("voc_phases")) stop("voc_phases is not loaded!")
  if(!exists("voc_phases_imputation_index")) stop("voc_phases_imputation_index is not loaded!")
  
  # extract information for observed VoCs
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
    voc_features_tmp[[j]] <- voc_features_inuse %>% dplyr::filter(voc_name == voc_phases_names[j])
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
  
  change_u <- voc_features_tmp %>% pull(change_u)
  change_y <- voc_features_tmp %>% pull(change_y)  
  change_ve <- voc_features_tmp %>% pull(change_ve)
  
  date_marker <- c(as.character(date_switch), as.character(lubridate::ymd(para$date0) + para$time1))
  expect_equal(length(date_switch), length(change_u))
  if(para$date0 < date_marker[1]) date_marker <- c(para$date0, date_marker)
  if(para$date0 > date_marker[1]) date_marker[1] <- para$date0
  t_range <- as.numeric(ymd(date_marker) - ymd(para$date0))
  
  # this is the table that will help us keep track of the names of things to be
  # changed in the schedule
  targets <- data.frame(
    scaler_label = c(
      "u_scaler",
      "uv_l_scaler",
      "uv_m_scaler",
      "uv_h_scaler",
      "ur_scaler",
      "uvr_l_scaler",
      "uvr_m_scaler",
      "uvr_h_scaler",
      "yv_l_scaler",
      "yv_m_scaler",
      "yv_h_scaler"
    ),
    
    variable_label = c("u",
                       "uv_l",
                       "uv_m",
                       "uv_h",
                       "ur",
                       "uvr_l",
                       "uvr_m",
                       "uvr_h",
                       "yv_l",
                       "yv_m",
                       "yv_h")
  )
  
  n_age_groups <- para$pop[[1]]$n_groups
  
  # assign initial estimates
  para$pop[[1]]$uv_l <- (1 - efficacy_baseline %>% dplyr::filter(protection_level_label == "l") %>%  pull(v_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uv_m <- (1 - efficacy_baseline %>% dplyr::filter(protection_level_label == "m") %>%  pull(v_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uv_h <- (1 - efficacy_baseline %>% dplyr::filter(protection_level_label == "h") %>%  pull(v_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$ur    <- (1 - efficacy_baseline$r_i_o[2]) * para$pop[[1]]$u
  para$pop[[1]]$uvr_l <- (1 - efficacy_baseline %>%  dplyr::filter(protection_level_label == "l") %>%  pull(vr_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uvr_m <- (1 - efficacy_baseline %>%  dplyr::filter(protection_level_label == "m") %>%  pull(vr_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uvr_h <- (1 - efficacy_baseline %>%  dplyr::filter(protection_level_label == "h") %>%  pull(vr_i_o)) * para$pop[[1]]$u
  
  para$pop[[1]]$yv_l <- para$pop[[1]]$y*(1 - efficacy_baseline$v_d_condition[1])
  para$pop[[1]]$yv_m <- para$pop[[1]]$y*(1 - efficacy_baseline$v_d_condition[2])
  para$pop[[1]]$yv_h <- para$pop[[1]]$y*(1 - efficacy_baseline$v_d_condition[3])
  
  # create modifier table
  data.frame(date = date_marker,
             phase = as.numeric(NA)) -> modifier
  
  for(j in 1:length(date_marker)){
    if(nrow(modifier) == (length(date_switch) + 1)){
      modifier[modifier$date == date_marker[j],"phase"] <- j
      modifier[modifier$date == date_marker[j],"diff_u"]  <- c(change_u)[j]
      modifier[modifier$date == date_marker[j],"diff_y"]  <- c(change_y)[j]
      modifier[modifier$date == date_marker[j],"diff_ve"] <- c(change_ve)[j]
    }
    
    if(nrow(modifier) == (length(date_switch) + 2)){
      modifier[modifier$date == date_marker[j],"phase"] <- j
      modifier[modifier$date == date_marker[j],"diff_u"]  <- c(1, change_u)[j]
      modifier[modifier$date == date_marker[j],"diff_y"]  <- c(1, change_y)[j]
      modifier[modifier$date == date_marker[j],"diff_ve"] <- c(1, change_ve)[j]
    }
  }
  
  modifier %<>% 
    mutate(diff_u = na_locf(diff_u),
           diff_y = na_locf(diff_y),
           diff_ve = na_locf(diff_ve))
  
  # VEs against infection and disease are implemented over "compartments"
  # VEs against severe, critical and mortality cases are implemented over "processes" 
  # Everything above infection in terms of outcome will need to use conditional
  # probability, because the infection step has already occurred to reach this 
  # endpoint
  
  # in the context of this model, disease preventing = clinical preventing
  modifier |> 
    mutate(u_scaler     = diff_u,
           uv_l_scaler  = diff_u*(1 - efficacy_baseline$v_i_o[1]*diff_ve)/(1 - efficacy_baseline$v_i_o[1]),
           uv_m_scaler  = diff_u*(1 - efficacy_baseline$v_i_o[2]*diff_ve)/(1 - efficacy_baseline$v_i_o[2]),
           uv_h_scaler  = diff_u*(1 - efficacy_baseline$v_i_o[3]*diff_ve)/(1 - efficacy_baseline$v_i_o[3]),
           ur_scaler    = diff_u,
           uvr_l_scaler = diff_u*(1 - efficacy_baseline$vr_i_o[1]*diff_ve)/(1 - efficacy_baseline$vr_i_o[1]),
           uvr_m_scaler = diff_u*(1 - efficacy_baseline$vr_i_o[2]*diff_ve)/(1 - efficacy_baseline$vr_i_o[2]),
           uvr_h_scaler = diff_u*(1 - efficacy_baseline$vr_i_o[3]*diff_ve)/(1 - efficacy_baseline$vr_i_o[3]),
           yv_l_scaler  = diff_y,
           yv_m_scaler  = diff_y,
           yv_h_scaler  = diff_y) -> modifier
  
  modifier |> 
    dplyr::select(ends_with("scaler")) |> 
    pivot_longer(targets$scaler_label) |> 
    group_by(name) |> summarise(value = min(value)) |> 
    pull(value) |> (function(y) y > 0)() |> all() -> test_range
  
  testthat::expect(test_range,
                   failure_message = "zero scaler values generated. 
                   please double check all change_xx variables.")
  
  for(j in seq_len(nrow(targets))){
    tmp <-   modifier |> 
      pull(targets$scaler_label[j]) |> 
      split(seq(nrow(modifier))) |> 
      map(rep, n_age_groups) |> 
      map(unname)
    
    para$schedule[[targets$scaler_label[j]]] <-  list(
      parameter = targets$variable_label[j],
      pops = numeric(),
      mode = "multiply",
      values = tmp,
      times = t_range
    )
    
    rm(tmp)
  }
  
  return(para)
}
