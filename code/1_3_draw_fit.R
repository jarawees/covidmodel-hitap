draw_fit <- function(input,
                     country = "Thailand",
                     draw_end = T,
                     voc_features = voc_features_test,
                     fit_vac_threshold = 0.1,
                     dt_tmp = 0.3){
  
  # debug
  # input <- c(2, 10, 0.1, 2)
  # country = "Thailand"
  # draw_end = T
  # voc_features = voc_features_test %>% 
  #   mutate(change_u = 1)
  # fit_vac_threshold = 0.3
  # dt_tmp = 0.3
  
  iso3c_tmp <- countrycode::countrycode(country, "country.name", "iso3c")
  if(country == "Kosovo") iso3c_tmp <- "XKX"
  params_tmp <- list()
  
  tmp <- owid_vac %>% 
    dplyr::filter(iso_code == iso3c_tmp,
                  people_fully_vaccinated_per_hundred > fit_vac_threshold*100) 
  
  if(nrow(tmp) > 0){
    params_tmp[["fit_end"]] <- tmp %>% 
      pull(date) %>% 
      min()
  }
  
  if(nrow(tmp) == 0){
    params_tmp[["fit_end"]] <- "2022-06-30"
  }
  
  rm(tmp)
  
  # params_tmp[["fit_start"]] <- epi %>% 
  #   mutate(cs = cumsum(new_death)) %>%  
  #   dplyr::filter(cs >= 1) %>% 
  #   pull(txn_date) %>% 
  #   min()
  
  params_tmp[["fit_start"]] <- "2021-02-15"
  
  params_tmp[["country"]] <- country
  
  params_tmp[["date_switch"]] <- voc_phases %>% 
    dplyr::filter(threshold == dt_tmp,
                  country_code == iso3c_tmp) %>% 
    pull(week_min)
  
  params_tmp[["voc_names"]] <- voc_phases %>% 
    dplyr::filter(threshold == dt_tmp,
                  country_code == iso3c_tmp) %>% 
    pull(voc_name)
  
  testthat::expect_equal(length(params_tmp$voc_names), length(params_tmp$date_switch))
  
  lapply(seq_along(params_tmp[["voc_names"]]), function(i){
    voc_features %>% dplyr::filter(voc_name ==   params_tmp$voc_names[i])
  }) %>% 
    bind_rows() -> tmp
  
  params_tmp[["change_u"]] <- tmp$change_u
  params_tmp[["change_y"]] <- tmp$change_y
  params_tmp[["change_ve"]] <- tmp$change_ve
  params_tmp[["change_severity"]] <- tmp$change_severity
  rm(tmp)
  
  params_tmp[["draw_end"]] <- if_else(draw_end == T,
                                      "2022-12-31",
                                      as.character(params_tmp[["fit_end"]]))
  
  suppressWarnings(
    fit_gen_country_basics(
      country_tmp = params_tmp[["country"]],
      country_code_tmp = country_list %>% dplyr::filter(country == params_tmp[["country"]]) %>% pull(country_code),
      date_start =  as.character(ymd(params_tmp[["fit_start"]]) - 30),
      date_end = params_tmp[["fit_end"]],
      period_wn = input[4]*365,
      R0_assumed = input[1],
      period_wv_m2l = 1*365, 
      processes_set = burden_processes_all,
      period_wv_h2m = 1*365, 
      prob_v_p_2l = 0.33,
      prob_v_p_2m = 0.33,
      prob_v_b_l2m = 0,
      deterministic = TRUE,
      seed = input[2]
    ) %>%
      update_u_y(
        para = ., 
        country_tmp = params_tmp[["country"]],
        country_code_tmp = country_list %>% dplyr::filter(country == params_tmp[["country"]]) %>% pull(country_code),
        efficacy_baseline = efficacy_all
      ) %>%
      emerge_voc_burden(
        para = ., 
        country_tmp = params_tmp[["country"]],
        country_code_tmp = country_list %>% dplyr::filter(country == params_tmp[["country"]]) %>% pull(country_code),
        efficacy_baseline = efficacy_all
      )  -> tmp
  )
  
  tmp %>%
    cm_simulate() %>%
    .[["dynamics"]] %>%
    filter(grepl("death", compartment)) %>%
    group_by(t, compartment) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(date = ymd(as.character(ymd(params_tmp[["fit_start"]]) - 30)) + t) %>%
    pivot_wider(names_from = compartment,
                values_from = value) -> sim_deaths
  
  sim_deaths %>% 
    left_join(data.frame(date = ymd(params_tmp[["date_switch"]]),
                         phase = params_tmp$voc_names),
              by = "date") %>% 
    mutate(date_min = min(date),
           phase = if_else(date == date_min, "wildtype", phase),
           phase = zoo::na.locf(phase)) %>% 
    dplyr::select(-date_min) %>% 
    pivot_longer(cols = starts_with("death_")) %>% 
    mutate(endpoint = name,
           endpoint = gsub("death_voc_", "", endpoint),
           endpoint = gsub("death_", "", endpoint),
           endpoint = gsub("_o", "", endpoint),
           endpoint = if_else(endpoint == "o", "wildtype",endpoint)) %>% 
    dplyr::filter(phase == endpoint) %>% 
    mutate(scaled = value*as.numeric(input[3])) %>% 
    arrange(date)-> predicted
  
  predicted %>%
    right_join(epi %>% 
                 dplyr::select(txn_date, new_death) %>% 
                 rename(deaths = new_death,
                        date = txn_date) %>% 
                 mutate(date = ymd(date)),
               by = "date") %>% 
    rename(unscaled = value) %>% 
    arrange(date) %>% 
    dplyr::select(date, scaled, unscaled, deaths) %>% 
    rename(observed = deaths,
           predicted = scaled,
           predicted_unscaled = unscaled) %>% 
    ggplot(., aes(x = date)) +
    geom_vline(xintercept = ymd(params_tmp[["fit_end"]]),
               linetype = 2) +
    geom_point(aes(y = predicted), color = "green") +
    geom_point(aes(y = observed), color = "purple") -> p
  
  return(p)
}

