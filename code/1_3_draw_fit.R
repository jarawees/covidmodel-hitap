draw_fit <- function(input,
                     country = "Thailand",
                     draw_end = T,
                     voc_features = voc_features_test %>% mutate(change_u = 1),
                     fit_vac_threshold = 0.1,
                     detection_threshold = 0.3){
  
  iso3c_tmp <- countrycode::countrycode(country, "country.name", "iso3c")
  if(country == "Kosovo") iso3c_tmp <- "XKX"
  params_tmp <- list()
  
  tmp <- owid_vac %>% 
    dplyr::filter(country_code == iso3c_tmp,
                  people_fully_vaccinated_per_hundred > 100*fit_vac_threshold)
  
  if(nrow(tmp) > 0){
    params_tmp[["fit_end"]] <- tmp %>% 
      pull(date) %>% 
      min()
  }
  
  if(nrow(tmp) == 0){
    params_tmp[["fit_end"]] <- "2022-06-30"
  }
  rm(tmp)
  
  params_tmp[["fit_start"]] <- lubridate::ymd("2020-03-01")
  
  params_tmp[["country"]] <- country
  
  params_tmp[["date_switch"]] <- voc_phases %>% 
    dplyr::filter(threshold == detection_threshold,
                  country_code == iso3c_tmp) %>% 
    pull(week_min)
  
  if(length(params_tmp[["date_switch"]]) == 0){
    params_tmp[["date_switch"]] <- impute_phase(country,
                                                detection_threshold,
                                                1) %>% 
      arrange(week) %>% 
      pull(week)
  }
  
  params_tmp[["voc_names"]] <- voc_phases %>% 
    dplyr::filter(threshold == detection_threshold,
                  country_code == iso3c_tmp) %>% 
    pull(voc_name)
  
  if(length(params_tmp[["voc_names"]]) == 0){
    params_tmp[["voc_names"]] <- impute_phase(country, 
                                              detection_threshold, 
                                              1) %>% 
      arrange(week) %>% 
      pull(voc_name)
  }
  
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
      date_start =  as.character(params_tmp[["fit_start"]] - 30),
      date_end = params_tmp[["fit_end"]],
      period_wn = 3*365,
      R0_assumed = input[1],
      # duration, waning of natural immunity
      # duration, waning from medium to low levels vaccine induced 
      period_wv_m2l = 1*365, 
      # this needs to be pre-calculated, generated from 
      # `gen_burden_processes` with special sets of 
      # vaccine efficacies
      processes_set = burden_processes_all,
      # duration, waning from medium to low levels vaccine induced 
      period_wv_h2m = 1*365, 
      prob_v_p_2l = 0.33,
      prob_v_p_2m = 0.33,
      prob_v_b_l2m = 0,
      # reduction in susceptibility among previously 
      # infected individuals
      deterministic = TRUE,
      seed = input[2]
    ) %>%
      update_u_y(
        para = ., 
        country_tmp = params_tmp[["country"]],
        country_code_tmp = country_list %>% dplyr::filter(country == params_tmp[["country"]]) %>% pull(country_code),
        detection_threshold = 0.3,
        efficacy_baseline = efficacy_all,
        voc_features_inuse = voc_features
      ) %>%
      emerge_voc_burden(
        para = ., 
        country_tmp = params_tmp[["country"]],
        country_code_tmp = country_list %>% dplyr::filter(country == params_tmp[["country"]]) %>% pull(country_code),
        detection_threshold = 0.3,
        split_E = T,
        voc_features_inuse = voc_features,
        efficacy_baseline = efficacy_all
      )  -> tmp
  )
  
  tmp %>%
    cm_simulate() %>%
    .[["dynamics"]] %>%
    filter(grepl("death", compartment)) %>%
    group_by(t, compartment) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(date = ymd(as.character(params_tmp[["fit_start"]] - 30)) + t) %>%
    pivot_wider(names_from = compartment,
                values_from = value) -> sim_deaths
  
  sim_deaths %>% 
    arrange(date) %>% 
    left_join(data.frame(date = ymd(params_tmp[["date_switch"]]),
                         phase = c(1:length(params_tmp[["date_switch"]]))),
              by = "date") %>% 
    mutate(date_min = min(date),
           phase = if_else(date == date_min, 0, phase),
           phase = na_locf(phase)) %>% 
    dplyr::select(-date_min) %>% 
    left_join(data.frame(phase_name = c("wildtype", params_tmp[["voc_names"]]),
                         phase = c(0:length(params_tmp[["date_switch"]]))),
              by = "phase") %>% 
    pivot_longer(cols = starts_with("death_")) %>% 
    separate(name, into = c("seg1", "seg2", "endpoint", "seg4"), sep = "_", fill = "right") %>% 
    dplyr::select(-seg1, -seg2, -seg4) %>% 
    mutate(endpoint = if_else(is.na(endpoint), "wildtype", endpoint)) %>% 
    dplyr::filter(phase_name == endpoint) %>% 
    mutate(scaled = value*as.numeric(input[3])) -> predicted
  
  predicted %>%
    right_join(epi %>% 
                 dplyr::select(txn_date, new_death) %>% 
                 rename(deaths = new_death,
                        date = txn_date) %>% 
                 mutate(date = ymd(date)),
               by = "date") %>% 
    rename(unscaled = value) %>% 
    arrange(date) %>% 
    # filter(date <= params_tmp[["fit_end"]]) %>% 
    dplyr::select(date, scaled, unscaled, deaths) %>% 
    rename(observed = deaths,
           predicted = scaled,
           predicted_unscaled = unscaled) %>% 
    ggplot(., aes(x = date)) +
    geom_point(aes(y = predicted), color = "green") +
    # geom_point(aes(y = predicted_unscaled), color = "red") +
    geom_point(aes(y = observed), color = "purple") -> p
  
  return(p)
}
