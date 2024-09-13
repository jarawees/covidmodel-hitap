gen_country_basics <- function(country_tmp = "Thailand",
                               country_code_tmp = "THA",
                               date_start = "2020-01-01",
                               date_end = "2023-12-31",
                               # duration, waning of natural immunity
                               # duration, waning from medium to low levels vaccine induced 
                               # this needs to be pre-calculated, generated from 
                               # `gen_burden_processes` with special sets of 
                               # vaccine efficacies
                               processes_set = burden_processes_all,
                               # duration, waning from medium to low levels vaccine induced 
                               prob_v_p_2l = 0.33,
                               prob_v_p_2m = 0.33,
                               prob_v_b_l2m = 0,
                               fitted_table_tmp = NULL, # sensitivity analysis 1
                               wn_lt = 3*365, # sensitivity analysis 2
                               period_wv_h2m = 1*365, # sensitivity analysis 3
                               period_wv_m2l = 1*365, # sensitivity analysis 3
                               # reduction in susceptibility among previously 
                               # infected individuals
                               deterministic = TRUE){

  # country_tmp = "Thailand"
  # country_code_tmp = "THA"
  # date_start = "2020-01-01"
  # date_end = "2025-12-31"
  # processes_set = burden_processes_all
  # prob_v_p_2l = 0.33
  # prob_v_p_2m = 0.33
  # prob_v_b_l2m = 0
  # fitted_table_tmp = out_all[1,]
  # period_wn  = 3*365 # duration, waning of natural immunity
  # period_wv_m2l = 1*365 # duration, waning from medium to low levels vaccine induced
  # period_wv_h2m = 1*365 # duration, waning from medium to low levels vaccine induced
  # deterministic = TRUE

  if(!exists("contact_schedule")){stop("contact_schedule has not been loaded yet.")}
  contact_tmp <- 
    contact_schedule %>% 
    filter(country_code == country_code_tmp) %>% 
    filter(date >= date_start,
           date <= date_end)
  
  # 
  processes_tmp <- 
    processes_set[[country_code_tmp]]
  
  #### birth rate modifications ####
  if(!exists("cbr")){stop("cbr has not been loaded yet.")}
  rate_birth <- 
    cbr %>% 
    dplyr::filter(country_code == country_code_tmp,
                  year == 2020) %>% 
    pull(cbr_daily)
  testthat::expect_length(rate_birth, 1)

  seq(as.numeric(substr(date_start, 1, 4)),
      as.numeric(substr(date_end, 1, 4))) %>% 
    paste0(., "-01-01") %>% 
    c(., date_start) %>% 
    lubridate::ymd() %>% 
    sort %>% 
    enframe(value = "date") %>% 
    dplyr::filter(date >= lubridate::ymd(date_start)) %>% 
    mutate(year = as.numeric(substr(as.character(date), 1, 4)),
           country_code = country_code_tmp) -> modification_dates
    
  modification_dates %>% 
    left_join(cbr,
              by = c("country_code", "year")) %>% 
    left_join(data.frame(date = seq(ymd(date_start),
                                    ymd(date_end),
                                    "day")) %>% 
                mutate(t = 1:n()) %>% 
                filter(date %in% .$date),
              by = "date") %>% 
    dplyr::select(-name) %>% 
    distinct() -> cbr_changes
  
  #### mortality rate modifications ####
  if(!exists("mu_weighted_16")){stop("mu_weighted_16 has not been loaded yet.")}
  rate_death <- 
    mu_weighted_16 %>% 
    dplyr::filter(country_code == country_code_tmp,
                  year == 2020) %>% 
    arrange(age_from) %>%
    pull(mu_daily)
  testthat::expect_length(rate_death, 16)
  
  modification_dates %>% 
    left_join(mu_weighted_16,
              by = c("country_code", "year"),
              relationship = "many-to-many") %>% 
    left_join(data.frame(date = seq(ymd(date_start),
                                    ymd(date_end),
                                    "day")) %>% 
                mutate(t = 1:n()) %>% 
                filter(date %in% .$date),
              by = "date") %>% 
    dplyr::select(-name) %>% 
    distinct() -> mu_changes

  #### generate the actual parameter set ####
  if(country_tmp != "United Kingdom"){
    para = cm_parameters_SEI3R(dem_locations = country_tmp, 
                               date_start = date_start, 
                               date_end = date_end,
                               A = rep(1/(365*5),16),
                               B = c(rate_birth, rep(0, 15)),
                               D = rate_death,
                               dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                               dEa = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                               dIp = cm_delay_gamma(1.5, 4.0, t_max = 15, t_step = 0.25)$p,
                               dIs = cm_delay_gamma(3.5, 4.0, t_max = 15, t_step = 0.25)$p,
                               dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p,
                               deterministic = deterministic)
  }
  
  if(country_tmp == "United Kingdom"){
    para = cm_parameters_SEI3R(dem_locations = "United Kingdom", 
                               mat_locations = "GBR",
                               date_start = date_start, 
                               date_end = date_end,
                               A = rep(1/(365*5),16),
                               B = c(rate_birth, rep(0, 15)),
                               D = rate_death,
                               dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                               dEa = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                               dIp = cm_delay_gamma(1.5, 4.0, t_max = 15, t_step = 0.25)$p,
                               dIs = cm_delay_gamma(3.5, 4.0, t_max = 15, t_step = 0.25)$p,
                               dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p,
                               deterministic = deterministic)
  }
  
  n_age_groups <- length(para$pop[[1]]$size)
  seed <- fitted_table_tmp %>% dplyr::filter(country_code == country_code_tmp) %>% pull(seed_20200101)
  seeds <- seed:(seed+14)
  
  for(i in 1:length(para$pop)){
    
    para$pop[[i]]$y <- cf
    para$pop[[i]]$u <- sus
    para$pop[[i]]$v_p_2l <- rep(prob_v_p_2l, 16)
    para$pop[[i]]$v_p_2m <- rep(prob_v_p_2m, 16)
    para$pop[[i]]$v_b_l2m <- rep(prob_v_b_l2m, 16)
    
    # scale u (susceptibility) to achieve desired R0
    current_R0 = cm_calc_R0(para, i); # calculate R0 in population i of params
    R0_assumed = fitted_table_tmp %>% dplyr::filter(country_code == country_code_tmp) %>% pull(R0_assumed_2)
    para$pop[[i]]$u = para$pop[[i]]$u * R0_assumed / current_R0
    
    # The purpose of this chunk of code is to update uv_l, uv_m, uv_h, ur, uvr_l
    # uvr_m, uvr_h to be consistent with u and to update yv_l, yv_m, and yv_h 
    # to be consistent with yv. we will not implement efficacy at this step 
    # just yet.
    
    para$pop[[i]]$uv_l  <- para$pop[[i]]$u
    para$pop[[i]]$uv_m  <- para$pop[[i]]$u
    para$pop[[i]]$uv_h  <- para$pop[[i]]$u
    para$pop[[i]]$uvr_l  <- para$pop[[i]]$u
    para$pop[[i]]$uvr_m  <- para$pop[[i]]$u
    para$pop[[i]]$uvr_h  <- para$pop[[i]]$u
    para$pop[[i]]$ur  <- para$pop[[i]]$u
    
    para$pop[[i]]$yv_l <- para$pop[[i]]$y
    para$pop[[i]]$yv_m <- para$pop[[i]]$y
    para$pop[[i]]$yv_h <- para$pop[[i]]$y
    
    ## Set seeds to control start of outbreak
    # infections start in individuals aged 20-50
    para$pop[[i]]$dist_seed_ages = 
      cm_age_coefficients(20, 
                          80, 
                          5 * (0:length(para$pop[[i]]$size))) 
    
    # 1 new infections each day for 14 days to see the outbreak
    para$pop[[i]]$seed_times <- seeds
  }
  
  para$processes = processes_tmp
  
  # NPI modified contacts
  para$schedule[["mobility"]] = list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    # values and times need to be the same length
    # values need to be a list
    # times need to be an array
    # this is true for all schedule objects
    values = split(contact_tmp[,c("home", "work", "school", "others")],
                   seq(nrow(contact_tmp))) %>%
      map(unlist) %>%
      map(as.vector) %>% 
      unname, #remove list structure, convert to vector, and remove item names
    times = 1:nrow(contact_tmp))
  
  # implement changing birth rates
  para$schedule[["birth"]] = list(
    parameter = "B",
    pops = numeric(),
    mode = "assign",
    values = cbr_changes %>% 
      split(., 1:nrow(.)) %>% 
      map(pull, cbr_daily_approximate) %>% 
      map(~c(., rep(0, 15))) %>% 
      map(unname) %>% 
      unname,
    times = cbr_changes$t)
  
  # implement changing death rateså
  para$schedule[["death"]] = list(
    parameter = "D",
    pops = numeric(),
    mode = "assign",
    values = mu_changes %>% 
      group_by(t) %>% 
      group_split() %>% 
      map(arrange, age_from) %>% 
      map(pull, mu_daily_approximate),
    times = sort(unique(mu_changes$t)))
  
  # waning vaccine-induced immunity
  # period_wn  <- (fitted_table_tmp %>% dplyr::filter(country_code == country_code_tmp) %>% pull(wn))*365
  period_wn <- 3*365
  para$pop[[1]]$wn     <- rep(1/period_wn, n_age_groups)
  para$pop[[1]]$wv_m2l <-  rep(1/period_wv_m2l, n_age_groups)
  para$pop[[1]]$wv_h2m <-  rep(1/period_wv_h2m, n_age_groups)
  
  # para$schedule[["wn_change"]] = list(
  #   parameter = "wn",
  #   pops = numeric(),
  #   mode = "assign",
  #   values = list(rep(1/(wn_lt), n_age_groups)),
  #   times = 911)
 # lubridate::ymd("2023-12-31") - lubridate::ymd("2020-01-01") # 1460 days
  
  return(para)
}
