# foundamental wrapper of the functions
gen_country_basics <- function(country = "Thailand",
                               # fitting
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2030-12-31",
                               processes = NULL,
                               contact = contact_schedule,
                               # duration, waning of natural immunity
                               period_wn  = 3*365, 
                               # duration, waning from medium to low levels vaccine induced 
                               period_wv_m2l = 1*365, 
                               # duration, waning from medium to low levels vaccine induced 
                               period_wv_h2m = 1*365, 
                               prob_v_p_2l = 0.5,
                               prob_v_p_2m = 0.3,
                               prob_v_b_l2m = 0.5,
                               # birth rate just from macrotrends.net
                               rate_birth = c((9.532/1000)/365, rep(0,15)), # rate per day
                               rate_death = mu_inuse$mu_mean_day, # rate per day
                               # reduction in susceptibility among previously 
                               # infected individuals
                               deterministic = TRUE){
  
  require(countrycode)
  
  iso3c_tmp = countrycode(country, "country.name", "iso3c")
  
  c_tmp = contact %>% 
    filter(country_region == country) %>% 
    filter(date >= date_start,
           date <= date_end)
  
  para = cm_parameters_SEI3R(dem_locations = as.character(country), 
                             date_start = date_start, 
                             date_end = date_end,
                             A = rep(1/(365*5),16),
                             B = rate_birth,
                             D = rate_death,
                             dE  = cm_delay_gamma(2.5, 2.5,
                                                  t_max = 15, t_step = 0.25)$p,
                             dEa = cm_delay_gamma(2.5, 2.5,
                                                  t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(1.5, 4.0,
                                                  t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(3.5, 4.0,
                                                  t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0,
                                                  t_max = 15, t_step = 0.25)$p,
                             deterministic = deterministic)
  

  n_age_groups <- length(para$pop[[1]]$size)
  
  for(i in 1:length(para$pop)){
    
    para$pop[[i]]$y <- cf
    para$pop[[i]]$u <- sus
    para$pop[[i]]$v_p_2l <- rep(prob_v_p_2l, 16)
    para$pop[[i]]$v_p_2m <- rep(prob_v_p_2m, 16)
    para$pop[[i]]$v_b_l2m <- rep(prob_v_b_l2m, 16)
    
    # scale u (susceptibility) to achieve desired R0
    current_R0 = cm_calc_R0(para, i); # calculate R0 in population i of params
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
    para$pop[[i]]$yvr_l <- para$pop[[i]]$y
    para$pop[[i]]$yvr_m <- para$pop[[i]]$y
    para$pop[[i]]$yvr_h <- para$pop[[i]]$y

    ## Set seeds to control start of outbreak
    # infections start in individuals aged 20-50
    para$pop[[i]]$dist_seed_ages = 
      cm_age_coefficients(20, 
                          80, 
                          5 * (0:length(para$pop[[i]]$size))) 
    
    # 1 new infections each day for 14 days to see the outbreak
    para$pop[[i]]$seed_times <- c(1:14)
  }
  
  para$processes = processes
  
  # this is a schedule oject
  para$schedule[["mobility"]] = list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    # values and times need to be the same length
    # values need to be a list
    # times need to be an array
    # this is true for all schedule objects
    values = split(c_tmp[,3:6],
                   seq(nrow(c_tmp))) %>%
      map(unlist) %>%
      map(as.vector) %>% 
      unname, #remove list structure, convert to vector, and remove item names
    times = 1:nrow(c_tmp))
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wn     <- rep(1/period_wn, n_age_groups)
  para$pop[[1]]$wv_m2l <-  rep(1/period_wv_m2l, n_age_groups)
  para$pop[[1]]$wv_h2m <-  rep(1/period_wv_h2m, n_age_groups)
  
  return(para)
}

# CJ(date = seq(ymd("2019-12-01"),
#               ymd("2030-12-31"),
#               by = "day")) %>%
#   mutate(weights_v_i_l = rnorm(nrow(.), 0.8, 0.05),
#          weights_v_i_m = rnorm(nrow(.), 0.8, 0.05),
#          weights_v_i_h = rnorm(nrow(.), 0.8, 0.05),
#          weights_v_d_l = rnorm(nrow(.), 0.8, 0.05),
#          weights_v_d_m = rnorm(nrow(.), 0.8, 0.05),
#          weights_v_d_h = rnorm(nrow(.), 0.8, 0.05)) -> efficacy_weights_test
# 
# write_rds(efficacy_weights_test,
#           paste0(data_path, 
#                  "intermediate/efficacy_weights_test.rds"))

efficacy_weights_test <- read_rds(paste0(data_path, 
                                         "intermediate/efficacy_weights_test.rds"))

# efficacy_weights_test |> 
#   mutate_at(vars(starts_with("weights")), function(x) x <- 1) -> efficacy_weights_one

# efficacy_weights_test |> 
#   ggplot(aes(x = date, y = weights_v_d_h)) +
#   geom_point() +
#   geom_line()

# this function will update u and y
# there's two sources for u and y changes
# (1) emerging VOCs, changes in susceptibility and clinical fraction, immune 
# evading characteristics
# (2) different vaccine efficacy values intput at different time steps due to 
# composition of vaccines
update_u_y <- function(para = NULL,
                       # group (1) changes
                       # date_switch marks the introduction of new VOCs
                       # check function to make sure that date_switch and rc_x
                       # have the same size
                       # rc - relative changes
                       date_switch = c("2021-01-15", "2021-04-15", "2021-12-15"),
                       rc_u = c(1, 1.5, 0.5), # relative changes in u
                       rc_y = c(1, 0.5, 0.5), # relative changes in y
                       # if ve can be specific to VOCs in relation to the 
                       # wildtype, rc_ve will be all 1s
                       rc_ve = c(1, 0.9, 0.7), # relative changes in evasiveness (infection part)
                       # group (2) changes 
                       efficacy_baseline = NULL, # vaccine efficacy
                       # will need something to check 
                       efficacy_weights = efficacy_weights_test # as a result of different vaccine composition
){
  # debug
  # date_switch = c("2021-01-15", "2021-04-15", "2021-12-15")
  # rc_u = c(1, 1.5, 0.5)
  # rc_y = c(1, 0.5, 0.5)
  # rc_ve = c(1, 0.9, 0.7)
  # efficacy_baseline = efficacy_all
  # efficacy_weights = efficacy_weights_test

  date_range <- c(lubridate::ymd(para$date0),
                  date_switch,
                  lubridate::ymd(para$date0) + para$time1)
  
  t_range <- c(para$time0:para$time1)
  
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
  para$pop[[1]]$ur    <- (1 - efficacy_baseline$r_i_o[1]) * para$pop[[1]]$u
  para$pop[[1]]$uvr_l <- (1 - efficacy_baseline |> 
                            filter(protection_level_label == "l") |> 
                            pull(vr_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uvr_m <- (1 - efficacy_baseline |> 
                            filter(protection_level_label == "m") |> 
                            pull(vr_i_o)) * para$pop[[1]]$u
  para$pop[[1]]$uvr_h <- (1 - efficacy_baseline |> 
                            filter(protection_level_label == "h") |> 
                            pull(vr_i_o)) * para$pop[[1]]$u
  
  # create modifier table
  CJ(date = seq(lubridate::ymd(para$date0),
                lubridate::ymd(para$date0) + para$time1,
                by = "day")) |> 
    mutate(phase = as.numeric(NA)) -> modifier
  
  for(i in 1:length(date_range)){
    modifier[modifier$date >= date_range[i] & 
               modifier$date <= date_range[i+1],"phase"] <- i
    modifier[modifier$date >= date_range[i] & 
               modifier$date <= date_range[i+1],"rc_u_prod"] <- prod(c(1,rc_u)[1:i])
    modifier[modifier$date >= date_range[i] & 
               modifier$date <= date_range[i+1],"rc_y_prod"] <- prod(c(1,rc_y)[1:i])
    modifier[modifier$date >= date_range[i] & 
               modifier$date <= date_range[i+1],"rc_ve_prod"] <- prod(c(1,rc_ve)[1:i])
  }
  

  
  # VEs against infection and disease are implemented over "compartments"
  # VEs against severe, critical and mortality cases are implemented over "processes" 
  # Everything above infection in terms of outcome will need to use conditional
  # probability, because the infection step has already occurred to reach this 
  # endpoint
  
  # in the context of this model, disease preventing = clinical preventing
  modifier |> 
    left_join(efficacy_weights,
              by = "date") |> 
    mutate(u_scaler     = rc_u_prod,
           uv_l_scaler  = rc_u_prod*(1 - efficacy_baseline$v_i_o[1]*weights_v_i_l*rc_ve_prod),
           uv_m_scaler  = rc_u_prod*(1 - efficacy_baseline$v_i_o[2]*weights_v_i_m*rc_ve_prod),
           uv_h_scaler  = rc_u_prod*(1 - efficacy_baseline$v_i_o[3]*weights_v_i_h*rc_ve_prod),
           ur_scaler    = rc_u_prod,
           uvr_l_scaler = rc_u_prod*(1 - efficacy_baseline$v_i_o[1]*weights_v_i_l*rc_ve_prod),
           uvr_m_scaler = rc_u_prod*(1 - efficacy_baseline$v_i_o[2]*weights_v_i_m*rc_ve_prod),
           uvr_h_scaler = rc_u_prod*(1 - efficacy_baseline$v_i_o[3]*weights_v_i_h*rc_ve_prod),
           yv_l_scaler  = rc_y_prod*(1 - efficacy_baseline$v_d_condition[1]*weights_v_d_l*rc_ve_prod),
           yv_m_scaler  = rc_y_prod*(1 - efficacy_baseline$v_d_condition[2]*weights_v_d_m*rc_ve_prod),
           yv_h_scaler  = rc_y_prod*(1 - efficacy_baseline$v_d_condition[3]*weights_v_d_h*rc_ve_prod)) -> modifier
  
  modifier |> 
    select(ends_with("scaler")) |> 
    pivot_longer(targets$scaler_label) |> 
    group_by(name) |> summarise(value = min(value)) |> 
    pull(value) |> (function(y) y > 0)() |> all() -> test_range
  
  testthat::expect(test_range,
                   failure_message = "zero scaler values generated. 
                   please double check all rc_xx variables.")
  
  for(i in seq_len(nrow(targets))){
    tmp <-   modifier |> 
      pull(targets$scaler_label[i]) |> 
      split(seq(nrow(modifier))) |> 
      map(rep, n_age_groups) |> 
      map(unname)
    
    para$schedule[[targets$scaler_label[i]]] <-  list(
      parameter = targets$variable_label[i],
      pops = numeric(),
      mode = "multiply",
      values = tmp,
      times = t_range
    )
    
    rm(tmp)
  }
  
  # return results
  return(para)
}

# this function is reduced now to only make changes to only modify the health
# system burden processes changed as a result of emerging VOC

emerge_VOC_burden <- function(
    para = NULL,
    # could add another variable for ve_against >= severe outcomes to change by
    # VOC stages, to-do
    rc_severity = c(1, 1.5,1.5), # relative change in ihr and ifr
    efficacy_baseline = NULL){
  
  # debug
  # para = params
  # rc_severity = c(1, 1.5,1.5)
  # efficacy_baseline = ve_all
  
  # set up
  # the number of VOC 
  # all rc_xx variables need to have the same length
  n_voc = length(rc_severity)
  compartments_E <- c("newE", "newEv_l", "newEv_m", "newEv_h")
  
  # generate death processes
  generate_death_processes <- function(source_compartment = NULL,
                                       voc_index = NULL){
    multiplier1 <- prod(rc_severity[1:voc_index])
    if(source_compartment == "newE") multiplier2 <- 1
    if(source_compartment == "newEv_l") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[1]
    if(source_compartment == "newEv_m") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[2]
    if(source_compartment == "newEv_h") multiplier2 <- 1 - efficacy_baseline$v_mort_condition[3]
    tmp_var <- paste0("death_voc", voc_index)
    
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
    multiplier1 <- prod(rc_severity[1:voc_index])
    
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
    
    tmp_var_to_severe <- paste0("to_severe_voc", voc_index)
    tmp_var_to_critical <- paste0("to_critical_voc", voc_index)
    tmp_var_severe <- paste0("severe_voc", voc_index)
    tmp_var_critical <- paste0("critical_voc", voc_index)
    
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
  
  for(i in 1:n_voc){
    tmp_var_to_severe <- paste0("to_severe_voc", i)
    tmp_var_to_critical <- paste0("to_critical_voc", i)
    tmp_var_severe <- paste0("severe_voc", i)
    tmp_var_critical <- paste0("critical_voc", i)
    
    list(  
      cm_multinom_process(tmp_var_to_severe, 
                          data.frame(rep(1,16)) |> setNames(tmp_var_severe),
                          delays = data.frame(delay_2hosp)|> 
                            setNames(tmp_var_severe),   report = "ip"),
      
      cm_multinom_process(tmp_var_to_critical,
                          data.frame(rep(1,16))|> setNames(tmp_var_critical),
                          delays = data.frame(delay_2hosp_critical)|> 
                            setNames(tmp_var_critical),   report = "ip")
    ) -> intermediate_processes_agg[[i]]
  }
  
  intermediate_processes_agg |> purrr::flatten() -> intermediate_processes_agg
  to_attach <- c(death_processes, 
                 intermediate_processes,
                 intermediate_processes_agg)
  
  para$processes <- c(para$processes, to_attach)

  return(para)
}

# helper code file for this function is 0_4_1 and 0_4_2
# these code files help you prepare for these input objects: owid_vac, 
# primary_allocation_plan
vaccinate_primary <- function(para = NULL,
                              vac_data = owid_vac,
                              values = primary_allocation_plan
                              ){
  
  # debug
  # para <- params
  # vac_data = owid_vac,
  # values = primary_allocation_plan
  # 
  require(lurbidate)
  n_age_groups <- length(para$pop[[1]]$size)
  date_start <- ymd(para$date0)
  date_end <- date_start + para$time1
  data.frame(date = seq(date_start, date_end, by = "day")) |> 
    mutate(t = 0:para$time1,
           empirical = date %in% (vac_data$date)) |> 
    filter(empirical == T) |> 
    pull(t) -> tmp_times
  
  c(0, tmp_times, max(tmp_times)+1) -> tmp_times
  c(list(rep(0,16)), values, list(rep(0,16))) -> tmp_allocation
  
  testthat::expect_equal(length(tmp_times), length(tmp_allocation))

  para$schedule[["primary_course"]] <- list(
    parameter = "v_p",
    pops = numeric(),
    mode = "assign",
    values = tmp_allocation,
    times = tmp_times
  )
  return(para)
}

# multiple booster campaigns is it a one time thing?
# duration of interval; the start of the first booster campaign; age prioritisation
# booster vaccine characteristics; 
# vaccinate_booster <- function(para = NULL,
#                               vac_data = owid_vac,
#                               # this is paused time
#                               program_interval = 30*6, #default set to 6 months
#                               # should this be age-specific as well?
#                               uptake_by_existing = 0.9, 
#                               # age-specific variables that defines the 
#                               # prioritisation, the numbers are essentially just
#                               # rankings; NA = not boosted
#                               # this is based on history, based on owid_vac
#                               prioritisation_initial = c(rep(NA, 4), rep(1,12)),
#                               # this is future policy
#                               prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
#                               boosters_daily = 300000
#                               ){
# 
#   require(lubridate)
#   # debug
#   # para <- params
#   # vac_data = owid_vac
#   # program_interval = 30*6
#   # uptake_by_existing = 0.9
#   # prioritisation_initial = c(rep(NA, 4), rep(1,12))
#   # prioritisation_followup = c(NA,rep(2,11),rep(1,4))
#   # boosters_daily = 300000
#   
#   if(length(uptake_by_existing) == 1) uptake_by_existing_tmp <- rep(uptake_by_existing, 16)
#   testthat::expect_length(uptake_by_existing, 16)
# 
#   time_range <- data.frame(date = seq(ymd(para$date0),
#                                       ymd(para$date0) + (para$time1),
#                                       by = "day")) |> 
#     rownames_to_column(var = "t") |> 
#     mutate(t = as.numeric(t))
#   
#   tmp_allocation <- tmp_times <- list()
#   
#   # initial boosting programmes
#   # with owid data
#   proportions_allocated_initial <- para$pop[[1]]$size/sum(para$pop[[1]]$size)
#   n_age_groups <- length(proportions_allocated_initial)
#   # we want the initial stage to not divide by stage and target all adults
#   proportions_allocated_initial_rescaled <- (prioritisation_initial*proportions_allocated_initial)/sum(prioritisation_initial*proportions_allocated_initial, na.rm = T)
#   proportions_allocated_initial_rescaled[is.na(proportions_allocated_initial_rescaled)] <- 0
#   
#   vac_data |> 
#     select(total_boosters_daily) %>%
#     split(seq(nrow(.))) |> 
#     map(unlist) |> 
#     map(.f = function(x) x*proportions_allocated_initial_rescaled) |> 
#     setNames(NULL) -> tmp_allocation[["initial"]]
#   
#   vac_data |> 
#     mutate(date = lubridate::date(date)) |> 
#     left_join(time_range, by = "date") |> 
#     dplyr::select(date, t) |> 
#     pull(t) -> tmp_times[["initial"]]
#   
#   # follow-up campaigns
#   # children coverage = 0.787; adolescent coverage = 0.812; adult coverage = 0.813
#   data.frame(
#     prioritisation_followup = prioritisation_followup,
#     pop = para$pop[[1]]$size,
#     cov_primary = c(NA, 0.787, rep(0.812, 2),
#                     rep(0.813, 12))
#   ) |>
#     mutate(
#       cov_followup = uptake_by_existing * cov_primary,
#       cov_followup_doses = cov_followup * pop,
#       cov_followup_doses_all = sum(cov_followup_doses, na.rm = T),
#       campaign_duration = round(cov_followup_doses_all / boosters_daily)
#     ) |>
#     group_by(prioritisation_followup) |>
#     mutate(cov_followup_doses_bygroup = sum(cov_followup_doses, na.rm = T)) |> ungroup() |>
#     mutate(campaign_duration_bygroup = round(cov_followup_doses_bygroup /
#                                                boosters_daily)) |>
#     dplyr::select(prioritisation_followup, campaign_duration_bygroup) |> unique() |>
#     filter(!is.na(prioritisation_followup)) |>
#     arrange(prioritisation_followup) -> follow_up_order
#   
#   follow_up_schedule <- c(program_interval, follow_up_order$campaign_duration_bygroup)
#   follow_up_unit <- sum(follow_up_schedule)
#   t_empirical_end <- (time_range |> filter(date ==  range(vac_data$date)[2]) |> pull(t))
#   date_booster_start <-  vac_data |> filter(boosters_daily > 0) |> pull(date) |> min()
# 
#   time_range |>
#     mutate(
#       cycle = (t - t_empirical_end) / follow_up_unit,
#       cycle = floor(cycle),
#       t_within = t - t_empirical_end - follow_up_unit * cycle
#     ) |> 
#     mutate(
#       vaccination_phase = case_when(
#         date >= date_booster_start & t <= t_empirical_end ~ "booster_initial",
#         cycle >= 0 &
#           t_within <= follow_up_schedule[1] ~ "pause",
#         cycle >= 0 &
#           t_within > follow_up_schedule[1] &
#           t_within <= sum(follow_up_schedule[1:2]) ~ "booster_OA",
#         cycle >= 0 &
#           t_within > sum(follow_up_schedule[1:2]) &
#           t_within <= sum(follow_up_schedule[1:3]) ~ "booster_all"
#       )
#     ) |> 
#     dplyr::filter(cycle >= 0) -> phase_def
# 
#   phase_def |> 
#     group_by(vaccination_phase, cycle) |> 
#     mutate(start = min(t_within)) |> 
#     filter(start == t_within,
#            vaccination_phase != "booster_initial") -> phase_list
#   
#   # proportions OA rescaled
#   proportions_allocation_OA_rescaled <- proportions_allocated_initial
#   proportions_allocation_OA_rescaled[1:12] <- 0
#   proportions_allocation_OA_rescaled <- proportions_allocation_OA_rescaled/sum(proportions_allocation_OA_rescaled)
#   
#   # proportions A rescaled
#   proportions_allocation_all_rescaled <- proportions_allocated_initial
#   proportions_allocation_all_rescaled[1] <- 0
#   proportions_allocation_all_rescaled[13:16] <- 0
#   proportions_allocation_all_rescaled <- proportions_allocation_all_rescaled/sum(proportions_allocation_all_rescaled)
#   
#   allocation_byphase <- list(pause = rep(0,16),
#                              booster_OA = c(rep(0,12), rep(1,4))*proportions_allocation_OA_rescaled,
#                              booster_all = c(rep(0,1),  rep(1,11), rep(0,4))*proportions_allocation_all_rescaled)
#   
#   tmp_times[["campaign"]] <- phase_list$t
#   tmp_allocation[["campaign"]] <- list()
#   for(j in 1:nrow(phase_list)){
#     if(phase_list$vaccination_phase[j] == "pause") tmp_allocation[["campaign"]][[j]] <- allocation_byphase$pause*boosters_daily
#     if(phase_list$vaccination_phase[j] == "booster_OA") tmp_allocation[["campaign"]][[j]] <- allocation_byphase$booster_OA*boosters_daily
#     if(phase_list$vaccination_phase[j] == "booster_all") tmp_allocation[["campaign"]][[j]] <- allocation_byphase$booster_all*boosters_daily
# 
#   }
#   
#   tmp_times_move <- c(0,unlist(tmp_times) |> array())
#   tmp_values_move <- c(list(rep(0,16)),
#                        tmp_allocation$initial,
#                        tmp_allocation$campaign) |> setNames(NULL)
# 
#   testthat::expect_equal(length(tmp_times_move),
#                          length(tmp_values_move))
# 
#   para$schedule[["booster"]] <- list(
#     parameter = "v_b",
#     pops = numeric(),
#     mode = "assign",
#     values = tmp_values_move,
#     times = tmp_times_move
#   )
#   
#   return(para)
# }

source("code/0_1_1_AnnualProgram.R")

cm_multinom_process <- function(
    src, outcomes, delays,
    report = ""
) {
  if ("null" %in% names(outcomes)) {
    if (length(report) != length(outcomes)) report <- rep(report, length(outcomes))
    report[which(names(outcomes)=="null")] <- ""
    if (!("null" %in% names(delays))) {
      delays$null <- c(1, rep(0, length(delays[[1]])-1))
    }
  } else if (!all(rowSums(outcomes)==1)) {
    report <- c(rep(report, length(outcomes)), "")
    outcomes$null <- 1-rowSums(outcomes)
    delays$null <- c(1, rep(0, length(delays[[1]])-1))
  }
  nrow <- length(outcomes)
  list(
    source = src, type="multinomial", names=names(outcomes), report = report,
    prob = t(as.matrix(outcomes)), delays = t(as.matrix(delays))
  )
}

check_vaccination_program <- function(type = "booster", # or primary_course
                                      para = NULL){
  # para <- params
  # type = "booster"
  # type = "primary_course"
  para$schedule[[type]]$values |>
    map(data.frame) |> map(t) |> map(data.frame) |>
    bind_rows() |> set_rownames(NULL) |>
    mutate(t = para$schedule[[type]]$t) |>
    full_join(data.frame(t = seq(para$time0, para$time1)) |>
                mutate(date = lubridate::ymd(para$date0) + as.numeric(t)),
              by = "t") |>
    arrange(date) |>
    mutate_at(vars(starts_with("X")),
              imputeTS::na_locf) |>
    pivot_longer(starts_with("X")) |>
    mutate(name = factor(name,
                         levels = paste0("X", 1:16))) -> p_table
  
  year_lims <- paste0(c(p_table$date |> lubridate::year() |> min(na.rm = T),
                 p_table$date |> lubridate::year() |> max(na.rm = T)),"-01-01") |> 
    lubridate::ymd()
  
  p_table |> filter(name == "X7", date >= "2023-01-01") |> pull(value) |> unique()
  
  p_table |> 
    ggplot(aes(x = date, y = value)) +
    geom_line() +
    facet_wrap(~name) +
    geom_vline(xintercept = seq(year_lims[1],
                                year_lims[2],
                                by = "year")) -> p
  
  return(p)
  
}
