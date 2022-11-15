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

gen_country_basics <- function(country,
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2022-12-31",
                               processes = NULL,
                               contact = contact_schedule,
                               period_wn  = 3*365, # duration, waning of natural immunity
                               period_wv_m2l = 1*365, # duration, waning from medium to low levels vaccine induced 
                               period_wv_h2m = 1*365, # duration, waning from medium to low levels vaccine induced 
                               prob_v_p_2l = 0.5,
                               prob_v_p_2m = 0.3,
                               prob_v_b_l2m = 0.5,
                               # this needs to be a three number array
                               ve_inf = 0.7, # probability reduction of breakthrough due to infection
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
    
    # update all copies of u and y to begin with
    # need to fix this as well
    
    # ve against infection by vaccine induced immunity only
    para$pop[[i]]$uv_l  <- para$pop[[i]]$u
    para$pop[[i]]$uv_m  <- para$pop[[i]]$u
    para$pop[[i]]$uv_h  <- para$pop[[i]]$u
    # ve against infection by infection only
    para$pop[[i]]$ur    <- (1 - ve_inf)*para$pop[[i]]$u
    # ve against infection by hybrid immunity
    para$pop[[i]]$uvr_l <- (1 - ve_inf)*para$pop[[i]]$u
    para$pop[[i]]$uvr_m <- (1 - ve_inf)*para$pop[[i]]$u
    para$pop[[i]]$uvr_h <- (1 - ve_inf)*para$pop[[i]]$u
    
    para$pop[[i]]$yv_l <- para$pop[[i]]$y
    para$pop[[i]]$yv_m <- para$pop[[i]]$y
    para$pop[[i]]$yv_h <- para$pop[[i]]$y

    # natural waning
    para$pop[[i]]$wn <- rep((1/period_wn), n_age_groups)
    
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
  
  para$schedule[["mobility"]] = list(
    parameter = "contact",
    pops = numeric(),
    mode = "assign",
    values = split(c_tmp[,3:6],
                   seq(nrow(c_tmp))) %>%
      map(unlist) %>%
      map(as.vector) %>%
      unname,
    times = 1:nrow(c_tmp))
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wn <- rep(1/period_wn, n_age_groups)
  para$pop[[1]]$wv_m2l <-  rep(1/period_wv_m2l, n_age_groups)
  para$pop[[1]]$wv_h2m <-  rep(1/period_wv_h2m, n_age_groups)
  
  return(para)
}


# # generate some test data of what efficacy weights should look like
# CJ(date = seq(ymd("2019-12-01"),
#               ymd("2024-01-01"),
#               by = "day")) %>%
#   mutate(weights_ve_i_l = rnorm(nrow(.),1,0.05),
#          weights_ve_i_m = rnorm(nrow(.),1,0.05),
#          weights_ve_d_l = rnorm(nrow(.),1,0.05),
#          weights_ve_d_m = rnorm(nrow(.),1,0.05)) -> efficacy_weights_test
# write_rds(efficacy_weights_test,
#           "data/intermediate/efficacy_weights_test.rds"
#           )
efficacy_weights_test <- read_rds("data/intermediate/efficacy_weights_test.rds")
efficacy_weights_test |> 
  mutate_at(vars(starts_with("weights")), function(x) x <- 1) -> efficacy_weights_one

# this function will update u and y
# there's two sources for u and y changes
# (1) VOC emergences related to changes in susceptibility and clinical fraction
# (2) different vaccine efficacy values intput at different time steps due to 
# composition of vaccines
update_u_y <- function(para = NULL,
                       date_switch = c("2021-01-15", "2021-04-15", "2021-12-15"),
                       rc_u = c(1, 1.5, 0.5), # relative changes in u
                       rc_y = c(1, 0.5, 0.5), # relative changes in y
                       rc_ve = c(1, 0.9, 0.7), # relative changes in evasiveness (infection part)
                       # relative chanves in 
                       efficacy_baseline = NULL, # vaccine efficacy 
                       efficacy_weights = efficacy_weights_test # as a result of different vaccine composition
){
  # debug
  # date_switch = c("2021-01-15", "2021-04-15", "2021-12-15")
  # rc_u = c(1, 1.5, 0.5)
  # rc_y = c(1, 0.5, 0.5)
  # rc_ve = c(1, 0.9, 0.7)
  # efficacy_baseline = ve_az
  # efficacy_weights = efficacy_weights_test

  date_range <- c(lubridate::ymd(para$date0),
                  date_switch,
                  lubridate::ymd(para$date0) + para$time1)
  t_range <- c(para$time0:para$time1)
  targets <- data.frame(scaler_label = c("uv_l_scaler", "uv_m_scaler",
                                         "uvr_l_scaler", "uvr_m_scaler",
                                         "yv_l_scaler", "yv_m_scaler"),
                        variable_label = c("uv_l", "uv_m",
                                           "uvr_l", "uvr_m",
                                           "yv_l", "yv_m"))
  n_age_groups <- para$pop[[1]]$n_groups
  
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
  
  modifier |> 
    left_join(efficacy_weights,
              by = "date") |> 
    mutate(uv_l_scaler = rc_u_prod*(1 - efficacy_baseline$ve_i_o[1]*weights_ve_i_l*rc_ve_prod),
           uv_m_scaler = rc_u_prod*(1 - efficacy_baseline$ve_i_o[2]*weights_ve_i_m*rc_ve_prod),
           uvr_l_scaler = rc_u_prod*(1 - efficacy_baseline$ve_i_o[1]*weights_ve_i_l*rc_ve_prod),
           uvr_m_scaler = rc_u_prod*(1 - efficacy_baseline$ve_i_o[2]*weights_ve_i_m*rc_ve_prod),
           yv_l_scaler = rc_y_prod*(1 - efficacy_baseline$ve_d[1]*weights_ve_d_l*rc_ve_prod),
           yv_m_scaler = rc_y_prod*(1 - efficacy_baseline$ve_d[2]*weights_ve_d_m*rc_ve_prod)) -> modifier

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
    rc_severity = c(1, 1.5,1.5), # relative change in ihr and ifr
    efficacy_baseline = NULL){
  
  # debug
  # para = params
  # rc_severity = c(1, 1.5,1.5)
  # efficacy_baseline = ve_az
  
  # set up
  n_voc = length(rc_severity)
  compartments_E <- c("newE", "newEv_l", "newEv_m", "newEv_h")
  
  # generate death processes
  generate_death_processes <- function(source_compartment = NULL,
                                       voc_index = NULL){
    multiplier1 <- prod(rc_severity[1:voc_index])
    if(source_compartment == "newE") multiplier2 <- 1
    if(source_compartment == "newEv_l") multiplier2 <- 1 - efficacy_baseline$ve_mort_condition[1]
    if(source_compartment == "newEv_m") multiplier2 <- 1 - efficacy_baseline$ve_mort_condition[2]
    if(source_compartment == "newEv_h") multiplier2 <- 1 - efficacy_baseline$ve_mort_condition[3]
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
      multiplier2_severe <- 1 - efficacy_baseline$ve_severe_condition[1]
      multiplier2_critical <- 1 - efficacy_baseline$ve_critical_condition[1]
      
    }
    
    if(source_compartment == "newEv_m") {
      multiplier2_severe <- 1 - efficacy_baseline$ve_severe_condition[2]
      multiplier2_critical <- 1 - efficacy_baseline$ve_critical_condition[2]
    }
    
    if(source_compartment == "newEv_h") {
      multiplier2_severe <- 1 - efficacy_baseline$ve_severe_condition[3]
      multiplier2_critical <- 1 - efficacy_baseline$ve_critical_condition[3]
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

vaccinate_primary <- function(para = NULL){
  n_age_groups <- length(para$pop[[1]]$size)
  
  tmp_time <- c(0,
                as.numeric(ymd("2021-06-15") - ymd(para$date0)),
                as.numeric(ymd("2021-11-15") - ymd(para$date0)),
                para$time1)
  rep_time <- diff(tmp_time)
  rep_time[1] <-   rep_time[1] + 1
  
  list_tmp <- list()
  for(i in 1:1096) {
    list_tmp[[i]] <- rep(0, n_age_groups)
    if(i >= tmp_time[2]+2 & i <= tmp_time[3]){
    list_tmp[[i]] <- rep(c(0,7500), c(4, n_age_groups - 4))
    }
  }
  
  para$schedule[["primary_course"]] <- list(
    parameter = "v_p",
    pops = numeric(),
    mode = "assign",
    values = list_tmp,
    times = 0:1095
  )
  return(para)
}

# multiple booster campaigns is it a one time thing?
# duration of interval; the start of the first booster campaign; age prioritisation
# booster vaccine characteristics; 
vaccinate_booster <- function(para = NULL,
                              program_start = NULL,
                              program_end = NULL){
  n_age_groups <- length(para$pop[[1]]$size)
  tmp_time <- c(0,
                as.numeric(ymd(program_start) - ymd(para$date0)),
                as.numeric(ymd(program_end) - ymd(para$date0)))
  
  para$schedule[["booster"]] <- list(
    parameter = "v_b",
    pops = numeric(),
    mode = "assign",
    values = list(rep(0, n_age_groups),
                  rep(c(0,10000), c(4, n_age_groups - 4)),
                  rep(0, n_age_groups)),
    times = tmp_time
  )
  return(para)
}
 