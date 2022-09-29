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
                               period_wv_ml = 1*365, # duration, waning from medium to low levels vaccine induced 
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
    
    # scale u (susceptibility) to achieve desired R0
    current_R0 = cm_calc_R0(para, i); # calculate R0 in population i of params
    para$pop[[i]]$u = para$pop[[i]]$u * R0_assumed / current_R0
    
    # natural waning
    para$pop[[i]]$wn <- rep((1/period_wn), 
                            n_age_groups)
    
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
  para$pop[[1]]$wv_ml <-  rep(1/period_wv_ml, n_age_groups)
  
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
update_u_y <- function(para,
                       date_switch = c("2021-01-15", "2021-04-15", "2021-12-15"),
                       rc_u = c(1, 1.5, 0.5), # relative changes in u
                       rc_y = c(1, 0.5, 0.5), # relative changes in y
                       rc_ve = c(1, 0.9, 0.7), # relative changes in ve
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
                                         "yv_l_scaler", "yv_m_scaler"),
                        variable_label = c("uv_l", "uv_m",
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

# emerge_VOC_burden <- function(
#     para = NULL,
#     rc_severity = c(1, 1.5,1.5), # relative change in ihr and ifr
#     VE = NULL,
# ){
#   
#   burden_processes_new <- 
#     # mortality cases
#     list(
#       cm_multinom_process("E",       
#                           data.frame(death_voc1 = P.death*rc_severity[1]),                   
#                           delays = data.frame(death_voc1 = delay_2death), 
#                           report = "o"),
#       cm_multinom_process("E",       
#                           data.frame(death_voc2 = P.death*prod(rc_severity[1:2])),                   
#                           delays = data.frame(death_voc2 = delay_2death), 
#                           report = "o"),
#       cm_multinom_process("E",       
#                           data.frame(death_voc3 = P.death*prod(rc_severity[1:3])),                   
#                           delays = data.frame(death_voc3 = delay_2death), 
#                           report = "o"),
#       
#       cm_multinom_process("Ev",      
#                           data.frame(death_voc1 = P.death*(1-VE$ve_mort_condition[1])*(rc_severity[1])), 
#                           delays = data.frame(death_voc1 = delay_2death), 
#                           report = "o"),
#       cm_multinom_process("Ev",      
#                           data.frame(death_voc2 = P.death*(1-VE$ve_mort_condition[1])*prod(rc_severity[1:2])), 
#                           delays = data.frame(death_voc2 = delay_2death), 
#                           report = "o"),
#       cm_multinom_process("Ev",      
#                           data.frame(death_voc3 = P.death*(1-VE$ve_mort_condition[1])*prod(rc_severity[1:3])), 
#                           delays = data.frame(death_voc3 = delay_2death), 
#                           report = "o"),
#       
#       cm_multinom_process("Ev2",     
#                           data.frame(death_voc1 = P.death*(1-VE$ve_mort_condition[2])*(rc_severity[1])), 
#                           delays = data.frame(death_voc1 = delay_2death), 
#                           report = "o"),
#       cm_multinom_process("Ev2",     
#                           data.frame(death_voc2 = P.death*(1-VE$ve_mort_condition[2])*prod(rc_severity[1:2])), 
#                           delays = data.frame(death_voc2 = delay_2death),
#                           report = "o"),
#       cm_multinom_process("Ev2",
#                           data.frame(death_voc3 = P.death*(1-VE$ve_mort_condition[2])*prod(rc_severity[1:3])), 
#                           delays = data.frame(death_voc3 = delay_2death), 
#                           report = "o"),
#       
#       # severe cases
#       cm_multinom_process("E",       
#                           data.frame(to_severe_voc1 = P.severe*(rc_severity[1]),
#                                      to_critical_voc1 = P.critical*(rc_severity[1])),                  
#                           delays = data.frame(to_severe_voc1 = delay_2severe,
#                                               to_critical_voc1 = delay_2severe)),
#       cm_multinom_process("E",       
#                           data.frame(to_severe_voc2 = P.severe*prod(rc_severity[1:2]),
#                                      to_critical_voc2 = P.critical*prod(rc_severity[1:2])),                  
#                           delays = data.frame(to_severe_voc2 = delay_2severe,
#                                               to_critical_voc2 = delay_2severe)),
#       cm_multinom_process("E",       
#                           data.frame(to_severe_voc3 = P.severe*prod(rc_severity[1:3]),
#                                      to_critical_voc3 = P.critical*prod(rc_severity[1:3])),                  
#                           delays = data.frame(to_severe_voc3 = delay_2severe,
#                                               to_critical_voc3 = delay_2severe)),
#       
#       cm_multinom_process("Ev",      
#                           data.frame(to_severe_voc1 = P.severe*(1-VE$ve_severe_condition[1])*(rc_severity[1]),
#                                      to_critical_voc1 = P.critical*(1-VE$ve_critical_condition[1])*(rc_severity[1])),   
#                           delays = data.frame(to_severe_voc1 = delay_2severe,
#                                               to_critical_voc1 = delay_2severe)),
#       cm_multinom_process("Ev",      
#                           data.frame(to_severe_voc2 = P.severe*(1-VE$ve_severe_condition[1])*prod(rc_severity[1:2]),
#                                      to_critical_voc2 = P.critical*(1-VE$ve_critical_condition[1])*prod(rc_severity[1:2])),   
#                           delays = data.frame(to_severe_voc2 = delay_2severe,
#                                               to_critical_voc2 = delay_2severe)),
#       cm_multinom_process("Ev",      
#                           data.frame(to_severe_voc3 = P.severe*(1-VE$ve_severe_condition[1])*prod(rc_severity[1:3]),
#                                      to_critical_voc3 = P.critical*(1-VE$ve_critical_condition[1])*prod(rc_severity[1:3])),   
#                           delays = data.frame(to_severe_voc3 = delay_2severe,
#                                               to_critical_voc3 = delay_2severe)),
#       
#       cm_multinom_process("Ev2",     
#                           data.frame(to_severe_voc1 = P.severe*(1-VE$ve_severe_condition[2])*(rc_severity[1]),
#                                      to_critical_voc1 = P.critical*(1-VE$ve_critical_condition[2])*(rc_severity[1])),   
#                           delays = data.frame(to_severe_voc1 = delay_2severe,
#                                               to_critical_voc1 = delay_2severe)),
#       cm_multinom_process("Ev2",     
#                           data.frame(to_severe_voc2 = P.severe*(1-VE$ve_severe_condition[2])*prod(rc_severity[1:2]),
#                                      to_critical_voc2 = P.critical*(1-VE$ve_critical_condition[2])*prod(rc_severity[1:2])),   
#                           delays = data.frame(to_severe_voc2 = delay_2severe,
#                                               to_critical_voc2 = delay_2severe)),
#       cm_multinom_process("Ev2",     
#                           data.frame(to_severe_voc3 = P.severe*(1-VE$ve_severe_condition[2])*prod(rc_severity[1:3]),
#                                      to_critical_voc3 = P.critical*(1-VE$ve_critical_condition[2])*prod(rc_severity[1:3])),   
#                           delays = data.frame(to_severe_voc3 = delay_2severe,
#                                               to_critical_voc3 = delay_2severe)),
#       
#       cm_multinom_process("to_severe_voc1", data.frame(severe_voc1 = rep(1,16)),                          
#                           delays = data.frame(severe_voc1 = delay_2hosp),   report = "ip"),
#       cm_multinom_process("to_severe_voc2", data.frame(severe_voc2 = rep(1,16)),                          
#                           delays = data.frame(severe_voc2 = delay_2hosp),   report = "ip"),
#       cm_multinom_process("to_severe_voc3", data.frame(severe_voc3 = rep(1,16)),                          
#                           delays = data.frame(severe_voc3 = delay_2hosp),   report = "ip"),
#       
#       cm_multinom_process("to_critical_voc1", data.frame(critical_voc1 = rep(1,16)),                          
#                           delays = data.frame(critical_voc1 = delay_2hosp_critical),   report = "ip"),
#       cm_multinom_process("to_critical_voc2", data.frame(critical_voc2 = rep(1,16)),                          
#                           delays = data.frame(critical_voc2 = delay_2hosp_critical),   report = "ip"),
#       cm_multinom_process("to_critical_voc3", data.frame(critical_voc3 = rep(1,16)),                          
#                           delays = data.frame(critical_voc3 = delay_2hosp_critical),   report = "ip")
#     )
#   
#   burden_updated <- c(para$processes, burden_processes_new)
#   para$processes <- burden_updated
#   
#   return(para)
# }

 