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
                               waning_nat = 52*7*3,
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2022-12-31",
                               processes = NULL,
                               deterministic = TRUE){
  
  require(countrycode)
  
  iso3c_tmp = countrycode(country, "country.name", "iso3c")
  
  c_tmp = contact_schedule %>% 
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
    para$pop[[i]]$wn <- rep((1/waning_nat), n_age_groups)
    
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
  
  return(para)
}

update_vac_char <- function(para,
                            ve_i_l  = NULL, # infection blocking VE among V_l
                            ve_i_m  = NULL, # infection blocking VE among V_m
                            ve_d_l  = NULL, # disease reducing VE among E_l
                            ve_d_m  = NULL, # disease reducing VE among E_m 
                            period_wn  = NULL, # duration, waning of natural immunity
                            period_wv_ml = NULL # duration, waning from medium to low levels vaccine induced 
){
  
  n_age <- length(para$pop[[1]]$size)
  
  # key parameters on infection and disease blocking mechanisms
  para$pop[[1]]$uv_l  <- para$pop[[1]]$u * (1 - ve_i_l)
  para$pop[[1]]$uv_m <- para$pop[[1]]$u * (ve_i_m)
  para$pop[[1]]$yv_l  <- para$pop[[1]]$y * (1 - ve_d_l)
  para$pop[[1]]$yv_m <- para$pop[[1]]$y * (1 - ve_d_m)
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wn <- rep(1/period_wn, n_age)
  para$pop[[1]]$wv_ml <-  rep(1/period_wv_ml, n_age)
  
  # return results
  return(para)
}
 