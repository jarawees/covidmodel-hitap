
critical2 <- 0
picu_cocin_func <- function(age) {
  x <- c(-0.1309118, 0, 17.2398874, 65.7016492, 100)
  y <- c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
  p <- splinefun(x, y)(age)
  exp(p) / (1 + exp(p))
}
picu_cocin <- picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin <- 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100
# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje <- exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85))
# Amalgamate probabilities
probabilities <- data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_levin, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities <- probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.critical <- probabilities[, ihr * picu]
P.severe <- probabilities[, ihr * (1 - picu)]
P.death <- probabilities[, ifr]
P.hosp <- P.critical + P.severe

delay_2death <- cm_delay_gamma(26, 5, 60, 0.25)$p
delay_2severe <- cm_delay_gamma(8.5, 5, 60, 0.25)$p
delay_2hosp <- cm_delay_gamma(14.6, 5, 60, 0.25)$p
delay_2hosp_critical <- cm_delay_gamma(15.6, 5, 60, 0.25)$p

gen_burden_processes <- function(VE){
  tmp <- list(
    cm_multinom_process("E",       
                        data.frame(death = P.death),                   
                        delays = data.frame(death = delay_2death), report = "o"),
    cm_multinom_process("Ev_l",      
                        data.frame(death = P.death*(1-VE$ve_mort_condition[1])), 
                        delays = data.frame(death = delay_2death), report = "o"),
    cm_multinom_process("Ev_m",     
                        data.frame(death = P.death*(1-VE$ve_mort_condition[2])), 
                        delays = data.frame(death = delay_2death), report = "o"),
    
    
    cm_multinom_process("E",
                        data.frame(to_severe = P.severe,
                                   to_critical = P.critical),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critival = delay_2severe)),
    
    cm_multinom_process("Ev_l",
                        data.frame(to_severe = P.severe*(1-VE$ve_severe_condition[1]),
                                   to_critical = P.critical*(1-VE$ve_critical_condition[1])),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critival = delay_2severe)),
    
    cm_multinom_process("Ev_m",
                        data.frame(to_severe = P.severe*(1-VE$ve_severe_condition[2]),
                                   to_critical = P.critical*(1-VE$ve_critical_condition[2])),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critival = delay_2severe)),
    
    cm_multinom_process("to_severe", 
                        data.frame(severe = rep(1,16)),                  
                        delays = data.frame(severe = delay_2hosp),   report = "ip"),
    
    cm_multinom_process("to_critical", 
                        data.frame(critical = rep(1,16)),                  
                        delays = data.frame(critical = delay_2hosp_critical),   report = "ip")
  )
  return(tmp)
}

burden_processes_pfizer <- gen_burden_processes(VE = ve_pfizer)
burden_processes_az <- gen_burden_processes(VE = ve_az)
burden_processes_pfizer_low <- gen_burden_processes(VE = ve_pfizer_low)
burden_processes_az_low <- gen_burden_processes(VE = ve_az_low)

