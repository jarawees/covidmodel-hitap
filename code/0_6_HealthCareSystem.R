# These parameters govern the "processes", and not "compartments"
# In this model in general, we have five levels of public health endpoints
# infection, disease, severe, critical (ICU), mortality
# infection and disease have already been taken care of as they are by product 
# of model run
# severe, critical, and mortality are taken care of using the parameters listed 
# below

# COVID-19 Clinical Information Network (COCIN)
critical2 <- 0
picu_cocin_func <- function(age) {
  x <- c(-0.1309118, 0, 17.2398874, 65.7016492, 100)
  y <- c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
  p <- splinefun(x, y)(age)
  exp(p) / (1 + exp(p))
}
picu_cocin <- picu_cocin_func(0:85)

# Infection fatality rate (derived from COVID-19 Forecasting Team, Lancet, age-specific global estimates)
ifr_c19forecast <- c(0.0054, 0.0054, 0.0040, 0.00320, 0.0027, 0.0024, 
               0.0023, 0.0023, 0.0023, 0.0025,0.0028,
               0.0031, 0.0036, 0.00420, 0.0050, 0.0060, 
               0.0071, 0.0085, 0.0100, 0.0118, 0.0138,
               0.0162, 0.0188, 0.02190, 0.0254, 0.293, 
               0.0337, 0.0386, 0.0442, 0.0504, 0.0573,
               0.0650, 0.0735, 0.08290, 0.0932, 0.1046, 
               0.1171, 0.1307, 0.1455, 0.1616, 0.1789,
               0.1976, 0.2177, 0.2391, 0.2620, 0.2863, 
               0.3119, 0.3389, 0.3672, 0.3968, 0.4278,
               0.4606, 0.4958, 0.5342, 0.5766, 0.6242, 
               0.6785, 0.7413, 0.8149, 0.9022, 1.0035,
               1.1162, 1.2413, 1.3803, 1.5346, 1.7058, 
               1.8957, 2.1064, 2.3399, 2.5986, 2.8851,
               3.2022, 3.5527, 3.9402, 4.3679, 4.8397, 
               5.3597, 5.9320, 6.5612, 7.2520, 8.0093,
               8.8381, 9.7437, 10.7311, 11.8054, 12.9717)/100
# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje <- exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85))
# Amalgamate probabilities
probabilities <- data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_c19forecast, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities <- probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
# critical and severe states are mutually exclusive
P.critical <- probabilities[, ihr * picu]
P.severe <- probabilities[, ihr * (1 - picu)]
P.death <- probabilities[, ifr]
P.hosp <- P.critical + P.severe

# the delay function for each outputs
delay_2death <- cm_delay_gamma(26, 5, 60, 0.25)$p
delay_2severe <- cm_delay_gamma(8.5, 5, 60, 0.25)$p
delay_2hosp <- cm_delay_gamma(14.6, 5, 60, 0.25)$p
delay_2hosp_critical <- cm_delay_gamma(15.6, 5, 60, 0.25)$p

# data.frame(delay_2death = delay_2death,
#            delay_2severe = delay_2severe) |> 
#   mutate(t = cm_delay_gamma(26, 5, 60, 0.25)$t) |> 
#   pivot_longer(starts_with("delay")) |> 
#   ggplot(aes(x = t, y = value, group = name, color = name)) +
#   geom_line()

gen_burden_processes <- function(VE){
  tmp <- list(
    # progressing to deaths
    # source names can be found in processes_spec.h
    cm_multinom_process(src = "newE",       
                        outcomes = data.frame(death = P.death),                   
                        delays = data.frame(death = delay_2death), report = "o"),
    cm_multinom_process("newEv_l",      
                        data.frame(death = P.death*(1-VE$v_mort_condition[1])), 
                        delays = data.frame(death = delay_2death), report = "o"),
    cm_multinom_process("newEv_m",     
                        data.frame(death = P.death*(1-VE$v_mort_condition[2])), 
                        delays = data.frame(death = delay_2death), report = "o"),
    cm_multinom_process("newEv_h",     
                        data.frame(death = P.death*(1-VE$v_mort_condition[3])), 
                        delays = data.frame(death = delay_2death), report = "o"),
    
    
    # progressing to severe and critical outcomes 
    cm_multinom_process(src = "newE",
                        outcomes = data.frame(to_severe = P.severe,
                                   to_critical = P.critical),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critical = delay_2severe)),
    
    cm_multinom_process("newEv_l",
                        data.frame(to_severe = P.severe*(1-VE$v_severe_condition[1]),
                                   to_critical = P.critical*(1-VE$v_critical_condition[1])),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critical = delay_2severe)),
    
    cm_multinom_process("newEv_m",
                        data.frame(to_severe = P.severe*(1-VE$v_severe_condition[2]),
                                   to_critical = P.critical*(1-VE$v_critical_condition[2])),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critical = delay_2severe)),
    cm_multinom_process("newEv_h",
                        data.frame(to_severe = P.severe*(1-VE$v_severe_condition[3]),
                                   to_critical = P.critical*(1-VE$v_critical_condition[3])),
                        delays = data.frame(to_severe = delay_2severe,
                                            to_critical = delay_2severe)),
    
    # channelling different together
    cm_multinom_process("to_severe", 
                        data.frame(severe = rep(1,16)),                  
                        delays = data.frame(severe = delay_2hosp),   report = "ip"),
    
    cm_multinom_process("to_critical", 
                        data.frame(critical = rep(1,16)),                  
                        delays = data.frame(critical = delay_2hosp_critical),   report = "ip")
  )
  return(tmp)
}

# the output is a 10 variable list object that define your health care system 
# output
# each list element include: 
# source: where's the process originating 
# type: multinomial
# names: 
# report: i - incidence; p - prevalence; o - observed
# probability: probability of progression
# delay: delay functions

burden_processes_all <- gen_burden_processes(VE = efficacy_all)

# burden_processes_az <- gen_burden_processes(VE = ve_az)

