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

# This processing code does not need to be ran every time
# source("code/0_6_1_IFR_Reed_Lancet.R")
ifr_reed_data <- read_rds(paste0(data_path,  "ifr_reed.rds")) %>% 
  rename(age = Age)

fread(paste0(data_path, "pop_str_2021.csv")) %>%
  gather(key = "sex", value = "pop", both) %>%
  mutate(pop = parse_number(pop)) %>% 
  dplyr::select(-male, -female, -sex) %>% 
  left_join(ifr_reed_data, by = "age") %>% 
  mutate(point = na_locf(point),
         UL = na_locf(UL),
         LL = na_locf(LL)) %>% 
  mutate(age = if_else(age > 85, 85, age),
         pop_weighted = pop*point) %>% 
  group_by(age) %>% 
  mutate(pop_tot = sum(pop),
         pop_weight = pop/pop_tot,
         point_weight = sum(pop_weight*point)) %>% 
  dplyr::select(age, point_weight) %>% 
  distinct %>% 
  ungroup %>% 
  pull(point_weight) -> ifr_reed

  # mutate(levin = ifr_levin) %>% 
  # ggplot(., aes(x = point_weight, y = ifr_levin)) +
  # geom_point() +
  # geom_abline(intercept = 0, slope = 1)
  # 
  # summarise(pop = sum(pop),
  #           point = sum(pop*point)/sum(pop)) %>% tail()


# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin <- 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100
# ifr_reed %>% 
#   head(86) %>% 
#   mutate(levin = ifr_levin) %>% 
#   ggplot(., aes(x = point, y = levin, color = Age)) +
#   geom_point() +
#   geom_abline(intercept = 0,
#               slope = 1) +
#   labs(x = "reed") +
#   scale_x_log10() +
#   scale_y_log10()

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje <- exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85))
# Amalgamate probabilities
probabilities <- data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_reed, picu = picu_cocin)
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

