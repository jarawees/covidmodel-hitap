# have fitting window to end in August 2021
# owid_vac %>% 
#   filter(people_fully_vaccinated_per_hundred < 10) %>% 
#   pull(date) %>% 
#   range()

# ymd("2021-04-01")-ymd("2021-02-01")
# ymd("2021-06-01")-ymd("2021-02-01")
# 
# 

fit_func <- function(input){

  suppressWarnings(
    gen_country_basics(
      country = "Thailand",
      R0_assumed = input[1],
      date_start = "2021-02-01",
      date_end = "2021-08-26",
      contact = contact_schedule,
      processes = gen_burden_processes(VE = efficacy_all),
      period_wn  = 3 * 365,
      # duration, waning of natural immunity
      period_wv_m2l = 1 * 365,
      # duration, waning from medium to low levels vaccine induced
      period_wv_h2m = 1 * 365,
      # duration, waning from medium to low levels vaccine induced
      prob_v_p_2l = 1,
      prob_v_p_2m = 0,
      prob_v_b_l2m = 0.5,
      deterministic = TRUE,
      scenario_primary = scenario3_primary,
      scenario_booster = scenario3_booster,
      seed = input[2]
    ) %>%
      update_u_y(
        para = .,
        date_switch = c("2021-01-15", "2021-07-05"),
        rc_u = c(1, 1.5),
        # relative changes in u
        rc_y = c(1, 1),
        # relative changes in y
        rc_ve = c(1, 0.9),
        # relative evasiveness
        efficacy_baseline = efficacy_all,
        efficacy_weights = efficacy_weights_one
      ) %>%
      emerge_VOC_burden(
        para = .,
        rc_severity = c(1, 1.5),
        # relative change in ihr and ifr
        efficacy_baseline = efficacy_all
      )  -> params_tmp
  )
  
  
  params_tmp %>%
      cm_simulate() %>%
      .[["dynamics"]] %>%
      filter(grepl("death", compartment)) %>%
      group_by(t, compartment) %>%
      summarise(value = sum(value), .groups = "drop") %>%
      mutate(date = ymd("2021-02-01") + t) %>%
      pivot_wider(names_from = compartment,
                  values_from = value) %>%
      mutate(deaths_sim = case_when(date <= "2021-01-15" ~ death_o,
                                    date > "2021-01-15" & date <= "2021-07-05" ~ death_voc1_o,
                                    date > "2021-07-05" & date <= "2022-01-01" ~ death_voc2_o),
             scaled = deaths_sim*as.numeric(input[3]),
             date = as.character(date)) -> predicted
    
  predicted %>% 
    mutate(date = ymd(date)) %>% 
      right_join(epi %>% 
                  mutate(txn_date = ymd(txn_date)),
                by = c("date" = "txn_date")) %>%
    mutate(scaled = if_else(is.na(scaled), 0.001, scaled),
           scaled = if_else((scaled) == 0, 0.001, scaled)) %>% 
    arrange(date) %>% 
    filter(date <= "2021-08-26") %>% 
    dplyr::select(date, scaled, new_death) %>% 
    rename(observed = new_death,
           predicted = scaled) %>% 
    mutate(observed = round(observed, 0),
           ll = dpois(observed, predicted, log = T)) %>%
      pull(ll) %>% sum -> a
  
  return(-a)
  
}

draw_fit <- function(input){
  suppressWarnings(
    gen_country_basics(
      country = "Thailand",
      R0_assumed = input[1],
      date_start = "2021-02-01",
      date_end = "2021-08-26",
      contact = contact_schedule,
      processes = gen_burden_processes(VE = efficacy_all),
      period_wn  = 3 * 365,
      # duration, waning of natural immunity
      period_wv_m2l = 1 * 365,
      # duration, waning from medium to low levels vaccine induced
      period_wv_h2m = 1 * 365,
      # duration, waning from medium to low levels vaccine induced
      prob_v_p_2l = 1,
      prob_v_p_2m = 0,
      prob_v_b_l2m = 0.5,
      deterministic = TRUE,
      seed = input[2]
    ) %>%
      update_u_y(
        para = .,
        date_switch = c("2021-01-15", "2021-07-05"),
        rc_u = c(1, 1.5),
        # relative changes in u
        rc_y = c(1, 1),
        # relative changes in y
        rc_ve = c(1, 0.9),
        # relative evasiveness
        efficacy_baseline = efficacy_all,
        efficacy_weights = efficacy_weights_one
      ) %>%
      emerge_VOC_burden(
        para = .,
        rc_severity = c(1, 1.5),
        # relative change in ihr and ifr
        efficacy_baseline = efficacy_all
      )  -> params_tmp
  )
  
  
  params_tmp %>%
    cm_simulate() %>%
    .[["dynamics"]] %>%
    filter(grepl("death", compartment)) %>%
    group_by(t, compartment) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(date = ymd("2021-02-01") + t) %>%
    pivot_wider(names_from = compartment,
                values_from = value) %>%
    mutate(deaths_sim = case_when(date <= "2021-01-15" ~ death_o,
                                  date > "2021-01-15" & date <= "2021-07-05" ~ death_voc1_o,
                                  date > "2021-07-05" & date <= "2022-01-01" ~ death_voc2_o),
           scaled = deaths_sim*as.numeric(input[3]),
           date = as.character(date)) -> predicted
  
  predicted %>% 
    mutate(date = ymd(date)) %>% 
    right_join(epi %>% 
                 mutate(txn_date = ymd(txn_date)),
               by = c("date" = "txn_date")) %>%
    mutate(scaled = if_else(is.na(scaled), 0.0001, scaled),
           deaths_sim = if_else(is.na(deaths_sim), 0.0001, deaths_sim)) %>% 
    arrange(date) %>% 
    filter(date <= "2021-08-26") %>% 
    dplyr::select(date, scaled, deaths_sim, new_death) %>% 
    rename(observed = new_death,
           predicted_unscaled = deaths_sim,
           predicted = scaled) %>% 
    ggplot(., aes(x = date)) +
    geom_point(aes(y = predicted), color = "green") +
    geom_point(aes(y = observed), color = "purple") -> p
  
  return(p)
}

controlDE <- list(reltol=1e-6, steptol=20, itermax = 400, trace = 10,
                  parallelType = 2)

DEoptim(fn = fit_func,
        lower = c(1, 1, 0.5),
        upper = c(5, 90, 0.9),
        control = controlDE) -> out

# write_rds(out, "data/out_20230328.rds")
# 
# 
# out <- read_rds("data/out.rds")
# out_20230327 <- read_rds("data/out_20230327.rds")
# out_20230328 <- read_rds("data/out_20230328.rds")
# input <- c(1.78887360, 0.05408726, 0.04184677 )
# draw_fit(out$optim$bestmem)
# draw_fit(out_20230327$optim$bestmem)
# draw_fit(out_20230328$optim$bestmem)
# draw_fit(input)
