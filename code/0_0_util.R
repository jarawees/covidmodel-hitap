# foundamental wrapper of the functions
source("code/0_0_1_gen_country_basics.R")
source("code/0_0_2_update_u_y.R")
source("code/0_0_3_emerge_voc_burden.R")
source("code/0_0_4_vaccinate_primary.R")
source("code/0_0_5_vaccinate_additional.R")

# helper code file for this function is 0_4_1 and 0_4_2
# these code files help you prepare for these input objects: owid_vac, 
# primary_allocation_plan


# vaccinate_boost_initial <- function(para = NULL,
#                               vac_data = owid_vac,
#                               values = booster_allocation_plan
# ){
#   
#   require(lurbidate)
#   n_age_groups <- length(para$pop[[1]]$size)
#   date_start <- ymd(para$date0)
#   date_end <- date_start + para$time1
#   data.frame(date = seq(date_start, date_end, by = "day")) |> 
#     mutate(t = 0:para$time1,
#            empirical = date %in% (vac_data$date)) |> 
#     filter(empirical == T) |> 
#     pull(t) -> tmp_times
#   
#   c(0, tmp_times, max(tmp_times)+1) -> tmp_times
#   c(list(rep(0,16)), values, list(rep(0,16))) -> tmp_allocation
#   
#   testthat::expect_equal(length(tmp_times), length(tmp_allocation))
#   
#   para$schedule[["booster_initial"]] <- list(
#     parameter = "v_b",
#     pops = numeric(),
#     mode = "assign",
#     values = tmp_allocation,
#     times = tmp_times
#   )
#   return(para)
# }

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

# source("code/0_1_1_AnnualProgram.R")

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

check_vaccination_program <- function(type = "booster_initial", # or primary_course
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

parameterise_setting <- function(f = 1,
                                 prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
                                 boosting_level = 0.3){
  
  para <- gen_country_basics(country = "Thailand",
                             R0_assumed = out$optim$bestmem[1],
                             date_start = "2021-02-01",
                             date_end = "2027-12-31",
                             contact = contact_schedule,
                             processes = gen_burden_processes(VE = efficacy_all),
                             period_wn  = 3*365, # duration, waning of natural immunity
                             period_wv_m2l = 1*365, # duration, waning from medium to low levels vaccine induced 
                             period_wv_h2m = 1*365, # duration, waning from medium to low levels vaccine induced 
                             prob_v_p_2l = 1,
                             prob_v_p_2m = 0,
                             prob_v_b_l2m = 0.5,
                             deterministic = TRUE,
                             scenario_primary = scenario3_primary,
                             scenario_booster = scenario3_booster,
                             seed = out$optim$bestmem[2]) %>% 
    update_u_y(para = .,
               date_switch = c("2021-01-15", "2021-07-05", "2021-12-31"),
               rc_u = c(1, 1.5, 1.1), # relative changes in u
               rc_y = c(1, 1, 1), # relative changes in y
               rc_ve = c(1, 0.9, 0.7), # relative evasiveness 
               efficacy_baseline = efficacy_all
    ) %>%
    emerge_VOC_burden(para = .,
                      rc_severity = c(1, 1.5, 0.7), # relative change in ihr and ifr
                      efficacy_baseline = efficacy_all) %>%
    vaccinate_primary(para = .,
                      vac_data = owid_vac,
                      values = primary_allocation_plan) %>%
    vaccinate_booster(para = .,
                      vac_data = owid_vac,
                      booster_plan =  booster_allocation_plan,
                      # this is paused time
                      # program_interval = 30*6, #default set to 6 months
                      # should this be age-specific as well?
                      uptake_by_existing = boosting_level, 
                      # age-specific variables that defines the 
                      # prioritisation, the numbers are essentially just
                      # rankings; NA = not boosted
                      # this is future policy
                      prioritisation_followup = prioritisation_followup,
                      campaign_month = c(5:9),
                      frequency = f)
  return(para)
}

vaccinate_booster <- function(para = NULL,
                              vac_data = owid_vac,
                              booster_plan =  booster_allocation_plan,
                              # this is paused time
                              # program_interval = 30*6, #default set to 6 months
                              # should this be age-specific as well?
                              uptake_by_existing = 0.9, 
                              # age-specific variables that defines the 
                              # prioritisation, the numbers are essentially just
                              # rankings; NA = not boosted
                              # this is future policy
                              # prioritisation_followup = c(NA, rep(1,15)),
                              prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
                              # campaign_month = c(10:12,1:2),
                              campaign_month = c(5:9),
                              frequency = 1
                              # boosters_daily = 300000
){
  require(lubridate)
  # debug
  # para <- params
  # vac_data = owid_vac
  # uptake_by_existing = 0.3
  # prioritisation_followup = c(NA,rep(1,15))
  # campaign_month = c(10:12,1:2)
  # frequency = 1
  
  uptake_by_existing_tmp <- uptake_by_existing
  if(length(uptake_by_existing_tmp) == 1) uptake_by_existing_tmp <- rep(uptake_by_existing_tmp, 16)
  testthat::expect_length(uptake_by_existing_tmp, 16)
  
  time_range <- data.frame(date = seq(# ymd("2021-05-08"),
    ymd(para$date0),
    ymd(para$date0) + (para$time1) + 365,
    by = "day")) |> 
    rownames_to_column(var = "t") |> 
    mutate(t = as.numeric(t),
           m = lubridate::month(date),
           d = lubridate::day(date),
           campaign_days = if_else(m %in% campaign_month, T, F),
           doy = lubridate::yday(date),
           year = lubridate::year(date),
           vaccination_phase = NA) |> 
    group_by(year) %>% 
    filter(date >= "2023-10-01")
  
  time_range |> ungroup() |> filter(m == campaign_month[1], d == 1) -> season_start 
  time_range |> filter(doy >= 365) |> pull(doy) -> season_size
  
  time_range |>
    ungroup() |> 
    mutate(t_within = dplyr::lead(doy, (season_size[1] - season_start$doy[1] + 2))) |> 
    filter(date <= ymd(para$date0) + (para$time1)) -> time_range
  
  time_range |> mutate(year = lubridate::year(date)) |> filter(year == 2025) |> pull(campaign_days) |> sum() -> campaign_durations
  
  tmp_allocation <- tmp_times <- list()
  
  # initial boosting programmes
  # with owid data
  proportions_allocated_initial <- para$pop[[1]]$size/sum(para$pop[[1]]$size)
  # n_age_groups <- length(proportions_allocated_initial)
  # # we want the initial stage to not divide by stage and target all adults
  # proportions_allocated_initial_rescaled <- (prioritisation_initial*proportions_allocated_initial)/sum(prioritisation_initial*proportions_allocated_initial, na.rm = T)
  # proportions_allocated_initial_rescaled[is.na(proportions_allocated_initial_rescaled)] <- 0
  # 
  # vac_data |> 
  #   select(total_boosters_daily) %>%
  #   split(seq(nrow(.))) |> 
  #   map(unlist) |> 
  #   map(.f = function(x) x*proportions_allocated_initial_rescaled) |> 
  #   setNames(NULL) -> tmp_allocation[["initial"]]
  # 
  # vac_data |> 
  #   mutate(date = lubridate::date(date)) |> 
  #   left_join(time_range, by = "date") |> 
  #   dplyr::select(date, t) |> 
  #   pull(t) -> tmp_times[["initial"]]
  
  # follow-up campaigns
  # children coverage = 0.787; adolescent coverage = 0.812; adult coverage = 0.813
  data.frame(
    prioritisation_followup = prioritisation_followup,
    pop = para$pop[[1]]$size,
    cov_primary = c(NA, 
                    0.524, 
                    rep(0.811, 2),
                    rep(0.832, 12))
  ) |>
    mutate(
      cov_followup = uptake_by_existing_tmp * cov_primary,
      cov_followup_doses = cov_followup * pop,
      cov_followup_doses_all = sum(cov_followup_doses, na.rm = T),
      campaign_daily_dose =  round(cov_followup_doses_all / campaign_durations)
      # campaign_duration = round(cov_followup_doses_all / boosters_daily)
    ) |>
    group_by(prioritisation_followup) |>
    mutate(cov_followup_doses_bygroup = sum(cov_followup_doses, na.rm = T)) |> ungroup() |>
    mutate(campaign_duration_bygroup = round(cov_followup_doses_bygroup / campaign_daily_dose)) |>
    dplyr::select(prioritisation_followup, campaign_duration_bygroup, campaign_daily_dose) |> unique() |>
    filter(!is.na(prioritisation_followup)) |>
    arrange(prioritisation_followup)  -> follow_up_order
  
  follow_up_t <- c(0, cumsum(follow_up_order$campaign_duration_bygroup), 366)
  
  for(i in 1:(nrow(follow_up_order))){
    time_range[time_range$t_within <= follow_up_t[i+1] & time_range$t_within > follow_up_t[i], "vaccination_phase"] <- i
  }
  
  time_range |> 
    mutate(vaccination_phase = if_else(is.na(vaccination_phase), 
                                       as.numeric(nrow(follow_up_order)+1),
                                       as.numeric(vaccination_phase))) -> time_range
  
  time_range |> 
    group_by(vaccination_phase, year) |> 
    mutate(start = min(t_within)) |> 
    filter(start == t_within,
           date >= max(vac_data$date)) %>% 
    mutate(year_switch = (year%%frequency == 0))  -> phase_list
  
  # proportions OA rescaled
  testthat::expect_equal(length(unique(follow_up_order$campaign_daily_dose)),1)
  proportions_allocation_rescaled <- list()
  # when we are vaccinating people
  for(i in 1:nrow(follow_up_order)){
    tmp <-   as.numeric(prioritisation_followup == i) * proportions_allocated_initial /
      (sum(
        as.numeric(prioritisation_followup == i) * proportions_allocated_initial,
        na.rm = T
      )) * unique(follow_up_order$campaign_daily_dose)
    tmp[is.na(tmp)] <- 0
    proportions_allocation_rescaled[[i]] <- tmp
  }
  
  # taking care of the pause phase
  # tmp_len <- length(proportions_allocation_rescaled)
  # proportions_allocation_rescaled[[tmp_len+1]] <- (rep(0,16))
  
  tmp_times[["campaign"]] <- phase_list$t
  tmp_allocation[["campaign"]] <- list()
  
  for (j in 1:nrow(phase_list)) {
    for (i in 1:nrow(follow_up_order)) {
      if (phase_list$vaccination_phase[j] == i) {
        tmp_allocation[["campaign"]][[j]] <-
          proportions_allocation_rescaled[phase_list$vaccination_phase[j]]
      }
      if (phase_list$vaccination_phase[j] == max(phase_list$vaccination_phase) |
          phase_list$year_switch[j] == F) {
        tmp_allocation[["campaign"]][[j]] <- list(rep(0,16))
      }
    }
  }
  
  tmp_times_move <- c(unlist(tmp_times) |> array())
  tmp_values_move <- c(tmp_allocation$campaign |> purrr::flatten()) |> setNames(NULL)
  testthat::expect_equal(length(tmp_times_move),
                         length(tmp_values_move))
  
  # booster activities already observed
  n_age_groups <- length(para$pop[[1]]$size)
  date_start <- ymd(para$date0)
  date_end <- date_start + para$time1
  data.frame(date = seq(date_start, date_end, by = "day")) |> 
    mutate(t = 0:para$time1,
           empirical = date %in% (vac_data$date)) |> 
    filter(empirical == T) |> 
    pull(t) -> tmp_times_initial
  c(0, tmp_times_initial, max(tmp_times_initial)+1) -> tmp_times_initial
  c(list(rep(0,16)), booster_plan, list(rep(0,16))) -> tmp_values_initial
  testthat::expect_equal(length(tmp_times_initial), length(tmp_values_initial))
  
  para$schedule[["booster"]] <- list(
    parameter = "v_b",
    pops = numeric(),
    mode = "assign",
    values = c(tmp_values_initial, tmp_values_move),
    times = c(tmp_times_initial, tmp_times_move)
  )
  
  return(para)
}
