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
                              prioritisation_followup = c(NA,rep(2,11),rep(1,4)),
                              campaign_month = c(10:12,1:2),
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
    filter(date >= "2023-09-30")
  
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
                    rep(0.703, 2),
                    rep(0.849, 12))
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
  
  follow_up_t <- c(0, follow_up_order$campaign_duration_bygroup, 366)
  
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
           date >= max(vac_data$date)) -> phase_list
  
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
  tmp_len <- length(proportions_allocation_rescaled)
  proportions_allocation_rescaled[[tmp_len+1]] <- (rep(0,16))
  
  tmp_times[["campaign"]] <- phase_list$t
  tmp_allocation[["campaign"]] <- list()
  
  for (j in 1:nrow(phase_list)) {
    for (i in 1:nrow(follow_up_order)) {
      if (phase_list$vaccination_phase[j] == i) {
        tmp_allocation[["campaign"]][[j]] <-
          proportions_allocation_rescaled[phase_list$vaccination_phase[j]]
      }
      if (phase_list$vaccination_phase[j] == nrow(follow_up_order) + 1) {
        tmp_allocation[["campaign"]][[j]] <-
          proportions_allocation_rescaled[phase_list$vaccination_phase[[phase_list$vaccination_phase[j]]]]
      }
    }
  }
  
  tmp_times_move <- c(unlist(tmp_times) |> array())
  tmp_values_move <- c(tmp_allocation$campaign |> purrr::flatten()) |> setNames(NULL)
  
  if(frequency == 2 & max(prioritisation_followup, na.rm = T) == 1){
    to_alter <- which(((1:length(tmp_values_move))%%4) == 3)
    for(a in to_alter){
      tmp_values_move[[a]] <- rep(0,16)
    }
  }
  
  if(frequency == 2 & max(prioritisation_followup, na.rm = T) == 2){
    to_alter <- c(4, 5, 10, 11, 16, 17, 22, 23)
    for(a in to_alter){
      tmp_values_move[[a]] <- rep(0,16)
    }
  }
  
  
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
