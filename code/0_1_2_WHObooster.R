### ADDITIONAL BOOSTER DOSE STRATEGY FROM 2024 ###

# This function considers additional booster dose campaigns from 2024
# Different (adult) age groups can be eligible for either annual or every 6m additional booster doses

# There are 3 phases:
      # Phase 1: before October 2022 - booster doses from owid data (booster_allocation_plan)
      # Phase 2: Oct 2022 - Dec 2023 - assume no further vaccination
      # Phase 3: January 2024 onwards - additional booster campaign in adults

# NOTE: booster campaigns must lie within the calendar year (i.e. end by December)
# NOTE: primary course coverage must be modified for vaccination of adolescents and/or children
# NOTE: the additional 6m campaign cannot overlap with the main annual campaign


vaccinate_additional <- function(para = NULL,
                          vac_data = owid_vac,
                          
                          booster_plan = booster_allocation_plan, 
                          #this is for existing boosters
                          
                          start_age_annual = 55,
                          # minimum age to be eligible for annual booster
                          
                          start_age_6m = 75,
                          # minimum age to be eligible for 6m booster
                          
                          cov_2024 = 0.5,
                          # percentage of people vaccinated with primary course who get additional booster
                          
                          month_annual = c(5:6),
                          # months to vaccinate everyone eligible for 12m and 6m booster
                          
                          month_6m = c(11:12)
                          # months to additionally vaccinate those eligible for 6m booster
){
  
  require(lubridate)
  
  # Create time list covering all phases
  date_start <- ymd(para$date0)
  date_additional_booster <- "2024-01-01" # Start of additional booster doses (phase 3)
  date_end <- date_start + para$time1
  
  data.frame(date = seq(date_start, date_end, by = "day")) |> 
    mutate(t = 0:para$time1) -> tmp_times_full
  
  
  ### PHASE 1 BOOSTER ###
  
  # Phase 1 time array
  tmp_times_phase1 <- tmp_times_full |> 
    mutate(empirical = date %in% (vac_data$date)) |> 
    filter(empirical == T) |> 
    pull(t) 
  tmp_times_phase1 <- c(0, tmp_times_phase1, max(tmp_times_phase1)+1) 
  #### Ja - I have zeros each end as it is in Yang's code. Can you check if it is still needed?
  
  # Phase 1 list of daily vaccine doses
  tmp_values_phase1 <- c(list(rep(0,16)), booster_plan, list(rep(0,16)))
  testthat::expect_equal(length(tmp_times_phase1), length(tmp_values_phase1))
  
  
  ### PHASE 2 NO BOOSTER ###
  
  # Phase 2 time array
  tmp_times_phase2 <- tmp_times_full |> 
    filter(date < date_additional_booster) |> 
    mutate(empirical2 = !(date %in% (vac_data$date))) |> 
    filter(empirical2 == T) |>
    pull(t) 
  
  # Phase 2 list of daily vaccine doses (assumes no vaccination)
  tmp_values_phase2 <- c(rep(list(rep(0,16)),length(tmp_times_phase2)))
  testthat::expect_equal(length(tmp_times_phase2), length(tmp_values_phase2))
  
  
  ### PHASE 3 ADDITIONAL BOOSTER DOSES ###
  
  # Phase 3 time array
  tmp_times_phase3 <- tmp_times_full |>  
    filter(date >= date_additional_booster) |> 
    pull(t) 
  
  # Phase 3 list of daily vaccine doses
  
  eligible_annual <- c(rep(NA,(start_age_annual/5)), rep(1,(16-start_age_annual/5))) 
    # eligible age groups for annual (12m + 6m) => eligible is 1, ineligible is NA
  eligible_6m <- c(rep(NA,(start_age_6m/5)), rep(1,(16-start_age_6m/5)))
    # eligible age groups for additional booster (6m only) => eligible is 1, ineligible is NA
  
  pop <- para$pop[[1]]$size 
    # population in each age group
  cov_primary <- as.numeric(covered_final[which(covered_final$age_group2 == "adult"),"p_covered"]) 
    # primary series coverage in adults
  
  doses_total_annual <- cov_2024 * pop * cov_primary * eligible_annual
    # total vaccine doses per age group in annual vaccination round
  doses_total_6m <- cov_2024 * pop * cov_primary * eligible_6m
    # total vaccine doses per age group in additional 6m vaccination round
  
  additional_booster <- tmp_times_full |> 
    filter(date >= date_additional_booster) |>
    mutate(m = lubridate::month(date),
           y = lubridate::year(date),
           campaign_annual = (m %in% month_annual),
           campaign_6m = (m %in% month_6m)
           )
  
  duration_annual <- additional_booster |> 
    filter(y == 2025) |> 
    pull(campaign_annual) |> sum() # number of days in annual campaign
  
  duration_6m <- additional_booster |> 
    filter(y == 2025) |> 
    pull(campaign_annual) |> sum() # number of days in supplementary 6m campaign
  
  additional_booster <- additional_booster |> 
    mutate(doses_per_day = case_when(
      campaign_annual == T ~ list(eligible_annual * doses_total_annual / duration_annual),
      campaign_6m == T ~ list(eligible_6m * doses_total_6m / duration_6m),
      TRUE ~ list(rep(0,16)))) # calculate doses per day for each age group
  
  tmp_values_phase3 <- c(additional_booster$doses_per_day)
  testthat::expect_equal(length(tmp_times_phase3), length(tmp_values_phase3))
  
  
  ### CREATE SCHEDULE OBJECT ###
  
  para$schedule[["booster"]] <- list(
    parameter = "v_b",
    pops = numeric(),
    mode = "assign",
    values = c(tmp_values_phase1, tmp_values_phase2, tmp_values_phase3),
    times = c(tmp_times_phase1, tmp_times_phase2, tmp_times_phase3)
  )
  
  return(para)
  
}
