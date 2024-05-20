vaccinate_primary <- function(para = NULL,
                              vac_data = owid_vac,
                              values = primary_allocation_plan
){

  n_age_groups <- length(para$pop[[1]]$size)
  date_start <- ymd(para$date0)
  date_end <- date_start + para$time1
  data.frame(date = seq(date_start, date_end, by = "day")) |> 
    mutate(t = 0:para$time1,
           empirical = date %in% (vac_data$date)) |> 
    filter(empirical == T) |> 
    pull(t) -> tmp_times
  
  c(0, tmp_times, max(tmp_times)+1) -> tmp_times
  c(list(rep(0,16)), values, list(rep(0,16))) -> tmp_allocation
  
  testthat::expect_equal(length(tmp_times), length(tmp_allocation))
  
  para$schedule[["primary_course"]] <- list(
    parameter = "v_p",
    pops = numeric(),
    mode = "assign",
    values = tmp_allocation,
    times = tmp_times
  )
  return(para)
}