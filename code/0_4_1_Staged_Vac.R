# impute people fully vaccinated and then convert to daily
model <- gam(people_fully_vaccinated ~ s(total_vaccinations), data = owid_vac)
summary(model)

fitted_values <- (predict.gam(model, newdata = data.frame(total_vaccinations = owid_vac$total_vaccinations))) %>%
  replace(., .<0,0)
#
owid_vac[,"people_fully_vaccinated_imputed"] <- fitted_values
owid_vac |>
  mutate(people_fully_vaccinated_daily = c(0,diff(people_fully_vaccinated_imputed))) -> owid_vac

# plot(owid_vac$people_fully_vaccinated_daily)
# lines(fitted_values)

# add in vaccination phases
phase1 <- seq(owid_vac |> filter(total_vaccinations >= 0) |> pull(date) |> min(),
            ymd("2021-09-30"), 
            by = "day")
phase2 <- seq(ymd("2021-10-01"),
              ymd("2022-01-31"),
              by = "day")
# phase3 <- seq(ymd("2022-02-01"),
#               ymd("2022-08-30"),
#               by = "day")       
phase3 <- seq(ymd("2022-02-01"),
              owid_vac |> filter(total_vaccinations > 0) |> pull(date) |> max(),
              by = "day")            

owid_vac %<>% 
  mutate(vaccination_phase = case_when(
    date %in% phase1 ~ 1,
    date %in% phase2 ~ 2,
    date %in% phase3 ~ 3
  ))
# 
# target_distribution <- list()
# target_distribution[["OA"]] <- rnorm(10000, 0.8, 0.02) %>% .[.<1]
# target_distribution[["A"]] <- rnorm(10000, 0.8, 0.02) %>% .[.<1]
# target_distribution[["UA"]] <- rnorm(10000, 0.5, 0.02) %>% .[.<1]
# 
# fit_beta <- list()
# fit_beta[["OA"]] <- EnvStats::ebeta(target_distribution[["OA"]])
# fit_beta[["A"]] <- EnvStats::ebeta(target_distribution[["A"]])
# fit_beta[["UA"]] <- EnvStats::ebeta(target_distribution[["UA"]])
# 
# calLogLik <- function(input){
#   # proportion of doses allocated for the corresponding age groups
#   coef_phase <- data.frame(adult = c(     1, input[1],   input[2]),
#                            adolescent = c(0, 1-input[1], input[3]),
#                            children = c(  0, 0,          1-input[2]-input[3]),
#                            vaccination_phase = c(1:3))
# 
#   owid_vac |>
#     left_join(coef_phase, by = "vaccination_phase") |>
#     mutate(fully_vaccinated_adult_daily = adult*people_fully_vaccinated_daily,
#            fully_vaccinated_adolescent_daily = adolescent*people_fully_vaccinated_daily,
#            fully_vaccinated_children_daily = children*people_fully_vaccinated_daily) |>
#     # the code here checks if things proceed according to the plan.
#     # mutate(fully_vaccinated_check = fully_vaccinated_adult + fully_vaccinated_adolescent + fully_vaccinated_children,
#     #        diff = people_fully_vaccinated - fully_vaccinated_check) |>
#     # pull(diff) |> table()
#     dplyr::select(starts_with("fully_vaccinated_")) |>
#     summarise_all(sum) |> t() |> data.frame() |>
#     rownames_to_column() |>
#     separate(rowname, into = paste0("v",1:4)) |>
#     dplyr::select(-v1, -v2, -v4) |>
#     setNames(c("age_group2", "covered")) -> covered_by_age
# 
#   # sum(covered_by_age$covered)/(popTH_cm |> filter(age_group > 1) |> pull(tot_age) |> sum())
# 
#   popTH_cm |>
#     mutate(age_group2 = case_when(age_group %in% 5:16 ~ "adult",
#                                   age_group %in% 3:4 ~ "adolescent",
#                                   age_group %in% 2 ~ "children")) |>
#     filter(!is.na(age_group2)) |>
#     group_by(age_group2) |> summarise(tot_age = sum(tot_age)) |>
#     left_join(covered_by_age, by = "age_group2") |>
#     mutate(p_covered = covered/tot_age) -> res
# 
#   return(data.frame(ll = sum(dbeta(res |> filter(age_group2 == "children") |> pull(p_covered),
#                                    fit_beta$UA$parameters[1],
#                                    fit_beta$UA$parameters[2],
#                                    log = T),
#                              dbeta(res |> filter(age_group2 == "adolescent") |> pull(p_covered),
#                                    fit_beta$A$parameters[1],
#                                    fit_beta$A$parameters[2],
#                                    log = T),
#                              dbeta(res |> filter(age_group2 == "adult") |> pull(p_covered),
#                                    fit_beta$OA$parameters[1],
#                                    fit_beta$OA$parameters[2],
#                                    log = T)),
#               cov_child = res$p_covered[3],
#               cov_adolescent = res$p_covered[1],
#               cov_adult = res$p_covered[2]))
# }
# 
# CJ(input1 = seq(0,1,0.05),
#    input2 = seq(0,1,0.05),
#    input3 = seq(0,1,0.05)) %>%
#   filter(input1 <= 1,
#          input2 + input3 <= 1) %>%
#   split(seq(nrow(.))) |>
#   map(unlist)  -> grid
# 
# res_final <- list()
# pb <- progress_bar$new(total = length(grid))
# for(i in 1:length(grid)) {
#   res_final[[i]] <- calLogLik(grid[[i]])
#   pb$tick()
# }
# 
# res_final |>
#   bind_rows(.id = "set") |>
#   filter(cov_adolescent < 1,
#          cov_adult < 1,
#          cov_child < 1,
#          ll != -Inf) |>
#   arrange(desc(ll)) |>
#   tibble() -> logLikRank
# 
# coef_phase_final  <- data.frame(adult = c(     1, grid[[logLikRank$set[1]]][1],   grid[[logLikRank$set[1]]][2]),
#                          adolescent = c(0, 1-grid[[logLikRank$set[1]]][1], grid[[logLikRank$set[1]]][3]),
#                          children = c(  0, 0,          1-grid[[logLikRank$set[1]]][2]-grid[[logLikRank$set[1]]][3]),
#                          vaccination_phase = c(1:3)) |>
#   set_rownames(NULL)
#  
# write_rds(coef_phase_final, "data/phased_introduction2.rds")

phased_introduction <- read_rds(paste0("data/phased_introduction.rds"))

owid_vac |> 
  left_join(phased_introduction, by = "vaccination_phase") |> 
  mutate(fully_vaccinated_adult_daily = adult*people_fully_vaccinated_daily,
         fully_vaccinated_adolescent_daily = adolescent*people_fully_vaccinated_daily,
         fully_vaccinated_children_daily = children*people_fully_vaccinated_daily) |> 
  # the code here checks if things proceed according to the plan. 
  # mutate(fully_vaccinated_check = fully_vaccinated_adult + fully_vaccinated_adolescent + fully_vaccinated_children,
  #        diff = people_fully_vaccinated - fully_vaccinated_check) |> 
  # pull(diff) |> table()
  dplyr::select(starts_with("fully_vaccinated_")) |> 
  summarise_all(sum) |> t() |> data.frame() |> 
  rownames_to_column() |> 
  separate(rowname, into = paste0("v",1:4)) |> 
  dplyr::select(-v1, -v2, -v4) |> 
  setNames(c("age_group2", "covered")) -> covered_by_age_final

# sum(covered_by_age$covered)/(popTH_cm |> filter(age_group > 1) |> pull(tot_age) |> sum())

popTH_cm |> 
  mutate(age_group2 = case_when(age_group %in% 5:16 ~ "adult",
                                age_group %in% 3:4 ~ "adolescent",
                                age_group %in% 2 ~ "children")) |> 
  filter(!is.na(age_group2)) |> 
  group_by(age_group2) |> summarise(tot_age = sum(tot_age)) |> 
  left_join(covered_by_age_final, by = "age_group2") |> 
  mutate(p_covered = covered/tot_age) -> covered_final

owid_vac %>%
  left_join(phased_introduction, by = "vaccination_phase") |> 
  mutate(fully_vaccinated_adult_daily = adult*people_fully_vaccinated_daily,
         fully_vaccinated_adolescent_daily = adolescent*people_fully_vaccinated_daily,
         fully_vaccinated_children_daily = children*people_fully_vaccinated_daily) -> owid_vac

# owid_vac |>
#   ggplot() +
#   geom_line(aes(x = date, y = fully_vaccinated_adult_daily), color = "green") +
#   geom_line(aes(x = date, y = fully_vaccinated_adolescent_daily), color = "orange") +
#   geom_line(aes(x = date, y = fully_vaccinated_children_daily), color = "blue")
