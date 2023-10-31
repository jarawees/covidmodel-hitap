# impute people fully vaccinated and then convert to daily
model <- gam(people_fully_vaccinated ~ s(total_vaccinations), data = owid_vac)

fitted_values <- (predict.gam(model, newdata = data.frame(total_vaccinations = owid_vac$total_vaccinations))) %>%
  replace(., .<0,0)

owid_vac[,"people_fully_vaccinated_imputed"] <- fitted_values

owid_vac |>
  mutate(people_fully_vaccinated_daily = c(0,diff(people_fully_vaccinated_imputed))) -> owid_vac

# add in vaccination phases
# these three stages represent different initial targets for vaccination during 
# the campaign phase
phase1 <- seq(owid_vac |> filter(total_vaccinations >= 0) |> pull(date) |> min(),
            ymd("2021-09-30"), 
            by = "day")
phase2 <- seq(ymd("2021-10-01"),
              ymd("2022-01-31"),
              by = "day")
phase3 <- seq(ymd("2022-02-01"),
              owid_vac |> filter(total_vaccinations > 0) |> pull(date) |> max(),
              by = "day")            

owid_vac %<>% 
  mutate(vaccination_phase = case_when(
    date %in% phase1 ~ 1,
    date %in% phase2 ~ 2,
    date %in% phase3 ~ 3
  ))

# The code below is for fitting and does not need to be ran every single time
# source("code/0_4_1_1_fit_rollout.R")
# this set of phase introduction is based on the assumption that the vaccine 
# coverage among adults and adolescents is 80% and that among children is 20%
phased_introduction <- read_rds(paste0("data/phased_introduction2.rds"))

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
