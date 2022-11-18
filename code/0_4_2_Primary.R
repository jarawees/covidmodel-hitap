# separate adult doses to prioritse older adults
popTH_cm |> filter(age_group > 12) |> pull(tot_age) |> sum() -> pop_OA
popTH_cm |> filter(age_group <= 12 & age_group > 4) |> pull(tot_age) |> sum() -> pop_A
covered_final |> filter(age_group2 == "adult") |> pull(p_covered) -> covered_tmp
pop_OA*covered_tmp -> target_OA 
pop_A*covered_tmp -> target_A 
sapply(1:nrow(owid_vac), function(x) sum(owid_vac$fully_vaccinated_adult_daily[1:x])) -> daily_all_adults
(which(daily_all_adults < target_OA)) -> days_OA
(max((which(daily_all_adults < target_OA))) + 2):nrow(owid_vac) -> days_A

owid_vac |>
  mutate(
    days = 1:n(),
    fully_vaccinated_OA_daily = case_when(
      days %in% days_OA ~ fully_vaccinated_adult_daily,
      days == max(days_OA) + 1 ~ target_OA - daily_all_adults[max(days_OA)],
      TRUE ~ 0
    ),
    fully_vaccinated_A_daily = case_when(
      days %in% days_A ~ fully_vaccinated_adult_daily,
      days == max(days_OA) + 1 ~ fully_vaccinated_adult_daily- fully_vaccinated_OA_daily,
      TRUE ~ 0
    )
  ) -> owid_vac

# get within group weights
popTH_cm |> 
  mutate(age_group = if_else(age_group > 16, 16, age_group)) |> 
  group_by(age_group) |> summarise(tot_cat = sum(tot_age)) |> 
  mutate(age_group3 = case_when(age_group == 2 ~ "children",
                                age_group == 1 ~ "toddlers",
                                age_group %in% 3:4 ~ "adolescents",
                                age_group %in% 5:12 ~ "adults",
                                age_group > 12 ~ "older adults")) |> 
  group_by(age_group3) |> 
  mutate(tot_cat3 = sum(tot_cat),
         age_group3 = factor(age_group3,
                             levels = c("toddlers",
                                        "children",
                                        "adolescents",
                                        "adults",
                                        "older adults")), 
         weights_within = tot_cat/tot_cat3) -> popTH_weights

to_allocate <- data.frame(date = owid_vac$date,
                          doses_infants = 0,
                          doses_children = owid_vac$fully_vaccinated_children_daily,
                          doses_adolescents = owid_vac$fully_vaccinated_adolescent_daily,
                          doses_adults = owid_vac$fully_vaccinated_A_daily,
                          doses_olderadults = owid_vac$fully_vaccinated_OA_daily)# ,
# code to check
# doses_all = owid_vac$people_fully_vaccinated_daily) |> 
# mutate(doses_check = doses_infants + doses_children + doses_adolescents + doses_adults + doses_olderadults)

lapply(1:nrow(to_allocate), function(m){
  map2(.x = popTH_weights |> group_by(age_group3) |> group_split() |> 
         map(pull, weights_within), 
       .y = to_allocate[m,2:6],
       function(x,y) x*y) |> unlist()
})  -> primary_allocation_plan

pre_tag_OA <- min(which(owid_vac$people_fully_vaccinated_daily > 0)) - 1
primary_allocation_plan[[pre_tag_OA]] <- c(rep(0,12), rep(0.1,4))

pre_tag_A <- min(which(owid_vac$fully_vaccinated_A_daily > 0)) - 1
primary_allocation_plan[[pre_tag_A]][5:12] <- 0.1
# 
pre_tag_adolescent <- min(which(owid_vac$fully_vaccinated_adolescent_daily > 0)) - 1
primary_allocation_plan[[pre_tag_adolescent]][3:4] <- 0.1

pre_tag_children<- min(which(owid_vac$fully_vaccinated_children_daily > 0)) - 1
primary_allocation_plan[[pre_tag_children]][2] <- 0.1

# primary_allocation_plan |>
#   map(data.frame) |> map(t) |> map(data.frame) |> bind_rows() |>
#   set_rownames(NULL) |>
#   setNames(paste0("y",1:16)) |>
#   mutate(date = owid_vac$date) |>
#   pivot_longer(starts_with("y")) |>
#   mutate(name = parse_number(name),
#          age_group3 = case_when(name == 2 ~ "children",
#                                 name == 1 ~ "toddlers",
#                                 name %in% 3:4 ~ "adolescents",
#                                 name %in% 5:11 ~ "adults",
#                                 name >= 12 ~ "older adults")) -> p_table
# 
# ggplot(data = p_table, aes(x = date, y = value, color = age_group3)) +
#   geom_line() +
#   facet_wrap(~name) -> p
# 
# ggsave("figs/diagnostics/age_specific_primary_allocation_plan.png",
#        plot = p,
#        width = 10, height = 8)
