# roughly 63% of all vaccinated received booster, assume children are not receiving boosters
# sum(owid_vac$total_boosters_daily)/sum(covered_final$covered[1:2])

covered_final %>% 
  mutate(b_covered = 0.63,
         combo_covered = p_covered*b_covered) -> covered_booster

popTH_cm |> filter(age_group <5 & age_group >2) |> pull(tot_age) |> sum() -> pop_adolescent
pop_OA*(covered_booster %>% filter(age_group2 == "adult") %>% pull(combo_covered)) -> target_OA 
pop_A*(covered_booster %>% filter(age_group2 == "adult") %>% pull(combo_covered)) -> target_A 
pop_adolescent*(covered_booster %>% filter(age_group2 == "adolescent") %>% pull(combo_covered)) -> target_adolescent


to_allocate_booster <- data.frame(date = owid_vac$date,
                          doses = owid_vac$total_boosters_daily) %>% 
  mutate(doses_cumsum = cumsum(doses))
which(to_allocate_booster$doses_cumsum <= target_OA) -> days_OA_booster
which(to_allocate_booster$doses_cumsum <= (target_OA + target_A) & to_allocate_booster$doses_cumsum > target_OA) -> days_A_booster
which(to_allocate_booster$doses_cumsum <= (target_OA + target_A + target_adolescent) & to_allocate_booster$doses_cumsum > (target_OA + target_A)) -> days_adolescent_booster

to_allocate_booster |>
  mutate(
    days = 1:n(),
    booster_OA_daily = case_when(
      days %in% days_OA_booster ~ doses,
      TRUE ~ 0
    ),
    booster_A_daily = case_when(
      days %in% days_A_booster ~ doses,
      TRUE ~ 0
    ),
    booster_adolescent_daily = case_when(
      days %in% days_adolescent_booster ~ doses,
      TRUE ~ 0
    ),
    booster_child_daily = 0,
    booster_infant_daily = 0
  ) %>% 
  dplyr::select(date,
                booster_infant_daily,
                booster_child_daily,
                booster_adolescent_daily,
                booster_A_daily,
                booster_OA_daily) -> to_allocate_booster

lapply(1:nrow(to_allocate_booster), function(m){
  map2(.x = popTH_weights |> group_by(age_group3) |> group_split() |> 
         map(pull, weights_within), 
       .y = to_allocate_booster[m,2:6],
       function(x,y) x*y) |> unlist()
})  -> booster_allocation_plan

pre_tag_OA <- min(which(to_allocate_booster$booster_OA_daily > 0)) - 1
booster_allocation_plan[[pre_tag_OA]] <- c(rep(0,12), rep(0.1,4))

pre_tag_A <- min(which(to_allocate_booster$booster_A_daily > 0)) - 1
booster_allocation_plan[[pre_tag_A]][5:12] <- 0.1
# 
pre_tag_adolescent <-  min(which(to_allocate_booster$booster_adolescent_daily > 0))  - 1
booster_allocation_plan[[pre_tag_adolescent]][3:4] <- 0.1

# booster_allocation_plan |>
#   map(data.frame) |> map(t) |> map(data.frame) |> bind_rows() |>
#   set_rownames(NULL) |>
#   setNames(paste0("y",1:16)) |>
#   mutate(date = owid_vac$date) |>
#   pivot_longer(starts_with("y")) |>
#   mutate(name = parse_number(name),
#          age_group3 = case_when(name == 2 ~ "children",
#                                 name == 1 ~ "toddlers",
#                                 name %in% 3:4 ~ "adolescents",
#                                 name %in% 5:11 ~ "adults",s
#                                 name >= 12 ~ "older adults")) -> p_table
# #
# ggplot(data = p_table, aes(x = date, y = value, color = age_group3)) +
#   geom_line() +
#   facet_wrap(~name) -> p
# #
# ggsave("figs/diagnostics/age_specific_boost_allocation_plan.png",
#        plot = p,
#        width = 10, height = 8)
