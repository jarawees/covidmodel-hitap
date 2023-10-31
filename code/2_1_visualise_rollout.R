i = 34

data.frame(t_within = setting_list[[i]]$schedule$booster$times) %>% 
  mutate(date = t_within + ymd(setting_list[[i]]$date0)) -> date_grid1 

date_grid1 %>%   
  right_join(data.frame(date = seq(ymd("2021-02-01"),
                                   ymd("2027-12-31"),
                                   by = "day")),
             by = "date") %>% 
  arrange(date) -> date_grid2

setting_list[[i]]$schedule$booster$values %>% bind_cols() -> tmp

paste0("f = ", panel[i,1], 
       "; boosting level = ", panel[i,2],
       "; prioritisation = ", panel[i,3]) -> tmp_title

tmp %>% 
  t %>% 
  data.table %>% 
  bind_cols(date_grid1) %>% 
  right_join(date_grid2, by = c("date", "t_within")) %>% 
  arrange(date) %>% 
  mutate_at(vars(starts_with("V")),
            na_locf) %>% 
  pivot_longer(cols = starts_with("V"),
               names_to = "age_group") %>% 
  mutate(age_group = parse_number(age_group)) %>% 
  left_join(pop_TH,
            by = "age_group") %>% 
  mutate(daily_cov = value/pop_age,
         age_group_broad = case_when(age_group == 1 ~ "infants/ toddlers",
                                     age_group %in% 2 ~ "children",
                                     age_group %in% c(3,4) ~ "adolescents",
                                     age_group %in% c(5:12) ~ "adults",
                                     age_group %in% c(13:16) ~ "older adults"),
         age_group_broad = factor(age_group_broad,
                                  levels = c("infants/ toddlers",
                                             "children",
                                             "adolescents",
                                             "adults",
                                             "older adults"))) %>% 
  ggplot(., aes(x = date, y = daily_cov, group = age_group, color = age_group_broad)) +
  geom_line() +
  facet_wrap(~age_group_broad, ncol = 1)

