# vaccination vs no vaccination, 1 vs 11
n = 1
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "no vaccine") -> a

n = 11
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "50% uptake, with prioritisation, annual") -> b
  
bind_rows(a, b) %>% 
  mutate(setting = factor(setting, levels = c("no vaccine",
                             "50% uptake, with prioritisation, annual"))) %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~setting, ncol = 1) +
  theme_hitap +
  labs(x = "Date", y = "Daily doses administered by age group",
       title = "Comparing uptake") -> dim_1

# with and without prioritisation
n = 11
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "50% uptake, with prioritisation, annual") -> a

n = 31
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "50% uptake, without prioritisation, annual") -> b

bind_rows(a, b) %>% 
  mutate(setting = factor(setting, levels = c("50% uptake, with prioritisation, annual",
                                              "50% uptake, without prioritisation, annual"))) %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~setting, ncol = 1) +
  theme_hitap +
  labs(x = "Date", y = "Daily doses administered by age group",
       title = "Comparing prioritisation strategy") -> dim_2

# annual vs bi annual, 11 compare to 12
n = 9
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "40% uptake, with prioritisation, annual") -> a

n = 10
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "40% uptake, with prioritisation, bi-annual") -> b

n = 17
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, para$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name == "X16") %>% 
  mutate(setting = "80% uptake, with prioritisation, bi-annual") -> c

bind_rows(a, b, c) %>% 
  mutate(setting = factor(setting, levels = c("40% uptake, with prioritisation, annual",
                                              "40% uptake, with prioritisation, bi-annual",
                                              "80% uptake, with prioritisation, bi-annual"))) %>% 
  ggplot(., aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~setting, ncol = 1) +
  theme_hitap +
  labs(x = "Date", y = "Daily doses administered by age group",
       title = "Comparing frequency") -> dim_3

plot_grid(dim_1 + geom_vline(xintercept = seq(ymd("2023-01-01"), ymd("2030-01-01"), by = "year"), linetype = 2), 
          dim_2 + geom_vline(xintercept = seq(ymd("2023-01-01"), ymd("2030-01-01"), by = "year"), linetype = 2), 
          dim_3 + geom_vline(xintercept = seq(ymd("2023-01-01"), ymd("2030-01-01"), by = "year"), linetype = 2), ncol = 1) -> p

ggsave("figs/fig3_v0.png", height = 15, width = 9)

