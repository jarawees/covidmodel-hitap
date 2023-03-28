n = 1
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "no vaccine",
         frequency = "annual") -> p_1

setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "no vaccine",
         frequency = "bi-annual") -> p_2

n = 11
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults only",
         frequency = "annual") -> p_3
n = 12
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults only",
         frequency = "bi-annual") -> p_4

n = 31
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults then adults",
         frequency = "annual") -> p_5

n = 32
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults then adults",
         frequency = "bi-annual") -> p_6


n = 51
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults and adults",
         frequency = "annual") -> p_7

n = 52
setting_list[[n]]$schedule$booster$values %>%
  map(data.frame) %>% map(t) %>% map(data.frame) %>% bind_rows() %>% 
  mutate(t = setting_list[[n]]$schedule[["booster"]]$t) |>
  full_join(data.frame(t = seq(setting_list[[n]]$time0, setting_list[[n]]$time1)) |>
              mutate(date = lubridate::ymd(setting_list[[n]]$date0) + as.numeric(t)),
            by = "t") |>
  arrange(date) |>
  mutate_at(vars(starts_with("X")),
            imputeTS::na_locf) |>
  pivot_longer(starts_with("X")) |>
  mutate(name = factor(name,
                       levels = paste0("X", 1:16))) %>% 
  filter(name %in% c("X16", "X5")) %>% 
  mutate(setting = "50% uptake, older adults and adults",
         frequency = "bi-annual") -> p_8


bind_rows(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8) %>%
  mutate(setting = factor(setting, 
                          levels = c("no vaccine",
                                     "50% uptake, older adults only",
                                     "50% uptake, older adults and adults",
                                     "50% uptake, older adults then adults"),
                          labels = c("No vaccine use",
                                     "50% uptake\nolder adults only",
                                     "50% uptake\nolder adults and adults",
                                     "50% uptake\nolder adults then adults")),
         name = factor(name,
                       levels = c("X5", "X16"),
                       labels = c("20-24", "75+"))) %>% 
  ggplot(., aes(x = date, y = value, group = name, color = name)) +
  geom_line() +
  facet_grid(setting~frequency) +
  theme_hitap +
  labs(x = "Date", y = "Daily doses administered by age group",
       title = "Comparison by target age group and prioritisation", color = "Age Group (examples)") +
  ggsci::scale_color_lancet() +
  geom_vline(xintercept = ymd("2022-09-30"), linetype = 2)

ggsave("figs/fig3_v1.png", height = 12, width = 12)

