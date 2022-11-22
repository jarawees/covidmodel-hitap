# load vaccine market share data

fn <- "data/Vaccination_market.xlsx"
sn <- excel_sheets(fn)
doc <- list()
for(i in c(1:3)) {doc[[i]] <- read_excel(fn, sheet = sn[i])}

# check variability between data source
doc[[1]] |> 
  dplyr::select(date, interval_second) |> 
  dplyr::filter(!is.na(interval_second)) |> 
  mutate(date = lubridate::ymd(date)) |> 
  arrange(date) |> 
  mutate(interval_second_cumsum = cumsum(interval_second)) |> 
  full_join(owid_vac |> 
              dplyr::select(date, people_fully_vaccinated) |> 
              dplyr::filter(!is.na(people_fully_vaccinated)) |> 
              mutate(date = lubridate::ymd(date)),
            by = "date") |> 
  mutate(diff = abs(people_fully_vaccinated - interval_second_cumsum),
         marker_100K = if_else(diff >= 100000, date, ymd(NA))) %>%
  # |> pull(diff) |> summary()
  ggplot() +
  geom_point(aes(x = interval_second_cumsum, y = people_fully_vaccinated)) + geom_abline(intercept = 0, slope = 1) +
  scale_y_log10() + scale_x_log10()
  # geom_point(aes(x = date, y = interval_second_cumsum)) +
  # geom_point(aes(x = date, y = people_fully_vaccinated), color = "red") +
  # geom_vline(aes(xintercept = marker_100K, color = diff),
  #            size = 1.5) 

# ggsave("figs/diagnostics/market_share_validation.png",
#        width = 10, height = 6)
  
doc[[2]]|> 
  dplyr::select(date, total_second) |> 
  full_join(owid_vac |> 
              dplyr::select(date, people_fully_vaccinated) |> 
              dplyr::filter(!is.na(people_fully_vaccinated)) |> 
              mutate(date = lubridate::ymd(date)),
            by = "date") |> 
  mutate(date = lubridate::ymd(date),
         diff = abs(people_fully_vaccinated - total_second),
         marker_100K = if_else(diff >= 100000, date, ymd(NA))) |> 
  ggplot() +
  geom_point(aes(x = date, y = total_second)) +
  geom_point(aes(x = date, y = people_fully_vaccinated), color = "red") +
    geom_vline(aes(xintercept = marker_100K, color = diff),
               size = 1.5) 

# check for first boosters
doc[[1]] |> 
  dplyr::select(date, interval_boost) |> 
  filter(interval_boost > 0) |> 
  mutate(date = lubridate::ymd(date),
         interval_boost_cumsum = cumsum(interval_boost )) |> 
  full_join(owid_vac |> 
              dplyr::select(date, total_boosters) |> 
              filter(total_boosters > 0) |> 
              mutate(date = ymd(date)),
            by = "date") |> 
  mutate(diff = total_boosters - interval_boost_cumsum ,
         diff = abs(diff),
         marker_100K = if_else(diff > 100000, date, ymd(NA))) |> 
  # arrange(date) |> head(200)
  ggplot() +
  geom_point(aes(x = date, y = total_boosters), color = "red")+
  geom_point(aes(x = date, y = interval_boost_cumsum)) + 
  geom_vline(aes(xintercept = marker_100K, color = diff),
             size = 1.5) 

doc[[1]] |> 
  dplyr::select(date, interval_boost) |> 
  filter(interval_boost > 0) |> 
  mutate(date = lubridate::ymd(date),
         interval_boost_cumsum = cumsum(interval_boost )) |> 
  full_join(owid_vac |> 
              dplyr::select(date, total_boosters) |> 
              filter(total_boosters > 0) |> 
              mutate(date = ymd(date)),
            by = "date") |> 
  mutate(diff = total_boosters - interval_boost_cumsum ,
         diff = abs(diff),
         marker_100K = if_else(diff > 100000, date, ymd(NA)))

# check consistency among sums
# what are these discrepancies?
doc[[1]] |> 
  dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose, interval_boost) |> 
  mutate(date = lubridate::ymd(date)) |> 
  drop_na() %>%
  mutate(rs = rowSums(.[,2:6]),
         diff1 = interval_dose - interval_first - interval_second - interval_boost,
         diff2 = rs - interval_second - interval_first) |> View()
  
doc[[1]] |> 
  dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose) |> 
  mutate(date = lubridate::ymd(date)) |> 
  drop_na() %>%
  mutate(rs = rowSums(.[,2:6])) %>%
  mutate_at(vars(starts_with("interval") & ends_with("primary")), function(x) x/.$rs) |> 
  dplyr::select(date, starts_with("interval") & ends_with("primary")) |> 
  drop_na() |> 
  # mutate_at(vars(starts_with("interval")), function(x) c(0,diff(x))) |> 
  pivot_longer(starts_with("interval")) |> 
  filter(date >= ymd("2021-04-01")) |> 
  separate(name, into = c("seg1","seg2", "seg3")) -> p_table
  
ggplot(data = p_table, 
         aes(x = date, y = value, color = seg2, fill = seg2)) +
  geom_bar(position = "stack", stat = "identity") 

doc[[1]] |> 
  dplyr::select(date, starts_with("interval") & ends_with("boost"), interval_boost) |> 
  mutate(date = lubridate::ymd(date)) |> 
  drop_na() %>%
  mutate(rs = rowSums(.[,3:7])) %>%
  mutate_at(vars(starts_with("interval") & ends_with("boost")), function(x) x/.$rs) |> 
  dplyr::select(date, starts_with("interval") & ends_with("boost")) |> 
  drop_na() |> 
  pivot_longer(starts_with("interval")) |> 
  separate(name, into = c("seg1","seg2", "seg3")) |> 
  filter(seg2 != "boost") |> 
  ggplot(aes(x = date, y = value, color = seg2, fill = seg2)) +
  geom_bar(stat = "identity", position = "stack")


# aggregate data
# boost, ignore sp and sv for boosting
doc[[1]] |> 
  dplyr::select(date, starts_with("interval") & ends_with("boost"), interval_boost) |> 
  mutate(date = lubridate::ymd(date)) |> 
  drop_na() %>%
  mutate(rs = rowSums(.[,c(3:6)])) %>%
  mutate_at(vars(starts_with("interval") & ends_with("boost")), function(x) x/.$rs) |> 
  dplyr::select(date, starts_with("interval") & ends_with("boost")) |> 
  drop_na() |> 
  pivot_longer(starts_with("interval")) |> 
  separate(name, into = c("seg1","seg2", "seg3")) |> 
  dplyr::filter(!seg2 %in% c("boost", "sp", "sv")) |> 
  dplyr::select(-seg1, -seg3) |> 
  rename(vac_type = seg2) |> 
  mutate(vac_type2 = "booster") -> seg_boost

doc[[1]] |> 
  dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose) |> 
  mutate(date = lubridate::ymd(date)) |> 
  drop_na() %>%
  mutate(rs = rowSums(.[,2:6])) %>%
  mutate_at(vars(starts_with("interval") & ends_with("primary")), function(x) x/.$rs) |> 
  dplyr::select(date, starts_with("interval") & ends_with("primary")) |> 
  drop_na() |> 
  #vmutate_at(vars(starts_with("interval")), function(x) c(0,diff(x))) |> 
  pivot_longer(starts_with("interval")) |> 
  separate(name, into = c("seg1","seg2", "seg3")) |> 
  dplyr::select(-seg1, -seg3) |> rename(vac_type = seg2) |> 
  mutate(vac_type2 = "primary") -> seg_primary
