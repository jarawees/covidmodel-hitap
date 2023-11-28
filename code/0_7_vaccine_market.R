# load vaccine market share data

fn <- paste0(data_path, "Vaccination_market.xlsx")
sn <- excel_sheets(fn)
doc <- list()
for(i in c(1:3)) {doc[[i]] <- read_excel(fn, sheet = sn[i])}

vac_type_list <- c("sv", "az", "sp", "pf", "md")

fn <- paste0(data_path, "vaccine combo possibilities.csv")
combo <- read_csv(fn) %>% 
  .[,1:5] %>% 
  setNames(c("combo_index",
             "dose1",
             "dose2",
             "dose3",
             "applicability"))

# clean for dose 1
doc[[2]] %>% 
  dplyr::select(date, ends_with("first")) %>% 
  .[,c(1,3:7)] %>% 
  mutate(date = ymd(date)) -> dose1
dose1[1,2:6] <- 0
dose1 %>% 
  right_join(data.frame(date = seq(ymd(min(dose1$date)),
                                   ymd(max(dose1$date)),
                                   by = "day")),
             by = "date") %>% 
  arrange(date) %>% 
  mutate(sv_first_imputed = na_interpolation(sv_first),
         az_first_imputed = na_interpolation(az_first),
         sp_first_imputed = na_interpolation(sp_first),
         pf_first_imputed = na_interpolation(pf_first),
         md_first_imputed = na_interpolation(moderna_first),
         sv_first_daily = c(0, diff(sv_first_imputed)),
         az_first_daily = c(0, diff(az_first_imputed)),
         sp_first_daily = c(0, diff(sp_first_imputed)),
         pf_first_daily = c(0, diff(pf_first_imputed)),
         md_first_daily = c(0, diff(md_first_imputed))) -> dose1

# clean for dose2
doc[[2]] %>% 
  dplyr::select(date, ends_with("second")) %>% 
  .[,c(1,3:7)] %>% 
  mutate(date = ymd(date)) -> dose2
dose2[1,2:6] <- 0

dose2 %>% 
  right_join(data.frame(date = seq(ymd(min(dose2$date)),
                                   ymd(max(dose2$date)),
                                   by = "day")),
             by = "date") %>% 
  arrange(date) %>% 
  mutate(sv_second_imputed = na_interpolation(sv_second),
         az_second_imputed = na_interpolation(az_second),
         sp_second_imputed = na_interpolation(sp_second),
         pf_second_imputed = na_interpolation(pf_second),
         md_second_imputed = na_interpolation(moderna_second),
         sv_second_daily = c(0, diff(sv_second_imputed)),
         az_second_daily = c(0, diff(az_second_imputed)),
         sp_second_daily = c(0, diff(sp_second_imputed)),
         pf_second_daily = c(0, diff(pf_second_imputed)),
         md_second_daily = c(0, diff(md_second_imputed)))  -> dose2

# clean for dose3
doc[[2]] %>% 
  dplyr::select(date, ends_with("boost")) %>% 
  .[,c(1,3:7)] %>% 
  mutate(date = ymd(date)) -> dose3
dose3[1,2:6] <- 0

dose3 %>% 
  right_join(data.frame(date = seq(ymd(min(dose1$date)),
                                   ymd(max(dose1$date)),
                                   by = "day")),
             by = "date") %>% 
  arrange(date) %>% 
  mutate(sv_boost_imputed = na_interpolation(sv_boost),
         az_boost_imputed = na_interpolation(az_boost),
         sp_boost_imputed = na_interpolation(sp_boost),
         pf_boost_imputed = na_interpolation(pf_boost),
         md_boost_imputed = na_interpolation(moderna_boost),
         sv_boost_daily = c(0, diff(sv_boost_imputed)),
         az_boost_daily = c(0, diff(az_boost_imputed)),
         sp_boost_daily = c(0, diff(sp_boost_imputed)),
         pf_boost_daily = c(0, diff(pf_boost_imputed)),
         md_boost_daily = c(0, diff(md_boost_imputed)),
         sp_boost_daily = if_else(sp_boost_daily != 0, 0, sp_boost_daily),
         sv_boost_daily = if_else(sv_boost_daily != 0, 0, sv_boost_daily))  -> dose3


# dose3 %>% 
#   dplyr::select(date, ends_with("daily")) %>% 
#   arrange(date) %>% 
#   mutate_at(vars(ends_with("daily")), cumsum) %>% 
#   pivot_longer(ends_with("daily")) %>% 
#   ggplot(., aes(x = date, y = value, group = name, color = name)) +
#   geom_line()

# method_1 <- read_rds(paste0(data_path,"vaccine_market_identical_results.rds"))
# method_2 <- read_rds(paste0(data_path,"vaccine_market_proportional_results.rds"))

# method_1 %>% 
#   pivot_longer(starts_with(c("2_", "3_")),
#                names_to = "combo") %>%
#   mutate(method = "identical") %>% 
#   bind_rows(method_2 %>% 
#               pivot_longer(starts_with(c("1_", "2_", "3_")),
#                            names_to = "combo") %>% 
#               mutate(method = "proportional")) %>% 
#   ggplot(., aes(x = date, y = value, color = method)) +
# geom_line() +
#   facet_wrap(~combo)

# tail(method_1)
# tail(method_2)

# dose3 %>% 
#   dplyr::select(-ends_with("daily")) %>% 
#   pivot_longer(ends_with(c("boost")),
#                names_to = "vac_type_observed",
#                values_to = "observed") %>% 
#   pivot_longer(ends_with("imputed"),
#                names_to = "vac_type_imputed",
#                values_to = "imputed") %>% 
#   mutate(vac_type_imputed = substr(vac_type_imputed, 1, 2),
#          vac_type_observed = substr(vac_type_observed, 1, 2),
#          vac_type_observed = if_else(vac_type_observed == "mo", "md", vac_type_observed)) %>% 
#   filter(vac_type_observed == vac_type_imputed) %>% 
#   ggplot(., aes(x = date)) +
#   geom_point(aes(y = observed)) +
#   geom_line(aes(y = imputed)) +
#   facet_wrap(~vac_type_imputed)

# combo %>% 
#   filter(applicability == "Y") %>% 
#   dplyr::select(dose1, dose2, dose3) %>% 
#   unique() -> combo_Y
# 
# c(paste0("1_", unique(combo_Y$dose1)),
#   paste0("2_", unique(
#     paste(combo_Y$dose1, combo_Y$dose2, sep = "_")
#   )),
#   paste0("3_", unique(
#     paste(combo_Y$dose1, combo_Y$dose2, combo_Y$dose3, sep = "_")
#   ))) -> vaccine_compartments

# check variability between data source
# doc[[1]] |> 
#   dplyr::select(date, interval_second) |> 
#   dplyr::filter(!is.na(interval_second)) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   arrange(date) |> 
#   mutate(interval_second_cumsum = cumsum(interval_second)) |> 
#   full_join(owid_vac |> 
#               dplyr::select(date, people_fully_vaccinated) |> 
#               dplyr::filter(!is.na(people_fully_vaccinated)) |> 
#               mutate(date = lubridate::ymd(date)),
#             by = "date") |> 
#   mutate(diff = abs(people_fully_vaccinated - interval_second_cumsum),
#          marker_100K = if_else(diff >= 100000, date, ymd(NA))) %>%
#   # |> pull(diff) |> summary()
#   ggplot() +
#   # geom_point(aes(x = interval_second_cumsum, y = people_fully_vaccinated)) + geom_abline(intercept = 0, slope = 1) +
#   # scale_y_log10() + scale_x_log10()
#   geom_point(aes(x = date, y = interval_second_cumsum)) +
#   geom_vline(aes(xintercept = ymd("2022-09-30"))) +
#   geom_point(aes(x = date, y = people_fully_vaccinated), color = "red")
  # geom_vline(aes(xintercept = marker_100K, color = diff),
  #            size = 1.5)

# ggsave("figs/diagnostics/market_share_validation.png",
#        width = 10, height = 6)
  
# doc[[2]]|> 
#   dplyr::select(date, total_second) |> 
#   full_join(owid_vac |> 
#               dplyr::select(date, people_fully_vaccinated) |> 
#               dplyr::filter(!is.na(people_fully_vaccinated)) |> 
#               mutate(date = lubridate::ymd(date)),
#             by = "date") |> 
#   mutate(date = lubridate::ymd(date),
#          diff = abs(people_fully_vaccinated - total_second),
#          marker_100K = if_else(diff >= 100000, date, ymd(NA))) |> 
#   ggplot() +
#   geom_point(aes(x = date, y = total_second)) +
#   geom_point(aes(x = date, y = people_fully_vaccinated), color = "red") +
#     geom_vline(aes(xintercept = marker_100K, color = diff),
#                size = 1.5) 

# check for first boosters
# doc[[1]] |> 
#   dplyr::select(date, interval_boost) |> 
#   filter(interval_boost > 0) |> 
#   mutate(date = lubridate::ymd(date),
#          interval_boost_cumsum = cumsum(interval_boost )) |> 
#   full_join(owid_vac |> 
#               dplyr::select(date, total_boosters) |> 
#               filter(total_boosters > 0) |> 
#               mutate(date = ymd(date)),
#             by = "date") |> 
#   mutate(diff = total_boosters - interval_boost_cumsum ,
#          diff = abs(diff),
#          marker_100K = if_else(diff > 100000, date, ymd(NA))) |> 
#   # arrange(date) |> head(200)
#   ggplot() +
#   geom_point(aes(x = date, y = total_boosters), color = "red")+
#   geom_point(aes(x = date, y = interval_boost_cumsum)) + 
#   geom_vline(aes(xintercept = marker_100K, color = diff),
#              size = 1.5) 

# doc[[1]] |> 
#   dplyr::select(date, interval_boost) |> 
#   filter(interval_boost > 0) |> 
#   mutate(date = lubridate::ymd(date),
#          interval_boost_cumsum = cumsum(interval_boost )) |> 
#   full_join(owid_vac |> 
#               dplyr::select(date, total_boosters) |> 
#               filter(total_boosters > 0) |> 
#               mutate(date = ymd(date)),
#             by = "date") |> 
#   mutate(diff = total_boosters - interval_boost_cumsum ,
#          diff = abs(diff),
#          marker_100K = if_else(diff > 100000, date, ymd(NA)))

# check consistency among sums
# what are these discrepancies?
# doc[[1]] |> 
#   dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose, interval_boost) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   drop_na() %>%
#   mutate(rs = rowSums(.[,2:6]),
#          diff1 = interval_dose - interval_first - interval_second - interval_boost,
#          diff2 = rs - interval_second - interval_first) |> View()
#   
# doc[[1]] |> 
#   dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   drop_na() %>%
#   mutate(rs = rowSums(.[,2:6])) %>%
#   mutate_at(vars(starts_with("interval") & ends_with("primary")), function(x) x/.$rs) |> 
#   dplyr::select(date, starts_with("interval") & ends_with("primary")) |> 
#   drop_na() |> 
#   # mutate_at(vars(starts_with("interval")), function(x) c(0,diff(x))) |> 
#   pivot_longer(starts_with("interval")) |> 
#   filter(date >= ymd("2021-04-01")) |> 
#   separate(name, into = c("seg1","seg2", "seg3")) -> p_table
  
# ggplot(data = p_table, 
#          aes(x = date, y = value, color = seg2, fill = seg2)) +
#   geom_bar(position = "stack", stat = "identity") 
# 
# doc[[1]] |> 
#   dplyr::select(date, starts_with("interval") & ends_with("boost"), interval_boost) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   drop_na() %>%
#   mutate(rs = rowSums(.[,3:7])) %>%
#   mutate_at(vars(starts_with("interval") & ends_with("boost")), function(x) x/.$rs) |> 
#   dplyr::select(date, starts_with("interval") & ends_with("boost")) |> 
#   drop_na() |> 
#   pivot_longer(starts_with("interval")) |> 
#   separate(name, into = c("seg1","seg2", "seg3")) |> 
#   filter(seg2 != "boost") |> 
#   ggplot(aes(x = date, y = value, color = seg2, fill = seg2)) +
#   geom_bar(stat = "identity", position = "stack"

# aggregate data
# boost, ignore sp and sv for boosting
# doc[[1]] |> 
#   dplyr::select(date, starts_with("interval") & ends_with("boost"), interval_boost) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   drop_na() %>%
#   mutate(rs = rowSums(.[,c(3:6)])) %>%
#   mutate_at(vars(starts_with("interval") & ends_with("boost")), function(x) x/.$rs) |> 
#   dplyr::select(date, starts_with("interval") & ends_with("boost")) |> 
#   drop_na() |> 
#   pivot_longer(starts_with("interval")) |> 
#   separate(name, into = c("seg1","seg2", "seg3")) |> 
#   dplyr::filter(!seg2 %in% c("boost", "sp", "sv")) |> 
#   dplyr::select(-seg1, -seg3) |> 
#   rename(vac_type = seg2) |> 
#   mutate(vac_type2 = "booster") -> seg_boost

# doc[[1]] |> 
#   dplyr::select(date, starts_with("interval") & ends_with("primary"), interval_first, interval_second, interval_dose) |> 
#   mutate(date = lubridate::ymd(date)) |> 
#   drop_na() %>%
#   mutate(rs = rowSums(.[,2:6])) %>% 
#   mutate_at(vars(starts_with("interval") & ends_with("primary")), function(x) x/.$rs) |> 
#   dplyr::select(date, starts_with("interval") & ends_with("primary")) |> 
#   drop_na() |> 
#   #vmutate_at(vars(starts_with("interval")), function(x) c(0,diff(x))) |> 
#   pivot_longer(starts_with("interval")) |> 
#   separate(name, into = c("seg1","seg2", "seg3")) |> 
#   dplyr::select(-seg1, -seg3) |> rename(vac_type = seg2) |> 
#   mutate(vac_type2 = "primary") -> seg_primary

dose1 %>% 
  left_join(dose2, by = "date") %>% 
  left_join(dose3, by = "date") %>% 
  dplyr::select(date, ends_with("daily")) %>% 
  pivot_longer(ends_with("daily")) %>% 
  separate(name, into = c("vac_type",
                          "dose_index",
                          "rm")) %>% 
  dplyr::select(-rm) -> vaccine_daily

# expand.grid(dose1 = vac_type_list,
#             vac_type = vac_type_list) %>% 
#   dplyr::filter(dose1 == vac_type) %>% 
#   mutate(Vl = c(0.821, 0.129, 0.821, 0, 0),
#          Vm = c(0.171, 0.868, 0.171, 0.88, 0.88),
#          Vh = c(0.008, 0.003, 0.008, 0.12, 0.12)) -> vp_levels

## Update vp_levels using our recent meta-analysis
expand.grid(dose1 = vac_type_list,
            vac_type = vac_type_list) %>% 
  dplyr::filter(dose1 == vac_type) %>% 
  mutate(Vl = c(0.570, 0.2225, 0.599, 0.130, 0.022),
         Vm = c(0.364, 0.523, 0.338, 0.455, 0.255),
         Vh = c(0.066, 0.2545, 0.063, 0.415, 0.723)) -> vp_levels

expand_grid(booster_type = vac_type_list) %>% 
  filter(booster_type %in% c("pf", "md", "az")) %>% 
  mutate(Vl2m = c(0.5,0,0),
         Vl2h = c(0.5,1,1)) -> bp_levels

source("code/0_7_3_consistent_primary_memoryless.R")
