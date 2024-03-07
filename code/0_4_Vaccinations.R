# download the data file from owid if it doesn't exist in your directory
# if(!file.exists(paste0("data/vaccinations.csv"))){
#   download.file("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv",
#                 paste0("data/vaccinations.csv"))
# }
# 
# update the data file from owid if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/vaccinations.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   download.file("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv",
#                 paste0("data/vaccinations.csv"))
# }
# ve_i_o: observed vaccine efficacy against infection
# ve_d_o: observed vaccine efficacy against disease
# ve_d: VE against disease conditioned on changes occurring RE infection
# ve_severe, ve_critical, ve_mort: observed VE against different outcomes
# ve_severe_condition, ve_critical_condition, ve_mort_condition: VE against different outcomes condition on infection

data.table(v_i_o = c(0.33, 0.6, 0.75),
           vr_i_o = c(0.42, 0.84, 0.99), # update vr_i_o using meta-analysis from Dec2023
           # r_i_o = c(0.7), # original r_i_o
           # r_i_o = c(0.85), # from https://www.thelancet.com/article/S0140-6736(22)02465-5/fulltext
           r_i_o = c(0.76), # update r_i_o using meta-analysis from Dec2023
           # v_d_o = c(0.44, 0.69, 0.85),
           v_d_o = c(0.33, 0.72, 0.94), # update v_d_o using meta-analysis from Dec2023
           v_severe_o = c(0.59, 0.95, 0.99), # update v_severe_o using meta-analysis from Dec2023
           v_critical_o = c(0.8, 0.95, 0.99),
           v_mort_o = c(0.85, 0.95, 0.99),
           protection_level_label = c("l", "m", "h")) %>% 
  # the following lines do not explicit reflect existing changes in infection
  # which has been explicitly modelled as changes in u
  # This equation is explained in Liu et al. 
  # https://www.medrxiv.org/content/10.1101/2022.05.09.22274846v1.supplementary-material
  # Supplemental material, p37, version 1
  mutate(v_d_condition = 1 - (1-v_d_o)/((1-v_i_o)),
         v_severe_condition = 1 - (1-v_severe_o)/((1-v_i_o)),
         v_critical_condition = 1 - (1-v_critical_o)/((1-v_i_o)),
         v_mort_condition = 1 - (1-v_mort_o)/((1-v_i_o))) -> efficacy_all

fread(paste0(data_path, "vaccinations.csv")) %>%
# read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv") |> 
  filter(location == "Thailand") |> 
  arrange(date) |> 
  mutate(people_fully_vaccinated = as.numeric(people_fully_vaccinated),
         people_fully_vaccinated = if_else(date <= (fread(paste0(data_path, "vaccinations.csv")) %>%
                                             filter(location == "Thailand", !is.na(people_fully_vaccinated)) |> 
                                             arrange(date) |> pull(date) |> min()),
                                  0,
                                  people_fully_vaccinated),
         people_fully_vaccinated = imputeTS::na_interpolation(people_fully_vaccinated),
         
         people_vaccinated = as.numeric(people_vaccinated),
         people_vaccinated = if_else(date <= (fread(paste0(data_path, "vaccinations.csv")) %>%
                                                      filter(location == "Thailand", !is.na(people_vaccinated)) |> 
                                                      arrange(date) |> pull(date) |> min()),
                                           0,
                                     people_vaccinated),
         people_vaccinated = imputeTS::na_interpolation(people_vaccinated),
         # people_fully_vaccinated = imputeTS::na_interpolation(people_fully_vaccinated),
         
         total_vaccinations = as.numeric(total_vaccinations),
         total_vaccinations = imputeTS::na_interpolation(total_vaccinations),
         
         total_boosters = as.numeric(total_boosters),
         total_boosters = if_else(date <= (fread(paste0(data_path, "vaccinations.csv")) %>%
                                             filter(location == "Thailand", !is.na(total_boosters)) |> 
                                             arrange(date) |> pull(date) |> min()),
                                  0,
                                  total_boosters),
         total_boosters = imputeTS::na_interpolation(total_boosters),
         total_boosters_daily = c(0, diff(total_boosters)),
         daily_vaccinations = as.numeric(daily_vaccinations),
         # daily_vaccinations = imputeTS::na_interpolation(daily_vaccinations),
         daily_vaccinations_per_million = as.numeric(daily_vaccinations_per_million),
         # daily_vaccinations_per_million = imputeTS::na_interpolation(daily_vaccinations_per_million),
         date_numeric = as.numeric(date)) -> owid_vac

source("code/0_4_1_Staged_Vac.R")
source("code/0_4_2_Primary.R")
source("code/0_4_3_Boost_initial.R")

# owid_date_total <- data.frame(date = seq(range(owid_vac$date)[1], range(owid_vac$date)[2], by = "day"))
# owid_vac |> full_join(owid_date_total, by = "date")
# diff(owid_vac$date_numeric)

#### look for reasonable daily vaccination rates #####
# vaccinations/ million-day
# could test 190000 for starters

# lapply(c(seq(ymd("2021-09-30"),
#            ymd("2021-10-17"),
#            by = "day"),ymd("2022-08-30")), function(x){
#              owid_vac |>
#                filter(date < x)
#            }) |>
#   map(pull, daily_vaccinations) |>
#   map(summary) |> bind_rows()

# owid_vac |>
#   mutate(date = ymd(date),
#          people_vaccinated = as.numeric(people_vaccinated)) |>
#   ggplot() +
#   # geom_line(aes(x = date, y = people_vaccinated)) +
#   # geom_line(aes(x = date, y = people_fully_vaccinated), color = "red") +
#   # geom_line(aes(x = date, y = total_boosters ), color = "blue")+
#   geom_line(aes(x = date, y = daily_vaccinations), color = "green") + geom_hline(yintercept = 190000)

# owid_vac |> 
#   dplyr::select(date, total_vaccinations, total_boosters,
#                 daily_vaccinations, people_fully_vaccinated_per_hundred) |> 
#   mutate(total_vaccinations = as.numeric(total_vaccinations),
#          total_boosters = as.numeric(total_boosters),
#          total_primary = total_vaccinations - total_boosters,
#          daily_boosters = c(0, diff(total_boosters)),
#          daily_primary = c(0, diff(total_primary))
#          ) |> 
#   pivot_longer(cols = c("total_vaccinations", "total_boosters", "daily_boosters","daily_primary",
#                         "daily_vaccinations", "people_fully_vaccinated_per_hundred")) |> 
#   filter(name %in% c("daily_boosters", "daily_primary")) |> 
#   ggplot(aes(x = date, y = value)) +
#   geom_point() +
#   facet_wrap(~name, scales = "free", ncol = 1)

#### look for the transition time between primary dose and booster campaign ####
# lapply(seq(0,1,0.1),
#        function(x){
#          smooth.spline(x = owid_vac$date_numeric,
#                        y = owid_vac$people_fully_vaccinated,
#                        spar = x)
#        }) |>
#   map(predict,
#       deriv = 2) |>
#   map(function(x) {x["y"]}) |>
#   bind_cols() |>
#   setNames(paste0("spar_",seq(0,1,0.1))) |>
#   bind_cols(owid_vac[,"date"]) |>
#   pivot_longer(starts_with("spar")) -> splines_all

# we can see that spar = 0-0.3 doesn't really make sense because of too much
# permutation, the transition between initial vaccination and booster dose
# occurred in the first two weeks of October
# 2021-09-30 ~ 2021-10-17

# splines_all %>%
#   data.table |>
#   # split(by = "name") |>
#   # map(filter, value < 0 & date >= "2021-08-01") |>
#   # map(filter, date == min(date, na.rm = T)) |>
#   # bind_rows() |> 
#   ggplot(aes(x = date, y = value, group = name, color = name)) +
#   geom_line() +
#   facet_wrap(~name)

#### daily new vaccines
# owid_vac |> 
#   filter(date <= "2022-06-01" & date >= "2021-01-15") |> 
#   ggplot(aes(x = date, y = people_fully_vaccinated)) +
#   geom_line() +
#   geom_smooth(method = "lm")

# lm(people_fully_vaccinated ~ date, data = owid_vac) |> summary()
# owid_vac$people_fully_vaccinated |> tail(1)
