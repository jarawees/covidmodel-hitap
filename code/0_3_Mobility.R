# download the data file from Google if it doesn't exist in your directory
# if(!file.exists(paste0("data/gm.csv"))){
#   download.file("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
#                 paste0("data/gm.csv"))
# }

# # update the data file from google if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/gm.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   download.file("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
#                 paste0("data/gm.csv"))
# }

# gm <- read_csv("data/gm.csv") %>%
#   filter(country_region == "Thailand")

# global <- read_csv("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")
# gm <- global |> filter(country_region == "Thailand")
# write_rds(gm, "data/MobilityReport.rds")

curves <- data.table(
  work_scaler = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133,
    0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271,
    0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41,
    0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552,
    0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701,
    0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856,
    0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029,
    1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188,
    1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361
  ),
  other_scaler = c(
    0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078,
    0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094,
    0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109,
    0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13,
    0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175,
    0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31,
    0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549,
    0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86,
    0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224,
    1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561,
    1.589, 1.617, 1.645, 1.673, 1.701
  ),
  perc = round(seq(0, 1.25, 0.01), 2)
)

curves |> 
  ggplot() +
  geom_point(aes(x = perc, y = other_scaler)) +
  geom_vline(xintercept = 1)

# load data
mobility_type_all <- c("retail", "grocery", "parks", 
                       "transit", "workplaces", "residential")
mobility_type_used <- c("retail", "grocery", "transit", "workplaces")
mobility <- read_rds("data/MobilityReport.rds") |> 
  dplyr::select(-c("sub_region_1", "sub_region_2", "metro_area", 
                   "iso_3166_2_code", "census_fips_code", 
                   "country_region_code", "place_id")) |> 
  group_by(country_region, date) |> 
  pivot_longer(ends_with("baseline"),
               names_to = "mobility_type",
               values_to = "mobility_level") |> 
  mutate(mobility_type = gsub("_percent_change_from_baseline",
                              "",
                              mobility_type)) |> 
  filter(grepl(pattern = paste0(mobility_type_used, collapse = "|"), 
              mobility_type))

# define date range for this analysis
mobility_date_max <- max(mobility$date)
stringency_date_max <- max(oxcgrt$Date)

# generate some covariates
mobility |> 
  mutate(month = factor(month(date), levels = 1:12),
         dow = factor(wday(date)),
         doy = yday(date),
         date_numeric = as.numeric(date),
         mobility_level = mobility_level/100) |> 
  separate(mobility_type, into = c("mt1","mt2","mt3")) |> 
  dplyr::select(-mt2, -mt3) |> rename(mobility_type = mt1) |> 
  left_join(oxcgrt,
            by = c("date" = "Date",
                   "country_region" = "CountryName")) |> 
  filter(date <= mobility_date_max,
         date <= stringency_date_max) -> fit_tab

# fit the gam model 
model_fit <- gam(formula = mobility_level ~  dow*mobility_type  + 
             s(ContainmentHealthIndex_Average) + month,
           data = fit_tab,
           na.action = na.omit)

summary(model_fit)

# generate prediction table
CJ(date = seq(range(fit_tab$date)[1],
              as.Date("2030-12-31"),
              by = 1),
   mobility_type = mobility_type_used) %>%
  .[, c("dow",
        "doy",
        "month",
        "date_numeric") := list(lubridate::wday(date) %>% factor,
                            lubridate::yday(date),
                            month(date),
                            if_else(date <= max(fit_tab$date),
                                    as.numeric(date),
                                    as.numeric(max(fit_tab$date))))]  %>%
  .[, month := factor(month, levels = 1:12)] %>% 
  .[, country_region := "Thailand"] %>%
  left_join(oxcgrt,
            by = c("date" = "Date",
                   "country_region" = "CountryName")) |> 
  split(by = "mobility_type") |> 
  map(arrange, date) |> 
  map(mutate, ContainmentHealthIndex_Average = imputeTS::na_locf(ContainmentHealthIndex_Average)) |> 
  bind_rows() -> pre_tab

# generate predictions
pre_tab[,"mobility_level_predicted"] <- predict.gam(model_fit, newdata = pre_tab)
pre_tab |> 
  left_join(fit_tab[,c("country_region",
                       "date",
                       "mobility_type",
                       "mobility_level")],
            by = c("country_region",
                   "date",
                   "mobility_type")) |> 
  mutate(source = if_else(is.na(mobility_level), "imputed", "observed"),
         mobility_level_complete = if_else(is.na(mobility_level),
                                           mobility_level_predicted,
                                           mobility_level)) |> 
  dplyr::select(country_region,
                "C1M_School closing",
                mobility_type,
                date,
                mobility_level_complete) |> 
  pivot_wider(names_from = mobility_type,
              values_from = mobility_level_complete) |> 
  rename(school = `C1M_School closing`) |> 
  arrange(date) |> 
  mutate(school = if_else(is.na(school), 
                          max(school, 
                              na.rm = T), school))  -> pre_tab_merged

pre_tab_merged |> 
  mutate_at(vars(mobility_type_used),
            .funs = function(x) x + 1) |> 
  mutate(workplaces = if_else(workplaces > 1.25, 1.25, workplaces),
         othx = 0.345*retail + 0.445*transit + 0.21*grocery,
         othx = if_else(othx > 1.25, 1.25, othx),
         workplaces = round(workplaces, 2),
         othx = round(othx, 2)) |> 
  left_join(curves[,c("perc","work_scaler")], 
            by = c("workplaces" = "perc"))  |> 
  left_join(curves[,c("perc", "other_scaler")], 
            by = c("othx" = "perc")) |> 
  dplyr::select(-c(grocery, retail, transit, workplaces, othx)) %>%
  rename(work = work_scaler,
         other = other_scaler) |> 
  mutate(school = school/max(school, na.rm = T),
         home = 1) %>%
  dplyr::select(country_region, date, home, work, school, other) -> gm_scaled

CJ(date = seq(as.Date("2019-12-01"), as.Date("2020-02-14"),1),
   country_region = "Thailand") %>%
  .[,c("home",
       "work",
       "school",
       "other",
       "date") :=
      list(1,1,1,1,(date))] -> schedule_before

# school term
gm_scaled |> 
  bind_rows(schedule_before) |> 
  arrange(date) |> 
  mutate(month = month(date),
         day = day(date),
         year = year(date)) %>%
  mutate(holiday = if_else(
    #winter holiday,
    (year >= 2022 & month == 10 & day >=  11) |
      (year >= 2022 & month == 10 & day <= 31) |
      # summer holiday
      (year > 2022 & month == 4) |
      (year > 2022 & month == 5 & day <= 15),
    T,
    F),
    school = if_else(holiday, 0, school)) %>%
  dplyr::select(-holiday, -month, -day, -year) -> contact_schedule

# contact_schedule |> 
#   pivot_longer(c("home","work","school","other")) |> 
#   ggplot(aes(x = date, y = value)) +
#   geom_point() +
#   facet_wrap(~name)

