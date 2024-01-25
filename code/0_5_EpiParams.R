cf <- c(
  0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
  0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
  0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705
)

###### susceptibility ####
# (based on Davies et al, Nature paper)
sus <- c(
  0.3956736, 0.3956736, 0.3815349, 0.3815349, 0.7859512,
  0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
  0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189
)

#### load contact matrices ####
# Two versions of contact matrices done by Prem et al. (2017) and (2021)
# (2021) has been validated by DHS surveys, and therefore quality has been 
# improved especially in LMIC settings. Here, baseline CovidM uses 
# Prem et al. (2017) - here we are manually change the contact data source
# to Prem et al. (2021), so that this can be implemented directly using 
# COVIDM infrastructure

for(x in 2:5){
  load(paste0(data_path, list.files(data_path, pattern = "contact"))[x])
}

# countrycode::countrycode("Thailand", "country.name", "iso3c")
# "Brazil" is used here just to obtained rownames and column names 
cm_matrices[["Thailand"]]$home <- as.matrix(contact_home$THA) 
dimnames(cm_matrices[["Thailand"]]$home) <- 
  list(rownames(cm_matrices[["Brazil"]]$home),
       colnames(cm_matrices[["Brazil"]]$home))

cm_matrices[["Thailand"]]$work <- as.matrix(contact_work$THA) 
dimnames(cm_matrices[["Thailand"]]$work) <- 
  list(rownames(cm_matrices[["Brazil"]]$work),
       colnames(cm_matrices[["Brazil"]]$work))

cm_matrices[["Thailand"]]$school <- as.matrix(contact_school$THA) 
dimnames(cm_matrices[["Thailand"]]$school) <- 
  list(rownames(cm_matrices[["Brazil"]]$school),
       colnames(cm_matrices[["Brazil"]]$school))

cm_matrices[["Thailand"]]$other <- as.matrix(contact_others$THA) 
dimnames(cm_matrices[["Thailand"]]$other) <- 
  list(rownames(cm_matrices[["Brazil"]]$other),
       colnames(cm_matrices[["Brazil"]]$other))

# https://apps.who.int/gho/data/view.searo.61640?lang=en
mu <- read_csv(paste0(data_path, "THAI_MortalityRateByAge_2019.csv")) %>%
  setNames(.[1,]) %>% .[-1,] |> dplyr::select(-Male, -Female)|> 
  separate(`Age Group`, into = c("age_LL", "age_UL", sep = "-")) |> 
  mutate(age_LL = if_else(age_LL == "", 0, as.numeric(age_LL)),
         age_UL = if_else(age_UL == "years", 120, as.numeric(age_UL)))|> 
  dplyr::select(-`-`, -"Indicator") |> 
  rename(mu = `Both sexes`) |> 
  right_join(data.frame(age_LL = 0:120),
             by = "age_LL") |> 
  arrange(age_LL) |> 
  mutate(mu = as.numeric(mu),
         mu = imputeTS::na_locf(mu)) |> 
  dplyr::select(-age_UL) |> 
  rename(age = age_LL) |> 
  mutate(age_group_covidm = case_when(age >= 0  & age <= 4 ~ 1,
                                      age >= 5  & age <= 9 ~ 2,
                                      age >= 10 & age <= 14 ~ 3,
                                      age >= 15 & age <= 19 ~ 4,
                                      age >= 20 & age <= 24 ~ 5,
                                      age >= 25 & age <= 29 ~ 6,
                                      age >= 30 & age <= 34 ~ 7,
                                      age >= 35 & age <= 39 ~ 8,
                                      age >= 40 & age <= 44 ~ 9,
                                      age >= 45 & age <= 49 ~ 10,
                                      age >= 50 & age <= 54 ~ 11,
                                      age >= 55 & age <= 59 ~ 12,
                                      age >= 60 & age <= 64 ~ 13,
                                      age >= 65 & age <= 69 ~ 14,
                                      age >= 70 & age <= 74 ~ 15,
                                      age >= 75 ~ 16))

pop_by1 <- fread(paste0(data_path, "pop_str_2021.csv")) %>%
  gather(key = "sex", value = "pop", both) %>%
  mutate(pop = parse_number(pop)) |> 
  dplyr::select(-male,
                -female,
                -sex)

pop_by1 |> 
  left_join(mu,
            by = "age") |> 
  group_by(age_group_covidm) |> 
  summarise(mu_mean_year = sum(pop*mu)/sum(pop),
            mu_mean_day = mu_mean_year/365) -> mu_inuse

source("code/0_5_1_VoCs.R")

# owid_epi <- read_csv(url("https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-data.csv?raw=true")) %>% 
#   filter(location == "Thailand")
# write_rds(owid_epi, paste0(data_path, "owid_epi.rds"))

owid <- read_rds(paste0(data_path, "owid_epi.rds"))

# is this consistent with what we know about the population proejct outlook?
# test <- cm_parameters_SEI3R(dem_locations = "Thailand", 
#                             date_start = "2020-01-01", 
#                             date_end = "2030-12-31",
#                             A = rep(1/(365*5),16),
#                             # birth rate just from macrotrends.net
#                             B = c((9.532/1000)/365,rep(0,15)),
#                             D = mu_inuse$mu_mean_day,
#                             deterministic = T)
# 
# res <- cm_simulate(test)
# dyna <- res$dynamics |> 
#   filter(!grepl("foi|clinical|cases", compartment))
# 
# dyna |> 
#   group_by(t, group) |> 
#   summarise(value = sum(value)) |> 
#   ggplot(aes(x = t, y = value)) +
#   geom_point() +
#   facet_wrap(~group, scales = "free")
# 
# dyna |> 
#   group_by(t) |> 
#   summarise(value = sum(value)) |> 
#   ggplot(aes(x = t, y = value)) +
#   geom_point()
