# let’s go with containment and health index or the stringency index

# download the data file from OxCGRT if it doesn't exist in your directory
# if(!file.exists(paste0("data/OxCGRT_latest.csv"))){
#   download.file("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_nat_latest.csv",
#                 paste0("data/OxCGRT_latest.csv"))
# }

# # update the data file from OxCGRT if the time difference is greater than a week
# if(as.numeric(abs(as.Date(file.info(paste0("data/OxCGRT_latest.csv"))$mtime) -
#                   as.Date(Sys.time()))) > 7){
#   download.file("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_nat_latest.csv",
#                 paste0("data/OxCGRT_latest.csv"))
# }

# school related contact data is located in this code file

oxcgrt <- fread(paste0(data_path, "OxCGRT_latest.csv")) %>%
  dplyr::select(-c(RegionName,RegionCode)) %>%
  filter(CountryName == "Thailand") %>%
  mutate_at(vars(Date), ~lubridate::ymd(.)) |> 
  dplyr::select(CountryName, Date, `C1M_School closing`, 
                "ContainmentHealthIndex_Average",
                "StringencyIndex_Average",
                "GovernmentResponseIndex_Average"
                ) |> 
  pivot_longer(ends_with("Average", ignore.case = F)) |> 
  mutate(value = if_else(is.na(value) & Date <= ymd("2020-02-01"),
                         0,
                         value)) |> 
  filter(!is.na(value)) |> 
  pivot_wider(names_from = name,
              values_from = value)
  
