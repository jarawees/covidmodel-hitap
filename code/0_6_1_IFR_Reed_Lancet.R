read_excel(paste0(data_path, "IFR_Reed.xlsx")) %>% 
  mutate(`IFR and 95% uncertainty interval` = gsub("·",".",`IFR and 95% uncertainty interval`)) %>% 
  separate(`IFR and 95% uncertainty interval`,
           into = c("point", "range"),
           sep = "%") %>% 
  mutate(point = as.numeric(point)/100) %>% 
  separate(range,
           into = c("LL", "UL"),
           sep = "–") %>% 
  mutate(LL = parse_number(LL)/100,
         UL = parse_number(UL)/100,
         Age = parse_number(Age)) -> tmp

bind_rows(tmp[1,],
          tmp) -> tmp

tmp$Age[1] <- 0

tmp

write_rds(tmp,
          paste0(data_path,
                 "ifr_reed.rds"))
