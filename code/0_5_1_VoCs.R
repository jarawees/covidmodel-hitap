voc <- read_csv("data/sars-cov-2-variants-dmsc.csv") %>% 
  .[,1:6] %>% 
  mutate(Date_Start = dmy(Date_Start), Date_End = dmy(Date_End)) %>% 
  .[complete.cases(.),] %>% 
  split(., 1:nrow(.)) %>% 
  map(mutate, range = list(seq(from = (Date_Start), to = (Date_End), by = "day"))) %>% 
  map(select, -Date_Start, -Date_End) %>% 
  map(unnest, cols = "range") %>% 
  bind_rows()