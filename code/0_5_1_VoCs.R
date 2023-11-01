voc <- read_csv("data/sars-cov-2-variants-dmsc.csv") %>% 
  .[,1:6] %>% 
  mutate(Date_Start = dmy(Date_Start), Date_End = dmy(Date_End)) %>% 
  .[complete.cases(.),] %>% 
  split(., 1:nrow(.)) %>% 
  map(mutate, range = list(seq(from = (Date_Start), to = (Date_End), by = "day"))) %>% 
  map(select, -Date_Start, -Date_End) %>% 
  map(unnest, cols = "range") %>% 
  bind_rows()

voc %>% 
  pivot_longer(starts_with("B", ignore.case = F)) %>% 
  group_by(range) %>% 
  mutate(tot = sum(value),
         value_p = value/tot) %>% 
  ggplot(., aes(x = range, y = value_p, color = name, fill = name)) +
  geom_bar(stat = "identity", position = "stack")

# Delta started in July 2021
voc %>% 
  pivot_longer(starts_with("B", ignore.case = F)) %>% 
  group_by(range) %>% 
  mutate(tot = sum(value),
         value_p = value/tot) %>% 
  filter(grepl("Delta",name),
         value_p > 0.5) 

# Omicron started in Jan 2022
voc %>% 
  pivot_longer(starts_with("B", ignore.case = F)) %>% 
  group_by(range) %>% 
  mutate(tot = sum(value),
         value_p = value/tot) %>% 
  filter(grepl("Omicron",name),
         value_p > 0.5) 
