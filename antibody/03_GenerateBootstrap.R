## Prepare efficacy data ####
source('code/01_LoadData.R')
source('code/02_NAbToEffectivenessFunctions.R')

run_boot = F

if(run_boot){
  boot_agg <- GenerateBootstrap(df_agg %>% filter(day_measure <= 60), 10000, max(df_agg$SEM))
  saveRDS(boot_agg, file = "data/Bootstrap aggregated_60days.RDS")
} else {
  boot_agg <- readRDS("data/Bootstrap aggregated_60days.RDS")
}

## Plot NAb distribution for individual data ####
boot_agg_long <- pivot_longer(boot_agg, cols = 1:(ncol(boot_agg)-1), names_to = "group") %>% 
  filter(!grepl("Convalescent_No_WT", group)) %>% 
  separate(group, into = c("study_id", "regimen", "hybrid", "variant"), sep = "__")

boot_agg_long_nab <- boot_agg_long %>% slice(1:as.numeric(nrow(boot_agg_long)/3))
boot_agg_long_inf <- boot_agg_long %>% slice(as.numeric(nrow(boot_agg_long)/3+1):as.numeric(nrow(boot_agg_long)*2/3))
boot_agg_long_sev <- boot_agg_long %>% slice(as.numeric(nrow(boot_agg_long)*2/3+1):as.numeric(nrow(boot_agg_long)))

boot_agg_long <- boot_agg_long_nab %>% 
  bind_cols(Infection = boot_agg_long_inf %>% pull(value)) %>% 
  bind_cols(Severe = boot_agg_long_sev %>% pull(value)) %>% 
  rename(NAb = value) %>% 
  select(-label) %>% 
  left_join(author_name, by = "study_id")

rm(boot_agg_long_nab, boot_agg_long_inf, boot_agg_long_sev)


## VE table
# Immunity from vaccine
boot_agg_long %>% 
  filter(!grepl("Conv", regimen) & variant == "WT" & regimen %in% c("Moderna2", "Pfizer2") & hybrid == "No") %>% 
  mutate(NAb_level = as.numeric(cut_number(NAb, 3)),
         NAb_level = factor(NAb_level, levels = 1:3, labels = c("low", "medium", "high"))) %>% 
  group_by(NAb_level) %>% 
  summarise(n = n(),
            mean_inf = mean(Infection),
            mean_sev = mean(Severe))

# Hybrid immunity
boot_agg_long %>% 
  filter(!grepl("Conv", regimen) & variant == "WT" & regimen %in% c("Moderna2", "Pfizer2") & hybrid != "No") %>% 
  mutate(NAb_level = as.numeric(cut_number(NAb, 3)),
         NAb_level = factor(NAb_level, levels = 1:3, labels = c("low", "medium", "high"))) %>% 
  group_by(NAb_level) %>% 
  summarise(n = n(),
            mean_inf = mean(Infection),
            mean_sev = mean(Severe))

# Natural immunity
boot_agg_long %>% 
  filter(variant == "WT" & regimen %in% c("Convalescent")) %>% 
  mutate(NAb_level = as.numeric(cut_number(NAb, 3)),
         NAb_level = factor(NAb_level, levels = 1:3, labels = c("low", "medium", "high"))) %>% 
  group_by(NAb_level) %>% 
  summarise(n = n(),
            mean_inf = mean(Infection),
            mean_sev = mean(Severe))
