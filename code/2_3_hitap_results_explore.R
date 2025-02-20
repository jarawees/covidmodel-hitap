# fn <- list.files("~/Dropbox/Github_Data/HITAP_CovidM/grid_results/")
# res <- paste0("~/Dropbox/Github_Data/HITAP_CovidM/grid_results/",fn) %>%
#   map(read_csv) %>%
#   setNames(fn) %>%
#   bind_rows(.id = "grid_table_index")
# 
# res %<>%
#   mutate(grid_table_index = gsub(".csv", "", grid_table_index)) %>%
#   left_join(grid_table, by = "grid_table_index") %>%
#   group_by(grid_table_index, year, compartment, fit_table_index, panel_final_index) %>%
#   summarise(incidence = sum(incidence))
# 
# write_rds(res,
#           "~/Dropbox/Github_Data/HITAP_CovidM/grid_results_combined.rds")

res <- read_rds("~/Dropbox/Github_Data/HITAP_CovidM/grid_results_combined.rds")

res %>% 
  dplyr::filter(compartment == "severe",
                grid_table_index == 121)

res %>% 
  left_join(out_all %>% 
              rowid_to_column(var = 'fit_table_index')) %>% 
  left_join(panel_final %>% 
              rowid_to_column(var = "panel_final_index")) %>% 
  dplyr::filter(fit_table_index == 10,
                compartment == "death") -> p_tab

p_tab %>% 
  ggplot(., aes(x = year, 
                y = incidence, 
                group = interaction(cov_2024, scenario),
                color = cov_2024)) +
  geom_line() +
  facet_wrap(~scenario)

#### 
res %>% 
  left_join(out_all %>% 
              rowid_to_column(var = 'fit_table_index')) %>% 
  left_join(panel_final %>% 
              rowid_to_column(var = "panel_final_index")) %>% 
  dplyr::filter(compartment == "severe") -> p_tab

p_tab %>% dplyr::filter(incidence == 0.4)

p_tab %>% 
  ggplot(., aes(x = year, 
                y = incidence, 
                group = interaction(cov_2024, fit_table_index),
                color = cov_2024)) +
  geom_line() +
  facet_wrap(~scenario)
