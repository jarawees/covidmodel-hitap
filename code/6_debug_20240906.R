source("code/0_LoadAll.R")

# para <- cm_parameters_SEI3R(dem_locations = "Thailand")
# para$processes = burden_processes_all$THA
# res <- cm_simulate(para)

# country_tmp = "Thailand"
# country_code_tmp = "THA"
# date_start = "2020-01-01"
# date_end = "2025-12-31"
# processes_set = burden_processes_all
# prob_v_p_2l = 0.33
# prob_v_p_2m = 0.33
# prob_v_b_l2m = 0
# fitted_table_tmp = out_all[1,]
# period_wn  = 3*365 # duration, waning of natural immunity
# period_wv_m2l = 1*365 # duration, waning from medium to low levels vaccine induced
# period_wv_h2m = 1*365 # duration, waning from medium to low levels vaccine induced
# deterministic = TRUE


date_start = "2021-02-15"
fit_results_dir <- paste0("fit/")

out_all <- paste0(fit_results_dir, list.files(fit_results_dir, pattern = "fit_0")) %>%
  map(read_rds) %>%
  map(., "optim") %>%
  map(., "bestmem") %>%
  bind_rows() %>%
  mutate(fit_end_threshold = list.files(fit_results_dir, pattern = "fit_0") %>%
           gsub("fit_", "", .) %>%
           gsub(".rds", "", .) %>%
           as.numeric(),
         seed_raw = ymd(date_start) - 30 + par2 - ymd("2020-01-01"),
         seed_20200101 = as.numeric(seed_raw),
         country = "Thailand",
         continent = "Asia",
         country_code = "THA") %>%
  rename(R0_assumed_2 = par1,
         wn = par4) 

para <- gen_country_basics(date_start = "2020-01-01",
                           date_end = "2030-12-31",
                           processes_set = burden_processes_all,
                           prob_v_p_2l = 0.33,
                           prob_v_p_2m = 0.33,
                           prob_v_b_l2m = 0,
                           fitted_table_tmp = out_all[1,], # this needs to vary over the possible fitted results
                           period_wv_h2m = 1*365,
                           period_wv_m2l = 1*365,
                           deterministic = TRUE) %>%
  update_u_y(para = .,
             voc_features_inuse = voc_features_test %>% mutate(change_u = 1), # we have been told that the original version of voc_features_test overestimate changes in transmissibility 
             future_severe = F,
             future_severe_lvl = "mean",
             efficacy_baseline = efficacy_all 
  ) %>% 
  emerge_voc_burden(para = .,
                    voc_features_inuse = voc_features_test %>% mutate(change_u = 1),
                    detection_threshold = 0.3,
                    efficacy_baseline = efficacy_all,
                    split_E = T)

cm_simulate(para)
