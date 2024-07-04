source("code/1_1_fit_gen_country_basics.R")
source("code/1_2_renew_fit_func.R")
source("code/1_3_draw_fit.R")

controlDE <- list(reltol=1e-4, steptol=20, itermax = 400, trace = 10,
                  parallelType = 2)
tmp_country <- "Thailand"
fvt <- 0.3 # seq(0.05, 0.3, 0.02)

for(i in 1:length(fvt)){
  model_to_fit <- renew_fit_func(country = tmp_country,
                                 dt_tmp = 0.3,
                                 fit_vac_threshold = fvt[i],
                                 voc_features = voc_features_test %>% 
                                   mutate(change_u = 1))
  DEoptim(fn = model_to_fit,
          lower = c(1.5, 1, 0.05, 1),
          upper = c(4, 90, 1, 3),
          control = controlDE) -> out
  
  write_rds(out, file =paste0("fit/fit_", fvt[i],".rds"))
  
  p_save <- draw_fit(input = out$optim$bestmem,
                     country = "Thailand",
                     draw_end = T,
                     voc_features = voc_features_test %>% 
                       mutate(change_u = 1),
                     fit_vac_threshold = fvt[i],
                     dt_tmp = 0.3) +
    labs(title = fvt[i])
  
  ggsave(plot = p_save,
         filename = paste0("fit/fit_figures/", fvt[i], "_2.png"))
}

fit_gen_country_basics(country_tmp = "Thailand",
                       country_code_tmp = "THA",
                       date_start = params_tmp$fit_start,
                       date_end = params_tmp$fit_end,
                       R0_assumed = 2,
                       period_wn = 3*365,
                       period_wv_m2l = 1*365, 
                       processes_set = burden_processes_all,
                       period_wv_h2m = 1*365, 
                       prob_v_p_2l = 0.33,
                       prob_v_p_2m = 0.33,
                       prob_v_b_l2m = 0,
                       seed = 10,
                       deterministic = TRUE)

input_tmp <- c(1.2, 90, 0.1)
