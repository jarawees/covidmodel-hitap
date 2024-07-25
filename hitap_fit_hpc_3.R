source("code/0_LoadAll.R")

source("code/1_1_fit_gen_country_basics.R")
source("code/1_2_renew_fit_func.R")
source("code/1_3_draw_fit.R")

controlDE <- list(reltol=1e-6, steptol=20, itermax = 400, trace = 10,
                  parallelType = 2)
tmp_country <- "Thailand"
fvt <- c(0.21, 0.23, 0.25, 0.27)

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
  
  write_rds(out, file =paste0("fit/fit_results/fit_", fvt[i],".rds"))
  
  p_save <- draw_fit(input = out$optim$bestmem,
                     country = "Thailand",
                     draw_end = T,
                     voc_features = voc_features_test %>% 
                       mutate(change_u = 1),
                     fit_vac_threshold = fvt[i],
                     dt_tmp = 0.3) +
    labs(title = fvt[i])
  
  ggsave(plot = p_save,
         filename = paste0("fit/fit_figures/", fvt[i], ".png"))
}

print("complete successfully!")
