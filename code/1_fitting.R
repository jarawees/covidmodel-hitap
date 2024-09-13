source("code/0_LoadAll.R")
source("code/1_1_fit_gen_country_basics.R")
source("code/1_2_renew_fit_func.R")
source("code/1_3_draw_fit.R")

cl <- makeCluster(4)
clusterEvalQ(cl, source("code/0_LoadAll.R"))
clusterEvalQ(cl, source("code/1_1_fit_gen_country_basics.R"))
clusterEvalQ(cl, source("code/1_2_renew_fit_func.R"))
clusterEvalQ(cl, source("code/1_3_draw_fit.R"))

controlDE <- list(reltol=1e-4, 
                  steptol=20, 
                  itermax = 400, 
                  trace = 10,
                  parallelType = 2,
                  cluster = cl)

fvt <- seq(0.05, 0.3, 0.02)

for(i in 9:length(fvt)){
  model_to_fit <- renew_fit_func(fit_vac_threshold = fvt[i])
  DEoptim(fn = model_to_fit,
          lower = c(1.5, 1, 0.1),
          upper = c(3, 90, 1),
          control = controlDE) -> out
  
  write_rds(out, file = paste0("fit/fit_", fvt[i],"_3.rds"))
  
  p_save <- draw_fit(input = out$optim$bestmem,
                     draw_end = T,
                     voc_features = voc_features_test %>% 
                       mutate(change_u = 1),
                     fit_vac_threshold = fvt[i],
                     detection_threshold = 0.3) +
    labs(title = fvt[i])
  
  ggsave(plot = p_save,
         filename = paste0("fit/fit_figures/", fvt[i], "_3.png"))
}

# input,
# country = "Thailand",
# draw_end = T,
# voc_features = voc_features_test %>% mutate(change_u = 1),
# fit_vac_threshold = 0.1,
# detection_threshold = 0.3

# draw_fit(input = `fit_0.3`$optim$bestmem,
#          country = "Thailand",
#          draw_end = T,
#          voc_features = voc_features_test %>% 
#            mutate(change_u = 1),
#          fit_vac_threshold = 0.1,
#          dt_tmp = 0.3)
