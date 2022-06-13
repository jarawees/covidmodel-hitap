# if(!require(pacman)) install.packages("pacman")
# library(pacman)
# p_load(tidyverse)

# require(tidyverse)

##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))

#### Simple example ####
params <- cm_parameters_SEI3R("Thailand") 
res <- cm_simulate(params)
