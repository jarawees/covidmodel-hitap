if(!require(pacman)) install.packages("pacman")
library(pacman)
p_load(tidyverse)

##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))

#### Simple example ####
cm_parameters_SEI3R("Thailand") %>% 
  cm_simulate -> test
