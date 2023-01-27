# COVID-19 neutralising antibody titre meta-analysis
# 26 January 2023
# Version 4

# Load libraries
if(!require(pacman)) install.packages("pacman")
library(pacman)
p_load(plyr, dplyr, stats, tidyverse, meta, purrr, ggridges)


### 1. Prepare data ####

## Variable description for CSV files
# NAb_meta1.csv files contain summary data of neutralising antibody titres and study details
# Type: either RCT (randomised controlled trial) or Obs (observational study)
# Group 1: vaccine regimen, SP – sinopharm, SV – sinovac, AZ – AstraZeneca, PZ – Pfizer, MD – Moderna, JJ – Johnson & Johnson
# Group 2: specific population groups, 18to65 – people aged between 18 and 65 years of age, HCW – healthcare worker
# LOD: limit of detection of the virus neutralisation assay
# Variant: 
# Day_measure: number of days after last vaccination or infection that NAb titre was measured

  # Read CSV files

    # Uncensored: all included studies
    csv_uncensored <- c("Anna_NAb_meta1.csv", "Dimple_NAb_meta1_v2.csv", "Dream_NAb_meta1_v2.csv", 
                        "Ja_NAb_meta1_PZ35.csv", "Pui_NAb_meta1.csv", "Siobhan_NAb_meta1.csv", 
                        "Saudamini_NAb_meta1.csv")
    csv_uncensored <- lapply(csv_uncensored, read.csv)
    csv_uncensored <- do.call("rbind", csv_uncensored)

    # Censored: subset of uncensored studies with individual data points
    csv_censored <- c("Dimple_NAb_meta2.csv", "Dream_NAb_meta2_v2.csv", "Ja_NAb_meta2_NODELTA.csv", 
                      "Pui_NAb_meta2.csv", "Siobhan_NAb_meta2_Pf5.csv")
    csv_censored <- lapply(csv_censored, read.csv) 
    csv_censored <- do.call("rbind", csv_censored)
    
    NAb_uncensored <- csv_uncensored %>%
      select(-DOI)
    NAb_censored_raw <- csv_censored %>%
      select(-DOI)
    
    NAb_uncensored$LOD <- as.numeric(NAb_uncensored$LOD)
    NAb_uncensored <- NAb_uncensored %>%
      mutate_at(vars(LOD), ~ ifelse(. > 2, log10(.), .)) %>%
      # Replace values below LOD with half of LOD
      mutate_at(vars(NAb_mean, NAb_SD, NAb_l95, NAb_u95, NAb_med, NAb_q1, NAb_q3, NAb_min, NAb_max), ~ ifelse(. == "<LOD", LOD/2, .)) %>%
      mutate_at(vars(Variant), ~ case_when(. == "Wuhan" ~ "WT",
                                           . == "D614" ~ "D614G",
                                           . == "kappa" ~ "Kappa",
                                           . == "Lambda " ~ "Lambda",
                                           . == "Epsilion" ~ "Epsilon",
                                           . == "theta" ~ "Theta",
                                           TRUE ~ .)) %>%
      filter(Variant != "SARS-COV-1") %>%
      mutate_at(vars(Group1), ~ str_replace_all(., c("Pf" = "PZ", "PF" = "PZ", "CV" = "SV")))
    
    #table(NAb_uncensored$Group1)

    NAb_censored_raw <- NAb_censored_raw %>%
      left_join(NAb_uncensored %>% select(Study_ID, LOD) %>% unique(), by = "Study_ID") %>%
      mutate_at(vars(LOD), ~ ifelse(. > 2, log10(.), .)) %>%
      mutate_at(vars(Titre_log), ~ ifelse(. == "<LOD", LOD/2, .)) %>%
      mutate_at(vars(Variant), ~ case_when(. == "XD" ~ "Deltacron",
                                           . == "Lambda (C.37)" ~ "Lambda",
                                           TRUE ~ .)) %>%
      mutate_at(vars(Group1), ~ str_replace_all(., "Pf", "PZ"))
      
    #table(NAb_censored_raw$Group1)
    
  # Ensure columns for calculations are numeric
    
    NAb_uncensored[,7:17] <- sapply(NAb_uncensored[,7:17], as.numeric)
    NAb_censored_raw$Titre_log <- as.numeric(NAb_censored_raw$Titre_log)
    NAb_censored_raw$Day_measure <- as.numeric(NAb_censored_raw$Day_measure)

  # Convert median and IQR to mean and SD 
  # Conversions from Wan 2014 https://doi.org/10.1186/1471-2288-14-135
    
    NAb_uncensored <-  NAb_uncensored %>%
      mutate_at(vars(NAb_mean), ~ ifelse(is.na(.), (NAb_q1 + NAb_med + NAb_q3) /
                                           3, .)) %>%
      mutate_at(vars(NAb_SD),
                ~ case_when(
                  !is.na(NAb_l95) ~ (NAb_u95 - NAb_l95) * sqrt(Num_ppl) / (2 * qnorm(0.975)),
                  !is.na(NAb_q1) ~ (NAb_q3 - NAb_q1) /
                    (2 * qnorm((0.75 * Num_ppl - 0.125) / (Num_ppl + 0.25))),
                  !is.na(NAb_min) ~ (NAb_max - NAb_min) /
                    (2 * qnorm((Num_ppl - 0.375) / (Num_ppl + 0.25))),
                  !is.na(.) ~ .
                )) %>%
      select(1:10, 18)
    
    
### 2. Censoring for studies with individual data ####

  # For grouping according to day_measured:
  # NAb_censored_raw_summary <- NAb_censored_raw %>% group_by(Study_ID, Day_measure) %>% summarise(n = n())
  # Set cut-off for ST/LT measurements within the same study at 80 days
    
    NAb_censored_raw <- NAb_censored_raw %>%
      mutate(Day_order = case_when(
        Day_measure > 300 ~ "vLT",
        Day_measure > 80 ~ "LT",
        TRUE ~ "ST"
      ))
    
    NAb_uncensored <- NAb_uncensored %>%
      mutate(Day_order = case_when(
        Day_measure > 300 ~ "vLT",
        Day_measure > 80 ~ "LT",
        TRUE ~ "ST"
      ))
      
    
  # Negative log likelihood of normal distribution model with censoring (adapted from Cromer & Khoury github)
  
    Likelihood = function(p,censT,data) {
      -sum(log(dnorm(data[data>censT],p[1],p[2]))) - sum(log(pnorm(data[data<=censT],p[1],p[2])))}
    
    NAb_censored <- unique(NAb_censored_raw[,c("Study_ID","Group1","Variant","Day_order")])

    NAb_censored$SD <- NA
    NAb_censored$EstimatedMean <- NA
    NAb_censored$Num_ppl <- NA
    
    for (i in 1:nrow(NAb_censored)){
      
      tempdata =  NAb_censored_raw$Titre_log[NAb_censored_raw$Study_ID==NAb_censored$Study_ID[i] & 
                                              NAb_censored_raw$Group1==NAb_censored$Group1[i] & 
                                              NAb_censored_raw$Variant==NAb_censored$Variant[i] &
                                              NAb_censored_raw$Day_order==NAb_censored$Day_order[i]]
      if(sd(tempdata) != 0 & length(tempdata) != 1){
        fitmdltemp <- nlm(function(p){Likelihood(p,NAb_uncensored$LOD[NAb_uncensored$Study_ID==NAb_censored$Study_ID[i] &
                                                                        NAb_uncensored$Group1==NAb_censored$Group1[i] &
                                                                        NAb_uncensored$Variant==NAb_censored$Variant[i] &
                                                                        NAb_uncensored$Day_order==NAb_censored$Day_order[i]],
                                                 tempdata)},c(mean(tempdata),sd(tempdata)))
      
      
      NAb_censored$SD[i] <- fitmdltemp$estimate[2]
      NAb_censored$EstimatedMean[i] <- fitmdltemp$estimate[1]
      NAb_censored$Num_ppl[i] <- length(tempdata)
      }
    }
    

### 3. Normalise to convalescent reference ####
  # Note: currently only studies for which convalescent data is WT are included (NR not included)
    
  # Function to calculate fold of reference, for uncensored data
    
    ## TODO: update to function that doesn't split into separate tables
    conv_fold_uncensored <- function(tbl, ref_tbl) {
      
      tbl <- tbl %>%
        inner_join(ref_tbl, by = "Study_ID") %>%  # only keeps rows that have corresponding reference
        mutate(NAb_fold = NAb_mean.x - NAb_mean.y) %>%
        mutate(SEM = sqrt((NAb_SD.x^2 / Num_ppl.x) + (NAb_SD.y^2 / Num_ppl.y))) %>%
        select(Study_ID, Type.x, Group1.x, Group2.x, Variant.x, Day_measure.x, Day_order.x, NAb_fold, SEM, Num_ppl.x) %>%
        rename_with(~str_remove(., '.x')) %>% # remove suffix .x from column names
        filter(!is.na(SEM)) # remove values with missing SEM

      return(tbl)
    }
    
  
  # Normalise uncensored data
    
    # create convalescent table
    conv_uncensored_all <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & Variant == 'WT')
    conv_uncensored_LOD <- conv_uncensored_all %>% filter(LOD_50 == "N") # only includes data with >50% points above LOD
    
    # Vaccination response, all data
    NAb_uncensored_fold_all <- NAb_uncensored %>% 
      filter(Group1 != "ConvWT" & Group1 != "ConvNR") %>% 
      conv_fold_uncensored(conv_uncensored_all) %>%
      arrange(Group1)
    
    # Vaccination response, only data with >50% points above LOD
    NAb_uncensored_fold_LOD <- NAb_uncensored %>% 
      filter(Group1 != "ConvWT" & Group1 != "ConvNR" & LOD_50 == "N") %>% 
      conv_fold_uncensored(conv_uncensored_LOD) %>%
      arrange(Group1)
    
    # Convalescent WT response to variants, all data
    NAb_uncensored_fold_ConvtoVar_all <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & Variant != 'WT') %>% 
      conv_fold_uncensored(conv_uncensored_all)%>%
      arrange(Group1)
    
    # Convalescent WT response to variants, only data with >50% points above LOD
    NAb_uncensored_fold_ConvtoVar_LOD <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & Variant != 'WT' & LOD_50 == "N") %>% 
      conv_fold_uncensored(conv_uncensored_LOD) %>%
      arrange(Group1)
    
    
  # Normalise censored data (subtraction, as per Khoury & Cromer)
    
    conv_censored <- NAb_censored %>% 
      filter(Group1 == "ConvWT" & Variant == 'WT')
    
    # Vaccination response, censored data
    NAb_censored_fold <- NAb_censored %>% 
      filter(Group1 != "ConvWT" & Group1 != "ConvNR") %>%
      inner_join(conv_censored, by = "Study_ID") %>% # only keeps rows that have corresponding reference
      mutate(NAb_fold = EstimatedMean.y - EstimatedMean.x) %>%
      mutate(SEM = sqrt((SD.x^2 / Num_ppl.x) + (SD.y^2 / Num_ppl.y))) %>%
      select(Study_ID, Group1.x, Variant.x, NAb_fold, SEM, Num_ppl.x, Day_order.x) %>%
      rename_with(~str_remove(., '.x')) %>% # remove suffix .x from column names
      arrange(Group1)
    
    # Convalescent WT response to variants, censored data
    NAb_censored_fold_ConvtoVar <- NAb_censored %>%
      filter(Group1 == "ConvWT" & Variant != 'WT') %>%
      inner_join(conv_censored, by = "Study_ID") %>% # only keeps rows that have corresponding reference
      mutate(NAb_fold = EstimatedMean.y - EstimatedMean.x) %>%
      mutate(SEM = sqrt((SD.x^2 / Num_ppl.x) + (SD.y^2 / Num_ppl.y))) %>%
      select(Study_ID, Group1.x, Variant.x, NAb_fold, SEM, Num_ppl.x, Day_order.x) %>%
      rename_with(~str_remove(., '.x')) %>% # remove suffix .x from column names
      arrange(Group1)


### 4. Meta-analysis #####
    
  # Pooled SD for censored data
  # !!! currently only for vaccination data - can update for convalescent variant responses #
    
    # Centre individual points at the reported mean
    
    NAb_censored_centred <- NAb_censored_raw %>%
      
      # fold of convalescent for individual points
      filter(Group1 != 'ConvWT') %>%
      inner_join((filter(NAb_censored_raw, Group1 == 'ConvWT' & Variant == 'WT')), by = "Study_ID") %>%
      mutate(fold_conv = Titre_log.y - Titre_log.x) %>%
      rename_with(~str_remove(., '.x')) %>%
      
      # centre individual points to normalised reported mean
      inner_join(NAb_uncensored_fold_all, by = c("Study_ID", "Group1", "Variant", "Day_measure")) %>%
      mutate(centred_data = fold_conv - NAb_fold) %>%
    
      # centre normalised LOD to normalised reported mean
      mutate(centred_LOD = LOD - Titre_log - NAb_fold) %>%
      select("Study_ID", "Group1", "Day_measure", "Variant", "centred_data","centred_LOD") 
      
    
    # Calculate pooled SD (function from Khoury & Cromer Github)
    
    PooledSDModelFit <- nlm(function(p) {Likelihood(p, NAb_censored_centred$centred_LOD, NAb_censored_centred$centred_data)},
                            c(mean(NAb_censored_centred$centred_data), sd(NAb_censored_centred$centred_data)))    
    SD_pooled <- PooledSDModelFit$estimate[2]
    SD_pooled
    
    
  # Summary NAb ratio stratified by vaccine combination
    # using uncensored data
    
    #table(NAb_uncensored_fold_LOD$Group1, NAb_uncensored_fold_LOD$Day_order)
    
    meta_uncensored <- function(NAb_uncensored_fold, type, variant = "WT", duration = "ST") {
      
      # NAb_uncensored_fold is ratio of neutralising antibody against convalescent antibody and removing studies with >50% points below LOD
      # hybrid is 1 (vaccine only), 2 (hybrid only), or 3 (overall)
      # variant 
      # duration is "ST" (< 80 days), "LT" (80-300 days), or "vLT" (> 300 days)
      
      if (type == 1) {
        tempdata <- NAb_uncensored_fold %>%
          filter(!str_starts(Group1, "Conv"))
      } else if (type == 2) {
        tempdata <- NAb_uncensored_fold %>%
          filter(str_starts(Group1, "Conv"))
      } else if (type == 3) {
        tempdata <- NAb_uncensored_fold
      }
      
      tempdata <- tempdata %>%
        separate(Group1, into = c("dose1","dose2","dose3","dose4"), sep = "_", remove = F) %>%
        mutate_at(c("dose1","dose2","dose3","dose4"), ~ str_replace_all(., "[:digit:]", "")) %>%
        filter((dose1 != "Conv" & !is.na(dose2)) | (dose1 == "Conv" & !is.na(dose3))) %>% # filter out single dose regimens
        unite(Group_meta, c(dose1, dose2, dose3, dose4), sep = "_") %>% # create variable for subgroup meta-analysis
        mutate(Group_meta = str_replace_all(Group_meta, "_NA", "")) %>%
        filter(!str_starts(Group1, "JJ") & # remove one study with mixed vaccination of JJ (zero market share in Thailand) and PZ
                 Variant == variant & Day_order == duration)
      
      meta <- metamean(
        n = tempdata$Num_ppl,
        mean = tempdata$NAb_fold,
        sd = tempdata$SEM,
        random = TRUE,
        subgroup = tempdata$Group_meta,
        studlab = tempdata$Group_meta
      )
      
      # overall_tbl <- data.frame(
      #   group = "Overall",
      #   mean_ab_ratio = unlist(meta['TE.random'], use.names = FALSE),
      #   se_ab_ratio = unlist(meta['seTE.random'], use.names = FALSE),
      #   n_study = sum(unlist(meta['k.w'], use.names = FALSE)),
      #   I2 = unlist(meta['I2'], use.names = FALSE))
      
      subgroup_tbl <- data.frame(
        group = unique(tempdata$Group_meta),
        mean_ab_ratio = unlist(meta['TE.random.w'], use.names = FALSE),
        se_ab_ratio = unlist(meta['seTE.random.w'], use.names = FALSE),
        n_study = unlist(meta['k.w'], use.names = FALSE),
        I2 = unlist(meta['I2.w'], use.names = FALSE))
      
      # result_tbl <- rbind(overall_tbl, subgroup_tbl)
      return(list(meta,subgroup_tbl))
    }
    
# short-term against WT only, removing points with >50% points below LOD
tbl_ab_ratio <- list()
tbl_ab_ratio$vac_wt_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 1, variant = "WT", duration = "ST")[[2]]
tbl_ab_ratio$hyb_wt_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 2, variant = "WT", duration = "ST")[[2]]
tbl_ab_ratio$all_wt_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "WT", duration = "ST")[[2]]
tbl_ab_ratio$vac_wt_lt <- meta_uncensored(NAb_uncensored_fold_LOD, type = 1, variant = "WT", duration = "LT")[[2]]
tbl_ab_ratio$hyb_wt_lt <- meta_uncensored(NAb_uncensored_fold_LOD, type = 2, variant = "WT", duration = "LT")[[2]]
tbl_ab_ratio$all_wt_lt <- meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "WT", duration = "LT")[[2]]

# for sanity check: short-term against alpha, beta, and omicron, removing points with >50% points below LOD
tbl_ab_ratio$all_alpha_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "Alpha", duration = "ST")[[2]]
tbl_ab_ratio$all_beta_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "Beta", duration = "ST")[[2]]
tbl_ab_ratio$all_omicron_st <- meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "Omicron", duration = "ST")[[2]]

# By group meta-analysis diagnostics
# tiff("Meta_NAb_uncensored_fold_LOD.tiff", height = 15, width = 12, res = 300, units = "in", compression = "lzw")
# forest(meta_uncensored(NAb_uncensored_fold_LOD, type = 3, variant = "WT", duration = "ST")[[1]])
# dev.off()


### 5. Convert NAb titres to vaccine efficacy ####
    # (covariate matrix from Deborah, email on 11 Oct 2022)
    
    n_vx <- 246 # number of samples from Khoury & Cromer paper
    mean_k_n50 <- c(1.1306607, -0.6966127)
    cov_k_n50 <- matrix(c(0.03106460, 0.010755914, 0.010755914, 0.005727749), nrow=2, ncol=2)
    
  # function to transform NAb data to effectiveness DETERMINISTIC
    vx_efficacy <- function(k, n, n50) {
      vx_eff <- 1/(1 + exp(-k*(n - n50)))
      return(vx_eff)
    }
    
  # function to transform NAb data to effectiveness PROBABILISTIC
    runs <- 1000 # number of Monte Carlo runs
    k_n50_samples <- mvtnorm::rmvnorm(runs, mean = mean_k_n50, sigma = cov_k_n50)
    colnames(k_n50_samples) <- c("k","IC50")
    
    tbl_vx_eff <- list()
    tbl_vx_eff$all_wt_st <- tbl_ab_ratio$all_wt_st
    tbl_vx_eff$all_alpha_st <- tbl_ab_ratio$all_alpha_st
    tbl_vx_eff$all_beta_st <- tbl_ab_ratio$all_beta_st
    tbl_vx_eff$all_omicron_st <- tbl_ab_ratio$all_omicron_st
    tbl_vx_eff$all_wt_lt <- tbl_ab_ratio$all_wt_lt
    
    tbl_vx_eff <- map(tbl_vx_eff, ~ .x %>% 
                    mutate(mean_vx_eff = NA,
                           sd_vx_eff = NA))
    
    n_loop <- 1 # index for storing Monte Carlo results

    tbl_vx_eff_monte_carlo <- list()
    tbl_ab_ratio_random <- list()

  # loop to transform NAb data to effectiveness PROBABILISTIC for each vaccine combination
    for (i in 1:length(tbl_vx_eff)){
      for(j in 1:nrow(tbl_vx_eff[[i]])){
        n_random <- rnorm(runs, mean = tbl_vx_eff[[i]][j,"mean_ab_ratio"], sd = tbl_vx_eff[[i]][j,"se_ab_ratio"])
        monte_carlo <- mapply(vx_efficacy, k = k_n50_samples[,1], n = n_random, n50 = k_n50_samples[,2])
        
        tbl_vx_eff[[i]][j,"mean_vx_eff"] <- mean(monte_carlo)
        tbl_vx_eff[[i]][j,"sd_vx_eff"] <- sd(monte_carlo)
        
        tbl_vx_eff_monte_carlo[[n_loop]] <- monte_carlo
        names(tbl_vx_eff_monte_carlo)[n_loop] <- paste("monte",names(tbl_vx_eff)[i],tbl_vx_eff[[i]][j, "group"], sep = "_")
        
        tbl_ab_ratio_random[[n_loop]] <- n_random
        names(tbl_ab_ratio_random)[n_loop] <- paste("ab",names(tbl_vx_eff)[i],tbl_vx_eff[[i]][j, "group"], sep = "_")
        
        n_loop <- n_loop + 1
        }
    }
    
    tbl_vx_eff # summary table of NAb titre ratio and vaccine efficacy for each vaccine combo
    
    df_ab_ratio_random <- ldply(tbl_ab_ratio_random, data.frame)
    names(df_ab_ratio_random) <- c("regimen","NAb_titre")
    df_ab_ratio_random <- df_ab_ratio_random %>% 
      mutate(id = 1:n())
    
    df_vx_eff_monte_carlo <- ldply(tbl_vx_eff_monte_carlo, data.frame)
    names(df_vx_eff_monte_carlo) <- c("regimen","vx_eff")
    df_vx_eff_monte_carlo <- df_vx_eff_monte_carlo %>% 
      mutate(id = 1:n())
    
    # tempdata_ab <- df_ab_ratio_random %>% 
    #   filter(regimen %in% c("all_wt_st_PZ_PZ", "all_wt_st_PZ_PZ_PZ"))
    # 
    #   ggplot(tempdata_ab, aes(x = NAb_titre, fill = regimen, col = regimen)) +
    #   geom_histogram(alpha = 0.7) +
    #     facet_wrap( ~ regimen)
    # 
    # t.test(tempdata_ab$NAb_titre ~ tempdata_ab$regimen)
    
  # Density figure of vaccine combo with antibody titre on the x-axes    
  # Sanity check: a. pz > everything else; b. vac + conv > vac
    #tiff(file = "NAb titre density.tiff", unit = "in", width = 8, height = 6, res = 300, compression = "lzw")
    df_ab_ratio_random %>%
      filter(startsWith(regimen, "ab_all_wt_st_")) %>%
      mutate(regimen = str_remove(regimen, "ab_all_wt_st_")) %>%
      ggplot(., aes(x = NAb_titre, fill = regimen, col = regimen)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = c("#FED976", "#FEB24C", "#FD8D3C", "#C6DBEF", "#9ECAE1", "#6BAED6", "#D5F591", "#8CDBA9", "#33AB5F", "#BDBDBD", "#969696"), 
                        breaks = c("AZ_AZ", "SP_SP", "PZ_PZ", "AZ_PZ", "AZ_SP", "SP_AZ", "AZ_AZ_PZ", "PZ_PZ_PZ",  "MD_MD_PZ", "Conv_PZ_PZ", "Conv_PZ_PZ_PZ")) +
      scale_color_manual(values = c("#FED976", "#FEB24C", "#FD8D3C",  "#C6DBEF", "#9ECAE1", "#6BAED6", "#D5F591", "#8CDBA9", "#33AB5F", "#BDBDBD", "#969696"), 
                        breaks = c("AZ_AZ", "SP_SP", "PZ_PZ", "AZ_PZ", "AZ_SP", "SP_AZ", "AZ_AZ_PZ", "PZ_PZ_PZ", "MD_MD_PZ", "Conv_PZ_PZ", "Conv_PZ_PZ_PZ"))
    #dev.off()
    
  # Sanity check: c. day_measure; short elapse > long elapse
    df_vx_eff <- ldply(tbl_vx_eff_monte_carlo, data.frame)
    names(df_vx_eff) <- c("regimen","vx_eff")
    
    #tiff(file = "vx_eff by day_order.tiff", unit = "in", width = 8, height = 6, res = 300, compression = "lzw")
    df_vx_eff %>%
      mutate(Day_order = if_else(str_detect(regimen, "lt"),"LT","ST")) %>%
      mutate(regimen = str_remove(regimen, "monte_all_wt_st_")) %>%
      mutate(regimen = str_remove(regimen, "monte_all_wt_lt_")) %>%
      filter(regimen %in% c("AZ_AZ", "Conv_PZ_PZ", "PZ_PZ", "SP_SP")) %>%
      ggplot(., aes(x = vx_eff, fill = Day_order, col = Day_order)) +
      geom_density(alpha = 0.7) +
      facet_wrap(~ regimen)
    #dev.off()

  # Sanity check: d. variant; short-term against WT vs Alpha, Beta, Omicron
    #tiff(file = "vx_eff by variant.tiff", unit = "in", width = 9, height = 6, res = 300, compression = "lzw")
    df_vx_eff %>%
      mutate(regimen = str_remove(regimen, "monte_all_")) %>%
      filter(regimen %in% c("wt_st_PZ_PZ", "wt_st_Conv_PZ_PZ", "alpha_st_PZ_PZ", "alpha_st_Conv_PZ_PZ", "beta_st_PZ_PZ", "beta_st_Conv_PZ_PZ", "omicron_st_PZ_PZ", "omicron_st_Conv_PZ_PZ")) %>%
      separate(regimen, into = c("Variant","Day_order", "regimen"), sep = "_", remove = T) %>%
      mutate(regimen = if_else(regimen == "PZ", "PZ_PZ", "Conv_PZ_PZ")) %>%
      ggplot(., aes(x = vx_eff, fill = Variant, col = Variant)) +
      geom_density(alpha = 0.7) +
      facet_wrap(~ regimen)
    #dev.off()
    
    
### 6. Set cut-point for NAb ####
    # TODO: 1. cut off for L, M, H; 2. percentage transition for v_b and v_p
    
  # set the cut-points for NAb ratio at 0 and 0.5
    df_ab_ratio_level <- df_ab_ratio_random %>%
      mutate(immune_level = case_when(
        NAb_titre < 0 ~ "Vl",
        NAb_titre <= 0.5 ~ "Vm",
        TRUE ~ "Vh")
      ) %>%
      mutate_at(vars(immune_level), ~factor(., levels = c("Vl","Vm","Vh")))
    
    table(df_ab_ratio_level$immune_level)
    
    tbl_immune_level <- df_ab_ratio_level %>%
      left_join(df_vx_eff_monte_carlo %>% select(-regimen), by = "id") %>%
      group_by(immune_level) %>%
      summarise(mean = mean(vx_eff),
                sd = sd(vx_eff))
    tbl_immune_level # summary table of vaccine efficacy in people with low, medium, and high protective efficacy
    
    #tiff(file = "vx_eff by level of protection.tiff", unit = "in", width = 9, height = 6, res = 300, compression = "lzw")
    df_ab_ratio_level %>%
      left_join(df_vx_eff_monte_carlo %>% select(-regimen), by = "id") %>%
      #filter(regimen == "monte_all_wt_st_Conv_PZ_PZ" | regimen == "monte_all_wt_st_PZ_PZ") %>%
      ggplot(., aes(x = vx_eff, fill = immune_level, col = immune_level)) +
      geom_density(alpha = 0.5)
    #dev.off()

    # prob of people going into vl, vm, vh for each vaccine combo
    tbl_prob_immune_level <- df_ab_ratio_level %>%
      left_join(df_vx_eff_monte_carlo %>% select(-regimen), by = "id") %>%
      filter(startsWith(regimen, "ab_all_wt_st_")) %>%
      mutate(regimen = str_remove(regimen, "ab_all_wt_st_")) %>%
      group_by(regimen, immune_level) %>%
      summarise(p = n()/1000)
    tbl_prob_immune_level
    #write_csv(tbl_prob_immune_level, file = "tbl_prob_immune_level.csv")
    
    # plot prob density of people going into vl, vm, vh for each vaccine combo
    #tiff(file = "vx_eff by combo.tiff", unit = "in", width = 12, height = 10, res = 300, compression = "lzw")
    df_ab_ratio_level %>%
      left_join(df_vx_eff_monte_carlo %>% select(-regimen), by = "id") %>%
      filter(startsWith(regimen, "ab_all_wt_st_")) %>%
      mutate(regimen = str_remove(regimen, "ab_all_wt_st_")) %>%
      mutate_at(vars(regimen), ~factor(., levels = c("AZ_AZ", "SP_SP", "PZ_PZ", "AZ_PZ", "AZ_SP", "SP_AZ", "AZ_AZ_PZ", "PZ_PZ_PZ", "MD_MD_PZ", "Conv_PZ_PZ", "Conv_PZ_PZ_PZ"))) %>%
      ggplot(., aes(x = vx_eff, y = fct_rev(regimen), fill = stat(x))) +
      geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
      scale_fill_viridis_c(name = "Vaccine efficacy", option = "D") +
      theme_bw(base_size = 18)
    #dev.off()
    

### ASSUMPTIONS ###
    
      # all NAb data is normally distributed


### OUTSTANDING TASKS ###
      
      # meta-regression for time component with exponential decay
      # comparison of censored/non-censored results
      # include >50% below LOD as Boolean variable
