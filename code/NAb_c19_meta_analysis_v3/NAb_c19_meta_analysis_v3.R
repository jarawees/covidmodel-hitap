# COVID-19 neutralising antibody titre meta-analysis
# 18 November 2022
# Version 4

# Load libraries
if(!require(pacman)) install.packages("pacman")
library(pacman)
p_load(dplyr, stats, tidyverse, meta)

# Set directory to same directory as the r-script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### 1. Prepare data ####

  # Read CSV files
  rm(list=ls())

    # Uncensored: all included studies
    csv_uncensored <- c("Anna_NAb_meta1.csv", "Dimple_NAb_meta1_v2.csv", "Dream_NAb_meta1_v2.csv", 
                        "Ja_NAb_meta1_PZ35.csv", "Pui_NAb_meta1.csv", "Siobhan_NAb_meta1.csv", 
                        "Saudamini_NAb_meta1.csv")
    csv_uncensored <- lapply(csv_uncensored, read.csv)
    csv_uncensored <- do.call("rbind", csv_uncensored)

    # Censored: subset of uncensored studies with individual data point    
    csv_censored <- c("Dimple_NAb_meta2.csv", "Dream_NAb_meta2_v2.csv", "Ja_NAb_meta2_NODELTA.csv", 
                      "Pui_NAb_meta2.csv", "Siobhan_NAb_meta2_Pf5.csv")
    csv_censored <- lapply(csv_censored, read.csv) 
    csv_censored <- do.call("rbind", csv_censored)
    
    NAb_uncensored <- csv_uncensored %>%
      select(-DOI)
    NAb_censored_raw <- csv_censored %>%
      select(-DOI)
    
  # Replace values below LOD with half of LOD
    
    NAb_uncensored$LOD <- as.numeric(NAb_uncensored$LOD)
    NAb_uncensored <- NAb_uncensored %>%
      mutate_at(vars(LOD), ~ ifelse(. > 2, log10(.), .)) %>%
      mutate_at(vars(NAb_mean, NAb_SD, NAb_l95, NAb_u95, NAb_med, NAb_q1, NAb_q3, NAb_min, NAb_max), ~ ifelse(. == "<LOD", LOD/2, .)) %>%
      mutate_at(vars(Variant), ~ case_when(. == "Wuhan" ~ "WT",
                                           . == "D614" ~ "D614G",
                                           . == "kappa" ~ "Kappa",
                                           . == "Lambda " ~ "Lambda",
                                           . == "Epsilion" ~ "Epsilon",
                                           . == "theta" ~ "Theta",
                                           TRUE ~ .)) %>%
      filter(Variant != "SARS-COV-1") %>%
      mutate_at(vars(Group1), ~ str_replace_all(., c("Pf" = "PZ", "PF" = "PZ",  "CV" = "SV")))
    
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
      mutate_at(vars(NAb_mean), ~ifelse(is.na(.), (NAb_q1 + NAb_med + NAb_q3)/3, .)) %>%
      mutate_at(vars(NAb_SD), ~case_when(
        !is.na(NAb_l95) ~ (NAb_u95 - NAb_l95)*sqrt(Num_ppl)/(2*qnorm(0.975)),
        !is.na(NAb_q1) ~ (NAb_q3 - NAb_q1)/(2*qnorm((0.75*Num_ppl - 0.125)/(Num_ppl + 0.25))),
        !is.na(NAb_min) ~ (NAb_max - NAb_min)/(2*qnorm((Num_ppl - 0.375)/(Num_ppl + 0.25))),
        !is.na(.) ~ .)) %>%
      select(1:10,18)
    
    
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
    # using uncensored data, short-term against WT only, removing points with >50% points below LOD
    
    #table(NAb_uncensored_fold_LOD$Group1, NAb_uncensored_fold_LOD$Day_order)
    
    meta_uncensored_WT_ST <- function(tbl, type = 3) {
      
      # hybrid is 1 (vaccine only), 2 (hybrid only), or 3 (overall)
      if (type == 1) {
        tempdata <- tbl %>%
          filter(!str_starts(Group1, "Conv")) %>%
          mutate(Group_meta = case_when(
            # create new group variable for subgroup meta-analysis
            str_length(Group1) == 2 ~ Group1,
            str_length(Group1) == 5 &
              (substr(Group1, 1, 2) == substr(Group1, 4, 5)) ~ Group1,
            str_length(Group1) == 7 &
              (substr(Group1, 1, 2) == substr(Group1, 4, 5)) ~ substr(Group1, 1, 5),
            str_length(Group1) > 7 ~ "BoosterPZ",
            TRUE ~ "MixedVac")
          )
        
      } else if (type == 2) {
        tempdata <- tbl %>%
          filter(str_starts(Group1, "Conv")) %>%
          mutate(Group_meta = case_when(
            # create new group variable for subgroup meta-analysis
            str_length(Group1) == 2+5 ~ Group1,
            str_length(Group1) == 5+5 &
              (substr(Group1, 1+5, 2+5) == substr(Group1, 4+5, 5+5)) ~ Group1,
            str_length(Group1) == 7+5 &
              (substr(Group1, 1+5, 2+5) == substr(Group1, 4+5, 5+5)) ~ substr(Group1, 1, 5+5),
            str_length(Group1) > 7+5 ~ "Conv_BoosterPZ",
            TRUE ~ "MixedVac")
          )
        
      } else if (type == 3) {
        tbl1 <- tbl %>%
          filter(!str_starts(Group1, "Conv")) %>%
          mutate(Group_meta = case_when(
            # create new group variable for subgroup meta-analysis
            str_length(Group1) == 2 ~ Group1,
            str_length(Group1) == 5 &
              (substr(Group1, 1, 2) == substr(Group1, 4, 5)) ~ Group1,
            str_length(Group1) == 7 &
              (substr(Group1, 1, 2) == substr(Group1, 4, 5)) ~ substr(Group1, 1, 5),
            str_length(Group1) > 7 ~ "BoosterPZ",
            TRUE ~ "MixedVac")
          )
        
        tbl2 <- tbl %>%
          filter(str_starts(Group1, "Conv")) %>%
          mutate(Group_meta = case_when(
            # create new group variable for subgroup meta-analysis
            str_length(Group1) == 2+5 ~ Group1,
            str_length(Group1) == 5+5 &
              (substr(Group1, 1+5, 2+5) == substr(Group1, 4+5, 5+5)) ~ Group1,
            str_length(Group1) == 7+5 &
              (substr(Group1, 1+5, 2+5) == substr(Group1, 4+5, 5+5)) ~ substr(Group1, 1, 5+5),
            str_length(Group1) > 7+5 ~ "Conv_BoosterPZ",
            TRUE ~ "MixedVac")
            )
        
        tempdata <- rbind(tbl1, tbl2)
      }
      
      tempdata <- tempdata %>%
        filter(!str_starts(Group1, "JJ") & # remove one study with mixed vaccination of JJ and PZ
                 Variant == "WT" & Day_order == "ST")
      
      meta <- metamean(
        n = tempdata$Num_ppl,
        mean = tempdata$NAb_fold,
        sd = tempdata$SEM,
        random = TRUE,
        subgroup = tempdata$Group_meta,
        studlab = tempdata$Group_meta
      )
      
      overall_tbl <- data.frame(
        group = "Overall",
        mean_ab_ratio = unlist(meta['TE.random'], use.names = FALSE),
        se_ab_ratio = unlist(meta['seTE.random'], use.names = FALSE)
        # l95 = unlist(meta['TE.random'], use.names = FALSE) - qnorm(0.975) *
        #   unlist(meta['seTE.random'], use.names = FALSE),
        # u95 = unlist(meta['TE.random'], use.names = FALSE) + qnorm(0.975) *
        #   unlist(meta['seTE.random'], use.names = FALSE)
      )
      
      subgroup_tbl <- data.frame(
          group = unique(tempdata$Group_meta),
          mean_ab_ratio = unlist(meta['TE.random.w'], use.names = FALSE),
          se_ab_ratio = unlist(meta['seTE.random.w'], use.names = FALSE)
          # l95 = unlist(meta['TE.random.w'], use.names = FALSE) - qnorm(0.975) *
          #   unlist(meta['seTE.random.w'], use.names = FALSE),
          # u95 = unlist(meta['TE.random.w'], use.names = FALSE) + qnorm(0.975) *
          #   unlist(meta['seTE.random.w'], use.names = FALSE)
        )
      
      result_tbl <- rbind(overall_tbl, subgroup_tbl)
      return(result_tbl)
    }
    
tbl_ab_ratio <- meta_uncensored_WT_ST(NAb_uncensored_fold_LOD,3)
    

### 5. Convert NAb titres to vaccine efficacy ####
    # (covariate matrix from Deborah, email on 11 Oct)
    
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
    
    tbl_vx_eff <- tbl_ab_ratio
    tbl_vx_eff$mean_vx_eff <- NA
    tbl_vx_eff$sd_vx_eff <- NA
    
    tbl_vx_eff_monte_carlo <-
      data.frame(
        Overall = rep(NA, 1000),
        AZ = rep(NA, 1000),
        AZ_AZ = rep(NA, 1000),
        BoosterPZ = rep(NA, 1000),
        MixedVac = rep(NA, 1000),
        PZ = rep(NA, 1000),
        PZ_PZ = rep(NA, 1000),
        SP_SP = rep(NA, 1000),
        Conv_PZ = rep(NA, 1000),
        Conv_PZ_PZ = rep(NA, 1000),
        Conv_BoosterPZ = rep(NA, 1000)
      )

  # loop to transform NAb data to effectiveness PROBABILISTIC for each vaccine combination
    for (i in 1:nrow(tbl_vx_eff)){
      n_random <- rnorm(runs, mean = tbl_vx_eff[i,2], sd = tbl_vx_eff[i,3])
      monte_carlo <- mapply(vx_efficacy, k = k_n50_samples[,1], n = n_random, n50 = k_n50_samples[,2])
      
      tbl_vx_eff$mean_vx_eff[i] <- mean(monte_carlo)
      tbl_vx_eff$sd_vx_eff[i] <- sd(monte_carlo)
      
      tbl_vx_eff_monte_carlo[,i] <- monte_carlo
    }
    
    tbl_vx_eff
    
    #tiff(file = "tbl_vx_eff_monte_carlo.tiff", unit = "in", width = 5, height = 12, res = 300)
    tbl_vx_eff_monte_carlo %>%
      gather(key = "group", value = "vx_eff") %>%
      mutate_at(vars(group), ~factor(., levels = c("AZ", "PZ", "AZ_AZ", "SP_SP", "PZ_PZ", "MixedVac", "BoosterPZ", "Conv_PZ", "Conv_PZ_PZ", "Conv_BoosterPZ", "Overall"))) %>%
      ggplot(., aes(x = vx_eff)) +
      geom_histogram(binwidth = 0.01) +
      facet_wrap(~ group, ncol = 1)
    #dev.off()
      
    
### ASSUMPTIONS ###
    
      # all NAb data is normally distributed


### OUTSTANDING TASKS ###
      
      # meta-regression for time component with exponential decay
      # comparison of censored/non-censored results
      # include >50% below LOD as Boolean variable
