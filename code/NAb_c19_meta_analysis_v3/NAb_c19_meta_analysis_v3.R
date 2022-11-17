# COVID-19 neutralising antibody titre meta-analysis
# 25 October 2022
# Version 2


### Load libraries

  list_pkg <- c("dplyr","stats","tidyverse","meta")
  for (i in seq.int (length (list_pkg) ) ) {
    if (! require (list_pkg[i], character.only = T )) {
      install.packages (list_pkg[i]) }
    library(list_pkg[i] , character.only = T)
  }
  


### 1. Prepare data ###

  # Read CSV files
  # Note: update with Saudamini files and revised files from Dream & Ja & Siobhan
    
    csv_uncensored <- c("Anna_NAb_meta1.csv", "Dimple_NAb_meta1_v2.csv",
                        "Dream_NAb_meta1.csv", "Ja_NAb_meta1_PZ35.csv",
                        "Pui_NAb_meta1.csv", "Siobhan_NAb_meta1.csv")
    csv_uncensored <- lapply(csv_uncensored, read.csv) 
    csv_uncensored <- do.call("rbind", csv_uncensored)
    
    csv_censored <- c("Dimple_NAb_meta2.csv","Ja_NAb_meta2_NODELTA.csv",
                        "Pui_NAb_meta2.csv", "Siobhan_NAb_meta2_Pf5.csv")
    csv_censored <- lapply(csv_censored, read.csv) 
    csv_censored <- do.call("rbind", csv_censored)
    
    NAb_uncensored <- csv_uncensored %>%
      select(-DOI)
    NAb_censor_data <- csv_censored %>%
      select(-DOI)
    
    
  # Replace values below LOD with half of LOD
    
    NAb_uncensored$LOD <- as.numeric(NAb_uncensored$LOD)
    NAb_uncensored <- NAb_uncensored %>%
      mutate_at(vars(LOD), ~ ifelse(. > 2, log10(.), .)) %>%
      mutate_at(vars(NAb_med,NAb_q1,NAb_q3,NAb_min), ~ ifelse(. == "<LOD", LOD/2, .))
    
    NAb_censor_data <- NAb_censor_data %>%
      left_join(NAb_uncensored %>% select(Study_ID, LOD) %>% unique(), by = "Study_ID") %>%
      mutate_at(vars(LOD), ~ ifelse(. > 2, log10(.), .)) %>%
      mutate_at(vars(Titre_log), ~ ifelse(. == "<LOD", LOD/2, .)) 
    
    
  # Ensure columns for calculations are numeric
    NAb_uncensored[,7:17] <- sapply(NAb_uncensored[,7:17], as.numeric)
    NAb_censor_data$Titre_log <- as.numeric(NAb_censor_data$Titre_log)
    NAb_censor_data$Day_measure <- as.numeric(NAb_censor_data$Day_measure)
    
    
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
    
    
    
### 2. Censoring for studies with individual data ###

  # For grouping according to day_measured:
  # NAb_censor_data_summary <- NAb_censor_data %>% group_by(Study_ID, Day_measure) %>% summarise(n = n())
  # Set cut-off for ST/LT measurements within the same study at 80 days
    
    NAb_censor_data <- NAb_censor_data %>%
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
  ### !!! ISSUE: when meta2 has different data to meta1 (e.g. delta for PZ39)
  
    Likelihood = function(p,censT,data) {
      -sum(log(dnorm(data[data>censT],p[1],p[2]))) - sum(log(pnorm(data[data<=censT],p[1],p[2])))}
    
    NAb_censored <- unique(NAb_censor_data[,c("Study_ID","Group1","Variant","Day_order")])

    NAb_censored$SD <- NA
    NAb_censored$EstimatedMean <- NA
    NAb_censored$Num_ppl <- NA
    
    for (i in 1:nrow(NAb_censored)){
      
      tempdata =  NAb_censor_data$Titre_log[NAb_censor_data$Study_ID==NAb_censored$Study_ID[i] & 
                                              NAb_censor_data$Group1==NAb_censored$Group1[i] & 
                                              NAb_censor_data$Variant==NAb_censored$Variant[i] &
                                              NAb_censor_data$Day_order==NAb_censored$Day_order[i]]
      if(sd(tempdata) != 0 & length(tempdata) != 1){
        fitmdltemp <- nlm(function(p){Likelihood(p,NAb_uncensored$LOD[NAb_uncensored$Study_ID==NAb_censored$Study_ID[i] &
                                                                        NAb_uncensored$Group1==NAb_censored$Group1[i] &
                                                                        NAb_uncensored$Variant==NAb_censored$Variant[i] &
                                                                        NAb_uncensored$Day_order==NAb_censored$Day_order[i]], #NAb_censor_data
                                                 tempdata)},c(mean(tempdata),sd(tempdata)))
      
      
      NAb_censored$SD[i] <- fitmdltemp$estimate[2]
      NAb_censored$EstimatedMean[i] <- fitmdltemp$estimate[1]
      NAb_censored$Num_ppl[i] <- length(tempdata)
      }
    }
    

### 3. Normalise to convalescent reference ###
  # Note: currently only studies for which convalescent data is WT are included (NR not included)
    
  # Function to calculate fold of reference, for uncensored data
    
    conv_fold_uncensored <- function(tbl, ref_tbl) {
      
      tbl <- tbl %>%
        inner_join(ref_tbl, by = "Study_ID") %>%
        # only keeps rows that have corresponding reference
        
        mutate(NAb_fold = NAb_mean.x - NAb_mean.y) %>%
        mutate(SEM = sqrt((NAb_SD.x^2 / Num_ppl.x) + (NAb_SD.y^2 / Num_ppl.y))) %>%
        select(Study_ID, Type.x, Group1.x, Group2.x, Variant.x, Day_measure.x, NAb_fold, SEM, Num_ppl.x) %>%
        filter(!is.na(SEM)) %>% # remove values with missing SEM
        mutate(Day_order = case_when(
          Day_measure.x > 80 ~ "LT",
          Day_measure.x <= 80 ~ "ST"
        )) # to be able to separate longand short term immunity for MA
        
        
      return(tbl)
      
    }
    
  
  # Normalise uncensored data
    
    conv_tbl_raw <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & (Variant == 'WT'))
    conv_tbl_LOD <- conv_tbl_raw %>% filter(LOD_50 == "N") # only includes data with >50% points above LOD
    
    # Vaccination response, raw data
    tbl_uncensored_raw <- NAb_uncensored %>% 
      filter(Group1 != "ConvWT") %>% 
      conv_fold_uncensored(conv_tbl_raw) %>%
      arrange(Group1.x)
    
    # Vaccination response, only data with >50% points above LOD
    tbl_uncensored_LOD <- NAb_uncensored %>% 
      filter(Group1 != "ConvWT" & LOD_50 == "N")%>% 
      conv_fold_uncensored(conv_tbl_LOD)%>%
      arrange(Group1.x)
    
    # Convalescent WT response to variants, raw data
    tbl_uncensored_conv_var_raw <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & Variant != 'WT') %>% 
      conv_fold_uncensored(conv_tbl_raw)%>%
      arrange(Group1.x)
    
    # Convalescent WT response to variants, only data with >50% points above LOD
    tbl_uncensored_conv_var_LOD <- NAb_uncensored %>% 
      filter(Group1 == "ConvWT" & Variant != 'WT' & LOD_50 == "N")%>% 
      conv_fold_uncensored(conv_tbl_LOD)%>%
      arrange(Group1.x)
    
    
  # Normalise censored data (subtraction, as per Khoury & Cromer)
    
    conv_censored <- NAb_censored %>% 
      filter(Group1 == "ConvWT" & (Variant == 'WT'))
    
    # Vaccination response, censored data
    tbl_censored <- NAb_censored %>% 
      filter(Group1 != "ConvWT") %>%
      inner_join(conv_censored, by = "Study_ID") %>% # only keeps rows that have corresponding reference
      mutate(NAb_fold = EstimatedMean.y - EstimatedMean.x) %>%
      mutate(SEM = sqrt((SD.x^2 / Num_ppl.x) + (SD.y^2 / Num_ppl.y))) %>%
      select(Study_ID, Group1.x, Variant.x, NAb_fold, SEM, Num_ppl.x)%>%
      arrange(Group1.x)
    
    # Convalescent WT response to variants, censored data
    tbl_censored_conv_var <- NAb_censored %>%
      filter(Group1 == "ConvWT" & Variant != 'WT') %>%
      inner_join(conv_censored, by = "Study_ID") %>% # only keeps rows that have corresponding reference
      mutate(NAb_fold = EstimatedMean.y - EstimatedMean.x) %>%
      mutate(SEM = sqrt((SD.x^2 / Num_ppl.x) + (SD.y^2 / Num_ppl.y))) %>%
      select(Study_ID, Group1.x, Variant.x, NAb_fold, SEM, Num_ppl.x)%>%
      arrange(Group1.x)

    
     
    
### 4. Meta-analysis ###
    
    
  # Pooled SD for censored data
  # !!! currently only for vaccination data - can update for convalescent variant responses #
    
    ## Centre individual points at the reported mean
    centred_data <- NAb_censor_data %>%
      
      # fold of convalescent for individual points
      filter(Group1 != 'ConvWT') %>%
      inner_join((filter(NAb_censor_data, Group1 == 'ConvWT' & Variant == 'WT')), by = "Study_ID") %>%
      mutate(fold_conv = Titre_log.y - Titre_log.x)%>%
      
      # centre individual points to normalised reported mean
      inner_join(tbl_uncensored_raw, by = c("Study_ID", "Group1.x", "Variant.x", "Day_measure.x")) %>%
      mutate(centred_data = fold_conv - NAb_fold) %>%
    
      # centre normalised LOD to normalised reported mean
      mutate(centred_LOD = LOD.x - Titre_log.x - NAb_fold) %>%
      select("Study_ID", "Group1.x", "Day_measure.x", "Variant.x", "centred_data","centred_LOD")
    
    
    ## Calculate pooled SD (function from Khoury & Cromer Github)
    
    PooledSDModelFit <- nlm(function(p) {Likelihood(p, centred_data$centred_LOD, centred_data$centred_data)},
                            c(mean(centred_data$centred_data), sd(centred_data$centred_data)))    
    SD_pooled <- PooledSDModelFit$estimate[2]
    
    
    
    # Analysis 1: Sinopharm 2 doses protection against WT up to 80 days
    
    SP_SP_wt <- function (tbl) {
      SP_SP_wt <- tbl %>%
        filter(Group1.x %in%  c("SP_SP21","SP_SP40") & 
                 Variant.x == "WT" & Day_order == "ST")
      
      SP_SP_wt_MA <- metamean(n = SP_SP_wt$Num_ppl.x,
                              mean = SP_SP_wt$NAb_fold,
                              sd = SP_SP_wt$SEM,
                              random = TRUE,
                              studlab = SP_SP_wt$Group1.x,
                              title = "Sinopharm 2 dose (WT)")
      
      return(SP_SP_wt_MA)
    }
    
    
    # forest.meta(SP_SP_wt_MA, sortvar = TE)
    # funnel(SP_SP_wt_MA)
    
    
    
    # Analysis 2: Sinopharm 2 doses protection against WT more than 80 days
    
    SP_SP_wt_LT <- function (tbl) {
      SP_SP_wt_LT <- tbl %>%
        filter(Group1.x %in%  c("SP_SP21","SP_SP40") & 
                 Variant.x == "WT" & Day_order == "LT")
      
      SP_SP_wt_LT_MA <- metamean(n = SP_SP_wt_LT$Num_ppl.x,
                                 mean = SP_SP_wt_LT$NAb_fold,
                                 sd = SP_SP_wt_LT$SEM,
                                 random = TRUE,
                                 studlab = SP_SP_wt_LT$Group1.x,
                                 title = "Sinopharm 2 dose (WT)")
      
      return(SP_SP_wt_LT_MA)
    }
    
    
    
    # Analysis 3: Astrazeneca two dose wild type, up to 80 days
    
    AZ_AZ_wt <- function (tbl) {
      AZ_AZ_wt <- tbl %>%
        filter(Group1.x %in%  c("AZ_AZ28","AZ_AZ56","AZ_AZ61") & 
                 Variant.x == "WT" & Day_order == "ST")
      
      AZ_AZ_wt_MA <- metamean(n = AZ_AZ_wt$Num_ppl.x,
                              mean = AZ_AZ_wt$NAb_fold,
                              sd = AZ_AZ_wt$SEM,
                              random = TRUE,
                              studlab = AZ_AZ_wt$Group1.x,
                              title = "Astrazeneca 2 dose (WT)")
      
      return(AZ_AZ_wt_MA)
    }
    
    
    # Analysis 4: Astrazeneca two dose wild type, more than 80 days
    
    AZ_AZ_wt_LT <- function (tbl) {
      AZ_AZ_wt_LT <- tbl %>%
        filter(Group1.x %in%  c("AZ_AZ28","AZ_AZ56","AZ_AZ61") & 
                 Variant.x == "WT" & Day_order == "LT")
      
      AZ_AZ_wt_LT_MA <- metamean(n = AZ_AZ_wt_LT$Num_ppl.x,
                                 mean = AZ_AZ_wt_LT$NAb_fold,
                                 sd = AZ_AZ_wt_LT$SEM,
                                 random = TRUE,
                                 studlab = AZ_AZ_wt_LT$Group1.x,
                                 title = "Astrazeneca 2 dose (WT)")
      
      return(AZ_AZ_wt_LT_MA)
    }
    
    
    # Analysis 5: Pfizer two dose wild type, up to 80 days
    
    Pf_Pf_wt <- function (tbl) {
      Pf_Pf_wt <- tbl %>%
        filter(Group1.x %in%  c("Pf_Pf","Pf_Pf21","PZ_PZ21") & 
                 Variant.x == "WT" & Day_order == "ST")
      
      Pf_Pf_wt_MA <- metamean(n = Pf_Pf_wt$Num_ppl.x,
                              mean = Pf_Pf_wt$NAb_fold,
                              sd = Pf_Pf_wt$SEM,
                              random = TRUE,
                              studlab = Pf_Pf_wt$Group1.x,
                              title = "Pfizer 2 dose (WT)")
      
      return(Pf_Pf_wt_MA)
    }
    
    # Analysis 6: Pfizer two dose wild type, more than 80 days
    
    Pf_Pf_LT_wt <- function (tbl) {
      Pf_Pf_LT_wt <- tbl %>%
        filter(Group1.x %in%  c("Pf_Pf","Pf_Pf21","PZ_PZ21") & 
                 Variant.x == "WT" & Day_order == "LT")
      
      Pf_Pf_LT_wt_MA <- metamean(n = Pf_Pf_LT_wt$Num_ppl.x,
                                 mean = Pf_Pf_LT_wt$NAb_fold,
                                 sd = Pf_Pf_LT_wt$SEM,
                                 random = TRUE,
                                 studlab = Pf_Pf_LT_wt$Group1.x,
                                 title = "Pfizer 2 dose (WT)")
      
      return(Pf_Pf_LT_wt_MA)
    }
    
    
    # Analysis 7: Pfizer three dose wild type
    
    Pf_Pf_Pf_wt <- function (tbl) {
      Pf_Pf_Pf_wt <- tbl %>%
        filter(Group1.x %in%  c("Pf_Pf21_Pf210","PZ_PZ21_PZ301") & 
                 Variant.x == "WT")
      
      Pf_Pf_Pf_wt_MA <- metamean(n = Pf_Pf_Pf_wt$Num_ppl.x,
                                 mean = Pf_Pf_Pf_wt$NAb_fold,
                                 sd = Pf_Pf_Pf_wt$SEM,
                                 random = TRUE,
                                 studlab = Pf_Pf_Pf_wt$Group1.x,
                                 title = "Pfizer 3 dose (WT)")
      
      return(Pf_Pf_Pf_wt_MA)
    }

    
    
    # Analysis 8: Convalescent and Pfizer two dose, up to 80 days
    
    Conv_Pf_Pf_wt <- function (tbl) {
      Conv_Pf_Pf_wt <- tbl %>%
        filter(Group1.x %in%  c("Conv_Pf_Pf","Conv_Pf_Pf21") & 
                 Variant.x == "WT" & Day_order == "ST")
      
      Conv_Pf_Pf_wt_MA <- metamean(n = Conv_Pf_Pf_wt$Num_ppl.x,
                                 mean = Conv_Pf_Pf_wt$NAb_fold,
                                 sd = Conv_Pf_Pf_wt$SEM,
                                 random = TRUE,
                                 studlab = Conv_Pf_Pf_wt$Group1.x,
                                 title = "Convalescent and Pfizer 2 dose (WT)")
      
      return(Conv_Pf_Pf_wt_MA)
    }
    
    
    # Analysis 9: Convalescent and Pfizer two dose, more than 80 days
    
    Conv_Pf_Pf_wt_LT <- function (tbl) {
      Conv_Pf_Pf_wt_LT <- tbl %>%
        filter(Group1.x %in%  c("Conv_Pf_Pf","Conv_Pf_Pf21") & 
                 Variant.x == "WT" & Day_order == "LT")
      
      Conv_Pf_Pf_LT_wt_MA <- metamean(n = Conv_Pf_Pf_wt_LT$Num_ppl.x,
                                 mean = Conv_Pf_Pf_wt_LT$NAb_fold,
                                 sd = Conv_Pf_Pf_wt_LT$SEM,
                                 random = TRUE,
                                 studlab = Conv_Pf_Pf_wt_LT$Group1.x,
                                 title = "Convalescent and Pfizer 2 dose (WT)")
      
      return(Conv_Pf_Pf_wt_LT_MA)
    }
    
    

### 5. Convert NAb titres to vaccine efficacy ###
    # (covariate matrix from Deborah, email on 11 Oct)
    
    
    n_vx <- 246 # number of samples from Khoury & Cromer paper
    mean_k_n50 <- c(1.1306607, -0.6966127)
    cov_k_n50 <- matrix(c(0.03106460, 0.010755914, 0.010755914, 0.005727749), nrow=2, ncol=2)
    
    
  # function to transform NAb data to effectiveness DETERMINISTIC
    vx_efficacy <- function(k,n,n50) {
      vx_eff <- 1/(1 + exp(-k*(n - n50)))
      return(vx_eff)
    }
    
    
  # function to transform NAb data to effectiveness PROBABILISTIC
    
    runs <- 1000 # number of Monte Carlo runs
    
    k_n50_samples <- mvtnorm::rmvnorm(runs, mean = mean_k_n50, sigma = cov_k_n50)
    colnames(k_n50_samples) <- c("k","IC50")
    
    
    efficacy_psa <- function(meta, censored) {
      
      # meta is meta object from metamean function
      # censored is 1 (censored) or 0 (not censored)
      
      n_mean <- unlist(meta['TE.random'], use.names = FALSE)
      
      if (censored == 1) {sd = SD_pooled} # or pooled_sd(centred_data, ...)
      else {sd = unlist(meta['seTE.random'], use.names = FALSE)}
      ### CHECK: using SEM from meta-analysis in data that is not censored
      
      n_random <- rnorm(runs, mean = n_mean, sd = sd)
      
      monte_carlo <- mapply(vx_efficacy, k = k_n50_samples[,1], n= n_random, n50 = k_n50_samples[,2])
      
      return(monte_carlo)
      
    }
    

### 6. Final output
  # (using uncensored data, short-term only, removing points with >50% points below LOD)
    
    # Sinopharm 2 doses protection against WT up to 80 days
    SP_SP_wt_LOD <- SP_SP_wt(tbl_uncensored_LOD)
    SP_SP_wt_LOD_psa <- efficacy_psa(SP_SP_wt_LOD,0)
    
    # Astrazeneca two dose wild type, up to 80 days
    AZ_AZ_wt_LOD <- AZ_AZ_wt(tbl_uncensored_LOD)
    AZ_AZ_wt_LOD_psa <- efficacy_psa(AZ_AZ_wt_LOD,0)
    
    # Pfizer two dose wild type, up to 80 days
    Pf_Pf_wt_LOD <- Pf_Pf_wt(tbl_uncensored_LOD)
    Pf_Pf_wt_LOD_psa <- efficacy_psa(Pf_Pf_wt_LOD,0)
    
    # Sinovac 2 doses wild type, up to 80 days
    SV_SV_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "PZ35" &
                                                                     tbl_uncensored_LOD$Variant.x == "WT" &
                                                                     tbl_uncensored_LOD$Group1.x == "SV_SV28"],
                          sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "PZ35" &
                                                             tbl_uncensored_LOD$Variant.x == "WT" &
                                                             tbl_uncensored_LOD$Group1.x == "SV_SV28"])
    SV_SV_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = SV_SV_wt_LOD, n50 = k_n50_samples[,2])
    
    # AZ first dose, Pfizer second dose
    AZ_PZ_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "PZ619" &
                                                                     tbl_uncensored_LOD$Variant.x == "WT" &
                                                                     tbl_uncensored_LOD$Group1.x == "AZ_PZ"],
                          sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "PZ619" &
                                                             tbl_uncensored_LOD$Variant.x == "WT" &
                                                             tbl_uncensored_LOD$Group1.x == "AZ_PZ"])
    AZ_PZ_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = AZ_PZ_wt_LOD, n50 = k_n50_samples[,2])
    
    # AZ first dose, Sinopharm second dose
    AZ_SP_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "SP4" &
                                                                     tbl_uncensored_LOD$Variant.x == "WT" &
                                                                     tbl_uncensored_LOD$Group1.x == "AZ_SP60"],
                          sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "SP4" &
                                                             tbl_uncensored_LOD$Variant.x == "WT" &
                                                             tbl_uncensored_LOD$Group1.x == "AZ_SP60"])
    AZ_SP_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = AZ_SP_wt_LOD, n50 = k_n50_samples[,2])
    
    # Sinopharm first dose, Astrazeneca second dose
    SP_AZ_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "SP4" &
                                                                     tbl_uncensored_LOD$Variant.x == "WT" &
                                                                     tbl_uncensored_LOD$Group1.x == "SP_AZ33"],
                          sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "SP4" &
                                                             tbl_uncensored_LOD$Variant.x == "WT" &
                                                             tbl_uncensored_LOD$Group1.x == "SP_AZ33"])
    SP_AZ_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = SP_AZ_wt_LOD, n50 = k_n50_samples[,2])
    
    
    # CONVALESCENT: two doses Pfizer
    Conv_Pf_Pf_wt_LOD <-  Conv_Pf_Pf_wt(tbl_uncensored_LOD)
    Conv_Pf_Pf_wt_LOD_psa <- efficacy_psa( Conv_Pf_Pf_wt_LOD,0)
    
    
    # BOOSTER: 2 doses Astrazeneca, booster Pfizer at 6 months
    AZ_AZ_Pf_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "AZ2" &
                                                                     tbl_uncensored_LOD$Variant.x == "WT" &
                                                                     tbl_uncensored_LOD$Group1.x == "AZ_AZ56_PZ175"],
                             sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "AZ2" &
                                                                tbl_uncensored_LOD$Variant.x == "WT" &
                                                                tbl_uncensored_LOD$Group1.x == "AZ_AZ56_PZ175"])
    AZ_AZ_Pf_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = AZ_AZ_Pf_wt_LOD, n50 = k_n50_samples[,2])
    
    # BOOSTER: 2 doses Pfizer, booster Pfizer at 7 months
    Pf_Pf_Pf_wt_LOD <- Pf_Pf_Pf_wt(tbl_uncensored_LOD)
    Pf_Pf_Pf_wt_LOD_psa <- efficacy_psa(Pf_Pf_Pf_wt_LOD,0)
    
    
    # CONVALESCENT & BOOSTER: 2 doses Pfizer, Pfizer booster at 8 months
    Conv_Pf3_wt_LOD <- rnorm(runs, mean = tbl_uncensored_LOD$NAb_fold[tbl_uncensored_LOD$Study_ID == "Pf3" &
                                                                        tbl_uncensored_LOD$Variant.x == "WT" &
                                                                        tbl_uncensored_LOD$Group1.x == "Conv_Pf_Pf21_Pf240"],
                             sd = tbl_uncensored_LOD$SEM[tbl_uncensored_LOD$Study_ID == "Pf3" &
                                                                tbl_uncensored_LOD$Variant.x == "WT" &
                                                                tbl_uncensored_LOD$Group1.x == "Conv_Pf_Pf21_Pf240"])
    Conv_Pf3_wt_LOD_psa <- vx_efficacy(k = k_n50_samples[,1], n = Conv_Pf3_wt_LOD, n50 = k_n50_samples[,2])

      
    
### ASSUMPTIONS ###
    
      # all NAb data is normally distributed


### OUTSTANDING TASKS ###
      
      # meta-regression for time component with exponential decay
      # comparison of censored/non-censored results
      # include >50% below LOD as Boolean variable
