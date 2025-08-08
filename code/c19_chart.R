library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
# library(wesanderson)


# Load CSV files
base_case <- read.csv("20250219_Result for CUA.csv") %>%
  mutate(group = case_when(group == '09-May' ~ '5-9',
                           group == '14-Oct' ~ '10-14',
                           TRUE ~ as.character (group)))

cua_base <- read.csv("cua_bc_2025.csv") %>%
  mutate(ref = "base")


# 1 Base case (40% coverage) ----------------------------------------------

# Figure 1: 2025-2030 for all outcomes and scenarios, assume 40% coverage
fig1 <- base_case %>% 
  filter((cov_2024 == 0.4 | scenario == "base_case") & year > 2022 &
           compartment %in% c("cases", "severe", "death") &
           fit_end_threshold == 0.2 &
           scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case")) %>%
  select(group_index, year, compartment, incidence, scenario) %>%
  pivot_wider(names_from = group_index, values_from = incidence) %>%
  mutate(total = rowSums(across(c(4:19)))) %>%
  select(year, compartment, scenario, total) %>%
  mutate(scenario = case_when(
    scenario == "WHO" ~ "60y+ w/75y+ booster",
    TRUE ~ scenario)) %>%
  group_by(scenario)

fig1a <- fig1 %>% filter(compartment == "cases") %>%
  ggplot() +
  geom_line(aes(x=year, y=total, group=scenario, color=scenario)) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  ggtitle("Cases") +
  xlab("") + ylab("")

fig1b <- fig1 %>% filter(compartment == "severe") %>%
  ggplot() +
  geom_line(aes(x=year, y=total, group=scenario, color=scenario)) +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
  ggtitle("Hospitalisations") + 
  xlab("Year") + ylab("")

fig1c <- fig1 %>% filter(compartment == "death") %>%
  ggplot() +
  geom_line(aes(x=year, y=total, group=scenario, color=scenario)) +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
  ggtitle("Deaths") +
  xlab("") + ylab("")

# Arrange charts for figure 1
p <- plot_grid(fig1a + theme(legend.position="none"),
                    fig1b + theme(legend.position="none"),
                    fig1c + theme(legend.position="none"),
                    align = "vh", nrow = 1)

legend <- get_legend(fig1c) # get fig 1c legend

fig1_p <- plot_grid(p, legend, rel_widths = c(3, .5))

fig1_p

# Table 1: 2025-2030 summary of averted cases, hospitalisations, and deaths (relative to base case)
table1 <- fig1 %>% pivot_wider(names_from = year, values_from = total) %>%
  mutate(Total = rowSums(across(c(2:8)))) %>%
  select(compartment, scenario, Total) %>%
  pivot_wider(names_from = scenario, values_from = Total) %>%
  rename("From10y" = "10 y+", "From20y" = "20 y+", "From30y" = "30 y+", "From40y" = "40 y+", "From50y" = "50 y+",
         "From60y" = "60 y+", "From60y_booster" = "60y+ w/75y+ booster") %>%
  mutate(From10y = base_case - From10y, From20y = base_case - From20y, From30y = base_case - From30y, From40y = base_case - From40y,
         From50y = base_case - From50y, From60y = base_case - From60y, From60y_booster = base_case - From60y_booster)

write.csv(table1, "Table1_SEIR.csv", row.names = FALSE)

# Figure 2: total cases, hospitalisations, deaths averted per strategy
table1a <- table1 %>% select(-base_case) %>% pivot_longer(-compartment, names_to = "scenario", values_to = "averted")

fig2a <- table1a %>% filter(compartment == "cases") %>%
  ggplot() +
  geom_bar(stat = "identity", aes(x = scenario, y = averted, group = scenario, fill = scenario), position=position_dodge()) +
  scale_fill_viridis_d(option = "G") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  ggtitle("Cases (2025-2029)") +
  xlab("Strategy") + ylab("") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


fig2b <- table1a %>% filter(compartment == "severe") %>%
  ggplot() +
  geom_bar(stat = "identity", aes(x = scenario, y = averted, group = scenario, fill = scenario), position=position_dodge()) +
  scale_fill_viridis_d(option = "G") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  ggtitle("Hospitalisations (2025-2029)") +
  xlab("Strategy") + ylab("") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


fig2c <- table1a %>% filter(compartment == "death") %>%
  ggplot() +
  geom_bar(stat = "identity", aes(x = scenario, y = averted, group = scenario, fill = scenario), position=position_dodge()) +
  scale_fill_viridis_d(option = "G") +
  scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
  ggtitle("Deaths (2025-2029)") +
  xlab("Strategy") + ylab("") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Arrange charts for figure 1
p2 <- plot_grid(fig2a + theme(legend.position="none"),
               fig2b + theme(legend.position="none"),
               fig2c + theme(legend.position="none"),
               align = "vh", nrow = 1)

legend <- get_legend(fig2c) # get fig 1c legend

fig2_p <- plot_grid(p2, legend, rel_widths = c(3, .5))




# 2 Age-standardised incidence --------------------------------------------

# Figure 2a: 2025-2030 cases for all scenarios, separated by age, assume 40% coverage
fig2 <- base_case %>% 
  filter((cov_2025 == 0.4 | scenario == "base_case") & year > 2023 & compartment == "cases" &
           fit_end_threshold == 0.2 &
           scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case")) %>%
  select(group, year, prop, scenario) %>%
  mutate(scenario = case_when(
    # scenario == "5 y+" ~ "05 y+",
    scenario == "WHO" ~ "60y+ w/75y+ booster",
    TRUE ~ scenario)) %>%
  group_by(scenario)

# chart for each age group
plot_list = list()
for(i in unique(fig2$group)) {
  p_temp <- fig2 %>% filter(group == i) %>%
    ggplot() +
    geom_line(aes(x = year, y = prop, group = scenario, colour = scenario)) +
    ylab("") +
    ggtitle(i) +
    ylim(0,6) +
    theme(legend.position="none")
  
  plot_list[[i]] = p_temp
}

grid.arrange(grobs = plot_list)


# Figure 2c: 2025-2030 deaths for all scenarios, separated by age, assume 40% coverage
fig2b <- base_case %>% 
  filter((cov_2025 == 0.4 | scenario == "base_case") & year > 2023 & 
           scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case") &
           fit_end_threshold == 0.2 & compartment == "death") %>%
  select(group, year, prop, scenario) %>%
  mutate(scenario = case_when(
    scenario == "5 y+" ~ "05 y+",
    scenario == "WHO" ~ "60y+ w/75y+ booster",
    TRUE ~ scenario)) %>%
  group_by(scenario)

# chart for each age group

plot_list2 = list()
for(i in unique(fig2b$group)) {
  p_temp <- fig2b %>% filter(group == i) %>%
    ggplot() +
    geom_line(aes(x = year, y = prop, group = scenario, colour = scenario)) +
    ylab("") +
    ggtitle(i) +
    #ylim(0,0.05) +
    theme(legend.position="none")
  
  plot_list2[[i]] = p_temp
}

grid.arrange(grobs = plot_list2)


# Figure 2b: 2025-2030 hospitalisations for all scenarios, separated by age, assume 40% coverage
fig2c <- base_case %>% 
  filter((cov_2025 == 0.4 | scenario == "base_case") & year > 2023 & 
           scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case") &
           fit_end_threshold == 0.2 & compartment == "severe") %>%
  select(group, year, prop, scenario) %>%
  mutate(scenario = case_when(
    scenario == "5 y+" ~ "05 y+",
    scenario == "WHO" ~ "60y+ w/75y+ booster",
    TRUE ~ scenario)) %>%
  group_by(scenario)

# chart for each age group

plot_list3 = list()
for(i in unique(fig2c$group)) {
  p_temp <- fig2c %>% filter(group == i) %>%
    ggplot() +
    geom_line(aes(x = year, y = prop, group = scenario, colour = scenario)) +
    ylab("") +
    ggtitle(i) +
    ylim(0,0.15) +
    theme(legend.position="none")
  
  plot_list3[[i]] = p_temp
}

grid.arrange(grobs = plot_list3)




# 3 Cost-utility analysis --------------------------------------------------

base_case_line <- data.frame(scenario = "base", deltacost_total = 0, delta_qalygain = 0)

cua_chart <- cua_base %>% 
  filter(coverage == 0.4 & scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case")) %>%
  select(scenario, deltacost_total, delta_qalygain) %>%
  mutate(scenario = recode(scenario, "WHO" = "60y+ w/75y+ booster")) 

cua_chart <- rbind(cua_chart,base_case_line)
  

ggplot(cua_chart, aes(x = delta_qalygain, y = deltacost_total, colour = scenario, group = 1)) + 
  geom_point(size = 2) + 
  theme_minimal() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_y_continuous("Cost (USD)", labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(-12500000000, 500000000)) +
  scale_x_continuous("QALY", labels = scales::unit_format(unit = ", 000", scale = 1e-3), limits = c(-50000,650000),
                     breaks = c(200000, 400000, 600000)) 


# 4 Uncertainty analysis of CUA --------------------------------------------------


cua_critical_low <- read.csv("cua_lowercritical_2025.csv")%>%
  mutate(ref = "critical_low")

cua_severe_low <- read.csv("cua_lowersevere_2025.csv")%>%
  mutate(ref = "severe_low")

cua_critical_high <- read.csv("cua_uppercritical_2025.csv")%>%
  mutate(ref = "critical_high")

cua_severe_high <- read.csv("cua_uppersevere_2025.csv")%>%
  mutate(ref = "severe_high")

cua_uncertainty <- rbind(cua_base, cua_critical_low, cua_severe_low, cua_critical_high, cua_severe_high) %>%
  filter(coverage == 0.4 & scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case")) %>%
  select(ref, scenario, deltacost_total, delta_qalygain) %>%
  mutate(scenario = recode(scenario, "WHO" = "60y+ w/75y+ booster")) 

### HOSPITAL COSTS ###
ggplot(cua_uncertainty, aes(x = delta_qalygain, y = deltacost_total, colour = scenario, shape = ref, group = 1)) + 
  geom_point() + 
  theme_minimal() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_y_continuous("Cost (USD)", labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(-13000000000, 500000000)) +
  scale_x_continuous("QALY", labels = scales::unit_format(unit = ", 000", scale = 1e-3), limits = c(-50000,650000),
                     breaks = c(200000, 400000, 600000)) 

### COVERAGE ###

cua_coverage <- cua_base %>%
  filter(coverage %in% c(0.2, 0.4, 0.6) & 
           scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO", "base_case")) %>%
  select(coverage, scenario, deltacost_total, delta_qalygain) %>%
  mutate(scenario = recode(scenario, "WHO" = "60y+ w/75y+ booster")) %>%
  mutate(coverage = case_when(
    coverage == 0.2 ~ "20%",
    coverage == 0.4 ~ "40%",
    .default = "60%"
  ))

cua_coverage$scenario = factor(cua_coverage$scenario)

ggplot(cua_coverage, aes(x = delta_qalygain, y = deltacost_total, colour = coverage, shape = scenario, group = 1)) + 
  geom_point(size = 3) + 
  theme_minimal() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  scale_y_continuous("Cost (USD)", labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(-13000000000, 500000000)) +
  scale_x_continuous("QALY", labels = scales::unit_format(unit = ", 000", scale = 1e-3), limits = c(-50000,650000),
                     breaks = c(200000, 400000, 600000)) +
  scale_shape_manual(values = c(10:18)) +
  scale_color_manual(values=wes_palette(n=3, name="Darjeeling1"))


# 5 BIA --------------------------------------------------

bia_base <- read.csv("bia2025.csv")

bia <- bia_base %>% 
  filter(coverage == 0.4 & scenario %in% c("10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "WHO")) %>%
  select(scenario, cost_vaccine, diffcosttx) %>%
  mutate(scenario = recode(scenario, "WHO" = "60y+ w/75y+ booster")) %>%
  rename("Vaccine" = "cost_vaccine", "Hospital" = "diffcosttx") %>%
  pivot_longer(cols = c(Vaccine, Hospital), names_to = "Budget", values_to = "cost")


ggplot(data = bia) +
  geom_bar(stat = "identity", aes(x = scenario, y = cost, group = Budget, fill = Budget), position=position_dodge()) +
  scale_fill_manual(values = c("midnightblue", "aquamarine3")) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  xlab("Scenario") + ylab("Budget (USD)") +
  theme_minimal()



# EXTRA Uncertainty analysis of SEIR --------------------------------------------------

# Figure 3: 2025-2030 for all outcomes with uncertainty ranges, restricted scenarios for comparison
fig3_unc <- base_case %>% 
  filter((cov_2025 %in% c(0.2, 0.4, 0.5)) & year > 2023 & compartment %in% c("cases", "severe", "death") &
           scenario %in% c("0 y+", "10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "base_case")) %>%
  select(group_index, year, compartment, incidence, scenario, cov_2025) %>%
  pivot_wider(names_from = group_index, values_from = incidence) %>%
  mutate(total = rowSums(across(c(4:19)))) %>%
  select(year, cov_2025, compartment, scenario, total) %>%
  # filter(scenario %in% c("0 y+", "20 y+", "40 y+", "60 y+", "base_case")) %>%
  pivot_wider(names_from = cov_2025, values_from = total) %>%
  rename(low_cov = "0.2", high_cov = "0.5", base = "0.4")

fig3_lowve <- low_ve %>% 
  filter((cov_2025 == 0.4) & year > 2023 & compartment %in% c("cases", "severe", "death") &
           scenario %in% c("0 y+", "10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "base_case")) %>%
  select(group_index, year, compartment, incidence, scenario, cov_2025) %>%
  pivot_wider(names_from = group_index, values_from = incidence) %>%
  mutate(total = rowSums(across(c(4:19)))) %>%
  select(year, cov_2025, compartment, scenario, total)# %>%
  # filter(scenario %in% c("0 y+", "20 y+", "40 y+", "60 y+", "base_case")) 


fig3_highve <- high_ve %>% 
  filter((cov_2025 == 0.4) & year > 2023 & compartment %in% c("cases", "severe", "death") &
           scenario %in% c("0 y+", "10 y+", "20 y+", "30 y+", "40 y+", "50 y+", "60 y+", "base_case")) %>%
  select(group_index, year, compartment, incidence, scenario, cov_2025) %>%
  pivot_wider(names_from = group_index, values_from = incidence) %>%
  mutate(total = rowSums(across(c(4:19)))) %>%
  select(year, cov_2025, compartment, scenario, total) #%>%
  # filter(scenario %in% c("0 y+", "20 y+", "40 y+", "60 y+", "base_case")) 

# main dataframe
fig3 <- fig3_unc %>%
  mutate(low_ve = fig3_lowve$total, high_ve = fig3_highve$total) %>%
  group_by(scenario) %>% 
  filter(compartment == "severe" & year == 2025) %>%
  pivot_longer(c(base, low_cov, high_cov, low_ve, high_ve))

# charts
fig3 %>% ggplot() +
  geom_bar(stat = "identity", aes(x = name, y = value, group = scenario, fill = scenario), position=position_dodge()) +
  scale_fill_viridis_d(option = "G") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  ggtitle("2025") +
  xlab("Uncertainty analysis") + ylab("Hospitalisations") +
  theme_minimal()


