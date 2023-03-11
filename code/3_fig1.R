p_load(ggsci)

theme_hitap <- theme_cowplot() +
  theme(legend.position = "top",
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14))

fig1_date_range <- data.frame(date = seq(ymd("2020-01-01"),
                       ymd("2022-09-30"),
                       by = "day"))

# Figure 1, date range = 2020-01-12 to 2022-09-20
epi %>% 
  mutate(used_in_fitting = if_else(txn_date <= "2021-08-26", "Included in fitting", "Excluded in fitting"),
         used_in_fitting = factor(used_in_fitting,
                                  levels = c("Included in fitting", "Excluded in fitting")),
         txn_date = ymd(txn_date)) %>% 
  right_join(fig1_date_range, by = c("txn_date" = "date")) %>% 
  ggplot(aes(x = txn_date,
             y = new_death,
             color = used_in_fitting)) +
  geom_point() +
  labs(x = "Date",
       y = "Daily Reported COVID-19 Deaths",
       color = "") +
  ggsci::scale_color_lancet(na.translate = F) +
  theme_hitap -> p1

owid_vac %>% 
  dplyr::select(date, people_fully_vaccinated_daily, total_boosters_daily) %>% 
  pivot_longer(ends_with("daily")) %>% 
  mutate(name = factor(name,
                       labels = c("Primary Course",
                                  "Booster Dose")),
         date = ymd(date)) %>% 
  right_join(fig1_date_range, by = "date") %>% 
  ggplot(.) +
  geom_point(aes(x = date, y = value, color = name))+
  labs(x = "Date", y = "Daily Completion", color = "") +
  theme_hitap +
  scale_color_manual(values = ggsci::pal_lancet()(4)[3:4], na.translate = F) +
  scale_y_continuous(labels = function(x) sprintf("%g", x)) -> p2

# c("retail", "grocery", "transit", "workplaces")
mobility %>% 
  filter(grepl("retail|grocery|transit|workplaces", mobility_type)) %>% 
  right_join(fig1_date_range, by = "date") %>% 
  mutate(mobility_type = factor(mobility_type,
                                levels = c("grocery_and_pharmacy",
                                           "retail_and_recreation",
                                           "transit_stations",
                                           "workplaces"),
                                labels = c("Grocery and pharmacy",
                                           "Retail and recreation",
                                           "Transit stations",
                                           "Workplaces"))) %>% 
  ggplot(.) +
  geom_point(aes(x = date, y = mobility_level, color = mobility_type)) +
  scale_color_manual(values = pal_lancet()(8)[5:8], na.translate = F) +
  theme_hitap +
  labs(x = "Date",
       y = "Mobility Level (% change)",
       color = "") -> p3

oxcgrt[,2:4] %>% 
  mutate(`C1M_School closing` = `C1M_School closing`*100/3) %>% 
  pivot_longer(cols = 2:3) %>% 
  rename(date = Date) %>% 
  right_join(fig1_date_range, by = "date") %>%
  mutate(name = factor(name,
                       labels = c("School closure",
                                  "Containment and health index"))) %>% 
  ggplot(., aes(x = date, y = value, color = name)) +
  geom_point() +
  scale_color_manual(values = pal_lancet()(9)[c(1,2)], na.translate = F) +
  labs(x = "Date",
       color = "",
       y = "Magnitude of implementaiton (%)") +
  theme_hitap -> p4

plot_grid(p1, p2, p3, p4, ncol = 1, align = "hv",
          labels = c("(a)","(b)","(c)","(d)"),
          label_size = 16) -> fig1

ggsave("figs/fig1_v0.png",
       width = 10, height = 20)
