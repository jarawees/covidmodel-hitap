dose2 %>% 
  dplyr::select(date, ends_with("daily")) %>% 
  pivot_longer(ends_with("daily")) %>% 
  mutate(name = factor(name,
                       labels = c("AstraZeneca",
                                  "Moderna",
                                  "Pfizer",
                                  "Sinopharm",
                                  "Sinovac"))) %>% 
  ggplot(., aes(x = date, y = value, color = name, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Date",
       y = "Number of individuals\ncompeting their primary course\nby vaccine product",
       color = "",
       fill = "") +
  theme_hitap +
  scale_color_lancet() +
  scale_fill_lancet() +
  scale_y_continuous(labels = function(x) sprintf("%g", x)) -> p1

dose3 %>% 
  dplyr::select(date, ends_with("daily")) %>% 
  pivot_longer(ends_with("daily")) %>% 
  mutate(name = factor(name,
                       labels = c("AstraZeneca",
                                  "Moderna",
                                  "Pfizer",
                                  "Sinopharm",
                                  "Sinovac"))) %>% 
  ggplot(., aes(x = date, y = value, color = name, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Date",
       y = "Number of individuals\nreceiving their booster dose\nby vaccine product",
       color = "",
       fill = "") +
  theme_hitap +
  scale_color_lancet() +
  scale_fill_lancet() +
  scale_y_continuous(labels = function(x) sprintf("%g", x)) -> p2

plot_grid(get_legend(p1),
          p1 + theme(legend.position = "none"), 
          p2 + theme(legend.position = "none"), 
          ncol = 1, align = "hv",
          rel_heights = c(1,10,10)) -> p

ggsave("figs/fig2_v0.png",
       plot = p,
       width = 10, height = 15)
