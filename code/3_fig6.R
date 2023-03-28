panel %>% 
  rownames_to_column() %>% 
  filter(prioritisation == "OA and A") %>% 
  pull(rowname) %>% 
  as.numeric -> to_pull

lapply(to_pull, function(x){
  res_all[[x]] %>% 
    summarise(cases = sum(cases),
              hospitalisation = sum(severe_all) + sum(critical_all),
              death = sum(death_all))
}) %>% 
  bind_rows() %>% 
  bind_cols(panel %>% 
              rownames_to_column() %>% 
              filter(prioritisation == "OA and A")) %>% 
  pivot_longer(cols = c("cases", "hospitalisation", "death")) %>% 
  dplyr::select(boosting_level, f, name, value) %>% 
  pivot_wider(names_from = f,
              values_from = value) %>% 
  mutate(r = `2`/`1` -1,
         name = factor(name,
                       levels = c("cases", "hospitalisation", "death"),
                       labels = c("Cases", "Hospitalisations", "Deaths"))) %>% 
  ggplot(., aes(x = boosting_level,
                y = r,
                color = name)) +
  geom_point() +
  geom_line(aes(group = name)) +
  theme_hitap +
  labs(x = "Annual boosting level",
       y = "Percentage health outcome increased\ndue to reduced programme frequency",
       color = "") +
  ggsci::scale_color_lancet() -> p

ggsave("figs/fig6_v1.png", plot = p,
       width = 10, height = 10)
