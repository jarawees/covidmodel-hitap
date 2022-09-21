# if(!require(pacman)) install.packages("pacman")
# library(pacman)
# p_load(tidyverse)

# require(tidyverse)



#### Simple example ####
params <- gen_country_basics("Thailand",
                             processes = burden_processes_az) %>% 
  update_vac_char(., 
                  ve_i = ve_az$ve_i_o[1],
                  )

params$pop[[1]]$ur <- rep(0, 16)
res <- cm_simulate(params)

res$dynamics |>
  filter(!compartment %in% c("cases", "cases_reported", "foi","foiv_l", "foiv_m", "subclinical")) |> 
  # filter(compartment %in% c("S","E", "Ia","Ip","Is","R")) |> 
  group_by(t, compartment) |> summarise(value = sum(value)) |> 
  # group_by(t) |> 
  # summarise(value = sum(value)) |> 
  # ggplot(aes(x = t, y = value)) + geom_line()
  ggplot(aes(x = t, y = value, group = compartment, color = compartment, fill = compartment)) +
  # geom_line() +
  geom_bar(position = "stack", stat = "identity")

unique(res$dynamics$compartment)
