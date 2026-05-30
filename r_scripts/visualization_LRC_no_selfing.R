#########################################
### visualization for LRC simulations ###

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

mean_data <- read.table("data/LRC_no_selfing/processed_data/LRC_no_selfing_means_ind.tsv", 
                             header = TRUE)

ess <- calculate_ess(mean_data, model = "LRC")

df.reaction.norm <- ess %>% 
  group_by(migration, param) %>%
  summarise(mean_value = mean(mean)) %>%
  pivot_wider(names_from = param, values_from = mean_value) %>%
  crossing(R = seq(0,1,0.1)) %>% 
  mutate(z = R*slope+intercept)

expected_sex_alloc <- function(alpha, delta, m) {
  return((m*(2-m)*(1-alpha*delta))/(1-alpha+m*(2-m)*(1-alpha*delta*(2-alpha))))
}

expected <- data.frame(R = 0.5, migration=seq(0.05,0.95, 0.1), exp = expected_sex_alloc(0,0,seq(0.05,0.95, 0.1)))

ggplot() +
  geom_line(data = df.reaction.norm, aes(x=R,y=z, color = factor(migration))) +
  geom_point(data = expected, aes(x=R, y=exp, colour = factor(migration)))+
  scale_color_brewer(palette = "RdYlGn") +
  labs(x="Resource budget", 
       y = "Z (sex allocation)", 
       color = "m", 
       title = "alpha = 0") +
  scale_y_continuous(limits = c(0,1)) +
  theme_light() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = c(0.9, 0.77),
        strip.text = element_text(size = 14),
        title = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 1
        ))

###
#cest cohérent de pas avoir ce qu'on attend quand on ajoute de la plasticité car
#la variation de resource disponible va influencer les quantités de chaque type 
#de gamete et donc les retours en fitness qui en découlent (en effet, si les 
#individus avec bcp de ressources sont plus mâles -> devient intéressant de produire des ovules )
###