#########################################
### visualization for LRC simulations ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# sampling random individuals through time

sampled_ind <- read.table("data/LRC_dominance/summary_data/LRC_dominance_sampled_ind.tsv", 
                          header = TRUE) %>%
  mutate(generation = generation/1000)

with(
  data = sampled_ind,
  sample_through_time(gen = generation,
                      value = Param_value,
                      param = param,
                      alpha = alpha,
                      color = c("orange", "darkblue"),
                      ylimit = c(-1.5, 1.8))
)

selected_seed <- unique(sampled_ind$seed)[1]

pheno <- ggplot(data = sampled_ind %>% filter(param != "dominance" & seed == selected_seed)) +
  geom_point(aes(x=generation, y=Param_value, color = param), alpha = 0.5) +
  scale_color_manual(values = c("orange", "darkblue"))+
  geom_abline(slope = 0, intercept = 0.5, color = "orange", linetype = "dashed", linewidth = 1.5) +
  geom_abline(slope = 0, intercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1.5) +
  labs(x = expression("Generation [x" * 10^3 * "]"), y = "Value") +
  theme_light() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = c(0.15, 0.9),
        strip.text = element_text(size = 14),
        title = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 1
        ))
pheno
  
sampled_hapl <- read.table("data/LRC_dominance/summary_data/LRC_dominance_sampled_hapl.tsv", 
                           header = TRUE)

sampled_hapl_long <- sampled_hapl %>% 
  pivot_wider(names_from = param, values_from = Param_value) %>%
  mutate(generation = generation/1000)

haplo <- ggplot(data = sampled_hapl_long %>% filter(seed == selected_seed)) +
  geom_point(aes(x = slope, y = intercept, color = generation), alpha = 0.5, size = 2.5) +
  annotate("point", x = 0, y = 0.5, shape = 3, size = 5, stroke = 1.5, color = "red") +
  #facet_grid(seed ~ .) +
  scale_color_viridis_c() +
  #scale_color_gradientn(colors = c("darkblue","indianred3", "gold")) +
  #scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-2, 2, 0.5)) +
  #scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-2, 2, 0.5)) +
  labs(x = "slope" , y = "intercept", color = expression("Generation [x" * 10^3 * "]")) +
  theme_light() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = c(0.2, 0.8),
        strip.text = element_text(size = 14),
        title = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 1
        ))
haplo

sampled_hapl_long_dominance <- sampled_hapl_long %>% 
  filter(seed == selected_seed & generation > 275) %>%
  pivot_longer(cols = c(intercept, slope), names_to = "param", values_to = "param_value")

dominance <- ggplot(data = sampled_hapl_long %>% 
         filter(seed == selected_seed & generation > 275) %>%
         pivot_longer(cols = c(intercept, slope), names_to = "param", values_to = "param_value")) +
  scale_color_manual(values = c("orange", "darkblue")) +
  geom_point(aes(x = dominance, y = param_value, color = param)) +
  labs(x = "Dominance coefficient", y = "Value") +
  theme_light() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = c(0.7, 0.8),
        strip.text = element_text(size = 14),
        title = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 1
        ))

(pheno | haplo)

dominance
