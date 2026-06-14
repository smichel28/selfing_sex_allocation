# figure 3 script

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

basepar <- par(no.readonly = TRUE)


# ---- fig 3A ----
sampled_ind <- read.table("data/simple_model/summary_data/simple_model_sampled_ind.tsv", 
                          header = TRUE)

sampled_ind_high_delta <- sampled_ind %>%
  filter(delta > 0.5 & alpha %in% c(0.0, 0.75)) %>%
  mutate(generation=generation/1000)

sel_seeds <- c()
for (a in unique(sampled_ind_high_delta$alpha)) {
  seeds <- unique(sampled_ind_high_delta$seed[sampled_ind_high_delta$alpha == a])
  sel_seeds <- c(sel_seeds, seeds[4])
}

sampled_ind_high_delta <- sampled_ind_high_delta %>%
  filter(seed %in% sel_seeds)

png("figures/figure3A.png", height = 1700, width = 1900, res = 300)

par(mfrow = c(2,2), mar = c(3,3,1,1), oma = c(3,3,3,0), cex.axis = 1.3)

colors.B <- c("intercept" = "cornflowerblue", "slope" = "orange2")
shapes.B <- c("intercept" = 1, "slope" = 5)

with(data = sampled_ind_high_delta,
     
     for (p in c("intercept", "slope")) {     
       for (a in unique(alpha)) {
         condition <- alpha == a & param == p
         plot(generation[condition],
              Param_value[condition],
              col = colors.B[param[condition]],
              pch = shapes.B[param[condition]],
              xlab = "",
              ylab = "",
              lwd = 2.5,
              cex = 0.8,
              ylim = c(-0.55, 1.2))
         text(70, -0.45, ifelse(p == "intercept", expression(italic(z)[0] * ", intercept"), 
                                expression(italic(b)* ", slope")), cex = 1.5)
         abline(h = ifelse(p == "intercept", 0.5, 0), 
                col = ifelse(p == "intercept", "royalblue4", "darkorange3"),
                lwd = 2, 
                lty = 2)
       }
     }
     
)

mtext(side = 1, text = expression("Generation [x" * 10^3 * "]"), outer = TRUE, cex = 1.8, line = 1.5)
mtext(side = 2, text = "Phenotypic value", outer = TRUE, cex = 1.8)
mtext("\u03B1 = 0", outer = TRUE, side = 3, at = 0.28, cex = 1.7)
mtext("\u03B1 = 0.75", outer = TRUE, side = 3, at = 0.78, cex = 1.7)
mtext("A", outer = TRUE, side = 3, line = 0.5, adj = 0, cex = 2, font = 2)
dev.off()

# ---- fig 2B ----

sampled_hapl <- read.table("data/simple_model/summary_data/simple_model_sampled_hapl.tsv", 
                           header = TRUE)

sampled_hapl_high_delta <- sampled_hapl %>%
  filter(delta > 0.5 & alpha == 0.75) %>%
  mutate(alpha = paste0("\u03B1 = ", alpha)) %>%
  mutate(generation=generation/1000) %>%
  pivot_wider(names_from = param, values_from = Param_value)

sampled_hapl_high_delta <- sampled_hapl_high_delta %>%
  filter(seed %in% sel_seeds)

my_theme <- theme(axis.title = element_text(size = 22),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  strip.text = element_text(size = 14),
                  title = element_text(size = 16),
                  legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
                  legend.margin = margin(6, 6, 6, 6),
                  panel.border = element_rect(color = "black",fill = NA, linewidth = 1),
                  axis.ticks = element_line(linewidth = 0.5, color = "black"),
                  plot.title = element_text(face = "bold", size = 24))

png("figures/figure3B.png", height = 1700, width = 1700, res = 300)
ggplot(data = sampled_hapl_high_delta) +
  geom_point(aes(x=slope, y=intercept, color = generation), size = 2.5) +
  annotate("text", x = -0.5, y = 1.1, label = unique(sampled_hapl_high_delta$alpha), size = 8) +
  annotate("point", x = 0, y = 0.5, shape = 3, size = 6, col = "red", stroke=2) +
  scale_color_viridis_c()+
  labs(x = expression(italic(b) * ", slope"),
       y = expression(italic(z)[0] * ", intercept"),
       color = expression("Generation [x" * 10^3 * "]"),
       title = "B") +
  theme_test() +
  my_theme +
  theme(legend.position = c(0.75,0.3))
dev.off()

