# figure 2 script

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

# theme 
my_theme <- theme(axis.title = element_text(size = 18),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 9),
                  legend.title = element_text(size = 13),
                  strip.text = element_text(size = 14),
                  title = element_text(size = 16),
                  legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
                  panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
                  axis.ticks = element_line(linewidth = 0.5, color = "black"))

# panel A

mean_data <- read.table("data/simple_model/summary_data/simple_model_means_ind.tsv", 
                        header = TRUE)
ess <- calculate_ess(mean_data, model = "simple")

df.reaction.norm <- ess %>% 
  filter(delta < 0.5 & alpha %in% c(0,0.3,0.7,0.9) == FALSE) %>%
  group_by(alpha, delta, param) %>%
  summarise(mean_value = mean(mean)) %>%
  pivot_wider(names_from = param, values_from = mean_value) %>%
  crossing(R = seq(0,1,0.1)) %>% 
  mutate(z = R*slope+intercept)

A <- ggplot(data = df.reaction.norm) +
  geom_line(aes(x=R,y=z, color = factor(alpha), linetype = factor(delta)), linewidth = 1.2) +
  geom_abline(slope = 0, intercept = 0.5, linetype = "dashed", color = "grey", linewidth = 2) +
  scale_color_brewer(palette = "RdYlGn") +
  labs(x=expression(italic(R) * ", resource budget"), 
       y = expression(italic(z) * ", sex allocation"), 
       color = "\u03B1",
       linetype = expression(delta)) +
  scale_y_continuous(limits = c(0.4,1)) +
  theme_grey()+
  my_theme
A

basepar <- par(no.readonly = TRUE)

df <- ess %>% 
  filter(delta < 0.5 & alpha %in% c(0,0.3,0.7,0.9) == FALSE) %>%
  group_by(alpha, delta, param) %>%
  summarise(mean_value = mean(mean)) %>% 
  ungroup() %>%
  pivot_wider(names_from = param, values_from = mean_value)

colors.A <- hcl.colors(10, palette = "RdYlGn", rev = TRUE)
names(colors.A) <- sort(unique(df$alpha))

par(mar = c(5, 4, 4, 6))
plot(0, 0,
     type = "n",
     xlim = c(0, 1),
     ylim = c(0.5, 1),
     xlab = expression(italic(R) * ", resource budget"),
     ylab = expression(italic(z) * ", sex allocation"))
abline(h = 0.5, lty = 1, lwd = 3, col = "lightgrey")
for (a in unique(df$alpha)) {
    condition <- df$alpha == a
    abline(a = df$intercept[condition & df$delta == 0], 
           b = df$slope[condition & df$delta == 0],
           col = colors.A[as.character(a)],
           lwd = 3)
    abline(a = df$intercept[condition & df$delta == 0.4], 
           b = df$slope[condition & df$delta == 0.4],
           col = colors.A[as.character(a)],
           lty = 2,
           lwd = 2)
}
legend(1.05, 1.02, legend = sort(unique(df$alpha)), 
       col = colors.A, lty = 1, lwd = 3, bg = "white", 
       title = "\u03B1", xpd = NA, bty = "n")
legend(1.05, 0.77, legend = sort(unique(df$delta)), 
       col = "black", lty = c(1,2), lwd = c(2,3), bg = "white",
       title = expression(delta), xpd = NA, bty = "n")


# panel B

sampled_ind <- read.table("data/simple_model/summary_data/simple_model_sampled_ind.tsv", 
                          header = TRUE)

sampled_ind_high_delta <- sampled_ind %>%
  filter(delta > 0.5 & alpha %in% c(0.0, 0.75)) %>%
  #mutate(alpha = paste0("\u03B1 = ", alpha)) %>%
  mutate(generation=generation/1000)

sel_seeds <- c()
for (a in unique(sampled_ind_high_delta$alpha)) {
  seeds <- unique(sampled_ind_high_delta$seed[sampled_ind_high_delta$alpha == a])
  sel_seeds <- c(sel_seeds, seeds[4])
}

sampled_ind_high_delta <- sampled_ind_high_delta %>%
  filter(seed %in% sel_seeds)

B <- ggplot(data = sampled_ind_high_delta) +
  geom_point(aes(x=generation, y=Param_value)) +
  facet_grid(param ~ alpha, labeller = labeller(
    param = as_labeller(c(
      "slope"     = "italic(b)",
      "intercept" = "italic(h)"
    ), label_parsed))) +
  labs(x = expression("Generation [x" * 10^3 * "]"),
       y = "Phenotypic value") +
  theme_grey() +
  my_theme
B

par(mfrow = c(2,2), mar = c(3,3,2,2), oma = c(3,3,0,0))

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
    text(70, -0.45, ifelse(p == "intercept", expression(italic(h) * " (intercept)"), 
                         expression(italic(b)* " (slope)")),
         cex = 1.5)
  }
}

)



mtext(side = 1, text = expression("Generation [x" * 10^3 * "]"), outer = TRUE)
mtext(side = 2, text = "Phenotypic value", outer = TRUE)

# panel C

sampled_hapl <- read.table("data/simple_model/summary_data/simple_model_sampled_hapl.tsv", 
                          header = TRUE)

sampled_hapl_high_delta <- sampled_hapl %>%
  filter(delta > 0.5 & alpha == 0.75) %>%
  mutate(alpha = paste0("\u03B1 = ", alpha)) %>%
  mutate(generation=generation/1000) %>%
  pivot_wider(names_from = param, values_from = Param_value)

sampled_hapl_high_delta <- sampled_hapl_high_delta %>%
  filter(seed %in% sel_seeds)

C <- ggplot(data = sampled_hapl_high_delta) +
  geom_point(aes(x=slope, y=intercept, color = generation)) +
  scale_color_viridis_c()+
  facet_grid(. ~ alpha) +
  labs(x = expression(italic(b) * ", slope"),
       y = expression(italic(h)* ", intercept"),
       color = expression("Generation [x" * 10^3 * "]")) +
  theme_grey() +
  my_theme +
  theme(legend.position = c(0.75,0.3))
C

png("figures/figure2.png", height = 2000, width = 2000, res = 300)
fig2 <- A + B + C +
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(1, 2, 2)) &
  theme(plot.tag = element_text(size = 18))
dev.off()

png("figures/figure2A.png", height = 1500, width = 1500, res = 300)
A
dev.off()
png("figures/figure2B.png", height = 1500, width = 1500, res = 300)
B
dev.off()
png("figures/figure2C.png", height = 1500, width = 1500, res = 300)
C
dev.off()

ggsave(filename = "figures/figure2A.png",
       plot = A,
       height = 5, width = 5,
       units = "cm",
       dpi = 300)



df <- ess %>% 
  filter(delta < 0.5 & alpha %in% c(0,0.3,0.7,0.9) == FALSE) %>%
  group_by(alpha, delta, param) %>%
  summarise(mean_value = mean(mean)) %>% 
  ungroup() %>%
  pivot_wider(names_from = param, values_from = mean_value)

colors.A <- hcl.colors(10, palette = "RdYlGn", rev = TRUE)
names(colors.A) <- sort(unique(df$alpha))

par(cex.lab = 1.8, cex.axis = 1.2)
par(mfrow = c(1,1), mar = c(5, 5, 4, 6) + 0.1, oma = c(0,0,0,0))

plot(0, 0,
     type = "n",
     xlim = c(0, 1),
     ylim = c(0.5, 1),
     xlab = expression(italic(R) * ", resource budget"),
     ylab = expression(italic(z) * ", sex allocation"))
abline(h = 0.5, lty = 1, lwd = 3, col = "lightgrey")
for (a in unique(df$alpha)) {
  condition <- df$alpha == a
  abline(a = df$intercept[condition & df$delta == 0], 
         b = df$slope[condition & df$delta == 0],
         col = colors.A[as.character(a)],
         lwd = 3)
  abline(a = df$intercept[condition & df$delta == 0.4], 
         b = df$slope[condition & df$delta == 0.4],
         col = colors.A[as.character(a)],
         lty = 2,
         lwd = 2)
}
legend(1.05, 1.02, legend = sort(unique(df$alpha)), 
       col = colors.A, lty = 1, lwd = 3, bg = "white", 
       title = "\u03B1", xpd = NA, bty = "n", title.cex = 1.3)
legend(1.05, 0.77, legend = sort(unique(df$delta)), 
       col = "black", lty = c(1,2), lwd = c(3,2), bg = "white",
       title = expression(delta), xpd = NA, bty = "n", title.cex = 1.3)




par(mfrow = c(2,2), mar = c(3,3,2,2), oma = c(3,3,0,0))

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
         text(70, -0.45, ifelse(p == "intercept", expression(italic(h) * " (intercept)"), 
                                expression(italic(b)* " (slope)")),
              cex = 1.5)
       }
     }
     
)

mtext(side = 1, text = expression("Generation [x" * 10^3 * "]"), outer = TRUE, cex = 1.8, line = 1.5)
mtext(side = 2, text = "Phenotypic value", outer = TRUE, cex = 1.8)

par(mfrow = c(1,1), mar = c(5, 5, 4, 2) + 0.1, oma = c(0,0,0,0))

palette <- viridis::viridis(100)
colors.C <- palette[as.numeric(cut(sampled_hapl_high_delta$generation, 100))]

plot(sampled_hapl_high_delta$slope, 
     sampled_hapl_high_delta$intercept,
     pch = 16, 
     col = colors.C)

legend("right",
       legend = c("0-9","10-19", "20-29", "30-39", "40-49", "50-59", "60-69"),
       fill = viridis(10),
       bty = "n",
       title = "generation", inset = 0.02)

C <- ggplot(data = sampled_hapl_high_delta) +
  geom_point(aes(x=slope, y=intercept, color = generation)) +
  scale_color_viridis_c()+
  facet_grid(. ~ alpha) +
  labs(x = expression(italic(b)),
       y = expression(italic(h)),
       color = expression("Generation [x" * 10^3 * "]")) +
  theme_test() +
  theme(legend.position = c(0.7,0.2),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 14),
        title = element_text(size = 16),
        legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
        axis.ticks = element_line(linewidth = 0.5, color = "black"))
C
