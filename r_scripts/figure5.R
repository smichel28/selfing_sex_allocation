# figure 5 script

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

basepar <- par(no.readonly = TRUE)


sampled_ind <- read.table("data/LRC_dominance/summary_data/LRC_dominance_sampled_ind.tsv", 
                          header = TRUE) %>% mutate(generation = generation/1000)

selected_seed <- unique(sampled_ind$seed)[1]

df <- sampled_ind %>% filter(seed == selected_seed)

# ---- figure 5A ----

png("figures/figure5A.png", height = 2000, width = 1900, res = 300)
par(mfrow = c(2,1), mar = c(3,3,1,1), oma = c(3,3,0,0), cex.axis = 1)

colors.B <- c("intercept" = "cornflowerblue", "slope" = "orange2")
shapes.B <- c("intercept" = 1, "slope" = 5)

with(data = df,
     
     for (p in c("intercept", "slope")) {     
        condition <- param == p
        plot(generation[condition],
            Param_value[condition],
            col = colors.B[param[condition]],
            pch = shapes.B[param[condition]],
            xlab = "",
            ylab = "",
            lwd = 2.5,
            cex = 0.8,
            ylim = c(-1, 1.8))
        text(35, 1.5, ifelse(p == "intercept", expression(italic(h) * ", intercept"), 
                              expression(italic(b)* ", slope")), cex = 1.5)
        abline(h = ifelse(p == "intercept", 0.5, 0), 
              col = ifelse(p == "intercept", "royalblue4", "darkorange3"),
              lwd = 2, 
              lty = 2)
     }
     
)

mtext(side = 1, text = expression("Generation [x" * 10^3 * "]"), outer = TRUE, cex = 1.8, line = 1.5)
mtext(side = 2, text = "Phenotypic value", outer = TRUE, cex = 1.8)
#mtext(expression("\u03B1 = " * 0.55 * ", " * italic(m) * " = " * 0.45 * ", " * delta * " = " *
#                   0.8), outer = TRUE, side = 3, cex = 1.7)
dev.off()

# ---- figure 5B ----

sampled_hapl <- read.table("data/LRC_dominance/summary_data/LRC_dominance_sampled_hapl.tsv", 
                           header = TRUE)

sampled_hapl_long <- sampled_hapl %>% 
  pivot_wider(names_from = param, values_from = Param_value) %>%
  mutate(generation = generation/1000) %>% 
  filter(seed == selected_seed)

my_theme <- theme(axis.title = element_text(size = 22),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  legend.title = element_text(size = 15),
                  strip.text = element_text(size = 14),
                  title = element_text(size = 16),
                  legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
                  legend.margin = margin(6, 6, 6, 6),
                  panel.border = element_rect(color = "black",fill = NA, linewidth = 1),
                  axis.ticks = element_line(linewidth = 0.5, color = "black"))

png("figures/figure5B.png", height = 2000, width = 2000, res = 300)
ggplot(data = sampled_hapl_long) +
  geom_point(aes(x=slope, y=intercept, color = generation), size = 2.5) +
  #annotate("text", x = -0.53, y = 1.8, label = expression("\u03B1 = " * 0.55 * ", " * italic(m) * " = " * 0.45 * ", " * delta * " = " *
  #                                                         0.8), size = 8) +
  annotate("point", x = 0, y = 0.5, shape = 3, size = 6, col = "red", stroke=2) +
  scale_color_viridis_c()+
  labs(x = expression(italic(b) * ", slope"),
       y = expression(italic(h)* ", intercept"),
       color = expression("Generation [x" * 10^3 * "]")) +
  theme_test() +
  my_theme +
  theme(legend.position = c(0.25,0.7))
dev.off()


# ---- figure 5C ----

png("figures/figure5C.png", height = 2000, width = 2000, res = 300)
par(cex.lab = 1.8, cex.axis = 1.5)
par(mfrow = c(1,1), mar = c(5, 5, 4, 6) + 0.1, oma = c(0,0,0,0))

df1 <- sampled_hapl_long %>%
  mutate(dominance = dominance/max(dominance)) %>%
  filter(generation>=290)

plot(df1$dominance,
     df1$intercept,
     xlab = "Scaled dominance coeff.",
     ylab = "Allelic value",
     xlim = c(-0.1,1.1),
     ylim = c(-1,2), 
     col = "cornflowerblue",
     pch = 1,
     lwd = 2.5)
points(df1$dominance,
     df1$slope, 
     col = "orange2",
     pch = 5,
     lwd = 2.5)
legend("topright", legend = c(expression(italic("h") * ", intercept"), 
                              expression(italic("b") * ", slope")), 
       col = c("cornflowerblue", "orange2"), pch = c(1,5), bg = "white",
       pt.lwd = 2.5, bty = "n", cex = 1.5)
dev.off()
