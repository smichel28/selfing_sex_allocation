# figure 2 script

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

basepar <- par(no.readonly = TRUE)

# ---- panel A ----

mean_data_A <- read.table("data/simple_model/summary_data/simple_model_means_ind.tsv", 
                          header = TRUE)
ess_A <- calculate_ess(mean_data_A, model = "simple")


df <- ess_A %>% 
  filter(delta < 0.5 & alpha %in% c(0,0.3,0.7,0.9) == FALSE) %>%
  group_by(alpha, delta, param) %>%
  summarise(mean_value = mean(mean)) %>% 
  ungroup() %>%
  pivot_wider(names_from = param, values_from = mean_value)

colors.A <- hcl.colors(10, palette = "RdYlGn", rev = TRUE)
names(colors.A) <- sort(unique(df$alpha))

png("figures/figure2.png", width = 2900, height = 2200, res = 300)
par(cex.lab = 1.8, cex.axis = 1.2)
#par(mfrow = c(1,2), mar = c(4, 2.5, 1, 6) + 0.1, oma = c(0,3.5,0,0))
par(mfrow = c(1,2), mar = c(5, 2.5, 3, 2) + 0.1, oma = c(0,3.5,0,0))

plot(0, 0,
     type = "n",
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = expression(italic(R) * ", resource budget"),
     ylab = "",
     xaxt = "n",
     yaxt = "n")
axis(1, at = seq(0,1,0.25))
axis(2, at = seq(0,1,0.25))
abline(h = 0.5, lty = 1, lwd = 3, col = "grey70")
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
legend("bottomright", legend = sort(unique(df$alpha)), 
       col = colors.A, lty = 1, lwd = 3, bg = "white", 
       title = "\u03B1", title.cex = 1.3)
legend("bottomleft", legend = sort(unique(df$delta)), 
       col = "black", lty = c(1,2), lwd = c(3,2), bg = "white",
       title = expression(delta), title.cex = 1.3)
#legend(1.05, 1.02, legend = sort(unique(df$alpha)), 
#       col = colors.A, lty = 1, lwd = 3, bg = "white", 
#       title = "\u03B1", xpd = NA, bty = "n", title.cex = 1.3)
#legend(1.05, 0.77, legend = sort(unique(df$delta)), 
#       col = "black", lty = c(1,2), lwd = c(3,2), bg = "white",
#       title = expression(delta), xpd = NA, bty = "n", title.cex = 1.3)
text(0.15, 1, expression(italic(m) * " = 1"), cex = 1.5)
mtext("A", side = 3, line = 0.5, adj = 0, cex = 2, font = 2)

# ---- panel B ----
mean_data_B <- read.table("data/LRC_no_selfing/processed_data/LRC_no_selfing_means_ind.tsv", 
                          header = TRUE)

ess_B <- calculate_ess(mean_data_B, model = "LRC")

dfB <- ess_B %>% 
  group_by(migration, param) %>%
  summarise(mean_value = mean(mean)) %>% 
  ungroup() %>%
  pivot_wider(names_from = param, values_from = mean_value)

colors.B <- hcl.colors(10, palette = "RdYlGn", rev = TRUE)
names(colors.B) <- sort(unique(dfB$migration))

#par(cex.lab = 1.8, cex.axis = 1.2)
#par(mfrow = c(1,1), mar = c(5, 5, 4, 6) + 0.1, oma = c(0,0,0,0))

plot(0, 0,
     type = "n",
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = expression(italic(R) * ", resource budget"),
     ylab = "",
     xaxt = "n",
     yaxt = "n")
axis(1, at = seq(0,1,0.25))
axis(2, at = seq(0,1,0.25))
abline(h = 0.5, lty = 1, lwd = 3, col = "grey70")
for (m in unique(dfB$migration)) {
  condition <- dfB$migration == m
  abline(a = dfB$intercept[condition], 
         b = dfB$slope[condition],
         col = colors.B[as.character(m)],
         lwd = 3)
}
legend("topright", legend = sort(unique(dfB$migration)), 
       col = colors.B, lty = 1, lwd = 3, bg = "white", 
       title = expression(italic(m)), title.cex = 1.3)
#legend(1.05, 1.02, legend = sort(unique(dfB$migration)), 
#       col = colors.B, lty = 1, lwd = 3, bg = "white", 
#       title = expression(italic(m)), xpd = NA, bty = "n", title.cex = 1.3)
text(0.15, 1, "\u03B1 = 0", cex = 1.5)
mtext("B", side = 3, line = 0.5, adj = 0, cex = 2, font = 2)
#mtext(side = 1, text = expression(italic(R) * ", resource budget"), outer = TRUE, cex = 1.8, line = 1.5)
mtext(side = 2, expression(italic(z) * ", sex allocation"), outer = TRUE, cex = 1.8, line = 1.5)

dev.off()


