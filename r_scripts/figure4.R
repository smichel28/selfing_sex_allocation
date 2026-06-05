# figure 4 script

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

basepar <- par(no.readonly = TRUE)

mean_data_0.05 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.05_means_ind.tsv", 
                             header = TRUE)
mean_data_0.15 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.15_means_ind.tsv", 
                             header = TRUE)
mean_data_0.25 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.25_means_ind.tsv", 
                             header = TRUE)
mean_data_0.35 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.35_means_ind.tsv", 
                             header = TRUE)
mean_data_0.45 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.45_means_ind.tsv", 
                             header = TRUE)
mean_data_0.55 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.55_means_ind.tsv", 
                             header = TRUE)
mean_data_0.65 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.65_means_ind.tsv", 
                             header = TRUE)
mean_data_0.75 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.75_means_ind.tsv", 
                             header = TRUE)
mean_data_0.85 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.85_means_ind.tsv", 
                             header = TRUE)
mean_data_0.95 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.95_means_ind.tsv", 
                             header = TRUE)

ess_0.05 <- calculate_ess(mean_data_0.05, model = "LRC")
ess_0.15 <- calculate_ess(mean_data_0.15, model = "LRC")
ess_0.25 <- calculate_ess(mean_data_0.25, model = "LRC")
ess_0.35 <- calculate_ess(mean_data_0.35, model = "LRC")
ess_0.45 <- calculate_ess(mean_data_0.45, model = "LRC")
ess_0.55 <- calculate_ess(mean_data_0.55, model = "LRC")
ess_0.65 <- calculate_ess(mean_data_0.65, model = "LRC")
ess_0.75 <- calculate_ess(mean_data_0.75, model = "LRC")
ess_0.85 <- calculate_ess(mean_data_0.85, model = "LRC")
ess_0.95 <- calculate_ess(mean_data_0.95, model = "LRC")

ess <- rbind(ess_0.05, ess_0.15, ess_0.25, ess_0.35, ess_0.45, ess_0.55, ess_0.65, ess_0.75, ess_0.85, ess_0.95)

# sampling random individuals through time

sampled_ind_0.05 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.05_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.15 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.15_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.25 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.25_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.35 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.35_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.45 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.45_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.55 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.55_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.65 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.65_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.75 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.75_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.85 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.85_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.95 <- read.table("data/LRC_parameter_sweep/summary_data/LRC_mig_0.95_sampled_ind.tsv", 
                               header = TRUE)

sampled_ind <- rbind(sampled_ind_0.05, 
                     sampled_ind_0.15, 
                     sampled_ind_0.25,
                     sampled_ind_0.35,
                     sampled_ind_0.45,
                     sampled_ind_0.55,
                     sampled_ind_0.65,
                     sampled_ind_0.75,
                     sampled_ind_0.85,
                     sampled_ind_0.95)

deltas <- unique(sampled_ind$delta)
alphas <- unique(sampled_ind$alpha)
migs <- unique(sampled_ind$mig)

#these replicates corresponds to a situation where there is no branching but lots
# of variation in the final phenotype distribution. It'll be used to determine 
# the threshold at which we consider that there is branching.
rep <- unique(sampled_ind$seed[sampled_ind$alpha==0.05 &
                                 sampled_ind$delta==0.6 & 
                                 sampled_ind$mig==0.85])

# visualise these simulations
par(mfrow=c(1,3))
for (r in rep) {
  
  d <- sampled_ind[sampled_ind$seed==r,]
  with(
    data = d,
    sample_through_time(gen = generation,
                        value = Param_value,
                        param = param,
                        alpha = alpha,
                        color = c("orange", "darkblue"),
                        ylimit = c(-1, 1.2))
  )
  
}

# calculate the maximal difference in phenotypes at the end of the simulation
# which will serve as a threshold to automatically determine if there was branching 
# in other simulations.
min.val <- c()
max.val <- c()
for (r in rep) {
  val <- sampled_ind$Param_value[sampled_ind$param=="intercept" &
                                   sampled_ind$alpha==0.05 &
                                   sampled_ind$delta==0.6 & 
                                   sampled_ind$mig==0.85 &
                                   sampled_ind$generation>80000 &
                                   sampled_ind$seed==r]
  
  min.val <- c(min.val, min(val))
  max.val <- c(max.val, max(val))
}
max.diff <- mean(max.val-min.val)

# the following code uses the the aforementionned threshold and goes through all 
# simulations and combinations of parameters to see if they branched.
is.branching.list <- list()
for (d in deltas) {
  
  is.branching <- matrix(FALSE, nrow = length(migs), ncol = length(alphas))
  colnames(is.branching) <- alphas
  rownames(is.branching) <- migs
  
  for (a in alphas) {
    for (m in migs) {
      reps <- unique(sampled_ind$seed[sampled_ind$alpha==a &
                                       sampled_ind$delta==d & 
                                       sampled_ind$mig==m])
      min.val <- c()
      max.val <- c()
      for (r in reps) {
        val <- sampled_ind$Param_value[sampled_ind$param=="intercept" &
                                         sampled_ind$alpha==a &
                                         sampled_ind$delta==d & 
                                         sampled_ind$mig==m &
                                         sampled_ind$generation>80000 &
                                         sampled_ind$seed==r]
        min.val <- c(min.val, min(val))
        max.val <- c(max.val, max(val))
        if (mean(max.val-min.val) > (max.diff+0.1))
          is.branching[as.character(m),as.character(a)] <- TRUE
      }
    }
  }
  
  is.branching.list[[as.character(d)]] <- is.branching
  
}

screen.sex.alloc <- function(r, slope, intercept) {
  return(pmin(pmax(r*slope+intercept,0),1))
}

r <- seq(0, 1, by = 0.01)

deltas <- unique(ess$delta)
alphas <- unique(ess$alpha)
migs <- unique(ess$migration)

# calculates for each combination of parameters, across all replicates, the 
# ess values of the phenotypes slope and intercept as well as the mean sex
# allocations in these simulations
# saves everything in some matrices
heatmaps <- list()
for (d in deltas) {
  
  means <- matrix(NA, nrow = length(migs), ncol = length(alphas))
  colnames(means) <- alphas
  rownames(means) <- migs
  
  vars <- slopes <- intercepts <- means
  
  for (a in alphas) {
    for (m in migs) {
      
      conditions <- ess$delta == d & ess$alpha == a & ess$migration == m
      
      mean.slope <- mean(ess$mean[ess$param == 'slope' & conditions])
      mean.intercept <- mean(ess$mean[ess$param == 'intercept' & conditions])
      
      z <- screen.sex.alloc(r, mean.slope, mean.intercept)
      
      m.z <- mean(z)
      var.z <- mean(z^2) - mean(z)^2
      
      means[as.character(m), as.character(a)] <- m.z
      vars[as.character(m), as.character(a)] <- var.z
      slopes[as.character(m), as.character(a)] <- mean.slope
      intercepts[as.character(m), as.character(a)] <- mean.intercept
      
    }
  }
  
  means[is.branching.list[[as.character(d)]]] <- NA
  vars[is.branching.list[[as.character(d)]]] <- NA
  slopes[is.branching.list[[as.character(d)]]] <- NA
  intercepts[is.branching.list[[as.character(d)]]] <- NA
  
  heatmaps[[as.character(d)]] <- list("means" = means,
                                      "vars" = vars,
                                      "slopes" = slopes,
                                      "intercepts" = intercepts)
}

# convert these matrices into long formats for visualisation
longs <- list()
for (d in names(heatmaps)) {
  for (mat in names(heatmaps[[d]])) {
    
    long <- as.data.frame(as.table(heatmaps[[d]][[mat]])) %>% 
      mutate(delta = d, param = mat)
    
    longs[[paste0(d,"_",mat)]] <- long
    
  }
}

final_long <- bind_rows(longs) 
colnames(final_long) <- c("m", "alpha", "ess", "delta", "param")
final_long <- final_long %>% 
  mutate(delta = paste0("\u03B4 = ", delta)) %>%
  mutate(label = ifelse(is.na(ess), "XY", ""))

my_theme <- theme(panel.border = element_rect(color = "black",fill = NA, linewidth = 1),
                  axis.ticks = element_line(linewidth = 0.5, color = "black"))

# actual plots
means <- ggplot(data = final_long %>% filter(param == "means"), aes(x = alpha, y = m, fill = ess)) + 
  geom_tile() +
  #geom_text(aes(label = label), size = 2)+
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "orange2", midpoint = 0.5) +
  facet_grid(delta ~ .) +
  labs(x = "\u03B1", y = "m", fill = "Mean SA") +
  theme_test() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 70, hjust = 1)) +
  my_theme
means
intercepts <- ggplot(data = final_long %>% filter(param == "intercepts"), aes(x = alpha, y = m, fill = ess)) + 
  geom_tile() +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  #scale_fill_gradient2(low = "chartreuse4", mid = "gold2", high = "darkmagenta", midpoint = 0.5) +
  facet_grid(delta ~ .) +
  #labs(x = "\u03B1", y = "m", fill = "\u2113") +
  labs(x = "\u03B1", y = "", fill = expression(italic(h))) +
  theme_test() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 70, hjust = 1)) +
  my_theme
intercepts
slopes <- ggplot(data = final_long %>% filter(param == "slopes"), aes(x = alpha, y = m, fill = ess)) + 
  geom_tile() +
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "firebrick1", midpoint = 0) +
  facet_grid(delta ~ .) +
  labs(x = "\u03B1", y = "", fill = expression(italic(b))) +
  theme_test() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 70, hjust = 1)) +
  my_theme
slopes

png("figures/figure4.png", width = 2000, height = 2400, res = 300)
(means | slopes | intercepts)
dev.off()
