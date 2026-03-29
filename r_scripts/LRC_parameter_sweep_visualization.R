#########################################
### visualization for LRC simulations ###

library(dplyr)
library(tidyr)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# evolution of mean phenotype through time

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

par(mfrow = c(4,3))
with(
  data = mean_data_0.05,
  mean_through_time(gen = generation,
                    value = means,
                    param = param,
                    alpha = alpha,
                    delta = delta,
                    migration = mig,
                    sim = seed,
                    ylimit = c(-0.5, 1)
  )
)

#.....................
# Boxplots ESS

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

for (m in unique(ess$migration)) {
  df <- ess[ess$migration==m,]
  with(data = df,
       boxplot_pheno(alpha, 
                     mean, 
                     param, 
                     delta, 
                     delta.threshold = 0.6,
                     color = c("orange", "darkblue"), 
                     ylimit = c(-0.75, 1),
                     legend.position = "topleft",
                     legend.size = 1.6))
  mtext(paste0("m = ", m), outer = TRUE, cex = 1.5, font = 2, line = 1.2)
}


#.....................
# sampling random individuals through time

sampled_ind_0.01 <- read.table("data/LRC/summary_data/LRC_mig_0.01_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.05 <- read.table("data/LRC/summary_data/LRC_mig_0.05_sampled_ind.tsv", 
                               header = TRUE)
sampled_ind_0.1 <- read.table("data/LRC/summary_data/LRC_mig_0.1_sampled_ind.tsv", 
                              header = TRUE)
sampled_ind_0.5 <- read.table("data/LRC/summary_data/LRC_mig_0.5_sampled_ind.tsv", 
                              header = TRUE)
sampled_ind_0.9 <- read.table("data/LRC/summary_data/LRC_mig_0.9_sampled_ind.tsv", 
                              header = TRUE)
sampled_ind_0.95 <- read.table("data/LRC/summary_data/LRC_mig_0.95_sampled_ind.tsv", 
                               header = TRUE)

sampled_ind <- rbind(sampled_ind_0.01, 
                     sampled_ind_0.05, 
                     sampled_ind_0.1, 
                     sampled_ind_0.5,
                     sampled_ind_0.9, 
                     sampled_ind_0.95)

for (m in c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95)) {
  
  for (d in c(0.6, 0.85)) {
    
    par(mfrow = c(2, 2), oma = c(0, 0, 4, 0))
    
    sampled_ind_high_delta <- sampled_ind %>%
      filter(delta == d & mig == m)
    
    sel_seeds <- c()
    for (a in unique(sampled_ind_high_delta$alpha)) {
      seeds <- unique(sampled_ind_high_delta$seed[sampled_ind_high_delta$alpha == a])
      sel_seeds <- c(sel_seeds, seeds[1])
    }
    
    sampled_ind_high_delta <- sampled_ind_high_delta %>%
      filter(seed %in% sel_seeds)
    
    if (m == 0.5) yl<- c(-1,1.2)
    
    with(
      data = sampled_ind_high_delta,
      sample_through_time(gen = generation,
                          value = Param_value,
                          param = param,
                          alpha = alpha,
                          color = c("orange", "darkblue"),
                          ylimit = c(-1, 1.2))
    )
    
    mtext(paste0("m = ", m, ", delta = ", d), outer = TRUE, cex = 1.5, font = 2)
    
  }
}

#.....................
# heatmaps

screen.sex.alloc <- function(r, slope, intercept) {
  return(pmin(pmax(r*slope+intercept,0),1))
}

r <- seq(0, 1, by = 0.01)

deltas <- unique(ess$delta)
alphas <- unique(ess$alpha)
migs <- unique(ess$migration)

colors.means <- colorRampPalette(c("darkblue", "beige", "orange2"))(100)
colors.vars <- colorRampPalette(c("white", "darkolivegreen2", "darkgreen"))(100)
colors.slopes <- colorRampPalette(c("cornflowerblue", "white", "indianred2"))(100)
colors.intercepts <- colorRampPalette(c("white", "indianred2"))(100)

breaks.means <- seq(0, 1, length.out = 101)
breaks.slopes <- seq(-max(abs(ess$mean[ess$param=="slope"])), 
                     max(abs(ess$mean[ess$param=="slope"])), 
                     length.out = 101)
breaks.intercepts <- seq(min(ess$mean[ess$param=="intercept"]),
                         max(ess$mean[ess$param=="intercept"]),
                         length.out = 101)
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
  
  p.means <- pheatmap::pheatmap(means, 
                                main = "mean",
                                color = colors.means, 
                                breaks = breaks.means,
                                cluster_cols = F, cluster_rows = F,
                                na_col = "grey40")
  p.vars <- pheatmap::pheatmap(vars, 
                               main = "variance",
                               color = colors.vars,
                               cluster_cols = F, cluster_rows = F,
                               na_col = "grey40")
  p.slopes <- pheatmap::pheatmap(slopes, 
                                 main = "slope",
                                 color = colors.slopes,
                                 breaks = breaks.slopes,
                                 cluster_cols = F, cluster_rows = F,
                                 na_col = "grey40")
  p.intercepts <- pheatmap::pheatmap(intercepts, 
                                     main = "intercept",
                                     color = colors.intercepts,
                                     breaks = breaks.intercepts,
                                     cluster_cols = F, cluster_rows = F,
                                     na_col = "grey40")
  
  gridExtra::grid.arrange(
    p.means$gtable,
    p.vars$gtable,
    p.slopes$gtable,
    p.intercepts$gtable,
    ncol = 2,
    top = paste0("delta = ", d)
  )
  
}

