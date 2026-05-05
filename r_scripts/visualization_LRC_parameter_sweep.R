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

rep <- unique(sampled_ind$seed[sampled_ind$alpha==0.05 &
                                 sampled_ind$delta==0.6 & 
                                 sampled_ind$mig==0.85])

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
  
  print( min(d$Param_value[d$param=="intercept"]))
  
}

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


is.branching.list <- list()
for (d in deltas) {
  
  is.branching <- matrix(FALSE, nrow = length(migs), ncol = length(alphas))
  colnames(is.branching) <- alphas
  rownames(is.branching) <- migs
  
  for (a in alphas) {
    for (m in migs) {
      rep <- unique(sampled_ind$seed[sampled_ind$alpha==a &
                                       sampled_ind$delta==d & 
                                       sampled_ind$mig==m])
      min.val <- c()
      max.val <- c()
      for (r in rep) {
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

for (replicat in 1:3) {
  for (m in migs) {
    
    png(paste0("~/GitHub/selfing_sex_allocation/figures/sampled_ind_mig_", m, "_rep_", replicat,".png"), width = 3000, height = 3000, res = 200)
    
    par(mfrow = c(4, 5), oma = c(0, 4, 4, 0))
    
    for (d in deltas[deltas>0.5]) {
      
      sampled_ind_high_delta <- sampled_ind %>%
        filter(delta == d & mig == m)
      
      sel_seeds <- c()
      for (a in unique(sampled_ind_high_delta$alpha)) {
        seeds <- unique(sampled_ind_high_delta$seed[sampled_ind_high_delta$alpha == a])
        sel_seeds <- c(sel_seeds, seeds[replicat])
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
      
      pos <- 0.25
      if (d == 0.6) pos <- 0.75
      
      mtext(paste0("delta = ", d), outer = TRUE, cex = 1.5, font = 2, side = 2, at = pos)
      
    }
    
    mtext(paste0("m = ", m), outer = TRUE, cex = 1.5, font = 2)
    dev.off()
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
  
  means[is.branching.list[[as.character(d)]]] <- NA
  vars[is.branching.list[[as.character(d)]]] <- NA
  slopes[is.branching.list[[as.character(d)]]] <- NA
  intercepts[is.branching.list[[as.character(d)]]] <- NA
  
  p.means <- pheatmap::pheatmap(means, 
                                main = "mean",
                                color = colors.means, 
                                breaks = breaks.means,
                                cluster_cols = F, cluster_rows = F,
                                na_col = "grey30",
                                border_color = NA,
                                fontsize = 8)
  p.vars <- pheatmap::pheatmap(vars, 
                               main = "variance",
                               color = colors.vars,
                               cluster_cols = F, cluster_rows = F,
                               na_col = "grey30",
                               border_color = NA,
                               fontsize = 8)
  p.slopes <- pheatmap::pheatmap(slopes, 
                                 main = "slope",
                                 color = colors.slopes,
                                 breaks = breaks.slopes,
                                 cluster_cols = F, cluster_rows = F,
                                 na_col = "grey30",
                                 border_color = NA,
                                 fontsize = 8)
  p.intercepts <- pheatmap::pheatmap(intercepts, 
                                     main = "intercept",
                                     color = colors.intercepts,
                                     breaks = breaks.intercepts,
                                     cluster_cols = F, cluster_rows = F,
                                     na_col = "grey30",
                                     border_color = NA,
                                     fontsize = 8)
  
  png(paste0("~/GitHub/selfing_sex_allocation/figures/heatmap_delta_", d, ".png"), width = 1700, height = 1500, res = 200)
  gridExtra::grid.arrange(
    p.means$gtable,
    p.vars$gtable,
    p.slopes$gtable,
    p.intercepts$gtable,
    ncol = 2,
    top = grid::textGrob(paste0("\u03B4 = ", d), gp = grid::gpar(fontsize = 14, fontface = "bold")),
    left = grid::textGrob("m", rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold")),
    bottom = grid::textGrob("\u03B1", gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  dev.off()
  
}

