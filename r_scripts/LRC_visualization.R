##################################################
### visualization for simple model simulations ###

library(dplyr)
library(tidyr)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# evolution of mean phenotype through time

mean_data_0.01 <- read.table("data/LRC/summary_data/LRC_mig_0.01_means_ind.tsv", 
                        header = TRUE)
mean_data_0.05 <- read.table("data/LRC/summary_data/LRC_mig_0.05_means_ind.tsv", 
                             header = TRUE)
mean_data_0.1 <- read.table("data/LRC/summary_data/LRC_mig_0.1_means_ind.tsv", 
                             header = TRUE)
mean_data_0.5 <- read.table("data/LRC/summary_data/LRC_mig_0.5_means_ind.tsv", 
                             header = TRUE)

par(mfrow = c(4,3))
with(
  data = mean_data_0.01,
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

ess_0.01 <- calculate_ess(mean_data_0.01, model = "LRC")
ess_0.05 <- calculate_ess(mean_data_0.05, model = "LRC")
ess_0.1 <- calculate_ess(mean_data_0.1, model = "LRC")
ess_0.5 <- calculate_ess(mean_data_0.5, model = "LRC")

ess <- rbind(ess_0.01, ess_0.05, ess_0.1, ess_0.5)

for (m in unique(ess$migration)) {
  df <- ess[ess$migration==m,]
  with(data = df,
       boxplot_pheno(alpha, mean, param, delta, color = c("orange", "darkblue"), ylimit = c(-1, 1)))
}
