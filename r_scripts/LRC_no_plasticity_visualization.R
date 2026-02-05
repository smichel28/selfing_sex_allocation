##################################################
### visualization for simple model simulations ###

library(dplyr)
library(tidyr)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# evolution of mean phenotype through time

mean_data <- read.table("data/LRC_no_plasticity/summary_data/LRC_no_plasticity_means_ind.tsv", 
                        header = TRUE)

par(mfrow = c(3,3))
with(
  data = mean_data,
  mean_through_time(gen = generation,
                    value = means,
                    param = param,
                    alpha = alpha,
                    delta = delta,
                    migration = mig,
                    sim = seed)
)

ess <- calculate_ess(mean_data, model = "LRC")

for (m in unique(ess$migration)) {
  df <- ess[ess$migration==m,]
  with(data = df,
       boxplot_pheno(alpha, mean, param, delta, color = c("orange", "darkblue")))
}

#.....................
# sampling random individuals through time
# 
# sampled_ind <- read.table("data/LRC_no_plasticity/summary_data/LRC_no_plasticity_sampled_ind.tsv", 
#                           header = TRUE)
# 
# par(mfrow = c(2,2))
# with(
#   data = sampled_ind,
#   sample_through_time(gen = generation,
#                       value = Param_value,
#                       param = param,
#                       alpha = alpha,
#                       color = c("orange", "darkblue"))
# )
# 

