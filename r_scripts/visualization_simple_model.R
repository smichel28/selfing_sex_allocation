##################################################
### visualization for simple model simulations ###

library(dplyr)
library(tidyr)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# evolution of mean phenotype through time

mean_data <- read.table("data/simple_model/summary_data/simple_model_means_ind.tsv", 
                        header = TRUE)

par(mfrow = c(3,1))
with(
  data = mean_data,
  mean_through_time(gen = generation,
                    value = means,
                    param = param,
                    alpha = alpha,
                    delta = delta,
                    sim = seed)
  )

#.....................
# sampling random individuals through time

sampled_ind <- read.table("data/simple_model/summary_data/simple_model_sampled_ind.tsv", 
                          header = TRUE)

sampled_ind_high_delta <- sampled_ind %>%
  filter(delta > 0.5)

sel_seeds <- c()
for (a in unique(sampled_ind_high_delta$alpha)) {
  seeds <- unique(sampled_ind_high_delta$seed[sampled_ind_high_delta$alpha == a])
  sel_seeds <- c(sel_seeds, seeds[2])
}

sampled_ind_high_delta <- sampled_ind_high_delta %>%
  filter(seed %in% sel_seeds)

par(mfrow = c(2,2))
with(
  data = sampled_ind_high_delta,
  sample_through_time(gen = generation,
                      value = Param_value,
                      param = param,
                      alpha = alpha,
                      color = c("orange", "darkblue"))
)

  
