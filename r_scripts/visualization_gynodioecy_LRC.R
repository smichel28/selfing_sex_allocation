#########################################
### visualization for LRC simulations ###

library(dplyr)
library(tidyr)

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_plot.R")

#=================================================

#.....................
# sampling random individuals through time

sampled_ind <- read.table("data/gynodioecy_LRC/summary_data/gynodioecy_LRC__sampled_ind.tsv", 
                          header = TRUE)

with(
  data = sampled_ind,
  sample_through_time(gen = generation,
                      value = Param_value,
                      param = param,
                      alpha = alpha,
                      color = c("orange", "darkblue"),
                      ylimit = c(-1.5, 1.8))
)
