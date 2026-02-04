##############################################################################
### data handling for LRC model without plastic sex allocation simulations ###

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_data.R")

# ================================================

read_wd <- "data/LRC_no_plasticity/raw_data/"
file_names <- get_file_names(read_wd)

extract_and_saves(files = file_names,
                  nind = 20,
                  nhapl = 20,
                  save.name = "LRC_no_plasticity",
                  read.wd = read_wd,
                  write.wd = "data/LRC_no_plasticity/summary_data")
