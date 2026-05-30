#################################################
### data processing for LRC model simulations ###

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_data.R")

# ================================================

read_wd <- "data/gynodioecy_LRC/raw_data/"
file_names <- get_file_names(read_wd)

extract_and_saves(files = file_names,
                  nind = 20,
                  nhapl = 20,
                  save.name = "gynodioecy_LRC_",
                  read.wd = read_wd,
                  write.wd = "data/gynodioecy_LRC/summary_data")
