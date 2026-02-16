############################################################
### data handling for constant selfing model simulations ###

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_data.R")

# ================================================

read_wd <- "data/constant_selfing/raw_data/"
file_names <- get_file_names(read_wd)

extract_and_saves(files = file_names,
                  nind = 20,
                  nhapl = 20,
                  save.name = "constant_selfing",
                  read.wd = read_wd,
                  write.wd = "data/constant_selfing/summary_data")
