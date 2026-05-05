##################################################
### data handling for simple model simulations ###

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_data.R")

# ================================================

read_wd <- "data/simple_model/raw_data/"
file_names <- get_file_names(read_wd)

extract_and_saves(files = file_names,
                  nind = 20,
                  nhapl = 20,
                  save.name = "simple_model",
                  read.wd = read_wd,
                  write.wd = "data/simple_model/summary_data")
