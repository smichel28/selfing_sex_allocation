###############################################
### data handling for LRC model simulations ###

setwd("~/GitHub/selfing_sex_allocation/")
source("r_scripts/functions/functions_data.R")

# ================================================

read_wd <- "data/LRC/raw_data/"
file_names <- get_file_names(read_wd)

filter_files <- strsplit(file_names, split = 'mig_')
selected_files <- c()

for (f in filter_files) {
  name <- f[2]
  selected_files <- c(selected_files, as.numeric(strsplit(name, "_")[[1]][1]))
}

for (m in c(0.01, 0.05, 0.1, 0.5)) {
  
  cat(paste0("Processing files migration = ", m, "\n"))
  
  files <- file_names[selected_files == m]
  extract_and_saves(files = files,
                    nind = 20,
                    nhapl = 20,
                    save.name = paste0("LRC_mig_", m),
                    read.wd = read_wd,
                    write.wd = "data/LRC/summary_data")
}
