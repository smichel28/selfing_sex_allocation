library(dplyr)
library(tidyr)
library(data.table)

get_file_names <- function(wd) {
  return(paste(wd, list.files(path = wd), sep = '/'))
}

read_file_name <- function(file.name, model = 'simple') {
  # reads the file's name and return the parameters of the simulations
  
  param <- strsplit(file.name, split = '_')[[1]] # splits the file's name
  col_names <- param[seq(1, length(param), 2)] # gathers the parameters' names
  values <- param[seq(2, length(param), 2)] # gathers the values of the parameters
  
  # for the seed, split the last element of the values to remove the ".txt"
  seed <- strsplit(values[length(values)], split = '.txt')[[1]][1]
  values[length(values)] <- seed
  
  return(list(columns = col_names, parameter_values = values))
}

extract_and_saves <- function(files, 
                              every.gen = 1000, 
                              model = 'simple',
                              nind = 10,
                              nhapl = 10,
                              save_name,
                              write.wd) {
  
  means_ind <- paste(write.wd, paste0(save_name, '_means_ind.tsv'), sep = '/')
  means_hapl <-  paste(write.wd, paste0(save_name, '_means_hapl.tsv'), sep = '/')
  
  sampled_ind <-  paste(write.wd, paste0(save_name, '_sampled_ind.tsv'), sep = '/')
  sampled_hapl <-  paste(write.wd, paste0(save_name, '_sampled_hapl.tsv'), sep = '/')
  
  file.create(means_ind)
  file.create(means_hapl)
  file.create(sampled_ind)
  file.create(sampled_hapl)
  
  write_col <- TRUE
  
  for (f in files) {
    
    # gets information about the simulation parameter for this specific file
    info <- read_file_name(f)
    
    # gets data from file and keeps only desired generations
    data <- fread(f, header = TRUE) %>%
      filter(generation %% every.gen == 0)
    
    n_row <- nrow(data)
    n_param <- length(info[['columns']])
    
    parameters_simulation <- NULL
    for (i in 1:n_param) {
      
      col_name <- info[['columns']][i]
      value <- info[['parameter_values']][i]
      
      if (col_name == 'nsampled') {n_sampled <- value}
      
      if (is.null(parameters_simulation)) {
        parameters_simulation <- data.frame(col_name = rep(value, times = n_row))
      } else {
        parameters_simulation <- cbind(parameters_simulation, 
                                       data.frame(col_name = rep(value, times = n_row)))
      }
    }
    
    # creates two datasets, one contains phenotypes and the other contains allele values
    data_ind <- cbind(parameters_simulation, data[,1:(n_sampled+2)])
    print(head(data_ind))
    data_hapl <- cbind(parameters_simulation, data[,c(1,2,(n_sampled+3):(n_sampled*2+2))])
    print(head(data_hapl))
    
    # saves the mean phenotype and allele every desired generation
    write.table(
      cbind(parameters_simulation, data.frame(means = rowMeans(data_ind[,-c(1:(n_param+2))]))),
      file = means_ind,
      sep = "\t",
      row.names = FALSE,
      col.names = write_col,
      append = TRUE
    )
    write.table(
      cbind(parameters_simulation, data.frame(means = rowMeans(data_hapl[,-c(1:(n_param+2))]))),
      file = means_hapl,
      sep = "\t",
      row.names = FALSE,
      col.names = write_col,
      append = TRUE
    )

    # saves a few phenotypes or alleles every desired generation
    write.table(
      data_ind[, 1:(n_param+2+nind)] %>% 
        pivot_longer(cols = starts_with("Ind"), 
                     names_to = "Individual", 
                     values_to = "Param_value"),
      file = sampled_ind,
      sep = "\t",
      row.names = FALSE,
      col.names = write_col,
      append = TRUE
    )
    write.table(
      data_hapl[, 1:(n_param+2+nhapl)] %>% 
        pivot_longer(cols = starts_with("Hapl"), 
                     names_to = "Haplotype", 
                     values_to = "Param_value"),
      file = sampled_hapl,
      sep = "\t",
      row.names = FALSE,
      col.names = write_col,
      append = TRUE
    )
    
    write_col <- FALSE
  }
}

