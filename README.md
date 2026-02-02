# selfing_sex_allocation/
Project investigating the effect of competing self-fertilization and local resource competition on sex allocation in hermaphrodites.
Uses individual-based simulations to follow the evolution of sex allocation under different situations.

## simulation_scripts/

  - **slim_models/**
    - *simple_model.slim*: base model which implements competing selfing and variable resource budget
    - *constant_selfing.slim*: similar to simple_model.slim but selfing rate is fixed (i.e. is independent of resource budget and sex allocation)
    - *LRC.slim*: implements, upon simple_model.slim, local ressource competition through limited seed dispersal (also includes a faster algorithm but does not avoid incidental selfing, which will happen with probability 1/POP_SIZE)
    - *LRC_no_plasticity.slim*: same as LRC.slim but no resource budget variation and no plasticity in sex allocation, used to check if LRC.slim provide expected results in simpler scenario

  - **simulation_launcher/**
    contains python and bash scripts to launch simulations as subprocess in parallel
 
## r_scripts/

  - **..._data.R**: scripts for processing raw outputs for each simulation model
  - **..._visualization.R**: script for visualizing model's outputs
  - **functions/**
    - *function_data.R*: contains function to collect and summarize results from simulations
    - *function_plot.R*: contains functions to create plots
  
## data/

Contains raw and processed data for each simulation model