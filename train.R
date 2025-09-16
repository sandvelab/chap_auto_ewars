# left side is the names used in the code, right side is the internal names in CHAP
# Cases = number of cases
# E = population
# month = month
# ID_year = year
# ID_spat = location
# rainsum = rainfall
# meantemperature = mean_temperature

options(warn=1)

train_chap <- function(train_fn, model_fn){
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  train_fn <- args[1]
  model_fn <- args[2]
  
  train_chap(train_fn, model_fn)
}


