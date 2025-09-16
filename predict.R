# left side is the names used in the code, right side is the internal names in CHAP
# Cases = number of cases
# E = population
# month = month
# ID_year = year
# ID_spat = location
# rainsum = rainfall
# meantemperature = mean_temperature

options(warn=1)

library(INLA)
library(dlnm)
library(dplyr)

#for spatial effects
library(sf)
library(spdep)

predict_chap <- function(model_fn, hist_fn, future_fn, preds_fn, graph_fn){
  #model <- readRDS(file = model_fn) #would normally load a model here
  
  df <- read.csv(future_fn)
  df$Cases <- rep(NA, nrow(df))
  df$disease_cases <- rep(NA, nrow(df)) #so we can rowbind it with historic
  
  historic_df = read.csv(hist_fn)
  df <- rbind(historic_df, df)
  df <- mutate(df, ID_year = ID_year - min(df$ID_year) + 1)
  
  #adding a counting variable for the months like 1, ..., 12, 13, ...
  #could also do years*12 + months, but fails for weeks
  df <-group_by(df, location) |>
    mutate(month_num = row_number())
  
  basis_meantemperature <- crossbasis(df$meantemperature, lag=3,
        argvar = list(fun = "ns", knots = equalknots(df$meantemperature, 2)),
        arglag = list(fun = "ns", knots = 3/2), group = df$ID_spat)
  colnames(basis_meantemperature) = paste0("basis_meantemperature.", colnames(basis_meantemperature))

  basis_rainsum <- crossbasis(df$rainsum, lag=3,
         argvar = list(fun = "ns", knots = equalknots(df$rainsum, 2)),
         arglag = list(fun = "ns", knots = 3/2), group = df$ID_spat)
  colnames(basis_rainsum) = paste0("basis_rainsum.", colnames(basis_rainsum))
  
  df$ID_spat <- as.factor(df$ID_spat)
  df$ID_spat_num <- as.numeric(as.factor(df$ID_spat))
  
  df <- cbind(df, basis_meantemperature, basis_rainsum)
  
  # the formula without spatial smoothing and with some specified PC priors
  #lagged_formula <- Cases ~ 1 + f(ID_spat, model='iid', hyper=list(prec = list(prior = "pc.prec",
  #  param = c(1, 0.01)))) + f(month_num, model = "rw1", scale.model = T,
  #  replicate = ID_spat_num, hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
  #  f(month, model='rw1', cyclic=T, scale.model=T, hyper=list(prec = list(prior = "pc.prec",
  #  param = c(1, 0.01)))) + basis_meantemperature + basis_rainsum

  lagged_formula <- Cases ~ 1 + f(ID_spat, model='iid') + f(month_num, model = "rw1", scale.model = T,
    replicate = ID_spat_num) +
    f(month, model='rw1', cyclic=T, scale.model=T) + basis_meantemperature + basis_rainsum
 
  if(graph_fn != ""){
    df$ID_spat_num2 <- df$ID_spat_num
    
    geojson <- st_read(graph_fn)
    geojson <- st_make_valid(geojson) #some geojson files will giver errors later if this is not cleaned here
    
    df_locs <- unique(df$location)
    geojson_red <- geojson[geojson$VARNAME_1 %in% df_locs, ] #only keeps the regions present in the dataframe
    #both df_locs and geosjon_red should be in alfabetical order, believe the files are alfabetical in CHAP, might need to sort them
    #if they are both alphabetical, and take care about which field in the geeojson, could be vietnamese or english (handle in CHAP?)
    #then the regions should harmonize correctly in the formulas below
    
    nb <- poly2nb(geojson_red, queen = TRUE) #adjacent polygons are neighbors
    adjacency <- nb2mat(nb, style = "B", zero.policy = TRUE) #converts it to an adjacency matrix which is passed to the bym2
    
    #lagged_formula <- Cases ~ 1 + f(ID_spat_num, model = "bym2", graph = adjacency) +
    #  f(month, model='rw1', cyclic=T, scale.model=T) + basis_meantemperature + basis_rainsum
    
    lagged_formula <- Cases ~ 1 + f(ID_spat_num, model = "besag", graph = adjacency, scale.model = T) +
      f(ID_spat_num2, model = "iid") + f(month_num, model = "rw1", scale.model = T) +
       f(month, model='rw1', cyclic=T, scale.model=T) + basis_meantemperature + basis_rainsum
  }
  
  model <- inla(formula = lagged_formula, data = df, family = "nbinomial", offset = log(E),
                control.inla = list(strategy = 'adaptive'),
                control.compute = list(config = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = TRUE, prec.intercept = 0.1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE),
                verbose = F, safe=FALSE)

  casestopred <- df$Cases # response variable
  
  # Predict only for the cases where the response variable is missing
  idx.pred <- which(is.na(casestopred)) #this then also predicts for historic values that are NA, not ideal
  mpred <- length(idx.pred)
  s <- 1000
  y.pred <- matrix(NA, mpred, s)
  # Sample parameters of the model
  xx <- inla.posterior.sample(s, model)  # This samples parameters of the model
  xx.s <- inla.posterior.sample.eval(function(idx.pred) c(theta[1], Predictor[idx.pred]), xx, idx.pred = idx.pred) # This extracts the expected value and hyperparameters from the samples
  
  # Sample predictions
  for (s.idx in 1:s){
    xx.sample <- xx.s[, s.idx]
    y.pred[, s.idx] <- rnbinom(mpred,  mu = exp(xx.sample[-1]), size = xx.sample[1])
  }
  
  # make a dataframe where first column is the time points, second column is the location, rest is the samples
  # rest of columns should be called sample_0, sample_1, etc
  new.df = data.frame(time_period = df$time_period[idx.pred], location = df$location[idx.pred], y.pred)
  colnames(new.df) = c('time_period', 'location', paste0('sample_', 0:(s-1)))
  
  # Write new dataframe to file
  write.csv(new.df, preds_fn, row.names = FALSE)
  #saveRDS(model, file = model_fn) # to evaluate the model
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  model_fn <- args[1]
  hist_fn <- args[2]
  future_fn <- args[3]
  preds_fn <- args[4]
  graph_fn <- args[5]
  
  predict_chap(model_fn, hist_fn, future_fn, preds_fn, graph_fn)
}


# testing
#library(dplyr)
# 
#  # testing with chap data 
# model_fn <- "model"
# hist_fn <- "example_data/historic_data.csv"
# future_fn <- "example_data/future_data.csv"
# preds_fn <- "example_data/predictions.csv"

# testing with vietnam data with polygons
# model_fn <- "model"
# hist_fn <- "example_data_Viet/historic_data.csv"
# future_fn <- "example_data_Viet/future_data.csv"
# preds_fn <- "example_data_Viet/predictions.csv"
# graph_fn <- "example_data_Viet/vietnam.json"
# 

# df$Cases <- rep(NA, nrow(df))
# df$disease_cases <- rep(NA, nrow(df)) #so we can rowbind it with historic
# 
# historic_df = read.csv(hist_fn)
# df <- rbind(historic_df, df)
# df <- offset_years_and_months(df)
# df$ID_year <- df$ID_year - min(df$ID_year) + 1
#
# # some data exploration with acf and ccf from timeseries analysis
#
# un <- unique(df$location)
# df_loc1 <- df[df$location == un[3], ][c("rainfall", "disease_cases")]
# 
# plot(model$summary.random$month_num$mean[800:1000])
# 
# library(ggplot2)
# 
# df_loc1$t <- 1:nrow(df_loc1)  # Create a time index from 1 to n
# df_long <- df_loc1 %>% pivot_longer(cols = -t, names_to = "variable", values_to = "value")
# 
# ggplot(df_long, aes(x = t, y = value, color = variable)) +
#   geom_line() +
#   labs(title = "title", 
#        x = "x", y = "y") +
#   theme(legend.title = element_blank())
# 
# library(zoo)
# df_loc1$rainfall <- na.locf(df_loc1$rainfall)
# df_loc1$disease_cases <- na.locf(df_loc1$disease_cases)
# 
# ccf(df_loc1$rainfall, df_loc1$disease_cases, 
#     main = "Cross-Correlation between Rainfall and Disease Cases")
# acf(df_loc1$disease_cases)

# 
# # testing for geospatial components
# geojson <- st_read(graph_fn)
# geojson <- st_make_valid(geojson)
# 
# #might also need to alphabetically sort both geojson after location and the
# # dataframe after location!
# df <- read.csv(future_fn)
# df_locs <- unique(df$location)
# geojson$VARNAME_1
# 
# nb <- poly2nb(geojson, queen = TRUE) #adjacent polygons are neighbors
# adjacency <- nb2mat(nb, style = "B", zero.policy = TRUE)
# 
# geojson_reduced <- geojson[geojson$VARNAME_1 %in% df_locs, ]
# nb_red <- poly2nb(geojson_reduced, queen = TRUE) #adjacent polygons are neighbors
# adjacency_red <- nb2mat(nb_red, style = "B", zero.policy = TRUE)
# row.names(adjacency_red) <- NULL
# image(adjacency_red)
# 
# #plot of Vietnam for all the regions
# plot(st_geometry(geojson),
#      col = "steelblue", lwd = 0.5, asp = 1)
# 
# #plot of the regions of Vietnam in the training data
# plot(st_geometry(geojson_reduced),
#      col = "steelblue", lwd = 0.5, asp = 1)
# 
# 
# #For Brazil
# 
# graph_fn <- "brazil_polygons.geojson"
# future_fn <- "brazil.csv"
# df <- read.csv(future_fn)
# df_locs <- unique(df$location)
# geojson <- st_read(graph_fn)
# geojson <- st_make_valid(geojson)
# geojson$id
# 
# nb <- poly2nb(geojson, queen = TRUE) #adjacent polygons are neighbors
# adjacency <- nb2mat(nb, style = "B", zero.policy = TRUE)
# 
# geojson_reduced <- geojson[geojson$id %in% df_locs, ]
# nb_red <- poly2nb(geojson_reduced, queen = TRUE) #adjacent polygons are neighbors
# adjacency_red <- nb2mat(nb_red, style = "B", zero.policy = TRUE)
# row.names(adjacency_red) <- NULL
# image(adjacency_red)
# 
# #plot of Vietnam for all the regions
# plot(st_geometry(geojson),
#      col = "steelblue", lwd = 0.5, asp = 1)
# 
# #plot of the regions of Vietnam in the training data
# plot(st_geometry(geojson_reduced),
#      col = "steelblue", lwd = 0.5, asp = 1)
# 
# df = read.csv(future_fn)
# #df <- mutate(df, ID_year = ID_year - min(df$ID_year) + 1)
# 
# #adding a counting variable for the months like 1, ..., 12, 13, ...
# #could also do years*12 + months, but fails for weeks
# df <-group_by(df, location) |>
#   mutate(month_num = row_number())
# 
# basis_meantemperature <- crossbasis(df$mean_temperature, lag=3,
#                                     argvar = list(fun = "ns", knots = equalknots(df$mean_temperature, 2)),
#                                     arglag = list(fun = "ns", knots = 3/2), group = df$ID_spat)
# colnames(basis_meantemperature) = paste0("basis_meantemperature.", colnames(basis_meantemperature))
# 
# basis_rainsum <- crossbasis(df$rainsum, lag=3,
#                             argvar = list(fun = "ns", knots = equalknots(df$rainsum, 2)),
#                             arglag = list(fun = "ns", knots = 3/2), group = df$ID_spat)
# colnames(basis_rainsum) = paste0("basis_rainsum.", colnames(basis_rainsum))
# 
# df$ID_spat <- as.factor(df$ID_spat)
# df$ID_spat_num <- as.numeric(as.factor(df$ID_spat))
# 
# df <- cbind(df, basis_meantemperature, basis_rainsum)
# 
# # the formula without spatial smoothing and with some specified PC priors
# #lagged_formula <- Cases ~ 1 + f(ID_spat, model='iid', hyper=list(prec = list(prior = "pc.prec",
# #  param = c(1, 0.01)))) + f(month_num, model = "rw1", scale.model = T,
# #  replicate = ID_spat_num, hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
# #  f(month, model='rw1', cyclic=T, scale.model=T, hyper=list(prec = list(prior = "pc.prec",
# #  param = c(1, 0.01)))) + basis_meantemperature + basis_rainsum
# 
# lagged_formula <- Cases ~ 1 + f(ID_spat, model='iid') + f(month_num, model = "rw1", scale.model = T,
#                                                           replicate = ID_spat_num) +
#   f(month, model='rw1', cyclic=T, scale.model=T) + basis_meantemperature + basis_rainsum
# 
# if(graph_fn != ""){
#   #df$ID_spat_num2 <- df$ID_spat_num
#   
#   geojson <- st_read(graph_fn)
#   geojson <- st_make_valid(geojson) #some geojson files will giver errors later if this is not cleaned here
#   
#   df_locs <- unique(df$location)
#   geojson_red <- geojson[geojson$VARNAME_1 %in% df_locs, ] #only keeps the regions present in the dataframe
#   #both df_locs and geosjon_red should be in alfabetical order, believe the files are alfabetical in CHAP, might need to sort them
#   #if they are both alphabetical, and take care about which field in the geeojson, could be vietnamese or english (handle in CHAP?)
#   #then the regions should harmonize correctly in the formulas below
#   
#   nb <- poly2nb(geojson_red, queen = TRUE) #adjacent polygons are neighbors
#   adjacency <- nb2mat(nb, style = "B", zero.policy = TRUE) #converts it to an adjacency matrix which is passed to the bym2
#   
#   lagged_formula <- Cases ~ 1 + f(ID_spat_num, model = "bym2", graph = adjacency, replicate = month_num) + 
#     f(month_num, model = "rw1", scale.model = T, replicate = ID_spat_num) +
#     f(month, model='rw2', cyclic=T, scale.model=T) + basis_meantemperature + basis_rainsum
# }
# 
# model <- inla(formula = lagged_formula, data = df, family = "nbinomial", offset = log(E),
#               control.inla = list(strategy = 'adaptive'),
#               control.compute = list(config = TRUE, return.marginals = FALSE),
#               control.fixed = list(correlation.matrix = TRUE, prec.intercept = 0.1, prec = 1),
#               control.predictor = list(link = 1, compute = TRUE),
#               verbose = F, safe=FALSE)
# 
