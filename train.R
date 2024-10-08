# Y = number of cases
# E = pop.var.dat
# T1 = week
# T2 = year
# S1 = district

library(INLA)
#
# # Read in command line args filenames
args = commandArgs(trailingOnly=TRUE)
if (FALSE){
  args = c('training_data.csv', 'tmp.csv')
}
data_filename = args[1] # 'training_data.csv'#  args[1]
output_model_filename = args[2] # 'tmp.csv'
source('lib.R')
df <- read.table(data_filename, sep=',', header=TRUE)
print(head(df))
#df$week = as.numeric(substr(df$time_period, 6, 8))
df = offset_years_and_months(df)
print(head(df))
selectedFormula=lagged_formula
save(selectedFormula, file=output_model_filename)
#basis_meantemperature = extra_fields(df)
#basis_rainfall = get_basis_rainfall(df)
#model = mymodel(basis_formula, df, config = TRUE)
#model2 = mymodel(lagged_formula, df, config = TRUE)

# GOF
#if(model$dic$dic < model2$dic$dic) {
#  selectedFormula = basis_formula
#}else {
#  selectedFormula = lagged_formula
#}
#save(selectedFormula, file=output_model_filename)
