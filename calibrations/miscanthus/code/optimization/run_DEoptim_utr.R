library(BioCro)
library(UTRMiscanthusBML)
library(DEoptim)
# Adapted from Yufeng He: Dec 12, 2022 by Ximin on Sep 24, 2025

fn.observed = "../../data/biomass_observation/il_observation_without_index.csv"
#Yufeng: the last DOY the obs has reduced STEM.
#since the model does not have STEM senescenece, it has a low weight on the last DOY
# Ximin: the UTR model includes senescence, so the weights are the same
observed <- read.csv(fn.observed)
observed$wt <- rep(0.2, 5)
print(observed)

#climate file for running the model
fn.climate <- "../../data/weather_data_NASA_POWER/BioCro_input_NASA/site_2_2002_2018_LowerTransmittance.csv"
# growing season weather from IL 
weather_all <- read.csv(fn.climate)
partial_gro_list = list()

# Modules required for miscanthus
source("../../data/utr_parameters/miscanthus_utr_modules.R")
source("../../data/utr_parameters/miscanthus_utr_initial_values.R")
source("../../data/utr_parameters/miscanthus_utr_params.R")

parameters_to_optimize <- c("Leaf_utilization_rate_constant", "Stem_utilization_rate_constant", "Root_utilization_rate_constant", "Rhizome_utilization_rate_constant",
                            "Leaf_utilization_km",            "Stem_utilization_km",            "Stem_utilization_km",            "Rhizome_utilization_km",
                            "Rhizome_storage_to_substrate_rate_max",     "Rhizome_storage_to_substrate_a",     "Rhizome_storage_to_substrate_b",
                            "Rhizome_substrate_to_storage_rate_max",     "Rhizome_substrate_to_storage_a",     "Rhizome_substrate_to_storage_a")

# to be updated
lower_bound_parameters <- c()
upper_bound_parameters <- c()

miscanthus_giganteus_initial_state$Rhizome = 24.1  #t/ha,https://encyclopedia.pub/entry/31968  
#obs data years 2006-2008 from doi: 10.1111/j.1757-1707.2011.01153.x
#plant year: 2002. We start simulation from the plant year for a continuous run
#This will be consistent with our multi-site validation
years = 2006:2008 
for (i in 1:length(years)){
  year_i   = years[i] 
  growing_season_weather = weather_all[weather_all$year == year_i, ]
  
  growing_season <- growing_season_weather[with(growing_season_weather,doy >=106 & doy <=350),] 
  
  # partial_gro_function <- partial_gro_solver(initial_state =  miscanthus_giganteus_initial_state,
  #                                            parameters = miscanthus_giganteus_logistic_parameters,
  #                                            varying_parameters = growing_season,
  #                                            steady_state_module_names = miscanthus_giganteus_ss_logistic_modules,
  #                                            derivative_module_names = miscanthus_giganteus_deriv_logistic_modules,
  #                                            arg_names = parameters_to_optimize,
  #                                            solver=solver,verbose = FALSE)
  partial_gro_list[[i]] = partial_gro_function
}

source("il_objfunlogistic.R")
cost_func <- function(x){
  il_objfunlogistic(x,partial_gro_list,observed)
}
# maximum number of iterations
max.iter <- 500

set.seed(1234)
# Call DEoptim function to run optimization
parVars <- c('il_objfunlogistic','partial_gro_list','observed')

cl <- makeCluster(12)
clusterExport(cl, parVars,envir=environment())
optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
                      control=list(VTR=10,itermax=max.iter,parallelType=1,packages=c('BioCroMis'),parVar=parVars,cl=cl))

#optim_result<-DEoptim(fn=cost_func, lower=lower_bound_parameters, upper = upper_bound_parameters, 
#                      control=list(VTR=10,itermax=max.iter,parallelType=0,packages=c('BioCroMis')))

opt_result <- data.frame(para=optim_result$par,MSE=optim_result$value)

saveRDS(opt_result,'opt_result_DEoptim_3year_run_r1.rds')
