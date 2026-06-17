library(BioCro)
library(UTRMiscanthusBML)
library(BioCroWater)
library(ggplot2)
library(tidyr)
library(lattice)
library(DEoptim)

# Ximin on Sep 24, 2025

# parameters and weather inputs
source('miscanthus_utr_params.R')
source('miscanthus_utr_initial_values.R')
source('miscanthus_utr_modules.R')

# Set up BioCroWater
source('set_up_BioCroWater.R')
initial_values       <- set_init_values(initial_values)
parameters           <- set_parameters(parameters)
parameters$kd        <- parameters$k_diffuse
direct_modules       <- set_direct_modules(miscanthus_giganteus_direct_utr_modules) 
differential_modules <- set_differential_modules(miscanthus_giganteus_differential_utr_modules) 


# Parameter Optimization
opt_params <- c('Leaf_utilization_rate_constant', # 1
                'Stem_utilization_rate_constant', 
                'Rhizome_utilization_rate_constant', 
                'Root_utilization_rate_constant',
                'Leaf_utilization_km', 
                'Stem_utilization_km', # 6
                'Rhizome_utilization_km', 
                'Root_utilization_km', 
                'Stem_respiration_factor', 
                'Rhizome_respiration_factor', 
                'Root_respiration_factor', # 11
                'Rhizome_storage_to_substrate_rate_max', 
                'Rhizome_substrate_to_storage_rate_max',
                'storage_release_threshold', 
                'storage_to_substrate_km', 
                'substrate_to_storage_km', # 16
                'substrate_conductance_Leaf_to_Stem', 
                'substrate_conductance_Stem_to_Rhizome',
                'substrate_conductance_Rhizome_to_Root',
                'Leaf_senescence_fraction_max', 
                'Stem_senescence_fraction_max', # 21
                'Rhizome_senescence_fraction_max', 
                'Root_senescence_fraction_max',
                'Leaf_senescence_alpha', 
                'Stem_senescence_alpha',
                'Rhizome_senescence_alpha', # 26 
                'Root_senescence_alpha',
                'Leaf_senescence_beta', 
                'Stem_senescence_beta', 
                'Rhizome_senescence_beta', 
                'Root_senescence_beta', # 31
                'Leaf_senescence_reuse_factor',  
                'Stem_senescence_reuse_factor',  
                'Root_senescence_reuse_factor' 
)

test_params <- as.numeric(unlist(parameters[opt_params]))
output <- "0.445600    0.456183    0.222439    0.404927    0.006914    0.001473    0.004175    0.008541    0.016783    0.027012    0.051935    0.012635    1.767666    0.169505    0.003992    0.000147   14.517889    2.256671    3.387821    0.001486    0.000286    0.000000    0.000000    2.850732    0.907102    0.256907    0.810725    2.077753    3.113324    1.571546    5.302402    2.078172    0.443028    0.782029"
test_params <- scan(text = output)
parameters[opt_params] <- test_params
upperlim <- test_params * 3
lowerlim <- test_params * 0.1

weather_all<-read.csv(
  "../data/weather_data_NASA_POWER/BioCro_input_NASA/site_2_2002_2018_LowerTransmittance.csv")
result <- list()
partial_biocro_list <- list()
years = 2006:2008
for (i in 1:length(years)){
  year_i = years[i]
  growing_season_weather = weather_all[weather_all$year == year_i,]
  growing_season_weather =  growing_season_weather[growing_season_weather$doy>=106 &
                                                    growing_season_weather$doy<=350,]
  

  result[[i]] <- run_biocro(initial_values =  initial_values,
                       parameters = parameters,
                       drivers = growing_season_weather,
                       direct_module_names = direct_modules,
                       differential_module_names = differential_modules,
                       ode_solver = BioCro::default_ode_solvers$boost_rkck54,
                       verbose = FALSE)
  
  partial_biocro_list[[i]] <- partial_run_biocro(
    initial_values =  initial_values,
    parameters = parameters,
    drivers = growing_season_weather,
    direct_module_names = direct_modules,
    differential_module_names = differential_modules,
    arg_names = opt_params,
    ode_solver= BioCro::default_ode_solvers$boost_rkck54,
    verbose = FALSE)
}


source("plot_single_site.R")
#avg of 2006-2008
predicted = 0
for (i in 1:length(years)){
  predicted = predicted+result[[i]][,c("doy","Stem","Leaf","Rhizome","Root")]
}
predicted = predicted/length(years)
predicted = reshape2::melt(predicted,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root","Rhizome"))

observed_biomass <- read.csv("../data/biomass_observation/il_observation_without_index.csv")

names(observed_biomass)  = c("doy", "Rhizome",     "Leaf"     ,"Stem", "Root")
observed_biomass$Root[1] = 0 # First data point of root should be zero. Observation does not make difference between alive and dead roots
observed_biomass = observed_biomass[, c("doy","Stem","Leaf","Rhizome","Root")]
observed = reshape2::melt(observed_biomass,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root","Rhizome"))
names(observed)  = c("x","varname","y")
names(predicted) = c("x","varname","y")
xlabtitle = expression(paste("Day of Year"))
ylabtitle = expression(paste("Dry Biomass (Mg/ha)"))
BiomassPartitioning <- 
  Compare_Observed_and_Predicted(observeddata = observed,predicteddata=predicted  ,xlabtitle,ylabtitle)
plot(BiomassPartitioning)
ggsave(filename="Biomass.png", plot=BiomassPartitioning, width=6, height =4, units = "in", dpi = 300)


optim_cost_func <- function (params, partial_biocro_list, observed){
  cat(params); flush.console()
  predicted = list()
  predicted$stem = 0
  predicted$root = 0
  predicted$rhizome = 0
  predicted$leaf = 0
  
  time_ind = c(1,1273,2713,4201,5665) # calculated from (doy-first doy)*24+1
  
  for (i in 1:length(partial_biocro_list))
  {
    partial_biocro_function = partial_biocro_list[[i]]
    # res <- partial_biocro_function(params)
    res <- tryCatch(
      partial_biocro_function(params),
      error = function(e) { message("Worker error: ", e$message); NULL }
    )
    if (is.null(res) || dim(res)[1] < 5880) {
      # print(paste0('Result is Null?: ', is.null(res), '; Result run till timestep ', dim(res)[1]))
      # flush.console()
      return(1e5)
    }
    predicted$stem    = predicted$stem    + res$Stem[time_ind]
    predicted$root    = predicted$root    + res$Root[time_ind]
    predicted$rhizome = predicted$rhizome + res$Rhizome[time_ind]
    predicted$leaf    = predicted$leaf    + res$Leaf[time_ind]
  }
  #3-year averages
  num_years <- length(partial_biocro_list)
  predicted$stem    = predicted$stem/num_years
  predicted$root    = predicted$root/num_years
  predicted$rhizome = predicted$rhizome/num_years
  predicted$leaf    = predicted$leaf/num_years
  
  # print(observed)
  # print(predicted)
  
  # First data point of root should be zero. Observation does not make difference between alive and dead roots
  observed$root[1] = 0
  
  #be careful of the variable names of the observation
  #if they don't match, the program won't stop!
  #you may see warning message like "stack imbalance". Very annoying!
  stemE    = observed$Stem - predicted$stem 
  rootE    = observed$Root - predicted$root 
  rhizomeE = observed$Rhizome - predicted$rhizome 
  leafE    = observed$Leaf - predicted$leaf 
  
  # deleted the original weights
  E = sum(stemE^2)  + sum(rootE^2) + sum(rhizomeE^2)  + sum(leafE^2)
  # print(paste0('Leaf Error: ', leafE))
  # print(paste0('Stem Error: ', stemE))
  # print(paste0('Rhizome Error: ', rhizomeE))
  # print(paste0('Root Error: ', rootE))
  
  if(is.na(E)) {
    # print(params)
    E=1e5
  }
  # print(E)
  return(E)
}

cost_func <- function(x){
  optim_cost_func(x , partial_biocro_list, observed_biomass)
}

cost_func(test_params)

# r <- partial_biocro_list[[1]](testing_params)
# xyplot(data=r, Leaf+Stem+Root+Pod~fractional_doy)

# maximum number of iterations
max.iter <- 500

set.seed(123)
# Call DEoptim function to run optimization
parVars <- c('cost_func', 'optim_cost_func','partial_biocro_list', 'observed_biomass')

cl <- makeCluster(8)
clusterExport(cl, parVars,envir=environment())
sink(paste0('Optmization_output_utr_', Sys.Date(), '.txt'))
optim_result<-DEoptim(fn=cost_func, lower=lowerlim, upper = upperlim, 
                      control=list(VTR=10,
                                   itermax=max.iter, 
                                   parallelType=1,
                                   packages=c('BioCro', 'UTRMiscanthusBML'),
                                   parVar=parVars,cl=cl))

opt_result <- data.frame(para=optim_result$par,MSE=optim_result$value)
print(opt_result)
sink()
sink()


# examining_plots = TRUE
# if (examining_plots){
#   # plot the substrate/storage C "concentrations"
#   concentrations <- with(result[[1]], data.frame(
#     doy = doy,
#     Leaf_substrate = Leaf_substrate_carbon/Leaf_structural_carbon,
#     Stem_substrate = Stem_substrate_carbon/Stem_structural_carbon,
#     Root_substrate = Root_substrate_carbon/Root_structural_carbon,
#     Rhizome_substrate = Rhizome_substrate_carbon/Rhizome_structural_carbon,
#     Rhizome_storage = Rhizome_storage_carbon/Rhizome_structural_carbon
#   ))
#   xyplot(data=concentrations[-c(1:24), ], 
#          Leaf_substrate + Stem_substrate + Root_substrate + Rhizome_substrate + Rhizome_storage ~ doy, 
#          ylab = "Ratio to Structural C",
#          pch = 19, cex = 0.5, alpha = 0.5,
#          auto.key = TRUE)
#   
#   xyplot(data=concentrations[1000:5880,], 
#          Leaf + Stem + Root + Rhizome ~ doy, 
#          ylab = "Substrate C: Structural C",
#          auto.key = TRUE)
#   
#   xyplot(data=result[[1]], 
#          Leaf_substrate_carbon + 
#            Stem_substrate_carbon + 
#            Root_substrate_carbon + 
#            Rhizome_substrate_carbon + 
#            Rhizome_storage_carbon ~ doy, 
#          ylab = "Substrate C (mol C / m^2)",
#          ylim = c(0,1),
#          auto.key = TRUE)
#   
#   xyplot(data=result[[1]], 
#          Leaf_structural_carbon + 
#            Stem_structural_carbon + 
#            Root_structural_carbon + 
#            Rhizome_structural_carbon ~ doy, 
#          ylab = "Structural C (mol C / m^2)",
#          auto.key = TRUE)
#   
#   # plot the utilization rate compared to transport to check 
#   # 1. if there is significant imbalances
#   # 2. if the rhizome exports first and then import C
#   # 3. if Stem_utilization_rate is high enough for the rapid growth
#   # 4. if not is it limited by the stem substrate concentration or the utilization parameters?
#   # 5. Is Rhizome growing at the end of the season?
#   # 6. Do leaf and stem utilization max out, so more can be transported to Rhizome?
#   xyplot(data=result[[1]], 
#          canopy_assimilation_rate * 3.33 +
#            Leaf_utilization_rate+
#            substrate_transport_Leaf_to_Stem+
#            Stem_utilization_rate+
#            substrate_transport_Stem_to_Rhizome+
#            Rhizome_utilization_rate
#          ~doy, 
#          pch = 19, cex = 0.4, alpha = 0.2,
#          ylab = 'Transport/Utilization rate (mol C/ m^2 / hr)', 
#          auto.key = list(
#            text = c("Canopy assimilation", 
#                     "Leaf utilization", 
#                     "Leaf to stem transport", 
#                     "Stem utilization", 
#                     "Stem to rhizome transport", 
#                     "Rhizome utilization")))
#   
#   # plot the utilization rate compared to transport to check if there is significant imbalances
#   transport <- with(result[[1]], data.frame(
#     doy = doy,
#     Leaf_Stem = substrate_transport_Leaf_to_Stem,
#     Stem_Root = substrate_transport_Rhizome_to_Root,
#     Stem_Rhizome = substrate_transport_Stem_to_Rhizome
#   ))
#   xyplot(data=transport, 
#          Leaf_Stem + Stem_Root + Stem_Rhizome ~ doy, 
#          ylab = "Transport rate (mol C/ m^2 / hr)" , 
#          pch = 19, cex = 0.3, alpha = 0.5,
#          auto.key = TRUE,)
#   
#   # Compare canopy_assimilation_rate with storage release rate.
#   xyplot(data=result[[1]],
#          Rhizome_storage_to_substrate_rate+
#            canopy_assimilation_rate~
#            doy, auto.key=TRUE)
#   
#   # Check Rhizome input and outputs
#   xyplot(data=result[[1]], 
#          Rhizome_substrate_carbon+
#            Rhizome_storage_carbon~
#            doy, auto.key=TRUE)
#   
#   xyplot(data=result[[1]][1500:5880,], 
#          substrate_transport_Stem_to_Rhizome+
#            Rhizome_storage_to_substrate_rate +
#            Rhizome_utilization_rate~
#            doy, auto.key=TRUE)
#   
#   xyplot(data=result[[1]][1300:5880,],
#            Rhizome_substrate_carbon~
#            doy, auto.key=TRUE)
#   
#   # Check why stem has low substrate C
#   xyplot(data=result[[1]],
#          substrate_transport_Leaf_to_Stem+
#            substrate_transport_Stem_to_Rhizome+
#            Stem_utilization_rate~doy, auto.key=TRUE)
#   # Conclusion: Stem used up the incoming C
#   # Two solutions: increase transport rate, decrease utilization/export rate
#   ### Check whether Stem maxed out on utilization
#   xyplot(data=result[[1]],
#          Stem_utilization_rate/Stem_structural_carbon
#          ~doy, auto.key=TRUE)
#   
#   xyplot(data=result[[1]],
#          Rhizome_utilization_rate/Rhizome_structural_carbon
#          ~doy, auto.key=TRUE)
#   
#   # Check Leaf balance
#   xyplot(data=result[[1]],
#          canopy_assimilation_rate+
#            Leaf_utilization_rate+
#            substrate_transport_Leaf_to_Stem
#          ~doy, auto.key=TRUE)
#   xyplot(data=result[[1]],
#          Leaf_utilization_rate/Leaf_structural_carbon
#          ~doy, auto.key=TRUE)
#   
#   xyplot(data=result[[1]], Rhizome_storage_to_substrate_rate~doy)
#   
#   # Check the relationship between the substrate concentration and utilization rate
#   xyplot(data=result[[1]], Leaf_substrate_carbon/Leaf_structural_carbon/10+
#            Leaf_utilization_rate/Leaf_structural_carbon~time, auto.key=TRUE)
#   xyplot(data=result[[1]], Stem_substrate_carbon/Stem_structural_carbon+
#            Stem_utilization_rate/Stem_structural_carbon~time, 
#          ylim=c(-0.05, 0.05),
#          auto.key=TRUE)
#   
#   source_minus_sink <- aggregate(data=result[[1]], cbind(Rhizome_storage_to_substrate_rate + canopy_assimilation_rate * 3.33 - (Leaf_utilization_rate+Stem_utilization_rate+Root_utilization_rate+ Rhizome_utilization_rate))~ doy, sum)
#   xyplot(data=source_minus_sink, V1~doy, auto.key = TRUE)
#   sum(result[[1]][which(result[[1]]$doy==106),"result[[1]]$Rhizome_storage_to_substrate_rate"])
#   
#   # xyplot(data=result[[1]], -Rhizome_storage_to_substrate_rate/Rhizome_utilization_rate~fractional_doy)
# }
# 



# Stem_Respiration_Factor <- c()
# Final_Yield_Average <- c()
# for (i in 1:10){
#   parameters$Stem_respiration_factor = 0.02 * i
#   final_yield = 0
#   for (i in 1:length(years)){
#     result[[i]] <- run_biocro(initial_values =  miscanthus_giganteus_initial_state,
#                               parameters = parameters,
#                               drivers = growing_season_weather,
#                               direct_module_names = miscanthus_giganteus_direct_utr_modules,
#                               differential_module_names = miscanthus_giganteus_differential_utr_modules,
#                               ode_solver = BioCro::default_ode_solvers$boost_rkck54,verbose = FALSE)
#     final_idx <- nrow(result[[i]])
#     final_yield <- final_yield + result[[i]][final_idx, 'Leaf'] + result[[i]][final_idx, 'Stem']
#   }
#   final_yield_average <- final_yield / length(years)
#   print(paste0(parameters$Stem_respiration_factor, " ", final_yield_average))
#   Stem_Respiration_Factor <- c(Stem_Respiration_Factor, parameters$Stem_respiration_factor)
#   Final_Yield_Average <- c(Final_Yield_Average, final_yield_average)
# }
# 
# # Create a data frame
# data <- data.frame(
#   Stem_Respiration_Factor = Stem_Respiration_Factor,
#   Final_Yield_Average = Final_Yield_Average
# )
# 
# # Create the scatter plot
# p <- ggplot(data, aes(x = Stem_Respiration_Factor, y = Final_Yield_Average)) +
#   geom_point(size = 3, color = "#117733") +
#   geom_line(color = "#117733", alpha = 0.5) +
#   labs(
#     title = "Effect of Stem Respiration Factor on Final Yield",
#     x = "Stem Respiration Factor",
#     y = "Final Yield Average (Mg/ha)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 10)
#   )
# print(p)
# 
# # Sink Strength
# Stem_Utilization_Rate_Change <- c()
# Final_Yield_Average <- c()
# for (i in 0:10){
#   percentage_change <- (i-5) * 0.1
#   parameters$Stem_utilization_rate_constant = 0.22 * (1 + percentage_change)
#   final_yield = 0
#   for (i in 1:length(years)){
#     result[[i]] <- run_biocro(initial_values =  miscanthus_giganteus_initial_state,
#                               parameters = parameters,
#                               drivers = growing_season_weather,
#                               direct_module_names = miscanthus_giganteus_direct_utr_modules,
#                               differential_module_names = miscanthus_giganteus_differential_utr_modules,
#                               ode_solver = BioCro::default_ode_solvers$boost_rkck54,verbose = FALSE)
#     final_idx <- nrow(result[[i]])
#     final_yield <- final_yield + result[[i]][final_idx, 'Leaf'] + result[[i]][final_idx, 'Stem']
#   }
#   final_yield_average <- final_yield / length(years)
#   print(paste0(percentage_change, " ", final_yield_average))
#   Stem_Utilization_Rate_Change <- c(Stem_Utilization_Rate_Change, percentage_change * 100)
#   Final_Yield_Average <- c(Final_Yield_Average, final_yield_average)
# }
# 
# # Create a data frame
# data <- data.frame(
#   Stem_Utilization_Rate_Change = Stem_Utilization_Rate_Change,
#   Final_Yield_Average = Final_Yield_Average
# )
# 
# # Create the scatter plot
# p <- ggplot(data, aes(x = Stem_Utilization_Rate_Change, y = Final_Yield_Average)) +
#   geom_point(size = 3, color = "#117733") +
#   geom_line(linewidth = 1, color = "#117733") +
#   labs(
#     title = "Effect of Stem Sink Strength on Yield",
#     x = "% Change of Max Stem Utilization Rate",
#     y = "Final Yield Average (Mg/ha)"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 10)
#   )
# print(p)
# 
# # Combine both
# calculate_final_yield_average <- function(params){
#   final_yield = 0
#   for (i in 1:length(years)){
#     r <- run_biocro(initial_values =  miscanthus_giganteus_initial_state,
#                     parameters = params,
#                     drivers = growing_season_weather,
#                     direct_module_names = miscanthus_giganteus_direct_utr_modules,
#                     differential_module_names = miscanthus_giganteus_differential_utr_modules,
#                     ode_solver = BioCro::default_ode_solvers$boost_rkck54,verbose = FALSE)
#     final_idx <- nrow(r)
#     final_yield <- final_yield + r[final_idx, 'Leaf'] + r[final_idx, 'Stem']
#   }
#   final_yield_average <- final_yield / length(years)
#   return(final_yield_average)
# }
# 
# Respiration_Scaling_Factor <- c()
# Utilization_Percentage_Change <- c()
# Respiration_Yield_Average <- c()
# Utilization_Yield_Average <- c()

# for (i in 1:10){
#   # For respiration: use scaling factor i
#   Respiration_Scaling_Factor <- c(Respiration_Scaling_Factor, i)
#   parameters$Stem_utilization_rate_constant = 0.22
#   parameters$Stem_respiration_factor = 0.02 * i
#   final_yield_average_res <- calculate_final_yield_average(parameters)
#   print(paste0("Respiration scaling: ", i, " Yield: ", final_yield_average_res))
#   Respiration_Yield_Average <- c(Respiration_Yield_Average, final_yield_average_res)
#   
#   # For utilization: use percentage change
#   percentage_change <- (i-5) * 0.1
#   Utilization_Percentage_Change <- c(Utilization_Percentage_Change, percentage_change * 100)
#   parameters$Stem_respiration_factor = 0.02
#   parameters$Stem_utilization_rate_constant = 0.22 * (1 + percentage_change)
#   final_yield_average_utl <- calculate_final_yield_average(parameters)
#   print(paste0("Utilization % change: ", percentage_change * 100, " Yield: ", final_yield_average_utl))
#   Utilization_Yield_Average <- c(Utilization_Yield_Average, final_yield_average_utl)
# }
# 
# # Create separate data frames for each parameter
# data_respiration <- data.frame(
#   x_value = Respiration_Scaling_Factor,
#   Yield = Respiration_Yield_Average,
#   Parameter = "Respiration"
# )
# 
# data_utilization <- data.frame(
#   x_value = Respiration_Scaling_Factor,
#   Yield = Utilization_Yield_Average,
#   Parameter = "Utilization"
# )
# 
# # Combine data
# data_combined <- rbind(data_respiration, data_utilization)
# 
# # Create the plot with two x-axes
# p <- ggplot(data_combined, aes(x = x_value, y = Yield, color = Parameter)) +
#   geom_point(size = 3) +
#   geom_line(linewidth = 1) +
#   scale_color_manual(
#     values = c("Respiration" = "#E69F00", "Utilization" = "#117733"),
#     labels = c("Respiration" = "Carbon Cost", "Utilization" = "Utilization")
#   ) +
#   scale_x_continuous(
#     name = "Carbon cost of stem biosynthesis",
#     breaks = seq(1, 10),
#     labels = 1:10,
#     sec.axis = sec_axis(~ (. - 5) * 10 , 
#                         name = "% Change in max stem utilization rate")
#   ) +
#   
#   labs(
#     y = "Yield (Mg/ha)"
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     # White background
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     
#     # Grid lines (subtle)
#     panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
#     panel.grid.minor = element_blank(),
#     
#     # Axes
#     axis.line = element_line(color = "black", linewidth = 0.7),
#     axis.ticks = element_line(color = "black", linewidth = 0.5),
#     axis.ticks.length = unit(0.2, "cm"),
#     axis.text = element_text(size = 11, color = "black"),
#     axis.title.x.bottom = element_text(size = 13, face = "bold", color = "#E69F00"),
#     axis.title.x.top = element_text(size = 13, face = "bold", color = "#117733"),
#     axis.title.y = element_text(size = 13, face = "bold", color = "black"),
#     axis.text.x.bottom = element_text(color = "#E69F00"),
#     axis.text.x.top = element_text(color = "#117733"),
#     
#     # Legend
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key = element_rect(fill = "white", color = NA),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 10),
#     legend.position = "bottom",
#     legend.key.width = unit(1.2, "cm"),
#     legend.key.height = unit(0.6, "cm"),
#     
#     # Plot margins
#     plot.margin = margin(15, 15, 15, 15)
#   ) +
#   guides(color = guide_legend(
#     override.aes = list(shape = 16, linetype = 1, linewidth = 1, size = 3)
#   ))
# 
# print(p)
# ggsave("yield_response.png", plot = p, width = 14/3, height = 11.23/3, units = "in", dpi = 300)

