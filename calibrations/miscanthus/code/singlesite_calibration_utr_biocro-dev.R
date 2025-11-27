library(BioCro)
library(UTRMiscanthusBML)
# Adapted from Yufeng He: Dec 12, 2022 by Ximin on Sep 24, 2025

# parameters and weather inputs
source('miscanthus_utr_params.R')
source('miscanthus_utr_initial_values.R')
source('miscanthus_utr_modules.R')

miscanthus_giganteus_utr_parameters <- parameters
miscanthus_giganteus_initial_state <- initial_state


weather_all<-read.csv(
  "../data/weather_data_NASA_POWER/BioCro_input_NASA/site_2_2002_2018_LowerTransmittance.csv")
result = list()
years = 2006:2008
for (i in 1:length(years)){
  year_i = years[i]
  growing_season_weather = weather_all[weather_all$year == year_i,]
  growing_season_weather =  growing_season_weather[growing_season_weather$doy>=106 &
                                                    growing_season_weather$doy<=350,]

  result[[i]] <- run_biocro(initial_values =  miscanthus_giganteus_initial_state,
                       parameters = miscanthus_giganteus_utr_parameters,
                       drivers = growing_season_weather,
                       direct_module_names = miscanthus_giganteus_direct_utr_modules,
                       differential_module_names = miscanthus_giganteus_differential_utr_modules,
                       ode_solver = BioCro::default_ode_solvers$boost_rkck54,verbose = FALSE)
}

library(lattice)
examining_plots = TRUE
if (examining_plots){
  # plot the substrate/storage C "concentrations"
  concentrations <- with(result[[1]], data.frame(
    doy = doy,
    Leaf = Leaf_substrate_carbon/Leaf_structural_carbon,
    Stem = Stem_substrate_carbon/Stem_structural_carbon,
    Root = Root_substrate_carbon/Root_structural_carbon,
    Rhizome = Rhizome_substrate_carbon/Rhizome_structural_carbon,
    Rhizome_storage = Rhizome_storage_carbon/Rhizome_structural_carbon
  ))
  xyplot(data=concentrations,# [1500:5880,], 
         Leaf + Stem + Root + Rhizome + Rhizome_storage ~ doy, 
         ylab = "Substrate C: Structural C",
         auto.key = TRUE)
  
  xyplot(data=result[[1]], 
         Leaf_substrate_carbon + 
           Stem_substrate_carbon + 
           Root_substrate_carbon + 
           Rhizome_substrate_carbon + 
           Rhizome_storage_carbon ~ doy, 
         ylab = "Substrate C (mol C / m^2)",
         auto.key = TRUE)
  
  # plot the utilization rate compared to transport to check 
  # 1. if there is significant imbalances
  # 2. if the rhizome exports first and then import C
  # 3. if Stem_utilization_rate is high enough for the rapid growth
  # 4. if not is it limited by the stem substrate concentration or the utilization parameters?
  # 5. Is Rhizome growing at the end of the season?
  # 6. Do leaf and stem utilization max out, so more can be transported to Rhizome?
  xyplot(data=result[[1]], 
         Leaf_utilization_rate+
           substrate_transport_Leaf_to_Stem+
           Stem_utilization_rate+
           substrate_transport_Stem_to_Rhizome+
           Rhizome_utilization_rate
         ~doy, alpha = 0.5,
         ylab = 'Transport/Utilization rate (mol C/ m^2 / hr)', auto.key = TRUE)
  
  # plot the utilization rate compared to transport to check if there is significant imbalances
  transport <- with(result[[1]], data.frame(
    doy = doy,
    Leaf_Stem = substrate_transport_Leaf_to_Stem,
    Stem_Root = substrate_transport_Stem_to_Root,
    Stem_Rhizome = substrate_transport_Stem_to_Rhizome
  ))
  xyplot(data=transport, 
         Leaf_Stem + Stem_Root + Stem_Rhizome ~ doy, 
         ylab = "Transport_rate (mol C/ m^2 / hr)" , auto.key = TRUE,)
  
  # Compare canopy_assimilation_rate with storage release rate.
  xyplot(data=result[[1]],
         Rhizome_storage_to_substrate_rate+
           canopy_assimilation_rate~
           doy, auto.key=TRUE)
  
  # Check Rhizome input and outputs
  xyplot(data=result[[1]], 
         Rhizome_substrate_carbon+
           Rhizome_storage_carbon~
           doy, auto.key=TRUE)
  
  xyplot(data=result[[1]][1500:5880,], 
         substrate_transport_Stem_to_Rhizome+
           Rhizome_storage_to_substrate_rate +
           Rhizome_utilization_rate~
           doy, auto.key=TRUE)
  
  xyplot(data=result[[1]][1300:5880,],
           Rhizome_substrate_carbon~
           doy, auto.key=TRUE)
  
  # Check why stem has low substrate C
  xyplot(data=result[[1]],
         substrate_transport_Leaf_to_Stem+
           substrate_transport_Stem_to_Rhizome+
           Stem_utilization_rate~doy, auto.key=TRUE)
  # Conclusion: Stem used up the incoming C
  # Two solutions: increase transport rate, decrease utilization/export rate
  ### Check whether Stem maxed out on utilization
  xyplot(data=result[[1]],
         Stem_utilization_rate/Stem_structural_carbon
         ~doy, auto.key=TRUE)
  
  xyplot(data=result[[1]],
         Rhizome_utilization_rate/Rhizome_structural_carbon
         ~doy, auto.key=TRUE)
  
  # Check Leaf balance
  xyplot(data=result[[1]],
         canopy_assimilation_rate+
           Leaf_utilization_rate+
           substrate_transport_Leaf_to_Stem
         ~doy, auto.key=TRUE)
  xyplot(data=result[[1]],
         Leaf_utilization_rate/Leaf_structural_carbon
         ~doy, auto.key=TRUE)
  
  xyplot(data=result[[1]], Rhizome_storage_to_substrate_rate~doy)
  

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
#observed <- observed_biomass[1:4, c("doy","Stem","Leaf","Rhizome","Root","lai")]
observed = observed_biomass[, c("doy","Stem","Leaf","Rhizome","Root")]
observed = reshape2::melt(observed,id.vars = c("doy"),measure.vars = c("Leaf","Stem","Root","Rhizome"))
names(observed)  = c("x","varname","y")
names(predicted) = c("x","varname","y")
xlabtitle = expression(paste("Day of Year"))
ylabtitle = expression(paste("Dry Biomass (Mg/ha)"))
BiomassPartitioning <- 
  Compare_Observed_and_Predicted(observeddata = observed,predicteddata=predicted  ,xlabtitle,ylabtitle)
plot(BiomassPartitioning)
# ggsave(filename="./BiomassPartitioning.png",dpi = 500, width =6.5, height =4.2)

