miscanthus_giganteus_direct_utr_modules <- c("BioCro:stomata_water_stress_linear",
                                             "BioCro:leaf_water_stress_exponential",
                                             "UTRMiscanthusBML:lai_from_structural_carbon",
                                             #"BioCro:parameter_calculator",
                                             "BioCro:solar_position_michalsky",
                                             # "BioCro:shortwave_atmospheric_scattering", # replaced by c4_canopy
                                             # "BioCro:incident_shortwave_from_ground_par", # replaced by c4_canopy
                                             # "BioCro:ten_layer_canopy_properties", # replaced by c4_canopy
                                             # "BioCro:ten_layer_c4_canopy", # replaced by c4_canopy
                                             canopy_photosynthesis = "BioCro:c4_canopy",
                                             "BioCro:stefan_boltzmann_longwave",
                                             "BioCro:canopy_gbw_thornley",
                                             "BioCro:height_from_lai",
                                             # "BioCro:ten_layer_canopy_integrator", # replaced by c4_canopy
                                             "BioCro:carbon_assimilation_to_biomass",
                                             # "BioCro:thermal_time_development_rate_calculator",
                                             "BioCro:development_index_from_thermal_time",
                                             "BioCro:format_time",
                                             "UTRMiscanthusBML:thornley_utilization_calculator_lsrr",
                                             "UTRMiscanthusBML:thornley_transport_calculator_lsrr",
                                             "UTRMiscanthusBML:thornley_biomass_calculator_lsrr")

miscanthus_giganteus_differential_utr_modules <- c("UTRMiscanthusBML:thornley_utilization_lsrr",
                                                   "UTRMiscanthusBML:thornley_transport_lsrr",
                                                   "BioCro:two_layer_soil_profile",
                                                   # "BioCro:development_index",
                                                   "BioCro:thermal_time_trilinear")