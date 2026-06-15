# The variables commented with "optimize" are UTR specific variables
parameters <- list(
  # Parameters unrelated to any module
  timestep = 1.0,
  
  default_carbon_to_mass_factor <- 0.3, # Mg/ha / (mol/m^2)
  base_utilization_rate_constant <- 0.2, # optimize 0.0075
  base_utilization_km <- 0.002, # optimize 0.02
  base_conductance <- 0.1, # optimize 0.1
  # Parameters related to the UTR model
  Leaf_carbon_to_mass_factor = default_carbon_to_mass_factor, 
  Leaf_utilization_rate_constant = 1.0 * base_utilization_rate_constant, # optimize 1.0
  Leaf_utilization_km = 2.2 * base_utilization_km, # optimize 1.0
  Leaf_respiration_factor = 0.0,				# Leaf respiration is accounted for
  # storage
  Leaf_storage_to_substrate_rate_max = 0.0,     
  Leaf_substrate_to_storage_rate_max = 0.0,     
  
  
  # by the canopy photosynthesis module
  Stem_carbon_to_mass_factor = default_carbon_to_mass_factor, 
  Stem_utilization_rate_constant = 1.1 * base_utilization_rate_constant, # optimize 1.05
  Stem_utilization_km = 0.5 * base_utilization_km, # optimize 0.01
  Stem_respiration_factor = 0.02, # optimize 0.02
  
  # storage
  Stem_storage_to_substrate_rate_max = 0.0,   
  Stem_substrate_to_storage_rate_max = 0.0,     
  
  Root_carbon_to_mass_factor = default_carbon_to_mass_factor, # optimize 1.0
  Root_utilization_rate_constant = 1 * base_utilization_rate_constant, # optimize 1.0
  Root_utilization_km = 2.0 * base_utilization_km, # optimize 0.5
  Root_respiration_factor = 0.03,
  
  # storage
  Root_storage_to_substrate_rate_max = 0.0,     
  Root_substrate_to_storage_rate_max = 0.0,  
  
  Rhizome_carbon_to_mass_factor = default_carbon_to_mass_factor,
  Rhizome_utilization_rate_constant = 0.5 * base_utilization_rate_constant, # optimize 0.9
  Rhizome_utilization_km = 2.0 * base_utilization_km, # optimize 10
  Rhizome_respiration_factor = 0.02,
  
  # storage
  Rhizome_storage_to_substrate_rate_max = 0.005,  # optimize 0.003
  Rhizome_substrate_to_storage_rate_max = 1.0,  # optimize 0.001   
  storage_release_threshold = 0.1,
  
  # storage-substrate conversion
  storage_to_substrate_hill_coefficient = 1,
  substrate_to_storage_hill_coefficient = 1.56,
  storage_to_substrate_km = 0.01,
  substrate_to_storage_km = 0.0001,
  
  substrate_conductance_Leaf_to_Stem = 50 * base_conductance, #10
  substrate_conductance_Stem_to_Rhizome = 30 * base_conductance, # optimize 3.0 # Stem_to_Root 10
  substrate_conductance_Rhizome_to_Root = 20 * base_conductance, # optimize 5.0 # Stem_to_Rhizome 20
  
  # early_growth_end_dvi =                   0.3, # optimize 0.3
  # substrate_use_refill_ratio =             0.1, # optimize 0.1
  
  transportation_beta_exponent = 1,
  
  # senescence_coefficient_logistic module
  Leaf_senescence_fraction_max =           0.0005, # optimize 0.0005
  Stem_senescence_fraction_max =           0.0001, # optimize 0.0001
  Root_senescence_fraction_max =           0, # 0.0005,  # optimize 0.0005/0
  Rhizome_senescence_fraction_max  =       0.000, # optimize 0.000
  Leaf_senescence_alpha    =               1.0,   # optimize 1.0
  Stem_senescence_alpha    =               0.5,  # optimize 0.5
  Root_senescence_alpha    =               0.5,    # optimize 0.5
  Rhizome_senescence_alpha     =           0.1, # optimize 0.1
  Leaf_senescence_beta     =               1.5, # optimize 1.5
  Stem_senescence_beta     =               1.8, # optimize 1.8
  Root_senescence_beta     =               2.0, # optimize 2.0
  Rhizome_senescence_beta      =           1, # optimize 1
  Leaf_senescence_reuse_factor =           0.9, # optimize 0.9
  Stem_senescence_reuse_factor =           0.9, # optimize 0.9
  Root_senescence_reuse_factor =           0.9, # optimize 0.9
  Rhizome_senescence_reuse_factor =        0, # optimize 0
  
  
  # Parameters related to the `parameter_calculator` module
  iSp = 1.7, # updated 
  Sp_thermal_time_decay = 0.00038, # updated
  LeafN_0 = 2.0,
  LeafN = 2.0,
  vmax_n_intercept = 0,
  Vcmax_at_25 = 39, # updated # Renamed from "vmax1"
  alphab1 = 0,
  alpha1 = 0.04, # 0 in biomass-crops-marginal-land, 0.04 in miscanthus_x_giganteous.R, 0.86 in Collatz et al. 1992 paper.
  
  # carbon_assimilation_to_biomass
  dry_biomass_per_carbon = 30.026, # # g CHO / mol C (glucose)
  
  # soil parameters (clay loam)
  soil_air_entry              = -2.6,
  soil_b_coefficient          = 5.2,
  soil_bulk_density           = 1.35,
  soil_clay_content           = 0.34,
  soil_field_capacity         = 0.32,
  soil_sand_content           = 0.32,
  soil_saturated_conductivity = 6.4e-05,
  soil_saturation_capacity    = 0.52,
  soil_silt_content           = 0.34,
  soil_wilting_point          = 0.2,
  
  # Parameters related to the `soil_evaporation` module
  rsec = 0.2,
  soil_clod_size = 0.04,
  soil_reflectance = 0.2,
  soil_transmission = 0.01,
  specific_heat_of_air=1010,
  stefan_boltzman = 5.67e-8,
  
  # thermal_time_development_rate_calculator module
  sowing_fractional_doy =                  0      , # from miscanthus_x_giganteous.R, non-existing in yufeng's biomass-crops-marginal-land
  TTemr =                                  400    ,  # Emergence, Tejera & Heaton 2019 ; 400 from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  TTveg =                                  1600   ,  # GGD6 = 2200, Assume it is stem elongation. Tejera & Heaton 2019; 900 from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  TTrep =                                  400   , # GGD6 = 2600, Leaf growth stop, Tejera & Heaton 2019 # from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  
  # thermal_time_trilinear
  topt_lower =                            28, # from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  topt_upper =                            31, # from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  tmax =                                  40, # from yufeng's biomass-crops-marginal-land, non existing in miscanthus_x_giganteous.R
  
  # incident_shortwave_from_ground_par module
  par_energy_fraction            =        0.5,
  par_energy_content             =        0.235,
  
  # height_from_lai module
  heightf                     = 3,           # Yufeng
  
  # canopy_gbw_thornley module
  min_gbw_canopy              = 0.005,       # m / s
  
  # stefan_boltzmann_longwave module
  emissivity_sky              = 1,
  
  # two_layer_soil_profile module
  soil_depth1=                             0.0,         # meters
  soil_depth2=                             2.5,         # meters
  soil_depth3=                             10.0,        # meters
  wsFun     =                              2,           # not used, but must be defined
  hydrDist  =                              0 ,          # same as in sorghum parameter file
  rfl       =                              0.2,         # same as in sorghum parameter file
  rsdf      =                              0.44,        # same as in sorghum parameter file
  phi1      =                              0.01,
  phi2      =                              1.5,          # from Sugarcane-BioCro, Jaiswal et al. 2017 (https://doi.org/10.1038/nclimate3410)
  
  # thermal_time_linear module
  tbase     =                              6,          # degrees C, Tejera & Heaton 2017.(Originally 10)
  
  # solar_position_michalsky module
  lat                         =            40,
  longitude                   =            -88,
  time_zone_offset            =            -6,
  
  # shortwave_atmospheric_scattering module
  atmospheric_pressure  =                  101325,
  atmospheric_transmittance =              0.6,       # Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Pg 173
  atmospheric_scattering  =                0.3, 
  
  # ten_layer_canopy_properties module
  chil                        = 1,           # 1 in miscanthus_x_giganteous.R,  0.81, Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Table 15.1, pg 253
  k_diffuse                   = 0.1,         # Estimated from Campbell and Norman, An Introduction to Environmental Biophysics, 2nd Edition, Figure 15.4, pg 254 # 0.1
  kpLN                        = 0.2,        
  leaf_reflectance_nir        = 0.38,        # Corn values from Table 7-1 from Norman & Arkebauer (1991) https://doi.org/10.2134/agronmonogr31.c7
  leaf_reflectance_par        = 0.09,        # Corn values from Table 7-1 from Norman & Arkebauer (1991) https://doi.org/10.2134/agronmonogr31.c7
  leaf_transmittance_nir      = 0.45,        # Corn values from Table 7-1 from Norman & Arkebauer (1991) https://doi.org/10.2134/agronmonogr31.c7
  leaf_transmittance_par      = 0.04,        # Corn values from Table 7-1 from Norman & Arkebauer (1991) https://doi.org/10.2134/agronmonogr31.c7
  lnfun                       = 0,           # not used in Soybean-BioCro
  
  # ten_layer_c4_canopy module
  # jmax                        = 195,         # Bernacchi et al. 2005 (https://doi.org/10.1007/s00425-004-1320-8), 2002 Seasonal average
  # jmax_mature                 = 195,         # Needed in the varying_Jmax25 module
  sf_jmax                     = 0.2,         # Scaling factor for jmax. Needed in the varying_Jmax25 module
  electrons_per_carboxylation = 4.5,         # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  electrons_per_oxygenation   = 5.25,        # Bernacchi et al. 2003 (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  tpu_rate_max                = 13,          # Fitted value based on the A-Ci data measured at UIUC in 2019-08 by Delgrado (unpublished data)
  RL_at_25                    = 1.28,        # Davey et al. 2004 (https://doi.org/10.1104/pp.103.030569), Table 3, cv Pana, co2 368 ppm # Renamed from `Rd`
  Catm                        = 400,      # micromol / mol, CO2 level in 2002
  O2                          = 210,         # millimol / mol
  b0                          = 0.08,        # miscanthus_x_giganteus.R
  b1                          = 3,           # miscanthus_x_giganteus.R
  Gs_min                      = 1e-3,
  theta                       = 0.83,        # miscanthus_x_giganteous.R
  windspeed_height            = 10,
  leafwidth                   = 0.04,       
  beta                        = 0.93,        # from miscanthus_x_giganteus.R in data file, no reference
  kparm                       = 0.7,         # from miscanthus_x_giganteus.R in data file, no reference
  lowerT                      = 3,           # from miscanthus_x_giganteus.R in data file, no reference, 10 in yufeng's biomass-crops-marginal-land setting
  upperT                      = 37.5,        # from miscanthus_x_giganteus.R in data file, no reference
  # ten_layer_canopy_integrator module
  
  growth_respiration_fraction = 0,
  # c4_canopy, all from miscanthus_x_giganteous.R
  nalphab0 = 0.02367,
  nalphab1 = 0.000488,
  nileafn = 85,
  nkln = 0.5,
  nkpLN = 0.17,
  nlayers = 10,
  nlnb0 = -5,
  nlnb1 = 18,
  nRdb0 = -4.5917,
  nRdb1 = 0.1247,
  nvmaxb0 = -16.25,
  nvmaxb1 = 0.6938
)