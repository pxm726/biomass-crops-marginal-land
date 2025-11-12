# The initial values commented with "adjust" require data validation and adjustment based on the simulation
cf <- 0.3   # [Mg/ha]/[mmol/m2] # from the optimized result
subs_frac <- 0.045 # adjust
stor_frac <- 0.005 # adjust
struc_frac <- 0.95 # adjust
# For the initial total seed mass per land area, we use the following equation:
# Number of seeds per meter * weight per seed * 1 / row spacing
#
# Number of seeds per meter = 20 (Morgan et al., 2004, https://doi.org/10.1104/pp.104.043968)
# weight per seed = 0.15 grams / seed (https://www.feedipedia.org/node/42, average of .12 to .18 grams)
# row spacing = 0.38 meters (Morgan et al., 2004)
#
# (20 seeds / meter) * (0.15 grams / seed) * (1 / 0.38 meter) = 7.89 g / m^2 = 0.0789 Mg / ha
# This value is used to determine the initial Leaf, Stem, and Root biomasses
Leaf <- 1E-4
Stem <- 1E-4
Root <- 1E-4
Rhizome <- 24.1

initial_state <- list(
  
  Leaf_respiration_loss = 0.0,
  Stem_respiration_loss = 0.0,
  Root_respiration_loss = 0.0,
  Rhizome_respiration_loss  = 0.0,
  
  # senescence related initial state
  Leaf_senescence_loss  = 0.0,
  Stem_senescence_loss  = 0.0,
  Root_senescence_loss  = 0.0,
  Rhizome_senescence_loss   = 0.0,
  
  # Other variables
  TTc    = 0.0,
  soil_water_content = 0.32,
  
  cws1                  = 0.32,          # dimensionless, current water status, soil layer 1
  cws2                  = 0.32,          # dimensionless, current water status, soil layer 2
  # DVI =                     -1,             # Sowing date: DVI=-1
  
  # Substrate and structural C
  Leaf_substrate_carbon = subs_frac * Leaf / cf,
  Leaf_storage_carbon = stor_frac * Leaf / cf,
  Leaf_structural_carbon = struc_frac * Leaf / cf,
  Stem_substrate_carbon = subs_frac *  Stem / cf,
  Stem_storage_carbon = stor_frac * Stem / cf,
  Stem_structural_carbon = struc_frac * Stem / cf,
  Root_substrate_carbon =  subs_frac * Root / cf,
  Root_storage_carbon = stor_frac * Root / cf,
  Root_structural_carbon = struc_frac * Root / cf,
  
  Rhizome_substrate_carbon = 0.25 * Rhizome / cf, # adjust
  Rhizome_storage_carbon = 0.25 * Rhizome / cf, # adjust
  Rhizome_structural_carbon = 0.5 * Rhizome / cf # adjust   
)

