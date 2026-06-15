# Ximin Piao: Jun 10, 2026

library(ncdf4)
# library(BioCroMis)
library(BioCro) # added
library(UTRMiscanthusBML) # added
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggmap)

run_cwrfsoilwater=FALSE

ttc_function<-function(temp_array){
  tbase = 6 # originally 10
  topt_lower = 28
  topt_upper = 31
  tmax = 40
  TTc = 0
  for (i in 1:length(temp_array)){
    temp = temp_array[i]
    if (temp <= tbase){
      gdd_rate = 0.0
    }else if (temp <= topt_lower){
      gdd_rate = temp - tbase
    }else if (temp > topt_lower & temp < topt_upper){
      gdd_rate = topt_lower-tbase
    }else if (temp >= topt_upper &  temp < tmax){
      gdd_rate = (tmax - temp) * (topt_lower - tbase) / (tmax - topt_upper)
    }else{
      gdd_rate = 0.0
    }
    #  Normalize to a rate per hour
    gdd_rate = gdd_rate/ 24.0
    
    TTc = TTc + gdd_rate
  }
  return(TTc)
}
rmse<-function(obs,pred){
  rmse = sqrt(sum((pred-obs)^2)/length(obs))
}
source("output_for_different_harvest_dates.R")

nasa_path = "../data/weather_data_NASA_POWER/BioCro_input_NASA/"
################################ COMMENT 1 STARTS ###############################
# parameters inputs
#################################################################################
source('miscanthus_utr_params.R')
source('miscanthus_utr_initial_values.R')
source('miscanthus_utr_modules.R')

miscanthus_giganteus_utr_parameters <- parameters
initial_state <- initial_state
#################################COMMENT 1 ENDS #################################
# experimental data
multisite <- read.csv("../data/biomass_observation/Miscanthus_Observation_20230529.csv")
multisite <- multisite[!is.na(multisite$biomass),] #remove NA observed biomass

Num_obs <- dim(multisite)[1]

multisite$predicted_stem_actual  = numeric(Num_obs)
multisite$predicted_stem_peak    = numeric(Num_obs)
multisite$predicted_stem_march   = numeric(Num_obs)

years = 2002:2018
TTc_urbana=data.frame(year = years,total_TTC=NA)
weather_urbana = read.csv(paste0(nasa_path,'site_',2,'_2002_2018_LowerTransmittance.csv'))
for (i in 1:length(years)){
  tmp = weather_urbana$temp[weather_urbana$year==years[i]]
  TTc_urbana$total_TTC[i] = ttc_function(tmp)
}

parameters0 = parameters

unique_IDs = unique(multisite$LocationID)
ttc_scaling_factor_all=c()

if(run_cwrfsoilwater){
  miscanthus_giganteus_differential_utr_modules = miscanthus_giganteus_differential_utr_modules[-3] #remove two_layer_soil_profile
  initial_state =
    initial_state[names(initial_state)!=c('cws1', 'cws2')]
}

rhizome_winter_loss = 0.34 # 34% rhizome dies during winter
rhizome_winter_loss = 0

#create a column for predicted stem (with winter loss)
multisite$predicted_stem_actual = NA

#go through each site
for (i in 1:length(unique_IDs)){
    unique_ID = unique_IDs[i]
    weather_all = read.csv(paste0(nasa_path,'site_',unique_ID,'_2002_2018_LowerTransmittance.csv'))
    
    print(c("processing site",unique_ID))
    
    site_i = multisite[multisite$LocationID==unique_ID,]
    
    #years/duration of simulations
    planted_year <- site_i$planted_year
    measure_year <- site_i$measure_year
    stand_age    <- site_i$age     #this is actual age of measurements
    
    #maximum number of years
    number_of_years <- max(measure_year) - planted_year[1] + 1
    age_all = 1:number_of_years
    
    #biomass for each age at site_i
    biomass_site_i = (1:number_of_years)*NA
    count_index = 1  #first time the loop age reaches the first value of the measured age
    
    #initial rhizome, this should fo inside age loop to account for multiple years of growth
    initialRhizome = 0.3
    #always start with the first year stand and do a multi-year continuous simulation
    for (age in age_all){
      
      year_current = planted_year[1] + age -1
      weather <- weather_all[weather_all$year==year_current,]
      growing_season_weather <- get_growing_season_climate(weather, threshold_temperature = 0)
      
      #read in CWRF soil water for the initial soil water content and soil type
      soil_data = 
        read.csv(paste0("../data/cwrf_soil_data_1m/site_",unique_ID,"/cwrf_soilwater_",year_current,".csv"))
      
      if(run_cwrfsoilwater){
        growing_season_weather$cws1 = soil_data$swc[soil_data$doy>=growing_season_weather$doy[1]
                                                    & soil_data$doy<=tail(growing_season_weather$doy,1)]
        growing_season_weather$cws2 = soil_data$swc[soil_data$doy>=growing_season_weather$doy[1]
                                                          & soil_data$doy<=tail(growing_season_weather$doy,1)]
      }

      total_ttc = ttc_function(growing_season_weather$temp)
      ttc_scaling_factor = total_ttc / TTc_urbana$total_TTC[TTc_urbana$year==year_current]
      parameters$TTemr = parameters0$TTemr * ttc_scaling_factor
      parameters$TTveg = parameters0$TTveg * ttc_scaling_factor
      parameters$TTrep = parameters0$TTrep * ttc_scaling_factor
      ttc_scaling_factor_all = c(ttc_scaling_factor_all,ttc_scaling_factor)
      
      initial_state$Rhizome_substrate_carbon  = 0.1 * initialRhizome / cf
      initial_state$Rhizome_storage_carbon    = 0.3 * initialRhizome / cf
      initial_state$Rhizome_structural_carbon = 0.6 * initialRhizome / cf 
      
      parameters$soil_depth  = site_i$bedrock[1]
      parameters$soil_depth2 = site_i$bedrock[1]/2
      parameters$soil_depth3 = site_i$bedrock[1]
      
      parameters$soil_type_indicator = site_i$biocro_soiltype[1] + 1 #  soil_data$soiltype[1] 
      soil_params <- BioCro::soil_parameters[[parameters$soil_type_indicator]]
      parameters[names(soil_params)] <- soil_params
  
      result <- run_biocro(initial_values =  initial_state,
                           parameters = parameters, 
                           drivers = growing_season_weather,
                           direct_module_names = miscanthus_giganteus_direct_utr_modules,
                           differential_module_names = miscanthus_giganteus_differential_utr_modules,
                           ode_solver = BioCro::default_ode_solvers$boost_rkck54,verbose = FALSE)
      
      initialRhizome = result$Rhizome[dim(result)[1]] *(1-rhizome_winter_loss)
      print(paste("year=", year_current))
      print(paste("age=", age, "ttc_scaling_factor=",ttc_scaling_factor,"initialRhizome=",initialRhizome))
      # print(paste("mean temp=", mean(growing_season_weather$temp)))
      # print(paste("mean solar=", mean(growing_season_weather$solar)))
      # print(paste("mean stomaWS=", mean(result$StomataWS)))
      # print(paste("mean swc=", mean(growing_season_weather$soil_water_content)))
      #calculate winter loss of stem biomass
      #this is only valid when the loop age is matching an actual measured age
      if(age %in% stand_age){
        
        
        # #derive from Fig. 10 Midwest https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12929
        # #it makes no difference here since the oldest stand is just 7-year old
        # #will use this in the multi-site validation script
        #   aging_loss = -0.7101*max(7,age) + 4.8145 #at least age 7. Negative values
        # 
        # if(age<=3){
        #   #for young stands, we apply a constant loss of 33%
        #   biomass_site_i[age] = result$Stem[dim(result)[1]] * 0.67
        # }else{
        # #for older stands, we apply a loss of 0.07 t/ha per day
        #   if(is.na(site_i$measure_month_n[count_index])){ #if no measure month, assume the March 1st (DOY=60)
        #     site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*(60+365-tail(growing_season_weather$doy,1))
        #   }else{
        #     if(site_i$measure_month_n[count_index]>6){ #if measure month is later than June, it suggests harvesting in Winter
        #       measure_doy = site_i$measure_doy[count_index]
        #       measure_doy = as.character(measure_doy)
        #       measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
        #       number_of_loss_days = measure_doy-tail(growing_season_weather$doy,1)+1
        #       site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
        #     }else{              #if measure month is less than June, it suggests harvesting in Spring
        #       measure_doy = site_i$measure_doy[count_index]
        #       measure_doy = as.character(measure_doy)
        #       measure_doy = as.numeric(format(as.Date(measure_doy,format = "%Y%m%d"),"%j"))
        #       number_of_loss_days = 365-tail(growing_season_weather$doy,1)+1+measure_doy
        #       site_i$predicted_stem_actual[count_index] = result$Stem[dim(result)[1]] - 0.07*number_of_loss_days
        #     }
        #     if(number_of_loss_days<0) stop(c('age is',age)) # number_of_loss_days never calculated in this bracket before...?
        #     site_i$predicted_stem_actual[count_index]  = site_i$predicted_stem_actual[count_index] # site_i$predicted_stem_actual[count_index] never defined in this bracket before...?
        #     #* site_i$correction_factor[count_index]
        #   }
        #   biomass_site_i[age] = site_i$predicted_stem_actual[count_index] + aging_loss
        # }
        
        biomass_site_i[age] = result$Stem[dim(result)[1]]
        count_index = count_index+1
      }else{
        biomass_site_i[age] = NaN
      }
      # if(age==2) stop()
      print(paste0('Predicted: ', result$Stem[dim(result)[1]]))
      if(age %in% stand_age) {
        print(paste0('Measured: ', site_i$biomass[which(site_i$age==age)]))
      } else {
        print('No measurement this year.')
      }
    }
    #bind the biomass back to the full multisite matrix for plotting later
    #note that the multisite matrix does not need all ages' records
    multisite$predicted_stem_actual[multisite$LocationID==unique_ID] = biomass_site_i[age_all%in%stand_age]
}
multisite_clean = multisite
uniqueID <- unique(multisite_clean$LocationID)
uniqueID = uniqueID[!is.na(uniqueID)]

observed = numeric(length(uniqueID))
obs_se   = numeric(length(uniqueID))
predicted =  numeric(length(uniqueID))
ii=1
for (i in uniqueID){
  print(i)
  tmp = multisite_clean[multisite_clean$LocationID==i,]
  tmp = tmp[!is.na(tmp$biomass),]
  observed[ii] = mean(tmp$biomass)
  obs_se[ii]   = sd(tmp$biomass)/sqrt(length(tmp$biomass))
  predicted[ii] = mean(tmp$predicted_stem_actual)
  ii=ii+1
}

# mean_rmse = rmse(observed,predicted) # 
# mean_ccc =  epi.ccc(observed,predicted)$rho.c$est  #

cc_utr = cor(observed,predicted)
print(cc_utr)
slope_parameters  = lm(predicted~observed+0) # 
meandataset = data.frame(obs=observed, pred = predicted, SE=obs_se, location=uniqueID)
actual_stem_comparison <- ggplot(data =  meandataset,aes(x = obs ,  y = pred))+ 
  geom_point() +
  stat_smooth(method="lm",formula=y~0+x, level = 0.9) +
  geom_errorbarh(aes(xmin=obs-SE, xmax=obs+SE), height=.2,position=position_dodge(.9)) +
  xlab ("Observed Yield (Mg/ha)") + ylab("Predicted Yield (Mg/ha)") +
  ylim(0,40) + xlim(0,40)  +
  geom_abline() +
  geom_text(aes(label=location),hjust=0,vjust=0)+
  geom_text(x=28 , y = 22 , label = "1:1 line" ) + 
  geom_text(x=10 , y = 30 , label = paste0("Predicted = ",round(slope_parameters$coefficients,2)," × Observed"), color="blue")

write.csv(meandataset, file="./multisite_validation_data_utr.csv")
plot(actual_stem_comparison)
ggsave(actual_stem_comparison, filename = "utr_multisite_miscanthus_validation.png", dpi = 500, width =6.5, height =4.2)



# Create another plot
plot_data <- multisite %>%
  select(LocationID, age, biomass, predicted_stem_actual, measure_month) %>%
  pivot_longer(cols = c(biomass, predicted_stem_actual),
               names_to = "type", values_to = "value") %>%
  mutate(
    type = recode(type,
                  "biomass" = "Measured",
                  "predicted_stem_actual" = "Predicted"),
    LocationID = factor(LocationID),
    measure_month = factor(measure_month, levels = 1:12)  # force all 12 levels
  )

month_colors <- c(
  "#2166AC", # Jan
  "#4393C3", # Feb
  "#74C476", # Mar
  "#41AB5D", # Apr
  "#006D2C", # May
  "#FFEDA0", # Jun
  "#FED976", # Jul
  "#FEB24C", # Aug
  "#F03B20", # Sep
  "#BD0026", # Oct
  "#6A0DAD", # Nov
  "#08306B"  # Dec
)
present_months <- as.character(sort(unique(multisite$measure_month[!is.na(multisite$measure_month)])))
# Compute correlation per LocationID
cor_labels <- multisite %>%
  group_by(LocationID) %>%
  summarise(cor = cor(biomass, predicted_stem_actual, use = "complete.obs")) %>%
  mutate(label = paste0("L", LocationID, " r=", round(cor, 2))) %>%
  pull(label, name = LocationID)  # named vector: names are LocationID, values are labels


ggplot(plot_data, aes(x = age, y = value, group = type, linetype = type)) +
  geom_line(color = "grey50") +
  geom_vline(xintercept = 7, color = "red", linetype = "solid")+
  geom_point(
    data = filter(plot_data, type == "Measured"),
    aes(color = measure_month), size = 3
  ) +
  geom_point(
    data = filter(plot_data, type == "Predicted"),
    color = "grey30", size = 2
  ) +
  scale_linetype_manual(values = c("Measured" = "dotted", "Predicted" = "solid"),
                        name = "Type") +
  scale_color_manual(
    values  = setNames(month_colors, as.character(1:12)),
    limits  = as.character(1:12),
    breaks  = present_months,
    labels  = month.abb[as.numeric(present_months)],  # only labels for present months
    name    = "Measure\nMonth"
  ) +
  scale_x_continuous(breaks = sort(unique(multisite$age)))+
  facet_wrap(~ LocationID, labeller = label_both, scales = "free_y") +
  # facet_wrap(~ LocationID, labeller = labeller(LocationID = cor_labels), scales = "free_y") +
  labs(x = "Stand Age (years)", y = "Stem Biomass (Mg/ha)") +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    strip.text      = element_text(size = 9, face = "bold"),
    legend.position = "right",
    axis.text       = element_text(size = 8),
    panel.spacing   = unit(0.5, "lines")
  )

# Some other statistical tests to figure out which factor contributes to the error the most
multisite <- multisite %>%
  mutate(
    error     = predicted_stem_actual - biomass,  # signed error
    abs_error = abs(error),
    rel_error = error / biomass                   # relative error
  )
lm_error <- lm(abs_error ~ age + lat + lon + measure_month, data = multisite)
summary(lm_error)
# Output:
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    7.88275   16.59303   0.475   0.6357  
# age            0.55945    0.24967   2.241   0.0272 *
#   lat           -0.04833    0.41322  -0.117   0.9071  
# lon            0.02918    0.10927   0.267   0.7900  
# measure_month  0.27385    0.15853   1.727   0.0871 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.987 on 104 degrees of freedom
# Multiple R-squared:    0.1,	Adjusted R-squared:  0.06539 
# F-statistic: 2.889 on 4 and 104 DF,  p-value: 0.02587


# Variance partitioning: how much does each factor explain?
library(relaimpo)
calc.relimp(lm_error, type = "lmg")  # R2 contribution of each predictor
# Output
# Response variable: abs_error 
# Total response variance: 38.35834 
# Analysis based on 109 observations 
# 
# 4 Regressors: 
#   age lat lon measure_month 
# Proportion of variance explained by model: 10%
# Metrics are not normalized (rela=FALSE). 
# 
# Relative importance metrics: 
#   
#   lmg
# age           0.0570380900
# lat           0.0019153240
# lon           0.0008616078
# measure_month 0.0401919391
# 
# Average coefficients for different model sizes: 
#   
#   1X         2Xs         3Xs         4Xs
# age            0.68624519  0.64189546  0.60023695  0.55945120
# lat           -0.27276268 -0.18106237 -0.10444453 -0.04832753
# lon            0.04004391  0.03047873  0.02875634  0.02918130
# measure_month  0.37382027  0.33994818  0.30701265  0.27385331


aov_soil   <- aov(abs_error ~ factor(soil_type), data = multisite)
aov_source <- aov(abs_error ~ Source,            data = multisite)
summary(aov_soil)
summary(aov_source)

# Post-hoc: which groups differ?
TukeyHSD(aov_soil)

lm_full <- lm(abs_error ~ age + lat + lon + factor(soil_type) + 
                measure_month + Source, data = multisite)
summary(lm_full)

# Fine, smectitic, mesic Typic Albaqualfs — +16 Mg/ha error (p=0.009), the worst-performing soil
# Fine-silty, mixed, mesic Oxyaquic Hapludalfs — +14.7 Mg/ha error (p=0.016), second worst

# Step-wise model selection
library(MASS)
lm_step <- stepAIC(lm_full, direction = "both")
summary(lm_step)

# Output
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10.4003  -4.0600  -0.8461   2.8769  15.2960 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     3.3298     1.7661   1.885   0.0621 .
# age             0.5616     0.2471   2.273   0.0251 *
#   measure_month   0.2772     0.1528   1.814   0.0726 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.934 on 106 degrees of freedom
# Multiple R-squared:  0.09901,	Adjusted R-squared:  0.08201 
# F-statistic: 5.824 on 2 and 106 DF,  p-value: 0.003983

# The model explains only ~10% of variance in prediction error (R² = 0.099), 
# meaning most of what drives prediction accuracy is not captured by age, 
# location, soil type, or harvest month. 
# The model is statistically significant (p=0.004) but weak.
