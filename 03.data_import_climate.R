#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 03.data_import_climate
# Extract climatic data from Batalla et al 20XX?
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 2 august 2021
#######################################################

#####
# 0. Defining project characteristics
#####
# Base Layer : Growing Degree Day GDD
gdd = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/GDD/GDD_0_Pyrenees.tif")

#####
# 1. Uploading variables
#####
# a. Potential Evapo Transpiration PET
pet.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/PET/PET_Annual/PET_Annual_Pyrenees.tif")
pet.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/PET/PET_Spring/PET_Spring_Pyrenees.tif")
pet.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/PET/PET_Summer/PET_Summer_Pyrenees.tif")
pet.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/PET/PET_Autumn/PET_Autumn_Pyrenees.tif")
pet.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/PET/PET_Winter/PET_Winter_Pyrenees.tif")
  
# b. Potential Solar Radiation PSR
psr.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Pot_solar_rad/Pot_solar_rad_Annual/Pot_solar_rad_Annual_Pyrenees.tif")
psr.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Pot_solar_rad/Pot_solar_rad_Spring/Pot_solar_rad_Spring_Pyrenees.tif")
psr.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Pot_solar_rad/Pot_solar_rad_Summer/Pot_solar_rad_Summer_Pyrenees.tif")
psr.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Pot_solar_rad/Pot_solar_rad_Autumn/Pot_solar_rad_Autumn_Pyrenees.tif")
psr.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Pot_solar_rad/Pot_solar_rad_Winter/Pot_solar_rad_Winter_Pyrenees.tif")
  
# c. Precipitation PREC
prec.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Precipitation/Precipitation_Annual/Precipitation_Annual_Pyrenees.tif")
prec.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Precipitation/Precipitation_Spring/Precipitation_Spring_Pyrenees.tif")
prec.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Precipitation/Precipitation_Summer/Precipitation_Summer_Pyrenees.tif")
prec.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Precipitation/Precipitation_Autumn/Precipitation_Autumn_Pyrenees.tif")
prec.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Precipitation/Precipitation_Winter/Precipitation_Winter_Pyrenees.tif")

# d. Maximum temperature TMAX
tmax.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_max/Temperature_max_Annual/Temperature_max_Annual_Pyrenees.tif")
tmax.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_max/Temperature_max_Spring/Temperature_max_Spring_Pyrenees.tif")
tmax.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_max/Temperature_max_Summer/Temperature_max_Summer_Pyrenees.tif")
tmax.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_max/Temperature_max_Autumn/Temperature_max_Autumn_Pyrenees.tif")
tmax.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_max/Temperature_max_Winter/Temperature_max_Winter_Pyrenees.tif")
  
# e. Mean temperature TMEAN
tmean.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_mean/Temperature_mean_Annual/Temperature_mean_Annual_Pyrenees.tif")
tmean.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_mean/Temperature_mean_Spring/Temperature_mean_Spring_Pyrenees.tif")
tmean.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_mean/Temperature_mean_Summer/Temperature_mean_Summer_Pyrenees.tif")
tmean.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_mean/Temperature_mean_Autumn/Temperature_mean_Autumn_Pyrenees.tif")
tmean.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_mean/Temperature_mean_Winter/Temperature_mean_Winter_Pyrenees.tif")

# f. Minimum temperature TMIN
tmin.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_min/Temperature_min_Annual/Temperature_min_Annual_Pyrenees.tif")
tmin.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_min/Temperature_min_Spring/Temperature_min_Spring_Pyrenees.tif")
tmin.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_min/Temperature_min_Summer/Temperature_min_Summer_Pyrenees.tif")
tmin.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_min/Temperature_min_Autumn/Temperature_min_Autumn_Pyrenees.tif")
tmin.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Temperature_min/Temperature_min_Winter/Temperature_min_Winter_Pyrenees.tif")

# g. Water availability WAT
wat.an = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Water_availability/Water_availability_Annual/Water_availability_Annual_Pyrenees.tif")
wat.spr = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Water_availability/Water_availability_Spring/Water_availability_Spring_Pyrenees.tif")
wat.sum = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Water_availability/Water_availability_Summer/Water_availability_Summer_Pyrenees.tif")
wat.aut = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Water_availability/Water_availability_Autumn/Water_availability_Autumn_Pyrenees.tif")
wat.win = raster ("/media/papuga/TOSHIBA EXT/02.spatial.data/03.environement/climat/pyrenees/Water_availability/Water_availability_Winter/Water_availability_Winter_Pyrenees.tif")
  
#####
# 2. Stack
#####
env.var = stack (gdd, # Growing degree day
                 pet.an, pet.spr, pet.sum, pet.aut, pet.win, # a. Potential Evapo Transpiration PET
                 psr.an, psr.spr, psr.sum, psr.aut, psr.win,    # b. Potential Solar Radiation PSR
                 prec.an, prec.spr, prec.sum, prec.aut, prec.win,   # c. Precipitation PREC
                 tmax.an, tmax.spr, tmax.sum, tmax.aut, tmax.win,  # d. Maximum temperature TMAX
                 tmean.an, tmean.spr, tmean.sum, tmean.aut, tmean.win,  # e. Mean temperature TMEAN
                 tmin.an, tmin.spr, tmin.sum, tmin.aut, tmin.win,  # f. Minimum temperature TMIN
                 wat.an, wat.spr, wat.sum,  wat.aut, wat.win)  # g. Water availability WAT

#####
# 3. Uploading locations
#####
# Load file
localisation = read.csv (here::here("data", "raw", "site.location.csv"), 
                         head = T, row.names = NULL, sep = ";", dec = ",")
# Convert to UTM
#xy <- data.frame(ID = 1:2, X = c(118, 119), Y = c(10, 50))
coordinates(localisation) <- c("longitude", "latitude")
proj4string(localisation) = CRS("+proj=longlat +datum=WGS84")  #

coord = projection(env.var)
res = spTransform(localisation, CRS(coord))

#####
# 4. Extract climatic data and save the dataset
#####
site = as.data.frame(res)
clim.data = raster::extract (env.var, res, df = T) %>%
  mutate (combe = site[,1]) %>% 
  dplyr::select(-ID, 
         -PET_Winter_Pyrenees) # equal 0

# save the dataset
write_csv( clim.data, 
           here::here ("data", "processed", "clim.data.csv"))



