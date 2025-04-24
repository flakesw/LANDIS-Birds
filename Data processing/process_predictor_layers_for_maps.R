# process downloaded GEE rasters and local data to make predictor layers
# for mapping

#For all this data, you'll first need to run "download_predictor_layers.R" and
#the script that processes Treemap data and the one that does NLCD data
# We also need the "high-resolution" rasters from the GEE data. You could do
# this without running that script, but you'd have to do some pre-processing
# to mask out water and urban areas.

#--------------------------------
# Import predictor data layers
#-----------------------------------------------
#TODO fix units!
library(tidyverse)
library(terra)
library(sf)

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_albers <- bcr %>% sf::st_transform(crs = "EPSG:5070")

template_raster <- terra::rasterize(vect(bcr_albers), rast(vect(bcr_albers), resolution = 1000))
plot(template_raster)

#-----------------------------------------------------------------------------
#aggregate to coarser grid for prediction
biomass_masked <- terra::rast("./predictor_layers/sdm/biomass_raw.tif")
NAflag(biomass_masked) <- -1
biomass_agg <- biomass_masked %>%
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::aggregate(fact = 40, fun = median, na.rm = TRUE) %>% 
  terra::project(template_raster, method = "bilinear")
plot(biomass_agg)

fhd_masked <- terra::rast("./predictor_layers/sdm/fhd_raw.tif")
NAflag(fhd_masked) <- -1
fhd_agg <- fhd_masked %>%
  terra::aggregate(fact = 40, median, na.rm = TRUE) %>% #aggregate to 1000m cells
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)

height_masked <- terra::rast("./predictor_layers/sdm/height_raw.tif")
NAflag(height_masked) <- -1
height_agg <- height_masked %>%
  terra::aggregate(fact = 40, median, na.rm = TRUE) %>% #aggregate to 1000m cells
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)
plot(height_agg)

open_masked <- terra::rast("./predictor_layers/sdm/prop_open_raw.tif")
NAflag(open_masked) <- -1
open_agg <- open_masked %>%
  terra::aggregate(fact = 40, median, na.rm = TRUE) %>% #aggregate to 1000m cells
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)
plot(open_agg)

understory_masked <- terra::rast("./predictor_layers/sdm/understory_raw.tif")
NAflag(understory_masked) <- -1
understory_agg <- understory_masked %>%
  terra::aggregate(fact = 40, median, na.rm = TRUE) %>% #aggregate to 1000m cells
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)
plot(understory_agg)

#topography variables
chili_rast <- list.files("predictor_layers", "chili", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
plot(chili_rast, range = c(150,250))
tpi_rast <- list.files("predictor_layers", "tpi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
plot(tpi_rast)
slope_rast <- list.files("predictor_layers", "slope", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
plot(slope_rast)

#climate variables
tmax_rast <- list.files("predictor_layers", "tmax", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(10)
plot(tmax_rast)
tmin_rast <- list.files("predictor_layers", "tmin", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(10)
pr_rast <- list.files("predictor_layers", "ppt", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
soil_rast <- list.files("predictor_layers", "soil", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(10)
aet_rast <- list.files("predictor_layers", "aet", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(10)
cwd_rast <- list.files("predictor_layers", "cwd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(10)
pet_rast <- aet_rast + cwd_rast
vpd_rast <- list.files("predictor_layers", "vpd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(100) 
pdsi_rast <- list.files("predictor_layers", "pdsi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster) %>%
  `/`(100)


#make NLCD layers

#already cropped and reprojected
nlcd_rast <- terra::rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif") %>% as.numeric()

decid_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1), 1, 2, byrow = TRUE), others = 0)
prop_decid <- terra::aggregate(decid_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_decid <- (prop_decid/900) %>%
  terra::project(template_raster)
plot(prop_decid)
# terra::writeRaster(prop_decid, "./predictor_layers/prop_decid.tiff", overwrite = TRUE)

con_rast <- terra::classify(nlcd_rast, rcl = matrix(c(42, 1), 1, 2, byrow = TRUE), others = 0)
prop_conifer <- terra::aggregate(con_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_conifer <- (prop_conifer/900) %>%
  terra::project(template_raster)
plot(prop_conifer)
# terra::writeRaster(prop_conifer, "./predictor_layers/prop_conifer.tiff", overwrite = TRUE)

forest_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1,
                                                         42,1,
                                                         43,1,
                                                         90, 1), 4, 2, byrow = TRUE), others = 0)
prop_forest <- terra::aggregate(forest_rast, 30, fun = sum)
prop_forest <- (prop_forest/900) %>%
  terra::project(template_raster)
plot(prop_forest)
# terra::writeRaster(prop_forest, "./predictor_layers/prop_forest.tiff", overwrite = TRUE)

grass_rast <- terra::classify(nlcd_rast, rcl = matrix(c( 71,1,
                                                         81,1), 2, 2, byrow = TRUE), others = 0)
prop_grass <- terra::aggregate(grass_rast, 30, fun = sum)
prop_grass<- (prop_grass/900) %>%
  terra::project(template_raster)
plot(prop_grass)
# terra::writeRaster(prop_grass, "./predictor_layers/prop_grass.tiff", overwrite = TRUE)

shrub_rast <- terra::classify(nlcd_rast, rcl = matrix(c(52, 1), 1, 2, byrow = TRUE), others = 0)
prop_shrub <- terra::aggregate(shrub_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_shrub <- (prop_shrub/900) %>%
  terra::project(template_raster)
plot(prop_shrub)
# terra::writeRaster(prop_shrub, "./predictor_layers/prop_shrub.tiff", overwrite = TRUE)

water_rast <- terra::classify(nlcd_rast, rcl = matrix(c(11, 1), 1, 2, byrow = TRUE), others = 0)
prop_water <- terra::aggregate(water_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_water <- (prop_water/900) %>%
  terra::project(template_raster)
plot(prop_water)
# terra::writeRaster(prop_water, "./predictor_layers/prop_water.tiff", overwrite = TRUE)

dev_light_rast <- terra::classify(nlcd_rast, rcl = matrix(c(21, 1,
                                                            22,1,
                                                            82,1), 3, 2, byrow = TRUE), others = 0)
prop_dev_light <- terra::aggregate(dev_light_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_dev_light <-  (prop_dev_light/900) %>%
  terra::project(template_raster)
plot(prop_dev_light)
# terra::writeRaster(prop_dev_light, "./predictor_layers/prop_dev_light.tiff", overwrite = TRUE)

dev_heavy_rast <- terra::classify(nlcd_rast, rcl = matrix(c(23, 1,
                                                            24, 1), 2, 2, byrow = TRUE), others = 0)
prop_dev_heavy <- terra::aggregate(dev_heavy_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_dev_heavy <- (prop_dev_heavy/900) %>%
  terra::project(template_raster)
plot(prop_dev_heavy)
# terra::writeRaster(prop_dev_heavy, "./predictor_layers/prop_dev_heavy.tiff", overwrite = TRUE)

open_rast <- terra::classify(nlcd_rast, rcl = matrix(c(12,21,31,51,52,71:74,81,82,95, rep(1, 12)),
                                                     12, 2), others = 0)
prop_open <- terra::aggregate(open_rast, 30, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_open <- (prop_open/900) %>%
  terra::project(template_raster)
plot(prop_open)

nlcd_rast <- c(prop_decid,
               prop_conifer,
               prop_forest,
               prop_grass,
               prop_shrub,
               prop_water,
               prop_dev_light,
               prop_dev_heavy,
               prop_open
)


#treemap variables
prop_spruce <- terra::rast("./predictor_layers/prop_spruce_fia.tiff") %>%
  terra::aggregate(30, fun = sum, na.rm = TRUE) %>%
  `/`(900) %>%
  terra::project(template_raster)
prop_oak <- terra::rast("./predictor_layers/prop_oak_fia.tiff")%>%
  terra::aggregate(30, fun = sum, na.rm = TRUE) %>%
  `/`(900) %>%
  terra::project(template_raster)
prop_hardwood <- terra::rast("./predictor_layers/prop_hardwood_fia.tiff") %>%
  terra::aggregate(30, fun = sum, na.rm = TRUE) %>%
  `/`(900) %>%
  terra::project(template_raster)


#------------------------
#layers for gbm predictions
predictor_stack <- c(biomass_agg,
                     understory_agg,
                     fhd_agg,
                     height_agg,
                     open_agg,
                     
                     aet_rast,
                     pet_rast,
                     cwd_rast,
                     soil_rast,
                     vpd_rast,
                     pdsi_rast,
                     tmax_rast,
                     tmin_rast,
                     pr_rast,
                     tpi_rast,
                     chili_rast,
                     slope_rast,
                     
                     prop_decid,
                     prop_conifer,
                     prop_forest,
                     prop_grass,
                     prop_shrub,
                     prop_water,
                     prop_dev_light,
                     prop_dev_heavy,
                     prop_open,
                     prop_oak,
                     prop_spruce,
                     prop_hardwood
)


names(predictor_stack) <- c("biomass",
                            "understory",
                            "fhd",
                            "height",
                            "area_short",
                            
                            
                            "aet",
                            "pet",
                            "cwd",
                            "soil",
                            "vpd",
                            "pdsi",
                            "tmmx",
                            "tmmn",
                            "pr",
                            "tpi",
                            "chili",
                            "slope",
                           
                            "prop_decid",
                            "prop_conifer",
                            "prop_forest",
                            "prop_grass",
                            "prop_shrub",
                            "prop_water",
                            "prop_dev_light",
                            "prop_dev_heavy",
                            "prop_open",
                            "prop_oak",
                            "prop_spruce",
                            "prop_hardwood")

terra::writeRaster(predictor_stack, 
                   "./predictor_layers/predictor_stack_bcr28.tif", 
                   overwrite = TRUE)
