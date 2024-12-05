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
  terra::aggregate(fact = 30, fun = median, na.rm = TRUE) %>% 
  # terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster, method = "bilinear")

fhd_masked <- terra::rast("./predictor_layers/sdm/biomass_raw.tif")
NAflag(fhd_masked) <- -1
fhd_agg <- fhd_masked %>%
  terra::aggregate(fact = 10, median, na.rm = TRUE) %>% #aggregate to 250m cells
  terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)


#topography variables
chili_rast <- list.files("predictor_layers", "chili", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
tpi_rast <- list.files("predictor_layers", "tpi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
slope_rast <- list.files("predictor_layers", "slope", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)

tmax_rast <- list.files("predictor_layers", "tmax", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
tmin_rast <- list.files("predictor_layers", "tmin", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pr_rast <- list.files("predictor_layers", "ppt", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
soil_rast <- list.files("predictor_layers", "soil", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
aet_rast <- list.files("predictor_layers", "aet", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
cwd_rast <- list.files("predictor_layers", "cwd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pet_rast <- aet_rast + cwd_rast
vpd_rast <- list.files("predictor_layers", "vpd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pdsi_rast <- list.files("predictor_layers", "pdsi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)

forest_rast <- terra::rast("./predictor_layers/prop_forest.tiff")%>%
  terra::project(template_raster)
decid_rast <- terra::rast("./predictor_layers/prop_decid.tiff")%>%
  terra::project(template_raster)
conifer_rast <- terra::rast("./predictor_layers/prop_conifer.tiff")%>%
  terra::project(template_raster)
grass_shrub_rast <- terra::rast("./predictor_layers/prop_grass_shrub.tiff")%>%
  terra::project(template_raster)
oak_rast <- terra::rast("./predictor_layers/prop_oak_fia.tiff")%>%
  terra::project(template_raster)
spruce_rast <- terra::rast("./predictor_layers/prop_spruce_fia.tiff")%>%
  terra::project(template_raster)
decid_rast_fia <- terra::rast("./predictor_layers/prop_hardwood_fia.tiff")%>%
  terra::project(template_raster)

#------------------------
#layers for gbm predictions

predictor_stack <- c(biomass,
                     understory,
                     fhd_rast,
                     height_rast,
                     open_area_rast,
                     
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
                     
                     forest_rast,
                     decid_rast,
                     conifer_rast,
                     grass_shrub_rast,
                     oak_rast,
                     spruce_rast,
                     decid_rast_fia
)


names(predictor_stack) <- c("biomass",
                            "understory_ratio",
                            "fhd_normal",
                            "height",
                            "open_area",
                            
                            
                            "aet",
                            "pet",
                            "def",
                            "soil",
                            "vpd",
                            "pdsi",
                            "tmmx",
                            "tmmn",
                            "pr",
                            "tpi",
                            "chili",
                            "slope",
                            
                            
                            "prop_forest",
                            "prop_decid",
                            "prop_conifer",
                            "prop_grass", #TODO change name
                            "prop_oak",
                            "prop_spruce",
                            "prop_decid_fia")

terra::writeRaster(predictor_stack, 
                   "./predictor_layers/predictor_stack_bcr28.tif", 
                   overwrite = TRUE)
