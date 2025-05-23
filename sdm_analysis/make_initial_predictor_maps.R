# Make initial SDM predictor maps
library(terra)
library(sf)
library(tidyverse)
library("dismo")
library("gbm")


## make maps ------------------------------------
ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_albers <- bcr %>% 
  st_transform("EPSG:5070")
# st_transform(crs(ecoregions)) #%>%
# st_crop(ecoregions) #

predictor_stack <- terra::rast("predictor_layers/predictor_stack_bcr28.tif") 
# predictor_stack[[c(2,3,4)]] <- predictor_stack[[c(2,3,4)]]/1000
# predictor_stack[[c(4)]] <- predictor_stack[[c(4)]] #TODO check on this?
# if(crs(predictor_stack) != crs(bcr_albers)) predictor_stack <- project(predictor_stack, bcr_albers)

names(predictor_stack)[8] <- "def"

# names(predictor_stack)[2] <- "understory_proportion"


mods <- list.files("./sdm_analysis/sdms/", full.names = TRUE) %>%
  # `[`(grepl("gwwa|cerw", .)) %>%
  # `[`(!grepl("all_scales", .))
  `[`(grepl("all_scales", .))


for(i in 1:length(mods)){
  mod_list <- readRDS(mods[i])
  
  for(j in 1:nrow(mod_list)){
    
    modname <- mod_list[[1]][j]
    mod1 <- mod_list[[2]][[j]]
    
    preds <- terra::predict(object = predictor_stack, model = mod1,
                            ext = st_bbox(bcr_albers),
                            const = data.frame(time_observations_started = 10,
                                               duration_minutes = 60,
                                               effort_radius_m = 1000)
    )
    
    values(preds) <- boot::inv.logit(values(preds))
    preds <- terra::crop(preds, vect(bcr_albers), mask = TRUE)
    plot(preds)
    writeRaster(preds, paste0("./sdm_analysis/outputs/prediction_maps/", modname, ".tiff"), overwrite = TRUE)
  }
}
