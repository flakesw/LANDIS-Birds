#make maps from LANDIS outputs
#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyverse")

species <- "cerw"
model_name <- "LowTLowV BAU"
model_dir <- paste0("D:/SApps LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 60

ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

full_model <- readRDS("./models_for_landis/cerw_dist_model_full_landis.RDS")


predictor_stack_init <- terra::rast(paste0("./landis_predictor_layers/pred_stack_", 
                                      model_name, "_", year, ".tif"))%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack$understory_ratio <- 1 - predictor_stack$understory_ratio

predictor_stack_init <- terra::rast(paste0("./landis_predictor_layers/pred_stack_", 
                                           model_name, "_", 0, ".tif"))%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack_init$understory_ratio <- 1 - predictor_stack_init$understory_ratio
# 
# pred_stack_orig <- terra::rast("predictor_layers/predictor_stack.grd") %>%
#   terra::project(ecoregions) %>%
#   mask(ecoregions, maskvalues = 1)
# # writeRaster(pred_stack_orig, "pred_stack_study_area.tif")

preds <- terra::predict(object = predictor_stack, 
                        model = full_model,
                        # ext = st_bbox(bcr_albers),
                        const = data.frame(time_observations_started = 7,
                                           duration_minutes = 60)
                        )

values(preds) <- boot::inv.logit(values(preds))
preds <- crop(preds, ecoregions) %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
# plot(preds)

states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)

library("tidyterra")
ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) + 
  geom_spatraster(data = preds) +
  scale_fill_terrain_c(alpha = 0.7) + 
  # geom_sf(data = bcr_albers, fill = NA) +
  labs(fill = "P(observation)")

writeRaster(preds, paste0("./landis_predictions/prediction_", species, "_", 
                          model_name, "_", year, ".tif"), overwrite = TRUE)

#--------------------

#climate variables held constant
clim_preds <- terra::predict(object = predictor_stack[[-c(4,5,6)]], 
                                      model = full_model,
                                      # ext = st_bbox(bcr_albers),
                                      const = data.frame(time_observations_started = 7,
                                                         duration_minutes = 60,
                                                         tmmn = mean(predictor_stack_init$tmmn[]),
                                                         tmmx = mean(predictor_stack_init$tmmx[]),
                                                         pr = mean(predictor_stack_init$pr[]))
        )

values(clim_preds) <- boot::inv.logit(values(clim_preds)) 

clim_preds <- clim_preds %>%
  crop(ecoregions) %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
  
ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) + 
  geom_spatraster(data = clim_preds) +
  scale_fill_terrain_c(alpha = 0.7) + 
  # geom_sf(data = bcr_albers, fill = NA) +
  labs(fill = "P(observation)")


#----------------
#vegetation variables held constant

veg_preds <- terra::predict(object = predictor_stack[[-c(1,2,3)]], 
                             model = full_model,
                             # ext = st_bbox(bcr_albers),
                             const = data.frame(time_observations_started = 7,
                                                duration_minutes = 60,
                                                understory_ratio = mean(predictor_stack_init$understory_ratio[]),
                                                biomass = mean(predictor_stack_init$biomass[]),
                                                prop_decid = mean(predictor_stack_init$prop_decid[]))
)

values(veg_preds) <- boot::inv.logit(values(veg_preds)) 

veg_preds <- veg_preds %>%
  crop(ecoregions) %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)

ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) + 
  geom_spatraster(data = veg_preds) +
  scale_fill_terrain_c(alpha = 0.7) + 
  # geom_sf(data = bcr_albers, fill = NA) +
  labs(fill = "P(observation)")

