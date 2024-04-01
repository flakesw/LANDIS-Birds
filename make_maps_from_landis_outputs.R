#make maps from LANDIS outputs
#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyverse")

ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

full_model <- readRDS("./models_for_landis/cerw_dist_model_full_landis.RDS")


predictor_stack_0 <- terra::rast("landis_predictor_layers/pred_stack_LowTLowV BAU_0.tif")%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack_0$understory_ratio <- predictor_stack_0$understory_ratio * 12.8

predictor_stack_60 <- terra::rast("landis_predictor_layers/pred_stack_LowTLowV BAU_60.tif")%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)

pred_stack_orig <- terra::rast("predictor_layers/predictor_stack.grd") %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)

preds <- terra::predict(object = predictor_stack_0, #pred_stack_orig, 
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


#just climate variables




#just habitat variables

