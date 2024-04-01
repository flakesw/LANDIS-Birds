#make maps from LANDIS outputs


#full model
library("dismo")

full_model <- readRDS("./models_for_landis/cerw_dist_model_full_landis.RDS")

predictor_stack <- terra::rast("landis_predictor_layers/pred_stack_LowTLowV BAU_0.tif")

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


#just climate variables




#just habitat variables

