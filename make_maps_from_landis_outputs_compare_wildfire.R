#make maps from LANDIS outputs
#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyterra")

species <- "gwwa"
model_name <- "HighTHighV BAU extra rx"
model_dir <- paste0("C:/Users/swflake/Documents/SApps-LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 60

ecoregions <- terra::rast("C:/Users/swflake/Documents/SApps-LANDIS/Inputs/Basic_inputs/Ecos11_NCLD.tif")

full_model <- readRDS("./models_for_landis/gwwa_dist_model_full_landis.RDS")


predictor_stack <- terra::rast(paste0("./landis_predictor_layers/pred_stack_", 
                                      model_name, "_", year, ".tif"))%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack$understory_ratio <- 1 - predictor_stack$understory_ratio

predictor_stack_init <- terra::rast(paste0("./landis_predictor_layers/pred_stack_",
                                           "HighTHighV BAU", "_", 60, ".tif"))%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack_init$understory_ratio <- 1 - predictor_stack_init$understory_ratio


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
plot(preds)

preds_init <- terra::predict(object = predictor_stack_init, 
                                      model = full_model,
                                      # ext = st_bbox(bcr_albers),
                                      const = data.frame(time_observations_started = 7,
                                                         duration_minutes = 60)
)

values(preds_init) <- boot::inv.logit(values(preds_init))
preds_init <- crop(preds_init, ecoregions) %>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
plot(preds_init)

# states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
#   sf::st_transform(crs = crs(preds)) %>%
#   st_crop(preds)
ggplot() +
  geom_spatraster(data = preds - preds_init) +
  scale_fill_distiller(palette = "RdBu", limits = c(-0.6, 0.6), 
                       direction = 1, na.value = "transparent")
hist((preds - preds_init)[])
mean((preds - preds_init)[], na.rm = TRUE)

