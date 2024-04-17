#make maps from LANDIS outputs
#full model
library("terra")
library("sf")
library("tidyverse")
library("dismo")
library("gbm")
library("tidyverse")
diverging_color_ramp <- function(ras){
  the_palette_fc <- leaflet::colorNumeric(palette = "RdBu", 
                                          # domain = c(-max(abs(ras[]), na.rm = TRUE), max(abs(ras[]), na.rm = TRUE)),
                                          domain = c(-1, 1),
                                          reverse = FALSE)
  # the_colors <- the_palette_fc(seq((-max(abs(ras[]), na.rm = TRUE)), max(abs(ras[]), na.rm = TRUE), length.out = 50))
  the_colors <- the_palette_fc(seq(-1, 1, length.out = 50))

}

theme_set(theme_bw())
theme_update(panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.position = "none")



states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)

species <- "gwwa"
model_name <- "LowTLowV BAU nofire"
model_dir <- paste0("D:/SApps LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 60

ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

full_model <- readRDS(paste0("./models_for_landis/", species, "_dist_model_full_landis.RDS"))


predictor_stack <- terra::rast(paste0("./landis_predictor_layers/pred_stack_", 
                                      model_name, "_", year, ".tif"))%>%
  terra::project(ecoregions) %>%
  mask(ecoregions, maskvalues = 1)
predictor_stack$understory_ratio <- 1 - predictor_stack$understory_ratio

predictor_stack_init <- terra::rast("./landis_predictor_layers/pred_stack_study_area.tif")

#making original predictor stack, onle need to do once
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


library("tidyterra")
ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) + 
  geom_spatraster(data = preds - preds_init) +
  scale_fill_terrain_c(alpha = 0.7, ) + 
  # geom_sf(data = bcr_albers, fill = NA) +
  labs(fill = "P(observation)")



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

test <- preds - preds_init
test[1] <- 0.6
plot(test, col = ramp)
ramp <- diverging_color_ramp(preds - preds_init)

writeRaster(test, "gwwa_habitat_change.tif")
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


plot(clim_preds - preds_init, col = diverging_color_ramp(clim_preds - preds_init))

ggplot() +
  geom_spatraster(data = clim_preds - preds_init) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), 
                       direction = 1, na.value = "transparent")

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


plot(veg_preds - preds_init, col = diverging_color_ramp(veg_preds - preds_init))

ggplot() +
  geom_spatraster(data = veg_preds - preds_init) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), 
                       direction = 1, na.value = "transparent")
