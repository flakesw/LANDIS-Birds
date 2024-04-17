### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
library("terra")
library("sf")
#--------------------------
species <- "gwwa"

###### BRT
data <- list.files("./environment_vars/", pattern = species, full.names = TRUE) %>%
  `[`(grep("combined", .))
combined <- read.csv(data[3])
combined <- combined[complete.cases(combined), ] 
# combined$species_observed <- ifelse(as.character(combined$species_observed) == "FOUND", TRUE, FALSE)
brt <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = c("tmmx", "tmmn", "pr", #"aet", "pet", "def", "soil","pdsi","vpd", 
                                 "tpi", "chili", "slope", 
                                 "biomass", "understory_ratio", #"fhd_normal", "height", 
                                "prop_decid", #"prop_conifer", "prop_grass", "prop_spruce", "prop_oak", "prop_forest",  
                                 "time_observations_started", "duration_minutes"), #, "effort_radius_m"
                       gbm.y = "species_observed",
                       interaction.depth = 2,
                       cv_folds = 10)
print(brt)
summary.gbm(brt)
gbm.interactions(brt)
gbm.plot(brt)
gbm.perspec(brt, x = 15, y = 17, theta = 210)
gbm.perspec(brt, x = 12, y = 13)
gbm.perf(brt)
brt$self.statistics$discrimination
brt$cv.statistics$discrimination.mean
roc(combined$species_observed, boot::inv.logit(predict(brt)),
    plot = TRUE)
saveRDS(brt, paste0("./models_for_landis/", species, "_dist_model_full_landis.RDS"))


brt_simp <- dismo::gbm.simplify(brt)
plot(brt_simp$deviance.summary$mean ~ c(1:length(brt_simp$deviance.summary$mean)))
brt2<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = brt_simp$pred.list[[2]], #TODO get this automatically from the simp model
                       gbm.y = "species_observed")

saveRDS(brt2, paste0("./models_for_landis/", species, "_dist_model_simplified.RDS"))


gbm.plot(brt2,
         n.plots = 9,
         common.scale = TRUE,
         show.contrib = TRUE,
         plot.layout = c(3,3))
gbm.perspec(brt, x = 7, y = 6)


brt_veg <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                           gbm.x = c("prop_decid", "biomass", "understory_ratio",
                                     "time_observations_started", "duration_minutes"),
                           gbm.y = "species_observed")
summary.gbm(brt_veg)
gbm.plot(brt_veg)
saveRDS(brt_veg, paste0("./models_for_landis/", species, "_dist_model_veg_landis.RDS"))

brt_clim <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("tmmn", "pr", "tmmx", "slope", "tpi", "chili",
                                      "time_observations_started", "duration_minutes"),
                            gbm.y = "species_observed")
summary.gbm(brt_clim)
gbm.plot(brt_clim)
saveRDS(brt_clim, paste0("./models_for_landis/", species, "_dist_model_clim_landis.RDS"))

brt_reduced <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("biomass", "understory_ratio", "prop_decid",
                                      "tmmx", "tmmn", "pr", "chili", "tpi",
                                      "time_observations_started", "duration_minutes"),
                            gbm.y = "species_observed")
summary.gbm(brt_reduced)
gbm.plot(brt_reduced)


# saveRDS(brt_reduced, "cerw_dist_test.RDS")

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

predictor_stack <- terra::rast("predictor_layers/predictor_stack.grd") 
# predictor_stack[[c(2,3,4)]] <- predictor_stack[[c(2,3,4)]]/1000
# predictor_stack[[c(4)]] <- predictor_stack[[c(4)]] #TODO check on this?
# if(crs(predictor_stack) != crs(bcr_albers)) predictor_stack <- project(predictor_stack, bcr_albers)

names(predictor_stack)

preds <- terra::predict(object = predictor_stack, model = brt,
                        ext = st_bbox(bcr_albers),
                        const = data.frame(time_observations_started = 10,
                                           duration_minutes = 60,
                                           effort_radius_m = 1000)
                        )
values(preds) <- boot::inv.logit(values(preds))
preds <- terra::crop(preds, vect(bcr_albers), mask = TRUE)
# plot(preds)

#--------------------------
# Make maps of predictions and validation layers
states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)

ggplot() +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "P(observation)")

atlas_crs <- sf::st_read("./atlas/gwwa/BirdAtlas_6420.kml", layer = "Golden-Winged Warbler") %>%
  st_crs()
atlas_current <- terra::rast("./atlas/gwwa/models/cur.png")[[1]]
ext(atlas_current) <- c(-112.993023, -58.990946, 19.066954, 54.427486)
crs(atlas_current) <- atlas_crs$wkt
atlas_current <- (255 - atlas_current)/255
atlas_current <- terra::project(atlas_current, preds) %>%
  terra::crop(bcr_albers, mask = TRUE)
# plot(atlas_current)

ggplot() +
  geom_spatraster(data = atlas_current) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "P(observation)")


ebd_st <- terra::rast("./status_and_trends/gowwar_abundance_seasonal_breeding_max_2022.tif") %>%
  terra::project(preds) %>%
  terra::crop(vect(bcr_albers), mask = TRUE) %>%
  terra::clamp(upper = 1)
  # `/`(max(.[], na.rm = TRUE))
  # terra::classify(rcl = c(0, 1, Inf))
  
# plot(ebd_st)
ggplot() +
  geom_spatraster(data = ebd_st*2.5) +
  scale_fill_terrain_c() +
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "P(observation)")

ebd_st <- ebd_st * (0.29/0.075)
atlas_current <- atlas_current * (mean(preds[], na.rm = TRUE) / mean(atlas_current[], na.rm = TRUE))

error <- ebd_st - preds
error <- atlas_current - preds
plot(error)
hist(error)
