### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
library("terra")
library("sf")
#--------------------------



###### BRT

combined <- read.csv("cerw_combined_data.csv")
combined <- combined[complete.cases(combined), ] %>%
  rename(biomass = MU,
         understory_ratio = understory_proportion)
combined$grass <- combined$prop_grass

# combined$species_observed <- ifelse(as.character(combined$species_observed) == "FOUND", TRUE, FALSE)

brt <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = c("aet",
                                 "tmmx", "vpd", #"pr",
                                 "tpi", "chili", "height",
                                 "prop_decid",  "biomass",
                                 "understory_ratio", "grass",
                                 "prop_oak", "prop_forest", 
                                 "time_observations_started", "duration_minutes"),
                       gbm.y = "species_observed",
                       interaction.depth = 2)
print(brt)
summary.gbm(brt)
gbm.interactions(brt)
gbm.plot(brt)
gbm.perspec(brt, x = 12, y = 13)
gbm.perspec(brt, x = 12, y = 13)
gbm.perf(brt)
brt$self.statistics$discrimination
brt$cv.statistics$discrimination.mean
roc(combined$species_observed, boot::inv.logit(predict(brt)),
    plot = TRUE)

brt_simp <- dismo::gbm.simplify(brt)
plot(brt_simp$deviance.summary$mean ~ c(1:length(brt_simp$deviance.summary$mean)))
brt2<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = brt_simp$pred.list[[1]],
                       gbm.y = "species_observed")

saveRDS(brt2, "cerw_dist_model.RDS")


gbm.plot(brt2,
         n.plots = 9,
         common.scale = TRUE,
         show.contrib = TRUE,
         plot.layout = c(3,3))
gbm.perspec(brt, x = 7, y = 6)


brt_veg <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                           gbm.x = c("deciduous", 
                                     "biomass", "understory_ratio"),
                           gbm.y = "species_observed")
summary.gbm(brt_veg)
gbm.plot(brt_veg)
saveRDS(brt_veg, "cerw_dist_veg.RDS")

brt_clim <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("tmmn", "pr", "tmmx", "tpi", "chili"),
                            gbm.y = "species_observed")
summary.gbm(brt_clim)
gbm.plot(brt_clim)
saveRDS(brt_clim, "cerw_dist_clim.RDS")

brt_reduced <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c(#"aet", #"pet",
                                      "tmmx", "tmmn", "pr",
                                      "tpi", "chili",
                                      "prop_decid", "biomass",
                                      "prop_oak",
                                      "understory_ratio",
                                      "time_observations_started", "duration_minutes"),
                            gbm.y = "species_observed")
summary.gbm(brt_reduced)
gbm.plot(brt_reduced)


saveRDS(brt_reduced, "cerw_dist_test.RDS")

## make maps ------------------------------------
# 
# bcr_conical <- bcr %>% st_transform(st_crs("PROJCS[\"NAD_1983_Albers\",
#     GEOGCS[\"NAD83\",
#         DATUM[\"North_American_Datum_1983\",
#             SPHEROID[\"GRS 1980\",6378137,298.257222101,
#                 AUTHORITY[\"EPSG\",\"7019\"]],
#             TOWGS84[0,0,0,0,0,0,0],
#             AUTHORITY[\"EPSG\",\"6269\"]],
#         PRIMEM[\"Greenwich\",0,
#             AUTHORITY[\"EPSG\",\"8901\"]],
#         UNIT[\"degree\",0.0174532925199433,
#             AUTHORITY[\"EPSG\",\"9108\"]],
#         AUTHORITY[\"EPSG\",\"4269\"]],
#     PROJECTION[\"Albers_Conic_Equal_Area\"],
#     PARAMETER[\"standard_parallel_1\",29.5],
#     PARAMETER[\"standard_parallel_2\",45.5],
#     PARAMETER[\"latitude_of_center\",23],
#     PARAMETER[\"longitude_of_center\",-96],
#     PARAMETER[\"false_easting\",0],
#     PARAMETER[\"false_northing\",0],
#     UNIT[\"meters\",1]]"))

ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_albers <- bcr %>% 
  st_transform(crs(ecoregions)) #%>%
  # st_crop(ecoregions) #st_transform("EPSG:5070")

predictor_stack <- terra::rast("predictor_layers/predictor_stack.grd") %>%
  project(crs(bcr_albers))
names(predictor_stack)

preds <- terra::predict(object = predictor_stack, model = brt,
                        ext = st_bbox(bcr_albers),
                        const = data.frame(time_observations_started = 700,
                                           duration_minutes = 60)
                        )
values(preds) <- boot::inv.logit(values(preds))
preds <- crop(preds, vect(bcr_albers), mask = TRUE)
plot(preds)

states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)

library("tidyterra")
ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "P(observation)")

