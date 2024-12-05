### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
#--------------------------

combined <- read.csv("gwwa_combined_data.csv")
combined <- combined[complete.cases(combined), ]

boxplot(height ~ species_observed, data = combined)
boxplot(deciduous ~ species_observed, data = combined)
boxplot(MU ~ species_observed, data = combined)
boxplot(understory_proportion ~ species_observed, data = combined)

mod <- glm(species_observed ~ soil + aet + def + tmmx + slope + tpi + chili + height*deciduous, data = combined,
           family = binomial(link = "logit"))

summary(mod)
plot(effects::allEffects(mod))

MuMIn::r.squaredGLMM(mod)
library("pROC")
roc <- roc(combined$species_observed, boot::inv.logit(predict(mod)))
pROC::plot.roc(roc)
pROC::auc(roc)

combined$species_observed <- as.factor(combined$species_observed)
levels(combined$species_observed) <- c("not_Found", "Found")

library("caret")
ctrl <- caret::trainControl(method="cv", 
                            number=10, 
                            savePredictions="all",
                            summaryFunction = twoClassSummary,
                            classProbs=TRUE)
glmFit <- caret::train(
  mod$formula,
  data = combined,
  method = "glm",
  preProc = c("center", "scale"),
  tuneLength = 15,
  trControl = ctrl,
  metric = "ROC"
)

summary(glmFit)




###### BRT

combined <- read.csv("gwwa_combined_data.csv")
combined <- combined[complete.cases(combined), ] %>%
  rename(biomass = MU,
         understory_ratio = understory_proportion)
combined$grass_shrub <- combined$grass + combined$scrub

# combined$species_observed <- ifelse(as.character(combined$species_observed) == "FOUND", TRUE, FALSE)

brt <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = c("soil", "aet", "def", 
                                 "pdsi", "tmmx", "vpd", "pr",
                                 "tpi", "chili", "height",
                                 "deciduous", "fhd_normal", #"biomass",
                                 # "understory_ratio", "grass_shrub",
                                 "time_observations_started", "duration_minutes"),
                       gbm.y = "species_observed",
                       interaction.depth = 2,
                       # var.monotone = c("decid"),
                       train.fraction = 0.7)
print(brt)
summary.gbm(brt)
gbm.interactions(brt)
gbm.plot(brt)
gbm.perspec(brt, x = 13, y = 12)
gbm.perf(brt)
brt$self.statistics$discrimination
brt$cv.statistics$discrimination.mean
roc(combined$species_observed, boot::inv.logit(predict(brt)),
    plot = TRUE)

brt_simp <- dismo::gbm.simplify(brt)
brt2<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                       gbm.x = brt_simp$pred.list[[9]],
                       gbm.y = "species_observed")
gbm.plot(brt2,
         # n.plots = 9,
         common.scale = TRUE,
         show.contrib = TRUE,
         plot.layout = c(3,3))
gbm.perspec(brt, x = 7, y = 6)


brt_veg <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                           gbm.x = c("height", "deciduous", "fhd_normal", 
                                     "biomass", "understory_ratio"),
                           gbm.y = "species_observed")
summary.gbm(brt_veg)
gbm.plot(brt_veg)

brt_clim <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("soil", "aet", "tmmx", "tpi", "chili"),
                            gbm.y = "species_observed")
summary.gbm(brt_clim)
gbm.plot(brt_clim)


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

bcr_albers <- bcr %>% st_transform("EPSG:5070")


predictor_stack <- terra::rast("predictor_layers/predictor_stack.grd") 
names(predictor_stack)

preds <- terra::predict(object = predictor_stack, model = brt2,
                         ext = st_bbox(bcr_albers))
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
  scale_fill_terrain_c(alpha = 0.7)
