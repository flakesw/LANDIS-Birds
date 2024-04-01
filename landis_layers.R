#bring in LANDIS layers to generate new prediction surfaces
library("terra")
library("sf")
library("tidyverse")

model_name <- "LowTLowV BAU"
model_dir <- paste0("D:/SApps LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 0

comm_output <- read.csv(paste0(model_dir, "/community-input-file-", year,".csv"))
comm_map <- terra::rast(paste0(model_dir, "/output-community-", year, ".img"))
plot(comm_map)

climate_future <- read.csv(paste0(model_dir, "Climate-future-input-log.csv"))
ecoregions <- terra::rast(paste0(input_dir, "Basic_inputs/Ecos11_NCLD.tif"))
plot(ecoregions)

predictor_stack <- terra::rast("predictor_layers/predictor_stack.grd")
predictor_stack2 <- predictor_stack %>%
  project(ecoregions, mask = TRUE, method = "near", align = FALSE) %>%
  mask(ecoregions, maskvalues = 1)
names(predictor_stack)

#biomass
comm_output[comm_output[] == 0] <- NA

total_biomass <- comm_output %>%
  group_by(MapCode) %>%
  summarise(biomass = sum(CohortBiomass))
biomass_rast <- terra::ifel(comm_map %in% total_biomass$MapCode, comm_map, 0)
biomass_rast <- terra::subst(x = biomass_rast, 
                             from = total_biomass$MapCode, 
                             to = total_biomass$biomass) / 100
plot(biomass_rast)
biomass2 <- ecoregions
biomass2[] <- biomass_rast[]
# saveRDS(biomass_rast, "biomass_raster_60.RDS")

#understory_ratio 
under_ratio <- comm_output %>%
  mutate(understory = CohortAge <= 20) %>%
  group_by(MapCode, understory) %>%
  summarise(biomass = sum(CohortBiomass)) %>%
  tidyr::pivot_wider(names_from = understory, values_from = biomass, values_fill = 0) %>%
  mutate(understory_ratio = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  dplyr::select(MapCode, understory_ratio)
under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
under_rast <-  terra::subst(x = under_rast, 
                            from = under_ratio$MapCode, 
                            to = under_ratio$understory_ratio)
plot(under_rast)
under2 <- ecoregions
under2[] <- under_rast[]

# saveRDS(under_rast, "understory_raster_60.RDS")


#deciduousness
decid_all <- comm_output
decid_all$decid = !grepl("Pinu|Tsug", comm_output$SpeciesName)
decid_sum <- decid_all %>%
  group_by(MapCode, decid) %>%
  summarise(biomass = sum(CohortBiomass))
deciduous <- decid_sum %>% 
  group_by(MapCode) %>%
  tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
  mutate(decid = `TRUE`/(`TRUE`+ `FALSE`)) %>%
  dplyr::select(MapCode, decid)
decid_rast <- terra::ifel(comm_map %in% deciduous$MapCode, comm_map, NA)
decid_rast <-  terra::subst(x = decid_rast, 
                            from = deciduous$MapCode, 
                            to = deciduous$decid)
decid_rast <- decid_rast*100
plot(decid_rast)
decid2 <- ecoregions
decid2[] <- decid_rast[]

#climate variables
if(year > 0){
  summer_clim <- climate_future %>%
    filter(Year > 2005 + year - 10 & Year <= 2005 + year, 
           Timestep > 182 & Timestep < 273) %>%
    group_by(EcoregionIndex) %>%
    summarise(tmax = median(max_airtemp)*100,#rescale to match GEE output
              tmin = median(min_airtemp)*100,
              precip = median(ppt)*1000)
  tmax_rast <- ecoregions %>%
    terra::subst(from = august_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = august_clim$tmax,
                 others = NA)
  tmin_rast <- ecoregions %>%
    terra::subst(from = august_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = august_clim$tmin,
                 others = NA)
  precip_rast <- ecoregions %>%
    terra::subst(from = august_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = august_clim$precip,
                 others = NA)
}



#stack all the predictors
pred_stack_60 <- c(biomass2, under2, decid2,
                      tmax_rast, tmin_rast, precip_rast, 
                   predictor_stack2$chili,
                   predictor_stack2$tpi)
names(pred_stack_60) <- c("biomass", "understory_ratio", "prop_decid",
                          "tmmx", "tmmn", "pr", "chili", "tpi")



#full model
library("dismo")

cerw_model <- readRDS("cerw_dist_test.RDS")
cerw_model$var.names

preds_bak <- preds

preds <- terra::predict(object = pred_stack_60, 
                        model = cerw_model,
                        # ext = st_bbox(bcr_albers),
                        const = data.frame(time_observations_started = 700,
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

