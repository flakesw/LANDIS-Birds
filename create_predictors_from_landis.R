#bring in LANDIS layers to generate new prediction surfaces
library("terra")
library("sf")
library("tidyverse")

model_name <- "HighTHighV BAU"
model_dir <- paste0("D:/SApps LANDIS/Model templates/", model_name, "/")
input_dir <- "D:/SApps LANDIS/Inputs/"
year <- 60

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
                             to = total_biomass$biomass) / 100 #convert to tonnes ha-1
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
decid_rast <- decid_rast
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
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$tmax,
                 others = NA)
  tmin_rast <- ecoregions %>%
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$tmin,
                 others = NA)
  precip_rast <- ecoregions %>%
    terra::subst(from = summer_clim$EcoregionIndex + 2, #add 2 to line up with ecoregions; different for each model!
                 to = summer_clim$precip,
                 others = NA)
}


#TODO estimate height, FHD!

#stack all the predictors
if(year == 0){
  pred_stack <- c(biomass2, 
                  under2, 
                  decid2,
                  # predictor_stack2$aet,
                  # predictor_stack2$pet,
                  # predictor_stack2$def,
                  # predictor_stack2$soil,
                  # predictor_stack2$vpd,
                  # predictor_stack2$pdsi,
                  predictor_stack2$tmmn/10,
                  predictor_stack2$tmmx/10,
                  predictor_stack2$pr,
                  predictor_stack2$slope,
                  predictor_stack2$chili,
                  predictor_stack2$tpi)
} else{
  pred_stack <- c(biomass2, 
                  under2, 
                  decid2,
                  # predictor_stack2$aet,
                  # predictor_stack2$pet,
                  # predictor_stack2$def,
                  # predictor_stack2$soil,
                  # predictor_stack2$vpd,
                  # predictor_stack2$pdsi,
                  tmin_rast,
                  tmax_rast,
                  precip_rast,
                  predictor_stack2$slope,
                  predictor_stack2$chili,
                  predictor_stack2$tpi)
}


names(pred_stack) <- c("biomass", "understory_ratio", "prop_decid",
                          "tmmn", "tmmx", "pr", "slope", "chili", "tpi")

writeRaster(pred_stack, 
            filename = paste0("./landis_predictor_layers/pred_stack_", 
                                          model_name, "_", year, ".tif"),
            overwrite = TRUE)
