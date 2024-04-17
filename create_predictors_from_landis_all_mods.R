#bring in LANDIS layers to generate new prediction surfaces
library("terra")
library("sf")
library("tidyverse")

#"HighTHighV BAU extra Rx"

model_names <- c("LowTLowV BAU", "LowTLowV BAU extra Rx", "LowTLowV BAU nofire",
                 "HighTHighV BAU nofire", "HighTHighV BAU extra Rx")

input_dir <- "C:/Users/swflake/Documents/SApps LANDIS/Inputs/"
years <- c(10,20,30,40,50, 60)

for(model in model_names){
  for(year in years){
    make_predictor_stack_from_landis(model_name = model, input_dir = input_dir, year)
  }
}

make_predictor_stack_from_landis <- function(model_name, input_dir, year)
{

  model_dir <- paste0("C:/Users/swflake/Documents/SApps LANDIS/Model templates/", model_name, "/")
  
  comm_output <- read.csv(paste0(model_dir, "/community-input-file-", year,".csv"))
  comm_map <- terra::rast(paste0(model_dir, "/output-community-", year, ".img"))
  plot(comm_map)
  
  climate_future <- read.csv(paste0(model_dir, "Climate-future-input-log.csv"))
  ecoregions <- terra::rast(paste0(input_dir, "Basic_inputs/Ecos11_NCLD.tif"))
  plot(ecoregions)
  
  predictor_stack <- terra::rast("./landis_predictor_layers/pred_stack_study_area.tif")
  predictor_stack2 <- predictor_stack 
  
  comm_rast <- terra::ifel(comm_map %in% comm_output$MapCode, comm_map, NA)
  
  #biomass
  comm_output[comm_output[] == 0] <- NA
  
  total_biomass <- comm_output %>%
    group_by(MapCode) %>%
    summarise(biomass = sum(CohortBiomass))
  biomass_rast <- comm_rast
  biomass_rast <- terra::subst(x = biomass_rast, 
                               from = total_biomass$MapCode, 
                               to = total_biomass$biomass) / 100 #convert to tonnes ha-1
  plot(biomass_rast)
  biomass2 <- ecoregions
  biomass2[] <- biomass_rast[]
  # saveRDS(biomass_rast, "biomass_raster_60.RDS")
  
  # lai <- terra::rast(paste0(model_dir, "/NECN/LAI-5.img"))
  
  #understory_ratio 
  # under_ratio <- comm_output %>%
  #   mutate(understory = as.integer(CohortAge <= 30)) %>%
  #   mutate(understory = ifelse(SpeciesName %in% c("CornFlor", "AmelArbo", "AcerPens", 
  #                                                 "HaleDipt", "OxydArbo", "SassAlid"),
  #                              1, understory)) %>%
  #   mutate(understory = ifelse(understory == 0, 0.5, understory)) %>%
  #   mutate(understory_biomass = understory*CohortBiomass) %>%
  #   group_by(MapCode) %>%
  #   summarise(understory_ratio = sum(understory_biomass) / sum(CohortBiomass)) 
  # under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
  # under_rast <-  terra::subst(x = under_rast, 
  #                             from = under_ratio$MapCode, 
  #                             to = under_ratio$understory_ratio)
  # plot(under_rast)
  # under2 <- ecoregions
  # under2[] <- under_rast[]
  
  #mean_age
  mean_age <- comm_output %>%
    group_by(MapCode) %>%
    summarise(mean_age = mean(CohortAge)) 
  age_rast <- comm_rast
  age_rast <-  terra::subst(x = age_rast, 
                              from = mean_age$MapCode, 
                              to = mean_age$mean_age)
  plot(age_rast)
  age2 <- ecoregions
  age2[] <- age_rast[]
  under2 <- age2/200
  
  # saveRDS(under_rast, "understory_raster_60.RDS")
  
  
  #foliage height diversity
  # fhd <- comm_output %>%
  #   group_by(MapCode) %>%
  #   summarize(age_diversity = )
  #   mutate(understory = as.integer(CohortAge <= 30)) %>%
  #   mutate(understory = ifelse(SpeciesName %in% c("CornFlor", "AmelArbo", "AcerPens", 
  #                                                 "HaleDipt", "OxydArbo", "SassAlid"),
  #                              1, understory)) %>%
  #   mutate(understory = ifelse(understory == 0, 0.5, understory)) %>%
  #   mutate(understory_biomass = understory*CohortBiomass) %>%
  #   group_by(MapCode) %>%
  #   summarise(understory_ratio = sum(understory_biomass) / sum(CohortBiomass)) 
  # under_rast <- terra::ifel(comm_map %in% under_ratio$MapCode, comm_map, NA)
  # under_rast <-  terra::subst(x = under_rast, 
  #                             from = under_ratio$MapCode, 
  #                             to = under_ratio$understory_ratio)
  # plot(under_rast)
  # under2 <- ecoregions
  # under2[] <- under_rast[]
  
  
  #deciduousness
  #areas dominated by trees generally greater than 5 meters tall, 
  #and greater than 20% of total vegetation cover. 
  #More than 75% of the tree species shed foliage simultaneously in response to seasonal change.
  #TODO convert LAI to total veg coverhttp://127.0.0.1:35167/graphics/plot_zoom_png?width=826&height=886
  
  decid_all <- comm_output
  decid_all$decid = !grepl("Pinu|Tsug", comm_output$SpeciesName)
  decid_sum <- decid_all %>%
    group_by(MapCode, decid) %>%
    summarise(biomass = sum(CohortBiomass))
  deciduous <- decid_sum %>% 
    group_by(MapCode) %>%
    tidyr::pivot_wider(names_from = decid, values_from = biomass, values_fill = 0) %>%
    mutate(decid_prop = `TRUE`/(`TRUE`+ `FALSE`)) %>%
    mutate(decid = ifelse(decid_prop > 0.75, 1, 0)) %>%
    dplyr::select(MapCode, decid)
  decid_rast <- comm_rast
  decid_rast <-  terra::subst(x = decid_rast, 
                              from = deciduous$MapCode, 
                              to = deciduous$decid)
  # decid_rast[is.na(decid_rast[])] <- 0
  decid_smooth <- terra::focal(decid_rast, 5, fun = mean, na.rm = TRUE)
  
  decid2 <- ecoregions
  decid2[] <- decid_smooth[]
  decid2 <-  mask(decid2, ecoregions, maskvalues = 1)
  
  #climate variables
  if(year > 0){
    summer_clim <- climate_future %>%
      filter(Year > 2005 + year - 10 & Year <= 2005 + year, 
             Timestep > 182 & Timestep < 273) %>%
      group_by(EcoregionIndex) %>%
      summarise(tmax = median(max_airtemp)*10,#rescale to match GEE output
                tmin = median(min_airtemp)*10,
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
                    predictor_stack2$tmmn,
                    predictor_stack2$tmmx,
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
}
