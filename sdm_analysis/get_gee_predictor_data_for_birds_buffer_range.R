###################
# Get predictor variables from GEE
# Using GEE for climate variables and remote sensing variables
#this is the newest version of script as of 12/3/2024
#------------------------------------------------------
library("tidyverse")
library("rgee")
library("sf")
library("dismo")
library("stars")
library("starsExtra")
library("terra")
#ee_install()
#ee_install_upgrade()
ee_Initialize(user = "swflake@ncsu.edu", drive = TRUE)

species_list <- c("bcch", "ewpw", "evgr", "wewa",
                 "nawa", "lowa", "prow",
                  "kewa", "praw", "ytwa",
                  "rugr", "nswo", "veer", "heth", "oven",
                   "woth")
# buffer_scale <- "effort"
# buffer_scale <- 10000


bufferBy100 = function(feature) {
  feature$buffer(feature)
}

bufferFeatureByEffort = function(f) {
  f = ee$Feature(f)
  if(buffer_scale == "effort"){
    buffer_size = f$get('effort_radius_m')
  } else{
    buffer_size = as.numeric(buffer_scale)
  }
  return(f$buffer(buffer_size))
}

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326") %>%
  filter(BCR == 28) %>%
  # filter(BCR %in% c(12,13,14,22,23,24,27,28,29,30)) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")
# plot(st_geometry(bcr))

bcr_ee <- sf_as_ee(bcr)
# bcr_ee$getInfo()
# Map$addLayer(bcr_ee)

#get all the processed bird data and extract the unique checklists
for(species in species_list){
  sp_start_time <- Sys.time()
  
  ebird_all_obs <- list.files("./ebird", full.names = TRUE) %>%
    `[`(grepl(species, .)) %>%
    `[`(grepl("csv", .)) %>%
    `[`(grepl("subsampled_balanced", .)) %>%
    `[`(grepl("BCR28", .)) %>%
    read_csv() %>% 
    dplyr::mutate(month = substr(observation_date, 6, 7)) %>%
    dplyr::filter((protocol_type == "Stationary" | effort_distance_km  < 5) &
                    month %in% c("06","07","08")) %>%
    sf::st_as_sf(coords = c("longitude", "latitude")) %>%
    sf::st_set_crs("epsg:4326") %>%
    mutate(effort_radius_m = ifelse(is.na(effort_distance_km), 0.5, 
                                    ifelse(effort_distance_km < 0.5 , 0.5, effort_distance_km)) * 1000 /2) 
  ebird_unique_checklists <- ebird_all_obs %>%
    dplyr::distinct(checklist_id, .keep_all = TRUE)
  
  ebird_ss <- ebird_unique_checklists %>%
    dplyr::select(checklist_id, geometry, effort_radius_m, species_observed, 
                  time_observations_started, duration_minutes, observation_date)
  
  scale_range <- c("effort", 100, 500, 1000, 5000, 10000)
  # scale_range <- c(100, 500, 1000, 5000, 10000)
  
  for(buffer_scale in scale_range){
  
    print(paste0("Running script for species ", species, " at scale ", buffer_scale))
    scale_start_time <- Sys.time()
    
    ebird_buff <- sf::st_buffer(ebird_ss, dist = ifelse(buffer_scale == "effort",
                                                        ebird_ss$effort_radius_m,
                                                        as.numeric(buffer_scale)))
    ebird_buff_albers <- ebird_buff %>% st_transform("epsg:5070")
    
    points_ee <-  ebird_ss %>% #ebird_ss_group %>% 
      dplyr::select(checklist_id, effort_radius_m) %>%
      sf_as_ee()
    
    # ee_print(points_ee)
    
    # make a GEE polygon buffer about the point
    aoi = points_ee$map(bufferFeatureByEffort)
    
    clim_buffer <- ifelse(is.numeric(buffer_scale), 
                          ifelse(buffer_scale > 2000, as.numeric(buffer_scale), 2000), 
                          2000) #minimum climate buffer of 2km
    aoi_clim = points_ee$map(function(poly) poly$buffer(clim_buffer))
    # aoi2 = ee$Feature(aoi$first())
    # ee_print(aoi2)
    # Map$addLayer(aoi2)
    
    #get climate data -------------------------------------------------
    
    #TODO reduce number of variables before exporting -- do filtering server-side instead of downloading
    # all the climate data for each point
    #TODO get year-of-observation weather (or week of observation?) as covariate
    
    terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    
    monthstart <- 6
    monthend <- 8
    years <- ee$List$sequence(2013, 2023)
    lag_years <- ee$Number(10) #TODO try different time lags!
    
    make_preceding_climate <- function(x) {
      endyear <- ee$Number(x)
      startyear <- endyear$subtract(lag_years)
      
      tc_filtered <- terraclimate$filter(ee$Filter$calendarRange(startyear, endyear, "year"))$
        filter(ee$Filter$calendarRange(ee$Number(monthstart), ee$Number(monthend), "month"))$
        filterBounds(aoi)
      
      tc_mean <- tc_filtered$mean()
      return(tc_mean)
      
    }
    
    tc_lagged <- years$map(ee_utils_pyfunc(make_preceding_climate)) %>%
      ee$ImageCollection()
    
    print("Extracting climate for points")
    mean_clim <- ee_extract(tc_lagged$select(c("pet", "aet", "tmmn", "tmmx",
                                               "soil", "pdsi", "vpd", "pr")), 
                            aoi_clim, 
                            fun = ee$Reducer$mean(), #TODO add more variables -- max or min for example
                            scale = 4000,
                            quiet = TRUE,
                            via = "drive")
    # tileScale = 4)
    
    climate_vars <- mean_clim %>%
      left_join(dplyr::select(ebird_buff_albers, checklist_id, observation_date) %>% st_drop_geometry(), by = "checklist_id") %>%
      pivot_longer(cols = !c(checklist_id, effort_radius_m, observation_date), 
                   names_to = c("year", "var"),
                   values_to = "value",
                   names_pattern = "X_?(.*)_(.*)") %>%
      mutate(obs_year = as.numeric(substr(observation_date, 1, 4)),
             clim_year = as.numeric(year) + 2013) %>%
      filter(clim_year == obs_year) %>%
      dplyr::select(checklist_id, var, value) %>%
      pivot_wider(names_from = var,
                  values_from = value) %>%
      mutate(def = pet - aet)
    
    
    #-------------------------------------------------------------------------------
    # MODIS NPP
    #-------------------------------------------------------------------------------
    
    #use similar method as climate data, above, to get time-lagged NPP
    
    modis <- ee$ImageCollection("MODIS/061/MOD17A3HGF")
    
    years <- ee$List$sequence(2013, 2023)
    lag_years <- ee$Number(10) #TODO try different time lags!
    
    make_preceding_npp <- function(x) {
      endyear <- ee$Number(x)
      startyear <- endyear$subtract(lag_years)
      
      modis_filtered <- modis$filter(ee$Filter$calendarRange(startyear, endyear, "year"))$
        filterBounds(aoi)
      
      modis_mean <- modis_filtered$mean()
      return(modis_mean)
      
    }
    
    modis_lagged <- years$map(ee_utils_pyfunc(make_preceding_npp)) %>%
      ee$ImageCollection()
    
    #units are kg*C/m^2 at scale 0.0001
    #convert to gC/m2 by multiplying by 1000*0.0001 = divide by 10
    print("Extracting NPP for points")
    mean_npp <- ee_extract(modis_lagged$select(c("Npp")), 
                            aoi_clim, 
                            fun = ee$Reducer$mean(), #TODO add more variables -- max or min for example
                            scale = 4000,
                            quiet = TRUE,
                            via = "drive")
    
    npp_vals <- mean_npp %>%
      left_join(dplyr::select(ebird_buff_albers, checklist_id, observation_date) %>% st_drop_geometry(), 
                by = "checklist_id") %>%
      pivot_longer(cols = !c(checklist_id, effort_radius_m, observation_date), 
                   names_to = c("year", "var"),
                   values_to = "value",
                   names_pattern = "X_?(.*)_(.*)") %>%
      mutate(obs_year = as.numeric(substr(observation_date, 1, 4)),
             clim_year = as.numeric(year) + 2013) %>%
      filter(clim_year == obs_year) %>%
      dplyr::select(checklist_id, var, value) %>%
      pivot_wider(names_from = var,
                  values_from = value)
    
    
    # topography -------------------------------------------------------------------
    
    #A digital elevation model.
    dem = ee$Image('NASA/NASADEM_HGT/001')$select('elevation')
    
    #Calculate slope. Units are degrees, range is [0,90).
    slope = ee$Terrain$slope(dem)
    
    #Calculate aspect. Units are degrees where 0=N, 90=E, 180=S, 270=W.
    # aspect = ee$Terrain$aspect(dem)
    print("Extracting slope for points")
    slope_extract <- ee_extract(slope, 
                                aoi, 
                                fun = ee$Reducer$mean(), 
                                scale = 100,
                                quiet = TRUE,
                                via = "drive",
                                tileScale = 16)%>%
      dplyr::select(c("checklist_id", "slope")) 
    
    #topographic position index
    # ALOS-derived mTPI ranging from negative (valleys) to positive (ridges) values
    tpi = ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI")
    print("Extracting TPI for points")
    tpi_extract <- ee_extract(tpi, 
                              aoi, 
                              fun = ee$Reducer$mean(), 
                              scale = 100,
                              quiet = TRUE,
                              via = "drive",
                              tileScale = 16) %>%
      dplyr::select(c("checklist_id", "AVE")) %>%
      mutate(tpi = AVE)
    
    #heat load index
    chili = ee$Image("CSP/ERGo/1_0/Global/ALOS_CHILI")
    print("Extracting CHILI for points")
    chili_extract = ee_extract(chili, 
                               aoi, 
                               fun = ee$Reducer$mean(), 
                               scale = 100,
                               quiet = TRUE,
                               via = "drive",
                               tileScale = 16) %>%
      dplyr::select(c("checklist_id", "constant")) %>%
      mutate(chili = constant) 
    
    
    
    #combine data ------------------------------------------------------------------
    combined <- left_join(dplyr::select(ebird_buff_albers %>% st_drop_geometry(), 
                                        checklist_id),
                          climate_vars, by = "checklist_id") %>%
      left_join(npp_vals, by = "checklist_id") %>%
      left_join(slope_extract, by = "checklist_id") %>%
      left_join(dplyr::select(tpi_extract, checklist_id, tpi), by = "checklist_id") %>%
      left_join(dplyr::select(chili_extract, checklist_id, chili), by = "checklist_id") %>%
      right_join(dplyr::select(ebird_all_obs, species, checklist_id, time_observations_started, duration_minutes, effort_radius_m, species_observed), by = "checklist_id") %>%
      ungroup() %>%
      st_as_sf() %>%
      dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                    lat = sf::st_coordinates(.)[,2]) %>%
      st_drop_geometry()
    combined %>%
      write_csv(., paste0("./environment_vars/", species, "_combined_data_bcr28_gee_", buffer_scale,".csv"))
    
    scale_elapsed <- system.time() - scale_start_time
      
    print(paste0("Scale ", buffer_scale, " done! It took ", scale_elapsed))
    
  }
  
  species_elapsed <- system.time() - species_start_time
  message(paste0("Species ", species, " done! It took ", species_elapsed))
}
