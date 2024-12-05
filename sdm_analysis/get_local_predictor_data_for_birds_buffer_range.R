# get local data for birds
library("tidyverse")
library("sf")
library("terra")

species_list <- c("bcch", "ewpw", "evgr", "wewa",
                  "nawa", "lowa", "prow",
                  "kewa", "praw", "ytwa",
                  "rugr", "nswo", "veer", "heth", "oven",
                  "woth")
# buffer_scale <- "effort"
# buffer_scale <- 10000


bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326") %>%
  filter(BCR == 28) %>%
  # filter(BCR %in% c(12,13,14,22,23,24,27,28,29,30)) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")
# plot(st_geometry(bcr))

for(species in species_list){
#get all the processed bird data and extract the unique checklists
  sp_start_time <- Sys.time()
  
  ebird_all_obs <- list.files("./ebird", full.names = TRUE) %>%
    `[`(grepl(species, .)) %>%
    `[`(grepl("csv", .)) %>%
    `[`(grepl("subsampled_balanced", .)) %>%
    `[`(grepl("BCR28", .)) %>% #this should hopefully only been one checklist!
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
  
  #process in groups of 500 at a time
  ebird_ss$group <- floor(seq(1, nrow(ebird_ss))/500) + 1
  
  predictor_stack <- terra::rast("./predictor_layers/predictor_stack_BCR28.tif")
  
  scale_range <- c("effort", 100, 500, 1000, 5000, 10000)
  
  for(buffer_scale in scale_range){
    print(paste0("Running script for species ", species, " at scale ", buffer_scale))
    scale_start_time <- Sys.time()
    
    ebird_buff <- sf::st_buffer(ebird_ss, dist = ifelse(buffer_scale == "effort",
                                                        ebird_ss$effort_radius_m,
                                                        as.numeric(buffer_scale)))
    ebird_buff_albers <- ebird_buff %>% st_transform("epsg:5070")
    
      ## GEDI Level 2A product-----------------------------------------------------
    #just grab one of the relative height bands
    ebird_buff_albers$height <- terra::extract(predictor_stack$height, 
                                               vect(ebird_buff_albers),
                                               ID = FALSE,
                                               fun = function(x) mean(x, na.rm = TRUE)
              )$height
    
      boxplot(ebird_buff_albers$height ~ ebird_buff_albers$species_observed)
    
    ebird_buff_albers$open_area<- terra::extract(predictor_stack$open_area, 
                                                 vect(ebird_buff_albers),
                                                 ID = FALSE,
                                                 fun = function(x) mean(x, na.rm = TRUE)
      )$open_area
    
      boxplot(ebird_buff_albers$open_area~ ebird_buff_albers$species_observed)
    
    ##GEDI Level 2B products--------------------------------------------------------
    ebird_buff_albers$fhd_normal <- terra::extract(predictor_stack$fhd_normal, 
                                                   vect(ebird_buff_albers),
                                                   ID = FALSE,
                                                   fun = function(x) mean(x, na.rm = TRUE)
      )$fhd_normal
    
    ebird_buff_albers$understory_ratio <- terra::extract(predictor_stack$understory_ratio, 
                                                         vect(ebird_buff_albers),
                                                         ID = FALSE,
                                                         fun = function(x) mean(x, na.rm = TRUE)
      )$understory_ratio
      boxplot(ebird_buff_albers$understory_ratio ~ ebird_buff_albers$species_observed)
    
    ebird_buff_albers$biomass <- terra::extract(predictor_stack$biomass, 
                                                vect(ebird_buff_albers),
                                                ID = FALSE,
                                                fun = function(x) mean(x, na.rm = TRUE)
      )$biomass
      boxplot(ebird_buff_albers$biomass ~ ebird_buff_albers$species_observed)
    
    
    #-----------------------------------------------------------
    # NLCD cover
    #-----------------------------------------------------------
    nlcd <- rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif")
    nlcd_points <- terra::extract(nlcd, vect(ebird_buff_albers), raw = TRUE) %>%
      as.data.frame() %>%
      group_by(ID) %>%
      summarise(total_cells = n(),
                prop_forest = sum(`NLCD Land Cover Class` %in% c(41,42,43,90))/total_cells,
                prop_decid = sum(`NLCD Land Cover Class` %in% c(41))/total_cells,
                prop_conifer = sum(`NLCD Land Cover Class` %in% c(42))/total_cells,
                prop_grass = sum(`NLCD Land Cover Class` %in% c(71, 81))/total_cells,
                prop_shrub = sum(`NLCD Land Cover Class` %in% c(52))/total_cells,
                prop_water = sum(`NLCD Land Cover Class` %in% c(11))/total_cells,
                prop_dev_light = sum(`NLCD Land Cover Class` %in% c(11))/total_cells,
                prop_dev_heavy = sum(`NLCD Land Cover Class` %in% c(11))/total_cells,
                prop_open = sum(`NLCD Land Cover Class` %in% c(12,21,31,51,52,71:74,81,82,95))/total_cells) %>%
      cbind(ebird_buff_albers)
    
    #--------------------------------------------------------------------------
    #treemap
    #--------------------------------------------------------------------------
    
    prop_oak <- terra::rast("./predictor_layers/prop_oak_fia.tiff")
    prop_spruce <- terra::rast("./predictor_layers/prop_spruce_fia.tiff")
    
    ebird_buff_treemap <- ebird_buff %>% sf::st_transform(crs(prop_oak))
    
    oak_extract <- vector("numeric", length = 0L)
    for(i in 1:length(unique(ebird_buff_treemap$group))){
      oak_extract_vals <- terra::extract(prop_oak, vect(ebird_buff_treemap[ebird_buff_treemap$group == i, ]), 
                                  fun = mean, na.rm = TRUE, ID = FALSE)
      oak_extract <- c(oak_extract, oak_extract_vals)
    }
    oak_extract <- unname(unlist(oak_extract))
    
    spruce_extract <- vector("numeric", length = 0L)
    for(i in 1:length(unique(ebird_buff_treemap$group))){
      spruce_extract_vals <- terra::extract(prop_spruce, vect(ebird_buff_treemap[ebird_buff_treemap$group == i, ]), 
                                            fun = mean, na.rm = TRUE, ID = FALSE)
      spruce_extract <- c(spruce_extract, spruce_extract_vals)
    }
    
    spruce_extract <- unname(unlist(spruce_extract))
    ebird_buff_treemap <- cbind(ebird_buff_treemap, oak_extract, spruce_extract)
    ebird_buff_treemap <- rename(ebird_buff_treemap,
                                 prop_oak = oak_extract,
                                 prop_spruce = spruce_extract) %>%
      st_drop_geometry()
    ebird_buff_albers <- st_drop_geometry(ebird_buff_albers)
    
    
    combined <- dplyr::select(nlcd_points, prop_forest, prop_decid, prop_conifer, prop_shrub, prop_grass, prop_water, 
                              prop_dev_light, prop_dev_heavy, prop_open, checklist_id) %>%
      left_join(dplyr::select(ebird_buff_albers, height, open_area, fhd_normal, understory_ratio, biomass, checklist_id), by = "checklist_id") %>%
      left_join(dplyr::select(ebird_buff_treemap, prop_oak, prop_spruce, checklist_id), by = "checklist_id") %>%
      right_join(dplyr::select(ebird_all_obs, species, checklist_id, time_observations_started, duration_minutes, effort_radius_m, species_observed), by = "checklist_id") %>%
      ungroup() %>%
      st_as_sf() %>%
      dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                    lat = sf::st_coordinates(.)[,2]) %>%
      st_drop_geometry()
    combined %>%
      write_csv(., paste0("./environment_vars/", species, "_combined_data_bcr28_local_", buffer_scale,".csv"))
  
    scale_elapsed <- system.time() - scale_start_time
    
    print(paste0("Scale ", buffer_scale, " done! It took ", scale_elapsed))
    
  }
  
  species_elapsed <- system.time() - species_start_time
  message(paste0("Species ", species, " done! It took ", species_elapsed))
}

