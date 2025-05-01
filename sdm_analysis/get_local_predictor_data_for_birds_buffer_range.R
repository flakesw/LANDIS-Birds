# get local data for birds
library("tidyverse")
library("sf")
library("terra")

#------------------

sample_and_average <- function(x, max_shots = 100){
  x <- x[!is.na(x)]
  x[x<0] <- 0
  n_shots <- min(max_shots, length(x)) #max of max_shots shots per point
  selected <- x[sample(n_shots)]
  mean <- mean(selected, na.rm = TRUE)
  return(mean)
}

extract_sample <- function(x, max_shots = 100){
  x <- x[!is.na(x)]
  n_shots <- min(max_shots, length(x)) #max of max_shots shots per point
  selected <- list(x[sample(n_shots)])
  return(selected)
}

#-------------------------

# species_list <- c(#"cerw", "gwwa", "recr", "swwa", "bwwa", 
#                   "blbw",
#                   "bcch", "ewpw", "evgr", "wewa",
#                   "nawa", "lowa", "prow",
#                   "kewa", "praw", "ytwa",
#                   "rugr", "nswo", "veer", "heth", "oven",
#                   "woth", "alfl", "bhnu", "osfl")

species_list <- c("acfl")

# "acfl" -- needs scale 1000 and 5000

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
  species_start_time <- Sys.time()
  
  ebird_all_obs <- list.files("./ebird", full.names = TRUE) %>%
    `[`(grepl(species, .)) %>%
    `[`(grepl("csv", .)) %>%
    `[`(grepl("subsampled_balanced", .)) %>%
    `[`(grepl("BCR28", .)) %>% #this should hopefully only been one checklist!
    read_csv(show_col_types = FALSE) %>% 
    dplyr::mutate(month = substr(observation_date, 6, 7)) %>%
    dplyr::filter((protocol_type == "Stationary" | effort_distance_km  < 5) &
                    month %in% c("06","07","08")) %>%
    sf::st_as_sf(coords = c("longitude", "latitude")) %>%
    sf::st_set_crs("epsg:4326") %>%
    mutate(effort_radius_m = ifelse(is.na(effort_distance_km), 0.5, 
                                    ifelse(effort_distance_km < 0.5 , 0.5, effort_distance_km)) * 1000 /2) 
  ebird_unique_checklists <- ebird_all_obs %>% #this only helps if you're doing multiple species at one time, which we're not actually doing at the moment. But it doesn't hurt anything either.
    dplyr::distinct(checklist_id, .keep_all = TRUE)
  
  
  ebird_ss <- ebird_unique_checklists %>%
    dplyr::select(checklist_id, geometry, effort_radius_m, species_observed, 
                  time_observations_started, duration_minutes, observation_date)
  
  #process in groups of 500 at a time
  ebird_ss$group <- floor(seq(0, (nrow(ebird_ss)-1))/500) + 1
  
  
  # scale_range <- c("effort", 100, 500, 1000, 5000, 10000)
  scale_range <- c(1000)
  
  for(buffer_scale in scale_range){
    print(paste0("Running script for species ", species, " at scale ", buffer_scale))
    scale_start_time <- Sys.time()
    
    ebird_buff <- sf::st_buffer(ebird_ss, dist = ifelse(buffer_scale == "effort",
                                                        ebird_ss$effort_radius_m,
                                                        as.numeric(buffer_scale)))
    ebird_buff_albers <- ebird_buff %>% st_transform("epsg:5070")
    
    gedi_stack <- rast("./predictor_layers/sdm/gedi_stack.tif")
    
    nlcd <- rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif")
    
    treemap <- rast("./predictor_layers/sdm/treemap_stack.tif")
    
    #---------------------------------------
    # GEDI stuff
    #--------------------------------------
    
        ## GEDI Level 2A product-----------------------------------------------------
      #just grab one of the relative height bands
    if(buffer_scale == 5000){
      gedi_stack <- rast("./predictor_layers/sdm/gedi_stack_5000.tif")
    } else if(buffer_scale == 10000){
      gedi_stack <- rast("./predictor_layers/sdm/gedi_stack_10000.tif")
    }
      
      ebird_buff_albers$height <- terra::extract(gedi_stack$height, 
                                                 vect(ebird_buff_albers),
                                                 ID = FALSE,
                                                 fun = sample_and_average
                )$height
      print("     Height done")
   
      ebird_buff_albers$area_short<- terra::extract(gedi_stack$area_short, 
                                                   vect(ebird_buff_albers),
                                                   ID = FALSE,
                                                   fun = sample_and_average
        )$area_short
      print("     Area short done")
      
      ##GEDI Level 2B products--------------------------------------------------------
      ebird_buff_albers$fhd <- terra::extract(gedi_stack$fhd, 
                                                     vect(ebird_buff_albers),
                                                     ID = FALSE,
                                                     fun = sample_and_average
        )$fhd
      print("     FHD done")
      
      ebird_buff_albers$understory <- terra::extract(gedi_stack$understory, 
                                                           vect(ebird_buff_albers),
                                                           ID = FALSE,
                                                           fun = sample_and_average
        )$understory
      
        print("     Understory done")
        
      ebird_buff_albers$biomass <- terra::extract(gedi_stack$biomass, 
                                                  vect(ebird_buff_albers),
                                                  ID = FALSE,
                                                  fun = sample_and_average
        )$biomass
      
        print("     Biomass done")
        
    #     gedi_data_extract[gedi_data_extract$group == group, ] <- ebird_buff_subset
    # }
    
    #-----------------------------------------------------------
    # NLCD cover
    #-----------------------------------------------------------
    if(buffer_scale == 5000){
      nlcd <- rast("./predictor_layers/sdm/nlcd_rescale_5000.tif")
    } else if(buffer_scale == 10000){
      nlcd <- rast("./predictor_layers/sdm/nlcd_rescale_10000.tif")
    }

    nlcd_points <- data.frame(checklist_id = character(length = 0),
                              total_cells = numeric(length = 0),
                              prop_forest = numeric(length = 0),
                              prop_decid = numeric(length = 0),
                              prop_conifer = numeric(length = 0),
                              prop_grass = numeric(length = 0),
                              prop_shrub = numeric(length = 0),
                              prop_water = numeric(length = 0),
                              prop_dev_light = numeric(length = 0),
                              prop_dev_heavy = numeric(length = 0),
                              prop_open = numeric(length = 0))
      
      for(group_select in unique(ebird_buff_albers$group)){
        ebird_subset <- ebird_buff_albers %>% 
          filter(group == group_select) %>%
          mutate(ID = seq(1:nrow(.)))
          
        nlcd_temp <- terra::extract(nlcd, vect(ebird_subset), fun = extract_sample) %>%
          as_tibble() %>%
          rename(NLCD.Land.Cover.Class = any_of("value")) %>% #rename values for the preprocessed raster stacks
          left_join(dplyr::select(ebird_subset, ID, checklist_id), by = "ID") %>%
          group_by(checklist_id) %>%
          summarise(total_cells = length(unlist(NLCD.Land.Cover.Class)),
                    prop_forest = sum(unlist(NLCD.Land.Cover.Class) %in% c(41,42,43,90))/total_cells,
                    prop_decid = sum(unlist(NLCD.Land.Cover.Class) %in% c(41))/total_cells,
                    prop_conifer = sum(unlist(NLCD.Land.Cover.Class) %in% c(42))/total_cells,
                    prop_grass = sum(unlist(NLCD.Land.Cover.Class) %in% c(71, 81))/total_cells,
                    prop_shrub = sum(unlist(NLCD.Land.Cover.Class) %in% c(52))/total_cells,
                    prop_water = sum(unlist(NLCD.Land.Cover.Class) %in% c(11))/total_cells,
                    prop_dev_light = sum(unlist(NLCD.Land.Cover.Class) %in% c(21, 22, 82))/total_cells,
                    prop_dev_heavy = sum(unlist(NLCD.Land.Cover.Class) %in% c(23, 24))/total_cells,
                    prop_open = sum(unlist(NLCD.Land.Cover.Class) %in% c(12,21,31,51,52,71:74,81,82,95))/total_cells)
        
        nlcd_points <- rbind(nlcd_points, nlcd_temp)
      }
    print("     NLCD done")
    
    #--------------------------------------------------------------------------
    #treemap

    if(buffer_scale == 5000){
      treemap <- rast("./predictor_layers/sdm/treemap_rescale_5000.tif")
    } else if(buffer_scale == 10000){
      treemap <- rast("./predictor_layers/sdm/treemap_rescale_10000.tif")
    }
    
    ebird_buff_treemap <- ebird_buff %>% sf::st_transform(crs(treemap))
    
    ebird_buff_treemap <- terra::extract(treemap, ebird_buff_treemap, fun = sample_and_average, ID = FALSE, bind = TRUE)
    
    print("     Treemap done")
    
    
    
    
    #----------------------------
    # combine data together
    ebird_buff_treemap <- st_as_sf(ebird_buff_treemap) %>% 
      st_drop_geometry(ebird_buff_treemap)
    
    ebird_buff_albers2 <- st_drop_geometry(ebird_buff_albers)
    
    
    combined <- dplyr::select(nlcd_points, prop_forest, prop_decid, prop_conifer, prop_shrub, prop_grass, prop_water, 
                              prop_dev_light, prop_dev_heavy, prop_open, checklist_id) %>%
      left_join(dplyr::select(ebird_buff_albers2, height, area_short, fhd, understory, biomass, checklist_id), by = "checklist_id") %>%
      left_join(dplyr::select(ebird_buff_treemap, prop_oak, prop_spruce, prop_hardwood, checklist_id), by = "checklist_id") %>%
      left_join(dplyr::select(ebird_all_obs, species, checklist_id, time_observations_started, duration_minutes, effort_radius_m, species_observed), by = "checklist_id") %>%
      ungroup() %>%
      st_as_sf() %>%
      dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                    lat = sf::st_coordinates(.)[,2]) %>%
      st_drop_geometry()


    combined %>%
      write_csv(., paste0("./environment_vars/", species, "_combined_data_bcr28_local_", buffer_scale,".csv"))

    scale_elapsed <- Sys.time() - scale_start_time
    
    message(paste0("Scale ", buffer_scale, " done! It took ",  round(scale_elapsed, 2), " ", attr(scale_elapsed, "units")))
    
  }
  
  species_elapsed <- Sys.time() - species_start_time
  message(paste0("Species ", species, " done! It took ",  round(species_elapsed, 2), " ", attr(species_elapsed, "units")))
  
}

