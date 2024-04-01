###################
# Get predictor variables from GEE and local sources
# Using GEE for climate variables and extracting others from local data
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

species <- "cerw"

bufferBy100 = function(feature) {
  feature$buffer(feature)
}

bufferFeatureByEffort = function(f) {
  f = ee$Feature(f)
  buffer_size = f$get('effort_radius_m')
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

ebird_ss <- read.csv(paste0("./ebird/", species, "_subsampled_balanced_bcr28.csv"))
ebird_ss$month <- substr(ebird_ss$observation_date, 6, 7)

ebird_ss <- ebird_ss %>% 
  dplyr::filter((protocol_type == "Stationary" | effort_distance_km  < 5) &
                  month %in% c("06","07","08")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  sf::st_set_crs("epsg:4326") %>%
  mutate(effort_radius_m = ifelse(is.na(effort_distance_km), 0.5, 
                                  ifelse(effort_distance_km < 0.5 , 0.5, effort_distance_km)) * 1000 /2) %>%
  dplyr::select(checklist_id, geometry, effort_radius_m, species_observed, 
                time_observations_started, duration_minutes, observation_date) 

# plot(ebird_ss["species_observed"])

ebird_buff <- sf::st_buffer(ebird_ss, dist = ebird_ss$effort_radius_m)
ebird_buff_albers <- ebird_buff %>% st_transform("epsg:5070")

points_ee <-  ebird_ss %>% 
  dplyr::select(checklist_id, effort_radius_m, observation_date) %>%
  sf_as_ee()

# ee_print(points_ee)

# make a GEE polygon buffer about the point
aoi = points_ee$map(bufferFeatureByEffort)
aoi_clim = points_ee$map(function(poly) poly$buffer(2000))
aoi2 = ee$Feature(aoi$first())
ee_print(aoi2)
Map$addLayer(aoi2)


predictor_stack <- terra::rast("./predictor_layers/predictor_stack.grd")

## GEDI Level 2A product-----------------------------------------------------
#just grab one of the relative height bands
ebird_buff_albers$height <- terra::extract(predictor_stack$height, 
                                            vect(ebird_buff_albers),
                                            ID = FALSE,
                                            fun = function(x) mean(x, na.rm = TRUE)
                                            )$height

boxplot(ebird_buff_albers$height ~ ebird_buff_albers$species_observed)

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

#get climate data -------------------------------------------------
#define study period

#TODO reduce number of variables before exporting
#TODO get year-of-observation weather (or week of observation?)

terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

monthstart <- 6
monthend <- 8
years <- ee$List$sequence(2013, 2023)
lag_years <- ee$Number(10)

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

visparams <- list(
  bands = c("def"),
  min = 0,
  max = 1000,
  gamma = 1.4
)

# Map$addLayer(tc_lagged$first(), visparams)+
#   Map$addLayer(aoi)

mean_clim <- ee_extract(tc_lagged$select(c("pet", "aet", "tmmn", "tmmx",
                                          "soil", "pdsi", "vpd", "pr")), 
                       aoi_clim, 
                       fun = ee$Reducer$mean(), 
                       scale = 4000,
                       via = "drive")
# tileScale = 4)

climate_vars <- mean_clim %>%
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

#---------------------------------------------------
# NLCD 

#landfire: var dataset = ee.ImageCollection('LANDFIRE/Vegetation/EVT/v1_4_0');
# 
# #nlcd: 41 is deciduous, 43 is mixed, 42 is evergreen, 
# nlcd = ee$ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')$
#   filterBounds(aoi)$
#   filter(ee$Filter$eq('system:index', '2019'))$first()$
#   select('landcover')
# 
# Map$addLayer(nlcd)
# 
# # lulc <- ee$Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")$
# #   select("discrete_classification")
# 

# TODO try something like this: https://code.earthengine.google.com/83dde965ee8862e4ff418f422cfcf82d

# # This is fairly slow! It doesn't leverage GEE server-side calculations very well.
# # From https://github.com/r-spatial/rgee/issues/199
# # TODO redo with better utilization of GEE
# 
# ee_area_nlcd <- function(img, region, scale = 300) {
#   lista_histo <- list()
#   
#   #process one at a time on GEE
#   
#   for (i in 1:nrow(region)) {
#     region_ee <- region[i, ] %>% sf_as_ee()
#     ee_histo <- img$reduceRegion(
#       reducer = ee$Reducer$frequencyHistogram(),
#       geometry = region_ee,
#       scale = scale
#     )
#     lista_histo[[i]] <- ee_histo$getInfo() %>%
#       map_df(., .f = as.data.frame) %>%
#       mutate(checklist_id = region[[i, 1]] %>% as.vector())
#   }
#   
#   #all local
#   histo_df <- map_df(lista_histo, .f = as.data.frame) %>%
#     mutate_if(is.numeric, .funs = function(x) {
#       x * scale * scale / 10000 #convert from m on side to ha
#     }) %>%
#     replace(is.na(.), 0) %>%
#     pivot_longer(., cols = contains("X"), names_prefix = "X", names_to = "Class", values_to = "Area")
#   return(histo_df)
# }
# 
# nlcd_points <- ee_area_nlcd(img = nlcd, region = ebird_buff, scale = 300)
# 
# nlcd_summary <- nlcd_points %>%
#   group_by(checklist_id) %>%
#   summarize(forest = sum(Area[Class %in% c(41,42,43)])/sum(Area),
#             deciduous = sum(Area[Class %in% c(41)])/sum(Area),
#             mixed = sum(Area[Class %in% c(43)])/sum(Area),
#             coniferous = sum(Area[Class %in% c(42)])/sum(Area),
#             scrub = sum(Area[Class %in% c(52)])/sum(Area),
#             grass = sum(Area[Class %in% c(71, 81)])/sum(Area))


# topography -------------------------------------------------------------------

#A digital elevation model.
dem = ee$Image('NASA/NASADEM_HGT/001')$select('elevation')

#Calculate slope. Units are degrees, range is [0,90).
slope = ee$Terrain$slope(dem)

#Calculate aspect. Units are degrees where 0=N, 90=E, 180=S, 270=W.
# aspect = ee$Terrain$aspect(dem)

slope_extract <- ee_extract(slope, 
                            aoi, 
                            fun = ee$Reducer$mean(), 
                            scale = 100,
                            via = "drive",
                            tileScale = 16)%>%
  dplyr::select(c("checklist_id", "slope")) 

#topographic position index
# ALOS-derived mTPI ranging from negative (valleys) to positive (ridges) values
tpi = ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI")
tpi_extract <- ee_extract(tpi, 
                          aoi, 
                          fun = ee$Reducer$mean(), 
                          scale = 100,
                          via = "drive",
                          tileScale = 16) %>%
  dplyr::select(c("checklist_id", "AVE")) %>%
  mutate(tpi = AVE)

#heat load index
chili = ee$Image("CSP/ERGo/1_0/Global/ALOS_CHILI")
chili_extract = ee_extract(chili, 
                           aoi, 
                           fun = ee$Reducer$mean(), 
                           scale = 100,
                           via = "drive",
                           tileScale = 16) %>%
  dplyr::select(c("checklist_id", "constant")) %>%
  mutate(chili = constant) 


#-------------------------------------------------------------------------------
#bring in local data
#NLCD
#-------------------------------------------------------------------------------

#this is actually much much faster with local NLCD data
#TODO update to match obs with closest NLCD date

nlcd <- rast("D:/Data/NLCD_landcover_2019_release_all_files_20210604/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img")
ebird_buff_nlcd <- sf::st_transform(ebird_buff, crs(nlcd))

nlcd_points <- terra::extract(nlcd, vect(ebird_buff_nlcd), raw = TRUE) %>%
  as.data.frame() %>%
  group_by(ID) %>%
  summarise(total_cells = n(),
            prop_forest = sum(`NLCD Land Cover Class` %in% c(41,42,43))/total_cells,
            prop_decid = sum(`NLCD Land Cover Class` %in% c(41))/total_cells,
            prop_conifer = sum(`NLCD Land Cover Class` %in% c(42))/total_cells,
            prop_grass = sum(`NLCD Land Cover Class` %in% c(52, 71, 81))/total_cells) %>%
  cbind(ebird_buff_nlcd)

#--------------------------------------------------------------------------
#treemap
prop_oak <- terra::rast("./predictor_layers/prop_oak_fia.tiff")
prop_spruce <- terra::rast("./predictor_layers/prop_spruce_fia.tiff")

ebird_buff_treemap <- ebird_buff %>% sf::st_transform(crs(prop_oak))
oak_extract <- terra::extract(prop_oak, vect(ebird_buff_treemap), raw = TRUE) %>%
  as.data.frame() %>%
  group_by(ID) %>%
  summarise(total_cells = n(),
            prop_oak = sum(`tl_id`, na.rm = TRUE)/total_cells) %>%
  ungroup() %>%
  cbind(ebird_buff_treemap)
spruce_fir_extract <- terra::extract(prop_spruce, vect(ebird_buff_treemap), raw = TRUE) %>%
  as.data.frame() %>%
  group_by(ID) %>%
  summarise(total_cells = n(),
            prop_spruce = sum(`tl_id`, na.rm = TRUE)/total_cells) %>%
  ungroup() %>%
  cbind(ebird_buff_treemap)


#combine data ------------------------------------------------------------------
combined <- left_join(dplyr::select(ebird_buff_albers %>% st_drop_geometry(), 
                                    checklist_id, 
                                    height, 
                                    fhd_normal, 
                                    understory_ratio,  
                                    biomass),
                      climate_vars, by = "checklist_id") %>%
  left_join(dplyr::select(nlcd_points, prop_forest, prop_decid, prop_conifer, prop_grass, checklist_id), by = "checklist_id") %>%
  left_join(dplyr::select(oak_extract, prop_oak, checklist_id), by = "checklist_id") %>%
  left_join(dplyr::select(spruce_fir_extract, prop_spruce, checklist_id), by = "checklist_id") %>%
  left_join(slope_extract, by = "checklist_id") %>%
  left_join(dplyr::select(tpi_extract, checklist_id, tpi), by = "checklist_id") %>%
  left_join(dplyr::select(chili_extract, checklist_id, chili), by = "checklist_id") %>%
  left_join(dplyr::select(ebird_ss, checklist_id, time_observations_started, duration_minutes, effort_radius_m, species_observed), by = "checklist_id") %>%
  ungroup() %>%
  st_as_sf() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

combined %>%
  write_csv(., paste0("./environment_vars/", species, "_combined_data_bcr28", Sys.Date(),".csv"))
