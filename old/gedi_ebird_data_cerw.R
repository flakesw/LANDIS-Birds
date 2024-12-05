
###################
# Get predictors from GEE ----------------
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


bufferBy100 = function(feature) {
  feature$buffer(feature)
}

bufferFeatureByEffort = function(f) {
  f = ee$Feature(f)
  buffer_size = f$get('effort_radius_m')
  return(f$buffer(buffer_size))
}

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")
plot(st_geometry(bcr))

bcr_ee <- sf_as_ee(bcr)
# bcr_ee$getInfo()
# Map$addLayer(bcr_ee)

ebird_ss <- read.csv("cerw_subsampled.csv")
ebird_ss$month <- substr(ebird_ss$observation_date, 6, 7)

ebird_ss <- ebird_ss %>% 
  dplyr::filter((protocol_type == "Stationary" | effort_distance_km  < 5) &
           month %in% c("06","07","08")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  sf::st_set_crs("epsg:4326") %>%
  mutate(effort_radius_m = ifelse(is.na(effort_distance_km), 0.5, 
                                  ifelse(effort_distance_km == 0, 0.5, effort_distance_km)) * 1000 /2) %>%
  dplyr::select(checklist_id, geometry, effort_radius_m, species_observed, time_observations_started, duration_minutes) 

plot(ebird_ss["species_observed"])

ebird_buff <- sf::st_buffer(ebird_ss, dist = ebird_ss$effort_radius_m)

points_ee <-  (sf_as_ee(ebird_ss))

# ee_print(points_ee)


# make a GEE polygon buffer about the point
aoi = points_ee$map(bufferFeatureByEffort)


## GEDI Level 2A product-----------------------------------------------------
#just grab one of the relative height bands

heightBands = paste0("rh", seq(0, 100, by = 5))

gedi = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(aoi)$
  map(function(image) {
    return(image$
             updateMask(image$select('quality_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select(heightBands))
  })$median() #collapses collection to image


#version with ee_extract which gets us the mean within each polygon

means <- ee_extract(gedi, 
                    aoi, 
                    fun = ee$Reducer$mean(), 
                    scale = 300,
                    via = "drive")
# sds <- ee_extract(gedi, aoi, fun = ee$Reducer$stdDev(), scale = 30)
# print(means)
# print(sds)

mean_heights <- tidyr::pivot_longer(means,
                                    cols = starts_with("rh"),
                                    names_to = "bin",
                                    names_prefix = "rh",
                                    values_to = "height") %>%
  mutate(bin = as.numeric(bin),
         height = as.numeric(height)) %>%
  group_by(species_observed, bin) #%>%
# summarise(height = mean(height, na.rm = TRUE))
# ggplot(data = mean_heights, aes(x = bin, y = height, col = species_observed)) +
#   geom_point(alpha = 0.3) +
#   geom_line(aes(group = checklist_id), alpha = 0.3) +
#   geom_hline(yintercept = 0) + 
#   geom_smooth(linewidth = 2)

# 
# sd_heights <- tidyr::pivot_longer(sds,
#                                   cols = starts_with("rh"),
#                                   names_to = "bin",
#                                   names_prefix = "rh",
#                                   values_to = "sd") %>%
#   mutate(bin = as.numeric(bin),
#          sd = as.numeric(sd)) %>%
#   left_join(mean_heights, by = c("checklist_id", "bin", "species_observed")) %>%
#   mutate(cv = sd/abs(height))
# ggplot(data = sd_heights, aes(x = bin, y = sd, col = species_observed)) +
#   geom_point(alpha = 0.3) +
#   geom_line(aes(group = checklist_id), alpha = 0.3) +
#   geom_hline(yintercept = 0) + 
#   geom_smooth(linewidth = 2)

mean_max_height <- mean_heights %>%
  filter(bin == 95)
# sd_max_height <- sd_heights %>%
  # filter(bin == 95)
# 


##GEDI Level 2B products--------------------------------------------------------

heightBands = c("fhd_normal", 
                paste0("pavd_z", seq(0, 25, by = 5)))

gedi2b = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(aoi)$
  map(function(image) {
    return(image$
             updateMask(image$select('algorithmrun_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select(heightBands))
  })$median() #collapses collection to image

all_pavd = gedi2b$
  select(c(paste0("pavd_z", seq(0, 25, by = 5))))$
  reduce(ee$Reducer$sum())

understory_prop = gedi2b$select("pavd_z0")$divide(all_pavd)

pavd_extract <- ee_extract(all_pavd, 
                       aoi, 
                       fun = ee$Reducer$mean(), 
                       scale = 100,
                       via = "drive",
                       tileScale = 16)

mean_under <- ee_extract(understory_prop, 
                         aoi, 
                         fun = ee$Reducer$mean(), 
                         scale = 100,
                         via = "drive",
                         tileScale = 16)
mean_fhd <- ee_extract(gedi2b$select("fhd_normal"), 
                       aoi, 
                       fun = ee$Reducer$mean(), 
                       scale = 100,
                       via = "drive",
                       tileScale = 16)
mean_under$understory_proportion <- mean_under$pavd_z0
# sds_under <- ee_extract(understory_prop, aoi, fun = ee$Reducer$stdDev(), scale = 30)
print(mean_under)
# print(sds_under)
# 
# ggplot(data = mean_under, aes(x = species_observed, y = pavd_z0)) +
#   geom_boxplot()

# ggplot(data = mean_pavd[!duplicated(mean_pavd$checklist_id), ], 
#        aes(x = fhd_normal, col = species_observed)) +
#   geom_density()


#-------------------------------------------------------------------
# gedi biomass layer
gedi4 <- ee$Image("LARSE/GEDI/GEDI04_B_002")
gedi4 <- ee$Image("LARSE/GEDI/GEDI04_B_002")$
          # updateMask(gedi4$select('QF')$eq(2))$
          select('MU')

mean_biomass <- ee_extract(gedi4, 
                           aoi, 
                           fun = ee$Reducer$mean(), 
                           scale = 1000,
                           via = "drive",
                           tileScale = 16)


#get climate data -------------------------------------------------
#define study period
startdate = ee$Date$fromYMD(2010,1,1);
enddate = ee$Date$fromYMD(2022,12,31);

terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
  filterBounds(aoi)

#function to aggregate each month
#need to wrap function in ee_utils_pyfunc if using Lists
# Create a number ee$List where each element represent a month
months <- ee$List$sequence(1, 12)

# Function to Calculate a monthly composite
monthly_terraclimate <- function(m) {
  terraclimate$filter(ee$Filter$calendarRange(m, m, "month")) %>%
    ee$ImageCollection$median() 
}

get_preceding_summers <- function(point) {
  obs_date <-  ee$Date$parse(point$get('observation_date'))$aside(print)
  startdate <- obs_date$advance(-2, "years")$aside(print)
  
  terraclimate$
    filter(ee$Filter$calendarRange(m, m, "month"))$
    filter(ee$Filter$date(startdate, enddate))$
    %>%
    ee$ImageCollection$median() 
}

l5_monthly <- months$map(ee_utils_pyfunc(monthly_terraclimate))

#get may and august climate normals
l5_mean_may <- ee$Image(l5_monthly$get(4))
l5_mean_aug <- ee$Image(l5_monthly$get(7))

## Vis parameters.
visparams <- list(
  bands = c("def"),
  min = 0,
  max = 300,
  gamma = 1.4
)

Map$centerObject(aoi, zoom = 10)
Map$addLayer(l5_mean_aug, visparams, name = "Aug") +
  Map$addLayer(l5_mean_may, visparams, name = "May")

terraclimate_may <- ee_extract(l5_mean_may, 
                               aoi, 
                               fun = ee$Reducer$mean(), 
                               scale = 1000,
                               via = "drive",
                               tileScale = 16)
# boxplot(terraclimate_may$def ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$aet ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$soil ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$tmmn ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$tmmx ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$vap ~ terraclimate_may$species_observed)
# boxplot(terraclimate_may$pdsi ~ terraclimate_may$species_observed)

terraclimate_aug <- ee_extract(l5_mean_aug, 
                               aoi, 
                               fun = ee$Reducer$mean(), 
                               scale = 1000,
                               via = "drive",
                               tileScale = 16)
# boxplot(terraclimate_aug$def ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$aet ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$soil ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$tmmn ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$tmmx ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$vap ~ terraclimate_aug$species_observed)
# boxplot(terraclimate_aug$pdsi ~ terraclimate_aug$species_observed)

#TODO aggregate to each month over last decade, to align with weather preceding
#ebird observations

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
combined <- left_join(dplyr::select(ungroup(mean_max_height), checklist_id, height),
                      dplyr::select(terraclimate_aug, !c(time_observations_started, duration_minutes, species_observed)), by = "checklist_id") %>%
  left_join(dplyr::select(nlcd_points, prop_forest, prop_decid, prop_conifer, prop_grass, checklist_id), by = "checklist_id") %>%
  left_join(dplyr::select(oak_extract, prop_oak, checklist_id), by = "checklist_id") %>%
  left_join(slope_extract, by = "checklist_id") %>%
  left_join(tpi_extract, by = "checklist_id") %>%
  left_join(chili_extract, by = "checklist_id") %>%
  left_join(dplyr::select(mean_biomass, checklist_id, MU), by = "checklist_id") %>%
  left_join(dplyr::select(mean_under, checklist_id, understory_proportion), by = "checklist_id") %>%
  left_join(dplyr::select(mean_fhd, checklist_id, fhd_normal), by = "checklist_id") %>%
  left_join(dplyr::select(ebird_ss, checklist_id, time_observations_started, duration_minutes, species_observed), by = "checklist_id") %>%
  ungroup() %>%
  st_as_sf() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

combined %>%
  write_csv(., 'cerw_combined_data.csv')

