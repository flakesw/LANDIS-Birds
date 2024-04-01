###################
# Get predictor variables from GEE and local sources-----
#------------------------------------------------------

#This script is old, use all_predictor_data.R instead!

#TODO combine the GEDI images into one collection and just do one extract step

library("tidyverse")
library("rgee")
library("sf")
library("dismo")
library("stars")
library("starsExtra")
library("terra")
#ee_install()
# ee_install_upgrade()
ee_Initialize(user = "swflake@ncsu.edu", drive = TRUE)

# write.csv(cerw, "cerw_test.csv")

species <- "gwwa"

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
  # filter(BCR == 28) %>%
  filter(BCR %in% c(12,13,14,22,23,24,27,28,29,30)) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")
# plot(st_geometry(bcr))

bcr_ee <- sf_as_ee(bcr)
# bcr_ee$getInfo()
# Map$addLayer(bcr_ee)
# cerw <- ebird_ss
ebird_ss <- read.csv(paste0("./ebird/", species, "_subsampled_balanced.csv"))
# ebird_ss <- read.csv(paste0("./ebird/", species, "_subsampled_balanced_orig.csv"))
ebird_ss$month <- substr(ebird_ss$observation_date, 6, 7)

ebird_ss <- ebird_ss %>% 
  dplyr::filter((protocol_type == "Stationary" | effort_distance_km  < 5) &
           month %in% c("06","07","08")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  sf::st_set_crs("epsg:4326") %>%
  mutate(effort_radius_m = ifelse(is.na(effort_distance_km), 1, 
                                  ifelse(effort_distance_km < 1 , 1, effort_distance_km)) * 1000 /2) %>%
  dplyr::select(checklist_id, geometry, effort_radius_m, species_observed, 
                time_observations_started, duration_minutes, observation_date, year) 

#TODO more filters: length of time > 10 min, distance minimum for traveling, cap n_observers,
#remove NAs from species_count

# plot(ebird_ss["species_observed"])

# ebird_ss2 <- ebird_ss

# bad <- ebird_ss[c(1,5,6,9,11), ] #1:5 and 6:12 does not work; 13:100 work fine!
# good <- ebird_ss[-c(1,5,6,9,11), ]
#1, 5,6,9,11 bad, 2,3,4,7,8,10,12 good, 
# ebird_ss <- rbind(cerw, ebird_ss)

ebird_buff <- sf::st_buffer(ebird_ss, dist = ebird_ss$effort_radius_m)

points_ee <-  ebird_ss %>% 
  dplyr::select(checklist_id, effort_radius_m, observation_date, year) %>%
  sf_as_ee()

# ee_print(points_ee)

# make a GEE polygon buffer about the point
aoi = points_ee$map(bufferFeatureByEffort)
aoi2 = ee$Feature(aoi$first())
# ee_print(aoi2)
# Map$addLayer(aoi)

## GEDI Level 2A product-----------------------------------------------------
#just grab one of the relative height bands

gedi = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(aoi)$
  map(function(image) {
    return(image$
             updateMask(image$select('quality_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select("rh95"))
  })$median() #collapses collection to image

# ee_print(gedi)
visparams <- list(
  bands = c("rh95"),
  min = 0,
  max = 50,
  gamma = 1.4
)

# Map$addLayer(gedi, visparams)

#version with ee_extract which gets us the mean within each polygon

means <- ee_extract(gedi, 
                    aoi, 
                    fun = ee$Reducer$mean(), 
                    scale = 100,
                    via = "drive",
                    tileScale = 32)
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
  left_join(ebird_ss, by = "checklist_id") %>%
  # group_by(species_observed, bin) %>%
  group_by(checklist_id, species_observed, bin) %>%
  summarise(height = mean(height, na.rm = TRUE))

ggplot(data = mean_heights, aes(x = bin, y = height, col = species_observed)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(group = checklist_id), alpha = 0.3) +
  geom_hline(yintercept = 0) +
  geom_smooth(linewidth = 2)

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
# gedi Level 4B layer
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

#TODO reduce number of variables before exporting
#TODO get year-of-observation weather (or week of observation?)

terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
  select(c("pet", "aet", "tmmn", "tmmx",
            "soil", "pdsi", "vpd"))

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


mean_clim <- ee_extract(tc_lagged, 
                           aoi, 
                           fun = ee$Reducer$mean(), 
                           scale = 4000,
                           via = "drive")

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
              values_from = value)

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
combined <- left_join(dplyr::select(ungroup(mean_max_height), checklist_id, height),
                      climate_vars, by = "checklist_id") %>%
  left_join(dplyr::select(nlcd_points, prop_forest, prop_decid, prop_conifer, prop_grass, checklist_id), by = "checklist_id") %>%
  left_join(dplyr::select(oak_extract, prop_oak, checklist_id), by = "checklist_id") %>%
  left_join(slope_extract, by = "checklist_id") %>%
  left_join(dplyr::select(tpi_extract, tpi), by = "checklist_id") %>%
  left_join(dplyr::select(chili_extract, chili), by = "checklist_id") %>%
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
  write_csv(., paste0("./environment_vars/", species, "_combined_data", Sys.Date(),".csv"))








#TRASH---------------------------------------------------



tc_list <- tc_lagged$toList(999) #convert to a list of images

# function zonalStats(ic, fc, params) {
#   # Initialize internal params dictionary.
#   _params=list(
#     reducer= ee.Reducer.mean(),
#     scale= null,
#     crs= null,
#     bands= null,
#     bandsRename= null,
#     imgProps= null,
#     imgPropsRename= null,
#     datetimeName= 'datetime',
#     datetimeFormat= 'YYYY-MM-dd HH:mm:ss'
#   );
#   
#   // Replace initialized params with provided params.
#   if (params){
#     for (var param in params){
#       _params[param]=params[param] || _params[param];
#     }
#   }
#   
#   // Set default parameters based on an image representative.
#   var imgRep=ic.first();
#   var nonSystemImgProps=ee.Feature(null)
#   .copyProperties(imgRep).propertyNames();
#   if (!_params.bands) _params.bands=imgRep.bandNames();
#   if (!_params.bandsRename) _params.bandsRename=_params.bands;
#   if (!_params.imgProps) _params.imgProps=nonSystemImgProps;
#   if (!_params.imgPropsRename) _params.imgPropsRename=_params
#   .imgProps;
#   
#   // Map the reduceRegions function over the image collection.
#   var results=ic.map(function(img){
#     // Select bands (optionally rename), set a datetime & timestamp property.
#     img=ee.Image(img.select(_params.bands, _params
#                             .bandsRename))
#     // Add datetime and timestamp features.
#     .set(_params.datetimeName, img.date().format(
#       _params.datetimeFormat))
#     .set('timestamp', img.get('system:time_start'));
#     
#     // Define final image property dictionary to set in output features.
#     var propsFrom=ee.List(_params.imgProps)
#     .cat(ee.List([_params.datetimeName,
#                   'timestamp']));
#     var propsTo= ee.List(_params.imgPropsRename)
#     .cat(ee.List([_params.datetimeName,
#                   'timestamp']));
#     var imgProps=img.toDictionary(propsFrom).rename(
#       propsFrom, propsTo);
#     
#     // Subset points that intersect the given image.
#     var fcSub=fc.filterBounds(img.geometry());
#     
#     // Reduce the image by regions.
#     return img.reduceRegions({
#       collection: fcSub,
#       reducer: _params.reducer,
#       scale: _params.scale,
#       crs: _params.crs
#     })
#     // Add metadata to each feature.
#     .map(function(f){
#       return f.set(imgProps);
#     });
#     
#     // Converts the feature collection of feature collections to a single
#     //feature collection.
#   }).flatten();
#   
#   return results;
# }

tc_extract <- function(poly){
  poly <- ee$Feature(poly)
  extract_year <- poly$get("year") %>% ee$Number()
  index <- years$indexOf(extract_year)
  
  clim_vars <- ee$Image(ee$List(tc_list)$get(index))$
    clip(poly)
  out <- clim_vars$reduceRegion(reducer = ee$Reducer$mean(),
                                geometry = poly$geometry(),
                                bestEffort = TRUE)
  
  poly2 <- poly$set(out)
  return(out)
}

Map$centerObject(poly) + Map$addLayer(poly)

test <- aoi$map(tc_extract)
ee_print(test)

test <- tc_list$map(function(im){
  im$reduceRegions(aoi, reducer = ee$Reducer$mean())
})


visparams <- list(
  bands = c("pet"),
  min = 0,
  max = 1000,
  gamma = 1.4
)
Map$addLayer(ee$Image(tc_lagged$select(2)), visparams)
