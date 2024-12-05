
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

ebird_ss <- read.csv("woth_subsampled.csv")
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

ee_print(gedi)

gediVis = list(
  min= 1,
  max= 60,
  palette= 'red,blue')

# gedi_clip <- gedi$clip(aoi)

Map$addLayer(points_ee, list(color = "yellow"), "Point") +
  Map$addLayer(aoi, list(color = "blue"), "Buffer") +
  Map$addLayer(gedi$select("rh0"), gediVis, 'rh0') +
  Map$addLayer(gedi$select("rh95"), gediVis, 'rh95')


# ## version using sampleRegions, which downloads every GEDI point
# sampled_pixels = gedi$sampleRegions(aoi, scale = 30)
# ee_print(sampled_pixels)
# Map$addLayer(sampled_pixels, list(color = "blue"))
# 
# #this takes a little bit of time
# test <- ee_as_sf(sampled_pixels,
#                  via = "drive")
# # calculate and plot the relative height profiles for each bird sample point
# test1 <- tidyr::pivot_longer(test,
#                              cols = starts_with("rh"),
#                              names_to = "bin",
#                              names_prefix = "rh",
#                              values_to = "height") %>%
#   mutate(bin = as.numeric(bin),
#          height = as.numeric(height)) %>%
#   group_by(checklist_id, bin) %>%
#   summarise(height = mean(height))
# test2 <- tidyr::pivot_longer(test,
#                              cols = starts_with("rh"),
#                              names_to = "bin",
#                              names_prefix = "rh",
#                              values_to = "height") %>%
#   mutate(bin = as.numeric(bin))
#   
# ggplot(test2, aes(y = bin, x = height, group = id)) +
#   geom_line()
# ggplot(test1, aes(x = bin, y = height, col = checklist_id)) +
#   geom_line() +
#   geom_hline(yintercept = 0)


#version with ee_extract which gets us the mean within each polygon

means <- ee_extract(gedi, aoi, fun = ee$Reducer$mean(), scale = 30)
sds <- ee_extract(gedi, aoi, fun = ee$Reducer$stdDev(), scale = 30)
print(means)
print(sds)

mean_heights <- tidyr::pivot_longer(means,
                                    cols = starts_with("rh"),
                                    names_to = "bin",
                                    names_prefix = "rh",
                                    values_to = "height") %>%
  mutate(bin = as.numeric(bin),
         height = as.numeric(height)) %>%
  group_by(species_observed, bin) #%>%
# summarise(height = mean(height, na.rm = TRUE))
ggplot(data = mean_heights, aes(x = bin, y = height, col = species_observed)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(group = checklist_id), alpha = 0.3) +
  geom_hline(yintercept = 0) + 
  geom_smooth(linewidth = 2)


sd_heights <- tidyr::pivot_longer(sds,
                                  cols = starts_with("rh"),
                                  names_to = "bin",
                                  names_prefix = "rh",
                                  values_to = "sd") %>%
  mutate(bin = as.numeric(bin),
         sd = as.numeric(sd)) %>%
  left_join(mean_heights, by = c("checklist_id", "bin", "species_observed")) %>%
  mutate(cv = sd/abs(height))
ggplot(data = sd_heights, aes(x = bin, y = sd, col = species_observed)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(group = checklist_id), alpha = 0.3) +
  geom_hline(yintercept = 0) + 
  geom_smooth(linewidth = 2)

mean_max_height <- mean_heights %>%
  filter(bin == 95)
sd_max_height <- sd_heights %>%
  filter(bin == 95)



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


mean_under <- ee_extract(understory_prop, aoi, fun = ee$Reducer$mean(), scale = 30)
mean_under$understory_proportion <- mean_under$pavd_z0
sds_under <- ee_extract(understory_prop, aoi, fun = ee$Reducer$stdDev(), scale = 30)
print(mean_under)
print(sds_under)

ggplot(data = mean_under, aes(x = species_observed, y = pavd_z0)) +
  geom_boxplot()

ggplot(data = mean_pavd[!duplicated(mean_pavd$checklist_id), ], 
       aes(x = fhd_normal, col = species_observed)) +
  geom_density()


#-------------------------------------------------------------------
# gedi biomass layer
gedi4 <- ee$Image("LARSE/GEDI/GEDI04_B_002")
gedi4 <- ee$Image("LARSE/GEDI/GEDI04_B_002")$
          # updateMask(gedi4$select('QF')$eq(2))$
          select('MU')

mean_biomass <- ee_extract(gedi4, aoi, fun = ee$Reducer$mean(), scale = 30)


#get climate data -------------------------------------------------
#define study period
startdate = ee$Date$fromYMD(2010,1,1);
enddate = ee$Date$fromYMD(2022,12,31);

terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
  filter(ee$Filter$date(startdate, enddate))$
  filterBounds(aoi)

#function to aggregate each month
#need to wrap function in ee_utils_pyfunc if using Lists
# Create a number ee$List where each element represent a month
months <- ee$List$sequence(1, 12)

# Function to Calculate a monthly composite
monthly_l5 <- function(m) {
  terraclimate$filter(ee$Filter$calendarRange(m, m, "month")) %>%
    ee$ImageCollection$median() 
}
l5_monthly <- months$map(ee_utils_pyfunc(monthly_l5))

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

terraclimate_may <- ee_extract(l5_mean_may, aoi, scale = 30)

boxplot(terraclimate_may$def ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$aet ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$soil ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$tmmn ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$tmmx ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$vap ~ terraclimate_may$species_observed)
boxplot(terraclimate_may$pdsi ~ terraclimate_may$species_observed)

terraclimate_aug <- ee_extract(l5_mean_aug, aoi, scale = 30)

boxplot(terraclimate_aug$def ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$aet ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$soil ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$tmmn ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$tmmx ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$vap ~ terraclimate_aug$species_observed)
boxplot(terraclimate_aug$pdsi ~ terraclimate_aug$species_observed)

#TODO aggregate to each month over last decade, to align with weather preceding
#ebird observations

#---------------------------------------------------
# NLCD 

#landfire: var dataset = ee.ImageCollection('LANDFIRE/Vegetation/EVT/v1_4_0');

#nlcd: 41 is deciduous, 43 is mixed, 42 is evergreen, 
nlcd = ee$ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')$
  filterBounds(aoi)$
  filter(ee$Filter$eq('system:index', '2019'))$first()$
  select('landcover')

Map$addLayer(nlcd)

lulc <- ee$Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019")$
  select("discrete_classification")

# This is fairly slow! It doesn't leverage GEE server-side calculations very well.
# From https://github.com/r-spatial/rgee/issues/199
# TODO redo with better utilization of GEE

ee_area_nlcd <- function(img, region, scale = 1000) {
  lista_histo <- list()
  for (i in 1:nrow(region)) {
    region_ee <- region[i, ] %>% sf_as_ee()
    ee_histo <- img$reduceRegion(
      reducer = ee$Reducer$frequencyHistogram(),
      geometry = region_ee,
      scale = scale
    )
    lista_histo[[i]] <- ee_histo$getInfo() %>%
      map_df(., .f = as.data.frame) %>%
      mutate(checklist_id = region[[i, 1]] %>% as.vector())
  }
  histo_df <- map_df(lista_histo, .f = as.data.frame) %>%
    mutate_if(is.numeric, .funs = function(x) {
      x * scale * scale / 10000 #convert from m on side to ha
    }) %>%
    replace(is.na(.), 0) %>%
    pivot_longer(., cols = contains("X"), names_prefix = "X", names_to = "Class", values_to = "Area")
  return(histo_df)
}

nlcd_points <- ee_area_nlcd(img = nlcd, region = ebird_buff, scale = 30)

nlcd_summary <- nlcd_points %>%
  group_by(checklist_id) %>%
  summarize(forest = sum(Area[Class %in% c(41,42,43)])/sum(Area),
            deciduous = sum(Area[Class %in% c(41)])/sum(Area),
            mixed = sum(Area[Class %in% c(43)])/sum(Area),
            coniferous = sum(Area[Class %in% c(42)])/sum(Area),
            scrub = sum(Area[Class %in% c(52)])/sum(Area),
            grass = sum(Area[Class %in% c(71, 81)])/sum(Area))


# topography -------------------------------------------------------------------

#A digital elevation model.
dem = ee$Image('NASA/NASADEM_HGT/001')$select('elevation')

#Calculate slope. Units are degrees, range is [0,90).
slope = ee$Terrain$slope(dem)

#Calculate aspect. Units are degrees where 0=N, 90=E, 180=S, 270=W.
# aspect = ee$Terrain$aspect(dem)

slope_extract <- ee_extract(slope, aoi, scale = 30)%>%
  dplyr::select(c("checklist_id", "slope")) 

#topographic position index
# ALOS-derived mTPI ranging from negative (valleys) to positive (ridges) values
tpi = ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI")
tpi_extract <- ee_extract(tpi, aoi, scale = 30) %>%
  dplyr::select(c("checklist_id", "AVE")) %>%
  mutate(tpi = AVE)

#heat load index
chili = ee$Image("CSP/ERGo/1_0/Global/ALOS_CHILI")
chili_extract = ee_extract(chili, aoi, scale = 30) %>%
  dplyr::select(c("checklist_id", "constant")) %>%
  mutate(chili = constant)

#combine data ------------------------------------------------------------------
combined <- left_join(mean_max_height, 
                      dplyr::select(terraclimate_may, !c(effort_radius_m, species_observed)), by = "checklist_id") %>%
  left_join(nlcd_summary, by = "checklist_id") %>%
  left_join(slope_extract, by = "checklist_id") %>%
  left_join(tpi_extract, by = "checklist_id") %>%
  left_join(chili_extract, by = "checklist_id") %>%
  left_join(dplyr::select(mean_biomass, checklist_id, MU), by = "checklist_id") %>%
  left_join(dplyr::select(mean_under, checklist_id, understory_proportion), by = "checklist_id") %>%
  left_join(dplyr::select(ebird_ss, checklist_id, time_observations_started, duration_minutes), by = "checklist_id") %>%
  ungroup() %>%
  st_as_sf() %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()

combined %>%
  write_csv(., 'woth_combined_data.csv')



###### BRT

combined <- read.csv("woth_combined_data.csv")
combined <- combined[complete.cases(combined), ] %>%
  rename(biomass = MU) %>%
  mutate(species_observed = as.logical(species_observed))
  
# combined$species_observed <- ifelse(as.character(combined$species_observed) == "FOUND", TRUE, FALSE)

library("gbm")
library("tidyterra")

brt <- dismo::gbm.step(data = combined, 
         gbm.x = c("soil", "aet", "def", "pdsi", 
                   "tmmx", "vpd", "tpi", "chili",
                   "biomass", "forest", "deciduous",
                   "understory_proportion", "time_observations_started",
                   "duration_minutes"),
         gbm.y = "species_observed")
summary.gbm(brt)
gbm.interactions(brt)
plot(brt)
gbm.plot(brt)
gbm.perspec(brt, x = 13, y = 12)

brt_simp <- dismo::gbm.simplify(brt)
brt <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                      gbm.x = brt_simp$pred.list[[3]],
                      gbm.y = "species_observed")
gbm.plot(brt,
         n.plots = 12,
         common.scale = TRUE,
         show.contrib = TRUE,
         plot.layout = c(3,4))
gbm.perspec(brt, x = 7, y = 6)


brt_veg <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                           gbm.x = c("height", "deciduous", "fhd_normal", "biomass", "understory_ratio"),
                           gbm.y = "species_observed")
summary.gbm(brt_veg)
gbm.plot(brt_veg)

brt_clim <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("soil", "aet", "tmmx", "tpi", "chili"),
                            gbm.y = "species_observed")
summary.gbm(brt_clim)
gbm.plot(brt_clim)




### predictors
preds <- raster::raster("./predictor_layers/predictor_stack.tiff")

preds <- raster::predict(object = predictor_stack, model = brt,
                         ext = st_bbox(bcr_conical),
                         const = data.frame(time_observations_started = 700,
                                            duration_minutes = 60
                                            ))
values(preds) <- boot::inv.logit(values(preds))
preds <- mask(rast(preds), vect(bcr_conical))
plot(preds)

states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)


ggplot() +
  geom_sf(data = states, colour = "black", fill = NA) +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c(alpha = 0.7)
