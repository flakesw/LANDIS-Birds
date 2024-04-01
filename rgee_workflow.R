library("tidyverse")
library("rgee")
library("sf")
ee_Initialize(user = "swflake@ncsu.edu", drive = TRUE)

bufferBy100 = function(feature) {
 feature$buffer(feature)
}

bufferFeatureByEffort = function(f) {
  f = ee$Feature(f)
  buffer_size = f$get('effort_distance_m')
  return(f$buffer(buffer_size))
}



cerw <- readRDS("./ebird_data/cerw.RDS")
cerw$month <- substr(cerw$observation_date, 6, 7)

cerw_filter <- cerw %>% 
  filter((protocol_type == "Stationary" | effort_distance_km  < 3) &
                                 locality_type == "P" &
                                 month %in% c("06","07","08")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  sf::st_set_crs("epsg:4326") %>%
  mutate(effort_distance_m = ifelse(is.na(effort_distance_km), 0.5, effort_distance_km) * 1000) %>%
  dplyr::select(global_unique_identifier, geometry, effort_distance_m) 

# sp <- sf::st_as_sf(cerw_filter, coords = c("longitude", "latitude"), crs = "EPSG:4326") %>%
#   sf::st_transform(crs = "EPSG:32119") #reproject to NC state plane

points <- cerw_filter[1:10, ]
# points$last_edited_date <- rdate_to_eedate(points$last_edited_date, timestamp = TRUE)

points_ee <-  (sf_as_ee(points))

ee_print(points_ee)

# GEE  -------------------------------------------------------------------------

# make a GEE polygon buffer about the point
#lon then lat for GEE
point_ee <- ee$Geometry$Point(list(cerw$longitude[2], cerw$latitude[2]))
aoi = points_ee$map(BufferFeature)

big_buffer <- points_ee$map(bufferFeatureByEffort)

#just grab one of the relative height bands

heightBands = paste0("rh", seq(0, 100, by = 5))

gedi = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(aoi)$
  map(function (image) {
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

gedi_clip <- gedi$clip(aoi)

Map$addLayer(points_ee, list(color = "yellow"), "Point") +
Map$addLayer(aoi, list(color = "blue"), "Buffer") +
Map$addLayer(gedi_clip$select("rh0"), gediVis, 'rh0') +
Map$addLayer(gedi_clip$select("rh100"), gediVis, 'rh100')


# ee_image_info(gedi_clip)


sampled_pixels = gedi$sampleRegions(aoi, scale = 30)
ee_print(sampled_pixels)
Map$addLayer(sampled_pixels, list(color = "blue"))

#this takes a little bit of time
test <- ee_as_sf(sampled_pixels,
                 via = "drive")
#calculate and plot the relative height profiles for each bird sample point
test1 <- tidyr::pivot_longer(test, 
                             cols = starts_with("rh"),
                             names_to = "bin",
                             names_prefix = "rh",
                             values_to = "height") %>%
  mutate(bin = as.numeric(bin),
         height = as.numeric(height)) %>%
  group_by(global_unique_identifier, bin) %>%
  summarise(height = mean(height))
ggplot(test1, aes(x = bin, y = height, col = global_unique_identifier)) + 
  geom_line() +
  geom_hline(yintercept = 0)


#--------------------
# LAI metrics


# make a GEE polygon buffer about the point
#lon then lat for GEE
aoi = points_ee$map(BufferFeature)

big_buffer <- points_ee$map(bufferFeatureByEffort)

#just grab one of the relative height bands

heightBands = c("cover",
                "pai",
                paste0("pai_z", seq(0, 29)))

gedi = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(aoi)$
  map(function (image) {
    return(image$
             updateMask(image$select('l2b_quality_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select(heightBands))
  })$median() #collapses collection to image


# ee_image_info(gedi_clip)


sampled_pixels = gedi$sampleRegions(aoi, scale = 30)
ee_print(sampled_pixels)
Map$addLayer(sampled_pixels, list(color = "blue"))

test <- ee_as_sf(sampled_pixels,
                 via = "drive")
#calculate and plot the vertical cover profiles
test1 <- tidyr::pivot_longer(test, 
                             cols = starts_with("pai_z"),
                             names_to = "bin",
                             names_prefix = "pai_z",
                             values_to = "pai_z") %>%
  mutate(bin = as.numeric(bin),
         cover_z = as.numeric(pai_z)) %>%
  mutate(height = bin * 5) %>%
  group_by(global_unique_identifier, height) %>%
  summarise(pai_z = mean(pai_z))
test1 <- rbind(test1 %>% mutate(height = height + 5), test1) #add the top of each bin


ggplot(test1, aes(y = height, x = pai_z, col = global_unique_identifier)) + 
  geom_line() +
  geom_hline(yintercept = 0)

#------------------------------------------------
# Other variables
#-----------------------------------------------

## Topographic heat load
var dataset = ee.Image('CSP/ERGo/1_0/US/CHILI');
var usChili = dataset.select('constant');
var usChiliVis = {
  min: 0.0,
  max: 255.0,
};
Map.setCenter(-105.8636, 40.3439, 11);
Map.addLayer(usChili, usChiliVis, 'US CHILI');


## LANDFIRE EVT
var evt = ee.ImageCollection("projects/sat-io/open-datasets/landfire/vegetation/EVT");

## NIghttime light (proxy for human use)
var dmsp = ee.ImageCollection("projects/sat-io/open-datasets/Harmonized_NTL/dmsp");
var viirs = ee.ImageCollection("projects/sat-io/open-datasets/Harmonized_NTL/viirs");

