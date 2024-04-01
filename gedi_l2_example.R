library("tidyverse")
library("rgee")
library("sf")
ee_Initialize(user = "swflake@ncsu.edu", drive = TRUE)

bufferBy100 = function(feature) {
  feature$buffer(100)
}

#coordinates for some arbitrary points
points <- data.frame(oid = c(0,1,2),
                     longitude = c(-75.44586, -81.86181, -92.85509),
                     latitude = c(39.94798, 40.55233, 37.72540))

#create sf object --> convert to ee object
points_ee <- sf::st_as_sf(points, coords = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  sf_as_ee() #convert to ee feature collection

ee_print(points_ee) #inspect ee featurecollection

# GEE  -------------------------------------------------------------------------

# make a GEE polygon buffer about each point
aoi = points_ee$map(bufferBy100)

#what bands do we want?

heightBands = paste0("rh", seq(0, 100, by = 5))

#get gedi shots within each polygon
gedi_rh = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(aoi)$ #filter to polygon collection
  map(function (image) {
    return(image$
             updateMask(image$select('quality_flag')$eq(1))$ #select high quality shots
             updateMask(image$select('degrade_flag')$eq(0))$ #select high quality shots
             select(heightBands)) #select specified bands
  })$median() #collapses collection to image; because the shots don't overlap, this shouldn't change values

ee_print(gedi_rh) #view the resulting object

gediVis = list(
  min= 1,
  max= 60,
  palette= 'red,blue')

#clip data to polygons. Not sure if this is necessary, and it might be slow for bigger datasets
gedi_clip <- gedi_rh$clip(aoi)

Map$addLayer(points_ee, list(color = "yellow"), "Point") +
  Map$addLayer(aoi, list(color = "blue"), "Buffer") +
  Map$addLayer(gedi_clip$select("rh50"), gediVis, 'rh50') + #height at which 50% of light is returned
  Map$addLayer(gedi_clip$select("rh95"), gediVis, 'rh95') #height at which 95% of light is returned (near max canopy height)

#add these layers to see all the gedi shots 
 # Map$addLayer(gedi_rh$select("rh50"), gediVis, 'rh50') + #height at which 50% of light is returned
  # Map$addLayer(gedi_rh$select("rh95"), gediVis, 'rh95') #height at which 95% of light is returned (near max canopy height)


# ee_image_info(gedi_clip)


sampled_pixels = gedi_rh$sampleRegions(aoi, scale = 30) #only select gedi pixels within buffered polygons
ee_print(sampled_pixels)
Map$addLayer(sampled_pixels, list(color = "blue"))

#export data to local environment, via Google Drive (EE --> Drive --> download to temp folder)
#This makes a table with a row for every shot; this could be very large if you have a lot of
# area being sampled.
#this takes a little bit of time to download
test <- ee_as_sf(sampled_pixels,
                 via = "drive")
nrow(test) #number of shots included
unique(test$oid) #which points are included? If no GEDI shots are within the polygon, those
   #points OIDs will not be included here

#calculate and plot the relative height profiles for each bird sample point
#This takes the mean of each height value within each rh bin. It might or might not
#be the right choice for the analysis you're doing
test1 <- tidyr::pivot_longer(test, 
                             cols = starts_with("rh"),
                             names_to = "bin",
                             names_prefix = "rh",
                             values_to = "height") %>%
  mutate(bin = as.numeric(bin),
         height = as.numeric(height),
         oid = as.factor(oid)) %>%
  group_by(oid, bin) %>%
  summarise(height = mean(height))

#plot relative height profile for each plot (with shots aggregated)
ggplot(test1, aes(x = bin, y = height, group = oid, col = oid)) + 
  geom_line() +
  geom_hline(yintercept = 0)


#--------------------
# PAI metrics
# Similar process, just using Level 2B products
# L2B has plant area index vertical profiles, canopy cover vertical profiles, and 
# plant volume vertical profiles

#choosing different bands
heightBands = c("cover",
                "pai",
                paste0("pai_z", seq(0, 29)))
#same process, just with a different ImageCollection (L2B product), and different bands
gedi_lai = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(aoi)$
  map(function (image) {
    return(image$
             updateMask(image$select('l2b_quality_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select(heightBands))
  })$median() #collapses collection to image

sampled_pixels = gedi_lai$sampleRegions(aoi, scale = 30)
ee_print(sampled_pixels)

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
  group_by(oid, height) %>%
  summarise(pai_z = mean(pai_z))

#plot LAI profile. LAI at each height from the ground (so, below all the canopy at z = 0, PAI = total plot PAI)
ggplot(test1, aes(y = height, x = pai_z, group = oid, col = oid)) + 
  geom_line() +
  geom_hline(yintercept = 0)
