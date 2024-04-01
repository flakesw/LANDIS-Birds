library(sf)
library(raster)
library(tidyverse)
library(basemaps)
library(ggmap)

#overall notes:

#GEDI is hosted on GEE, https://developers.google.com/earth-engine/datasets/catalog/LARSE_GEDI_GEDI02_B_002_MONTHLY


# install.packages("BiocManager")
# BiocManager::install("rhdf5")
# devtools::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = FALSE)

# resources
# https://cornelllabofornithology.github.io/ebird-best-practices/index.html

square_buffer <- function(x, dist){
  sf::st_buffer(x, dist)%>%
    sf::st_transform(crs = "EPSG:4326") %>%
    sf::st_bbox(x) %>%
    return(.)
}


cerw <- readRDS("./ebird_data/cerw.RDS")
cerw$month <- substr(cerw$observation_date, 6, 7)

cerw_filter <- cerw %>% filter((protocol_type == "Stationary" | effort_distance_km  < 3) &
                                 locality_type == "P" &
                                 month %in% c("06","07","08"))

ext <- draw_ext()
saveRDS(ext, file = "test_box.RDS")
str(ext$geometry)

sp <- sf::st_as_sf(cerw_filter, coords = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  sf::st_transform(crs = "EPSG:32119") #reproject to NC state plane

test_box <- st_bbox(ext$geometry)
test <- test_box
# crop <- crop(sp, test_box)
# 


test <- square_buffer(sp[1, ],
                      dist = 1000)

basemap_plot(ext)
points(test)



nc <- ggmap::get_map(location = sf::st_coordinates(sf::st_transform(sp[1, ], crs = "EPSG:4326")))

ggmap::ggmap(nc) + 
  geom_sf(data = sf::st_transform(sp[1, ], crs = "EPSG:4326"), inherit.aes = FALSE, )


# Specifying the date range
daterange=c("2019-07-01","2023-01-01")


library("rGEDI")
library("rhdf5")
# gLevel2A<- gedifinder(product="GEDI02_A", ul_lat = test[4], ul_lon = test[1], test[2], test[3], version="002",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B", ul_lat = test[4], ul_lon = test[1], test[2], test[3], version="002",daterange=daterange)

outdir <- "./gedi test/"
  
# download files (they're big! 600 mb for each level 2b file, 2.5gb for each level 2a)
# gediDownload(filepath=gLevel2B[1], outdir=outdir)

#what files do we have to work with?
gedi_files <- list.files("./gedi test/", full.names = TRUE)

#get info on the h5 file
gedi_test <- h5ls(gedi_files[1])
gedi_test

#check out an attribute -- Hansen tree cover 
h5readAttributes(gedi_files[1], "/BEAM0000/land_cover_data/landsat_treecover")

#plant area index
h5readAttributes(gedi_files[1], "/BEAM0000/pai")

#plant area vertical distribution
h5readAttributes(gedi_files[1], "/BEAM0000/pai_z")



#read some data
testSubset <- h5read(file = gedi_files[1], 
                     name = "/BEAM0000/pai_z")

str(testSubset)


#using rgedi
gedilevel2b <- readLevel2B(level2Bpath = gedi_files[1])
level2BVPM<-getLevel2BVPM(gedilevel2b)
head(level2BVPM[,c("beam","shot_number","pai","fhd_normal","omega","pgap_theta","cover")])

level2BVPM$shot_number<-paste0(level2BVPM$shot_number)

# Converting GEDI Vegetation Profile Biophysical Variables as data.table to SpatialPointsDataFrame
level2BVPM_spdf<-SpatialPointsDataFrame(cbind(level2BVPM$longitude_lastbin,level2BVPM$latitude_lastbin),data=level2BVPM)

level2BPAIProfile<-getLevel2BPAIProfile(gedilevel2b)
head(level2BPAIProfile[,c("beam","shot_number","pai_z0_5m","pai_z5_10m")])




# google earth engine tutorial -------------------------------------------------

library(rgee)
# ee_install(py_env = "rgee") #only run once

# ee_install_upgrade() #run once to fix the version installed


ee_check()#make sure the environment is right

# authenticate and initialize Earth Engine
ee_Initialize(user = "swflake@ncsu.edu")

ee_get_earthengine_path()

#example 1

srtm <- ee$Image("USGS/SRTMGL1_003")

## example 2
# Load an image.
image <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")

# Display the image.
Map$centerObject(image)
Map$addLayer(image, name = "Landsat 8 original image")

# Define visualization parameters in an object literal.
vizParams <- list(
  bands = c("B5", "B4", "B3"),
  min = 5000, max = 15000, gamma = 1.3
)

Map$addLayer(image, vizParams, "Landsat 8 False color")

# Use Map to add features and feature collections to the map. For example,
counties <- ee$FeatureCollection("TIGER/2016/Counties")

Map$addLayer(
  eeObject = counties,
  visParams = vizParams,
  name = "counties"
)

#example 3
collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1")

point <- ee$Geometry$Point(-122.262, 37.8719)
start <- ee$Date("2014-06-01")
finish <- ee$Date("2014-10-01")

filteredCollection <- ee$ImageCollection("LANDSAT/LC08/C01/T1")$
  filterBounds(point)$
  filterDate(start, finish)$
  sort("CLOUD_COVER", TRUE)

first <- filteredCollection$first()

# Define visualization parameters in an object literal.
vizParams <- list(
  bands = c("B5", "B4", "B3"),
  min = 5000,
  max = 15000,
  gamma = 1.3
)

Map$addLayer(first, vizParams, "Landsat 8 image")

# Load a feature collection.
featureCollection <- ee$FeatureCollection("TIGER/2016/States")

# Filter the collection.
filteredFC <- featureCollection$filter(ee$Filter$eq("NAME", "California"))

# Display the collection.
Map$addLayer(
  eeObject = filteredFC,
  visParams = list(palette = "red"),
  name = "California"
)


##example 4


# This function gets NDVI from Landsat 5 imagery.
getNDVI <- function(image) {
  return(image$normalizedDifference(c("B4", "B3")))
}

# Load two Landsat 5 images, 20 years apart.
image1 <- ee$Image("LANDSAT/LT05/C01/T1_TOA/LT05_044034_19900604")
image2 <- ee$Image("LANDSAT/LT05/C01/T1_TOA/LT05_044034_20100611")

# Compute NDVI from the scenes.
ndvi1 <- getNDVI(image1)
ndvi2 <- getNDVI(image2)

# Compute the difference in NDVI.
ndviDifference <- ndvi2$subtract(ndvi1)

ndviParams <- list(palette = c(
  "#d73027", "#f46d43", "#fdae61", "#fee08b",
  "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850"
))
ndwiParams <- list(min = -0.5, max = 0.5, palette = c("FF0000", "FFFFFF", "0000FF"))

Map$centerObject(ndvi1)
Map$addLayer(ndvi1, ndviParams, "NDVI 1") +
  Map$addLayer(ndvi2, ndviParams, "NDVI 2") +
  Map$addLayer(ndviDifference, ndwiParams, "NDVI difference")


#example 5
# This function gets NDVI from Landsat 8 imagery.
addNDVI <- function(image) {
  return(image$addBands(image$normalizedDifference(c("B5", "B4"))))
}

# Load the Landsat 8 raw data, filter by location and date.
collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1")$
  filterBounds(ee$Geometry$Point(-122.262, 37.8719))$
  filterDate("2014-06-01", "2014-10-01")

# Map the function over the collection.
ndviCollection <- collection$map(addNDVI)

first <- ndviCollection$first()
print(first$getInfo())

bandNames <- first$bandNames()
print(bandNames$getInfo())

#example 6

# Load a Landsat 8 collection.
collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1")$
  filterBounds(ee$Geometry$Point(-122.262, 37.8719))$
  filterDate("2014-01-01", "2014-12-31")$
  sort("CLOUD_COVER")

# Compute the median of each pixel for each band of the 5 least cloudy scenes.
median <- collection$limit(5)$reduce(ee$Reducer$median())

# Define visualization parameters in an object literal.
vizParams <- list(
  bands = c("B5_median", "B4_median", "B3_median"),
  min = 5000, max = 15000, gamma = 1.3
)

Map$addLayer(
  eeObject = median,
  visParams = vizParams,
  name = "Median image"
)


#example 7


# Load and display a Landsat TOA image.
image <- ee$Image("LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140318")
Map$addLayer(
  eeObject = image,
  visParams = list(bands = c("B4", "B3", "B2"), max = 30000),
  name = "Landsat 8"
)

Map$setCenter(-122.262, 37.8719)
# Create an arbitrary rectangle as a region and display it.
region <- ee$Geometry$Rectangle(-122.2806, 37.1209, -122.0554, 37.2413)
Map$addLayer(
  eeObject = region,
  name = "Region"
)

# Get a dictionary of means in the region.  Keys are bandnames.
mean <- image$reduceRegion(
  reducer = ee$Reducer$mean(),
  geometry = region,
  scale = 30
)

print(mean$getInfo())


#example 8--------------------------------------------------------------------

# This function gets NDVI from Landsat 8 imagery.
addNDVI <- function(image) {
  return(image$addBands(image$normalizedDifference(c("B5", "B4"))))
}

# This function masks cloudy pixels.
cloudMask <- function(image) {
  clouds <- ee$Algorithms$Landsat$simpleCloudScore(image)$select("cloud")
  return(image$updateMask(clouds$lt(10)))
}

# Load a Landsat collection, map the NDVI and cloud masking functions over it.
collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1_TOA")$
  filterBounds(ee$Geometry$Point(c(-122.262, 37.8719)))$
  filterDate("2014-03-01", "2014-05-31")$
  map(addNDVI)$
  map(cloudMask)

# Reduce the collection to the mean of each pixel and display.
meanImage <- collection$reduce(ee$Reducer$mean())
vizParams <- list(
  bands = c("B5_mean", "B4_mean", "B3_mean"),
  min = 0,
  max = 0.5
)

Map$addLayer(
  eeObject = meanImage,
  visParams = vizParams,
  name = "mean"
)

# Load a region in which to compute the mean and display it.
counties <- ee$FeatureCollection("TIGER/2016/Counties")
santaClara <- ee$Feature(counties$filter(
  ee$Filter$eq("NAME", "Santa Clara")
)$first())

Map$addLayer(
  eeObject = santaClara,
  visParams = list(palette = "yellow"),
  name = "Santa Clara"
)

# Get the mean of NDVI in the region.
mean <- meanImage$select("nd_mean")$reduceRegion(
  reducer = ee$Reducer$mean(),
  geometry = santaClara$geometry(),
  scale = 30
)

# Print mean NDVI for the region.
cat("Santa Clara spring mean NDVI:", mean$get("nd_mean")$getInfo())


# example ------------------------------------------------------------------
# Load a Landsat 8 ImageCollection for a single path-row.
collection <- ee$ImageCollection("LANDSAT/LC08/C01/T1_TOA")$
  filter(ee$Filter$eq("WRS_PATH", 44))$
  filter(ee$Filter$eq("WRS_ROW", 34))$
  filterDate("2014-03-01", "2014-09-01")
ee_print(collection)

# Get the number of images.
count <- collection$size()
cat("Count: ", count$getInfo())

# Get the date range of images in the collection.
range <- collection$reduceColumns(
  ee$Reducer$minMax(),
  list("system:time_start")
)

col_min <- eedate_to_rdate(range$get("min"))
col_max <- eedate_to_rdate(range$get("max"))
cat("Date range: ", as.character(col_min), as.character(col_max))

# Get statistics for a property of the images in the collection.
sunStats <- collection$aggregate_stats("SUN_ELEVATION")
cat("Sun elevation statistics: ")
sunStats$getInfo()

# Sort by a cloud cover property, get the least cloudy image.
image <- ee$Image(collection$sort("CLOUD_COVER")$first())
cat("Least cloudy image: ")
image$getInfo()

# Limit the collection to the 10 most recent images.
recent <- collection$sort("system:time_start", FALSE)$limit(10)
cat("Recent images: ")
recent$getInfo()


#gedi example--------------------------------------------------------------

#make polygon
polygon <- ee$Geometry$Polygon(
  list(
    c(test$xmin, test$ymin), c(test$xmin, test$ymax), 
    c(test$xmax, test$ymax), c(test$xmax, test$ymin)
  )
)

# Create a planar polygon.
planarPolygon <- ee$Geometry(polygon, {}, FALSE)
polygon <- ee$FeatureCollection(polygon)
planarPolygon <- ee$FeatureCollection(planarPolygon)

# Display the polygons by adding them to the map.
Map$centerObject(polygon, zoom = 13)
Map$addLayer(polygon, list(color = "FF0000"), "geodesic polygon") +
  Map$addLayer(planarPolygon, list(color = "000000"), "planar polygon")



qualityMask = function(im) {
  return(im$updateMask(im$select('l2b_quality_flag')$eq(1))$
           updateMask(im$select('degrade_flag')$eq(0)))
}

#load GEDI data
dataset = ee$FeatureCollection('LARSE/GEDI/GEDI02_B_002')$
  ee_print()
  

ee_print(dataset)


gediVis = list(
  min= 1,
  max= 60,
  palette= 'red, green, blue')
Map$setCenter(12.60033, 51.01051, 12)
Map$addLayer(dataset, gediVis, 'Solar Elevation')












#-------------------------------------------------------------------------------
# example using gee; convert to rgee    ----------------------------------------

# example from stackoverflow ---------------------------------------------------
#https://gis.stackexchange.com/questions/426742/error-in-random-forest-classification-for-the-gedi-datasets-in-gee
#https://code.earthengine.google.com/37d6104d659f4dc63f36baccde372a16

Paracou = 
  ee$Geometry$Polygon(
    list(c(-52.992228448173776, 5.3639624330339855),
      c(-52.992228448173776, 5.216962579558439),
      c(-52.898158013603464, 5.216962579558439),
      c(-52.898158013603464, 5.3639624330339855)))

qualityMask = function(im) {
  return (im$updateMask(im$select('quality_flag')$eq(1))$
  updateMask(im$select('degrade_flag')$eq(0)))
}

img = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(Paracou)$
  map(qualityMask)$
  select('rh13','rh99')$
  median()
  
  
gediVis = list(
  min= 1,
  max= 4)

img = img$clip(Paracou)
Map$centerObject(Paracou,12)
Map$addLayer(img, gediVis, 'rh13-99')
print(img$getInfo());

### sample() doesn't work for some reason!
points = img$sample(
  list(
  'region' = Paracou,
  'scale' = 30,
  'numPixels' = 200000,
  'seed' = 0,
  'geometries'= TRUE))

Map.addLayer(points,{},'training',true);
print(points.size().getInfo());
print(points.first().getInfo());
//SPLIT AND TRAINING
// Use these bands for prediction
var bands = ['rh13','rh99'];

var label = 'Relative Height';

var sample = img.select(bands).sampleRegions({
  collection: points,
  properties : [label],
  scale : 30
});
var sample = sample.randomColumn();
var split = 0.7;
var training = sample.filter(ee.Filter.lt('random',split));
var validation = sample.filter(ee.Filter.gte('random',split));
print(training.first().getInfo());
print(validation.first().getInfo());
var classifier = ee.Classifier.smileRandomForest(10).train(training,label,bands);
// Classification
var result = img.select(bands).classify(classifier);
Map.addLayer(result.randomVisualizer(),{},'classified');


#-------------------------------------------------------------------------------
# GEDI example from stackoverflow -- in answers to question
#https://gis.stackexchange.com/questions/426742/error-in-random-forest-classification-for-the-gedi-datasets-in-gee

aoi = ee$Geometry$Polygon(list(
  c(-52.992228448173776, 5.3639624330339855),
  c(-52.992228448173776, 5.216962579558439),
  c(-52.898158013603464, 5.216962579558439),
  c(-52.898158013603464, 5.3639624330339855)))$
    buffer(1e4)

heightBand = 'rh95'

gedi = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(aoi)$
  map(function (image) {
  return(image$
    updateMask(image$select('quality_flag')$eq(1))$
    updateMask(image$select('degrade_flag')$eq(0))$
    select(heightBand))
  })$median()


gediVis = list(
  min= 1,
  max= 60,
  palette= 'red,blue')

Map$centerObject(aoi,12)
Map$addLayer(gedi, gediVis, 'rh95')
print(img$getInfo());












referenceData = gedi$
  map(toReferenceData)$  
  filterMetadata('count', 'not_equals', 0)

trainingData = referenceData$
  map(toTrainingData)$
  flatten()

composite = toComposite('2021-01-01', '2022-01-01')

classifier = ee$Classifier$smileRandomForest(25)$
  train(trainingData, heightBand, composite$bandNames())$
  setOutputMode('REGRESSION')

classification = composite$classify(classifier)


print('trainingData', trainingData)
print('classifier', classifier$explain())
Map$addLayer(referenceData.flatten(), null, 'reference data', false)
Map$addLayer(gedi$mosaic()$clip(aoi), {bands: heightBand, min: 0, max: 40, palette: '#fffdcd, #e1cd73, #aaac20, #5f920c, #187328, #144b2a, #172313'}, heightBand)
Map$addLayer(composite, {bands: 'red,green,blue', min: 0, max: 2000, gamma: 1.6}, 'optical')
Map$addLayer(composite, {bands: 'VV,VH,ratio_VV_VH', min: [-20, -25, 3], max: [0, -5, 14]}, 'sar')
Map$addLayer(classification, {min: 0, max: 40, palette: '#fffdcd, #e1cd73, #aaac20, #5f920c, #187328, #144b2a, #172313'}, 'height')
Map$centerObject(aoi, 12)


toReferenceData <- function(gedi) {
  referenceData = gedi$
    sample(list(
    region= aoi,
    scale = 25,
    numPixels = 200000,
    geometries= TRUE
  ))$
    map(function (feature) {
    #Can only use ints as classes
    return (feature$set('rh99', feature$getNumber('rh99')$round()))
  })
  
  return (referenceData$
    set(list(
    year= gedi$getNumber('year'),
    month= gedi$getNumber('month'),
    count= referenceData$size()
  )))
}


toTrainingData <- function (referenceData) {
  year = referenceData$getNumber('year')
  month = referenceData$getNumber('month')
  fromDate = ee$Date$fromYMD(year, month, 1)
  toDate = fromDate$advance(1, 'month')
  composite = toComposite(fromDate, toDate)
  return (composite$sampleRegions(list(
    collection= referenceData,
    properties = "heightBand",
    scale = 30)))
}


toComposite <- function (fromDate, toDate) {
  return (toOpticalComposite(fromDate, toDate)$
            addBands(toSarComposite(fromDate, toDate)))
}

toSarComposite <- function (fromDate, toDate) {  
  composite = ee$ImageCollection('COPERNICUS/S1_GRD_FLOAT')$
    filter(ee$Filter$and(
    ee$Filter$bounds(aoi),
    ee$Filter$date(fromDate, toDate),
    ee$Filter$eq('instrumentMode', 'IW'),
    ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'),
    ee$Filter$listContains('transmitterReceiverPolarisation', 'VH')
  ))$
    map(function (image) {
    return (maskBorder(
      toDb(
        toGamma0(image)
      )
    )$select(c('VV', 'VH')))
  })$median()
  
  return (composite$addBands(
    composite$select('VV')$divide(composite.select('VH'))$rename('ratio_VV_VH')
  ))
}
  
  
toGamma0  <- function(image) {
    gamma0 = image$expression('i/(cos(angle * pi / 180))', list(
      'i'= image$select(c('VV', 'VH')),
      'angle'= image$select('angle'),
      'pi'= Math$PI
    ))
    return (image$addBands(gamma0, null, TRUE))    
}
  
toDb <- function(image) {
    return (image$addBands(
      image$select(c('VV', 'VH'))$log10()$multiply(10), null, TRUE
    )    )
}
  
maskBorder <-  function(image) {
    angle = image$select('angle')
    return (image$updateMask(
      angle$gt(31)$and(angle.lt(45))
    )    )
}


toOpticalComposite <- function (fromDate, toDate) {
  return (ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')$
            filterDate(fromDate, toDate)$
            filterBounds(aoi)$
            map(function (image) {
                qa = image.select('QA_PIXEL')
                cloudShadow = bitwiseExtract(qa, 4)
                snow = bitwiseExtract(qa, 5)$rename('snow')
                cloud = bitwiseExtract(qa, 6)$not()$rename('cloud')
                return (image$
                          updateMask(cloudShadow$not()$
                                       and(snow$not())$
                                       and(cloud$not())
                  )$
                    select(list(
                      c('SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'),
                      c('blue', 'green', 'red', 'nir', 'swir1', 'swir2')))$
                     multiply(2.75e-05)$
                  add(-0.2)$
                  multiply(10000)$
                  int16()
                $median()
                )
            }
            )
  )
}

bitwiseExtract = function (value, fromBit, toBit) {
  if (is.null(toBit))
    toBit = fromBit
  maskSize = eeNumber(1)$add(toBit)$subtract(fromBit)
  mask = ee$Number(1)$leftShift(maskSize)$subtract(1)
  return (value$rightShift(fromBit)$bitwiseAnd(mask))
}





#------------------------------

#feature operations


# Create two circular geometries.
poly1 <- ee$Geometry$Point(c(-50, 30))$buffer(1e6)
poly2 <- ee$Geometry$Point(c(-40, 30))$buffer(1e6)

# Display polygon 1 in red and polygon 2 in blue.
Map$setCenter(-45, 30)
Map$addLayer(poly1, list(color = "FF0000"), "poly1") +
  Map$addLayer(poly2, list(color = "0000FF"), "poly2")

# Compute the intersection, display it in blue.
intersection <- poly1$intersection(poly2, ee$ErrorMargin(1))
Map$addLayer(intersection, list(color = "00FF00"), "intersection")

# Compute the union, display it in magenta.
union <- poly1$union(poly2, ee$ErrorMargin(1))
Map$addLayer(union, list(color = "FF00FF"), "union")

# Compute the difference, display in yellow.
diff1 <- poly1$difference(poly2, ee$ErrorMargin(1))
Map$addLayer(diff1, list(color = "FFFF00"), "diff1")

# Compute symmetric difference, display in black.
symDiff <- poly1$symmetricDifference(poly2, ee$ErrorMargin(1))
Map$addLayer(symDiff, list(color = "000000"), "symmetric difference")


#----------------------------------------------------------------------------


fc <- ee$FeatureCollection(list(
  ee$Feature(
    ee$Geometry$Polygon(
      list(
        c(-109.05, 41),
        c(-109.05, 37),
        c(-102.05, 37),
        c(-102.05, 41)
      )
    ),
    list(name = "Colorado", fill = 1)
  ),
  ee$Feature(
    ee$Geometry$Polygon(
      list(
        c(-114.05, 37.0),
        c(-109.05, 37.0),
        c(-109.05, 41.0),
        c(-111.05, 41.0),
        c(-111.05, 42.0),
        c(-114.05, 42.0)
      )
    ),
    list(name = "Utah", fill = 2)
  )
))

# Fill, then outline the polygons into a blank image.
image1 <- ee$Image(0)$mask(0)$toByte()
image2 <- image1$paint(fc, "fill") # Get color from property named 'fill'
image3 <- image2$paint(fc, 3, 5) # Outline using color 3, width 5.

Map$setCenter(-107, 41, 4)
Map$addLayer(
  eeObject = image3,
  visParams = list(
    palette = c("000000", "FF0000", "00FF00", "0000FF"),
    max = 3,
    opacity = 0.5
  ),
  name = "Colorado & Utah"
)
