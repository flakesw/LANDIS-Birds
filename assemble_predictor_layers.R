# #-------------------------------------------------------------------------------
# # Get predictor map surfaces
# #-------------------------------------------------------------------------------
# #-----
# # getting surfaces for prediction

library("tidyverse")
library("rgee")
library("sf")
library("terra")
ee_Initialize(user = "swflake@ncsu.edu", drive = TRUE)

output_scale <- 100

#region to extract prediction layers for
bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:5070") %>%
  filter(BCR == 28) %>%
  # filter(BCR %in% c(12,13,14,22,23,24,27,28,29,30)) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")
# plot(st_geometry(bcr))

bcr_ee <- sf_as_ee(bcr)

proj = ee$Projection('EPSG:5070')

gedi_proj = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  select('quality_flag')$
  first()$
  projection()
gedi_proj$getInfo()

#---------------------------------------------------
# rh95 
height_red = ee$ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')$
  filterBounds(bcr_ee$geometry())$
  map(function(image) {
    return(image$select("rh95")$multiply(1000)$toShort()$
             updateMask(image$select('quality_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))
             )
  })$
  median()$
  unmask(-1) #set the NA value to -1 so we can fix it later

height_download <- ee_image_to_drive(height_red,
                                     description = "rh95_raster",
                                     fileNamePrefix = "rh95_raster",
                                     folder = "prediction_layers",
                                     region = bcr_ee$geometry(),
                                     scale = 25,
                                     crs = proj,
                                     maxPixels = 4097260648,
                                     fileFormat = 'GeoTIFF')
height_download$start()
ee_monitoring(height_download)

# #----------
# # FHD

fhd_red = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(bcr_ee$geometry())$
  map(function(image) {
    return(image$select("fhd_normal")$multiply(1000)$toShort()$
             updateMask(image$select('algorithmrun_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))
             )
  })$
  median()$
  unmask(-1)

fhd_download <- ee_image_to_drive(fhd_red,
                                     description = "fhd_raster",
                                     fileNamePrefix = "fhd_raster",
                                     folder = "prediction_layers",
                                     region = bcr_ee$geometry(),
                                     scale = 25,
                                     crs = proj,
                                     maxPixels = 4097260648,
                                     fileFormat = 'GeoTIFF')

fhd_download$start()
ee_monitoring(fhd_download)


### pavd
# heightBands = c(paste0("pavd_z", seq(0, 25, by = 5)))
# 
# pavd_red = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
#   filterBounds(bcr_ee$geometry())$
#   map(function(image) {
#     return(image$
#              updateMask(image$select('algorithmrun_flag')$eq(1))$
#              updateMask(image$select('degrade_flag')$eq(0))$
#              select(heightBands))
#   })$
#   median()$ #collapses collection to image
#   reproject(gedi_proj)$
#   reduceResolution(reducer = ee$Reducer$mean(),
#                    maxPixels = 65536)
# 
# #make a map!
# pavd_download <- ee_image_to_drive(pavd_red,
#                                    description = "pavd_raster_with_reducer",
#                                    folder = "prediction_layers",
#                                    region = bcr_ee$geometry(),
#                                    scale = output_scale,
#                                    crs = proj,
#                                    maxPixels = 2048630324)
# pavd_download$start()
# ee_monitoring(pavd_download)

# #---------------------------
# biomass raster
# gedi4b <- ee$Image("LARSE/GEDI/GEDI04_B_002")$
#   # updateMask(gedi4$select('QF')$eq(2))$
#   select('MU')
# 
# biomass_download_4b <- ee_image_to_drive(gedi4b,
#                                       description = "biomass_raster_4b",
#                                       folder = "prediction_layers",
#                                       region = bcr_ee$geometry(),
#                                       scale = 1000,
#                                       crs = proj)
# biomass_download_4b$start()
# ee_monitoring(biomass_download)

gedi4a <- ee$ImageCollection("LARSE/GEDI/GEDI04_A_002_MONTHLY")$
  filterBounds(bcr_ee$geometry())$
  map(function(image) {
    return(image$select('agbd')$
            toShort()$
            updateMask(image$select('algorithm_run_flag')$eq(1))$
            updateMask(image$select('degrade_flag')$eq(0)))
  })$
  median()$
  unmask(-1)

biomass_download_4a <- ee_image_to_drive(gedi4a,
                                      description = "biomass_raster_4a",
                                      fileNamePrefix = "biomass_raster_4a",
                                      folder = "prediction_layers",
                                      region = bcr_ee$geometry(),
                                      scale = 25,
                                      crs = proj,
                                      maxPixels = 4097260648,
                                      fileFormat = 'GeoTIFF'
                                      )
biomass_download_4a$start() #TODO - with degrade flag
ee_monitoring(biomass_download_4a)

#--------------
# proportion understory
heightBands = c(paste0("pavd_z", seq(0, 25, by = 5)))

gedi2b = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(bcr_ee$geometry())$
  map(function(image) {
    return(image$
             updateMask(image$select('algorithmrun_flag')$eq(1))$
             updateMask(image$select('degrade_flag')$eq(0))$
             select(heightBands))
  })$
  median()

all_pavd = gedi2b$
  select(c(paste0("pavd_z", seq(0, 25, by = 5))))$
  reduce(ee$Reducer$sum())

understory_prop = gedi2b$select("pavd_z0")$
  divide(all_pavd)$
  multiply(1000)$
  toShort()$
  unmask(-1)

under_download <- ee_image_to_drive(understory_prop, 
                                    description = "understory_raster",
                                    fileNamePrefix = "understory_raster",
                                    folder = "prediction_layers",
                                    region = bcr_ee$geometry(),
                                    scale = 25,
                                    crs = proj,
                                    maxPixels = 4097260648,
                                    fileFormat = 'GeoTIFF')

under_download$start()
ee_monitoring(under_download)



# #---------------------------
# #topographic position index
# ALOS-derived mTPI ranging from negative (valleys) to positive (ridges) values
tpi = ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI")
#make a map!
tpi_download <- ee_image_to_drive(tpi, 
                                  description = "tpi_raster",
                                  fileNamePrefix = "tpi_raster",
                                  folder = "prediction_layers",
                                  region = bcr_ee$geometry(),
                                  scale = output_scale,
                                  crs = proj,
                                  maxPixels = 4097260648,
                                  fileFormat = 'GeoTIFF')
tpi_download$start()
# ee_monitoring(tpi_download)

#heat load index
chili = ee$Image("CSP/ERGo/1_0/Global/ALOS_CHILI")
chili_download <- ee_image_to_drive(chili,
                                  description = "chili_raster",
                                  fileNamePrefix = "chili_raster",
                                  folder = "prediction_layers",
                                  region = bcr_ee$geometry(),
                                  scale = output_scale,
                                  crs = proj,
                                  maxPixels = 4097260648,
                                  fileFormat = 'GeoTIFF')
chili_download$start()
ee_monitoring(chili_download)


#Calculate slope. Units are degrees, range is [0,90).
dem = ee$Image('NASA/NASADEM_HGT/001')$select('elevation')
slope = ee$Terrain$slope(dem)

slope_download <- ee_image_to_drive(slope, 
                            description = "slope_raster",
                            fileNamePrefix = "slope_raster",
                            folder = "prediction_layers",
                            region = bcr_ee$geometry(),
                            scale = output_scale,
                            crs = proj,
                            maxPixels = 4097260648,
                            fileFormat = 'GeoTIFF')
slope_download$start()
ee_monitoring(slope_download)

# #----------------
# # climate layers
#TODO reduce number of variables before exporting
#TODO get year-of-observation weather (or week of observation?)

terraclimate <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

monthstart <- 6
monthend <- 8
years <- ee$Number(2023)
lag_years <- ee$Number(10)

endyear <- ee$Number(2023)
startyear <- endyear$subtract(lag_years)

tc_mean <- terraclimate$filter(ee$Filter$calendarRange(startyear, endyear, "year"))$
  filter(ee$Filter$calendarRange(ee$Number(monthstart), ee$Number(monthend), "month"))$
  filterBounds(bcr_ee)$
  mean()


#soil moisture
soil_download <- ee_image_to_drive(tc_mean$select('soil'),
                                    description = "soil_raster",
                                   fileNamePrefix = "soil_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = output_scale,
                                   crs = proj,
                                   maxPixels = 4097260648,
                                   fileFormat = 'GeoTIFF')
soil_download$start()
# ee_monitoring(soil_download)


#aet
aet_download <- ee_image_to_drive(tc_mean$select('aet'),
                                   description = "aet_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = output_scale,
                                   crs = proj,
                                  maxPixels = 2048630324)
aet_download$start()
# ee_monitoring(aet_download)

#tmax
tmax_download <- ee_image_to_drive(tc_mean$select('tmmx'),
                                  description = "tmax_raster",
                                  folder = "prediction_layers",
                                  region = bcr_ee$geometry(),
                                  scale = output_scale,
                                  crs = proj,
                                  maxPixels = 2048630324)
tmax_download$start()
# ee_monitoring(tmax_download)

#tmin
tmin_download <- ee_image_to_drive(tc_mean$select('tmmn'),
                                   description = "tmin_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = output_scale,
                                   crs = proj,
                                   maxPixels = 2048630324)
tmin_download$start()
ee_monitoring(tmin_download)

#ppt
ppt_download <- ee_image_to_drive(tc_mean$select('pr'),
                                   description = "ppt_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = output_scale,
                                   crs = proj,
                                  maxPixels = 2048630324)
ppt_download$start()
# ee_monitoring(ppt_download)


#vpd
vpd_download <- ee_image_to_drive(tc_mean$select('vap'),
                                  description = "vpd_raster",
                                  folder = "prediction_layers",
                                  region = bcr_ee$geometry(),
                                  scale = output_scale,
                                  crs = proj,
                                  maxPixels = 2048630324)
vpd_download$start()
# ee_monitoring(vpd_download)

#cwd
cwd_download <- ee_image_to_drive(tc_mean$select('def'),
                                  description = "cwd_raster",
                                  folder = "prediction_layers",
                                  region = bcr_ee$geometry(),
                                  scale = output_scale,
                                  crs = proj,
                                  maxPixels = 2048630324)
cwd_download$start()
# ee_monitoring(cwd_download)

#pdsi
pdsi_download <- ee_image_to_drive(tc_mean$select('pdsi'),
                                   description = "pdsi_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = output_scale,
                                   crs = proj,
                                   maxPixels = 2048630324)
pdsi_download$start()
# ee_monitoring(pdsi_download)


#------------------------------
#Download all the rasters
#---------------------------------

#retrieve files from Google Drive. Check if they're ready first!

# ee_drive_to_local(height_download, dsn = "./predictor_layers/rh95.tif", overwrite = TRUE)
# ee_drive_to_local(fhd_download, dsn = "./predictor_layers/fhd.tif", overwrite = TRUE) #TODO
# ee_drive_to_local(biomass_download_4a, dsn = "./predictor_layers/biomass_degrade_flag.tif", overwrite = TRUE)
# ee_drive_to_local(under_download, dsn = "./predictor_layers/understory.tif", overwrite = TRUE)
# 
# ee_drive_to_local(soil_download, dsn = "./predictor_layers/soil_moisture.tif", overwrite = TRUE)
# ee_drive_to_local(aet_download, dsn = "./predictor_layers/aet.tif", overwrite = TRUE)
# ee_drive_to_local(tmax_download, dsn = "./predictor_layers/tmax.tif", overwrite = TRUE)
# ee_drive_to_local(tmin_download, dsn = "./predictor_layers/tmin.tif", overwrite = TRUE)
# ee_drive_to_local(ppt_download, dsn = "./predictor_layers/ppt.tif", overwrite = TRUE)
# ee_drive_to_local(vpd_download, dsn = "./predictor_layers/vpd.tif", overwrite = TRUE)
# ee_drive_to_local(cwd_download, dsn = "./predictor_layers/cwd.tif", overwrite = TRUE)
# ee_drive_to_local(pdsi_download, dsn = "./predictor_layers/pdsi.tif", overwrite = TRUE)
# 
# ee_drive_to_local(chili_download, dsn = "./predictor_layers/chili.tif", overwrite = TRUE)
# ee_drive_to_local(tpi_download, dsn = "./predictor_layers/tpi.tif", overwrite = TRUE)
# ee_drive_to_local(slope_download, dsn = "./predictor_layers/slope.tif", overwrite = TRUE)


# #---------------------------------------
# # nlcd landcover
# #---------------------------------------
nlcd_rast <- rast("D:/Data/NLCD_landcover_2019_release_all_files_20210604/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img")
bcr_conical <- bcr %>% st_transform(crs(nlcd_rast))
#nlcd is in a non-EPSG CRS: https://spatialreference.org/ref/sr-org/6630/

nlcd_rast <- terra::crop(nlcd_rast, bcr_conical) %>%
  terra::project("EPSG:5070")
# plot(nlcd_rast)

pclass <- function(x, y=c(41)) {
  return( length(which(x %in% y)) / length(x) )
}

decid_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1), 1, 2), others = 0)

prop_decid <- terra::aggregate(decid_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_decid <- prop_decid/16
plot(prop_decid)

terra::writeRaster(prop_decid, "./predictor_layers/prop_decid.tiff", overwrite = TRUE)

con_rast <- terra::classify(nlcd_rast, rcl = matrix(c(43, 1), 1, 2), others = 0)

prop_conifer <- terra::aggregate(con_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_conifer <- prop_conifer/16
plot(prop_conifer)

terra::writeRaster(prop_conifer, "./predictor_layers/prop_conifer.tiff", overwrite = TRUE)

forest_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1,
                                                         42,1,
                                                         43,1), 3, 2), others = 0)
prop_forest <- terra::aggregate(forest_rast, 4, fun = sum)
prop_forest <- prop_forest/16
plot(prop_forest)

terra::writeRaster(prop_forest, "./predictor_layers/prop_forest.tiff", overwrite = TRUE)

grass_shrub_rast <- terra::classify(nlcd_rast, rcl = matrix(c(52, 1,
                                                              71,1,
                                                              81,1), 3, 2), others = 0)
prop_grass_shrub <- terra::aggregate(grass_shrub_rast, 4, fun = sum)
prop_grass_shrub <- prop_grass_shrub/16
plot(prop_grass_shrub)

terra::writeRaster(prop_grass_shrub, "./predictor_layers/prop_grass_shrub.tiff", overwrite = TRUE)



#--------------------------------
# Import predictor data layers
#-----------------------------------------------
#TODO fix units!

bcr_albers <- bcr %>% sf::st_transform(crs = "EPSG:5070")

template_raster <- terra::rasterize(vect(bcr_albers), rast(vect(bcr_albers), resolution = 500))
plot(template_raster)

biomass <- list.files("predictor_layers", "biomass", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic()  %>%
  `NAflag<-`(., value = -1) %>%
  terra::aggregate(x = biomass2, fact = 10, fun = median, na.rm = TRUE) %>% 
  terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)
fhd_rast <- list.files("predictor_layers", "fhd", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1) %>%
  terra::aggregate(fact = 10, median, na.rm = TRUE) %>% #aggregate to 250m cells
  terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)
height_rast <- list.files("predictor_layers", "rh95", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1) %>%
  terra::aggregate(fact = 10, median, na.rm = TRUE) %>% #aggregate to 250m cells
  terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster) %>%
  `/`(1000)
understory <- list.files("predictor_layers", "understory", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1) %>%
  terra::aggregate(fact = 10, median, na.rm = TRUE) %>% #aggregate to 250m cells
  terra::focal(w=15, fun=median, na.policy="only", na.rm=T) %>%
  terra::project(template_raster)%>%
  `/`(1000)

chili_rast <- list.files("predictor_layers", "chili", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
tpi_rast <- list.files("predictor_layers", "tpi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
slope_rast <- list.files("predictor_layers", "slope", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)

tmax_rast <- list.files("predictor_layers", "tmax", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
tmin_rast <- list.files("predictor_layers", "tmin", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pr_rast <- list.files("predictor_layers", "ppt", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
soil_rast <- list.files("predictor_layers", "soil", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
aet_rast <- list.files("predictor_layers", "aet", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
cwd_rast <- list.files("predictor_layers", "cwd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pet_rast <- aet_rast + cwd_rast
vpd_rast <- list.files("predictor_layers", "vpd", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)
pdsi_rast <- list.files("predictor_layers", "pdsi", full.names = TRUE)%>%
  terra::rast() %>%
  terra::project(template_raster)

forest_rast <- terra::rast("./predictor_layers/prop_forest.tiff")%>%
  terra::project(template_raster)
decid_rast <- terra::rast("./predictor_layers/prop_decid.tiff")%>%
  terra::project(template_raster)
grass_shrub_rast <- terra::rast("./predictor_layers/prop_grass_shrub.tiff")%>%
  terra::project(template_raster)
oak_rast <- terra::rast("./predictor_layers/prop_oak_fia.tiff")%>%
  terra::project(template_raster)
spruce_rast <- terra::rast("./predictor_layers/prop_spruce_fia.tiff")%>%
  terra::project(template_raster)
decid_rast_fia <- terra::rast("./predictor_layers/prop_hardwood_fia.tiff")%>%
  terra::project(template_raster)

#------------------------
#layers for gbm predictions

predictor_stack <- c(biomass,
                      understory,
                      fhd_rast,
                      height_rast,
                      
                      aet_rast,
                      pet_rast,
                      cwd_rast,
                      soil_rast,
                      vpd_rast,
                      pdsi_rast,
                      tmax_rast,
                      tmin_rast,
                      pr_rast,
                      tpi_rast,
                      chili_rast,
                      slope_rast,
                      
                      forest_rast,
                      decid_rast,
                      grass_shrub_rast,
                      oak_rast,
                      spruce_rast,
                      decid_rast_fia
)


names(predictor_stack) <- c("biomass",
                            "understory_ratio",
                            "fhd_normal",
                            "height",
                            
                            
                            "aet",
                            "pet",
                            "def",
                            "soil",
                            "vpd",
                            "pdsi",
                            "tmmx",
                            "tmmn",
                            "pr",
                            "tpi",
                            "chili",
                            "slope",
                            
                            
                            "prop_decid",
                            "prop_forest",
                            "prop_grass", #TODO change name
                            "prop_oak",
                            "prop_spruce",
                            "prop_decid_fia")

terra::writeRaster(predictor_stack, 
                   "./predictor_layers/predictor_stack_bcr28.grd", 
                   overwrite = TRUE)
