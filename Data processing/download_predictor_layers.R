# #-------------------------------------------------------------------------------
# # Get predictor map surfaces
# #-------------------------------------------------------------------------------
# #-----
# # getting surfaces for prediction

# updated 2024-11-19 to just download, will process in a different script

# TODO try doing everything as INTs -- the understory raster is a short and it's 
# only 2 tiles instead of 6 when downloading

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
    return(image$select("rh98")$multiply(1000)$toShort()$
             updateMask(image$select('quality_flag')$eq(1))#$
             # updateMask(image$select('degrade_flag')$eq(0))
             )
  })$
  median()$
  unmask(-1) #set the NA value to -1 so we can fix it later

height_download <- ee_image_to_drive(height_red,
                                     description = "rh98_raster",
                                     fileNamePrefix = "rh98_raster",
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
             updateMask(image$select('l2b_quality_flag')$eq(1))
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
            updateMask(image$select('algorithm_run_flag')$eq(1))#$
            # updateMask(image$select('degrade_flag')$eq(0))
           )
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
biomass_download_4a$start() 
ee_monitoring(biomass_download_4a)

#--------------
# proportion understory
heightBands = c(paste0("pavd_z", seq(0, 25, by = 5)))

gedi2b = ee$ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')$
  filterBounds(bcr_ee$geometry())$
  map(function(image) {
    return(image$
             updateMask(image$select('algorithmrun_flag')$eq(1))$
             updateMask(image$select('l2b_quality_flag')$eq(1))$
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


# #----------------
# # npp layer
#TODO reduce number of variables before exporting
#TODO get year-of-observation weather (or week of observation?)

modis <- ee$ImageCollection("MODIS/061/MOD17A3HGF")

endyear <- ee$Number(2023)
lag_years <- ee$Number(10) #TODO try different time lags!
startyear <- endyear$subtract(lag_years)

npp_mean <- modis$filter(ee$Filter$calendarRange(startyear, endyear, "year"))$
  select('Npp')$
  filterBounds(bcr_ee)$
  mean()

Map$addLayer(npp_mean)

npp_download <- ee_image_to_drive(npp_mean,
                                   description = "NPP_raster",
                                   fileNamePrefix = "NPP_raster",
                                   folder = "prediction_layers",
                                   region = bcr_ee$geometry(),
                                   scale = 500,
                                   crs = proj,
                                   maxPixels = 4097260648,
                                   fileFormat = 'GeoTIFF')
npp_download$start()
ee_monitoring(npp_download)


#------------------------------
#Download all the rasters
#---------------------------------

#retrieve files from Google Drive. Check if they're ready first!

# ee_drive_to_local(height_download, dsn = "./predictor_layers/rh98.tif", overwrite = TRUE)
# ee_drive_to_local(fhd_download, dsn = "./predictor_layers/fhd.tif", overwrite = TRUE) #TODO
# ee_drive_to_local(biomass_download_4a, dsn = "./predictor_layers/biomass.tif", overwrite = TRUE)
# ee_drive_to_local(under_download, dsn = "./predictor_layers/understory.tif", overwrite = TRUE)
 
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

# ee_drive_to_local(npp_download, dsn = "./predictor_layers/npp.tif", overwrite = TRUE)
