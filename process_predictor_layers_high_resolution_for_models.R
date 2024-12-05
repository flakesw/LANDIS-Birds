#TODO fix units!
library(tidyverse)
library(terra)
library(sf)


bcr_albers <- bcr %>% sf::st_transform(crs = "EPSG:5070")
template <- list.files("predictor_layers", "biomass", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic()


#make urban and water mask
nlcd_rast <- terra::rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif") %>% as.numeric()
water_urban_mask <- terra::classify(nlcd_rast, rcl = matrix(c(23, 1,
                                                              24,1,
                                                              11,1), 3, 2, byrow = TRUE), others = 0)
water_urban_mask <- terra::project(water_urban_mask, template, mask = TRUE)

#aboveground biomass
biomass_masked <- list.files("predictor_layers", "biomass", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic()  %>%
  `NAflag<-`(., value = -1) %>% 
  terra::mask(water_urban_mask,  maskvalue = 1, updatevalue = -1)
#output for getting bird point-level data
writeRaster(biomass_masked, "./predictor_layers/sdm/biomass_raw.tif", overwrite = TRUE)

fhd_rast_masked <- list.files("predictor_layers", "fhd", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1)  %>%
  terra::mask(water_urban_mask,  maskvalue = 1, updatevalue = -1)
#output for getting bird point-level data
writeRaster(fhd_masked, "./predictor_layers/sdm/fhd_raw.tif", overwrite = TRUE)

height_rast_masked <- list.files("predictor_layers", "rh98", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1) %>%
  terra::mask(water_urban_mask,  maskvalue = 1, updatevalue = -1)
#output for getting bird point-level data
writeRaster(height_rast_masked, "./predictor_layers/sdm/height_raw.tif", overwrite = TRUE)

open_area_rast_masked <- list.files("predictor_layers", "rh98", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1)%>%
  `/`(1000) %>%
  terra::focal(w = 9, fun = function(x){ x <- x[!is.na(x)]
                                          ncells <- length(x)
                                          if(ncells == 0 | is.na(ncells)){
                                            return(NA)
                                          }else{
                                            cells_less_than_5 <- sum(x < 5)
                                            return(cells_less_than_5 / ncells)}
                                          }, 
                na.policy = "all")  %>%
  terra::mask(water_urban_mask,  maskvalue = 1, updatevalue = -1)
#output for getting bird point-level data
writeRaster(open_area_rast_masked, "./predictor_layers/sdm/prop_open_raw.tif", overwrite = TRUE)

understory_masked <- list.files("predictor_layers", "understory", full.names = TRUE) %>%
  lapply(terra::rast) %>%
  terra::sprc() %>%
  terra::mosaic() %>%
  `NAflag<-`(., value = -1) %>%
  terra::mask(water_urban_mask,  maskvalue = 1, updatevalue = -1)
#output for getting bird point-level data
writeRaster(understory_masked, "./predictor_layers/sdm/understory_raw.tif", overwrite = TRUE)
  
