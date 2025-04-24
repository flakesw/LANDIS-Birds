# process NLCD landcover data to make map predictors
# #---------------------------------------
# # nlcd landcover
# #---------------------------------------

# nlcd_rast <- rast("D:/Data/NLCD_landcover_2019_release_all_files_20210604/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img")
# bcr_conical <- bcr %>% st_transform(crs(nlcd_rast))
# #nlcd is in a non-EPSG CRS: https://spatialreference.org/ref/sr-org/6630/
# 
# nlcd_rast <- terra::crop(nlcd_rast, bcr_conical) %>%
#   terra::project("EPSG:5070")
# # plot(nlcd_rast)

#already cropped and reprojected
nlcd_rast <- terra::rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif") %>% as.numeric()

decid_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1), 1, 2, byrow = TRUE), others = 0)
prop_decid <- terra::aggregate(decid_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_decid <- prop_decid/16
plot(prop_decid)
terra::writeRaster(prop_decid, "./predictor_layers/prop_decid.tiff", overwrite = TRUE)

con_rast <- terra::classify(nlcd_rast, rcl = matrix(c(42, 1), 1, 2, byrow = TRUE), others = 0)
prop_conifer <- terra::aggregate(con_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_conifer <- prop_conifer/16
plot(prop_conifer)
terra::writeRaster(prop_conifer, "./predictor_layers/prop_conifer.tiff", overwrite = TRUE)

forest_rast <- terra::classify(nlcd_rast, rcl = matrix(c(41, 1,
                                                         42,1,
                                                         43,1,
                                                         90, 1), 4, 2, byrow = TRUE), others = 0)
prop_forest <- terra::aggregate(forest_rast, 4, fun = sum)
prop_forest <- prop_forest/16
plot(prop_forest)

terra::writeRaster(prop_forest, "./predictor_layers/prop_forest.tiff", overwrite = TRUE)

grass_rast <- terra::classify(nlcd_rast, rcl = matrix(c( 71,1,
                                                         81,1), 2, 2, byrow = TRUE), others = 0)
prop_grass <- terra::aggregate(grass_rast, 4, fun = sum)
prop_grass<- prop_grass/16
plot(prop_grass)

terra::writeRaster(prop_grass, "./predictor_layers/prop_grass.tiff", overwrite = TRUE)

shrub_rast <- terra::classify(nlcd_rast, rcl = matrix(c(52, 1), 1, 2, byrow = TRUE), others = 0)
prop_shrub <- terra::aggregate(shrub_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_shrub <- prop_shrub/16
plot(prop_shrub)
terra::writeRaster(prop_shrub, "./predictor_layers/prop_shrub.tiff", overwrite = TRUE)

water_rast <- terra::classify(nlcd_rast, rcl = matrix(c(11, 1), 1, 2, byrow = TRUE), others = 0)
prop_water <- terra::aggregate(water_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_water <- prop_water/16
plot(prop_water)
terra::writeRaster(prop_water, "./predictor_layers/prop_water.tiff", overwrite = TRUE)

dev_light_rast <- terra::classify(nlcd_rast, rcl = matrix(c(21, 1,
                                                            22,1,
                                                            82,1), 3, 2, byrow = TRUE), others = 0)
prop_dev_light <- terra::aggregate(dev_light_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_dev_light <- prop_dev_light/16
plot(prop_dev_light)
terra::writeRaster(prop_dev_light, "./predictor_layers/prop_dev_light.tiff", overwrite = TRUE)

dev_heavy_rast <- terra::classify(nlcd_rast, rcl = matrix(c(23, 1,
                                                            24, 1), 2, 2, byrow = TRUE), others = 0)
prop_dev_heavy <- terra::aggregate(dev_heavy_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_dev_heavy <- prop_dev_heavy/16
plot(prop_dev_heavy)
terra::writeRaster(prop_dev_heavy, "./predictor_layers/prop_dev_heavy.tiff", overwrite = TRUE)

open_rast <- terra::classify(nlcd_rast, rcl = matrix(c(12,21,31,51,52,71:74,81,82,95, rep(1, 12)),
                                                     12, 2), others = 0)
prop_open <- terra::aggregate(open_rast, 4, fun = sum) #aggregate to 120m; TODO change this to a specified scale
prop_open <- prop_open/16
plot(prop_open)
terra::writeRaster(prop_open, "./predictor_layers/prop_open.tiff", overwrite = TRUE)

