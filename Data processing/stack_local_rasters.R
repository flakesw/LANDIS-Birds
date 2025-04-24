buffer_scale <- 10000
if(buffer_scale == "effort"){ 
  aggregate_scale <- 1
} else{
  aggregate_scale <- buffer_scale/1000
  aggregate_scale <- ifelse(aggregate_scale < 1, 1, aggregate_scale)
}

gedi_stack_orig <- terra::rast(c("./predictor_layers/sdm/height_raw.tif",
                                 "./predictor_layers/sdm/biomass_raw.tif",
                                 "./predictor_layers/sdm/fhd_raw.tif",
                                 "./predictor_layers/sdm/prop_open_raw.tif",
                                 "./predictor_layers/sdm/understory_raw.tif"))
NAflag(gedi_stack_orig) <- -1
names(gedi_stack_orig) <- c("height", "biomass", "fhd", "area_short", "understory")

# gedi_stack <- terra::aggregate(gedi_stack_orig, fact = aggregate_scale, fun = "mean", na.rm = TRUE)
# writeRaster(gedi_stack_orig, "./predictor_layers/sdm/gedi_stack.tif", overwrite = TRUE)


nlcd <- rast("./predictor_layers/nlcd_cropped_bcr28_epsg5070.tif") %>%
  aggregate(fact = aggregate_scale, fun = "modal", na.rm = TRUE)
plot(nlcd)
nlcd_int <- as.numeric(nlcd)
writeRaster(nlcd_int, 
            paste0("./predictor_layers/sdm/nlcd_rescale_", buffer_scale, ".tif"),
            datatype = "INT8U",
            overwrite = TRUE)


prop_oak <- terra::rast("./predictor_layers/prop_oak_fia.tiff") %>%
  subst(from = NA, to = 0)
prop_spruce <- terra::rast("./predictor_layers/prop_spruce_fia.tiff")%>%
  subst(from = NA, to = 0)
prop_hardwood <- terra::rast("./predictor_layers/prop_hardwood_fia.tiff")%>%
  subst(from = NA, to = 0)

treemap_stack <- c(prop_oak, prop_spruce, prop_hardwood)
names(treemap_stack) <- c("prop_oak", "prop_spruce", "prop_hardwood")

treemap_scale <- treemap_stack %>%
  aggregate(fact = aggregate_scale, fun = "mean", na.rm = TRUE)
plot(treemap_scale)
writeRaster(treemap_scale, paste0("./predictor_layers/sdm/treemap_rescale_", buffer_scale, ".tif"), overwrite = TRUE)
# writeRaster(treemap_stack, "./predictor_layers/sdm/treemap_stack.tif", overwrite = TRUE)
