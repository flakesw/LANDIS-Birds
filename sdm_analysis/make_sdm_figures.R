#make figures for sdms
library("tidyverse")

## make maps ------------------------------------
ecoregions <- terra::rast("C:/Users/Sam/Documents/Research/Project-Southern-Appalachians-2018/Models/LANDIS_Sapps_Ecosystem_Papers/Ecos11_NCLD.tif")

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_albers <- bcr %>% 
  st_transform("EPSG:5070")
# st_transform(crs(ecoregions)) #%>%
# st_crop(ecoregions) #

predictor_stack <- terra::rast("predictor_layers/predictor_stack_bcr28.tif") 
# predictor_stack[[c(2,3,4)]] <- predictor_stack[[c(2,3,4)]]/1000
# predictor_stack[[c(4)]] <- predictor_stack[[c(4)]] #TODO check on this?
# if(crs(predictor_stack) != crs(bcr_albers)) predictor_stack <- project(predictor_stack, bcr_albers)

preds <- rast("./sdm_analysis/outputs/prediction_maps/", species, "_gedi_only.tiff")
# writeRaster(preds, paste0("./sdm_analysis/outputs/prediction_maps/", species, "_gedi_only.tiff"))


#---------------------------------------------------
# partial effects
#---------------------------------------------------

mod1 <- readRDS("./sdm_analysis/sdms/cerw_dist_model_gedi_only.RDS")
mod2 <- readRDS("./sdm_analysis/sdms/gwwa_dist_model_full.RDS")

gbm.plot(mod1, n.plots = 12, plot.layout = c(4,3))
gbm.plot(mod2, n.plots = 12, plot.layout = c(4,3))

#---------------------------------------------------
# Variable importance plot
#---------------------------------------------------
var_tab <- read.csv("./sdm_analysis/variable_table.csv")
relinf1 <- summary(mod1)
cBars <- nrow(relinf1)
relinf2 <- summary(mod2)

par(mar=c(5,8,4,2) + 0.1)
barplot(relinf1[cBars:1, "rel.inf"], 
        horiz = TRUE,
        las = 2,
        col = rainbow(cBars, start = 3/6, end = 4/6), 
        names = relinf1$var[cBars:1],
        xlab = "Relative influence")


relinf_combined <- left_join(relinf1, relinf2, by = c("var")) %>%
 left_join(var_tab) %>%
  rename(CERW = rel.inf.x,
         GWWA = rel.inf.y)

#TODO fix the colors!
par(mar=c(5,7,4,2) + 0.1)
barplot(formula = cbind(relinf_combined$GWWA[cBars:1], relinf_combined$CERW[cBars:1]) ~ relinf_combined$short_name[cBars:1],
        col = relinf_combined$Color[cBars:1],
        horiz = TRUE,
        las = 2,
        beside = TRUE,
        ylab = "",
        xlab = "Relative influence",
        legend.text = c("GWWA", "CERW"),
        args.legend = list(x = "bottomright",
                           inset = c(0.08, 0.08)))

