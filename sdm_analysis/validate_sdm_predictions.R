#compare results to BBS
library("sf")
library("tidyverse")
library("terra")
library("tidyterra")
library("gbm")


bbs_bird_list <- read.csv("validation data/bbs/2024Release_Nor/SpeciesList.csv") %>%
  mutate(English_Common_Name = tolower(English_Common_Name))
route_list <- read.csv("validation data/bbs/2024Release_Nor/Routes.csv") %>%
  mutate(Route = sprintf("%03d", Route),
         StateNum = sprintf("%02d", StateNum),
         unique_id = paste0(StateNum, Route)) %>%
  filter(BCR == 28)

bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326")%>%
  filter(BCR == 28) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_albers <- bcr %>% 
  st_transform("EPSG:5070")


#can find routes here: https://www.mbr-pwrc.usgs.gov/bbs/geographic_information/Instructions_trend_route.htm
# nabbs02_mis_alb/nabbs02_mis_alb.shp
# route_lines/vy474dv5024.shp
routes <- sf::st_read("./validation data/bbs/nabbs02_mis_alb/nabbs02_mis_alb.shp") %>%
  sf::st_set_crs("EPSG:5070") %>%
  st_transform(st_crs(bcr_albers)) %>%
  st_join(bcr_albers, join = st_within) %>%
  mutate(RTENO = sprintf("%05d", RTENO)) %>%
  filter(BCR == 28)
# plot(st_geometry(routes))

routes_buffer <- st_buffer(routes, 100) %>% 
  group_by(RTENO) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup() %>%
  st_as_sf()

bbs_data_all <- list.files("validation data/bbs/2024Release_Nor/50-StopData", full.names = TRUE) %>% 
  read_csv()
bbs_n_years <- bbs_data_all %>%
  mutate(unique_id = paste0(StateNum, Route)) %>%
  filter(Year %in% c(2000:2020),
         unique_id %in% route_list$unique_id) %>%
  group_by(Year, unique_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(unique_id) %>%
  summarize(n_years = n())



#---------------------------------------------
# Extract data for species
species_name_long <- "cerulean warbler"
species_abbr <- "cerw"


AOU_num <- bbs_bird_list[bbs_bird_list$English_Common_Name == species_name_long, ]$AOU

bbs_data_sp <- bbs_data_all %>%
  mutate(AOU = str_remove(AOU, "^0+")) %>%
  mutate(unique_id = paste0(StateNum, Route)) %>%
  filter(unique_id %in% routes_buffer$RTENO) %>%
  filter(Year %in% c(2010:2020)) %>%
  filter(AOU %in% AOU_num) %>%
  rowwise() %>% 
  mutate(across(.cols = Stop1:Stop50, .fns = function(x) x > 0)) %>%
  mutate(route_year_total = sum(c_across(Stop1:Stop50), na.rm = T)) %>%
  ungroup()

bbs_data_summary <- bbs_data_sp %>%
  group_by(unique_id) %>%
  summarize(route_total = sum(route_year_total),
            route_mean = route_total / 50) %>%
  left_join(bbs_n_years, by = "unique_id") %>% 
  mutate(prop_overall = route_mean / n_years)


preds_list <- list.files("sdm_analysis/outputs/prediction_maps/", full.names = TRUE) %>%
  `[`(grepl(species_abbr, .))
mod_name <- NA
cor_val  <- NA
auc <- NA
sensitivity <- NA
specificity <- NA
atlas_cor <- NA
ebd_cor <- NA


random_pts <- sf::st_sample(st_buffer(bcr_albers, -10000), 10000) %>% st_as_sf()

atlas_crs <- sf::st_read("./validation data/atlas/gwwa/BirdAtlas_6420.kml", layer = "Golden-Winged Warbler") %>%
  st_crs()
atlas_current <- terra::rast("./validation data/atlas/gwwa/models/cur.png")[[1]]
ext(atlas_current) <- c(-112.993023, -58.990946, 19.066954, 54.427486)
crs(atlas_current) <- atlas_crs$wkt
atlas_current <- (255 - atlas_current)/255
atlas_current <- terra::project(atlas_current, preds) %>%
  terra::crop(bcr_albers, mask = TRUE)
random_pts$atlas <- terra::extract(atlas_current, vect(random_pts))$cur_1

ebd_st <- terra::rast("./validation data/status_and_trends/gowwar_abundance_seasonal_breeding_max_2022.tif") %>%
  terra::project(preds) %>%
  terra::crop(vect(bcr_albers), mask = TRUE) %>%
  terra::clamp(upper = 1)
random_pts$ebd_st <- terra::extract(ebd_st, vect(random_pts))$breeding

bbs_threshold <- quantile(routes_preds$prop_overall, 0.75)
preds_threshold <- quantile(routes_preds$lyr1, 0.75)

for(i in 1:length(preds_list)){
  preds <- rast(preds_list[i])
  routes_preds <- terra::extract(preds, vect(routes_buffer), fun = mean, ID = FALSE, bind = TRUE) %>% 
    st_as_sf() %>%
    left_join(bbs_data_summary, by = c("RTENO" = "unique_id"))
  routes_preds <- routes_preds %>%
    mutate(prop_overall = ifelse(is.na(prop_overall), 0, prop_overall))
  
  mod_name[i] <- preds_list[i]
  cor_val[i]  <- cor.test(routes_preds$prop_overall, routes_preds$lyr1, method = "spearman")$estimate
  auc[i] <- pROC::auc(pROC::roc(response = routes_preds$prop_overall > bbs_threshold, predictor = routes_preds$lyr1, auc = TRUE))
  
  conMat <- caret::confusionMatrix(as.factor(routes_preds$lyr1 > preds_threshold),
                         as.factor(routes_preds$prop_overall > bbs_threshold),
                         positive = "TRUE")
  sensitivity[i] <- conMat$byClass["Sensitivity"]
  specificity[i] <- conMat$byClass["Specificity"]
  
  
  random_pts$preds <- terra::extract(preds, vect(random_pts))$lyr1
  atlas_cor[i] <- cor.test(random_pts$preds, random_pts$atlas, method = "spearman")$estimate
  ebd_cor[i] <- cor.test(random_pts$preds, random_pts$ebd_st, method = "spearman")$estimate
  
}

validation_stats <- data.frame(mod_name, cor_val, auc, sensitivity, specificity, atlas_cor, ebd_cor)
validation_stats$tss <- validation_stats$sensitivity + validation_stats$specificity - 1
validation_stats



#--------------------------
# Make maps of predictions and validation layers
states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)



preds_list <- list.files("sdm_analysis/outputs/prediction_maps/", full.names = TRUE) %>%
  `[`(grepl("cerw", .))
preds <- rast(preds_list[3])
pred_plot <- ggplot() +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index") 

atlas_crs <- sf::st_read("./validation data/atlas/cerw/BirdAtlas_6580.kml", layer = "Cerulean Warbler") %>%
  st_crs()
atlas_current <- terra::rast("./validation data/atlas/cerw/models/cur.png")[[1]]
ext(atlas_current) <- c(-112.993023, -58.990946, 19.066954, 54.427486)
crs(atlas_current) <- atlas_crs$wkt
atlas_current <- (255 - atlas_current)/255
atlas_current <- terra::project(atlas_current, preds) %>%
  terra::crop(bcr_albers) %>%
  terra::mask(bcr_albers)

atlas_plot <- ggplot() +
  geom_spatraster(data = atlas_current) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")


ebd_st <- terra::rast("./validation data/status_and_trends/cerwar_abundance_seasonal_breeding_max_2022.tif") %>%
  terra::project(preds) %>%
  terra::crop(vect(bcr_albers), mask = TRUE) %>%
  terra::clamp(upper = 1, values = TRUE)


# plot(ebd_st)
ebd_st_plot <- ggplot() +
  geom_spatraster(data = ebd_st) +
  scale_fill_terrain_c() +
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")


title <- ggdraw() + 
  draw_label(
    "Cerulean warbler",
    fontface = 'bold',
    x = 0,
    hjust = 0)
  # theme(
  #   # add margin on the left of the drawing canvas,
  #   # so title is aligned with left edge of first plot
  #   plot.margin = margin(0, 0, 0, 7)
  # )

#combine plots
all_products_cerw <- plot_grid(pred_plot + theme(legend.position="none"), 
                          atlas_plot + theme(legend.position="none"), 
                          ebd_st_plot + theme(legend.position="none"), 
                  nrow = 1, ncol = 3,
                  labels = c("Current study", "CC Bird Atlas", "eBird S&T"),
                  label_size = 12) %>%
  plot_grid(title, ., ncol = 1, rel_heights = c(0.1,1))


#-----------
# same for gwwa
#------------
preds_list <- list.files("sdm_analysis/outputs/prediction_maps/", full.names = TRUE) %>%
  `[`(grepl("gwwa", .))
preds <- rast(preds_list[3])
pred_plot <- ggplot() +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index") 
# geom_sf(ebd, mapping = aes(color = species_observed)) #TODO add to appendix!

atlas_crs <- sf::st_read("./validation data/atlas/gwwa/BirdAtlas_6420.kml", layer = "Golden-Winged Warbler") %>%
  st_crs()
atlas_current <- terra::rast("./validation data/atlas/gwwa/models/cur.png")[[1]]
ext(atlas_current) <- c(-112.993023, -58.990946, 19.066954, 54.427486)
crs(atlas_current) <- atlas_crs$wkt
atlas_current <- (255 - atlas_current)/255
atlas_current <- terra::project(atlas_current, preds) %>%
  terra::crop(bcr_albers) %>%
  terra::mask(bcr_albers)

atlas_plot <- ggplot() +
  geom_spatraster(data = atlas_current) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")


ebd_st <- terra::rast("./validation data/status_and_trends/gowwar_abundance_seasonal_breeding_max_2022.tif") %>%
  terra::project(preds) %>%
  terra::crop(vect(bcr_albers), mask = TRUE) %>%
  terra::clamp(upper = 1, values = TRUE)


# plot(ebd_st)
ebd_st_plot <- ggplot() +
  geom_spatraster(data = ebd_st) +
  scale_fill_terrain_c() +
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")

ebd_legend <- ggplot() +
  geom_spatraster(data = ebd_st) +
  scale_fill_terrain_c()+
  # theme(legend.position = "bottom") +
  # theme(legend.direction = "horizontal") +
  guides(fill = guide_colourbar(nrow = 1, position = "bottom", direction = "horizontal")) +
  labs(fill = "Habitat index")
plot(ebd_legend)


#combine plots
all_products <- plot_grid(pred_plot + theme(legend.position="none"), 
                          atlas_plot + theme(legend.position="none"), 
                          ebd_st_plot + theme(legend.position="none"), 
                          nrow = 1, ncol = 3,
                          labels = c("Current study", "CC Bird Atlas", "eBird S&T"),
                          label_size = 12)

# extract a legend that is laid out horizontally
legend_b <- get_plot_component(ebd_legend, 'guide-box-bottom', return_all = TRUE)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
title <- ggdraw() + 
  draw_label(
  "Golden-winged warbler",
  fontface = 'bold',
  x = 0,
  hjust = 0)

all_products_gwwa <- plot_grid(title, all_products, ncol = 1, rel_heights = c(0.1, 1))
plot(all_products_gwwa)


all_products_combined <- plot_grid(all_products_cerw, all_products_gwwa, legend_b, 
                                   ncol = 1, rel_heights = c(1,1,0.15))

plot(all_products_combined)

ggsave(all_products_combined, filename = "./sdm_analysis/outputs/figure_compare_model_products.svg", 
       dpi = 600, width = 7, height = 7)
