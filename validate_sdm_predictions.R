#compare results to BBS
library("sf")
library("tidyverse")
library("terra")
library("tidyterra")

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

routes_buffer <- st_buffer(routes, 400) %>% 
  group_by(RTENO) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup() %>%
  st_as_sf()

bbs_data_all <- list.files("validation data/bbs/2024Release_Nor/50-StopData", full.names = TRUE) %>% 
  read_csv()
bbs_n_years <- bbs_data_all %>%
  mutate(unique_id = paste0(StateNum, Route)) %>%
  filter(Year %in% c(2000:2019),
         unique_id %in% route_list$unique_id) %>%
  group_by(Year, unique_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(unique_id) %>%
  summarize(n_years = n())



#---------------------------------------------
# Extract data for species

AOU_num <- bbs_bird_list[bbs_bird_list$English_Common_Name == "cerulean warbler", ]$AOU

bbs_data_sp <- bbs_data_all %>%
  mutate(AOU = str_remove(AOU, "^0+")) %>%
  mutate(unique_id = paste0(StateNum, Route)) %>%
  filter(unique_id %in% routes_buffer$RTENO) %>%
  filter(Year %in% c(2000:2020)) %>%
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


preds <- rast("sdm_analysis/outputs/prediction_maps/cerw_gedi_only.tiff") #TODO make automatic
routes_preds <- terra::extract(preds, vect(routes_buffer), fun = mean, ID = FALSE, bind = TRUE) %>% 
  st_as_sf() %>%
  left_join(bbs_data_summary, by = c("RTENO" = "unique_id"))
routes_preds <- routes_preds %>%
  mutate(prop_overall = ifelse(is.na(prop_overall), 0, prop_overall))

ggplot() +
    geom_spatraster(data = preds) +
    geom_sf(data = routes_preds, aes(col = prop_overall))


plot(routes_preds$prop_overall ~ routes_preds$lyr1,
     xlab = "BRT predicted encounter rate",
     ylab = "BBS presence")       
boxplot(routes_preds$prop_overall ~ routes_preds$lyr1 > 0.5,
        xlab = "Predicted > 50% encounter rate",
        ylab = "BBS encounter rate",
        main = "CERW")
cor.test(routes_preds$prop_overall, routes_preds$lyr1, method = "spearman")
plot(pROC::roc(response = routes_preds$prop_overall > 0.5, predictor = routes_preds$lyr1, auc = TRUE))
auc(pROC::roc(response = routes_preds$prop_overall > 0.5, predictor = routes_preds$lyr1, auc = TRUE))


#--------------------------
# Make maps of predictions and validation layers
states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)
ebd <- read.csv(paste0("./ebird/", species, "_subsampled_balanced_bcr28.csv")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = crs(preds))
ebd$predicted <- terra::extract(preds, vect(ebd))$lyr1
boxplot(ebd$predicted ~ ebd$species_observed)

ggplot() +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index") #+
# geom_sf(ebd, mapping = aes(color = species_observed))


atlas_crs <- sf::st_read("./validation data/atlas/cerw/BirdAtlas_6580.kml", layer = "Cerulean Warbler") %>%
  st_crs()
atlas_current <- terra::rast("./validation data/atlas/cerw/models/cur.png")[[1]]
ext(atlas_current) <- c(-112.993023, -58.990946, 19.066954, 54.427486)
crs(atlas_current) <- atlas_crs$wkt
atlas_current <- (255 - atlas_current)/255
atlas_current <- terra::project(atlas_current, preds) %>%
  terra::crop(bcr_albers, mask = TRUE)
# plot(atlas_current)

ggplot() +
  geom_spatraster(data = atlas_current) +
  scale_fill_terrain_c() + 
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")


ebd_st <- terra::rast("./validation data/status_and_trends/cerwar_abundance_seasonal_breeding_max_2022.tif") %>%
  terra::project(preds) %>%
  terra::crop(vect(bcr_albers), mask = TRUE) %>%
  terra::clamp(upper = 1)
# `/`(max(.[], na.rm = TRUE))
# terra::classify(rcl = c(0, 1, Inf))

# plot(ebd_st)
ggplot() +
  geom_spatraster(data = ebd_st*2.5) +
  scale_fill_terrain_c() +
  geom_sf(data = bcr_albers, fill = NA) +
  geom_sf(data = states, colour = alpha("black",0.5), fill = NA) +
  labs(fill = "Habitat index")

# ebd_st <- ebd_st * (0.29/0.075)
atlas_current <- atlas_current * (mean(preds[], na.rm = TRUE) / mean(atlas_current[], na.rm = TRUE))

error <- ebd_st - preds
error <- atlas_current - preds
plot(error)
hist(error)


