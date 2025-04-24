#make treemap surface for getting predictors and making predictions to

library("rFIA")
library("tidyverse")
library("terra")
treemap <- terra::rast("D:/Data/treemap/RDS-2019-0026_Data/Data/national_c2014_tree_list.tif") 



bcr <- sf::st_read("./maps/BCR_NA.shp") %>%
  sf::st_transform("EPSG:4326") %>%
  filter(BCR == 28) %>%
  # filter(BCR %in% c(12,13,14,22,23,24,27,28,29,30)) %>%
  sf::st_make_valid() %>%
  dplyr::select("BCR")

bcr_treemap <- bcr %>% sf::st_transform(crs(treemap))
tm2 <- treemap %>%
  terra::crop(vect(bcr_treemap), mask = TRUE)

sp_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_SPECIES.csv")

treemap_plots <- unique(tm2, na.rm = TRUE)$tl_id

tl_plots <- read.csv("D:/Data/treemap/RDS-2019-0026_Data/Data/TL_CN_Lookup.txt") %>%
  # filter(tl_id %in% ebird_treemap_extract$tl_id)
  filter(tl_id %in% treemap_plots)
tl_trees <- read.csv("D:/Data/treemap/RDS-2019-0026_Data/Data/Tree_table_CONUS.txt") %>%
  # filter(tl_id %in% ebird_treemap_extract$tl_id)
  filter(tl_id %in% treemap_plots)
sp_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_SPECIES.csv")

fia_trees <- purrr::map(.x = unique(tl_trees$State_Abbreviation),
            .f = function(x) {
                  readFIA(dir = 'D:/Data/fia/rFIA_downloads/',
                          tables = c('TREE'),
                          states = x) %>%
                  '[['('TREE') %>%
                  dplyr::filter(PLT_CN %in% tl_plots$CN) %>%
                  dplyr::select(PLT_CN, SPCD, DRYBIO_AG)},
            .progress = TRUE) %>%
  dplyr::bind_rows()

tree_hardwood <- fia_trees %>%
  left_join(sp_ref, by = "SPCD") %>%
  group_by(PLT_CN) %>%
  mutate(site_bio = sum(DRYBIO_AG, na.rm = TRUE)) %>%
  group_by(PLT_CN, SFTWD_HRDWD) %>%
  summarize(bio = sum(DRYBIO_AG, na.rm = TRUE),
            site_bio = site_bio[1]) %>%
  ungroup() %>%
  mutate(prop = bio/site_bio) %>%
  dplyr::select(PLT_CN, SFTWD_HRDWD, prop, site_bio) %>%
  pivot_wider(names_from = SFTWD_HRDWD, values_from = prop, values_fill = 0) %>%
  mutate(H = ifelse(site_bio == 0 | is.na(site_bio), NA, H),
         S = ifelse(site_bio == 0 | is.na(site_bio), NA, S)) #if there's no biomass, 
                                                            # make proportion = NA


tree_oak_hickory <- fia_trees %>%
  left_join(sp_ref, by = "SPCD") %>%
  mutate(oak_hickory = ifelse(GENUS %in% c("Quercus", "Carya"), "OH", "Other")) %>%
  group_by(PLT_CN) %>%
  mutate(site_bio = sum(DRYBIO_AG, na.rm = TRUE)) %>%
  group_by(PLT_CN, oak_hickory) %>%
  summarize(bio = sum(DRYBIO_AG, na.rm = TRUE),
            site_bio = site_bio[1]) %>%
  group_by(PLT_CN) %>%
  mutate(prop = bio/site_bio) %>%
  dplyr::select(PLT_CN, oak_hickory, prop, site_bio) %>%
  pivot_wider(names_from = oak_hickory, values_from = prop, values_fill = 0) %>%
  mutate(OH = ifelse(site_bio == 0, NA, OH)) #if there's no biomass, make proportion oak/hickory = NA

tree_spruce_fir <- fia_trees %>%
  left_join(sp_ref, by = "SPCD") %>%
  mutate(spruce_fir = ifelse(GENUS %in% c("Abies", "Picea"), "SF", "Other")) %>%
  group_by(PLT_CN) %>%
  mutate(site_bio = sum(DRYBIO_AG, na.rm = TRUE)) %>%
  group_by(PLT_CN, spruce_fir) %>%
  summarize(bio = sum(DRYBIO_AG, na.rm = TRUE),
            site_bio = site_bio[1]) %>%
  group_by(PLT_CN) %>%
  mutate(prop = bio/site_bio) %>%
  dplyr::select(PLT_CN, spruce_fir, prop, site_bio) %>%
  pivot_wider(names_from = spruce_fir, values_from = prop, values_fill = 0) %>%
  mutate(SF = ifelse(site_bio == 0, NA, SF)) #if there's no biomass, make proportion oak/hickory = NA

tl_plots_add <- tl_plots %>%
  left_join(tree_hardwood, by = c("CN" = "PLT_CN")) %>%
  left_join(tree_oak_hickory, by = c("CN" = "PLT_CN"))  %>%
  left_join(tree_spruce_fir, by = c("CN" = "PLT_CN"))

#this takes a really really long time
fia_prop_hardwood_rast <- terra::classify(tm2, tl_plots_add[, c("tl_id", "H")])
# fia_prop_softwood_rast <- terra::classify(tm2, tl_plots_add[, c("tl_id", "S")])
fia_prop_oak_rast <- terra::classify(tm2, tl_plots_add[, c("tl_id", "OH")])
fia_prop_spruce_rast <- terra::classify(tm2, tl_plots_add[, c("tl_id", "SF")])


terra::writeRaster(fia_prop_hardwood_rast, "./predictor_layers/prop_hardwood_fia.tiff")
# terra::writeRaster(fia_prop_softwood_rast, "./predictor_layers/prop_softwood_fia.tiff")
terra::writeRaster(fia_prop_oak_rast, "./predictor_layers/prop_oak_fia.tiff")
terra::writeRaster(fia_prop_spruce_rast, "./predictor_layers/prop_spruce_fia.tiff")
