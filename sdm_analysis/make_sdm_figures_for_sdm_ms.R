#make figures for sdms
library("tidyverse")
library("sf")
library("terra")
library("dismo")
library("gbm")
library("lemon")
library("cowplot")
library("gridExtra")
library("tidyterra")
library("colorspace")
source("./sdm_analysis/r_functions.R")

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank(),
             strip.background = element_rect(fill = "white"))

## import some stuff for mapping ------------------------------------
species <- "cerw"

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

#import predictions just to use as a template
preds <- rast(paste0("./sdm_analysis/outputs/prediction_maps/", species, "_dist_model_full.tiff"))
states <- sf::st_read("C:/Users/Sam/Documents/Maps/Basic maps/state boundaries/cb_2018_us_state_5m/cb_2018_us_state_5m.shp") %>%
  sf::st_transform(crs = crs(preds)) %>%
  st_crop(preds)

#---------------------------------------------------
# partial effects
#---------------------------------------------------

mod1 <- readRDS("./sdm_analysis/sdms/cerw_dist_model_full.RDS")
mod2 <- readRDS("./sdm_analysis/sdms/gwwa_dist_model_full.RDS")

gbm.plot(mod1, n.plots = 24, plot.layout = c(4,3), common.scale = TRUE, write.title = TRUE)
gbm.plot(mod2, n.plots = 24, plot.layout = c(4,3), common.scale = TRUE, write.title = TRUE)

#---------------------------------------------------
# Variable importance plot -- Fig. 3
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

svg(filename = "./sdm_analysis/outputs/fig_3 relative importance.svg")
par(mar=c(5,7,4,2) + 0.1)
barplot(formula = cbind(relinf_combined$GWWA[cBars:1], relinf_combined$CERW[cBars:1]) ~ relinf_combined$short_name[cBars:1],
        col = c("#ebca14", "#00a6e2"),
        horiz = TRUE,
        las = 2,
        beside = TRUE,
        ylab = "",
        xlab = "Relative influence",
        legend.text = c("GWWA", "CERW"),
        args.legend = list(x = "bottomright",
                           inset = c(0.08, 0.08)))

dev.off()

#----------------------------------------------------------------------------
# SDM partial effects and predictions -- Fig. 2
#---------------------------------------------------------------------------
var_tab <- read.csv("./sdm_analysis/variable_table.csv", 
                    encoding = "latin1", 
                    check.names = FALSE)

mod1 <- readRDS("./sdm_analysis/sdms/cerw_dist_model_full.RDS")

vars_plot_data <- as_tibble(summary(mod1)[1:9, ]) %>% #this gets variables in order of rel.inf
  left_join(var_tab, by = c("var" = "var")) #%>%
  # mutate(short_name_units = gsub("/n", "\n", short_name_units))
#get values to make partial effects plots
vars_plot_data$plot_vals <- map(vars_plot_data$var, function(vars) plot.gbm(mod1, vars, return.grid = TRUE))


plot_val_data <- vars_plot_data$plot_vals[c(1:9)] %>%
    bind_rows() %>%
    pivot_longer(cols = !y,
                 values_drop_na = TRUE) %>%
    left_join(var_tab, by = c("name" = "var")) %>%
    mutate(short_name = factor(short_name_units, levels = vars_plot_data$short_name_units)) #reorder by variable importance

# parse_units <- function(x) {x <- gsub("/n", "\n", x)
#                             bquote(x)}


effects_plot_cerw <-   ggplot(plot_val_data, aes(x = value, y = boot::inv.logit(y))) +
                       geom_line() +
                       ylab("Habitat index") +
                       ggtitle("Cerulean warbler") +
                       scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits = c(0.2, 0.8)) + 
                       scale_x_continuous(breaks=scales::pretty_breaks(n=3)) + 
                       facet_wrap(facets = "short_name", 
                                  scales = "free_x", 
                                  strip.position = "bottom",
                                  labeller = as_labeller(function(x) str_wrap(string = x, width = 10)))+ 
                       theme(strip.background = element_blank(),
                              strip.placement='outside',
                             strip.text = element_text(margin = margin(t = 0, b = 0)),
                             axis.title.x=element_blank(),
                             axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                             plot.title = element_text(size = 15, hjust = 0.5, margin = margin(t = 4, b = 4)),
                             plot.margin=unit(c(0,2,4,2), 'points')) 
plot(effects_plot_cerw)
  
  preds <- rast("./sdm_analysis/outputs/prediction_maps/cerw_dist_model_full.tiff")
map_cerw <-   ggplot() +
    geom_spatraster(data = preds) +
    # scale_fill_continuous_divergingx(palette = 'Geyser', mid = 0, rev = TRUE) +
    # scale_fill_whitebox_c(palette = "bl_yl_rd", direction = -1, na.value = "transparent")+
    # scale_fill_viridis_c(option = "magma") +
    # scale_fill_gradientn(colours = terrain.colors(7)) +
    # scale_fill_gradient2(low = "white", mid = "yellow", high = "darkgreen", midpoint = .5) +
    scale_fill_distiller(palette="YlGn", direction = 1, na.value = "transparent") + #distiller for continuous palette
  # scale_fill_terrain_c() +   
  geom_sf(data = bcr_albers, fill = NA, linewidth = 1) +
    geom_sf(data = states, fill = NA, alpha = 0.5) +
  theme(plot.margin=unit(c(0,0,0,0), 'lines'))+
  labs(fill = "Habitat index")


#GWWA
mod1 <- readRDS("./sdm_analysis/sdms/gwwa_dist_model_full.RDS")

vars_plot_data <- as_tibble(summary(mod1)[1:9, ]) %>% #this gets variables in order of rel.inf
  left_join(var_tab, by = c("var" = "var")) #%>%
# mutate(short_name_units = gsub("/n", "\n", short_name_units))
#get values to make partial effects plots
vars_plot_data$plot_vals <- map(vars_plot_data$var, function(vars) plot.gbm(mod1, vars, return.grid = TRUE))


plot_val_data <- vars_plot_data$plot_vals[c(1:9)] %>%
  bind_rows() %>%
  pivot_longer(cols = !y,
               values_drop_na = TRUE) %>%
  left_join(var_tab, by = c("name" = "var")) %>%
  mutate(short_name = factor(short_name_units, levels = vars_plot_data$short_name_units)) #reorder by variable importance

effects_plot_gwwa <-   ggplot(plot_val_data, aes(x = value, y = boot::inv.logit(y))) +
  geom_line() +
  ylab("Habitat index") +
  ggtitle(label = "Golden-winged warbler") +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits = c(0.2, 0.8)) + 
  facet_wrap(facets = "short_name", 
             scales = "free_x", 
             strip.position = "bottom",
             labeller = as_labeller(function(x) str_wrap(string = x, width = 10)))+
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme(strip.background = element_blank(),
        strip.placement='outside',
        strip.text = element_text(margin = margin(t = 0, b = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size = 15, hjust = 0.5, margin = margin(t = 4, b = 4)),
        plot.margin=unit(c(0,2,4,2), 'points'))
plot(effects_plot_gwwa)

preds <- rast("./sdm_analysis/outputs/prediction_maps/gwwa_dist_model_full.tiff")
map_gwwa <-   ggplot() +
  geom_spatraster(data = preds) +
  # scale_fill_terrain_c() + 
  scale_fill_distiller(palette = "YlGn", direction = 1, na.value = NA,
                       limits = c(0,1)) +
  guides(fill = guide_colourbar(nrow = 1, position = "bottom", 
                                direction = "horizontal")) +
  labs(fill = "Habitat index") +
  geom_sf(data = bcr_albers, fill = NA, linewidth = 1) +
  geom_sf(data = states, fill = NA, alpha = 0.5)
  
fill_legend <- get_plot_component(map_gwwa, 'guide-box-bottom', return_all = TRUE)

#combine plots
combined_plot <- plot_grid(effects_plot_cerw,
                  effects_plot_gwwa, 
                  map_cerw+ theme(legend.position="none"), 
                  map_gwwa+ theme(legend.position="none")) %>%
  plot_grid(fill_legend, ncol = 1, rel_heights = c(1, 0.06))
plot(combined_plot)
ggsave(combined_plot, 
       filename = "./sdm_analysis/outputs/fig_2_maps_and_effects.svg", 
       width = 6, height = 8)

#----------------------------------------------
# variable importance ~ spatial scale. Fig. 5
#----------------------------------------------
species_list <- c("cerw", "gwwa")

model_list <- list(map(paste0("./sdm_analysis/sdms/", species_list, "_all_scales_full_model.RDS"),
                               readRDS)) %>%
  bind_rows() %>%
  mutate(species = map(strsplit(names, "_"), pluck(1)) %>% unlist(),
         scale = map(strsplit(names, "_"), pluck(2)) %>% unlist())

vars_to_plot <- c("height", "biomass", "fhd", "understory", "area_short")

# walk(model_list, print)
# walk(model_list, summary)
# walk(model_list$models, function(x) gbm.plot(x, n.plots = 9, plot.layout = c(3,3)))

model_list$model_cv_stats_list <- map(model_list$models, .f = function(x) pluck(x) %>% 
                                      pluck("cv.statistics") %>% 
                                      pluck("discrimination.mean")) %>%
                                  unlist() 

model_list$var_importance_list <- map(model_list$models, .f = function(x) pluck(x) %>% 
                             pluck("contributions"))

model_list$gedi_vars_preds_list <- map(model_list$models, .f = function(mod){map(vars_to_plot, function(vars) 
                                                                  plot.gbm(mod, vars, return.grid = TRUE))})

varimp_list <- model_list %>%
  dplyr::select(species, scale, var_importance_list) %>%
  unnest(cols = var_importance_list) %>%
  filter(var %in% vars_to_plot) %>%
  mutate(scale = factor(scale, levels = c("100", "500", "1000", "5000", "10000", "effort"))) %>%
  mutate(var = factor(var)) %>%
  mutate(var = forcats::fct_recode(var, FHD = "fhd", RatioUnder = "understory", Biomass = "biomass", Ht = "height", PropShort = "area_short"))

varimp_scale <- ggplot(varimp_list, aes(x = scale, y = rel.inf, col = species, group = species)) +
  geom_point() +
  geom_line(data = varimp_list %>% filter(!(scale %in% "effort"))) +
  labs(y = expression(paste("Relative variable influence (percent)")), x = "Scale (m)") +
  scale_color_manual(labels = c("CERW", "GWWA"), values = c("#00a6e2", "#ebca14")) +
  theme(legend.title = element_blank()) +
  facet_wrap(facets = "var") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
varimp_scale <- tag_facet(varimp_scale)
varimp_scale <- shift_legend2(varimp_scale)

plot(varimp_scale)
ggsave(file="./sdm_analysis/outputs/fig_5_variable_importance.svg", plot=varimp_scale, width=5, height=3)



#---------------------------------------------------
# Fig. xx sum of importance of different variable types
#TODO finish this part?
#---------------------------------------------------


model_list$model_cv_stats_list <- map(model_list$models, .f = function(x) pluck(x) %>% 
                                        pluck("cv.statistics") %>% 
                                        pluck("discrimination.mean")) %>%
  unlist() 

# model_list$var_importance_list <- map(model_list$models, .f = function(x) pluck(x) %>% 
#                                         pluck("contributions")) %>%
#   map(left_join())


varimp_list <- model_list %>%
  dplyr::select(species, scale, var_importance_list) %>%
  unnest(cols = var_importance_list) %>%
  # filter(var %in% vars_to_plot) %>%
  mutate(scale = factor(scale, levels = c("100", "500", "1000", "5000", "10000", "effort"))) %>%
  left_join(var_tab, by = c("var" = "var")) %>%
  group_by(species, scale, Type) %>%
  dplyr::summarize(total_importance = sum(rel.inf))

varimp_scale <- ggplot(varimp_list, aes(x = scale, y = total_importance, col = species, group = species)) +
  geom_point() +
  geom_line(data = varimp_list %>% filter(!(scale %in% "effort"))) +
  labs(y = expression(paste("Relative variable influence (percent)")), x = "Scale (m)") +
  scale_color_manual(labels = c("CERW", "GWWA"), values = c("#00a6e2", "#ebca14")) +
  theme(legend.title = element_blank()) +
  facet_wrap(facets = "Type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
varimp_scale <- tag_facet(varimp_scale)
varimp_scale <- shift_legend2(varimp_scale)

plot(varimp_scale)
ggsave(file="./sdm_analysis/outputs/fig_xx_relative influence of predictor types.svg", plot=varimp_scale, width=5, height=3)




#-----------------------------------------------------
# Fig. 4 -- predicted effects with different variable subsets
#-----------------------------------------------------
species <- "gwwa"
model_names  <- list.files("./sdm_analysis/sdms", pattern = species, full.names = TRUE) %>%
  `[`(grep("dist", .))
model_list <- tibble(names = basename(model_names),
                    models = map(model_names, readRDS)) %>%
  mutate(species = map(strsplit(names, "_"), pluck(1)) %>% unlist(),
         type = map(basename(model_names), function(x){str_match(pattern = "model_\\s*(.*?)\\s*.RDS", string = x) %>%
                                                       unlist() %>% 
                                                       pluck(2)}) %>%
                    unlist())

vars_to_plot

#extract predictions for each variable
model_list$gedi_vars_preds_list <- map(model_list$models, 
                                       .f = function(mod){map(vars_to_plot, function(vars){ 
                                                              if(vars %in% mod$gbm.call$predictor.names)
                                                              return(plot.gbm(mod, vars, return.grid = TRUE)) else{
                                                                return(NA)
                                                              }})})

preds_list <- model_list %>%
  dplyr::select(species, type, gedi_vars_preds_list) %>%
  unnest(cols = gedi_vars_preds_list) %>%
  unnest(cols = gedi_vars_preds_list, keep_empty = TRUE) %>%
  pivot_longer(cols = c(height, biomass:area_short),
               names_to = "variable",
               values_to = "x") %>%
  mutate(variable = forcats::fct_recode(variable, FHD = "fhd", RatioUnder = "understory", 
                                        Biomass = "biomass", Ht = "height", PropShort = "area_short")) %>%
  mutate(type = forcats::fct_recode(type, Topoclimate = "clim", "Topoclimate + GEDI" = "clim_gedi", 
                                    All = "full", "GEDI only" = "gedi_only", "All - GEDI" = "no_gedi",
                                    "Landcover + Inventory" = "veg", 
                                    "Landcover + \nInventory + GEDI" = "veg_gedi")) %>%
  filter(type %in% c("Topoclimate + GEDI", "All", "GEDI only", "Landcover + \nInventory + GEDI"))




gedi_effects <- ggplot(preds_list, aes(y= boot::inv.logit(y), x = x, group = type, col = type)) +
  geom_line() +
  facet_wrap(facets = c("variable"), strip.position = "bottom", scales = "free_x") + 
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside',
        axis.title.x = element_blank()) +
  ylab("Habitat Index") +
  ggtitle("Golden-winged warbler")
gedi_effects <- tag_facet(gedi_effects)
gedi_effects <- shift_legend2(gedi_effects)
# plot(gedi_effects)
gedi_effects_gwwa <- gedi_effects


species <- "cerw"
model_names  <- list.files("./sdm_analysis/sdms", pattern = species, full.names = TRUE) %>%
  `[`(grep("dist", .))
model_list <- tibble(names = basename(model_names),
                     models = map(model_names, readRDS)) %>%
  mutate(species = map(strsplit(names, "_"), pluck(1)) %>% unlist(),
         type = map(basename(model_names), function(x){str_match(pattern = "model_\\s*(.*?)\\s*.RDS", string = x) %>%
             unlist() %>% 
             pluck(2)}) %>%
           unlist())

vars_to_plot

#extract predictions for each variable
model_list$gedi_vars_preds_list <- map(model_list$models, 
                                       .f = function(mod){map(vars_to_plot, function(vars){ 
                                         if(vars %in% mod$gbm.call$predictor.names)
                                           return(plot.gbm(mod, vars, return.grid = TRUE)) else{
                                             return(NA)
                                           }})})

preds_list <- model_list %>%
  dplyr::select(species, type, gedi_vars_preds_list) %>%
  unnest(cols = gedi_vars_preds_list) %>%
  unnest(cols = gedi_vars_preds_list, keep_empty = TRUE) %>%
  pivot_longer(cols = c(height, biomass:area_short),
               names_to = "variable",
               values_to = "x") %>%
  mutate(variable = forcats::fct_recode(variable, FHD = "fhd", RatioUnder = "understory", 
                                        Biomass = "biomass", Ht = "height", PropShort = "area_short")) %>%
  mutate(type = forcats::fct_recode(type, Topoclimate = "clim", "Topoclimate + GEDI" = "clim_gedi", 
                                    All = "full", "GEDI only" = "gedi_only", "All - GEDI" = "no_gedi",
                                    "Landcover + Inventory" = "veg", 
                                    "Landcover + \nInventory + GEDI" = "veg_gedi")) %>%
  filter(type %in% c("Topoclimate + GEDI", "All", "GEDI only", "Landcover + \nInventory + GEDI"))




gedi_effects <- ggplot(preds_list, aes(y= boot::inv.logit(y), x = x, group = type, col = type)) +
  geom_line() +
  facet_wrap(facets = c("variable"), strip.position = "bottom", scales = "free_x") + 
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside',
        axis.title.x = element_blank()) +
  ylab("Habitat Index") +
  ggtitle("Cerulean warbler")
gedi_effects <- tag_facet(gedi_effects, tag_pool = letters[-c(1:5)])
gedi_effects <- shift_legend2(gedi_effects)
plot(gedi_effects)
gedi_effects_cerw <- gedi_effects

test <- plot_grid(gedi_effects_gwwa,
                  gedi_effects_cerw,
                  ncol = 1)
plot(test)
ggsave(test, filename = "./sdm_analysis/outputs/fig_4_gedi_effects.svg", dpi = 600, width = 6, height = 6)


#----------------------------------------------------------------------------
# supplementary figure all vars
# TODO start here
#---------------------------------------------------------------------------
var_tab <- read.csv("./sdm_analysis/variable_table.csv")

mod1 <- readRDS("./sdm_analysis/sdms/cerw_dist_model_full.RDS")

vars_plot_data <- as_tibble(summary(mod1)[1:10, ]) %>% #this gets variables in order of rel.inf
  left_join(var_tab, by = c("var" = "var"))
#get values to make partial effects plots
vars_plot_data$plot_vals <- map(vars_plot_data$var, function(vars) plot.gbm(mod1, vars, return.grid = TRUE))


plot_val_data <- vars_plot_data$plot_vals[c(1:9)] %>%
  bind_rows() %>%
  pivot_longer(cols = !y,
               values_drop_na = TRUE) %>%
  left_join(var_tab, by = c("name" = "var")) %>%
  mutate(short_name = factor(short_name, levels = vars_plot_data$short_name)) #reorder by variable importance

effects_plot_cerw <-   ggplot(plot_val_data, aes(x = value, y = boot::inv.logit(y))) +
  geom_line() +
  ylab("Habitat index") +
  facet_wrap(facets = "short_name", scales = "free_x", strip.position = "bottom")+ 
  theme(strip.background = element_blank(),
        strip.placement='outside',
        axis.text.x=element_blank())

preds <- rast("./sdm_analysis/outputs/prediction_maps/cerw_dist_model_full.tiff")
map_cerw <-   ggplot() +
  geom_spatraster(data = preds) +
  # scale_fill_continuous_divergingx(palette = 'Geyser', mid = 0, rev = TRUE) +
  # scale_fill_whitebox_c(palette = "bl_yl_rd", direction = -1, na.value = "transparent")+
  # scale_fill_viridis_c(option = "magma") +
  # scale_fill_gradientn(colours = terrain.colors(7)) +
  # scale_fill_gradient2(low = "white", mid = "yellow", high = "darkgreen", midpoint = .5) +
  # scale_fill_distiller(palette="YlGn", direction = 1,na.value = "transparent") + #distiller for continuous palette
  scale_fill_terrain_c() +   
  geom_sf(data = bcr_albers, fill = NA, linewidth = 1.5) +
  geom_sf(data = states, fill = NA, alpha = 0.5) +
  labs(fill = "Habitat index")


#GWWA
var_tab <- read.csv("./sdm_analysis/variable_table.csv")

mod1 <- readRDS("./sdm_analysis/sdms/gwwa_dist_model_full.RDS")

vars_plot_data <- as_tibble(summary(mod1)[1:10, ]) %>%
  left_join(var_tab, by = c("var" = "var"))

vars_plot_data$plot_vals <- map(vars_plot_data$var, function(vars) plot.gbm(mod1, vars, return.grid = TRUE))

plot_val_data <- vars_plot_data$plot_vals[c(1:9)] %>%
  bind_rows() %>%
  pivot_longer(cols = !y,
               values_drop_na = TRUE) %>%
  left_join(var_tab, by = c("name" = "var")) %>%
  mutate(short_name = factor(short_name, levels = vars_plot_data$short_name)) #reorder by variable importance

effects_plot_gwwa <-   ggplot(plot_val_data, aes(x = value, y = boot::inv.logit(y))) +
  geom_line() +
  ylab("Habitat index") +
  facet_wrap(facets = "short_name", scales = "free_x", strip.position = "bottom")+ 
  theme(strip.background = element_blank(),
        strip.placement='outside',
        axis.text.x=element_blank())

preds <- rast("./sdm_analysis/outputs/prediction_maps/gwwa_dist_model_full.tiff")
map_gwwa <-   ggplot() +
  geom_spatraster(data = preds) +
  scale_fill_terrain_c() + 
  guides(fill = guide_colourbar(nrow = 1, position = "bottom", direction = "horizontal")) +
  labs(fill = "Habitat index") +
  geom_sf(data = bcr_albers, fill = NA, linewidth = 1.5) +
  geom_sf(data = states, fill = NA, alpha = 0.5)


fill_legend <- get_plot_component(map_gwwa, 'guide-box-bottom', return_all = TRUE)

#combine plots
test <- plot_grid(effects_plot_cerw,
                  effects_plot_gwwa, 
                  map_cerw+ theme(legend.position="none"), 
                  map_gwwa+ theme(legend.position="none")) %>%
  plot_grid(fill_legend, ncol = 1, rel_heights = c(1, 0.1))
plot(test)
ggsave(test, filename = "./sdm_analysis/outputs/test.svg", dpi = 600, width = 6, height = 6)

