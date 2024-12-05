### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
library("terra")
library("sf")
library("tidyterra")
#--------------------------

species_list <- c("cerw", "gwwa")


#some helper functions
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(file_name = basename(flnm))
}

models_all <- tibble(names = character(length = 0L),
                    models = vector(mode = "list", length = 0L))

for(species in species_list){
  data_files <- list.files("./environment_vars/", pattern = species, full.names = TRUE) %>%
    `[`(grep("combined", .)) 
  data_type <- data.frame(file_name = basename(data_files),
                          file_location = data_files) %>%
    mutate(species = map(strsplit(file_name, "[[:punct:]]"), pluck(1)) %>% unlist(),
           data_source = map(strsplit(file_name, "[[:punct:]]"), pluck(5)) %>% unlist(),
           scale = map(strsplit(file_name, "[[:punct:]]"), pluck(6)) %>% unlist(),
           unique = paste(species, scale, sep = "_"))
  
  #are there any duplicate files or species/scale combos? If either of these is FALSE,
  #there's a problem to fix
  if(sum(duplicated(data_type$file_name)) != 0) break()
  if(sum(duplicated(data_type$unique)) != nrow(data_type) / 2) break()
  
  data_local <- map_df(data_files %>% `[`(grepl("local", .)), read_plus) %>%
    left_join(dplyr::select(data_type, -species), by = c("file_name"))
  data_gee <- map_df(data_files %>% `[`(grepl("gee", .)), read_plus) %>%
    left_join(dplyr::select(data_type, -species), by = c("file_name"))
  
  #check that things match! Otherwise we've got to fix something
  if(!identical(data_local$unique, data_gee$unique)) break()
  
  #combine data sets together
  combined <- cbind(data_local, dplyr::select(data_gee, aet:chili))
  
  combined <- combined %>%
    mutate(prop_spruce = replace_na(prop_spruce, 0),
           prop_oak = replace_na(prop_oak, 0))
  
  combos <- unique(combined$unique)
  model_list <- vector("list", length = length(combos))
  
  #fit all the models for each spatial scale
  for(i in 1:length(combos)){
    data_subset <- combined %>%
      filter(unique == combos[i])
    brt_all<- dismo::gbm.step(data = as.data.frame(data_subset[complete.cases(data_subset), ]), 
                              gbm.x = c("aet", "pdsi", "pet", "pr", "soil", "tmmx", "tmmn", "vpd", "def",
                                        #"Npp",
                                        "tpi", "chili", "slope",
                                        "height", "biomass", "fhd_normal", "understory_ratio", "open_area",
                                        "prop_forest", "prop_decid", "prop_conifer","prop_shrub", "prop_grass",
                                        "prop_water", "prop_dev_light", "prop_dev_heavy", "prop_open",
                                        #"prop_spruce", "prop_oak",
                                        "time_observations_started", "duration_minutes"), 
                              gbm.y = "species_observed",
                              tree.complexity = 5,
                              learning.rate = 0.01,
                              n.folds = 10)
    model_list[i] <- list(brt_all)
  }
  
  names(model_list) <- combos
  
  #wrangle models and add them to the overall list
  model_tib <- tibble(names = combos,
                      models = model_list)
  models_all <- bind_rows(models_all, model_tib)
  
  
  }


vars_to_plot <- c("height", "biomass", "fhd_normal", "understory_ratio", "open_area")

walk(model_list, print)
walk(model_list, summary)
walk(model_list, function(x) gbm.plot(x,n.plots = 9, plot.layout = c(3,3)))

model_cv_stats_list <- map(model_list, .f = function(x) pluck(x) %>% 
                             pluck("cv.statistics") %>% 
                             pluck("discrimination.mean"))

var_importance_list <- map(model_list, .f = function(x) pluck(x) %>% 
                             pluck("contributions"))

gedi_vars_preds_list <- map(model_list, .f = function(mod){map(vars_to_plot, function(vars) plot.gbm(mod, vars, return.grid = TRUE))})


test <- tibble(model = combos, 
          cv = model_cv_stats_list, 
          var_imp = var_importance_list, 
          preds = gedi_vars_preds_list)




test <- model_list[[1]]
fits<- gbm.plot(test, n.plots = 9, plot.layout = c(3,3))

