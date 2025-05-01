### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
library("terra")
library("sf")
library("tidyterra")
#--------------------------

#for this, you're going to want the combined data from "combine_local_and_gee_data.R"


# species_list <- c("cerw", "gwwa")

species_list <- c("bcch", "blbw",  "cerw", "gwwa", "heth", "kewa", "lowa",
                   "praw", "prow", "recr", "swwa", "veer", "wewa","ytwa",
                  "ewpw", "rugr", "nswo",  "acfl", "alfl", "bhnu")

scale_select <- 500


for(species in species_list){
  models_all <- tibble(names = character(length = 0L),
                       models = vector(mode = "list", length = 0L))
  
  combined <- list.files("./environment_vars/combined/", pattern = species, full.names = TRUE) %>%
    `[`(grep("combined", .)) %>%
    read_csv() %>%
    filter(scale == scale_select)
  
  combos <- unique(combined$unique)
  model_list <- vector("list", length = length(combos))
  
  #fit all the models for each spatial scale
  for(i in 1:length(combos)){
    data_subset <- combined %>%
      filter(unique == combos[i])
    brt_all<- dismo::gbm.step(data = as.data.frame(data_subset[complete.cases(data_subset), ]), 
                              gbm.x = c("aet", "pet", "pr", "soil", "tmmx", "tmmn", # "vpd", "def", "pdsi",
                                        "tpi", "chili", "slope", #"Npp",
                                        "height", "biomass", "fhd", "area_short", "understory", 
                                        "prop_forest", #"prop_decid", "prop_conifer", "prop_shrub", "prop_grass",
                                        #"prop_water", "prop_dev_heavy", "prop_open", "prop_dev_light"
                                        "prop_spruce", "prop_oak", "prop_hardwood",
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
  
  saveRDS(models_all, file = paste0("./sdm_analysis/sdms/", species, "full_model_landis_v2.RDS"))
  
}
