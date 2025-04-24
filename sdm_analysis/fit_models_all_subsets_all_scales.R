### fit models
library("tidyverse")
library("pROC")
library("gbm")
library("dismo")
library("terra")
library("sf")
library("tidyterra")
#--------------------------

species <- "cerw"
scale_list <- c("100", "500", "1000", "5000", "10000", "effort")


for(i in 1:length(scale_list)){
  
  scale <- scale_list[i]
  ###### BRT
  combined <- list.files("./environment_vars/combined", pattern = species, full.names = TRUE) %>%
    `[`(grep("all_scales", .)) %>%
    read_csv() %>%
    filter(scale == scale)
  
  
  # combined$species_observed <- ifelse(as.character(combined$species_observed) == "FOUND", TRUE, FALSE)
  brt_all<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                         gbm.x = c("aet", "pdsi", "pet", "pr", "soil", "tmmx", "tmmn", "vpd", "def",
                                           "tpi", "chili", "slope",
                                           "height", "biomass", "fhd", "understory", "area_short",
                                           "prop_forest", "prop_decid", "prop_conifer","prop_shrub", "prop_grass",
                                           "prop_water", "prop_dev_light", "prop_dev_heavy", "prop_open",
                                           "prop_hardwood", "prop_spruce", "prop_oak",
                                           "time_observations_started", "duration_minutes"), 
                         gbm.y = "species_observed",
                         tree.complexity = 5,
                         learning.rate = 0.01,
                         n.folds = 10)
  
  # print(brt_all)
  # summary.gbm(brt_all)
  # brt_int <- gbm.interactions(brt_all)
  # gbm.plot(brt_all, n.plots = 9, plot.layout = c(3,3))
  # gbm.perspec(brt_all, x = brt_int$rank.list$var1.index[1], y = brt_int$rank.list$var2.index[1], theta = 210)
  # 
  # brt_all$self.statistics$discrimination
  # brt_all$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_all$cv.statistics$correlation.mean
  # sqrt(mean(brt_all$residuals^2))
  # min(brt_all$train.error)
  # 
  # roc(combined$species_observed, boot::inv.logit(predict(brt_all)),
  #     plot = TRUE) #BRTs often overfit
  saveRDS(brt_all, paste0("./sdm_analysis/sdms/", species, "_dist_model_full_", scale, "_v2.RDS"))
  
  
  # brt_simp <- dismo::gbm.simplify(brt_all)
  # plot(brt_simp$deviance.summary$mean ~ c(1:length(brt_simp$deviance.summary$mean)))
  # preds_selected <- brt_simp$pred.list[which(brt_simp$deviance.summary$mean == min(brt_simp$deviance.summary$mean))]
  # brt2<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
  #                        gbm.x = preds_selected, #TODO get this automatically from the simp model
  #                        gbm.y = "species_observed")
  # 
  # saveRDS(brt2, paste0("./sdm_analysis/sdms/", species, "_dist_model_simplified.RDS"))
  # 
  # 
  # gbm.plot(brt2,
  #          n.plots = 9,
  #          common.scale = TRUE,
  #          show.contrib = TRUE,
  #          plot.layout = c(3,3))
  # gbm.perspec(brt, x = 7, y = 6)
  
  #----------------------
  # Just NLCD + Treemap
  #----------------------
  brt_veg <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x = c("prop_forest", "prop_decid", "prop_conifer","prop_shrub", "prop_grass",
                                     "prop_water","prop_open", "prop_dev_light", "prop_dev_heavy", 
                                     "prop_hardwood", #"prop_spruce", "prop_oak",
                                     "time_observations_started", "duration_minutes"),
                             gbm.y = "species_observed",
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             n.folds = 10)
  # summary.gbm(brt_veg)
  # gbm.plot(brt_veg)
  # brt_veg$self.statistics$discrimination
  # brt_veg$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_veg$cv.statistics$correlation.mean
  # sqrt(mean(brt_veg$residuals^2))
  # min(brt_veg$train.error)
  saveRDS(brt_veg, paste0("./sdm_analysis/sdms/", species, "_dist_model_veg_", scale, "_v2.RDS"))
  
  #----------------------
  # NLCD + Treemap + GEDI
  #----------------------
  
  brt_veg_gedi <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                             gbm.x = c("height", "biomass", "fhd", "understory", "area_short",
                                       "prop_forest", "prop_decid", "prop_conifer","prop_shrub", "prop_grass",
                                       "prop_water", "prop_open", #"prop_dev_light", "prop_dev_heavy", 
                                       "prop_hardwood", #"prop_spruce", "prop_oak",
                                       "time_observations_started", "duration_minutes"),
                             gbm.y = "species_observed",
                             tree.complexity = 5,
                             learning.rate = 0.01,
                             n.folds = 10)
  # summary.gbm(brt_veg_gedi)
  # gbm.plot(brt_veg_gedi)
  # brt_veg_gedi$self.statistics$discrimination
  # brt_veg_gedi$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_veg_gedi$cv.statistics$correlation.mean
  # sqrt(mean(brt_veg_gedi$residuals^2))
  # min(brt_veg_gedi$train.error)
  saveRDS(brt_veg_gedi, paste0("./sdm_analysis/sdms/", species, "_dist_model_veg_gedi_", scale, "_v2.RDS"))
  
  #-------------------------------------
  # Topography and climate
  #---------------------------------------
  
  brt_clim <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                              gbm.x = c("pdsi", "pet", "pr", "soil", "tmmn", 
                                        "tmmx", "vpd", "def", "tpi", "chili", "slope",
                                        "time_observations_started", "duration_minutes"),
                              gbm.y = "species_observed",
                              tree.complexity = 5,
                              learning.rate = 0.01,
                              n.folds = 10)
  # summary.gbm(brt_clim)
  # gbm.plot(brt_clim)
  # brt_clim$self.statistics$discrimination
  # brt_clim$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_clim$cv.statistics$correlation.mean
  # sqrt(mean(brt_clim$residuals^2))
  # min(brt_clim$train.error)
  saveRDS(brt_clim, paste0("./sdm_analysis/sdms/", species, "_dist_model_clim_", scale, "_v2.RDS"))
  
  #-------------------------------------
  # Topography and climate + GEDI
  #---------------------------------------
  
  brt_clim_gedi <- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                              gbm.x = c("height", "biomass", "fhd", "understory", "area_short",
                                        "pdsi", "pet", "pr", "soil", "tmmn", 
                                        "tmmx", "vpd", "def", "tpi", "chili", "slope",
                                        "time_observations_started", "duration_minutes"),
                              gbm.y = "species_observed",
                              tree.complexity = 5,
                              learning.rate = 0.01,
                              n.folds = 10)
  # summary.gbm(brt_clim_gedi)
  # gbm.plot(brt_clim_gedi)
  # brt_clim_gedi$self.statistics$discrimination
  # brt_clim_gedi$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_clim_gedi$cv.statistics$correlation.mean
  # sqrt(mean(brt_clim_gedi$residuals^2))
  # min(brt_clim_gedi$train.error)
  saveRDS(brt_clim_gedi, paste0("./sdm_analysis/sdms/", species, "_dist_model_clim_gedi_", scale, "_v2.RDS"))
  
  
  #-------------------------------------------------------
  # All without GEDI
  #--------------------------------------------------------
  brt_all_no_gedi<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                            gbm.x =  c("aet", "pdsi", "pet", "pr", "soil", "tmmx", "tmmn", "vpd", "def",
                                       "tpi", "chili", "slope",
                                       "prop_forest", "prop_decid", "prop_conifer","prop_shrub", "prop_grass",
                                       "prop_water", "prop_open",  #"prop_dev_light", "prop_dev_heavy",
                                       "prop_hardwood", #"prop_spruce", "prop_oak", 
                                       "time_observations_started", "duration_minutes"), 
                            gbm.y = "species_observed",
                            tree.complexity = 5,
                            learning.rate = 0.01,
                            n.folds = 10)
  # 
  # summary.gbm(brt_all_no_gedi)
  # gbm.plot(brt_all_no_gedi)
  # brt_all_no_gedi$self.statistics$discrimination
  # brt_all_no_gedi$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_all_no_gedi$cv.statistics$correlation.mean
  # sqrt(mean(brt_all_no_gedi$residuals^2))
  # min(brt_all_no_gedi$train.error)
  
  saveRDS(brt_all_no_gedi, paste0("./sdm_analysis/sdms/", species, "_dist_model_no_gedi_", scale, "_v2.RDS"))
  
  #-------------------------------------------------------
  # GEDI only
  #--------------------------------------------------------
  brt_gedi<- dismo::gbm.step(data = as.data.frame(combined[complete.cases(combined), ]), 
                                    gbm.x = c("height", "biomass", "fhd", "understory", "area_short",
                                              "time_observations_started", "duration_minutes"), 
                                    gbm.y = "species_observed",
                                    tree.complexity = 5,
                                    learning.rate = 0.01,
                                    n.folds = 10)
  
  # summary.gbm(brt_gedi)
  # gbm.plot(brt_gedi)
  # brt_gedi$self.statistics$discrimination
  # brt_gedi$cv.statistics$discrimination.mean #mean(brt_all$cv.roc.matrix)
  # brt_gedi$cv.statistics$correlation.mean
  # sqrt(mean(brt_gedi$residuals^2))
  # min(brt_gedi$train.error)
  
  saveRDS(brt_gedi, paste0("./sdm_analysis/sdms/", species, "_dist_model_gedi_only_", scale, "_v2.RDS"))

}