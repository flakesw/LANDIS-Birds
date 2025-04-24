model_list <- readRDS(TODO)

vars_to_plot <- c("height", "biomass", "fhd_normal", "understory_ratio", "open_area")

walk(model_list, print)
walk(model_list, summary)
walk(model_list, function(x) gbm.plot(x,n.plots = 9, plot.layout = c(3,3)))

model_cv_stats_list <- map(model_list, .f = function(x) pluck(x) %>% 
                             pluck("cv.statistics") %>% 
                             pluck("discrimination.mean"))

var_importance_list <- map(model_list, .f = function(x) pluck(x) %>% 
                             pluck("contributions"))

gedi_vars_preds_list <- map(model_list, .f = function(mod){map(vars_to_plot, function(vars) 
  plot.gbm(mod, vars, return.grid = TRUE))})


test <- tibble(model = combos, 
               cv = model_cv_stats_list, 
               var_imp = var_importance_list, 
               preds = gedi_vars_preds_list)

