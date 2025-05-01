library("tidyverse")


# species_list <- c("cerw", "gwwa")

#full species list
# species_list <- c("bcch", "blbw", "bwwa", "cerw",
#                   "evgr", "ewpw", "gwwa", "heth",
#                   "kewa", "lowa", "nawa", "nswo",
#                   "oven", "praw", "prow", "recr",
#                   "rugr", "swwa", "veer", "wewa",
#                   "woth", "ytwa")

species_list <- c("acfl") #, "alfl", "bhnu", "osfl")

#some helper functions
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(file_name = basename(flnm))
}

for(species in species_list){

  local_data_files <- list.files("./environment_vars/local/", pattern = species, full.names = TRUE) %>%
    `[`(grep("combined", .)) 
  data_type_local <- data.frame(file_name = basename(local_data_files),
                          file_location = local_data_files) %>%
    mutate(species = map(strsplit(file_name, "[[:punct:]]"), pluck(1)) %>% unlist(),
           data_source = map(strsplit(file_name, "[[:punct:]]"), pluck(5)) %>% unlist(),
           scale = map(strsplit(file_name, "[[:punct:]]"), pluck(6)) %>% unlist(),
           unique = paste(species, scale, sep = "_"))
  
  gee_data_files <- list.files("./environment_vars/gee_temp/", pattern = species, full.names = TRUE) %>%
    `[`(grep("combined", .)) 
  data_type_gee <- data.frame(file_name = basename(gee_data_files),
                                file_location = gee_data_files) %>%
    mutate(species = map(strsplit(file_name, "[[:punct:]]"), pluck(1)) %>% unlist(),
           data_source = map(strsplit(file_name, "[[:punct:]]"), pluck(5)) %>% unlist(),
           scale =   map(strsplit(file_name, "[[:punct:]]|group"), pluck(6)) %>%
                     unlist(),
           unique = paste(species, scale, sep = "_"))
  
  #are there any duplicate files or species/scale combos? If either of these is FALSE,
  #there's a problem to fix
  # if(sum(duplicated(data_type$file_name)) != 0) break()
  # if(sum(duplicated(data_type$unique)) != nrow(data_type)) break()
  
  data_local <- map_df(local_data_files, read_plus) %>%
    left_join(dplyr::select(data_type_local, -species), by = c("file_name")) %>%
    arrange(scale, checklist_id)
  data_gee <- map_df(gee_data_files, read_plus) %>%
    left_join(dplyr::select(data_type_gee, -species), by = c("file_name")) %>%
    arrange(scale, checklist_id)
  
  #check that things match! Otherwise we've got to fix something
  if(!identical(data_local$unique, data_gee$unique)) break()
  
  data_local <- mutate(data_local,
                       fhd = fhd/1000,
                       understory = understory/1000,
                       height = height/1000) %>%
    dplyr::arrange(species, scale, species_observed, checklist_id)
  data_gee <- mutate(data_gee,
                     tmmn = tmmn/10, #convert to C
                     tmmx = tmmx/10,
                     aet = aet/10, #convert to mean mm/month
                     pet = pet/10,
                     def = def/10,
                     vpd = vpd/100, #convert to mean kPa
                     pdsi = pdsi/100) %>%
    dplyr::arrange(species, scale, species_observed, checklist_id)
  
  #this needs to be true or the data gets all messed up
  identical(paste0(data_local$unique, data_local$checklist_id), paste0(data_gee$unique, data_gee$checklist_id))
  
  #combine data sets together
  combined <- cbind(data_local, dplyr::select(data_gee, aet:chili))
  
  combined <- combined %>%
    mutate(prop_spruce = replace_na(prop_spruce, 0),
           prop_oak = replace_na(prop_oak, 0))
  
  write.csv(combined, paste0("./environment_vars/combined/", species, "_all_scales.csv"))
  
}
