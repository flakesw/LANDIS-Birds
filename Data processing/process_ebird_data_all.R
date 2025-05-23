#import and wrangle ebird data

library("auk")
library("tidyverse")
library("lubridate")
library("sf")
library("dggridR")

# ebird_taxonomy <- read.csv("E:/ebird_data/Clements-v2023-October-2023.csv")
banding_codes <- read.csv("./bird_info/banding_codes.csv")
species_data <- data.frame(species = character(4))

# #which species to use?
# species_data$species <- c("bcch", "ewpw", "evgr", "recr", "wewa",
#                           "swwa", "nawa", "lowa", "prow", "cerw",
#                           "kewa", "praw", "ytwa", "blbw", "gwwa",
#                           "rugr", "nswo", "veer", "heth", "oven",
#                           "bwwa", "woth")
# species_data$species_common <- banding_codes[match(species_data$species, tolower(banding_codes$Alpha.Code)), "Common.Name"]
# species_data$species_scientific <- banding_codes[match(species_data$species, tolower(banding_codes$Alpha.Code)), "Scientific.Name"]
# species_data$species_common[16] <- "Ruffed Grouse"
# species_data$species_scientific[16] <- "Bonasa umbellus"

species_data$species <- c("acfl", "alfl", "bhnu", "osfl")
species_data$species_common <- banding_codes[match(species_data$species, tolower(banding_codes$Alpha.Code)), "Common.Name"]
species_data$species_scientific <- banding_codes[match(species_data$species, tolower(banding_codes$Alpha.Code)), "Scientific.Name"]

season_start <- "*-06-15"
season_end <- "*-07-30"
region <- "BCR28"
# region <- "all_regions"
max_n <- 5000 #maximum number of presences to extract


time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}


species <- species_data$species
species_common <- species_data$species_common

output_file <- paste0("./ebird/ebd_filtered_all_species.txt")
output_sampling <- paste0("./ebird/ebd_sampling_filtered_all_species.txt")
f_ebd <- paste0("E:/ebird_data/ebd_US_relJan-2024.txt") 
f_smp <- "E:/ebird_data/ebd_US_relJan-2024_sampling.txt"
ebd_filters <-  auk_ebd(f_ebd, file_sampling = f_smp) %>% 
  auk_species(species_common, taxonomy_version = 2023) %>%
  auk_bcr(if(region == "BCR28"){28} else{ c(12,13,14,22,23,24,27,28,29,30)}) %>% #for just BCR28 or surrounding BCRs too?
  auk_complete() %>%
  auk_year(year = 2013:2022) %>%
  auk_date(date = c(season_start, season_end)) %>% 
  auk_distance(distance = c(0,3.1)) %>%
  auk_duration(c(0,180)) %>%
  auk_filter(file = output_file,
             file_sampling = output_sampling, #TODO figure out if we can do the output separately
             overwrite = TRUE)
ebd_zf <- auk_zerofill(ebd_filters)
ebd_zf_df <- collapse_zerofill(ebd_zf)

# clean up variables
ebd_zf_df <- ebd_zf_df %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

ebird <- ebd_zf_df %>% 
  dplyr::select(checklist_id, observer_id, sampling_event_identifier,
                scientific_name,
                observation_count, species_observed, 
                state_code, locality_id, latitude, longitude,
                protocol_type, all_species_reported,
                observation_date, year, day_of_year,
                time_observations_started, 
                duration_minutes, effort_distance_km,
                number_observers) %>%
  left_join(species_data, by = c("scientific_name" = "species_scientific"))

write_csv(ebird, paste0("./ebird/all_species_filtered.csv"), na = "")

#-spatial and temporal thin----------------------------------------------------
##spatially thin to hex grid

# generate hexagonal grid with 3.1 km between cells (maximum travel distance)
# or spatial scale if it's larger than 3.1
dggs <- dgconstruct(spacing = ifelse(scale > 3.1, scale, 3.1))
dggs <-  dgconstruct(spacing = 3.1)
# dggs_10 <- dgconstruct(spacing = 10)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird %>%
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         month = month(observation_date),
         week = week(observation_date))


for(i in 1:nrow(species_data)){

  # sample one checklist per grid cell per month
  # sample detection/non-detection independently
  ebird_ss <- checklist_cell %>%
    dplyr::filter(species == species_data$species[i]) %>%
    group_by(species_observed, year, month, cell) %>%
    sample_n(size = 1) %>%
    ungroup()
  
  write.csv(ebird_ss, paste0("./ebird/", species_data$species[i], "_subsampled.csv"))
  # ebird_ss <- read.csv(paste0("./ebird/", species, "_subsampled.csv"))
}

for(i in 1:nrow(species_data)){
  
  #------balance majority
  ebird_ss <- read.csv(paste0("./ebird/", species_data$species[i], "_subsampled.csv"))
  ebird_ss_bal <- ebird_ss %>%
    group_by(species_observed) %>%
    sample_n(size = min(max_n, sum(.$species_observed == TRUE)))
  
  write.csv(ebird_ss_bal, paste0("./ebird/", species_data$species[i], "_subsampled_balanced_", region, ".csv"))

}



# old function

# process_ebird_data <- function(species, species_common, season_start, season_end, region, scale, max_n){
#   #Data pre-treatment -------------------------------------------
#   
#   species <- species_data$species
#   species_common <- species_data$species_common
#   
#   output_file <- paste0("./ebird/ebd_filtered_", species, ".txt")
#   output_sampling <- paste0("./ebird/ebd_sampling_filtered_", species, ".txt")
#   f_ebd <- paste0("E:/ebird_data/ebd_US_relJan-2024.txt") 
#   f_smp <- "E:/ebird_data/ebd_US_relJan-2024_sampling.txt"
#   ebd_filters <-  auk_ebd(f_ebd, file_sampling = f_smp) %>% 
#     auk_species(species_common, taxonomy_version = 2023) %>%
#     auk_bcr(ifelse(region == "BCR28", 28, c(12,13,14,22,23,24,27,28,29,30))) %>% #for just BCR28 or surrounding BCRs too?
#     auk_complete() %>%
#     auk_year(year = 2013:2022) %>%
#     auk_date(date = c(season_start, season_end)) %>% 
#     auk_distance(distance = c(0,3.1)) %>%
#     auk_duration(c(0,180)) %>%
#     auk_filter(file = output_file,
#                file_sampling = output_sampling, #TODO figure out if we can do the output separately
#                overwrite = TRUE)
#   ebd_zf <- auk_zerofill(ebd_filters)
#   ebd_zf_df <- collapse_zerofill(ebd_zf)
#   
#   # clean up variables
#   ebd_zf_df <- ebd_zf_df %>% 
#     mutate(
#       # convert X to NA
#       observation_count = if_else(observation_count == "X", 
#                                   NA_character_, observation_count),
#       observation_count = as.integer(observation_count),
#       # effort_distance_km to 0 for non-travelling counts
#       effort_distance_km = if_else(protocol_type != "Traveling", 
#                                    0, effort_distance_km),
#       # convert time to decimal hours since midnight
#       time_observations_started = time_to_decimal(time_observations_started),
#       # split date into year and day of year
#       year = year(observation_date),
#       day_of_year = yday(observation_date)
#     )
#   
#   ebird <- ebd_zf_df %>% 
#     dplyr::select(checklist_id, observer_id, sampling_event_identifier,
#                   scientific_name,
#                   observation_count, species_observed, 
#                   state_code, locality_id, latitude, longitude,
#                   protocol_type, all_species_reported,
#                   observation_date, year, day_of_year,
#                   time_observations_started, 
#                   duration_minutes, effort_distance_km,
#                   number_observers)
#   write_csv(ebird, paste0("./ebird/", species, "_filtered.csv"), na = "")
# }


# 
# subset_and_balance_data <- function(species, species_common, season_start, season_end, region, scale, max_n)  {
#   #-----------------------set up for GEE
#   
#   ebird <- read.csv(paste0("./ebird/", species, "_filtered.csv"))
#   
#   # ebird_sf <- ebird %>%
#   #   # convert to spatial points
#   #   dplyr::rename_with(tolower) %>%
#   #   sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
#   #   sf::st_transform(crs = 5070) %>%
#   #   dplyr::select(species_observed, checklist_id)
#   
#   
#   #------------spatial thin
#   library(dggridR)
#   ##spatially thin to hex grid
#   # generate hexagonal grid with 1.5 km between cells (half the maximum travel distance)
#   
#   #3.1 km to be conservative
#   dggs <- dgconstruct(spacing = 3.1)
#   # dggs_10 <- dgconstruct(spacing = 10)
#   # get hexagonal cell id and week number for each checklist
#   checklist_cell <- ebird %>%
#     mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
#            year = year(observation_date),
#            month = month(observation_date),
#            week = week(observation_date))
#   # checklist_cell_absence <- filter(ebird, !species_observed) %>%
#   #   mutate(cell = dgGEO_to_SEQNUM(dggs_absences, longitude, latitude)$seqnum,
#   #          year = year(observation_date),
#   #          month = month(observation_date),
#   #          week = week(observation_date))
#   
#   # sample one checklist per grid cell per month
#   # sample detection/non-detection independently
#   ebird_ss <- checklist_cell %>%
#     group_by(species_observed, year, month, cell) %>%
#     sample_n(size = 1) %>%
#     ungroup()
#   
#   write.csv(ebird_ss, paste0("./ebird/", species, "_subsampled.csv"))
#   # ebird_ss <- read.csv(paste0("./ebird/", species, "_subsampled.csv"))
#   
#   #------balance majority
#   
#   ebird_ss_bal <- ebird_ss %>%
#     group_by(species_observed) %>%
#     sample_n(size = min(max_n, sum(.$species_observed == TRUE)))
#   # ebird_ss_sf <- ebird_sf %>%
#   #   filter(checklist_id %in% ebird_ss_bal$checklist_id)
#   # plot(ebird_ss_sf["species_observed"])
#   
#   write.csv(ebird_ss_bal, paste0("./ebird/", species, "_subsampled_balanced_", region,"_", Sys.Date(), ".csv"))
#   
#   ebird_ss_true <- ebird_ss %>%
#     filter(species_observed == TRUE)
#   ebird_ss_bal2 <- ebird_ss %>%
#     filter(species_observed == FALSE) %>%
#     sample_n(ifelse(nrow(.) > nrow(ebird_ss_true) * 10,
#                     nrow(ebird_ss_true) * 10,
#                     nrow(.))) %>%
#     bind_rows(ebird_ss_true)
#   
#   # ebird_ss_sf <- ebird_sf %>%
#   #   filter(checklist_id %in% ebird_ss_bal$checklist_id)
#   # plot(ebird_ss_sf["species_observed"])
#   
#   write.csv(ebird_ss_bal2, paste0("./ebird/", species, "_subsampled_balanced_extra_", region, ".csv"))
#   
#   # #figure -------------------------------------------------------
#   # 
#   # all_pts <- ebird %>%  
#   #   st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
#   #   st_transform(crs = 4326) %>% 
#   #   dplyr::select(species_observed)
#   # ss_pts <- ebird_ss %>%  
#   #   st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
#   #   st_transform(crs = 4326) %>% 
#   #   dplyr::select(species_observed)
#   # both_pts <- list(before_ss = all_pts, after_ss = ss_pts)
#   # 
#   # 
#   # 
#   # p <- par(mfrow = c(2, 1))
#   # for (i in seq_along(both_pts)) {
#   #   par(mar = c(0.25, 0.25, 0.25, 0.25))
#   #   # set up plot area
#   #   plot(st_geometry(both_pts[[i]]), col = NA)
#   #   # contextual gis data
#   #   # plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
#   #   # plot(bcr, col = "#cccccc", border = NA, add = TRUE)
#   #   # plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
#   #   # plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#   #   # ebird observations
#   #   # not observed
#   #   plot(st_geometry(both_pts[[i]]),
#   #        pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
#   #        add = TRUE)
#   #   # observed
#   #   plot(filter(both_pts[[i]], species_observed) %>% st_geometry(),
#   #        pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
#   #        add = TRUE)
#   #   # legend
#   #   legend("bottomright", bty = "n",
#   #          col = c("#555555", "#4daf4a"),
#   #          legend = c("Non-detection", "Detection"),
#   #          pch = 19)
#   #   box()
#   #   par(new = TRUE, mar = c(0, 0, 3, 0))
#   #   if (names(both_pts)[i] == "before_ss") {
#   #     title("cerw eBird Observations\nBefore subsampling")
#   #   } else {
#   #     title("After subsampling")
#   #   }
#   # }
#   # par(p)
#   
#   
# }
# 
# #need to restart at 5
# for(i in 6:22){
#   process_ebird_data(species = species_data$species[i], species_common = species_data$species_common[i], 
#                      season_start = season_start, season_end = season_end, region = region, 
#                      scale = scale, max_n = max_n)
#   subset_and_balance_data(species = species_data$species[i], species_common = species_data$species_common[i], 
#                           season_start = season_start, season_end = season_end, region = region, 
#                           scale = scale, max_n = max_n)
# }
