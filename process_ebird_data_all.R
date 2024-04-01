#import and wrangle ebird data

library("auk")
library("tidyverse")
library("lubridate")
library("sf")

#Make sure to change species *everywhere*!
species <- "gwwa"
species_long <- "gowwar"
species_common <- "Golden-winged warbler"
season_start <- "*-06-15"
season_end <- "*-07-30"
scale <- NA #TODO add different spatial scales for aggregation

#Data pre-treatment -------------------------------------------

output_file <- paste0("./ebird/ebd_filtered_", species, ".txt")
output_sampling <- paste0("./ebird/ebd_sampling_filtered", species, ".txt")
f_ebd <- paste0("E:/ebird_data/ebd_US_relJan-2024.txt") 
f_smp <- "E:/ebird_data/ebd_US_relJan-2024_sampling.txt"
ebd_filters <-  auk_ebd(f_ebd, file_sampling = f_smp) %>% 
  auk_species(species_common, taxonomy_version = 2023) %>%
  # auk_bcr(c(12,13,14,22,23,24,27,28,29,30)) %>%
  auk_bcr(28) %>% #Just for Appalachian BCR
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


time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

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
         number_observers)
# write_csv(ebird, paste0("./ebird/", species, "_filtered.csv"), na = "")

#-----------------------set up for GEE

# ebird <- read.csv(paste0("./ebird/", species, "_filtered.csv"))

ebird_sf <- ebird %>%
  # convert to spatial points
  dplyr::rename_with(tolower) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = 5070) %>%
  dplyr::select(species_observed, checklist_id)


#------------spatial thin
library(dggridR)
##spatially thin to hex grid
# generate hexagonal grid with 1.5 km between cells (half the maximum travel distance)
dggs <- dgconstruct(spacing = 1.55)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird %>%
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))
# sample one checklist per grid cell per week
# sample detection/non-detection independently
ebird_ss <- checklist_cell %>%
  group_by(species_observed, year, week, cell) %>%
  sample_n(size = 1) %>%
  ungroup()

# write.csv(ebird_ss, paste0("./ebird/", species, "_subsampled.csv"))

#------balance majority

ebird_ss_bal <- ebird_ss %>%
  group_by(species_observed) %>%
  sample_n(size = sum(.$species_observed == TRUE))
ebird_ss_sf <- ebird_sf %>%
  filter(checklist_id %in% ebird_ss_bal$checklist_id)
plot(ebird_ss_sf["species_observed"])

write.csv(ebird_ss_bal, paste0("./ebird/", species, "_subsampled_balanced_bcr28.csv"))

# #figure -------------------------------------------------------
# 
# all_pts <- ebird %>%  
#   st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
#   st_transform(crs = 4326) %>% 
#   dplyr::select(species_observed)
# ss_pts <- ebird_ss %>%  
#   st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
#   st_transform(crs = 4326) %>% 
#   dplyr::select(species_observed)
# both_pts <- list(before_ss = all_pts, after_ss = ss_pts)
# 
# 
# 
# p <- par(mfrow = c(2, 1))
# for (i in seq_along(both_pts)) {
#   par(mar = c(0.25, 0.25, 0.25, 0.25))
#   # set up plot area
#   plot(st_geometry(both_pts[[i]]), col = NA)
#   # contextual gis data
#   # plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
#   # plot(bcr, col = "#cccccc", border = NA, add = TRUE)
#   # plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
#   # plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#   # ebird observations
#   # not observed
#   plot(st_geometry(both_pts[[i]]),
#        pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
#        add = TRUE)
#   # observed
#   plot(filter(both_pts[[i]], species_observed) %>% st_geometry(),
#        pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
#        add = TRUE)
#   # legend
#   legend("bottomright", bty = "n",
#          col = c("#555555", "#4daf4a"),
#          legend = c("Non-detection", "Detection"),
#          pch = 19)
#   box()
#   par(new = TRUE, mar = c(0, 0, 3, 0))
#   if (names(both_pts)[i] == "before_ss") {
#     title("cerw eBird Observations\nBefore subsampling")
#   } else {
#     title("After subsampling")
#   }
# }
# par(p)



