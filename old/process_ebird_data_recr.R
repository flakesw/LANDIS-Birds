#analyze CERW detection via GEE

#import and wrangle ebird data

library("auk")
library("tidyverse")
library("lubridate")
library("sf")
library("dismo")

#Data pre-treatment -------------------------------------------

output_file <- "ebd_filtered_recr.txt"
output_sampling <- "ebd_sampling_filtered.txt"
f_ebd <- "D:/Data/ebird_data/ebd_US_redcro_201001_202301_relJul-2023/ebd_US_redcro_201001_202301_relJul-2023.txt"
f_smp <- "D:/Data/ebird_data/ebd_sampling_relJul-2023.txt"
ebd_filters <-  auk_ebd(f_ebd, file_sampling = f_smp) %>%  
  auk_bcr(c(12,13,14,22,23,24,27,28,29,30)) %>% 
  auk_complete() %>%
  auk_filter(file = output_file,
             file_sampling = output_sampling, 
             overwrite = TRUE)
ebd_zf <- auk_zerofill(ebd_filters)
ebd_zf_df <- collapse_zerofill(ebd_zf)

# ebd_thin <- ebd_zf_df %>%
#   dplyr::group_by(species_observed) %>%
#   dplyr::slice_sample(n = 10000)
# ebd_thin <- ebd_zf_df
# rm(ebd_zf, ebd_zf_df)
rm(ebd_zf)

# ebd_zf_df <- ebd_thin

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


ebd_zf_filtered <- ebd_zf_df %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year >= 2010,
    # 10 or fewer observers
    number_observers <= 10)

ebird <- ebd_zf_filtered %>% 
  dplyr::select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebird, "recr_filtered.csv", na = "")


#-----------------------set up for GEE

ebird <- read.csv("recr_filtered.csv")
ebird_sf <- read.csv("recr_filtered.csv") %>% 
  # convert to spatial points
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  dplyr::select(species_observed)


#----- figure
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
# plot(st_geometry(ebird_sf), col = NA)
# contextual gis data
# plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
# plot(bcr, col = "#cccccc", border = NA, add = TRUE)
# plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
# plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(st_geometry(ebird_sf),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)
# observed
plot(filter(ebird_sf, species_observed) %>% st_geometry(),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklists", "CERW sightings"),
       pch = 19)
box()
par(new = TRUE, mar = c(0, 0, 3, 0))
title("Cerulean warbler eBird Observations\nJune 2010-2023, BCR 28")


ext <- sf::st_bbox(ebird_sf) %>%
  as.vector()%>%
  as.data.frame() %>%
  rename(left, bottom, right, top)

library(ggmap)
ggmap::register_google(key = "AIzaSyCIhXCi5LQjdbGF4MFzaOTg4Nq9CT1SIwQ")

nc_map <- get_map(location = ext, zoom = 10)

ggmap(nc_map) +
  geom_sf(data = ebird_sf, aes(color = species_observed),
          show.legend = "point", inherit.aes = FALSE) +
  theme_minimal()

install.packages("simplevis")
library("simplevis")
leaf_sf_col(ebird_sf, 
               col_var = species_observed)



#------------spatial thin
library(dggridR)
##spatially thin to hex grid
# generate hexagonal grid with 2.5 km between cells (half the maximum travel distance)
dggs <- dgconstruct(spacing = 2500)
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

write.csv(ebird_ss, "recr_subsampled.csv")

#figure -------------------------------------------------------

all_pts <- ebird %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
  st_transform(crs = 4326) %>% 
  dplyr::select(species_observed)
ss_pts <- ebird_ss %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
  st_transform(crs = 4326) %>% 
  dplyr::select(species_observed)
both_pts <- list(before_ss = all_pts, after_ss = ss_pts)



p <- par(mfrow = c(2, 1))
for (i in seq_along(both_pts)) {
  par(mar = c(0.25, 0.25, 0.25, 0.25))
  # set up plot area
  plot(st_geometry(both_pts[[i]]), col = NA)
  # contextual gis data
  # plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  # plot(bcr, col = "#cccccc", border = NA, add = TRUE)
  # plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  # plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  # ebird observations
  # not observed
  plot(st_geometry(both_pts[[i]]),
       pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
       add = TRUE)
  # observed
  plot(filter(both_pts[[i]], species_observed) %>% st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
       add = TRUE)
  # legend
  legend("bottomright", bty = "n",
         col = c("#555555", "#4daf4a"),
         legend = c("Non-detection", "Detection"),
         pch = 19)
  box()
  par(new = TRUE, mar = c(0, 0, 3, 0))
  if (names(both_pts)[i] == "before_ss") {
    title("CERW eBird Observations\nBefore subsampling")
  } else {
    title("After subsampling")
  }
}
par(p)



