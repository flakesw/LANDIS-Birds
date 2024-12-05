#data nicely organized here: https://www1.usgs.gov/csas/swap/

info_dir <- "./bird_info"
sgcn_files <- list.files(info_dir, full.names = TRUE)
sgcn_data <- read_tsv(sgcn_files)

sgcn_data <- sgcn_data %>%
  filter(provided_taxonomy_group %in% c("Bird", "Birds")) %>%
  filter(!duplicated(`Common Name`))
