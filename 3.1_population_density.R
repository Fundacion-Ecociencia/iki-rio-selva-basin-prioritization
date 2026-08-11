# ==============================================================================
# Step 3.1 - Population density by micro-watershed
#
# Aggregates GHS-POP rasters (population count per pixel) to micro-watershed
# level for each available year and derives population density.
#
# Population is summed with coverage weighting, so pixels partially inside a
# micro-watershed contribute proportionally to their overlapping area.
# Areas are taken from the micro-watershed geometries in EPSG:32718 (metric).
#
# Output: population_by_microbasin.csv with columns
#   comid | year | population | basin_area_km2 | pop_density
# ==============================================================================

library(terra)
library(sf)
library(exactextractr)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

pop_dir     <- "inputs/ghs_pop"
basins_path <- "inputs/geographical_data/geoglows_basins.gpkg"
output_path <- "workspace/08_population_by_basin.csv"

# Regular expression capturing the year in the raster file names.
year_pattern <- "(\\d{4})"

# ------------------------------------------------------------------------------
# 2. Load micro-watersheds and compute their area
# ------------------------------------------------------------------------------

basins <- st_read(basins_path, quiet = TRUE) %>%
  select(comid = linkno) %>%
  mutate(basin_area_km2 = as.numeric(st_area(.)) / 1e6)

# ------------------------------------------------------------------------------
# 3. Extract total population per micro-watershed for one year
# ------------------------------------------------------------------------------

population_for_year <- function(file_path) {

  year <- as.integer(sub(paste0(".*", year_pattern, ".*"), "\\1",
                         basename(file_path)))
  message("Processing ", year)

  r <- rast(file_path)

  # Align the polygons to the raster CRS instead of resampling population counts
  basins_proj <- st_transform(basins, crs(r))

  tibble(
    comid      = basins$comid,
    year       = year,
    population = exact_extract(r, basins_proj, "sum", progress = FALSE)
  )
}

# ------------------------------------------------------------------------------
# 4. Process all years and derive density
# ------------------------------------------------------------------------------

raster_files <- list.files(pop_dir, pattern = "\\.tif$", full.names = TRUE)

population <- lapply(raster_files, population_for_year) %>%
  bind_rows() %>%
  left_join(st_drop_geometry(basins), by = "comid") %>%
  mutate(pop_density = population / basin_area_km2) %>%   # inhabitants / km2
  select(comid, year, population, basin_area_km2, pop_density) %>%
  arrange(comid, year)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(population, output_path, row.names = FALSE)

message("Done: ", nrow(population), " rows written to ", output_path)
