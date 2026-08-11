# ==============================================================================
# Step 1.1 - Land cover composition by microbasin
#
# Computes, for each microbasin and each year (1985-2024), the area and
# percentage covered by each functional land cover group, based on MapBiomas
# annual rasters.
#
# CRS handling:
#   - MapBiomas rasters are in EPSG:4326 (geographic)
#   - Microbasin are in EPSG:32718 (UTM 18S, metric)
#   Areas are computed in EPSG:32718; zonal extraction is done in EPSG:4326
#   to avoid resampling the rasters.
#
# Output (workspace/01_lulc_by_microbasin.csv), one row per microbasin x
# functional group x year:
#   comid          - Microbasin identifier (linkno from geoglows_basins.gpkg).
#   group          - Functional land cover group name (see reclass_table
#                    below for the MapBiomas class codes making up each
#                    group). Every microbasin/year has one row per group,
#                    even if its coverage is 0 (fraction filled with 0 by
#                    complete()).
#   year           - Year of the MapBiomas annual raster (1985-2024).
#   area_km2       - Area covered by this group within the microbasin, in
#                    km2 (= fraction of covered cells x basin_area_km2).
#   basin_area_km2 - Total area of the microbasin, in km2 (constant across
#                    groups/years for a given comid).
#   pct            - Percentage of the microbasin covered by this group
#                    (= fraction x 100). Percentages across all groups for
#                    a given comid/year sum to ~100% minus any NA cells
#                    (e.g. no-data/cloud) excluded by exact_extract's "frac".
#
# ==============================================================================

library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

raster_dir  <- "inputs/map_biomas"
basins_path <- "inputs/geographical_data/geoglows_basins.gpkg"
output_path <- "workspace/01_lulc_by_microbasin.csv"

# ------------------------------------------------------------------------------
# 2. Reclassification table: MapBiomas class code -> functional group
# ------------------------------------------------------------------------------

reclass_table <- tribble(
  ~code, ~group_id, ~group,
  # A. Natural forest
  3,  1, "Bosque natural",
  4,  1, "Bosque natural",
  6,  1, "Bosque natural",

  # B. High-Andean vegetation (paramo / bofedal)
  81, 2, "Vegetacion altoandina",
  82, 2, "Vegetacion altoandina",
  11, 2, "Vegetacion altoandina",

  # C. Other non-forest natural formation
  12, 3, "Otra formacion natural",
  13, 3, "Otra formacion natural",
  29, 3, "Otra formacion natural",

  # D. Agriculture, silviculture and livestock
  9,  4, "Agricultura, silvicultura y ganaderia",
  21, 4, "Agricultura, silvicultura y ganaderia",
  74, 4, "Agricultura, silvicultura y ganaderia",

  # E. Anthropic infrastructure
  24, 5, "Infraestructura antropica",
  25, 5, "Infraestructura antropica",

  # F. Mining
  30, 6, "Mineria",

  # G. Natural non-vegetated area
  23, 7, "Area no vegetada natural",
  68, 7, "Area no vegetada natural",

  # H. Natural water bodies
  33, 8, "Cuerpos de agua",
  34, 8, "Cuerpos de agua",

  # I. Aquaculture
  31, 9, "Acuicultura"
)

# Lookup of group_id -> group name, used to label extraction results
group_labels <- reclass_table %>%
  distinct(group_id, group)

# Two-column matrix (from -> to) required by terra::classify
reclass_matrix <- as.matrix(reclass_table[, c("code", "group_id")])

# ------------------------------------------------------------------------------
# 3. Load micro-watersheds and compute their area in the metric CRS
# ------------------------------------------------------------------------------

basins <- st_read(basins_path, quiet = TRUE) %>%
  select(comid = linkno) %>%
  mutate(basin_area_km2 = as.numeric(st_area(.)) / 1e6)

# Geographic version used for zonal extraction against the rasters
basins_geo <- st_transform(basins, 4326)

# ------------------------------------------------------------------------------
# 4. Process each annual raster
# ------------------------------------------------------------------------------

raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)

extract_year <- function(file_path) {

  year <- as.integer(tools::file_path_sans_ext(basename(file_path)))
  message("Processing ", year)

  # Reclassify MapBiomas codes into functional groups; unlisted codes -> NA
  r <- classify(rast(file_path), reclass_matrix, others = NA)

  # "frac" returns the coverage-weighted fraction of each group per polygon,
  # as columns named frac_<group_id>. Fractions exclude NA cells.
  fractions <- exact_extract(r, basins_geo, "frac", progress = FALSE)

  bind_cols(st_drop_geometry(basins_geo), fractions) %>%
    pivot_longer(
      cols         = starts_with("frac_"),
      names_to     = "group_id",
      names_prefix = "frac_",
      values_to    = "fraction"
    ) %>%
    mutate(
      group_id = as.integer(group_id),
      year     = year
    )
}

composition <- lapply(raster_files, extract_year) %>%
  bind_rows()

# ------------------------------------------------------------------------------
# 5. Complete missing combinations, derive areas and format output
# ------------------------------------------------------------------------------

composition <- composition %>%
  select(comid, group_id, year, fraction) %>%
  complete(comid, group_id, year, fill = list(fraction = 0)) %>%
  left_join(group_labels, by = "group_id") %>%
  left_join(st_drop_geometry(basins), by = "comid") %>%
  mutate(
    area_km2 = fraction * basin_area_km2,
    pct      = fraction * 100
  ) %>%
  select(comid, group, year, area_km2, basin_area_km2, pct) %>%
  arrange(comid, year, group)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(composition, output_path, row.names = FALSE)

message("Done: ", nrow(composition), " rows written to ", output_path)
