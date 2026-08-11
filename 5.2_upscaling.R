# ==============================================================================
# Step 5.2 - Upscaling of the pressure index to N6 hydrographic units
#
# Re-aggregates the microbasin-level results of step 4.2 (1747 GEOGloWS
# microbasins) to the 62 level-6 hydrographic units of the national Pfafstetter
# classification (cuencas_N6_50k_2020.gpkg), using an AREA-WEIGHTED MEAN.
#
# Method
# ------------------------------------------------------------------------------
# The two layers are not nested: a microbasin can straddle the boundary between
# two N6 units. The layers are therefore intersected, and each resulting piece
# contributes to its N6 unit in proportion to the area of that piece:
#
#   value_N6 = sum(area_piece * value_microbasin) / sum(area_piece)
#
# so a microbasin split 70/30 between two N6 units contributes 70% of its
# weight to one and 30% to the other, and large microbasins weigh more than
# small ones. Both layers cover the same 32 000 km2 study area, so in practice
# every N6 unit is fully covered (see pct_area_analizada below).
#
# NA handling: NA COUNTS AS 0. The tau_* columns arrive with the NA of steps
# 1.2-3.2 (constant or too-short series), which means "no trend measurable" and
# therefore no pressure. Those pieces enter the numerator as 0 while their area
# still counts in the denominator, so every column is averaged over the full
# unit area — the same NA -> 0 substitution step 4.2 already applied to the
# ori_*, tema_* and indice_presion columns, now extended to tau_* so all
# columns are weighted identically.
#
# Read the tau_* columns with that in mind: a value near 0 can mean either
# "weak trend across the whole unit" or "strong trend in a small part of it and
# no measurable trend elsewhere". n_subind_validos tells the two apart.
#
# CRS: both layers are reprojected to EPSG:32718, the metric CRS the rest of
# the pipeline computes areas in. The output is written in that CRS, so it
# overlays 10_/11_pressure_index*.gpkg directly.
#
# Output (workspace/12_pressure_index_v2_by_n6.gpkg), one feature per N6 unit
# (polygon geometry from cuencas_N6_50k_2020.gpkg):
#   N6_name            - Name of the level-6 hydrographic unit (kept from the
#                         source layer).
#   N5_name            - Name of the level-5 unit it belongs to (kept from the
#                         source layer).
#   area_km2           - Area of the N6 unit in EPSG:32718, in km2.
#   area_analizada_km2 - Area of the unit actually covered by microbasins from
#                         the index layer, in km2. [diagnostic]
#   pct_area_analizada - area_analizada_km2 / area_km2 * 100. Below 100% the
#                         unit extends beyond the study area and its averages
#                         describe only the covered part. [diagnostic]
#   n_microcuencas     - Number of microbasins contributing to the unit.
#                         [diagnostic]
#   tau_<indicator>    - Area-weighted mean of the raw mk_tau of each of the 15
#                         sub-indicators over the full unit area, counting NA
#                         as 0. A unit whose microbasins are all NA for an
#                         indicator therefore gets 0, not NA.
#   ori_<indicator>    - Area-weighted mean of the same value re-oriented
#                         towards pressure (+1 = more pressure).
#   tema_*             - Area-weighted mean of the five theme scores.
#   indice_presion     - Area-weighted mean of the final index. Range [-1, 1];
#                         +1 = highest pressure.
#   n_subind_validos   - Area-weighted mean of the number of sub-indicators
#                         with real data (out of 15), so a non-integer value
#                         is expected.
#
# The diagnostic columns are derived here, not taken from the source layer;
# set include_diagnostics <- FALSE to keep only N6_name, N5_name and the
# weighted averages.
# ==============================================================================

library(sf)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

# Point index_path at workspace/10_pressure_index_by_basin.gpkg to upscale the
# v1 index instead; the column set is discovered from the file, not hard-coded.
index_path  <- "workspace/11_pressure_index_v2_by_basin.gpkg"
n6_path     <- "inputs/geographical_data/cuencas_N6_50k_2020.gpkg"
output_path <- "workspace/12_pressure_index_v2_by_n6.gpkg"

# Metric CRS used for every area computation in this pipeline
target_crs <- "EPSG:32718"

# Fields kept from the N6 layer; everything else is dropped
n6_fields <- c("N6_name", "N5_name")

include_diagnostics <- TRUE

# ------------------------------------------------------------------------------
# 2. Load both layers in the common metric CRS
# ------------------------------------------------------------------------------

index <- st_read(index_path, quiet = TRUE) %>%
  st_transform(target_crs)

n6 <- st_read(n6_path, quiet = TRUE) %>%
  st_transform(target_crs) %>%
  select(all_of(n6_fields)) %>%
  mutate(n6_id = row_number(), area_km2 = as.numeric(st_area(.)) / 1e6)

# Every numeric column of the index except the identifier is upscaled. Keeping
# this data-driven means a sub-indicator added to step 4.2 flows through here
# without editing this script.
value_cols <- index %>%
  st_drop_geometry() %>%
  select(where(is.numeric), -comid) %>%
  names()

message("Upscaling ", length(value_cols), " columns from ", nrow(index),
        " microbasins to ", nrow(n6), " N6 units")

# ------------------------------------------------------------------------------
# 3. Intersect and weigh each piece by its area
# ------------------------------------------------------------------------------
# Attributes are constant over each polygon, so declaring the agr silences the
# "assumed to be spatially constant" warning of st_intersection().

st_agr(n6)    <- "constant"
st_agr(index) <- "constant"

pieces <- st_intersection(n6, index) %>%
  mutate(piece_area = as.numeric(st_area(.))) %>%
  st_drop_geometry() %>%
  # Shared borders yield zero-area line/point pieces; they carry no weight
  filter(piece_area > 0)

# ------------------------------------------------------------------------------
# 4. Area-weighted mean, ignoring pieces where the value is missing
# ------------------------------------------------------------------------------

weighted_mean_area <- function(x, w) {

  keep <- is.finite(w) & w > 0

  if (!any(keep)) return(NA_real_)

  x <- x[keep]
  w <- w[keep]

  # NA (and the odd non-finite tau) means "no trend measurable", which is 0
  # pressure, so it enters the numerator as 0 while its area still counts in
  # the denominator. Same substitution step 4.2 applies to the ori_* columns.
  x[!is.finite(x)] <- 0

  sum(x * w) / sum(w)
}

upscaled <- pieces %>%
  group_by(n6_id) %>%
  summarise(
    area_analizada_km2 = sum(piece_area) / 1e6,
    n_microcuencas     = n_distinct(comid),
    across(all_of(value_cols), ~ weighted_mean_area(.x, piece_area)),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 5. Attach to the N6 geometries and write
# ------------------------------------------------------------------------------
# The N6 layer drives the join, so all 62 units are kept even if one had no
# overlap at all (its averages would then be NA).

result <- n6 %>%
  left_join(upscaled, by = "n6_id") %>%
  mutate(
    area_analizada_km2 = coalesce(area_analizada_km2, 0),
    n_microcuencas     = coalesce(n_microcuencas, 0L),
    pct_area_analizada = area_analizada_km2 / area_km2 * 100
  ) %>%
  select(
    all_of(n6_fields),
    area_km2, area_analizada_km2, pct_area_analizada, n_microcuencas,
    all_of(value_cols)
  )

coverage <- result$pct_area_analizada   # kept for the report below

if (!include_diagnostics) {
  result <- select(result, -area_km2, -area_analizada_km2,
                   -pct_area_analizada, -n_microcuencas)
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
st_write(result, output_path, delete_dsn = TRUE, quiet = TRUE)

message("Done: ", nrow(result), " N6 units written to ", output_path)
message("  area covered by microbasins: ",
        sprintf("%.1f%% to %.1f%%", min(coverage), max(coverage)))
message("  indice_presion range: ",
        sprintf("[%.4f, %.4f]", min(result$indice_presion, na.rm = TRUE),
                max(result$indice_presion, na.rm = TRUE)))
