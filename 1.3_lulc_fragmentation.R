# ==============================================================================
# Step 1.3 - Landscape fragmentation of natural cover by micro-watershed
#
# Computes class-level fragmentation and connectivity metrics for natural cover
# in each micro-watershed and year, and evaluates their temporal trend.
#
# Land cover is reduced to a binary landscape:
#   1 = natural cover      (groups A, B, C, G, H)
#   2 = anthropic use      (groups D, E, F, I)
# Both classes are kept so that metrics are computed on a complete landscape;
# results are then filtered to class 1.
#
# Metrics reported (class level, natural cover):
#   np       - number of patches
#   pd       - patch density (patches / 100 ha)
#   area_mn  - mean patch area (ha)
#   lpi      - largest patch index (%)
#   ed       - edge density (m / ha)
#   cohesion - patch cohesion index (0-100), physical connectedness
#   mesh     - effective mesh size (ha), inverse measure of fragmentation
#
# Outputs:
#   fragmentation_by_microbasin.csv        -> one row per comid x year x metric
#   fragmentation_trends_by_microbasin.csv -> one row per comid x metric
#
# Output columns
# ------------------------------------------------------------------------------
# 04_fragmentation_by_basin.csv (one row per comid x year x metric, class 1
# "natural cover" only):
#   comid   - Microbasin identifier (linkno).
#   year    - Year of the MapBiomas annual raster (1985-2024).
#   metric  - Metric code, one of np/pd/area_mn/lpi/ed/cohesion/mesh (see
#             "Metrics reported" above for the definition of each).
#   value   - Metric value for that comid/year, in the unit stated above for
#             the metric (patches, patches/100 ha, ha, %, m/ha, 0-100 index).
#             Can be NA (natural cover absent that year) or Inf (e.g. ratio
#             metrics with a zero denominator when a single patch remains);
#             these are excluded before computing trends in step 5.
#
# 05_fragmentation_trends_by_basin.csv (one row per comid x metric, summarizing
# the full 1985-2024 series of the metric above):
#   comid            - Microbasin identifier (linkno).
#   metric           - Metric code (np/pd/area_mn/lpi/ed/cohesion/mesh).
#   n_years          - Number of years with a finite (non-NA/non-Inf) value,
#                       used in the trend test.
#   value_mean       - Mean of the finite values across the series; NA if no
#                       finite value exists.
#   mk_tau           - Mann-Kendall tau statistic (non-parametric trend
#                       strength/direction, -1 to 1); NA if n_years < 5 or the
#                       series is constant.
#   mk_pvalue        - p-value of the Mann-Kendall test; NA under the same
#                       conditions as mk_tau.
#   sen_slope        - Sen's slope estimator (robust, non-parametric trend),
#                       in metric units per year; NA under the same
#                       conditions as mk_tau.
#   trend_direction  - "incremento" / "descenso" when the Mann-Kendall test is
#                       significant (mk_pvalue < 0.05) and sen_slope is
#                       positive/negative respectively; "sin tendencia"
#                       otherwise (including series with mk_pvalue = NA).
# ==============================================================================

library(terra)
library(sf)
library(landscapemetrics)
library(dplyr)
library(tidyr)
library(trend)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

raster_dir  <- "inputs/map_biomas"
basins_path <- "inputs/geographical_data/geoglows_basins.gpkg"

metrics_output <- "workspace/04_fragmentation_by_basin.csv"
trends_output  <- "workspace/05_fragmentation_trends_by_basin.csv"

# Metric CRS required by landscapemetrics (square cells, metric units)
target_crs <- "EPSG:32718"
target_crs <- "+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs"

selected_metrics <- c("np", "pd", "area_mn", "lpi", "ed", "cohesion", "mesh")

# ------------------------------------------------------------------------------
# 2. Reclassification: MapBiomas code -> binary landscape
# ------------------------------------------------------------------------------
# 1 = natural cover, 2 = anthropic use

reclass_matrix <- as.matrix(tribble(
  ~code, ~class,
  3,  1,   # Bosque natural
  4,  1,
  6,  1,
  81, 1,   # Vegetacion altoandina
  82, 1,
  11, 1,
  12, 1,   # Otra formacion natural
  13, 1,
  29, 1,
  23, 1,   # Area no vegetada natural
  68, 1,
  33, 1,   # Cuerpos de agua
  34, 1,
  9,  2,   # Agricultura y silvicultura
  21, 2,
  74, 2,
  24, 2,   # Infraestructura antropica
  25, 2,
  30, 2,   # Mineria
  31, 2    # Acuicultura
))

# ------------------------------------------------------------------------------
# 3. Load micro-watersheds
# ------------------------------------------------------------------------------

basins <- st_read(basins_path, quiet = TRUE) %>%
  select(comid = linkno) %>%
  st_transform(target_crs)

# sample_lsm() identifies polygons by their row position
basin_ids <- tibble(plot_id = seq_len(nrow(basins)), comid = basins$comid)

# ------------------------------------------------------------------------------
# 4. Compute metrics for each year
# ------------------------------------------------------------------------------

raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)

metrics_for_year <- function(file_path) {

  year <- as.integer(tools::file_path_sans_ext(basename(file_path)))
  message("Processing ", year)

  # Reclassify first, then reproject: nearest-neighbour on two classes
  # preserves patch structure better than reprojecting 20 original classes
  r <- classify(rast(file_path), reclass_matrix, others = NA) %>%
    project(target_crs, method = "near")

  sample_lsm(
    r,
    y      = basins,
    level  = "class",
    metric = selected_metrics,
    progress = FALSE
  ) %>%
    filter(class == 1) %>%           # natural cover only
    select(plot_id, metric, value) %>%
    mutate(year = year)
}

fragmentation <- lapply(raster_files, metrics_for_year) %>%
  bind_rows() %>%
  left_join(basin_ids, by = "plot_id") %>%
  select(comid, year, metric, value) %>%
  arrange(comid, metric, year)

write.csv(fragmentation, metrics_output, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Temporal trend of each metric (Mann-Kendall + Sen's slope)
# ------------------------------------------------------------------------------

trend_stats <- function(df) {

  # Keep only finite values: metrics can be NA (class absent) or Inf
  # (e.g. ratios with a zero denominator when a single patch remains)
  values <- df$value[is.finite(df$value)]

  # A trend test requires a minimum series length and some variability
  if (length(values) < 5 || sd(values) == 0) {
    return(tibble(
      n_years    = length(values),
      value_mean = if (length(values) > 0) mean(values) else NA_real_,
      mk_tau     = NA_real_,
      mk_pvalue  = NA_real_,
      sen_slope  = NA_real_
    ))
  }

  mk  <- mk.test(values)
  sen <- sens.slope(values)

  tibble(
    n_years    = length(values),
    value_mean = mean(values),
    mk_tau     = as.numeric(mk$estimates["tau"]),
    mk_pvalue  = mk$p.value,
    sen_slope  = as.numeric(sen$estimates)
  )
}




trends <- fragmentation %>%
  group_by(comid, metric) %>%
  group_modify(~ trend_stats(.x)) %>%
  ungroup() %>%
  mutate(
    trend_direction = case_when(
      is.na(mk_pvalue) | mk_pvalue >= 0.05 ~ "sin tendencia",
      sen_slope > 0                        ~ "incremento",
      TRUE                                 ~ "descenso"
    )
  ) %>%
  arrange(comid, metric)

write.csv(trends, trends_output, row.names = FALSE)

message("Done: ", nrow(fragmentation), " metric rows and ",
        nrow(trends), " trend rows written")
