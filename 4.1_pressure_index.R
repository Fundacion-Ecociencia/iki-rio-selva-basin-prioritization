# ==============================================================================
# Step 4.1 - Aquatic ecosystem pressure index by micro-watershed
#
# Combines the 15 Mann-Kendall trend signals mapped in 5.1_plots.R into a
# single synthetic pressure index per microbasin, so basins can be ranked
# and prioritised.
#
# Sign orientation
# ------------------------------------------------------------------------------
# Each sub-indicator is re-oriented so that +1 ALWAYS means more pressure on
# aquatic ecosystems:
#
#   ori_i = direction_i * tau_i          (tau_i = NA  ->  ori_i = 0)
#
#   direction = -1  -> plotted with `colorbar`  in 5.1_plots.R, where tau
#                      near -1 means more pressure (e.g. forest loss), so the
#                      sign is INVERTED.
#   direction = +1  -> plotted with `colorbar2` in 5.1_plots.R, where tau
#                      near +1 means more pressure (e.g. mining expansion),
#                      so the sign is KEPT.
#
# NA handling: mk_tau is NA when the series is constant (sd == 0) or too
# short. That is semantically "no trend", so it contributes 0 pressure. This
# keeps the denominator fixed and the index comparable across basins (NA
# rates are high and very uneven: 93.5% for mining, 0.1% for discharge).
# The raw tau_* columns keep their original NA so the substitution stays
# auditable.
#
# Invalid tau: Kendall's tau is bounded to [-1, 1], but trend::mk.test
# returns +-Inf or |tau| > 1 for some heavily tied series in step 1.3 (the
# "pd" metric is the main offender). Those are numerical artefacts rather
# than trends, so they are demoted to NA before anything else and then
# treated like any other constant series.
#
# Aggregation: the mean is taken WITHIN each theme first, then across the
# five themes, so each theme weighs 20% regardless of how many metrics it
# contributes (fragmentation has 4 strongly correlated metrics, population
# only 1).
#
#   tema_t         = mean(ori_i for i in theme t)
#   indice_presion = mean(tema_lulc, tema_fragmentacion, tema_caudal,
#                         tema_deficit, tema_poblacion)
#
# Since every ori_i is in [-1, 1], each tema_t and the final index are also
# bounded by [-1, 1], with +1 = maximum pressure.
#
# Output (workspace/10_pressure_index_by_basin.gpkg), one feature per
# microbasin (polygon geometry from geoglows_basins.gpkg):
#   comid             - Microbasin identifier (linkno).
#   tau_<indicator>   - Raw mk_tau of each of the 15 sub-indicators, as
#                       produced by steps 1.2, 1.3, 2.1, 2.2 and 3.2. NA is
#                       preserved here (constant/short series, or a basin
#                       absent from that source table).
#   ori_<indicator>   - Same value re-oriented towards pressure
#                       (direction * tau) with NA replaced by 0. This is
#                       what actually feeds the index.
#   tema_lulc         - Mean of the 5 land-cover ori_* values.
#   tema_fragmentacion- Mean of the 4 fragmentation ori_* values.
#   tema_caudal       - Mean of the 3 discharge ori_* values.
#   tema_deficit      - Mean of the 2 water-deficit ori_* values.
#   tema_poblacion    - The single population ori_* value.
#   indice_presion    - Final index, mean of the five tema_* columns.
#                       Range [-1, 1]; +1 = highest pressure.
#   n_subind_validos  - Number of sub-indicators with a non-NA raw tau (out
#                       of 15). Indicates how much real data supports the
#                       index in that basin.
# ==============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(tibble)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

basins_path            <- "inputs/geographical_data/geoglows_basins.gpkg"
lulc_trends_path       <- "workspace/02_lulc_trends_by_basin.csv"
frag_trends_path       <- "workspace/05_fragmentation_trends_by_basin.csv"
streamflow_trends_path <- "workspace/06_streamflow_trends_by_basin.csv"
deficit_trends_path    <- "workspace/07_streamflow_deficit_trends.csv"
population_trends_path <- "workspace/09_population_trends_by_basin.csv"

output_path <- "workspace/10_pressure_index_by_basin.gpkg"

# ------------------------------------------------------------------------------
# 2. Sub-indicator registry
# ------------------------------------------------------------------------------
# Mirrors exactly the 15 maps produced at the end of 5.1_plots.R. This table
# is the only place to edit when adding or removing a sub-indicator: the
# extraction, orientation and theme aggregation below are all driven by it.

subindicators <- tribble(
  ~tema,            ~indicador,               ~source,   ~filter_col,  ~filter_value,                 ~direction,
  # --- Land cover (step 1.2) -------------------------------------------------
  "lulc",           "bosque_natural",         "lulc",    "group_code", "A",                                   -1,
  "lulc",           "veg_altoandina",         "lulc",    "group_code", "B",                                   -1,
  "lulc",           "agricultura",            "lulc",    "group_code", "D",                                    1,
  "lulc",           "mineria",                "lulc",    "group_code", "F",                                    1,
  "lulc",           "cuerpos_agua",           "lulc",    "group_code", "H",                                   -1,

  # --- Landscape fragmentation (step 1.3) ------------------------------------
  "fragmentacion",  "num_parches",            "frag",    "metric",     "np",                                   1,
  "fragmentacion",  "densidad_parches",       "frag",    "metric",     "pd",                                   1,
  "fragmentacion",  "parche_mayor",           "frag",    "metric",     "lpi",                                 -1,
  "fragmentacion",  "cohesion",               "frag",    "metric",     "cohesion",                            -1,

  # --- Discharge (step 2.1) --------------------------------------------------
  "caudal",         "caudal_medio",           "flow",    "variable",   "qavg",                                -1,
  "caudal",         "caudal_minimo",          "flow",    "variable",   "qmin",                                -1,
  "caudal",         "caudal_maximo",          "flow",    "variable",   "qmax",                                -1,

  # --- Water deficit (step 2.2) ----------------------------------------------
  "deficit",        "duracion_max_evento",    "deficit", "variable",   "duracion_maxima_evento_dias",          1,
  "deficit",        "frecuencia_eventos",     "deficit", "variable",   "frecuencia_eventos",                   1,

  # --- Population (step 3.2) -------------------------------------------------
  "poblacion",      "densidad_poblacional",   "pop",     "variable",   "pop_density",                          1
)

# ------------------------------------------------------------------------------
# 3. Load the microbasins and the five trend tables
# ------------------------------------------------------------------------------

basins <- st_read(basins_path, quiet = TRUE) %>%
  select(comid = linkno)

# The discharge/deficit tables cover the whole GEOGloWS domain (~12 400
# comid); the join against `basins` at the end keeps only the study area.
sources <- list(
  lulc    = read.csv(lulc_trends_path),
  frag    = read.csv(frag_trends_path),
  flow    = read.csv(streamflow_trends_path),
  deficit = read.csv(deficit_trends_path),
  pop     = read.csv(population_trends_path)
)

# ------------------------------------------------------------------------------
# 4. Extract and orient each sub-indicator
# ------------------------------------------------------------------------------

extract_subindicator <- function(indicador, source, filter_col, filter_value,
                                 direction) {

  df <- sources[[source]]

  if (is.null(df)) {
    stop("Unknown source '", source, "' for sub-indicator '", indicador, "'")
  }
  if (!filter_col %in% names(df)) {
    stop("Column '", filter_col, "' not found in source '", source,
         "' (sub-indicator '", indicador, "')")
  }
  if (!"mk_tau" %in% names(df)) {
    stop("Column 'mk_tau' not found in source '", source, "'")
  }

  rows <- df[df[[filter_col]] == filter_value, , drop = FALSE]

  if (nrow(rows) == 0) {
    stop("No rows found for ", filter_col, " = '", filter_value,
         "' in source '", source, "' (sub-indicator '", indicador, "')")
  }

  tau <- suppressWarnings(as.numeric(rows$mk_tau))

  # Kendall's tau is bounded to [-1, 1] by definition. Some upstream series
  # are degenerate enough (heavily tied values) that trend::mk.test returns
  # +-Inf or |tau| > 1 from its tie-corrected denominator - notably the
  # fragmentation "pd" metric. Those are numerical artefacts, not trends, so
  # they are demoted to NA and treated exactly like a constant series.
  invalid <- !is.na(tau) & (!is.finite(tau) | tau < -1 | tau > 1)
  if (any(invalid)) {
    message("  ", indicador, ": ", sum(invalid),
            " invalid tau values (non-finite or outside [-1, 1]) set to NA")
    tau[invalid] <- NA_real_
  }

  tibble(
    comid     = rows$comid,
    indicador = indicador,
    tau       = tau,
    # Orient towards pressure; a constant/short series (NA) adds no pressure
    ori       = direction * coalesce(tau, 0)
  )
}

long <- lapply(seq_len(nrow(subindicators)), function(i) {
  s <- subindicators[i, ]
  extract_subindicator(s$indicador, s$source, s$filter_col, s$filter_value,
                       s$direction)
}) %>%
  bind_rows()

# ------------------------------------------------------------------------------
# 5. Reshape to one row per microbasin
# ------------------------------------------------------------------------------

wide <- long %>%
  pivot_wider(
    id_cols     = comid,
    names_from  = indicador,
    values_from = c(tau, ori),
    names_glue  = "{.value}_{indicador}"
  )

# Anchor on the full basin set so every microbasin gets a row, then fill the
# ori_* gaps left by basins missing from a source table (5 lack fragmentation
# metrics, 12 lack discharge/deficit series) with 0 = no pressure trend.
# The tau_* columns keep their NA so the substitution remains visible.
wide <- tibble(comid = basins$comid) %>%
  left_join(wide, by = "comid") %>%
  mutate(across(starts_with("ori_"), ~ coalesce(.x, 0)))

# ------------------------------------------------------------------------------
# 6. Aggregate: within each theme first, then across themes
# ------------------------------------------------------------------------------

temas <- unique(subindicators$tema)

for (tm in temas) {
  ori_cols <- paste0("ori_",
                     subindicators$indicador[subindicators$tema == tm])
  wide[[paste0("tema_", tm)]] <-
    rowMeans(as.matrix(wide[, ori_cols, drop = FALSE]))
}

tema_cols <- paste0("tema_", temas)
tau_cols  <- paste0("tau_", subindicators$indicador)
ori_cols  <- paste0("ori_", subindicators$indicador)

wide$indice_presion <- rowMeans(as.matrix(wide[, tema_cols, drop = FALSE]))

# Data support: how many of the 15 sub-indicators had a real (non-NA) trend
wide$n_subind_validos <-
  rowSums(!is.na(as.matrix(wide[, tau_cols, drop = FALSE])))

# ------------------------------------------------------------------------------
# 7. Write the geopackage
# ------------------------------------------------------------------------------

index_sf <- basins %>%
  left_join(
    wide %>%
      select(comid, all_of(tau_cols), all_of(ori_cols), all_of(tema_cols),
             indice_presion, n_subind_validos),
    by = "comid"
  )

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
st_write(index_sf, output_path, delete_dsn = TRUE, quiet = TRUE)

message("Done: ", nrow(index_sf), " microbasins written to ", output_path)
message("  indice_presion range: [",
        round(min(index_sf$indice_presion, na.rm = TRUE), 4), ", ",
        round(max(index_sf$indice_presion, na.rm = TRUE), 4), "]")
