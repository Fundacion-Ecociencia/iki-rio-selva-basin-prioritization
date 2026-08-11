# ==============================================================================
# Step 5.1 - Mapping of trend indicators and the pressure index by microbasin
#
# Provides six mapping functions to visualize any column of the trend tables
# produced in steps 1.2, 1.3, 2.1, 2.2 and 3.2, plus the synthetic pressure
# index built in steps 4.1 / 4.2. Each one resolves a filter and a title and
# then delegates the drawing to map_choropleth():
#
#   plot_lulc_trend_map()          -> maps a column of
#                                      02_lulc_trends_by_basin.csv, filtered
#                                      to one functional group (group_code:
#                                      A-I, or NAT/ANT if step 1.2 was run
#                                      with include_aggregates = TRUE).
#
#   plot_fragmentation_trend_map() -> maps a column of
#                                      05_fragmentation_trends_by_basin.csv,
#                                      filtered to one landscape metric
#                                      (metric: np, pd, area_mn, lpi, ed,
#                                      cohesion or mesh).
#
#   plot_streamflow_trend_map()    -> maps a column of
#                                      06_streamflow_trends_by_basin.csv,
#                                      filtered to one streamflow variable
#                                      (variable: qavg = "Caudal medio",
#                                      qmin = "Caudal mínimo", qmax =
#                                      "Caudal máximo").
#
#   plot_deficit_trend_map()       -> maps a column of
#                                      07_streamflow_deficit_trends.csv,
#                                      filtered to one deficit variable
#                                      (variable: dias_bajo_umbral,
#                                      duracion_maxima_evento_dias,
#                                      duracion_media_eventos_dias,
#                                      frecuencia_eventos,
#                                      porcentaje_dias_bajo_umbral).
#
#   plot_population_trend_map()    -> maps a column of
#                                      09_population_trends_by_basin.csv,
#                                      filtered to one population variable
#                                      (variable: population = "Poblacion
#                                      total", pop_density = "Densidad
#                                      poblacional").
#
#   plot_pressure_index_map()      -> maps a column of the geopackage
#                                      10_pressure_index_by_basin.gpkg built
#                                      in step 4.1: the synthetic pressure
#                                      index ("indice_presion", the default),
#                                      any theme score ("tema_*") or any
#                                      individual sub-indicator ("tau_*" /
#                                      "ori_*"). No filter needed, since that
#                                      table already holds one row per
#                                      microbasin.
#
# Both functions join the requested trend table to geoglows_basins.gpkg by
# comid (= linkno) and render a choropleth map with ggplot2:
#   - Microbasin fill AND border are both colored by `value_col`, using the
#     same color ramp (so borders blend into their own polygon rather than
#     standing out in a fixed color).
#   - The color ramp is passed in via the `colors` argument (a vector of
#     colors, in the order they should appear along the scale). It defaults
#     to a diverging blue -> green -> white -> orange -> red ramp, centered
#     on 0 by default (blue = negative/loss, white = no change,
#     red = positive/increase), matching the sign convention of ols_slope,
#     sen_slope and mk_tau in steps 1.2 and 1.3. Columns without a
#     meaningful zero (e.g. pct_mean, value_mean, n_years) should be mapped
#     with center_zero = FALSE. Non-numeric columns (e.g. trend_direction)
#     use `colors` as a discrete palette instead.
#   - The outline of the overall study basin (pastaza.gpkg) is drawn on top
#     as a fixed-color reference boundary.
#   - Province boundaries (provinces.gpkg) are drawn on top of everything
#     else, with transparent fill and a thin (linewidth 0.5) grey outline;
#     toggle with `show_provinces`.
#   - The map extent (zoom) is fixed to the bounding box of
#     geoglows_basins.gpkg, so overlaying the (typically larger) study basin
#     or province layers never zooms the map out to fit them. `padding`
#     (map units, meters) is added to each side of that bounding box before
#     fixing the extent, to leave breathing room around the basins.
#   - The legend is a horizontal color bar (or horizontal discrete legend)
#     placed below the map, and the plot panel has a black border frame.
#   - If `output_path` is given, the figure is also saved to that PNG file,
#     at the requested size in centimeters.
# ==============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

basins_path      <- "inputs/geographical_data/geoglows_basins.gpkg"
study_basin_path <- "inputs/geographical_data/pastaza.gpkg"
provinces_path   <- "inputs/geographical_data/provinces.gpkg"
lulc_trends_path <- "workspace/02_lulc_trends_by_basin.csv"
frag_trends_path <- "workspace/05_fragmentation_trends_by_basin.csv"
streamflow_trends_path <- "workspace/06_streamflow_trends_by_basin.csv"
deficit_trends_path <- "workspace/07_streamflow_deficit_trends.csv"
population_trends_path <- "workspace/09_population_trends_by_basin.csv"
pressure_index_path <- "workspace/10_pressure_index_by_basin.gpkg"
pressure_index_v2_path <- "workspace/11_pressure_index_v2_by_basin.gpkg"

# Spanish labels for the streamflow "variable" column values, used to build
# the map title in plot_streamflow_trend_map().
streamflow_variable_labels <- c(
  qavg = "Caudal medio",
  qmin = "Caudal mínimo",
  qmax = "Caudal máximo"
)

# Spanish labels for the population "variable" column values, used to build
# the map title in plot_population_trend_map().
population_variable_labels <- c(
  population  = "Poblacion total",
  pop_density = "Densidad poblacional"
)

basins <- st_read(basins_path, quiet = TRUE) %>%
  select(comid = linkno)

# Outline of the overall study basin, dissolved to a single boundary in case
# the source layer holds more than one feature, reprojected to match basins.
study_basin <- st_read(study_basin_path, quiet = TRUE) %>%
  st_transform(st_crs(basins)) %>%
  st_union() %>%
  st_as_sf()

# Fixed map extent, so that adding the province reference layer never zooms
# the map out to fit it. `padding` (in map_choropleth) adds breathing room
# around this box.
basins_bbox <- st_bbox(study_basin)

# Province boundaries, drawn on top of everything as a reference layer
# (transparent fill, thin grey outline).
provinces <- st_read(provinces_path, quiet = TRUE) %>%
  st_transform(st_crs(basins))

# ------------------------------------------------------------------------------
# 2. Shared helpers
# ------------------------------------------------------------------------------

# Default diverging blue -> green -> white -> orange -> red ramp
diverging_palette <- c("#08519c", "#74c476", "#ffffff", "#fd8d3c", "#a50f15")

# Renders the choropleth: `value_col` drives both fill and border color,
# using `colors` as the ramp (continuous gradient for numeric columns,
# discrete palette for categorical ones such as trend_direction). The study
# basin outline is drawn on top, and the figure is optionally saved to
# `output_path` at `width_cm` x `height_cm`.
map_choropleth <- function(data, value_col, title, legend_title, colors,
                           center_zero = TRUE,
                           limits = NULL,
                           study_boundary_color = "black",
                           show_provinces = TRUE,
                           padding = 1000,
                           output_path = NULL,
                           width_cm = 20,
                           height_cm = 15) {

  values <- data[[value_col]]

  base_map <- ggplot(data) +
    geom_sf(aes(fill = .data[[value_col]], color = .data[[value_col]]),
            linewidth = 0.15) +
    geom_sf(data = study_basin, fill = NA,
            color = study_boundary_color, linewidth = 0.7)

  if (show_provinces) {
    base_map <- base_map +
      geom_sf(data = provinces, fill = NA, color = "black", linewidth = 0.3)
  }

  base_map <- base_map +
    coord_sf(
      xlim   = c(basins_bbox["xmin"] - padding, basins_bbox["xmax"] + padding),
      ylim   = c(basins_bbox["ymin"] - padding, basins_bbox["ymax"] + padding),
      expand = FALSE
    ) +
    labs(title = title, fill = legend_title) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5)
    )

  if (is.numeric(values)) {

    if (is.null(limits)) {
      if (center_zero) {
        max_abs <- max(abs(values), na.rm = TRUE)
        limits  <- c(-max_abs, max_abs)
      } else {
        limits <- range(values, na.rm = TRUE)
      }
    }

    map <- base_map +
      scale_fill_gradientn(
        colours  = colors,
        limits   = limits,
        oob      = scales::squish,
        na.value = "grey85",
        guide    = guide_colorbar(
          direction      = "horizontal",
          title.position = "top",
          title.hjust    = 0.5,
          barwidth       = unit(6, "cm"),
          barheight      = unit(0.4, "cm")
        )
      ) +
      scale_color_gradientn(
        colours  = colors,
        limits   = limits,
        oob      = scales::squish,
        na.value = "grey85",
        guide    = "none"
      )

  } else {

    map <- base_map +
      scale_fill_manual(
        values = colors, na.value = "grey85",
        guide  = guide_legend(
          direction      = "horizontal",
          title.position = "top",
          title.hjust    = 0.5,
          nrow            = 1
        )
      ) +
      scale_color_manual(values = colors, na.value = "grey85", guide = "none")
  }

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(output_path, plot = map, width = width_cm, height = height_cm,
           units = "cm", dpi = 300)
  }

  map
}

# ------------------------------------------------------------------------------
# 3. Map function - LULC trends (step 1.2 output)
# ------------------------------------------------------------------------------

#' Map a column of 02_lulc_trends_by_basin.csv for one functional group
#'
#' @param group_code_filter Group code to filter to: "A".."I" (see
#'   group_lookup in step 1.2), or "NAT"/"ANT" if step 1.2 was run with
#'   include_aggregates = TRUE.
#' @param value_col   Name of the column to map, e.g. "sen_slope",
#'   "ols_slope", "mk_tau", "pct_mean" or "trend_direction".
#' @param legend_label Text shown as the legend title, e.g. "Pendiente tau"
#'   for value_col = "mk_tau". If NULL (default), `value_col` itself is used.
#' @param colors      Vector of colors used to build the scale (continuous
#'   gradient for numeric columns, discrete palette for categorical ones).
#'   Defaults to a diverging blue-green-white-orange-red ramp.
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Set to FALSE for columns without a meaningful zero,
#'   such as "pct_mean".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range, e.g. c(-1, 1). If NULL (default), the range is
#'   derived automatically from the data via `center_zero`. Values outside
#'   `limits` are squished to the nearest end of the scale rather than
#'   turned into NA.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param trends_path Path to the trends CSV (default: step 1.2 output).
#' @return A ggplot object.
plot_lulc_trend_map <- function(group_code_filter,
                                value_col,
                                legend_label = NULL,
                                colors = diverging_palette,
                                center_zero = TRUE,
                                limits = NULL,
                                show_provinces = TRUE,
                                padding = 1000,
                                output_path = NULL,
                                width_cm = 20,
                                height_cm = 15,
                                trends_path = lulc_trends_path) {

  trends <- read.csv(trends_path) %>%
    filter(group_code == group_code_filter)

  if (nrow(trends) == 0) {
    stop("No rows found for group_code = '", group_code_filter,
         "' in ", trends_path)
  }
  if (!value_col %in% names(trends)) {
    stop("Column '", value_col, "' not found in ", trends_path)
  }

  map_data <- basins %>%
    left_join(trends, by = "comid")

  group_label <- unique(trends$group)[1]

  if (is.null(legend_label)) legend_label <- value_col

  map_choropleth(
    data         = map_data,
    value_col    = value_col,
    title        = paste0(
                          group_label,
                          " (", group_code_filter, ") - ",
                          legend_label),
    legend_title = legend_label,
    colors         = colors,
    center_zero    = center_zero,
    limits         = limits,
    show_provinces = show_provinces,
    padding        = padding,
    output_path    = output_path,
    width_cm       = width_cm,
    height_cm      = height_cm
  )
}

# ------------------------------------------------------------------------------
# 4. Map function - fragmentation trends (step 1.3 output)
# ------------------------------------------------------------------------------

#' Map a column of 05_fragmentation_trends_by_basin.csv for one metric
#'
#' @param metric_filter Landscape metric to filter to: "np", "pd",
#'   "area_mn", "lpi", "ed", "cohesion" or "mesh" (see step 1.3 for the
#'   definition of each).
#' @param value_col   Name of the column to map, e.g. "sen_slope",
#'   "mk_tau", "value_mean" or "trend_direction".
#' @param legend_label Text shown as the legend title, e.g. "Pendiente tau"
#'   for value_col = "mk_tau". If NULL (default), `value_col` itself is used.
#' @param colors      Vector of colors used to build the scale (continuous
#'   gradient for numeric columns, discrete palette for categorical ones).
#'   Defaults to a diverging blue-green-white-orange-red ramp.
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Set to FALSE for columns without a meaningful zero,
#'   such as "value_mean" or "n_years".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range, e.g. c(-1, 1). If NULL (default), the range is
#'   derived automatically from the data via `center_zero`. Values outside
#'   `limits` are squished to the nearest end of the scale rather than
#'   turned into NA.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param trends_path Path to the trends CSV (default: step 1.3 output).
#' @return A ggplot object.
plot_fragmentation_trend_map <- function(metric_filter,
                                         value_col,
                                         legend_label = NULL,
                                         colors = diverging_palette,
                                         center_zero = TRUE,
                                         limits = NULL,
                                         show_provinces = TRUE,
                                         padding = 1000,
                                         output_path = NULL,
                                         width_cm = 20,
                                         height_cm = 15,
                                         trends_path = frag_trends_path) {

  trends <- read.csv(trends_path) %>%
    filter(metric == metric_filter)

  if (nrow(trends) == 0) {
    stop("No rows found for metric = '", metric_filter, "' in ", trends_path)
  }
  if (!value_col %in% names(trends)) {
    stop("Column '", value_col, "' not found in ", trends_path)
  }

  map_data <- basins %>%
    left_join(trends, by = "comid")

  if (is.null(legend_label)) legend_label <- value_col

  map_choropleth(
    data = map_data,
    value_col = value_col,
    title = paste0("Fragmentacion - ", metric_filter, " - ", legend_label),
    legend_title = legend_label,
    colors = colors,
    center_zero = center_zero,
    limits = limits,
    show_provinces = show_provinces,
    padding = padding,
    output_path = output_path,
    width_cm = width_cm,
    height_cm = height_cm
  )
}

# ------------------------------------------------------------------------------
# 5. Map function - streamflow trends (step 2.1 output)
# ------------------------------------------------------------------------------

#' Map a column of 06_streamflow_trends_by_basin.csv for one streamflow
#' variable
#'
#' @param variable_filter Streamflow variable to filter to: "qavg" (mapped
#'   to the label "Caudal medio"), "qmin" ("Caudal mínimo") or "qmax"
#'   ("Caudal máximo").
#' @param value_col   Name of the column to map, e.g. "sen_slope", "mk_tau",
#'   "q_mean", "sen_slope_rel" or "trend_direction".
#' @param legend_label Text shown as the legend title, e.g. "Pendiente tau"
#'   for value_col = "mk_tau". If NULL (default), `value_col` itself is used.
#' @param colors      Vector of colors used to build the scale (continuous
#'   gradient for numeric columns, discrete palette for categorical ones).
#'   Defaults to a diverging blue-green-white-orange-red ramp.
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Set to FALSE for columns without a meaningful zero,
#'   such as "q_mean" or "n_years".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range, e.g. c(-1, 1). If NULL (default), the range is
#'   derived automatically from the data via `center_zero`. Values outside
#'   `limits` are squished to the nearest end of the scale rather than
#'   turned into NA.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param trends_path Path to the trends CSV (default: step 2.1 output).
#' @return A ggplot object.
plot_streamflow_trend_map <- function(variable_filter,
                                      value_col,
                                      legend_label = NULL,
                                      colors = diverging_palette,
                                      center_zero = TRUE,
                                      limits = NULL,
                                      show_provinces = TRUE,
                                      padding = 1000,
                                      output_path = NULL,
                                      width_cm = 20,
                                      height_cm = 15,
                                      trends_path = streamflow_trends_path) {

  trends <- read.csv(trends_path) %>%
    filter(variable == variable_filter)

  if (nrow(trends) == 0) {
    stop("No rows found for variable = '", variable_filter, "' in ", trends_path)
  }
  if (!value_col %in% names(trends)) {
    stop("Column '", value_col, "' not found in ", trends_path)
  }

  map_data <- basins %>%
    left_join(trends, by = "comid")

  variable_label <- unname(streamflow_variable_labels[variable_filter])
  if (is.na(variable_label)) variable_label <- variable_filter

  if (is.null(legend_label)) legend_label <- value_col

  map_choropleth(
    data = map_data,
    value_col = value_col,
    title = paste0(variable_label, " - ", legend_label),
    legend_title = legend_label,
    colors = colors,
    center_zero = center_zero,
    limits = limits,
    show_provinces = show_provinces,
    padding = padding,
    output_path = output_path,
    width_cm = width_cm,
    height_cm = height_cm
  )
}

# ------------------------------------------------------------------------------
# 6. Map function - streamflow deficit trends (step 2.2 output)
# ------------------------------------------------------------------------------

#' Map a column of 07_streamflow_deficit_trends.csv for one deficit variable
#'
#' @param variable_filter Deficit variable to filter to: "dias_bajo_umbral",
#'   "duracion_maxima_evento_dias", "duracion_media_eventos_dias",
#'   "frecuencia_eventos" or "porcentaje_dias_bajo_umbral".
#' @param value_col   Name of the column to map, e.g. "sen_slope", "mk_tau",
#'   "value_mean", "prop_zeros" or "trend_direction".
#' @param legend_label Text shown as the legend title, e.g. "Pendiente tau"
#'   for value_col = "mk_tau". If NULL (default), `value_col` itself is used.
#' @param colors      Vector of colors used to build the scale (continuous
#'   gradient for numeric columns, discrete palette for categorical ones).
#'   Defaults to a diverging blue-green-white-orange-red ramp.
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Set to FALSE for columns without a meaningful zero,
#'   such as "value_mean", "prop_zeros" or "n_years".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range, e.g. c(-1, 1). If NULL (default), the range is
#'   derived automatically from the data via `center_zero`. Values outside
#'   `limits` are squished to the nearest end of the scale rather than
#'   turned into NA.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param trends_path Path to the trends CSV (default: step 2.2 output).
#' @return A ggplot object.
plot_deficit_trend_map <- function(variable_filter,
                                   value_col,
                                   legend_label = NULL,
                                   colors = diverging_palette,
                                   center_zero = TRUE,
                                   limits = NULL,
                                   show_provinces = TRUE,
                                   padding = 1000,
                                   output_path = NULL,
                                   width_cm = 20,
                                   height_cm = 15,
                                   trends_path = deficit_trends_path) {

  trends <- read.csv(trends_path) %>%
    filter(variable == variable_filter)

  if (nrow(trends) == 0) {
    stop("No rows found for variable = '", variable_filter, "' in ", trends_path)
  }
  if (!value_col %in% names(trends)) {
    stop("Column '", value_col, "' not found in ", trends_path)
  }

  map_data <- basins %>%
    left_join(trends, by = "comid")

  if (is.null(legend_label)) legend_label <- value_col

  map_choropleth(
    data = map_data,
    value_col = value_col,
    title = paste0(variable_filter, " - ", legend_label),
    legend_title = legend_label,
    colors = colors,
    center_zero = center_zero,
    limits = limits,
    show_provinces = show_provinces,
    padding = padding,
    output_path = output_path,
    width_cm = width_cm,
    height_cm = height_cm
  )
}

# ------------------------------------------------------------------------------
# 7. Map function - population trends (step 3.2 output)
# ------------------------------------------------------------------------------

#' Map a column of 09_population_trends_by_basin.csv for one population
#' variable
#'
#' @param variable_filter Population variable to filter to: "population"
#'   (mapped to the label "Poblacion total") or "pop_density" ("Densidad
#'   poblacional").
#' @param value_col   Name of the column to map, e.g. "sen_slope", "mk_tau",
#'   "value_mean", "sen_slope_rel" or "trend_direction".
#' @param legend_label Text shown as the legend title, e.g. "Pendiente tau"
#'   for value_col = "mk_tau". If NULL (default), `value_col` itself is used.
#' @param colors      Vector of colors used to build the scale (continuous
#'   gradient for numeric columns, discrete palette for categorical ones).
#'   Defaults to a diverging blue-green-white-orange-red ramp.
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Set to FALSE for columns without a meaningful zero,
#'   such as "value_mean" or "n_years".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range, e.g. c(-1, 1). If NULL (default), the range is
#'   derived automatically from the data via `center_zero`. Values outside
#'   `limits` are squished to the nearest end of the scale rather than
#'   turned into NA.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param trends_path Path to the trends CSV (default: step 3.2 output).
#' @return A ggplot object.
plot_population_trend_map <- function(variable_filter,
                                      value_col,
                                      legend_label = NULL,
                                      colors = diverging_palette,
                                      center_zero = TRUE,
                                      limits = NULL,
                                      show_provinces = TRUE,
                                      padding = 1000,
                                      output_path = NULL,
                                      width_cm = 20,
                                      height_cm = 15,
                                      trends_path = population_trends_path) {

  trends <- read.csv(trends_path) %>%
    filter(variable == variable_filter)

  if (nrow(trends) == 0) {
    stop("No rows found for variable = '", variable_filter, "' in ", trends_path)
  }
  if (!value_col %in% names(trends)) {
    stop("Column '", value_col, "' not found in ", trends_path)
  }

  map_data <- basins %>%
    left_join(trends, by = "comid")

  variable_label <- unname(population_variable_labels[variable_filter])
  if (is.na(variable_label)) variable_label <- variable_filter

  if (is.null(legend_label)) legend_label <- value_col

  map_choropleth(
    data = map_data,
    value_col = value_col,
    title = paste0(variable_label, " - ", legend_label),
    legend_title = legend_label,
    colors = colors,
    center_zero = center_zero,
    limits = limits,
    show_provinces = show_provinces,
    padding = padding,
    output_path = output_path,
    width_cm = width_cm,
    height_cm = height_cm
  )
}

# ------------------------------------------------------------------------------
# 8. Map function - synthetic pressure index (step 4.1 / 4.2 output)
# ------------------------------------------------------------------------------

#' Map a column of 10_pressure_index_by_basin.gpkg
#'
#' Unlike the functions above, this one reads a geopackage that already holds
#' one row per microbasin, so it takes no filter argument. Its attributes are
#' joined onto `basins` (rather than plotted from the gpkg geometry directly)
#' so the map keeps exactly the same extent and layer stack as every other
#' figure in this script.
#'
#' @param value_col   Name of the column to map. Defaults to
#'   "indice_presion" (the final index). Also accepts any theme score
#'   ("tema_lulc", "tema_fragmentacion", "tema_caudal", "tema_deficit",
#'   "tema_poblacion"), any pressure-oriented sub-indicator ("ori_*"), any
#'   raw Mann-Kendall tau ("tau_*") or "n_subind_validos".
#' @param legend_label Text shown as the legend title, e.g. "Indice de
#'   presion". If NULL (default), `value_col` itself is used.
#' @param title       Map title. If NULL (default), `legend_label` is used.
#' @param colors      Vector of colors used to build the scale. Note the
#'   index is already oriented so that +1 = highest pressure, which matches
#'   `colorbar2` (blue -> white -> red).
#' @param center_zero Center the color scale at 0 (default TRUE). Ignored if
#'   `limits` is given. Keep TRUE for the index and any tau/ori/tema column;
#'   set FALSE for "n_subind_validos".
#' @param limits      Optional length-2 numeric vector fixing the legend/
#'   color scale range. Unlike the raw mk_tau maps, do NOT force c(-1, 1)
#'   here: averaging 15 signals (many of them 0) pulls the index towards 0
#'   and the map would come out almost entirely white. Leaving it NULL makes
#'   the scale span the actual data range, symmetric around 0.
#' @param show_provinces Draw province boundaries (provinces.gpkg) on top of
#'   the map, as a transparent-fill layer with a thin grey outline (default
#'   TRUE).
#' @param padding     Padding added around the microbasins bounding box, in
#'   map units (meters, since basins are in a UTM CRS), applied to each of
#'   xmin/xmax/ymin/ymax before fixing the map extent (default 1000).
#' @param output_path Optional path to save the figure as PNG. If NULL
#'   (default), the figure is not saved to disk.
#' @param width_cm    Figure width in centimeters when saving (default 20).
#' @param height_cm   Figure height in centimeters when saving (default 15).
#' @param index_path  Path to the pressure index geopackage (default: step
#'   4.1 output).
#' @return A ggplot object.
plot_pressure_index_map <- function(value_col = "indice_presion",
                                    legend_label = NULL,
                                    title = NULL,
                                    colors = diverging_palette,
                                    center_zero = TRUE,
                                    limits = NULL,
                                    show_provinces = TRUE,
                                    padding = 1000,
                                    output_path = NULL,
                                    width_cm = 20,
                                    height_cm = 15,
                                    index_path = pressure_index_path) {

  if (!file.exists(index_path)) {
    stop("Pressure index not found at '", index_path,
         "'. Run 4.1_pressure_index.R first.")
  }

  index <- st_read(index_path, quiet = TRUE) %>%
    st_drop_geometry()

  if (!value_col %in% names(index)) {
    stop("Column '", value_col, "' not found in ", index_path)
  }

  map_data <- basins %>%
    left_join(index, by = "comid")

  if (is.null(legend_label)) legend_label <- value_col
  if (is.null(title)) title <- legend_label

  map_choropleth(
    data = map_data,
    value_col = value_col,
    title = title,
    legend_title = legend_label,
    colors = colors,
    center_zero = center_zero,
    limits = limits,
    show_provinces = show_provinces,
    padding = padding,
    output_path = output_path,
    width_cm = width_cm,
    height_cm = height_cm
  )
}

# ------------------------------------------------------------------------------
# 9. Example usage
# ------------------------------------------------------------------------------

colorbar <- c(
  "#f4662e", "#f7f7f7", "#2166ac"
)

colorbar2 <- c(
  "#2166ac", "#f7f7f7", "#f4662e"
)

colorbar3 <- c(
  "#2166ac", "#1c9c71", "#f7f7f7", "#f4d32e", "#f4662e"
)

plot_lulc_trend_map(
  "A", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/01_mk_tau_bosque_natural.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_lulc_trend_map(
  "B", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/02_mk_tau_vegetacion_altoandina.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)


plot_lulc_trend_map(
  "D", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/03_mk_tau_agricultura.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_lulc_trend_map(
  "F", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/04_mk_tau_mineria.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_lulc_trend_map(
  "H", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/05_mk_tau_agua.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_fragmentation_trend_map(
  "np", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/06_mk_tau_number_patch.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_fragmentation_trend_map(
  "pd", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/07_mk_tau_patch_density.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_fragmentation_trend_map(
  "lpi", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/08_mk_tau_patch_area.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_fragmentation_trend_map(
  "cohesion", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/09_mk_tau_cohesion.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)


plot_streamflow_trend_map(
  "qavg", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/10_mk_tau_qavg.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_streamflow_trend_map(
  "qmin", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/11_mk_tau_qmin.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_streamflow_trend_map(
  "qmax", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar,
  limits       = c(-1, 1),
  output_path  = "outputs/12_mk_tau_qmax.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_deficit_trend_map(
  "duracion_maxima_evento_dias", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/13_mk_tau_es_duration.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_deficit_trend_map(
  "frecuencia_eventos", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/14_mk_tau_es_frecuency.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

plot_population_trend_map(
  "pop_density", "mk_tau",
  legend_label = "Pendiente tau",
  colors       = colorbar2,
  limits       = c(-1, 1),
  output_path  = "outputs/15_mk_tau_pop_density.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

# Synthetic pressure index (requires 4.1_pressure_index.R to have been run).
# It is already oriented so that +1 = highest pressure, hence colorbar2.
# `limits` is intentionally left NULL: averaging 15 signals pulls the index
# towards 0, so forcing c(-1, 1) would wash the map out to white.
plot_pressure_index_map(
  legend_label = "Indice de presion",
  colors       = colorbar3,
  output_path  = "outputs/16_indice_presion.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)

# The same function maps any theme score or individual sub-indicator, e.g.:
# plot_pressure_index_map("tema_lulc", legend_label = "Presion - uso del suelo",
#                         colors = colorbar2)
# plot_pressure_index_map("n_subind_validos", legend_label = "Sub-indicadores validos",
#                         colors = colorbar2, center_zero = FALSE)

# Version 2 of the index (step 4.2): plain mean of the 15 sub-indicators,
# each weighing 1/15, instead of the 20%-per-theme weighting of step 4.1.
# Same function, just pointed at the other geopackage via `index_path`.
plot_pressure_index_map(
  legend_label = "Indice de presion (v2)",
  colors       = colorbar3,
  index_path   = pressure_index_v2_path,
  output_path  = "outputs/17_indice_presion_v2.png",
  width_cm     = 18,
  height_cm    = 14,
  padding      = 15000
)