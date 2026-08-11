# ==============================================================================
# Step 2.2 - Trend analysis of ecological flow deficit indicators
#
# Evaluates whether the occurrence of flows below the Tennant ecological
# threshold is becoming more or less frequent, longer or shorter, in each
# micro-watershed.
#
# Variables analysed:
#   porcentaje_dias_bajo_umbral    share of the year below the threshold (%)
#   dias_bajo_umbral               days per year below the threshold
#   frecuencia_eventos             number of separate deficit events per year
#   duracion_media_eventos_dias    mean length of deficit events (days)
#   duracion_maxima_evento_dias    longest deficit event of the year (days)
#
# As in step 2.1, lag-1 autocorrelation is tested and the Hamed-Rao variance
# correction is applied where it is significant.
#
# Note: these series contain many zero values in micro-watersheds where the
# threshold is rarely crossed. The share of zeros is reported so that weakly
# informative series can be identified.
#
# Output: ecoflow_deficit_trends.csv, one row per comid x variable
# ==============================================================================

library(dplyr)
library(tidyr)
library(trend)
library(modifiedmk)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

input_path  <- "inputs/geoglows/ecological_discharge.csv"
output_path <- "workspace/07_streamflow_deficit_trends.csv"

# Analysis window, aligned with the MapBiomas series
year_min <- 1985
year_max <- 2024

# Minimum daily records for a year to be considered complete
min_records <- 365

deficit_variables <- c(
  "porcentaje_dias_bajo_umbral",
  "dias_bajo_umbral",
  "frecuencia_eventos",
  "duracion_media_eventos_dias",
  "duracion_maxima_evento_dias"
)

# ------------------------------------------------------------------------------
# 2. Load and reshape
# ------------------------------------------------------------------------------

deficit <- read.csv(input_path) %>%
  rename(year = anio) %>%
  filter(year >= year_min, year <= year_max, dias_observados >= min_records) %>%
  pivot_longer(
    cols      = all_of(deficit_variables),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  arrange(comid, variable, year)

# ------------------------------------------------------------------------------
# 3. Trend statistics for one series
# ------------------------------------------------------------------------------

trend_stats <- function(df) {

  x <- df$value[is.finite(df$value)]
  n <- length(x)

  empty <- tibble(
    n_years = n, value_mean = if (n > 0) mean(x) else NA_real_,
    prop_zeros = if (n > 0) mean(x == 0) else NA_real_,
    lag1_acf = NA_real_, autocorrelated = NA,
    mk_tau = NA_real_, mk_pvalue = NA_real_, mk_pvalue_hr = NA_real_,
    pvalue_final = NA_real_, sen_slope = NA_real_
  )

  # A constant series (typically all zeros) carries no trend information
  if (n < 10 || sd(x) == 0) return(empty)

  # Lag-1 autocorrelation, significant beyond the 95% white-noise bounds
  lag1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
  autocorrelated <- abs(lag1) > 1.96 / sqrt(n)

  mk  <- mk.test(x)
  sen <- sens.slope(x)
  hr  <- mmkh(x)   # Hamed-Rao (1998) correction for serial correlation

  tibble(
    n_years        = n,
    value_mean     = mean(x),
    prop_zeros     = mean(x == 0),
    lag1_acf       = lag1,
    autocorrelated = autocorrelated,
    mk_tau         = as.numeric(mk$estimates["tau"]),
    mk_pvalue      = mk$p.value,
    mk_pvalue_hr   = as.numeric(hr["new P-value"]),
    pvalue_final   = if (autocorrelated) as.numeric(hr["new P-value"])
                     else mk$p.value,
    sen_slope      = as.numeric(sen$estimates)   # units of the variable per year
  )
}

# ------------------------------------------------------------------------------
# 4. Apply to every micro-watershed and variable
# ------------------------------------------------------------------------------

trends <- deficit %>%
  group_by(comid, variable) %>%
  group_modify(~ trend_stats(.x)) %>%
  ungroup() %>%
  mutate(
    trend_direction = case_when(
      is.na(pvalue_final) | pvalue_final >= 0.05 ~ "sin tendencia",
      sen_slope > 0                              ~ "incremento",
      TRUE                                       ~ "descenso"
    ),
    # For deficit variables an increase means worsening conditions
    interpretacion = case_when(
      trend_direction == "sin tendencia" ~ "sin cambio detectable",
      trend_direction == "incremento"    ~ "mayor estres hidrico",
      TRUE                               ~ "menor estres hidrico"
    )
  ) %>%
  arrange(comid, variable)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(trends, output_path, row.names = FALSE)

message("Done: ", nrow(trends), " rows written to ", output_path)
