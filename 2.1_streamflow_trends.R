# ==============================================================================
# Step 2.1 - Trend analysis of annual discharge series
#
# Applies non-parametric trend tests to the annual discharge series of each
# micro-watershed, for three variables:
#   qavg - mean annual discharge
#   qmin - annual minimum discharge (water security relevance)
#   qmax - annual maximum discharge (flood regime)
#
# Hydrological series are commonly autocorrelated, which inflates the
# significance of the classic Mann-Kendall test. Lag-1 autocorrelation is
# therefore tested first, and the Hamed-Rao variance correction is applied
# where it is significant.
#
# Output: discharge_trends_by_microbasin.csv, one row per comid x variable
# ==============================================================================

library(dplyr)
library(tidyr)
library(trend)
library(modifiedmk)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

input_path  <- "inputs/geoglows/anual_discharge.csv"
output_path <- "workspace/06_streamflow_trends_by_basin.csv"

# Analysis window, aligned with the MapBiomas series
year_min <- 1985
year_max <- 2024

# Minimum number of daily records for a year to be considered complete
min_records <- 365

# ------------------------------------------------------------------------------
# 2. Load and reshape discharge data
# ------------------------------------------------------------------------------

discharge <- read.csv(input_path) %>%
  filter(year >= year_min, year <= year_max, records >= min_records) %>%
  pivot_longer(
    cols      = c(qavg, qmin, qmax),
    names_to  = "variable",
    values_to = "q"
  ) %>%
  arrange(comid, variable, year)

# ------------------------------------------------------------------------------
# 3. Trend statistics for one series
# ------------------------------------------------------------------------------

trend_stats <- function(df) {

  q <- df$q[is.finite(df$q)]
  n <- length(q)

  if (n < 10 || sd(q) == 0) {
    return(tibble(
      n_years = n, q_mean = if (n > 0) mean(q) else NA_real_,
      lag1_acf = NA_real_, autocorrelated = NA,
      mk_tau = NA_real_, mk_pvalue = NA_real_, mk_pvalue_hr = NA_real_,
      pvalue_final = NA_real_, sen_slope = NA_real_, sen_slope_rel = NA_real_
    ))
  }

  # Lag-1 autocorrelation, significant beyond the 95% white-noise bounds
  lag1 <- acf(q, lag.max = 1, plot = FALSE)$acf[2]
  autocorrelated <- abs(lag1) > 1.96 / sqrt(n)

  mk  <- mk.test(q)
  sen <- sens.slope(q)

  # Hamed-Rao (1998) variance correction for serially correlated data
  hr <- mmkh(q)

  q_mean    <- mean(q)
  sen_slope <- as.numeric(sen$estimates)

  tibble(
    n_years        = n,
    q_mean         = q_mean,
    lag1_acf       = lag1,
    autocorrelated = autocorrelated,
    mk_tau         = as.numeric(mk$estimates["tau"]),
    mk_pvalue      = mk$p.value,
    mk_pvalue_hr   = as.numeric(hr["new P-value"]),
    pvalue_final   = if (autocorrelated) as.numeric(hr["new P-value"])
                     else mk$p.value,
    sen_slope      = sen_slope,                        # m3/s per year
    sen_slope_rel  = sen_slope / q_mean * 100          # % of mean per year
  )
}

# ------------------------------------------------------------------------------
# 4. Apply to every micro-watershed and variable
# ------------------------------------------------------------------------------

trends <- discharge %>%
  group_by(comid, variable) %>%
  group_modify(~ trend_stats(.x)) %>%
  ungroup() %>%
  mutate(
    # Direction based on the p-value appropriate to each series
    trend_direction = case_when(
      is.na(pvalue_final) | pvalue_final >= 0.05 ~ "sin tendencia",
      sen_slope > 0                              ~ "incremento",
      TRUE                                       ~ "descenso"
    )
  ) %>%
  arrange(comid, variable)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(trends, output_path, row.names = FALSE)

message("Done: ", nrow(trends), " rows written to ", output_path)
