# ==============================================================================
# Step 3.2 - Population trend analysis by micro-watershed
#
# Applies the same non-parametric trend analysis used in steps 1.2/1.3/2.1
# (Mann-Kendall test + Sen's slope estimator) to the population time series
# produced in step 3.1, for two variables:
#   population   - total inhabitants in the microbasin
#   pop_density  - inhabitants per km2
#
# The GHS-POP series only has 9 time points, five years apart (1985-2025).
# That is too short to reliably estimate lag-1 autocorrelation, so unlike
# step 2.1's daily discharge series, no Hamed-Rao correction is applied here;
# a classic Mann-Kendall test / Sen's slope is used instead, as in steps 1.2
# and 1.3.
#
# Output (workspace/09_population_trends_by_basin.csv), one row per comid x
# variable:
#   comid            - Microbasin identifier (linkno).
#   variable         - "population" (total inhabitants) or "pop_density"
#                       (inhabitants / km2).
#   n_years          - Number of finite values in the series (out of 9).
#   value_mean       - Mean of the series over 1985-2025, in the variable's
#                       own unit.
#   ols_slope        - OLS linear trend of value vs. year (unit/year). 0 for
#                       constant series (e.g. a microbasin with population 0
#                       throughout).
#   ols_pvalue       - p-value of the OLS slope; NA for constant series.
#   mk_tau           - Mann-Kendall tau statistic (non-parametric trend
#                       strength/direction, -1 to 1); NA for constant or
#                       too-short series.
#   mk_pvalue        - p-value of the Mann-Kendall test; NA under the same
#                       conditions as mk_tau.
#   sen_slope        - Sen's slope estimator (robust, non-parametric trend),
#                       in the variable's unit per year; 0 for constant
#                       series.
#   sen_slope_rel    - sen_slope as a percentage of value_mean (%/year),
#                       comparable across microbasins regardless of their
#                       population size. NA if value_mean is 0.
#   trend_direction  - "aumento" / "disminucion" when the Mann-Kendall test
#                       is significant (mk_pvalue < 0.05) and sen_slope is
#                       positive/negative respectively; "sin tendencia"
#                       otherwise (including constant series, where
#                       mk_pvalue is NA).
# ==============================================================================

library(dplyr)
library(tidyr)
library(trend)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

input_path  <- "workspace/08_population_by_basin.csv"
output_path <- "workspace/09_population_trends_by_basin.csv"

# Minimum number of time points required to attempt a trend test
min_years <- 5

# ------------------------------------------------------------------------------
# 2. Load and reshape population data
# ------------------------------------------------------------------------------

population <- read.csv(input_path) %>%
  pivot_longer(
    cols      = c(population, pop_density),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  arrange(comid, variable, year)

# ------------------------------------------------------------------------------
# 3. Trend statistics for one series
# ------------------------------------------------------------------------------

trend_stats <- function(df) {

  finite <- is.finite(df$value)
  value  <- df$value[finite]
  n      <- length(value)

  # Constant series (e.g. a microbasin with no population throughout) or too
  # few points yield no meaningful trend
  if (n < min_years || sd(value) == 0) {
    return(tibble(
      n_years       = n,
      value_mean    = if (n > 0) mean(value) else NA_real_,
      ols_slope     = 0,
      ols_pvalue    = NA_real_,
      mk_tau        = NA_real_,
      mk_pvalue     = NA_real_,
      sen_slope     = 0,
      sen_slope_rel = NA_real_
    ))
  }

  fit <- lm(value ~ year, data = df[finite, ])
  mk  <- mk.test(value)
  sen <- sens.slope(value)

  value_mean <- mean(value)
  sen_slope  <- as.numeric(sen$estimates)

  tibble(
    n_years       = n,
    value_mean    = value_mean,
    ols_slope     = coef(fit)[["year"]],
    ols_pvalue    = summary(fit)$coefficients["year", "Pr(>|t|)"],
    mk_tau        = as.numeric(mk$estimates["tau"]),
    mk_pvalue     = mk$p.value,
    sen_slope     = sen_slope,
    sen_slope_rel = if (value_mean > 0) sen_slope / value_mean * 100
                    else NA_real_
  )
}

# ------------------------------------------------------------------------------
# 4. Apply to every micro-watershed and variable
# ------------------------------------------------------------------------------

trends <- population %>%
  group_by(comid, variable) %>%
  group_modify(~ trend_stats(.x)) %>%
  ungroup() %>%
  mutate(
    trend_direction = case_when(
      is.na(mk_pvalue) | mk_pvalue >= 0.05 ~ "sin tendencia",
      sen_slope > 0                        ~ "aumento",
      TRUE                                 ~ "disminucion"
    )
  ) %>%
  arrange(comid, variable)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(trends, output_path, row.names = FALSE)

message("Done: ", nrow(trends), " rows written to ", output_path)
