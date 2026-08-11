# ==============================================================================
# Step 1.2 - Land cover change rates by micro-watershed
#
# Takes the long-format composition table produced in step 1.1 and computes,
# for each micro-watershed and each functional group (A to I):
#
#   (a) Full-series linear trend (OLS slope, in percentage points per year)
#   (b) Non-parametric trend: Mann-Kendall test and Sen's slope estimator
#   (c) Change rates for the standard analysis periods: the general period
#       1985-2024, plus 1985-1995, 1995-2005, 2005-2015, 2015-2024
#
# Outputs:
#   trends_by_microbasin.csv        -> one row per comid x group
#   period_rates_by_microbasin.csv  -> one row per comid x group x period
#
# Output columns
# ------------------------------------------------------------------------------
# 02_lulc_trends_by_basin.csv (one row per comid x functional group, using the
# full 1985-2024 series):
#   comid            - Microbasin identifier (linkno).
#   group_code       - Short code for the group: A-I for the nine functional
#                       groups (see group_lookup), or NAT/ANT for the optional
#                       aggregates (total natural cover / total anthropic use,
#                       only present if include_aggregates = TRUE).
#   group            - Full group name (or "Cobertura natural" / "Uso
#                       antropico" for the aggregates).
#   pct_mean         - Mean coverage of the group over 1985-2024, in % of the
#                       microbasin area.
#   ols_slope        - OLS linear trend of pct vs. year, in percentage points
#                       per year (pp/yr). 0 for constant series (no coverage
#                       change across the whole period).
#   ols_pvalue       - p-value of the OLS slope; NA for constant series.
#   mk_tau           - Mann-Kendall tau statistic (non-parametric trend
#                       strength/direction, -1 to 1); NA for constant series.
#   mk_pvalue        - p-value of the Mann-Kendall test; NA for constant
#                       series.
#   sen_slope        - Sen's slope estimator (robust, non-parametric trend),
#                       in percentage points per year; 0 for constant series.
#   trend_direction  - "incremento" / "perdida" when the Mann-Kendall test is
#                       significant (mk_pvalue < 0.05) and sen_slope is
#                       positive/negative respectively; "sin tendencia"
#                       otherwise (including constant series, where mk_pvalue
#                       is NA).
#
# 03_lulc_period_trends_by_basin.csv (one row per comid x functional group x
# analysis period, from the `periods` table: 1985-2024, 1985-1995, 1995-2005,
# 2005-2015, 2015-2024):
#   comid            - Microbasin identifier (linkno).
#   group_code       - Same coding as above (A-I, plus NAT/ANT if enabled).
#   group            - Full group name.
#   period           - Analysis period label (e.g. "1995-2005").
#   pct_start        - Coverage (%) in the period's start year.
#   pct_end          - Coverage (%) in the period's end year.
#   area_start_km2   - Area (km2) covered by the group in the start year.
#   area_end_km2     - Area (km2) covered by the group in the end year.
#   abs_rate_pp_yr   - Absolute change rate within the period: OLS slope of
#                       pct vs. year, in percentage points/year. 0 if pct is
#                       constant across the period's years.
#   rel_rate_pct_yr  - Relative (compound annual) change rate, per Puyravaud
#                       (2003): log(area_end/area_start) / (n years) * 100,
#                       in %/year. Comparable across groups regardless of
#                       their size. NA if either area_start_km2 or
#                       area_end_km2 is 0 (log undefined / rate not
#                       meaningful for a group that starts or ends absent).
# ==============================================================================

library(dplyr)
library(tidyr)
library(trend)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

input_path          <- "workspace/01_lulc_by_basin.csv"
trends_output       <- "workspace/02_lulc_trends_by_basin.csv"
period_rates_output <- "workspace/03_lulc_period_trends_by_basin.csv"

# Set to FALSE to report only the nine functional groups
include_aggregates <- TRUE

# Functional groups, in reporting order. Names must match those written by
# step 1.1; group_code is added for traceability.
group_lookup <- tribble(
  ~group,                                  ~group_code, ~is_natural,
  "Bosque natural",                        "A",         TRUE,
  "Vegetacion altoandina",                 "B",         TRUE,
  "Otra formacion natural",                "C",         TRUE,
  "Agricultura, silvicultura y ganaderia", "D",         FALSE,
  "Infraestructura antropica",             "E",         FALSE,
  "Mineria",                               "F",         FALSE,
  "Area no vegetada natural",              "G",         TRUE,
  "Cuerpos de agua",                       "H",         TRUE,
  "Acuicultura",                           "I",         FALSE
)

# Standard analysis periods: the general period plus four sub-periods used
# to detect acceleration or deceleration of processes.
periods <- tribble(
  ~period,      ~year_start, ~year_end,
  "1985-2024",  1985,        2024,
  "1985-1995",  1985,        1995,
  "1995-2005",  1995,        2005,
  "2005-2015",  2005,        2015,
  "2015-2024",  2015,        2024
)

# ------------------------------------------------------------------------------
# 2. Load composition table and assemble the series to analyse
# ------------------------------------------------------------------------------

composition <- read.csv(input_path) %>%
  left_join(group_lookup, by = "group")

# Series for the nine functional groups
series <- composition %>%
  select(comid, group_code, group, year, area_km2, basin_area_km2, pct)

# Optional aggregates: total natural cover and total anthropic use
if (include_aggregates) {

  aggregated <- composition %>%
    mutate(
      group_code = if_else(is_natural, "NAT", "ANT"),
      group      = if_else(is_natural, "Cobertura natural", "Uso antropico")
    ) %>%
    group_by(comid, group_code, group, year, basin_area_km2) %>%
    summarise(
      area_km2 = sum(area_km2),
      pct      = sum(pct),
      .groups  = "drop"
    )

  series <- bind_rows(series, aggregated)
}

series <- arrange(series, comid, group_code, year)

# ------------------------------------------------------------------------------
# 3. Full-series trends: OLS, Mann-Kendall and Sen's slope
# ------------------------------------------------------------------------------

trend_stats <- function(df) {

  # Constant series (e.g. a group absent throughout) yield no meaningful trend
  if (sd(df$pct) == 0) {
    return(tibble(
      pct_mean = mean(df$pct),
      ols_slope = 0, ols_pvalue = NA_real_,
      mk_tau = NA_real_, mk_pvalue = NA_real_, sen_slope = 0
    ))
  }

  fit <- lm(pct ~ year, data = df)
  mk  <- mk.test(df$pct)
  sen <- sens.slope(df$pct)

  tibble(
    pct_mean   = mean(df$pct),
    ols_slope  = coef(fit)[["year"]],
    ols_pvalue = summary(fit)$coefficients["year", "Pr(>|t|)"],
    mk_tau     = as.numeric(mk$estimates["tau"]),
    mk_pvalue  = mk$p.value,
    sen_slope  = as.numeric(sen$estimates)
  )
}

trends <- series %>%
  group_by(comid, group_code, group) %>%
  group_modify(~ trend_stats(.x)) %>%
  ungroup() %>%
  mutate(
    # Direction of change, reported only where Mann-Kendall is significant
    trend_direction = case_when(
      is.na(mk_pvalue) | mk_pvalue >= 0.05 ~ "sin tendencia",
      sen_slope > 0                        ~ "incremento",
      TRUE                                 ~ "perdida"
    )
  ) %>%
  arrange(comid, group_code)

write.csv(trends, trends_output, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Period-based change rates
# ------------------------------------------------------------------------------
# Two complementary rates are reported:
#   abs_rate_pp_yr  -> OLS slope within the period, in percentage points/year
#   rel_rate_pct_yr -> compound annual rate (Puyravaud 2003), relative to the
#                      initial area; comparable across groups of different size

period_rate <- function(df, y_start, y_end) {

  window <- df %>% filter(year >= y_start, year <= y_end)

  area_start <- window$area_km2[window$year == y_start]
  area_end   <- window$area_km2[window$year == y_end]

  tibble(
    pct_start       = window$pct[window$year == y_start],
    pct_end         = window$pct[window$year == y_end],
    area_start_km2  = area_start,
    area_end_km2    = area_end,
    abs_rate_pp_yr  = if (sd(window$pct) == 0) 0
                      else coef(lm(pct ~ year, data = window))[["year"]],
    rel_rate_pct_yr = if (area_start > 0 && area_end > 0)
                        log(area_end / area_start) / (y_end - y_start) * 100
                      else NA_real_
  )
}

period_rates <- series %>%
  group_by(comid, group_code, group) %>%
  group_modify(~ {
    df <- .x
    periods %>%
      group_by(period) %>%
      group_modify(~ period_rate(df, .x$year_start, .x$year_end)) %>%
      ungroup()
  }) %>%
  ungroup() %>%
  arrange(comid, group_code, period)

write.csv(period_rates, period_rates_output, row.names = FALSE)

message("Done: ", nrow(trends), " trend rows and ",
        nrow(period_rates), " period rows written")
