# Basin Prioritization — Pastaza River Basin

An R analysis pipeline that ranks the ~1,747 microbasins of the Pastaza river basin (Ecuador) by
pressure on aquatic ecosystems. Each stage derives a per-basin time series from raster or tabular
sources, fits a Mann-Kendall trend to it, and the final stages combine all trends into a single
synthetic pressure index and a set of choropleth maps.

## Contents

- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Input data](#input-data)
- [Running the pipeline](#running-the-pipeline)
- [Methodology](#methodology)
- [Outputs](#outputs)
- [Known limitations](#known-limitations)

## Requirements

- R ≥ 4.5 (developed and tested on R 4.5.1)
- R packages: `sf`, `terra`, `exactextractr`, `landscapemetrics`, `dplyr`, `tidyr`, `tibble`,
  `trend`, `modifiedmk`, `ggplot2`, `scales`

```r
install.packages(c(
  "sf", "terra", "exactextractr", "landscapemetrics",
  "dplyr", "tidyr", "tibble", "trend", "modifiedmk",
  "ggplot2", "scales"
))
```

No package manifest (`renv.lock`, `DESCRIPTION`) is tracked yet; the list above is exhaustive as of
this pipeline's current scripts.

## Repository structure

```
basin_prioritization/
├── 1.1_lulc_composition.R        Land cover composition by microbasin (MapBiomas)
├── 1.2_lulc_change_rates.R       Land cover trend analysis (Mann-Kendall + Sen's slope)
├── 1.3_lulc_fragmentation.R      Landscape fragmentation metrics + trends
├── 2.1_streamflow_trends.R       Annual discharge trend analysis
├── 2.2_streamflow_ecological.R   Ecological flow deficit trend analysis
├── 3.1_population_density.R      Population density by microbasin (GHS-POP)
├── 3.2_populations_trends.R      Population trend analysis
├── 4.1_pressure_index.R          Synthetic pressure index, v1 (20% per theme)
├── 4.2_pressure_index_v2.R       Synthetic pressure index, v2 (equal weight per sub-indicator)
├── 5.1_plots.R                   All choropleth maps (trend indicators + pressure index)
├── inputs/                       Read-only source data (see below)
│   ├── map_biomas/               Annual land-cover rasters, 1985–2024
│   ├── ghs_pop/                  Population rasters, 1985–2025 (5-year steps)
│   ├── geoglows/                 Discharge and ecological-flow time series (CSV)
│   └── geographical_data/        Microbasin, study-area and province polygons (GeoPackage)
├── workspace/                    Numbered intermediate tables/geopackages (01–11)
├── outputs/                      Numbered PNG figures (01–17)
└── CLAUDE.md                     Internal guidance for AI coding agents working in this repo
```

Each script is numbered `<stage>.<step>_<name>.R` and opens with a comment block documenting every
output column — that header is the authoritative reference for a table's schema; this README
summarizes it.

## Input data

Not tracked in this repository (`inputs/` is ~5.5 GB). You must source these yourself and place
them under `inputs/` following the layout above:

| Path | Content | Used by |
|---|---|---|
| `inputs/map_biomas/<year>.tif` | Annual MapBiomas land-cover classification, EPSG:4326, one GeoTIFF per year 1985–2024 (40 files) | 1.1, 1.3 |
| `inputs/ghs_pop/<year>.tif` | GHS-POP population-count rasters, one per 5-year step 1985–2025 | 3.1 |
| `inputs/geoglows/anual_discharge.csv` | Annual `qavg`/`qmin`/`qmax` discharge per `comid` (GEOGloWS) | 2.1 |
| `inputs/geoglows/ecological_discharge.csv` | Annual Tennant ecological-flow deficit metrics per `comid` | 2.2 |
| `inputs/geographical_data/geoglows_basins.gpkg` | The ~1,747 study microbasins, keyed by `linkno` | all scripts (as `comid`) |
| `inputs/geographical_data/pastaza.gpkg` | Study-area outline, used as a map reference boundary | 5.1 |
| `inputs/geographical_data/provinces.gpkg` | Province boundaries, used as a map reference layer | 5.1 |

`inputs/geoglows/monthly_discharge.csv` and `inputs/geographical_data/pastaza_basin.gpkg` are also
present but are not read by any script in this pipeline.

## Running the pipeline

Run from the repository root — every path in the scripts is relative to it.

```bash
Rscript 1.1_lulc_composition.R      # MapBiomas rasters      -> workspace/01_lulc_by_basin.csv
Rscript 1.2_lulc_change_rates.R     # 01                     -> workspace/02, 03
Rscript 1.3_lulc_fragmentation.R    # MapBiomas rasters      -> workspace/04, 05
Rscript 2.1_streamflow_trends.R     # geoglows/anual_discharge.csv       -> workspace/06
Rscript 2.2_streamflow_ecological.R # geoglows/ecological_discharge.csv -> workspace/07
Rscript 3.1_population_density.R    # GHS-POP rasters        -> workspace/08
Rscript 3.2_populations_trends.R    # 08                     -> workspace/09
Rscript 4.1_pressure_index.R        # 02, 05, 06, 07, 09     -> workspace/10 (.gpkg)
Rscript 4.2_pressure_index_v2.R     # 02, 05, 06, 07, 09     -> workspace/11 (.gpkg)
Rscript 5.1_plots.R                 # everything above       -> outputs/*.png
```

Stages `1.x`, `2.x` and `3.x` are independent of each other and can run in any order (or in
parallel); only `4.x` and `5.1` need the earlier stages' outputs already in `workspace/`. The raster
stages (`1.1`, `1.3`, `3.1`) are the slowest, since they iterate over the 40 MapBiomas and 9 GHS-POP
rasters. `5.1_plots.R` also takes a few minutes, as it regenerates all 17 figures every time it
runs — it is not meant to be `source()`d just to reuse its plotting functions.

`comid` is the join key across every table in this pipeline, always derived from the basins layer's
`linkno` field.

## Methodology

### Trend detection

All trend tables (`workspace/02`, `05`, `06`, `07`, `09`) are computed the same way: for each
microbasin (and, where relevant, each functional group / metric / variable), a non-parametric
**Mann-Kendall test** is applied to the annual time series, together with a **Sen's slope**
estimator (via the `trend` package). Every trend table exposes an `mk_tau` column bounded to
[-1, 1] — the common currency the pressure index is built from — plus the corresponding p-value and
a Spanish `trend_direction` label.

The two discharge-based stages (`2.1`, `2.2`) additionally test lag-1 autocorrelation and, where it
is significant, apply the **Hamed-Rao (1998) variance correction** (`modifiedmk::mmkh`) before
judging significance. The shorter land-cover, fragmentation and population series use the plain
Mann-Kendall p-value instead.

### Synthetic pressure index

`4.1_pressure_index.R` and `4.2_pressure_index_v2.R` combine 15 `mk_tau` signals — 5 land-cover
groups, 4 fragmentation metrics, 3 discharge variables, 2 deficit variables and 1 population
metric — into a single index per microbasin, bounded to [-1, 1], where **+1 means the highest
pressure**.

1. **Sign orientation.** Each sub-indicator is reoriented so a positive value always means more
   pressure: `ori = direction × tau`, where `direction` is `-1` for indicators where a *falling*
   tau signals pressure (e.g. forest loss, discharge decline) and `+1` where a *rising* tau does
   (e.g. mining expansion, fragmentation). This mirrors the two color scales used in `5.1_plots.R`
   (`colorbar` = red at −1, `colorbar2` = red at +1).
2. **Missing trends.** `mk_tau = NA` (a constant or too-short series) is treated as "no trend" and
   contributes 0 pressure, keeping the denominator fixed across basins. A handful of numerically
   invalid `tau` values from `trend::mk.test` (`|tau| > 1` or non-finite, mostly from the
   fragmentation `pd` metric) are demoted to `NA` before aggregation.
3. **Aggregation** differs between the two versions:
   - **v1** (`4.1`) averages within each of the 5 themes first, then averages the 5 theme scores —
     every theme weighs 20%, regardless of how many sub-indicators it contributes.
   - **v2** (`4.2`) is a plain average of all 15 sub-indicators, each weighing 1/15 (implied theme
     weights: land cover 33%, fragmentation 27%, discharge 20%, deficit 13%, population 7%).

Both outputs keep the raw `tau_*` (with its original `NA`) alongside the oriented `ori_*`, the
theme scores `tema_*`, the final `indice_presion` and `n_subind_validos` (how many of the 15
sub-indicators had real, non-`NA` data for that basin), so the aggregation is fully auditable from
the geopackage alone.

## Outputs

### `workspace/` — intermediate tables and geopackages

| File | Produced by | Content |
|---|---|---|
| `01_lulc_by_basin.csv` | 1.1 | Land-cover composition (% by functional group), one row per basin × group × year |
| `02_lulc_trends_by_basin.csv` | 1.2 | Mann-Kendall trend of land-cover % over time, per basin × group |
| `03_lulc_period_trends_by_basin.csv` | 1.2 | Change rates over fixed sub-periods (1985–1995, …, 2015–2024) |
| `04_fragmentation_by_basin.csv` | 1.3 | Landscape fragmentation metrics (`landscapemetrics`), per basin × year |
| `05_fragmentation_trends_by_basin.csv` | 1.3 | Mann-Kendall trend of each fragmentation metric |
| `06_streamflow_trends_by_basin.csv` | 2.1 | Mann-Kendall trend of annual `qavg`/`qmin`/`qmax` discharge |
| `07_streamflow_deficit_trends.csv` | 2.2 | Mann-Kendall trend of ecological-flow deficit indicators |
| `08_population_by_basin.csv` | 3.1 | Population count and density per basin × year |
| `09_population_trends_by_basin.csv` | 3.2 | Mann-Kendall trend of population and population density |
| `10_pressure_index_by_basin.gpkg` | 4.1 | Pressure index v1 (20% per theme), with all raw/oriented sub-indicators |
| `11_pressure_index_v2_by_basin.gpkg` | 4.2 | Pressure index v2 (equal weight per sub-indicator) |

Full column-by-column definitions live in the header comment of the script that produces each
table.

### `outputs/` — figures

`5.1_plots.R` produces 17 choropleth maps of microbasins, all sharing the same visual language: fill
and border colored by the mapped variable, the study-area outline and province boundaries overlaid
as reference layers, and a fixed map extent.

- `01`–`05`: Mann-Kendall tau of land-cover trends (natural forest, high-Andean vegetation,
  agriculture, mining, water bodies)
- `06`–`09`: Mann-Kendall tau of fragmentation metrics (patch count, patch density, patch area,
  cohesion)
- `10`–`12`: Mann-Kendall tau of discharge trends (mean, minimum, maximum flow)
- `13`–`14`: Mann-Kendall tau of ecological-flow deficit trends (event duration, event frequency)
- `15`: Mann-Kendall tau of population density trend
- `16`: Synthetic pressure index, version 1
- `17`: Synthetic pressure index, version 2
