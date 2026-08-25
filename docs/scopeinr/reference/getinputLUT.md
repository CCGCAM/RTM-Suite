# Get main input for SCOPE

`getinputLUT` extracts one group of SCOPE model input parameters (a
single row, or a set of rows, of a LUT/input table) into the named list
structure expected by the rest of SCOPEinR (leaf biology, meteo, canopy
structure, geometry, soil, coordinates, time series, or multilayer leaf
biochemistry profile). For `dataset == 'soil'`, it also derives the soil
thermal inertia (`GAM`, via [`Soil_Inertia0`](Soil_Inertia0.md) or
[`Soil_Inertia1`](Soil_Inertia1.md)) and, optionally, the soil/canopy
resistances `rss`/`rbs` (via [`calc_rssrbs`](calc_rssrbs.md)).

## Usage

``` r
getinputLUT(inputLUT, dataset, calc.heat, calc.rss_rbs)
```

## Arguments

- inputLUT:

  data.frame or matrix. LUT/input table with one column per named
  parameter (e.g. `N`, `Cab`, `LAI`, `tts`, `Ta`, ...) and one row per
  simulation.

- dataset:

  character. Which group of parameters to extract; one of `'leafbio'`,
  `'meteo'`, `'canopy'`, `'angles'`, `'soil'`, `'Coordinates'`, `'mly'`
  or `'timeseries'`.

- calc.heat:

  numeric. Soil heat flux method used only when `dataset == 'soil'`: `0`
  estimates `GAM` with [`Soil_Inertia0`](Soil_Inertia0.md) (from `cs`,
  `rhos`, `lambdas`), `1` estimates `GAM` with
  [`Soil_Inertia1`](Soil_Inertia1.md) (from `SMC`), any other value
  (default `2`) leaves `GAM` as `NA` (soil heat flux computed elsewhere
  as `0.35*Rn`).

- calc.rss_rbs:

  numeric. Flag used only when `dataset == 'soil'`: `0` (default) uses
  `rss`/`rbs` directly from `inputLUT`; any other value recomputes them
  from `SMC`/`LAI`/`rbs` via [`calc_rssrbs`](calc_rssrbs.md).

## Value

A named list with the extracted parameters for the requested `dataset`
(e.g. `leafbio`, `meteo`, `canopy`, `angles`, `soil`, `coord`,
`timeseries`, or `mly`).

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
leafbio <- getinputLUT(inputLUT, dataset = 'leafbio')
soil <- getinputLUT(inputLUT, dataset = 'soil', calc.heat = 2, calc.rss_rbs = 0)
} # }
```
