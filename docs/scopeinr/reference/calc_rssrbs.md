# Function for calculating the rss and rbs values based on the input parameters (SMC, LAI and rbs)

`calc_rssrbs` calculates the soil surface resistance to evaporation
(`rss`) as an empirical function of soil moisture content, and scales
the boundary-layer resistance of the soil (`rbs`) with leaf area index.

## Usage

``` r
calc_rssrbs(SMC, LAI, rbs)
```

## Arguments

- SMC:

  numeric. Soil moisture content (m^3 m^-3).

- LAI:

  numeric. Leaf area index of the canopy (m^2 m^-2).

- rbs:

  numeric. Boundary-layer resistance of the soil at reference LAI = 3.3
  (s m^-1); overwritten in the returned list with the LAI-scaled value.

## Value

A list with:

- rss:

  numeric. Soil surface resistance to evaporation (s m^-1), computed as
  `11.2 * exp(42 * (0.22 - SMC))`.

- rbs:

  numeric. Soil boundary-layer resistance scaled by LAI (s m^-1),
  computed as `rbs * LAI / 3.3`.

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
rss_rbs <- calc_rssrbs(SMC, LAI, rbs)
} # }
```
