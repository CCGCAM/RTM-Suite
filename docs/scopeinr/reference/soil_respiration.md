# soil respiration

soil respiration

## Usage

``` r
soil_respiration(Ts)
```

## Arguments

- Ts:

  numeric. Soil temperature.

## Value

R: soil respiration rate, in umol m-2 s-1. This R port always returns 0
(soil respiration is not simulated) - kept for interface compatibility
with the original SCOPE model.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
soil_respiration(20)
#> [1] 0
```
