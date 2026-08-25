# get.Monin.Obukhov function

`get.Monin.Obukhov` computes the Monin-Obukhov length, a stability
parameter used to characterize atmospheric stability and to correct the
aerodynamic resistances above the canopy for non-neutral conditions.

## Usage

``` r
get.Monin.Obukhov(data.meteo, H)
```

## Arguments

- data.meteo:

  list or data.frame. Meteorological forcing data; must contain `ustar`
  (friction velocity, m/s) and `Ta` (air temperature, deg C).

- H:

  numeric vector. Sensible heat flux from the canopy/soil to the
  atmosphere (W m^-2).

## Value

numeric vector. Monin-Obukhov length `L` (m). Values that are `NA` (e.g.
when `H` is zero) are set to `-1e6`, representing near-neutral
stability.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
L <- get.Monin.Obukhov(data.meteo, H)
} # }
```
