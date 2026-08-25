# get.gsFun `get.gsFun` get stomatal conductance

get.gsFun `get.gsFun` get stomatal conductance

## Usage

``` r
get.gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
```

## Arguments

- Cs:

  numeric. CO2 concentration at the leaf surface.

- RH:

  numeric. Relative humidity (0-1).

- A:

  numeric. Net assimilation rate, in the same CO2 units as Cs, per m2
  per s.

- BallBerrySlope:

  numeric. Slope parameter of the Ball-Berry stomatal conductance model.

- BallBerry0:

  numeric. Intercept (minimum conductance) parameter of the Ball-Berry
  model.

## Value

gs: stomatal conductance.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.gsFun(Cs = 400, RH = 0.7, A = 15, BallBerrySlope = 9, BallBerry0 = 0.01)
#> [1] 0.24625
```
