# subfunction ph for stability correction (eg. Paulson, 1970)

subfunction ph for stability correction (eg. Paulson, 1970)

## Usage

``` r
get.psih(z, L, unstable, stable, x)
```

## Arguments

- z:

  numeric. Height above the surface, in meters.

- L:

  numeric. Monin-Obukhov length, in meters (atmospheric stability
  scale).

- unstable:

  logical. TRUE if the atmosphere is in an unstable stratification
  regime.

- stable:

  logical. TRUE if the atmosphere is in a stable stratification regime.

- x:

  numeric. Stability correction intermediate variable, typically (1 -
  16\*z/L)^0.25.

## Value

ph vector

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.psih(z = 2, L = -50, unstable = TRUE, stable = FALSE, x = 1.2)
#> [1] 0.3977017
```
