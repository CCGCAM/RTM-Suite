# subfunction phs for stability correction (eg. Paulson, 1970)

subfunction phs for stability correction (eg. Paulson, 1970)

## Usage

``` r
get.phstar(z, zR, d, L, stable, unstable, x)
```

## Arguments

- z:

  numeric. Height above the surface, in meters.

- zR:

  numeric. Reference height, in meters.

- d:

  numeric. Zero-plane displacement height, in meters.

- L:

  numeric. Monin-Obukhov length, in meters (atmospheric stability
  scale).

- stable:

  logical. TRUE if the atmosphere is in a stable stratification regime.

- unstable:

  logical. TRUE if the atmosphere is in an unstable stratification
  regime.

- x:

  numeric. Stability correction intermediate variable, typically (1 -
  16\*z/L)^0.25.

## Value

phs vector

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.phstar(z = 2, zR = 10, d = 0.7, L = -50, stable = FALSE, unstable = TRUE, x = 1.2)
#> [1] 0.02520712
```
