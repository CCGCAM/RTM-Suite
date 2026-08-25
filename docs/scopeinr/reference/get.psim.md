# subfunction pm for stability correction (eg. Paulson, 1970)

`get.psim` computes the stability correction function for momentum
transfer, used to correct the logarithmic wind profile for non-neutral
atmospheric conditions.

## Usage

``` r
get.psim(z, L, unstable, stable, x)
```

## Arguments

- z:

  numeric. Height above the zero-plane displacement height, in meters.

- L:

  numeric. Monin-Obukhov length, in meters (atmospheric stability
  scale).

- unstable:

  logical. TRUE if the atmosphere is in an unstable stratification
  regime.

- stable:

  logical. TRUE if the atmosphere is in a stable stratification regime.

- x:

  numeric. Stability correction intermediate variable, typically
  `(1 - 16*z/L)^0.25`.

## Value

numeric. Momentum stability correction function `pm`; `0` under neutral
conditions.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.psim(z = 2, L = -50, unstable = TRUE, stable = FALSE, x = 1.2)
#> [1] 0.2081514
```
