# MD12 algorithm for the computation of fluorescence yield

MD12 algorithm for the computation of fluorescence yield

## Usage

``` r
get.MD12(ps, Ja, Jms, kps, kf, kds, kDs)
```

## Arguments

- ps:

  numeric. Photochemical yield (fraction of absorbed light used in
  photochemistry).

- Ja:

  numeric. Actual electron transport rate.

- Jms:

  numeric. Maximum (light-saturated) electron transport rate.

- kps:

  numeric. Rate constant for photochemistry.

- kf:

  numeric. Rate constant for fluorescence.

- kds:

  numeric. Rate constant for (light-independent) thermal dissipation.

- kDs:

  numeric. Rate constant for (light-dependent, sustained) thermal
  dissipation.

## Value

fs: PSII fluorescence yield.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
get.MD12(ps = 0.3, Ja = 40, Jms = 120, kps = 4, kf = 0.05, kds = 0.95, kDs = 0)
#> [1] 0.005625
```
