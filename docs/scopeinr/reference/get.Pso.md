# get.Pso function

`get.Pso` calculates the bi-directional gap probability `Pso`, i.e. the
joint probability of a sunlit and simultaneously viewed leaf (or soil
background) at a given normalized canopy depth, accounting for the
hot-spot effect between the solar and viewing directions.

## Usage

``` r
get.Pso(K, k, LAI, q, dso, xl)
```

## Arguments

- K:

  numeric value. Canopy extinction coefficient in the direction of the
  observer.

- k:

  numeric value. Canopy extinction coefficient in the direction of the
  sun.

- LAI:

  numeric value. Total (one-sided) leaf area index of the canopy (m2
  m-2).

- q:

  numeric value. Hot-spot size parameter (leaf width to canopy height
  ratio, dimensionless).

- dso:

  numeric value. Normalized distance between the sun and observer beams
  at the canopy level, derived from the solar/viewing zenith and
  relative azimuth angles (dimensionless).

- xl:

  numeric value. Normalized cumulative canopy depth at which `Pso` is
  evaluated (dimensionless, 0 at the top of the canopy to -1 at the
  bottom, expressed as a fraction of `LAI`).

## Value

numeric value. Bi-directional gap probability `Pso` at depth `xl`
(dimensionless, 0-1).

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
pso <- get.Pso(K, k, LAI, q, dso, xl)
} # }
```
