# get.zo_and_d model `get.zo_and_d` Calculates roughness length for momentum and zero plane displacement from vegetation height and LAI

get.zo_and_d model `get.zo_and_d` Calculates roughness length for
momentum and zero plane displacement from vegetation height and LAI

## Usage

``` r
get.zo_and_d(inputLUT, constants, calc.heat, calc.rss_rbs)
```

## Arguments

- inputLUT:

  list. Canopy structure inputs (vegetation height, LAI, etc.) for the
  roughness-length calculation.

- constants:

  list. Physical constants used in the aerodynamic roughness
  formulation.

- calc.heat:

  logical. Whether to include the heat-flux-related roughness adjustment
  terms.

- calc.rss_rbs:

  logical. Whether to also compute soil/boundary-layer resistance terms
  alongside roughness length and displacement height.

## Value

zo_and_d: a list with roughness length for momentum and zero-plane
displacement height.

## References

Verhoef, McNaughton & Jacobs (1997), HESS 1, 81-91

## Author

    A. Verhoef (Original version in Matlab)

Carlos Camino (Ported version into R)

last updates:

- 17 November 2008

## Examples

``` r
if (FALSE) { # \dontrun{
get.zo_and_d(inputLUT, constants, calc.heat = TRUE, calc.rss_rbs = TRUE)
} # }
```
