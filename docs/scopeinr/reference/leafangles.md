# Subroutine FluorSail_dladgen (Version 2.3)

`leafangles` computes the leaf inclination distribution function (LIDF)
frequencies for the 13 standard leaf angle classes, from the two Verhoef
ellipsoidal LIDF shape parameters, by numerically inverting the
cumulative distribution via [`dcum`](dcum.md).

## Usage

``` r
leafangles(a, b)
```

## Arguments

- a:

  numeric. LIDFa parameter, controls the average leaf slope (bimodality
  parameter of the ellipsoidal distribution).

- b:

  numeric. LIDFb parameter, controls the distribution's
  bimodality/spread.

## Value

A list with: `lidf` (numeric vector of length 13, the fraction of leaf
area in each leaf angle class) and `litab` (numeric vector of length 13,
the leaf angle class centers, degrees: `5,15,...,89`).

## References

For more information look to page 128 of "theory of radiative transfer
models applied in optical remote sensing of vegetation canopies"

## Author

    Joris Timmermans, Christiaan van der Tol  (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
LeafDistribution <- leafangles(a = -0.35, b = -0.15)
} # }
```
