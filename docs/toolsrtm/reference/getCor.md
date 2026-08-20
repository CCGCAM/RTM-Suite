# Generate multiple mutually-correlated variables at a target correlation level

Builds `n_inputs` variables, each `nLUT` values long, all
pairwise-correlated at (approximately) `rho`, then rescales each to its
own `[MinRange, MaxRange]` bound. Useful for LUTs where several traits
should co-vary (e.g. LAI and canopy height) rather than being drawn
independently.

## Usage

``` r
getCor(
  n_inputs = NULL,
  nLUT = 100,
  distribution = "Uniform",
  setseed = 123,
  rho = NULL,
  Varnames = NULL,
  MinRange = NULL,
  MaxRange = NULL
)
```

## Arguments

- n_inputs:

  integer. How many variables to generate.

- nLUT:

  integer. How many values per variable (LUT rows). Default 100.

- distribution:

  character. `"Uniform"` or `"Normal"`. Default `"Uniform"`.

- setseed:

  integer. Random seed. Default 123.

- rho:

  numeric (-1 to 1). Target pairwise correlation between every pair of
  variables.

- Varnames:

  character vector, length `n_inputs`. Column names for the output; if
  `NULL`, columns are named `Var_1`, `Var_2`, ...

- MinRange:

  numeric vector, length `n_inputs`. Minimum of each variable's output
  range.

- MaxRange:

  numeric vector, length `n_inputs`. Maximum of each variable's output
  range.

## Value

A list: `LUT` (data.frame, `nLUT` rows x `n_inputs` columns) and
`Covarianza` (the realized correlation matrix of `LUT`, for checking how
close the sample came to the requested `rho`).

## Examples

``` r
out <- getCor(n_inputs = 2, nLUT = 200, distribution = "Uniform", rho = 0.7,
               Varnames = c("LAI", "Height"), MinRange = c(0.5, 2), MaxRange = c(7, 30))
#> Generating a Uniform distribution for all correlated inputs ...
cor(out$LUT$LAI, out$LUT$Height)  # close to 0.7
#> [1] 0.7273228
```
