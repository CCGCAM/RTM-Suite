# count_k

`count_k` increments a multi-digit "odometer" counter, used to enumerate
all combinations of a set of discretized variables (e.g. when building a
full-factorial input LUT). Starting from digit `id`, every digit that is
already at its maximum value is reset to 1 and carried over to the next
digit (cyclically, wrapping from `nvars` back to 1); the first digit
found that is not at its maximum is incremented by one.

## Usage

``` r
count_k(nvars, v, vmax, id)
```

## Arguments

- nvars:

  integer. Total number of digits (variables) in the counter.

- v:

  integer vector of length `nvars`. Current value of each digit.

- vmax:

  integer vector of length `nvars`. Maximum value allowed for each
  digit.

- id:

  integer. Index of the digit at which to start incrementing.

## Value

integer vector of length `nvars`. Updated digit vector `vnew` after
incrementing.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
vnew <- count_k(nvars, v, vmax, id)
} # }
```
