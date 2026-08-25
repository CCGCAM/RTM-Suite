# Simpson-like trapezoidal integration

`Sint` numerically integrates the vector `y` over `x` using the
trapezoidal rule, i.e. the sum of the areas of the trapezoids formed by
consecutive point pairs.

## Usage

``` r
Sint(y, x)
```

## Arguments

- y:

  numeric vector. Values of the integrand, same length as `x`.

- x:

  numeric vector. Values of the integration variable (e.g. wavelength),
  same length as `y`; must be a monotonically increasing series.

## Value

numeric value. The integral of `y` with respect to `x`.

## References

WV Jan. 2013, for SCOPE 1.40

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
int <- Sint(y, x)
} # }
```
