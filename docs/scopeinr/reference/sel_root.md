# quadratic formula, root of least magnitude

quadratic formula, root of least magnitude

## Usage

``` r
sel_root(a, b, c, dsign)
```

## Arguments

- a:

  numeric. Quadratic coefficient (a in a*x^2 + b*x + c = 0).

- b:

  numeric. Linear coefficient.

- c:

  numeric. Constant term.

- dsign:

  numeric. Sign of the discriminant term to select which root to return
  (+1 or -1).

## Value

The selected root of least magnitude.

## Examples

``` r
sel_root(a = 1, b = -3, c = 2, dsign = 1)
#> [1] 2
```
