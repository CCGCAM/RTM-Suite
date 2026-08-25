# Generate a variable correlated with an existing one

Draws a new variable `y = r*x + noise`, with noise scaled so that
`cor(x, y)` is approximately `r`. Negative results are clipped to 0
(this assumes `x`/`y` are non-negative quantities, e.g. pigment
concentrations – see [`get_distributionLUT`](get_distributionLUT.md)'s
`DepCab` option, which uses this to correlate Car with Cab).

## Usage

``` r
correlatedValue(x, r)
```

## Arguments

- x:

  numeric vector. The variable to correlate against.

- r:

  numeric (-1 to 1). Target correlation coefficient between `x` and the
  result.

## Value

A numeric vector the same length as `x`, correlated with it at
approximately `r`.

## Examples

``` r
Cab <- runif(200, 10, 80)
Car <- correlatedValue(x = Cab / 4, r = 0.8)
cor(Cab, Car)  # close to 0.8
#> [1] 0.9883911
```
