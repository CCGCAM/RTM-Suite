# Gradient scale Map a numeric vector to a color gradient

Gradient scale Map a numeric vector to a color gradient

## Usage

``` r
color.gradient(x, colors = c("red", "yellow", "green"), colsteps = 100)
```

## Arguments

- x:

  numeric vector. Values to map onto the color scale.

- colors:

  character vector. Colors defining the gradient, interpolated in order.
  Default `c("red","yellow","green")`.

- colsteps:

  integer. Number of discrete color steps in the gradient. Default 100.

## Value

a scale
