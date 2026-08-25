# calculates the saturated vapour pressure at temperature T (degrees C) and the derivative of es to temperature s (kPa/C)

calculates the saturated vapour pressure at temperature T (degrees C)
and the derivative of es to temperature s (kPa/C)

## Usage

``` r
satvap(Temp)
```

## Arguments

- Temp:

  numeric. Air/leaf temperature, in degrees Celsius.

## Value

es: the saturated vapour pressure at Temp, in hPa/mbar.

## Details

Date: 2003.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
satvap(20)
#> [1] 23.37787
```
