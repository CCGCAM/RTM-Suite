# calculates the saturated vapour pressure at temperature T (degrees C) and the derivative of es to temperature s (kPa/C)

calculates the saturated vapour pressure at temperature T (degrees C)
and the derivative of es to temperature s (kPa/C)

## Usage

``` r
slope_satvap(Temp)
```

## Arguments

- Temp:

  numeric. Air/leaf temperature, in degrees Celsius.

## Value

The derivative of saturated vapour pressure with respect to temperature
(slope), in kPa/degree C.

## Details

Date: 2003.

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
slope_satvap(20)
#> $es
#> [1] 23.37787
#> 
#> $s
#> [1] 1.447115
#> 
```
