# Get atmospheric pressure (in hPa) as a function of altitude

This function calculates the atmospheric pressure at a given altitude
using the barometric formula.

## Usage

``` r
get.Altitude2Pa(altitude, Pa0 = 1013.25, Ta0 = 288.15)
```

## Arguments

- altitude:

  Altitude in meters.

- Pa0:

  Air pressure at sea level in hPa (default is 1013.25 hPa).

- Ta0:

  Air temperature at sea level in Kelvin (default is 288.15 K).

## Value

The atmospheric pressure at the given altitude in hPa.

## Examples

``` r
# Atmospheric pressure at 1000 meters altitude
get.Altitude2Pa(altitude = 1000)
#> [1] 899.9491

# Atmospheric pressure at 5000 meters altitude with custom sea-level pressure
get.Altitude2Pa(altitude = 5000, Pa0 = 1010)
#> [1] 558.2492
```
