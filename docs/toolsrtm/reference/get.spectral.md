# Get spectral bandset for interpolation

This function retrieves the spectral bandset for various sensors to be
used in interpolation or analysis.

## Usage

``` r
get.spectral(sensor = "Sentinel2a")
```

## Arguments

- sensor:

  A character string with the name of the sensor; available options are:
  "ALI", "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8",
  "MODIS", "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b",
  "WorldView2-4", "WorldView2-8".

## Value

A data frame containing the spectral bands for the specified sensor.

## Examples

``` r
df.sentinel2a <- get.spectral(sensor='Sentinel2a')
```
