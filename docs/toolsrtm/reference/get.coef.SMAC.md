# Get the coefficients for the SMAC model

This function retrieves the atmospheric correction coefficients for
different satellite sensors including Landsat, Sentinel-2, Sentinel-3,
and MODIS (Aqua and Terra).

## Usage

``` r
get.coef.SMAC(sensor)
```

## Arguments

- sensor:

  A character string specifying the sensor name. Options include
  "Landsat", "Sentinel-2", "Sentinel-3", "MODIS_Aqua", and
  "MODIS_Terra".

## Value

A list of coefficients corresponding to the specified sensor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get coefficients for Landsat
coef_landsat <- get.coef.SMAC("Landsat")

# Get coefficients for Sentinel-2
coef_sentinel2 <- get.coef.SMAC("Sentinel-2")
} # }
```
