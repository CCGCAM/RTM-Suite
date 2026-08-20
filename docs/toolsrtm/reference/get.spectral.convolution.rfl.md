# Get Spectral Convolution for Reflectance

This function performs spectral convolution on reflectance data using a
specified sensor's spectral response function. The convolution is
essential for simulating how well the sensor captures the true
reflectance values based on its spectral characteristics.

## Usage

``` r
get.spectral.convolution.rfl(df, sensor.i, get.plots = T)
```

## Arguments

- df:

  A data frame containing the high-resolution reflectance data to be
  convoluted.

- sensor.i:

  A character string with sensor information for which the calculations
  are performed. Supported sensors should have predefined spectral
  response functions. Options include: "LANDSAT4.TM", "LANDSAT5.TM",
  "LANDSAT7.ETM", "LANDSAT8.OLI", "Sentinel2A.MSI", "Sentinel2B.MSI",
  "Sentinel3A.OLCI", "Sentinel3B.OLCI", and "TerraAqua.MODIS".

- get.plots:

  A boolean value indicating whether to generate plots of the spectral
  convolution. Default is TRUE.

## Value

A data frame containing the convoluted reflectance values for the
specified sensor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example data frame with reflectance data
df <- data.frame(Wavelength = seq(400, 2500, by = 10),
                 Reflectance = runif(211)) # Simulated reflectance values

# Perform spectral convolution for a specific sensor
convoluted_results <- get.spectral.convolution.rfl(df, sensor.i = "LANDSAT8.OLI", get.plots = TRUE)
} # }
```
