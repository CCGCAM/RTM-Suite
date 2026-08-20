# get.spectral.convolution

This function performs spectral convolution on the input irradiance (or
radiance) dataframe using the spectral response function of the given
sensor.

## Usage

``` r
get.spectral.convolution(df.irradiance, sensor.i, get.plots = T)
```

## Arguments

- df.irradiance:

  A dataframe containing high-resolution irradiance (or radiance) data
  with wavelengths and corresponding values.

- sensor.i:

  A dataframe or list containing the spectral response function
  (wavelengths and sensitivity) of the sensor.

- get.plots:

  A boolean indicating whether to generate a plot showing the
  convolution results. Default is TRUE.

## Value

A dataframe of convolved values for each sensor band, optionally with a
plot if `get.plots` is TRUE.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage with sample irradiance data and sensor response function:
df.irradiance <- data.frame(wavelength = seq(400, 700, 1), irradiance = rnorm(301))  # Example irradiance data
sensor.i <- data.frame(wavelength = seq(400, 700, 1), response = runif(301))        # Example sensor response
get.spectral.convolution(df.irradiance, sensor.i, get.plots = TRUE)

} # }
```
