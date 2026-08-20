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
  convoluted. Must have a `wave` column (wavelength, nm) and an `rfl`
  column (reflectance) – e.g.
  `data.frame(wave = dataSpec_PDB[, 1], rfl = my_reflectance)`.

- sensor.i:

  The sensor to convolve onto – one of the package's bundled sensor
  objects, e.g. `ToolsRTM::Sentinel2A.MSI`, `ToolsRTM::LANDSAT8.OLI`,
  `ToolsRTM::TerraAqua.MODIS` (also accepts the bare sensor name as a
  character string, e.g. `"LANDSAT8.OLI"`, which is resolved to the
  matching bundled object automatically). Supported sensors:
  "LANDSAT4.TM", "LANDSAT5.TM", "LANDSAT7.ETM", "LANDSAT8.OLI",
  "Sentinel2A.MSI", "Sentinel2B.MSI", "Sentinel3A.OLCI",
  "Sentinel3B.OLCI", and "TerraAqua.MODIS".

- get.plots:

  A boolean value indicating whether to generate plots of the spectral
  convolution. Default is TRUE.

## Value

A data frame containing the convoluted reflectance values for the
specified sensor, with columns `wave` (band center wavelength, nm) and
`RFL` (convolved reflectance).

## Examples

``` r
df <- data.frame(wave = ToolsRTM::dataSpec_PDB[, 1],
                  rfl = 0.05 + 0.3 * ToolsRTM::dataSpec_PDB[, 1] / 2500) # toy reflectance
convoluted_results <- get.spectral.convolution.rfl(df, sensor.i = "LANDSAT8.OLI", get.plots = FALSE)
```
