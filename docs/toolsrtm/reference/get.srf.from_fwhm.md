# Generate Spectral Response Function from Full Width at Half Maximum (FWHM)

This function generates a spectral response function based on a Gaussian
distribution, which is defined by the center wavelength and standard
deviation (sd). A normalization based on the maximum (max) value is
applied. The sensors considered in this function include:

- **Sentinel2a** and **Sentinel2b**: For these sensors, the function
  uses the spectral response data stored in `ToolsRTM::str.sentinel2a`
  and `ToolsRTM::str.sentinel2b`, respectively. These datasets provide
  the necessary characteristics to accurately model their spectral
  responses.

- **PRISMA**: The function utilizes the original resampling function for
  PRISMA, with the spectral response characteristics stored in
  `ToolsRTM::srf.prisma`. This allows for precise calculations based on
  the unique properties of the PRISMA sensor.

- Other sensors (e.g., ALI, Hyperion, Landsat4-8, MODIS, Quickbird,
  RapidEye, WorldView2-4, WorldView2-8) use the Full Width at Half
  Maximum (FWHM) function for their spectral response calculations,
  based on their respective characteristics.

## Usage

``` r
get.srf.from_fwhm(sensor = "", save = T, path.out, get.plot = FALSE)
```

## Arguments

- sensor:

  A string specifying the sensor name. Options include: "ALI",
  "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8", "MODIS",
  "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b", "WorldView2-4",
  "WorldView2-8".

- save:

  Logical indicating whether to save SRF file (default is TRUE).

- path.out:

  A string specifying the output path for the CSV files.

- get.plot:

  Logical indicating whether to plot the spectral response curves
  (default is FALSE).

## Value

A data frame containing the spectral response functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate spectral response function for Sentinel-2b
srf.matrix <- get.srf.from_fwhm(sensor = "Sentinel2b", save=T, path.out = "Tables", get.plot = TRUE)

} # }
```
