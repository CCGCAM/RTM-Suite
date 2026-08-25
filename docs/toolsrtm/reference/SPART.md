# Compute Top-of-Canopy and Top-of-Atmosphere Reflectance and Radiance

This function computes the top-of-canopy (TOC) reflectance,
top-of-atmosphere (TOA) reflectance, and TOA radiance using the
Soil-Plant-Atmosphere Radiative Transfer (SPART) model based on the
provided input parameters and models.

## Usage

``` r
SPART(
  inputLUT,
  optipar = NULL,
  CanopyModel = "fourSAIL",
  LeafModel = "PROSPECT-PRO",
  sensor.i = NULL,
  df.irradiance = NULL,
  rsoil = NULL,
  get.plots = T
)
```

## Arguments

- inputLUT:

  A lookup table (LUT) containing the necessary parameters for the
  radiative transfer model.

- optipar:

  Optional parameters for optical calculations (default is NULL).

- CanopyModel:

  A character string specifying the canopy model to use. Default is
  'fourSAIL'.

- LeafModel:

  A character string specifying the leaf model to use. Default is
  'PROSPECT-PRO'.

- sensor.i:

  A dataset with sensor information for which the calculations are
  performed. Options include "LANDSAT4.TM", "LANDSAT5.TM", "LANDSAT5.TM,
  "LANDSAT7.ETM", "LANDSAT8.OLI","Sentinel2A.MSI", "Sentinel2B.MSI",
  "Sentinel3A.OLCI", "Sentinel3B.OLCI"and "TerraAqua.MODIS"

- df.irradiance:

  A data frame containing irradiance data for the calculations.

- rsoil:

  Optional numeric vector, length 2001 (400-2400nm, 1nm step – matching
  `get.spectra.spart(getSpectral = TRUE)$reg1`), of soil reflectance.
  When supplied, this is used directly as the soil spectrum instead of
  computing one from the built-in BSM (Brightness-Shape-Moisture) soil
  model – e.g. to feed in a spectrum from
  [`get.marmit.rsoil`](get.marmit.rsoil.md) (pass `wl.out = 400:2400`
  there to get a matching-length vector). The thermal region (2400nm+)
  is still padded the same way as the BSM path (a flat `0.06`
  reflectance). Default `NULL` (use BSM).

- get.plots:

  A boolean indicating whether to generate plots of the results. Default
  is TRUE.

## Value

A list containing TOC reflectance, TOA reflectance, and TOA radiance.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires a real LUT (see get.LUTfromRanges() / getLUT()) and real
# irradiance data — not run automatically as no example dataset ships
# with the package for this.
LUT <- get.LUTfromRanges(...)          # build/load a real LUT here
df.irradiance <- data.frame(...)        # real irradiance data here
results <- SPART(inputLUT = LUT, optipar = NULL,
                  CanopyModel = 'fourSAIL', LeafModel = 'PROSPECT-PRO',
                  sensor.i = ToolsRTM::TerraAqua.MODIS, df.irradiance = df.irradiance,
                  get.plots = TRUE)
} # }
```
