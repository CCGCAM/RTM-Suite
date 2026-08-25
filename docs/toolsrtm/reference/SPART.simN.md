# Compute Top-of-Canopy and Top-of-Atmosphere Reflectance and Radiance with Simulations

This function computes the top-of-canopy (TOC) reflectance,
top-of-atmosphere (TOA) reflectance, and TOA radiance using the
Soil-Plant-Atmosphere Radiative Transfer (SPART) model based on the
provided input parameters and models. It also supports multiple
simulations to assess variability.

## Usage

``` r
SPART.simN(
  inputLUT,
  optipar = NULL,
  CanopyModel = "fourSAIL",
  LeafModel = "PROSPECT-PRO",
  sensor.i = NULL,
  df.irradiance = NULL,
  rsoil = NULL
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

  Optional numeric vector, length 2001 (400-2400nm, 1nm step), of soil
  reflectance to use instead of the built-in BSM soil model – e.g. from
  [`get.marmit.rsoil`](get.marmit.rsoil.md). See [`SPART`](SPART.md) for
  details; behaves identically here. Default `NULL` (use BSM).

## Value

A list containing TOC reflectance, TOA reflectance, and TOA radiance for
each simulation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires a real LUT (see get.LUTfromRanges() / getLUT()) and real
# irradiance data — not run automatically as no example dataset ships
# with the package for this.
LUT <- get.LUTfromRanges(...)          # build/load a real LUT here
df.irradiance <- data.frame(...)        # real irradiance data here
results <- SPART.simN(inputLUT = LUT, optipar = NULL,
                       CanopyModel = 'fourSAIL', LeafModel = 'PROSPECT-PRO',
                       sensor.i = ToolsRTM::TerraAqua.MODIS, df.irradiance = df.irradiance,
                       N = 10)
} # }
```
