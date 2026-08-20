# Apply an atmospheric correction using SMAC method (Rahman and Dedieu, 1994)

This function applies the Simplified Method for Atmospheric Correction
(SMAC) to satellite imagery based on the provided LUT (Look-Up Table)
and sensor characteristics.

## Usage

``` r
get.smac(inputLUT = LUT, sensor)
```

## Arguments

- inputLUT:

  A data frame containing the LUT used for atmospheric correction.

- sensor:

  A dataset with sensor information for which the calculations are
  performed. Options include "LANDSAT4.TM", "LANDSAT5.TM", "LANDSAT5.TM,
  "LANDSAT7.ETM", "LANDSAT8.OLI","Sentinel2A.MSI", "Sentinel2B.MSI",
  "Sentinel3A.OLCI", "Sentinel3B.OLCI"and "TerraAqua.MODIS"

## Value

A data frame with atmospheric corrected values.

## Examples

``` r
if (FALSE) { # \dontrun{
inputLUT = ToolsRTM::inputsSPART
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = ToolsRTM::inputsSPART, nLUT = 1, setseed = 1234))
corrected_data <- get.smac(inputLUT, sensor = ToolsRTM::TerraAqua.MODIS)

} # }
```
