# get.RTM - Run Radiative Transfer Model (RTM)

This function runs a Radiative Transfer Model (RTM) using specified
models for the leaf and canopy, and allows for optional BRDF correction
and plot generation.

## Usage

``` r
get.RTM(
  LUT.table,
  leaf.model = "fluspect-CX",
  canopy.model = "fourSAIL",
  option = "Foward-Mode",
  BRDF = T,
  get.plots = T
)
```

## Arguments

- LUT.table:

  A lookup table with input variables for the RTM model (e.g., leaf
  optical properties, canopy structure).

- leaf.model:

  The leaf model to be used (default is 'fluspect-CX').

- canopy.model:

  The canopy model to be used (default is 'fourSAIL').

- option:

  A string specifying the option or mode (e.g., 'Forward-Mode',
  'Inverse-Mode'). Default is 'Forward-Mode'.

- BRDF:

  A boolean to indicate whether Bidirectional Reflectance Distribution
  Function (BRDF) correction should be applied. Default is TRUE.

- get.plots:

  A boolean to indicate whether to generate plots of the results.
  Default is TRUE.

## Value

A list with results of the RTM, potentially including reflectance,
transmittance, and optional plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run RTM with default settings
get.RTM(LUT.table)

# Run RTM without BRDF correction and plots
get.RTM(LUT.table, BRDF = FALSE, get.plots = FALSE)

} # }
```
