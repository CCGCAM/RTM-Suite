# Get.merge.Output generate the LUT table adding the apparent reflectance, radiance or fluorescence emission

Get.merge.Output generate the LUT table adding the apparent reflectance,
radiance or fluorescence emission

## Usage

``` r
get.merge.SCOPE(
  paths.sims = NULL,
  inputs,
  options = c("vegetation", "fluxes", "fluorescence", "aPAR"),
  spectra.data = "rfl"
)
```

## Arguments

- paths.sims:

  character. Path to the folder of a SCOPE simulation run (containing
  the `Parameters/` subfolder and the output CSVs to merge).

- inputs:

  A vector specifying the input parameters to include in the final
  lookup table. If missing or NULL all values will be returned.

- options:

  A character vector specifying which additional data tables to merge
  into the LUT.

- spectra.data:

  A string specifying the type of spectral data to merge ("rfl.app",
  "rfl", "Lo", "Lo.sif", "sif", "sif.leaves", "sif.shaded",
  "sif.sunlit", "sif.scattered", "sif.hemis").

## Value

The merged lookup table.

## Examples

``` r
if (FALSE) { # \dontrun{
paths.sims <- "path/to/sims/directory/"
inputs <- c('Cab','Car','Anth','LMA','EWT','Vcmax25','CBC','Prot','Cs',
             'Cbrown','Cx','LIDFa','LIDFb','LAI','Rin','Rli')
merged_LUT <- get.merge.SCOPE(paths.sims = paths.sims, inputs = inputs,
                               options = c("vegetation", "fluxes", "fluorescence", "aPAR"),
                               spectra.data = "rfl.app")
} # }
```
