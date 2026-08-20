# Convolve reflectance onto a sensor using a plain per-band SRF table

[`get.spectral.convolution.rfl()`](get.spectral.convolution.rfl.md)
needs a sensor object with SMAC atmospheric-correction coefficients
bundled (see its own docs for the 9 supported sensors) – fine for
Landsat/OLCI/MODIS, but PRISMA has no SMAC coefficients at all, and
Sentinel-2A/B additionally ship a second, plain, publisher-original
spectral response function table alongside their SMAC bundle. This
function convolves onto any of those "just a per-band SRF table, nothing
else" sensors instead – one shared implementation for
`ToolsRTM::srf.prisma`, `ToolsRTM::srf.sentinel2a`, and
`ToolsRTM::srf.sentinel2b` (generalizes what used to be two separate,
app-only helper functions, `convolve_prisma()`/`convolve_smac_sensor()`,
in the AEO-Course PROSAIL Shiny app's own app.R).

## Usage

``` r
get.spectral.convolution.srf(df, srf, fwhm = NULL, get.plots = FALSE)
```

## Arguments

- df:

  A data frame with a `wave` column (wavelength, nm) and an `rfl` column
  (reflectance) – same convention as
  [`get.spectral.convolution.rfl()`](get.spectral.convolution.rfl.md).

- srf:

  A plain per-band SRF table: its FIRST column is wavelength (nm, any
  column name), every other column is one sensor band's SRF weight at
  that wavelength (e.g. `ToolsRTM::srf.prisma`,
  `ToolsRTM::srf.sentinel2a`, `ToolsRTM::srf.sentinel2b`).

- fwhm:

  Optional data frame with a `fwhm` column, one row per band in `srf`
  (same order), e.g. `ToolsRTM::fwhm.prisma` – when supplied, its
  (precisely interpolated) FWHM values are used instead of the coarser
  half-max-crossing estimate this function would otherwise derive
  directly from `srf`'s own sampled weight profile. Sentinel-2A/B have
  no equivalent bundled FWHM table, so leave this `NULL` for them.

- get.plots:

  logical, plot the convolved spectrum? Default `FALSE`.

## Value

A data frame with one row per SRF column: `band` (index), `wl`
(SRF-weighted mean center wavelength, nm), `fwhm` (full width at half
maximum, nm), `RFL` (convolved reflectance).

## Examples

``` r
df <- data.frame(wave = ToolsRTM::dataSpec_PDB[, 1],
                  rfl = 0.05 + 0.3 * ToolsRTM::dataSpec_PDB[, 1] / 2500)
prisma_bands <- get.spectral.convolution.srf(df, ToolsRTM::srf.prisma,
                                              fwhm = ToolsRTM::fwhm.prisma)
s2a_bands <- get.spectral.convolution.srf(df, ToolsRTM::srf.sentinel2a)
```
