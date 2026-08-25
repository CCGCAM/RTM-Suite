# Cost function for fitting Fluspect leaf parameters to measured reflectance/transmittance

`get.COST_4Fluspect.for.SCOPE` builds the residual vector between
simulated and measured leaf reflectance/transmittance for a given
parameter vector, intended for use inside a non-linear least-squares
optimizer when calibrating leaf biochemistry/structure parameters
against measured spectra. Note: this function calls an unqualified
`fluspect_B_CX_PSI_PSII_combined()`, which is neither defined in this
package nor exported by ToolsRTM - it will fail with "could not find
function" until that dependency is resolved. Not exported; needs review
before use.

## Usage

``` r
get.COST_4Fluspect.for.SCOPE(params, measurement, input)
```

## Arguments

- params:

  numeric vector of length 8, in order: Cab, Cdm, Cw, Cs, Cca, Cant, Cx,
  N - candidate leaf biochemistry/structure parameters.

- measurement:

  list. Measured spectra, with elements `refl`, `tran` (measured
  reflectance/transmittance) and `std` (measurement uncertainty), all as
  a function of wavelength.

- input:

  list of length 6: `leafbio` (baseline leaf parameters), `optipar`
  (Fluspect optical parameters), `spectral` (spectral configuration,
  including `wlP` model wavelengths and `wlM` measurement wavelengths),
  `include` (list of 0/1 flags selecting which parameters in `params`
  are actually varied), `target` (character, `"0"` for combined
  reflectance+transmittance residuals, `"1"` for reflectance only,
  otherwise transmittance only), and `range` (list with `wlmin`/`wlmax`
  defining the wavelength range used to compute residuals).

## Value

A list with: `er` (residual vector, or two-column matrix of
reflectance/transmittance residuals when `target == "0"`), `refl`/`tran`
(simulated reflectance/transmittance interpolated to `spectral$wlM`),
and `leafopt` (full Fluspect output).

## Author

Carlos Camino

## Examples

``` r
if (FALSE) { # \dontrun{
cost <- get.COST_4Fluspect.for.SCOPE(params, measurement, input)
} # }
```
