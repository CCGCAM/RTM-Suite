# Leaf FLUSPECT-B model for SCOPE

`getFluspect.B.SCOPE` calculates reflectance and transmittance spectra
of a leaf using FLUSPECT-B, calculates reflectance and transmittance
spectra of a leaf using FLUSPECT-B, plus four excitation-fluorescence
matrices

## Usage

``` r
getFluspect.B.SCOPE(inputsLeaf, inputsOptipar, version = "D")
```

## Arguments

- inputsLeaf:

  a LUT for FLUXspect

- inputsOptipar:

  internal parameters

- version:

  Leaf model uses for Kcal: 'D' refers as PROSPECT-D; 'PRO', refers as
  PROSPECT-PRO

## Value

a list which contains:

- refl (reflectance)

- tran (transmittance)

- Mb (backward scattering fluorescence matrix, I for PSI and II for
  PSII)

- Mf (forward scattering fluorescence matrix, I for PSI and II for PSII)

## Author

Wout Verhoef, Christiaan van der Tol, Joris Timmermans, Nastassia Vilfan
(Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
inputs <- ToolsRTM::inputsRTMs
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
sim <- getFluspect.B.SCOPE(inputsLeaf = LUT, inputsOptipar = SCOPEinR::optipar, version = 'D')
} # }
```
