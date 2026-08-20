# Leaf FLUSPECT-B-Cx model for SCOPE

`getFluspect.Cx.SCOPE` calculates reflectance and transmittance spectra
of a leaf using FLUSPECT-B, calculates reflectance and transmittance
spectra of a leaf using FLUSPECT-B, plus four excitation-fluorescence
matrices

## Usage

``` r
getFluspect.Cx.SCOPE(inputsLeaf, inputsOptipar, version = "Cx", step)
```

## Arguments

- inputsLeaf:

  a LUT with main plnt traits for Fluspect-B Cx

- inputsOptipar:

  internal parameters

- version:

  Leaf model uses for Kcal: 'Cx'

- step:

  a step for wavelengths for getting a matrix of Mf and Mb, by default
  uses a step of 1 nm, This step return matrix of 211x351 the SCOPE
  model uses a step of 5 for getting matrix of 53x71

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

Matlab version: Date: 2007-2020

Ported in R: Jun-2023

## Examples

``` r
if (FALSE) { # \dontrun{
inputs <- ToolsRTM::inputsRTMs
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = 1, setseed = 1234))
sim <- getFluspect.Cx.SCOPE(inputsLeaf = LUT, inputsOptipar = SCOPEinR::optipar, version = 'Cx')
} # }
```
