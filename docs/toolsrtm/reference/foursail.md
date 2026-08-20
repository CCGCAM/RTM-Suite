# fourSAIL model coupled with several leaf models

`foursail` fourSAIL simulation based on a set of combinations of input
parameters

## Usage

``` r
foursail(inputLUT, rsoil, LeafModel = "PROSPECT-PRO", spectrum.all = T)
```

## Arguments

- inputLUT:

  LUT table with distribution of biophysical parameters used as input
  parameters in the model

- rsoil:

  numeric. Soil reflectance

- LeafModel:

  Version of PROSPECT model, Liberty model and fluspect model. For
  PROSPECT model: 'PROSPECT-PRO' or 'PROSPECT-D' is accepted. By default
  'PROSPECT-PRO' is used. fluspect model is valid in two options:
  'Fluspect-B' and 'Fluspect-B-Cx'. Liberty as 'Liberty'

- spectrum.all:

  a boolean value, False is for SPART and Fluspect Models (400-2400 nm),
  True for the PROSPECT and Liberty models (400-2500 nm) Liberty model
  is leaf radiative transfer model designed for conifer needles. Uses
  'Liberty' Fluspect-B model is leaf radiative transfer model designed
  for adding fluorescence emission. Uses 'Fluspect-B' or
  'Fluspect-B-Cx'.

## Value

list. rdot,rsot,rddt,rsdt

rdot: hemispherical-directional reflectance factor in viewing direction
rsot: bi-directional reflectance factor rsdt: directional-hemispherical
reflectance factor for solar incident flux rddt: bi-hemispherical
reflectance factor

## References

Verhoef W & Bach H, 2007. Coupled soil–leaf-canopy and atmosphere
radiative transfer modeling to simulate hyperspectral multi-angular
surface reflectance and TOA radiance data. Remote Sensing of
Environment, 109:166-182. doi:10.1016/j.rse.2006.12.013

Verhoef W, Jia L, Xiao Q & Su Z, 2007. Unified optical-thermal
four-stream radiative transfer theory for homogeneous vegetation
canopies. IEEE Transactions in Geosciences and Remote Sensing,
45:1808–1822. https://doi.org/10.1109/TGRS.2007.895844

Jacquemoud S, Verhoef W, Baret F, Bacour C, Zarco-Tejada PJ, Asner GP,
François C & Ustin SL, 2009. PROSPECT+ SAIL models: A review of use for
vegetation characterization. Remote Sensing of Environment, 113:S56–S66.
https://doi.org/doi:10.1016/j.rse.2008.01.026

Berger K, Atzberger C, Danner M, D’Urso G, Mauser W, Vuolo F & Hank T
2018. Evaluation of the PROSAIL Model Capabilities for Future
Hyperspectral Model Environments: A Review Study. Remote Sensing, 10:85.
https://doi.org/10.3390/rs10010085

Authors:

## Author

Wout Verhoef, Bach H. , JJean-Baptiste Feret (Original version in
Matlab) This version is also included in prospect package

Carlos Camino (Ported version into R wit modification from fourSAIL
model in prospect package)

The fourSAIL model is based on a version provided by Wout Verhoef et al.
(2007)

original version downloadable at
http://teledetection.ipgp.jussieu.fr/prosail/

Improved and extended version of SAILH model that avoids numerical
singularities and works more efficiently if only few parameters change.
