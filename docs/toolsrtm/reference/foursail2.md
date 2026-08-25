# Performs fourSAIL2 + PROSPECT simulation based on a set of combinations of input parameters

Performs fourSAIL2 + PROSPECT simulation based on a set of combinations
of input parameters

## Usage

``` r
foursail2(
  LUT_GB = NULL,
  inputLUT,
  rsoil,
  PROSPECTversion = "PRO",
  FieldObserv = NULL,
  LeafModel = NULL
)
```

## Arguments

- LUT_GB:

  dataframe Includes distribution of biophysical parameters used as
  green vegetation (first column). Includes distribution of biophysical
  parameters used (second column)

- inputLUT:

  LUT table with distribution of biophysical parameters used as input
  parameters in the model

- rsoil:

  numeric. Soil reflectance

- PROSPECTversion:

  Legacy switch, kept for backward compatibility: 'PRO' or anything else
  ('D'). Ignored if `LeafModel` is given.

- FieldObserv:

  logical/NULL. If `NULL` (default), the "green vegetation" leaf optics
  are taken from `inputLUT` itself (i.e. the same leaf spectrum used for
  the main canopy), overriding whatever `LUT_GB[1,]` would otherwise
  give. If not `NULL`, the green/brown vegetation leaf optics are taken
  as computed from `LUT_GB` instead – for matching a field-observed
  green/brown fraction defined independently of `inputLUT`.

- LeafModel:

  character. One of 'PROSPECT-PRO', 'PROSPECT-D', 'Liberty',
  'Fluspect-B', 'Fluspect-B-Cx' – same 5 leaf models
  [`foursail`](foursail.md) and [`inform`](inform.md) support. Defaults
  to `PROSPECTversion`-derived value if not given.

## Value

list. rdot,rsot,rddt,rsdt rdot: hemispherical-directional reflectance
factor in viewing direction rsot: bi-directional reflectance factor
rsdt: directional-hemispherical reflectance factor for solar incident
flux rddt: bi-hemispherical reflectance factor alfast: canopy
absorptance for direct solar incident flux alfadt: canopy absorptance
for hemispherical diffuse incident flux

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

Verhoef W.

Bach H.

Authors of the R version:

Jean-Baptiste Feret

The fourSAIL model is based on a version provided by Wout Verhoef et al.
(2007)

original version downloadable at
http://teledetection.ipgp.jussieu.fr/prosail/

Improved and extended version of SAILH model that avoids numerical
singularities and works more efficiently if only few parameters change.
