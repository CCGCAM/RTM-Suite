# `getRTMo` Calculates the spectra of hemisperical and directional observed visible and thermal radiation (fluxes E and radiances L), as well as the single and bi-directional gap probabilities

updates:

- 10 Sep 2007 (CvdT) - calculation of Rn

- 5 Nov 2007 - included observation direction

- 12 Nov 2007 - included abs. PAR spectrum output - improved calculation
  efficiency

- 13 Nov 2007 - written readme lines

- 11 Feb 2008 (WV&JT) - changed Volscat - small change in calculation
  Po,Ps,Pso (author:JT) - introduced parameter 'lazitab' - changed
  nomenclature - Appendix IV: cosine rule

- 04 Aug 2008 (JT) - Corrections for Hotspot effect in the probabilities

- 05 Nov 2008 (CvdT) - Changed layout

- 04 Jan 2011 - Included Pso function (Appendix IV) (JT&CvdT) - removed
  the analytical function (for checking)(JT&CvdT)

- 02 Oct 2012 (CvdT) - included incident PAR in output

- Jan/Feb 2013 (WV) - Major revision towards SCOPE version 1.40: -
  Parameters passed using structures - Improved interface with MODTRAN
  atmospheric data - Now also calculates 4-stream - reflectances rso,
  rdo, rsd and rdd analytically

- Apri 2013 (CvT) - improvements in variable names and descriptions

- Dec 2019 CvdT mSCOPE representation, lite option

## Usage

``` r
getRTMo(
  data.spectral,
  atmo,
  data.soil,
  data.leafopt,
  data.canopy,
  data.leafbio,
  data.angles,
  data.meteo,
  data.opts,
  canopy.model,
  get.plots = T
)
```

## Arguments

- data.spectral:

  information about wavelengths and resolutions

- atmo:

  MODTRAN atmospheric parameters

- data.soil:

  soil properties

- data.leafopt:

  leaf optical properties

- data.canopy:

  canopy properties (such as LAI and height)

- data.leafbio:

  leaf biochemical parameters (Cab, Car...)

- data.angles:

  viewing and observation angles

- data.meteo:

  has the meteorological variables. Is only used to correct, the total
  irradiance if a specific value is provided instead of the usual
  Modtran output.

- data.opts:

  simulation options. Here, the option

- canopy.model:

  Selection of canopy model. THe canopy models available are 'fourSAIL',
  'fourSAIL' and 'INFORM'. By default fourSAIL model will be used.

- get.plots:

  is true plot the intermediate plots

## References

Verhoef (1998), 'Theory of radiative transfer models applied in optical
remote sensing of vegetation canopies'. PhD Thesis Univ. Wageninegn.
Verhoef, W., Jia, L., Xiao, Q. and Su, Z. (2007) Unified optical -
thermal four - stream radiative transfer theory for homogeneous
vegetation canopies. IEEE Transactions on geoscience and remote sensing,
45,6. Verhoef (1985), 'Earth Observation Modeling based on Layer
Scattering Matrices', Remote sensing of Environment, 17:167-175.

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
rad <- getRTMo(data.spectral, atmo, data.soil, data.leafopt, data.canopy,
                data.leafbio, data.angles, data.meteo, data.opts,
                canopy.model = "fourSAIL", get.plots = FALSE)
} # }
```
