# get.RTMt.planck `get.RTMt.planck`analogue to get.RTMt.sb, this function calculates total outgoing radiation in hemispherical direction and total absorbed radiation per leaf and soil component. Radiation is integrated over the whole thermal spectrum with Stefan-Boltzman's equation. This function is a simplified version of 'get.RTMt.planck', and is less time consuming since it does not do the calculation for each wavelength separately.

Date: 05 Nov 2007 Update:

- 13 Nov 2007

- 16 Nov 2007 CvdT: improved calculation of net radiation.

- 27 Mar 2008 JT: added directional calculation of radiation.

- 24 Apr 2008 JT: introduced dx as thickness of layer (see parameters).

- 31 Oct 2008 JT: introduced optional directional calculation.

- 31 Oct 2008 JT: changed initialisation of F1 and F2 -\> zeros.

- 07 Nov 2008 CvdT: changed layout.

- 16 Mar 2009 CvdT: removed Tbright calculation.

- Feb 2013 WV: introduces structures for version 1.40.

- 04 Dec 2019 CvdT: adapted for SCOPE-lite.

- 17 Mar 2020 CvdT: mSCOPE representation.

## Usage

``` r
get.RTMt.planck(
  data.spectral,
  data.rad,
  data.soil,
  data.leafbio,
  data.leafopt,
  data.canopy,
  data.gap,
  Tcu,
  Tch,
  Tsu,
  Tsh,
  get.plots = T
)
```

## Arguments

- data.spectral:

  information about wavelengths and resolutions

- data.rad:

  a large number of radiative fluxes: spectrally distributed and
  integrated, and canopy radiative transfer coefficients.

- data.soil:

  soil properties

- data.leafbio:

  leaf properties (Cab....)

- data.leafopt:

  leaf optical properties

- data.canopy:

  canopy properties (such as LAI and height)

- data.gap:

  probabilities of direct light penetration and viewing

- Tcu:

  Temperature of sunlit leaves (oC), (13x36x60)

- Tch:

  Temperature of shaded leaves (oC), (13x36x60)

- Tsu:

  Temperature of sunlit soil (oC), (1)

- Tsh:

  Temperature of shaded soil (oC), (1)

- get.plots:

  is true plot the intermediate plots

## Value

a large number of radiative fluxes: spectrally distributed and
integrated, and canopy radiative transfer coefficients., Here thermal
fluxes are added.

## Author

Wout Verhoef and Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
data.rad <- get.RTMt.planck(data.spectral, data.rad, data.soil, data.leafbio,
                             data.leafopt, data.canopy, data.gap,
                             Tcu, Tch, Tsu, Tsh, get.plots = FALSE)
} # }
```
