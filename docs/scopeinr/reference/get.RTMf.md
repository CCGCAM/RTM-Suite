# get.RTMf `get.RTMf` Calculates the spectrum of fluorescent radiance in the observer's direction and also the TOC spectral hemispherical upward Fs flux.

get.RTMf `get.RTMf` Calculates the spectrum of fluorescent radiance in
the observer's direction and also the TOC spectral hemispherical upward
Fs flux.

## Usage

``` r
get.RTMf(
  data.spectral,
  data.rad,
  data.soil,
  data.leafopt,
  data.canopy,
  data.gap,
  data.angles,
  data.etau,
  data.etah,
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

- data.leafopt:

  leaf optical properties

- data.canopy:

  canopy properties (such as LAI and height)

- data.gap:

  probabilities of direct light penetration and viewing

- data.angles:

  viewing and observation angles

- data.etau:

  relative fluorescence emission efficiency for sunlit leaves

- data.etah:

  relative fluorescence emission efficiency for shaded leaves

- get.plots:

  logical. If `TRUE`, intermediate diagnostic plots are produced.
  Default `TRUE`.

## Value

a rad object with a large number of radiative fluxes: spectrally
distributed and integrated, and canopy radiative transfer coefficients.

## Details

Date: 12 Dec 2007. Update history:

- 26 Aug 2008 CvdT small correction to matrices

- 07 Nov 2008 CvdT changed layout

- 19 Mar 2009 CvdT major corrections: lines 95-96,101-107, and 119-120.

- 07 Apr 2009 WV & CvdT major correction: lines 89-90, azimuth
  dependence was not there in previous verions (implicit assumption of
  azimuth(solar-viewing) = 0). This has been corrected

- May-June 2012 WV & CvdT Add calculation of hemispherical Fs fluxes

- Jan-Feb 2013 WV Inputs and outputs via structures for SCOPE Version
  1.40

- Aug-Oct 2016 PY Re-write the calculation of emitted SIF of each layer.
  It doesnt use loop at all. with the function bsxfun, the calculation
  is much faster

- Oct 2017-Feb 2018 PY Re-write the RTM of fluorescence

- Jan 2020 CvdT Modified to include 'lite' option, mSCOPE representation

- 25 Jun 2020 PY Po, Ps, Pso. fix the problem we have with the oblique
  angles above 80 degrees

## Author

    Wout Verhoef and Christiaan van der Tol  (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
get.RTMf(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
         data.gap, data.angles, data.etau, data.etah)
} # }
```
