# get.RTMz `get.RTMz` Calculates the small modification of TOC outgoing radiance due to the conversion of Violaxanthin into Zeaxanthin in leaves

Date: 08 Dec 2016 Update:

- 17 Mar 2020 CvdT: added cluming, mSCOPE representation

- 25 Jun 2020 CvdT: Po, Ps, Pso. fix the problem we have with the
  oblique angles above 80 degrees

## Usage

``` r
get.RTMz(
  data.spectral,
  data.rad,
  data.soil,
  data.leafopt,
  data.canopy,
  data.gap,
  data.angles,
  data.Knu,
  data.Knh,
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

- data.Knu:

  ...kn for sunlit leaves (data.bcu\$Kn)

- data.Knh:

  .. kn for shaded leavess (data.bch\$Kn)

- get.plots:

  is true plot the intermediate plots

## Value

a rad object with a large number of radiative fluxes: spectrally
distributed and integrated, and canopy radiative transfer coefficients.

## Author

Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
data.rad <- get.RTMz(data.spectral, data.rad, data.soil, data.leafopt,
                      data.canopy, data.gap, data.angles, data.Knu, data.Knh,
                      get.plots = FALSE)
} # }
```
