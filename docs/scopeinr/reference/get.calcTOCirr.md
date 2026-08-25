# get.calcTOCirr

`get.calcTOCirr` calculates the top-of-canopy (TOC) direct solar
(`Esun_`) and diffuse sky (`Esky_`) irradiance spectra. If `atmo`
already contains pre-computed `Esun_`/`Esky_` these are returned
directly; otherwise they are derived from the MODTRAN atmospheric
transmittance/radiance functions (T1-T16) evaluated at the SCOPE
wavelengths, optionally rescaled to match the observed broadband
incoming shortwave (`Rin`) and longwave (`Rli`) radiation.

## Usage

``` r
get.calcTOCirr(atmo, meteo, rdd, rsd, wl, nwl)
```

## Arguments

- atmo:

  list or data.frame. Atmospheric input, either pre-computed
  `Esun_`/`Esky_` spectra, or MODTRAN output containing wavenumber
  (`WN`) and transmittance/radiance functions
  `T1, T2, T3, T4, T5, T12, T16`.

- meteo:

  list. Meteorological data; must contain `Ta` (air temperature, deg C)
  and, when rescaling is required, `Rin` (incoming shortwave radiation,
  W m-2) and `Rli` (incoming longwave radiation, W m-2). Use `-999` for
  `Rin` to skip rescaling.

- rdd:

  numeric vector of length `nwl`. Top-of-canopy
  hemispherical-hemispherical reflectance spectrum.

- rsd:

  numeric vector of length `nwl`. Top-of-canopy
  directional-hemispherical reflectance spectrum.

- wl:

  numeric vector of length `nwl`. Wavelengths of the spectral domain
  (nm).

- nwl:

  integer. Number of wavelengths in the spectral domain.

## Value

A list with:

- Esun\_:

  numeric vector of length `nwl`. Top-of-canopy direct solar irradiance
  spectrum.

- Esky\_:

  numeric vector of length `nwl`. Top-of-canopy diffuse sky irradiance
  spectrum.

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
TOCirr <- get.calcTOCirr(atmo, meteo, rdd, rsd, wl, nwl)
} # }
```
