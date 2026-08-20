# Brightness-Shape-Moisture soil model

`getBSM` Calculates the spectrum of fluorescent radiance in the
observer's direction and also the TOC spectral hemispherical upward Fs
flux.

## Usage

``` r
getBSM(soilpar, spec, emp)
```

## Arguments

- soilpar:

  list. Soil parameters: `BSMBrightness` (overall soil brightness),
  `BSMlat` (spectral shape "latitude", typical range 20-40 deg),
  `BSMlon` (spectral shape "longitude", typical range 45-65 deg).

- spec:

  list. Spectral inputs: `GSV` (Global Soil Vectors basis spectra,
  matrix nwl x 3), `Kw` (water absorption coefficient spectrum, nwl x
  1), `nw` (refraction index of water spectrum, nwl x 1).

- emp:

  list. Empirical wetting parameters: `SMp` (soil moisture volume
  percentage, range 5-55), `SMC` (soil moisture capacity, recommended
  0.25), `film` (effective optical thickness of a single water film,
  recommended 0.015).

## Value

    Wet soil reflectance spectra across 400 nm to 2400 nm

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
save.soil <- getBSM(soilpar, spec, emp)
} # }
```
