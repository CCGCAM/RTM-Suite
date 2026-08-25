# Brightness-Shape-Moisture soil model

`getBSM.toolsRTM` Calculates the spectrum of fluorescent radiance in the
observer's direction and also the TOC spectral hemispherical upward Fs
flux.

## Usage

``` r
getBSM.toolsRTM(soilpar, spec, emp)
```

## Arguments

- soilpar:

  list. Soil parameters: BSMBrightness (soil brightness), BSMlat
  (spectral shape latitude, 20-40 deg), BSMlon (spectral shape
  longitude, 45-65 deg).

- spec:

  list. Spectral inputs: GSV (Global Soil Vectors spectra), Kw (water
  absorption spectrum), nw (water refraction index spectrum).

- emp:

  list. Empirical parameters, including SMp (soil moisture volume
  percentage, 5-55).

## Value

    Wet soil reflectance spectra across 400 nm to 2400 nm

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
save.soil <- getBSM.toolsRTM(soilpar, spec, emp)
} # }
```
