# Generate a normalized Gaussian spectral response function

Generate a normalized Gaussian spectral response function

## Usage

``` r
resample_fun(center, wl, fwhm)
```

## Arguments

- center:

  numeric. Center wavelength (nm) of the spectral band.

- wl:

  numeric vector. Wavelengths (nm) at which to evaluate the response
  function.

- fwhm:

  numeric. Full width at half maximum (nm) of the spectral band,
  controlling the width of the Gaussian.

## Value

A numeric vector of relative spectral response values, normalized to the
range 0-1, one per element of `wl`.
