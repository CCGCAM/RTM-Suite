# Planck function

`get.Planck` estimates the spectral radiance emitted by a blackbody (or
greybody, if `em` is supplied) at given wavelength(s) and
temperature(s), following Planck's radiation law.

## Usage

``` r
get.Planck(wl, Tb, em = NULL)
```

## Arguments

- wl:

  numeric vector. Wavelength(s) of interest, in nanometers (nm).

- Tb:

  numeric vector. Temperature(s) of the emitting object, in Kelvin (K).

- em:

  numeric vector, optional. Emissivity of the object at each
  wavelength/temperature (dimensionless, 0-1). If missing, emissivity is
  set to 1 (ideal blackbody) for every element of `Tb`.

## Value

numeric vector. Spectral radiance `Lb` emitted at each
wavelength/temperature pair (W m^-2 sr^-1 um^-1).

## Author

    Christiaan van der Tol (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
Lb <- get.Planck(wl, Tb)
} # }
```
