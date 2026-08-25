# get.e2phot calculates the number of moles of photons

`get.e2phot` calculates the number of moles of photons corresponding to
E Joules of energy of wavelength lambda (m)

## Usage

``` r
get.e2phot(lambda, E, constants)
```

## Arguments

- lambda:

  numeric vector. Wavelength(s) of the photons (m).

- E:

  numeric vector. Radiant energy at each wavelength (J), same length as
  `lambda`.

- constants:

  data.frame. Physical constants table (as used throughout SCOPEinR),
  must contain rows for `h` (Planck's constant, J s), `c` (speed of
  light, m s^-1) and `A` (Avogadro's number, mol^-1) in columns
  `constant`/`value`.

## Value

numeric vector. Number of moles of photons (mol) corresponding to `E`
Joules of energy at each wavelength `lambda`.

## Author

     Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)

Carlos Camino (Ported version into R)

## Examples

``` r
if (FALSE) { # \dontrun{
molphotons <- get.e2phot(lambda, E, constants)
} # }
```
